import argparse
import yaml
import json
import os
import sys
import glob
import shutil

# Import shared folding engine utilities
sys.path.append(os.path.dirname(__file__))
from folding_engine import get_chain_sequence, build_hybrid_msa, build_unpaired_complex_a3m, predict_structure, calculate_ca_rmsd

def main():
    parser = argparse.ArgumentParser(description="Module 3: Fail-Fast Filtering (A' Validation via AF2 / AF3)")
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--design_dir', default='results/02_rupture_design')
    parser.add_argument('--analysis_dir', default='results/01_analysis')
    parser.add_argument('--out_dir', default='results/03_fail_fast')
    args = parser.parse_args()

    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    os.makedirs(args.out_dir, exist_ok=True)
    execution_mode = config.get('pipeline', {}).get('execution_mode', 'mock')
    chain_A_id = config.get('pipeline', {}).get('chain_A', 'A')
    chain_B_id = config.get('pipeline', {}).get('chain_B', 'B')
    wt_pdb = config.get('pipeline', {}).get('input_pdb', 'data/inputs/complex_S1_S2.pdb')
    wt_a3m_A = config.get('pipeline', {}).get('input_msa_A', config.get('pipeline', {}).get('input_msa', 'data/inputs/S1_chain_A.a3m'))
    wt_a3m_B = config.get('pipeline', {}).get('input_msa_B', 'data/inputs/S2_chain_B.a3m')

    plddt_min = config['thresholds']['monomer_stability']['plddt_min']
    rmsd_max = config.get('thresholds', {}).get('monomer_stability', {}).get('rmsd_max', 2.0)
    iptm_rupture_max = config['thresholds']['negative_design_rupture']['iptm_rupture_max']
    folding_engine = config.get('folding', {}).get('engine', 'af3')

    print(f"Module 3: Initializing Fail-Fast Filtering using folding engine: {folding_engine.upper()}")

    mapping_file = os.path.join(args.analysis_dir, 'index_mapping.json')
    modified_indices = []
    fixed_ids = []
    if os.path.exists(mapping_file):
        with open(mapping_file, 'r') as f:
            mapping_data = json.load(f)
            modified_indices = mapping_data.get('neighborhood_ids', [])
            fixed_ids = mapping_data.get('fixed_ids_A', [])

    # Get WT sequence of Chain B
    seq_B_wt = get_chain_sequence(wt_pdb, chain_id=chain_B_id)
    if not seq_B_wt:
        seq_B_wt = "M" * 50  # Fallback dummy sequence if file absent in mock

    design_pdbs = sorted(glob.glob(os.path.join(args.design_dir, "A_prime_candidate_*.pdb")))
    if not design_pdbs:
        design_pdbs = sorted(glob.glob(os.path.join(args.design_dir, "*_mpnn.pdb")))
    if not design_pdbs:
        fallback_pdb = os.path.join(args.design_dir, "A_prime_rfd_mpnn.pdb")
        if execution_mode == 'mock':
            if not os.path.exists(fallback_pdb):
                with open(fallback_pdb, 'w') as f:
                    f.write("DUMMY PDB CONTENT\n")
            design_pdbs = [fallback_pdb]
        else:
            print(f"Error: No design PDBs found in {args.design_dir}")
            sys.exit(1)

    print(f"Module 3: Found {len(design_pdbs)} candidate(s) to evaluate.")
    passed_candidates = []

    for idx, design_pdb in enumerate(design_pdbs):
        basename = os.path.basename(design_pdb)
        cand_name = os.path.splitext(basename)[0]
        print(f"\n[{idx+1}/{len(design_pdbs)}] Module 3: Evaluating candidate {basename}...")

        # Extract actual sequence of redesigned Chain A from PDB
        seq_A_prime = get_chain_sequence(design_pdb, chain_id=chain_A_id)
        if not seq_A_prime:
            seq_A_prime = "M" * 50

        # 1. Monomer Hybrid MSA & Stability Filter (A')
        monomer_out = os.path.join(args.out_dir, "monomer", cand_name)
        if os.path.exists(monomer_out):
            shutil.rmtree(monomer_out)
        os.makedirs(monomer_out, exist_ok=True)
        monomer_a3m = os.path.join(monomer_out, f"{basename}.a3m")
        build_hybrid_msa(wt_a3m_A, seq_A_prime, modified_indices, monomer_a3m)

        metrics_monomer = predict_structure(
            input_pdb=design_pdb,
            fasta_sequences=[('A_prime', seq_A_prime)],
            out_dir=monomer_out,
            msa_path=monomer_a3m,
            config=config
        )

        plddt = metrics_monomer.get("plddt", 0.0)
        rmsd_monomer = calculate_ca_rmsd(wt_pdb, monomer_out, chain_ref=chain_A_id)
        rmsd_framework = calculate_ca_rmsd(wt_pdb, monomer_out, chain_ref=chain_A_id, subset_res_ids=fixed_ids)
        print(f"  Monomer pLDDT: {plddt:.2f} (Threshold min: {plddt_min})")
        print(f"  Monomer RMSD vs WT: {rmsd_monomer:.3f} Å (Framework: {rmsd_framework:.3f} Å, Threshold max: {rmsd_max:.1f} Å)")

        if plddt < plddt_min:
            print(f"  FAIL-FAST: Monomer stability too low ({plddt:.2f} < {plddt_min}). Candidate rejected.")
            continue
        if rmsd_monomer > rmsd_max:
            print(f"  FAIL-FAST: Monomer RMSD too high ({rmsd_monomer:.2f} > {rmsd_max}). Candidate rejected.")
            continue

        # 2. Complex Rupture Test (A' + WT B heterodimer)
        rupture_out = os.path.join(args.out_dir, "complex_rupture", cand_name)
        if os.path.exists(rupture_out):
            shutil.rmtree(rupture_out)
        os.makedirs(rupture_out, exist_ok=True)
        complex_a3m = os.path.join(rupture_out, f"{cand_name}_complex.a3m")
        build_unpaired_complex_a3m(wt_a3m_A, wt_a3m_B, seq_A_prime, seq_B_wt, modified_indices, complex_a3m)

        metrics_rupture = predict_structure(
            input_pdb=design_pdb,
            fasta_sequences=[('A_prime', seq_A_prime), ('B_wt', seq_B_wt)],
            out_dir=rupture_out,
            msa_path=complex_a3m,
            config=config
        )

        iptm_rupture = metrics_rupture.get("iptm", 1.0)
        print(f"  Rupture iPTM: {iptm_rupture:.2f} (Threshold max: {iptm_rupture_max})")

        cand_metrics_path = os.path.join(args.out_dir, f"metrics_{cand_name}.json")
        with open(cand_metrics_path, 'w') as f:
            json.dump({
                "candidate": cand_name,
                "pdb": design_pdb,
                "plddt_monomer": plddt,
                "rmsd_monomer": rmsd_monomer,
                "rmsd_framework": rmsd_framework,
                "iptm_rupture": iptm_rupture,
                "passed": bool(iptm_rupture <= iptm_rupture_max and rmsd_monomer <= rmsd_max)
            }, f, indent=2)

        if iptm_rupture > iptm_rupture_max:
            print(f"  FAIL-FAST: Rupture failed. A' binds WT B too strongly ({iptm_rupture:.2f} > {iptm_rupture_max}). Candidate rejected.")
            continue

        print(f"  SUCCESS: Candidate {basename} passed Fail-Fast validation!")
        passed_candidates.append(design_pdb)

    passed_file = os.path.join(args.out_dir, "passed_candidates.txt")
    with open(passed_file, 'w') as f:
        for candidate in passed_candidates:
            f.write(candidate + "\n")

    print(f"Module 3 Complete: {len(passed_candidates)} candidate(s) passed and written to {passed_file}.")

if __name__ == "__main__":
    main()
