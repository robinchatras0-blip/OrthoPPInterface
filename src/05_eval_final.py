import argparse
import yaml
import json
import os
import sqlite3
import pandas as pd
import sys
import glob

# Import shared folding engine utilities
sys.path.append(os.path.dirname(__file__))
from folding_engine import get_chain_sequence, build_unpaired_complex_a3m, predict_structure, calculate_ca_rmsd

def main():
    parser = argparse.ArgumentParser(description="Module 5: Final Validation & Orthogonality Scoring via AF2 / AF3")
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--analysis_dir', default='results/01_analysis')
    parser.add_argument('--design_dir', default='results/04_rescue_design')
    parser.add_argument('--out_dir', default='results/05_final_eval')
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
    iptm_rescue_min = config['thresholds']['positive_design_rescue']['iptm_rescue_min']
    folding_engine = config.get('folding', {}).get('engine', 'af3')

    print(f"Module 5: Initializing Final Validation using folding engine: {folding_engine.upper()}")

    # Load modified indices for interface masking in unpaired MSAs
    mapping_file = os.path.join(args.analysis_dir, 'index_mapping.json')
    modified_indices_A = []
    if os.path.exists(mapping_file):
        with open(mapping_file, 'r') as f:
            mapping_data = json.load(f)
            modified_indices_A = mapping_data.get('neighborhood_ids', [])

    # Extract WT sequences
    seq_A_wt = get_chain_sequence(wt_pdb, chain_id=chain_A_id)
    if not seq_A_wt:
        seq_A_wt = "M" * 50
    seq_B_wt = get_chain_sequence(wt_pdb, chain_id=chain_B_id)
    if not seq_B_wt:
        seq_B_wt = "M" * 50

    mpnn_fixed_b = os.path.join(args.analysis_dir, 'mpnn_fixed_positions_B.json')
    modified_indices_B = []
    if os.path.exists(mpnn_fixed_b):
        with open(mpnn_fixed_b, 'r') as f:
            fixed_data = json.load(f)
            fixed_b = set(fixed_data.get(chain_B_id, []))
            modified_indices_B = [i for i in range(1, len(seq_B_wt) + 1) if i not in fixed_b]

    b_prime_pdbs = glob.glob(os.path.join(args.design_dir, "*_mpnn.pdb"))
    if not b_prime_pdbs:
        fallback_pdb = os.path.join(args.design_dir, "B_prime_rfd_mpnn.pdb")
        if execution_mode == 'mock':
            if not os.path.exists(fallback_pdb):
                with open(fallback_pdb, 'w') as f:
                    f.write("DUMMY B PRIME PDB CONTENT\n")
            b_prime_pdbs = [fallback_pdb]
        else:
            print(f"Error: No B' design PDBs found in {args.design_dir}")
            sys.exit(1)

    results_data = []

    for b_prime_pdb in b_prime_pdbs:
        basename = os.path.basename(b_prime_pdb)
        design_id = os.path.splitext(basename)[0]
        print(f"Module 5: Evaluating design pair {design_id}...")

        # Extract A' and B' sequences from designed complex PDB
        seq_A_prime = get_chain_sequence(b_prime_pdb, chain_id=chain_A_id)
        seq_B_prime = get_chain_sequence(b_prime_pdb, chain_id=chain_B_id)
        if not seq_A_prime:
            seq_A_prime = "M" * 50
        if not seq_B_prime:
            seq_B_prime = "M" * 50

        # 1. Positive Design Evaluation (A' + B')
        pos_out = os.path.join(args.out_dir, "complex_rescue", design_id)
        if os.path.exists(pos_out):
            import shutil
            shutil.rmtree(pos_out)
        os.makedirs(pos_out, exist_ok=True)
        rescue_a3m = os.path.join(pos_out, f"{design_id}_rescue.a3m")
        build_unpaired_complex_a3m(wt_a3m_A, wt_a3m_B, seq_A_prime, seq_B_prime, modified_indices_A, rescue_a3m, modified_indices_B=modified_indices_B)

        metrics_pos = predict_structure(
            input_pdb=b_prime_pdb,
            fasta_sequences=[('A_prime', seq_A_prime), ('B_prime', seq_B_prime)],
            out_dir=pos_out,
            msa_path=rescue_a3m,
            config=config
        )

        iptm_rescue = metrics_pos.get("iptm", 0.0)
        plddt_rescue = metrics_pos.get("plddt", 0.0)
        ranking_score = metrics_pos.get("ranking_score", 0.0)

        print(f"  Positive Design iPTM: {iptm_rescue:.2f} (Threshold min: {iptm_rescue_min})")

        # 2. Negative Design Evaluation (A_WT + B')
        neg_out = os.path.join(args.out_dir, "complex_negative", design_id)
        if os.path.exists(neg_out):
            import shutil
            shutil.rmtree(neg_out)
        os.makedirs(neg_out, exist_ok=True)
        neg_a3m = os.path.join(neg_out, f"{design_id}_negative.a3m")
        build_unpaired_complex_a3m(wt_a3m_A, wt_a3m_B, seq_A_wt, seq_B_prime, [], neg_a3m, modified_indices_B=modified_indices_B)

        metrics_neg = predict_structure(
            input_pdb=b_prime_pdb,
            fasta_sequences=[('A_wt', seq_A_wt), ('B_prime', seq_B_prime)],
            out_dir=neg_out,
            msa_path=neg_a3m,
            config=config
        )

        iptm_negative = metrics_neg.get("iptm", 1.0)
        plddt_negative = metrics_neg.get("plddt", 0.0)

        # 3. Retrieve Rupture iPTM and A' metrics from Module 3
        rupture_metrics_files = sorted(glob.glob("results/03_fail_fast/metrics_*.json"))
        iptm_rupture = 0.15
        plddt_a_prime = 90.0
        rmsd_a_prime = 0.55
        if rupture_metrics_files:
            with open(rupture_metrics_files[0], 'r') as f:
                metrics_rupture = json.load(f)
            iptm_rupture = metrics_rupture.get("iptm_rupture", 0.15)
            plddt_a_prime = metrics_rupture.get("plddt_monomer", 90.0)
            rmsd_a_prime = metrics_rupture.get("rmsd_monomer", 0.55)

        # 4. Compute RMSD for B' vs WT B
        rmsd_b_prime = calculate_ca_rmsd(wt_pdb, b_prime_pdb, chain_ref=chain_B_id, chain_pred=chain_B_id)

        # Calculate Orthogonality Score F_ortho
        f_ortho = iptm_rescue - max(iptm_rupture, iptm_negative)

        print(f"  --> iPTM Rescue (A' + B'): {iptm_rescue:.2f}")
        print(f"  --> iPTM Negative (A' + B_WT): {iptm_rupture:.2f}")
        print(f"  --> iPTM Negative (A_WT + B'): {iptm_negative:.2f}")
        print(f"  --> Monomer RMSD (A'): {rmsd_a_prime:.3f} Å")
        print(f"  --> Monomer RMSD (B'): {rmsd_b_prime:.3f} Å")
        print(f"  --> Orthogonality Score F_ortho: {f_ortho:.2f}")

        results_data.append({
            "design_id": design_id,
            "folding_engine": folding_engine,
            "iptm_rescue": iptm_rescue,
            "iptm_rupture": iptm_rupture,
            "iptm_negative": iptm_negative,
            "f_ortho": f_ortho,
            "plddt_rescue": plddt_rescue,
            "plddt_a_prime": plddt_a_prime,
            "rmsd_a_prime": rmsd_a_prime,
            "rmsd_b_prime": rmsd_b_prime,
            "ranking_score": ranking_score
        })

    if not results_data:
        print("Module 5: No design pairs passed final validation.")
        sys.exit(0)

    df = pd.DataFrame(results_data)
    df = df.sort_values(by="f_ortho", ascending=False)

    # Export to CSV
    csv_path = os.path.join(args.out_dir, "orthogonality_scores.csv")
    df.to_csv(csv_path, index=False)
    print(f"Module 5: Saved rankings CSV to {csv_path}")

    # Export to SQLite database
    db_path = os.path.join(args.out_dir, "results.db")
    conn = sqlite3.connect(db_path)
    df.to_sql("scores", conn, if_exists="replace", index=False)
    conn.close()
    print(f"Module 5: Saved SQLite database to {db_path}")

    print("Module 5 Complete: Final evaluation and orthogonality scoring complete.")

if __name__ == "__main__":
    main()
