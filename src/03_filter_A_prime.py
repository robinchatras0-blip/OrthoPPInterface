import argparse
import glob
import json
import os
import sys

sys.path.append(os.path.dirname(__file__))
from design_utils import interface_columns, load_chain, load_config, std_residues
from folding_engine import get_chain_sequence, prepare_msa, predict_structures, calculate_ca_rmsd


def main():
    parser = argparse.ArgumentParser(description="Module 3: Fail-Fast Filtering (A' validation via RoseTTAFold-3)")
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--design_dir', default='results/02_rupture_design')
    parser.add_argument('--analysis_dir', default='results/01_analysis')
    parser.add_argument('--out_dir', default='results/03_fail_fast')
    args = parser.parse_args()

    config = load_config(args.config)
    os.makedirs(args.out_dir, exist_ok=True)
    scope = config.get('folding', {}).get('msa_mask_scope', 'interface')   # recorded so that Module 5 can reuse results

    pcfg = config['pipeline']
    chain_A_id, chain_B_id = pcfg.get('chain_A', 'A'), pcfg.get('chain_B', 'B')
    wt_pdb = pcfg['input_pdb']
    wt_a3m_A, wt_a3m_B = pcfg['input_msa_A'], pcfg['input_msa_B']

    th = config['thresholds']
    plddt_min = th['monomer_stability']['plddt_min']
    rmsd_max = th['monomer_stability'].get('rmsd_max', 2.0)
    iptm_rupture_max = th['negative_design_rupture']['iptm_rupture_max']
    max_passing = th.get('fail_fast', {}).get('max_passing_candidates')

    with open(os.path.join(args.analysis_dir, 'index_mapping.json')) as f:
        mapping = json.load(f)
    fixed_ids = mapping.get('fixed_ids_A', [])

    seq_A_wt = get_chain_sequence(wt_pdb, chain_A_id)
    seq_B_wt = get_chain_sequence(wt_pdb, chain_B_id)
    if not seq_A_wt or not seq_B_wt:
        sys.exit(f"Error: could not read WT sequences from {wt_pdb}")

    # Residue ids -> 0-based alignment columns (robust to PDB numbering that does not start at 1)
    neighborhood = set(mapping.get('neighborhood_ids', []))
    modified_cols = [i for i, r in enumerate(std_residues(load_chain(wt_pdb, chain_A_id))) if r.id[1] in neighborhood]

    with open(os.path.join(args.analysis_dir, 'mpnn_fixed_positions_B.json')) as f:
        fixed_b_ids = json.load(f).get(chain_B_id, [])
    cols_A, cols_B = interface_columns(wt_pdb, chain_A_id, chain_B_id, mapping, fixed_b_ids)

    # Sanity checks of the predictor: the native complex must be recognised with full MSAs, and we record the
    # ceiling of the configured MSA regime (WT/WT with the interface masked), the best any design can reach there.
    reg_dir = os.path.join(args.out_dir, "wt_control_regime")
    os.makedirs(reg_dir, exist_ok=True)
    msa_A_reg, msa_B_reg = os.path.join(reg_dir, "A.a3m"), os.path.join(reg_dir, "B.a3m")
    prepare_msa(wt_a3m_A, seq_A_wt, msa_A_reg, config, interface_columns=cols_A)
    prepare_msa(wt_a3m_B, seq_B_wt, msa_B_reg, config, interface_columns=cols_B)
    m_ctrl, m_reg = predict_structures([([('A_wt', seq_A_wt, wt_a3m_A), ('B_wt', seq_B_wt, wt_a3m_B)],
                                          os.path.join(args.out_dir, "wt_control_full_msa")),
                                         ([('A_wt', seq_A_wt, msa_A_reg), ('B_wt', seq_B_wt, msa_B_reg)], reg_dir)], config)
    with open(os.path.join(args.out_dir, "wt_control.json"), 'w') as f:
        json.dump({"iptm_full_msa": m_ctrl.get("iptm"), "iptm_regime": m_reg.get("iptm"), "ipsae_regime": m_reg.get("ipsae"),
                   "plddt": m_ctrl.get("plddt"), "msa_scope": scope}, f, indent=2)
    print(f"  WT control (A_wt + B_wt): iPTM = {m_ctrl.get('iptm', 0.0):.3f} with full MSAs, "
          f"{m_reg.get('iptm', 0.0):.3f} in the configured MSA regime (= ceiling for designs)")
    if m_ctrl.get("iptm", 0.0) < 0.6:
        print("  WARNING: the predictor does not recognise the native complex even with full MSAs "
              "-> check MSA files / chain order / RF3 install before trusting any downstream score.")
    if m_reg.get("iptm", 0.0) < 0.4:
        print("  WARNING: the MSA regime leaves too little signal (WT ceiling < 0.4): iPTM cannot discriminate designs.")

    design_pdbs = sorted(glob.glob(os.path.join(args.design_dir, "A_prime_candidate_*.pdb")))
    if not design_pdbs:
        sys.exit(f"Error: No design PDBs found in {args.design_dir}")
    print(f"Module 3: Found {len(design_pdbs)} candidate(s) to evaluate.")
    passed_candidates = []

    batch = int(config.get('folding', {}).get('batch_size', 16))    # candidates per RF3 process (early stop between chunks)
    for start in range(0, len(design_pdbs), batch):
        chunk = design_pdbs[start:start + batch]
        # ---- monomer stability of A' (MSA masked like the design), one RF3 process for the chunk
        cands = []
        for design_pdb in chunk:
            basename = os.path.basename(design_pdb)
            cand_name = os.path.splitext(basename)[0]
            seq_A_prime = get_chain_sequence(design_pdb, chain_A_id)
            if not seq_A_prime:
                print(f"  Skipping {basename}: no chain {chain_A_id} sequence.")
                continue
            monomer_out = os.path.join(args.out_dir, "monomer", cand_name)
            os.makedirs(monomer_out, exist_ok=True)
            monomer_a3m = os.path.join(monomer_out, f"{basename}.a3m")
            msa_stats = prepare_msa(wt_a3m_A, seq_A_prime, monomer_a3m, config, modified_indices=modified_cols,
                                    interface_columns=cols_A)
            print(f"  {cand_name}: {msa_stats['n_masked']} columns re-designed ({msa_stats['frac_masked']:.0%}) "
                  f"mode={msa_stats['mode']}, rows={msa_stats['n_rows']}")
            cands.append({"pdb": design_pdb, "name": cand_name, "basename": basename, "seq": seq_A_prime,
                          "monomer_out": monomer_out, "a3m": monomer_a3m, "msa_stats": msa_stats})
        print(f"\nModule 3: monomer folds for candidates {start + 1}-{start + len(chunk)} of {len(design_pdbs)}...")
        monomers = predict_structures([([('A_prime', c["seq"], c["a3m"])], c["monomer_out"]) for c in cands], config)
        survivors = []
        for c, m in zip(cands, monomers):
            plddt = m.get("plddt", 0.0)
            rmsd_monomer = calculate_ca_rmsd(wt_pdb, c["monomer_out"], chain_ref=chain_A_id)
            rmsd_framework = calculate_ca_rmsd(wt_pdb, c["monomer_out"], chain_ref=chain_A_id, subset_res_ids=fixed_ids)
            print(f"  {c['name']}: monomer pLDDT {plddt:.2f} (min {plddt_min}), RMSD vs WT {rmsd_monomer:.3f} A "
                  f"(framework {rmsd_framework:.3f} A, max {rmsd_max:.1f} A)")
            if plddt < plddt_min:
                print(f"    FAIL-FAST: monomer stability too low ({plddt:.2f} < {plddt_min}). Candidate rejected.")
                continue
            if rmsd_monomer > rmsd_max:
                print(f"    FAIL-FAST: monomer RMSD too high ({rmsd_monomer:.2f} > {rmsd_max}). Candidate rejected.")
                continue
            c.update(plddt=plddt, rmsd_monomer=rmsd_monomer, rmsd_framework=rmsd_framework)
            survivors.append(c)

        # ---- rupture A'.B_WT for the survivors, one RF3 process
        reqs = []
        for c in survivors:
            c["rupture_out"] = os.path.join(args.out_dir, "complex_rupture", c["name"])
            os.makedirs(c["rupture_out"], exist_ok=True)
            msa_B_wt = os.path.join(c["rupture_out"], "Bwt.a3m")
            prepare_msa(wt_a3m_B, seq_B_wt, msa_B_wt, config, interface_columns=cols_B)
            reqs.append(([('A_prime', c["seq"], c["a3m"]), ('B_wt', seq_B_wt, msa_B_wt)], c["rupture_out"]))
        if reqs:
            print(f"Module 3: rupture folds for {len(reqs)} stable candidates...")
        for c, m in zip(survivors, predict_structures(reqs, config) if reqs else []):
            iptm_rupture = m.get("iptm", 1.0)
            print(f"  {c['name']}: rupture iPTM {iptm_rupture:.2f} (max {iptm_rupture_max})")
            with open(os.path.join(args.out_dir, f"metrics_{c['name']}.json"), 'w') as f:
                json.dump({
                    "candidate": c["name"], "pdb": c["pdb"], "plddt_monomer": c["plddt"], "rmsd_monomer": c["rmsd_monomer"],
                    "rmsd_framework": c["rmsd_framework"], "iptm_rupture": iptm_rupture,
                    "iptm_rupture_mean": m.get("iptm_mean"), "ipsae_rupture": m.get("ipsae"),
                    "n_masked_columns": c["msa_stats"]["n_masked"],
                    "msa_scope": scope, "passed": bool(iptm_rupture <= iptm_rupture_max and c["rmsd_monomer"] <= rmsd_max)
                }, f, indent=2)
            if iptm_rupture > iptm_rupture_max:
                print(f"    FAIL-FAST: rupture failed, A' binds WT B ({iptm_rupture:.2f} > {iptm_rupture_max}). Candidate rejected.")
                continue
            print(f"    SUCCESS: {c['basename']} passed fail-fast validation.")
            passed_candidates.append(c["pdb"])
            if max_passing and len(passed_candidates) >= int(max_passing):
                break
        if max_passing and len(passed_candidates) >= int(max_passing):
            print(f"\nModule 3: target of {max_passing} passing candidate(s) reached, stopping early.")
            break

    passed_file = os.path.join(args.out_dir, "passed_candidates.txt")
    with open(passed_file, 'w') as f:
        f.writelines(c + "\n" for c in passed_candidates)
    print(f"Module 3 Complete: {len(passed_candidates)} candidate(s) passed, written to {passed_file}.")


if __name__ == "__main__":
    main()
