import argparse
import json
import os
import sqlite3
import pandas as pd
import sys
import glob

# Import shared folding engine utilities
sys.path.append(os.path.dirname(__file__))
from folding_engine import get_chain_sequence, prepare_msa, predict_structure, calculate_ca_rmsd
from energy_engine import enabled as energy_enabled, score_interfaces
from scoring import combined_f_ortho, coherence_report, energy_fields, f_iptm_rel, iptm_margin, ok
from dockq import calculate_dockq
from design_utils import (diff_positions, interface_columns, load_chain, load_config, match_residues, interface_stats,
                          ligand_rmsd_after_receptor_fit, residue_columns, std_residues, write_complex)

NAN = float('nan')


def self_consistency(design_pdb, pred_pdb, chain_A, chain_B):
    """How well RF3's predicted complex reproduces the designed complex geometry.
    Returns (ligand-RMSD of B after fitting on A, #interface contacts in the prediction)."""
    if not pred_pdb or not os.path.exists(pred_pdb):
        return NAN, NAN
    try:
        dA, dB = load_chain(design_pdb, chain_A), load_chain(design_pdb, chain_B)
        pA = load_chain(pred_pdb, chain_A)
        pB = load_chain(pred_pdb, chain_B)
        l_rms = ligand_rmsd_after_receptor_fit(match_residues(dA, pA), match_residues(dB, pB))
        contacts, _ = interface_stats(pA, pB)
        return (NAN if l_rms is None else l_rms), float(contacts)
    except Exception as e:
        print(f"  Warning: self-consistency failed: {e}")
        return NAN, NAN


def main():
    parser = argparse.ArgumentParser(description="Module 5: Final Validation & Orthogonality Scoring via RoseTTAFold-3")
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--analysis_dir', default='results/01_analysis')
    parser.add_argument('--design_dir', default='results/04_rescue_design')
    parser.add_argument('--filter_dir', default=None, help="Directory containing Module 3 fail-fast outputs")
    parser.add_argument('--out_dir', default='results/05_final_eval')
    args = parser.parse_args()

    config = load_config(args.config)

    os.makedirs(args.out_dir, exist_ok=True)
    pcfg = config.get('pipeline', {})
    chain_A_id, chain_B_id = pcfg.get('chain_A', 'A'), pcfg.get('chain_B', 'B')
    wt_pdb = pcfg['input_pdb']
    wt_a3m_A, wt_a3m_B = pcfg['input_msa_A'], pcfg['input_msa_B']

    th = config['thresholds']
    iptm_rescue_min = th['positive_design_rescue']['iptm_rescue_min']
    iptm_rupture_max = th['negative_design_rupture']['iptm_rupture_max']
    iptm_negative_max = th.get('negative_design_orthogonality', {}).get('iptm_negative_max', iptm_rupture_max)
    rescue_rel = th['positive_design_rescue'].get('relative_to_ceiling')       # e.g. 0.85
    f_ortho_min = th.get('orthogonality', {}).get('f_ortho_min', 0.30)
    early_exit = float(th.get('early_exit', {}).get('rescue_iptm_floor', 0.0))
    fcfg = config.get('folding', {})
    folding_engine = fcfg.get('engine', 'rf3')
    neg_regime = fcfg.get('negative_msa_regime', 'matched')      # 'matched' | 'native'
    do_ceiling = bool(fcfg.get('wt_ceiling_control', True))
    ecfg = config.get('energy', {})
    energy_on = energy_enabled(config)
    w_iptm = float(ecfg.get('weight_iptm', 0.5))
    min_binding = float(ecfg.get('min_binding_fraction', 0.5))
    f_energy_min = float(ecfg.get('f_energy_min', 0.0))
    energy_jobs = [{"name": "native", "pdb": wt_pdb}] if energy_on else []

    print(f"Module 5: Final validation with {folding_engine.upper()} "
          f"(mask_mode={fcfg.get('msa_mask_mode', 'gap')}, scope={fcfg.get('msa_mask_scope', 'diff')}, "
          f"negative regime={neg_regime}, ceiling control={do_ceiling}, PyRosetta energy={energy_on})")

    # Redesign-allowed residue ids (used only when msa_mask_scope == 'mutable')
    with open(os.path.join(args.analysis_dir, 'index_mapping.json')) as f:
        mapping = json.load(f)
    fixed_b_path = os.path.join(args.analysis_dir, 'mpnn_fixed_positions_B.json')
    fixed_b = set(json.load(open(fixed_b_path)).get(chain_B_id, [])) if os.path.exists(fixed_b_path) else set()
    iface_A, iface_B = interface_columns(wt_pdb, chain_A_id, chain_B_id, mapping, fixed_b)
    # legacy 'mutable' scope: every position that was left redesignable
    cols_A = residue_columns(wt_pdb, chain_A_id, mapping.get('neighborhood_ids', []))
    cols_B = [i for i, r in enumerate(std_residues(load_chain(wt_pdb, chain_B_id))) if r.id[1] not in fixed_b]

    seq_A_wt = get_chain_sequence(wt_pdb, chain_id=chain_A_id)
    seq_B_wt = get_chain_sequence(wt_pdb, chain_id=chain_B_id)
    if not seq_A_wt or not seq_B_wt:
        print(f"Error: cannot read WT sequences from {wt_pdb}")
        sys.exit(1)

    # 1. Designed complexes
    b_prime_pdbs = sorted(glob.glob(os.path.join(args.design_dir, "*_B_cand_*.pdb")))
    if not b_prime_pdbs:
        sys.exit(f"Error: No B' design PDBs found in {args.design_dir}")

    meta_pairs = {}
    meta_path = os.path.join(args.design_dir, "diversity_pairs_metadata.json")
    if os.path.exists(meta_path):
        meta_pairs = json.load(open(meta_path))

    filter_dir = args.filter_dir or os.path.join(os.path.dirname(os.path.abspath(args.out_dir)), "03_fail_fast")
    if not os.path.exists(filter_dir):
        filter_dir = "results/03_fail_fast"

    results_data = []
    for b_prime_pdb in b_prime_pdbs:
        design_id = os.path.splitext(os.path.basename(b_prime_pdb))[0]
        meta = meta_pairs.get(design_id, {})
        parent_a = meta.get("parent_a", "") or (design_id.split("_B_")[0] if "_B_" in design_id else "")
        print(f"\nModule 5: Evaluating {design_id} (parent A' motif: {parent_a or 'primary'})...")

        seq_A_prime = get_chain_sequence(b_prime_pdb, chain_id=chain_A_id)
        seq_B_prime = get_chain_sequence(b_prime_pdb, chain_id=chain_B_id)
        if not seq_A_prime or not seq_B_prime:
            print(f"  Skipping {design_id}: missing chain sequence(s) in design PDB.")
            continue
        if len(seq_A_prime) != len(seq_A_wt) or len(seq_B_prime) != len(seq_B_wt):
            print(f"  WARNING: length differs from WT (A {len(seq_A_prime)}/{len(seq_A_wt)}, "
                  f"B {len(seq_B_prime)}/{len(seq_B_wt)}); MSA regime may fall back to single-sequence.")
        n_mut_A = len(diff_positions(seq_A_wt, seq_A_prime)) if len(seq_A_wt) == len(seq_A_prime) else -1
        n_mut_B = len(diff_positions(seq_B_wt, seq_B_prime)) if len(seq_B_wt) == len(seq_B_prime) else -1

        pos_out = os.path.join(args.out_dir, "complex_rescue", design_id)
        neg_out = os.path.join(args.out_dir, "complex_negative", design_id)
        rup_out = os.path.join(args.out_dir, "complex_rupture", design_id if neg_regime == 'matched' else (parent_a or design_id))
        ctl_out = os.path.join(args.out_dir, "complex_wt_ceiling", design_id)
        for d in (pos_out, neg_out, rup_out, ctl_out):
            os.makedirs(d, exist_ok=True)

        # MSAs: each designed chain is masked where IT differs from WT
        msa_A_prime = os.path.join(pos_out, f"{design_id}_A.a3m")
        msa_B_prime = os.path.join(pos_out, f"{design_id}_B.a3m")
        st_a = prepare_msa(wt_a3m_A, seq_A_prime, msa_A_prime, config, modified_indices=cols_A, interface_columns=iface_A)
        st_b = prepare_msa(wt_a3m_B, seq_B_prime, msa_B_prime, config, modified_indices=cols_B, interface_columns=iface_B)
        print(f"  Mutations vs WT: A'={n_mut_A}, B'={n_mut_B} | masked cols: A={st_a['n_masked']}, B={st_b['n_masked']}")

        # 2. Positive design: A' + B'
        m_pos = predict_structure([('A_prime', seq_A_prime, msa_A_prime), ('B_prime', seq_B_prime, msa_B_prime)], pos_out, config)
        iptm_rescue = m_pos.get("iptm", NAN)
        print(f"  Positive (A'+B') iPTM: {iptm_rescue:.3f}  (min {iptm_rescue_min})")

        iptm_negative = iptm_rupture = iptm_ceiling = NAN
        status = "ok"
        if iptm_rescue < early_exit:
            status = "early_exit_low_rescue"
            print(f"  Early exit: rescue iPTM < {early_exit} -> negative/rupture tests skipped.")
        else:
            # 3. Negative design A_WT + B' (same information regime as the positive test when 'matched')
            if neg_regime == 'matched':
                msa_A_wt = os.path.join(neg_out, f"{design_id}_Awt_matched.a3m")
                prepare_msa(wt_a3m_A, seq_A_wt, msa_A_wt, config, reference_sequence=seq_A_prime, interface_columns=iface_A)
            else:
                msa_A_wt = wt_a3m_A
            m_neg = predict_structure([('A_wt', seq_A_wt, msa_A_wt), ('B_prime', seq_B_prime, msa_B_prime)], neg_out, config)
            iptm_negative = m_neg.get("iptm", NAN)

            # 4. Rupture A' + B_WT (recomputed here so it shares the regime of the other tests)
            if neg_regime == 'matched':
                msa_B_wt = os.path.join(rup_out, f"{design_id}_Bwt_matched.a3m")
                prepare_msa(wt_a3m_B, seq_B_wt, msa_B_wt, config, reference_sequence=seq_B_prime, interface_columns=iface_B)
            else:
                msa_B_wt = wt_a3m_B
            m_rup = predict_structure([('A_prime', seq_A_prime, msa_A_prime), ('B_wt', seq_B_wt, msa_B_wt)], rup_out, config)
            iptm_rupture = m_rup.get("iptm", NAN)

            # 5. Ceiling: WT complex under the very same masked columns = best achievable iPTM in this regime
            if do_ceiling:
                msa_A_c = os.path.join(ctl_out, f"{design_id}_Awt_ceiling.a3m")
                msa_B_c = os.path.join(ctl_out, f"{design_id}_Bwt_ceiling.a3m")
                prepare_msa(wt_a3m_A, seq_A_wt, msa_A_c, config, reference_sequence=seq_A_prime, interface_columns=iface_A)
                prepare_msa(wt_a3m_B, seq_B_wt, msa_B_c, config, reference_sequence=seq_B_prime, interface_columns=iface_B)
                m_ctl = predict_structure([('A_wt', seq_A_wt, msa_A_c), ('B_wt', seq_B_wt, msa_B_c)], ctl_out, config)
                iptm_ceiling = m_ctl.get("iptm", NAN)

        # 6. Module 3 monomer metrics for the parent A' (NaN when unavailable - never invented)
        plddt_a_prime = rmsd_a_prime_m3 = NAN
        pm = os.path.join(filter_dir, f"metrics_{parent_a}.json") if parent_a else None
        if pm and os.path.exists(pm):
            md = json.load(open(pm))
            plddt_a_prime = md.get("plddt_monomer", NAN)
            rmsd_a_prime_m3 = md.get("rmsd_monomer", NAN)

        # 7. Geometry
        rmsd_a_prime = calculate_ca_rmsd(wt_pdb, b_prime_pdb, chain_ref=chain_A_id, chain_pred=chain_A_id)
        rmsd_b_prime = calculate_ca_rmsd(wt_pdb, b_prime_pdb, chain_ref=chain_B_id, chain_pred=chain_B_id)
        sc_lrms, pred_contacts = self_consistency(b_prime_pdb, m_pos.get("model_pdb"), chain_A_id, chain_B_id)

        dockq_res = calculate_dockq(wt_pdb, pos_out, chain_A=chain_A_id, chain_B=chain_B_id)
        dockq_des = calculate_dockq(b_prime_pdb, pos_out, chain_A=chain_A_id, chain_B=chain_B_id)

        # 8. Complexes for the interface energy (same frame: A' was fitted onto the WT frame in Module 4)
        if energy_on:
            edir = os.path.join(args.out_dir, "energy", design_id)
            os.makedirs(edir, exist_ok=True)
            p_awt_bp, p_ap_bwt = os.path.join(edir, "AWT.Bp.pdb"), os.path.join(edir, "Ap.BWT.pdb")
            write_complex(wt_pdb, b_prime_pdb, p_awt_bp, chain_A_id, chain_B_id)
            write_complex(b_prime_pdb, wt_pdb, p_ap_bwt, chain_A_id, chain_B_id)
            energy_jobs += [{"name": f"{design_id}__Ap.Bp", "pdb": b_prime_pdb},
                            {"name": f"{design_id}__AWT.Bp", "pdb": p_awt_bp},
                            {"name": f"{design_id}__Ap.BWT", "pdb": p_ap_bwt}]

        f_raw = iptm_margin(iptm_rescue, iptm_rupture, iptm_negative)
        f_rel = f_iptm_rel(iptm_rescue, iptm_rupture, iptm_negative, iptm_ceiling)
        rel = (iptm_rescue / iptm_ceiling) if iptm_ceiling and iptm_ceiling == iptm_ceiling and iptm_ceiling > 0 else NAN
        rescue_thr = iptm_rescue_min
        if rescue_rel and iptm_ceiling == iptm_ceiling:
            rescue_thr = min(iptm_rescue_min, rescue_rel * iptm_ceiling)
        pass_rescue = bool(iptm_rescue >= rescue_thr)
        pass_rupture = bool(iptm_rupture <= iptm_rupture_max) if iptm_rupture == iptm_rupture else False
        pass_negative = bool(iptm_negative <= iptm_negative_max) if iptm_negative == iptm_negative else False

        print(f"  --> RF3: rescue {iptm_rescue:.3f} | rupture {iptm_rupture:.3f} | negative {iptm_negative:.3f} | "
              f"ceiling {iptm_ceiling:.3f} | F_iptm(rel) {f_rel:.3f}")
        print(f"  --> B' RMSD vs WT {rmsd_b_prime:.2f} A | self-consistency L-RMSD {sc_lrms:.2f} A | "
              f"DockQ(WT) {dockq_res['dockq']:.3f} ({dockq_res['quality']}) DockQ(design) {dockq_des['dockq']:.3f}")

        results_data.append({
            "design_id": design_id, "parent_a_motif": parent_a or "primary", "folding_engine": folding_engine,
            "status": status,
            "iptm_rescue": iptm_rescue, "iptm_rupture": iptm_rupture, "iptm_negative": iptm_negative,
            "iptm_ceiling": iptm_ceiling, "iptm_rescue_rel": rel, "f_ortho_iptm": f_raw, "f_iptm_rel": f_rel,
            "pass_rescue": pass_rescue, "pass_rupture": pass_rupture, "pass_negative": pass_negative,
            "iptm_rescue_mean": m_pos.get("iptm_mean", NAN), "iptm_rescue_std": m_pos.get("iptm_std", NAN),
            "pae_min_rescue": m_pos.get("chain_pair_pae_min") if m_pos.get("chain_pair_pae_min") is not None else NAN,
            "has_clash": m_pos.get("has_clash", False),
            "dockq": dockq_res["dockq"], "dockq_quality": dockq_res["quality"], "fnat": dockq_res["fnat"],
            "irms": dockq_res["irms"], "lrms": dockq_res["lrms"], "dockq_vs_design": dockq_des["dockq"],
            "selfcons_lrms": sc_lrms, "pred_interface_contacts": pred_contacts,
            "plddt_rescue": m_pos.get("plddt", NAN), "plddt_a_prime": plddt_a_prime,
            "rmsd_a_prime": rmsd_a_prime, "rmsd_b_prime": rmsd_b_prime, "rmsd_a_prime_monomer": rmsd_a_prime_m3,
            "n_mut_A": n_mut_A, "n_mut_B": n_mut_B,
            "scaffold_idx": meta.get("scaffold_idx", NAN), "scaffold_flex_rmsd": meta.get("scaffold_flex_rmsd", NAN),
            "proxy_contacts_a_wt": meta.get("contacts_a_wt", NAN), "proxy_clashes_a_wt": meta.get("clashes_a_wt", NAN),
            "ranking_score": m_pos.get("ranking_score", NAN),
        })

    if not results_data:
        print("Module 5: No design pairs could be evaluated.")
        sys.exit(0)

    # ---- PyRosetta interface energies (one parallel batch for all designs) and the combined score
    energy_cols = ("dG_rescue", "dG_negative", "dG_rupture", "dG_native", "dsasa_rescue", "unsat_hb_rescue", "b_rescue",
                   "b_negative", "b_rupture", "gap_negative", "gap_rupture", "f_energy")
    scores = score_interfaces(energy_jobs, os.path.join(args.out_dir, "energy_scores.json"), config) if energy_on else {}

    def metric(name, key="dG"):
        return scores.get(name, {}).get(key, NAN)

    dG_native = metric("native")
    if energy_on and not (ok(dG_native) and dG_native < 0):
        print("  WARNING: the native complex has no binding energy under this protocol -> energy scores are NaN.")
    for row in results_data:
        did = row["design_id"]
        for k in energy_cols:
            row[k] = NAN
        row["pass_energy"] = NAN
        if energy_on:
            row.update(dG_rescue=metric(f"{did}__Ap.Bp"), dG_negative=metric(f"{did}__AWT.Bp"), dG_rupture=metric(f"{did}__Ap.BWT"),
                       dG_native=dG_native, dsasa_rescue=metric(f"{did}__Ap.Bp", "dSASA"), unsat_hb_rescue=metric(f"{did}__Ap.Bp", "unsat_hb"))
            row.update(energy_fields(row["dG_rescue"], row["dG_negative"], row["dG_rupture"], dG_native))
            row["pass_energy"] = bool(ok(row["b_rescue"], row["f_energy"]) and row["b_rescue"] >= min_binding and row["f_energy"] >= f_energy_min)
        f_i = row["f_iptm_rel"] if ok(row["f_iptm_rel"]) else row["f_ortho_iptm"]
        row["f_ortho"] = combined_f_ortho(f_i, row["f_energy"], w_iptm, energy_on)
        energy_ok = row["pass_energy"] if energy_on else True
        row["passes"] = bool(row["status"] == "ok" and row["pass_rescue"] and row["pass_rupture"] and row["pass_negative"]
                             and energy_ok and ok(row["f_ortho"]) and row["f_ortho"] >= f_ortho_min)
        print(f"  {did}: F_iptm(rel) {row['f_iptm_rel']:.3f} | F_energy {row['f_energy']:.3f} | F_ortho {row['f_ortho']:.3f} | PASS={row['passes']}")

    df = pd.DataFrame(results_data).sort_values(by=["passes", "f_ortho"], ascending=[False, False], na_position='last')

    ceilings = df["iptm_ceiling"].dropna()
    if len(ceilings) and ceilings.median() < 0.5:
        print(f"\n!! DIAGNOSTIC: median WT ceiling iPTM is {ceilings.median():.2f}. Even the NATIVE complex is not "
              f"recognised under this MSA regime -> absolute iPTM thresholds are unreachable; "
              f"try msa_mask_scope: diff / msa_mask_mode: substitute, or rely on rescue_rel & structural metrics.")

    if energy_on:
        coh = coherence_report(df)
        with open(os.path.join(args.out_dir, "coherence.json"), 'w') as f:
            json.dump(coh, f, indent=2)
        print("\nRF3 vs PyRosetta coherence over %d designs (Spearman): rescue %s | cross %s | rupture %s | margins %s | sign agreement %s" % (
            coh["n_designs"], *(("%.2f" % coh[k]) if coh[k] is not None else "n/a" for k in (
                "spearman_rescue_iptm_vs_binding", "spearman_negative_iptm_vs_binding", "spearman_rupture_iptm_vs_binding",
                "spearman_margins", "margin_sign_agreement"))))
        if coh["disagreements"]:
            print("  designs where the two measures disagree most:", ", ".join(coh["disagreements"]))

    csv_path = os.path.join(args.out_dir, "orthogonality_scores.csv")
    df.to_csv(csv_path, index=False)
    print(f"Module 5: Saved rankings CSV to {csv_path}")

    conn = sqlite3.connect(os.path.join(args.out_dir, "results.db"))
    df.to_sql("scores", conn, if_exists="replace", index=False)
    conn.close()

    xlsx_path = os.path.join(args.out_dir, "orthogonality_scores.xlsx")
    html_path = os.path.join(args.out_dir, "orthogonality_scores.html")
    try:
        from export_styled_reports import generate_styled_excel, generate_interactive_html
        generate_styled_excel(csv_path, xlsx_path)
        generate_interactive_html(csv_path, html_path)
        parent_run = os.path.dirname(os.path.abspath(args.out_dir))
        if os.path.exists(parent_run) and os.path.basename(args.out_dir) == "05_final_eval":
            generate_styled_excel(csv_path, os.path.join(parent_run, "SUMMARY.xlsx"))
            generate_interactive_html(csv_path, os.path.join(parent_run, "SUMMARY.html"))
    except Exception as e:
        print(f'Warning: Failed to generate styled reports: {e}')
    print(f"Module 5 Complete: {int(df['passes'].sum())}/{len(df)} designs pass all orthogonality criteria.")


if __name__ == "__main__":
    main()
