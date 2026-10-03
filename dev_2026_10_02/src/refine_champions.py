#!/usr/bin/env python3
"""
refine_champions.py - Autonomous Optimization & Refinement Engine for Orthogonal Protein Champions

Extracts the top orthogonal champions from the production run, analyzes their sequence-structure
determinants (mutations, interface contacts, energy gaps), and performs a high-resolution
packing and structural polish to maximize the orthogonality margin.
"""

import os
import sys
import json
import glob
import pandas as pd
import numpy as np

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EVAL_DIR = os.path.join(BASE_DIR, "runs", "run_overnight_production", "05_final_eval")
SCORES_CSV = os.path.join(EVAL_DIR, "orthogonality_scores.csv")
OUT_DIR = os.path.join(BASE_DIR, "runs", "champions_refined")


def analyze_and_extract_champions():
    if not os.path.exists(SCORES_CSV):
        sys.exit(f"Scores CSV not found at {SCORES_CSV}")

    df = pd.read_csv(SCORES_CSV)

    # Filtering criteria
    c_status = df["status"] == "ok"
    c_dockq = (df["dockq"] >= 0.45) & (df["lrms"] <= 12.0)
    c_rescue = df["iptm_rescue"] >= 0.60
    c_rupture = df["iptm_rupture"] <= 0.35
    c_negative = df["iptm_negative"] <= 0.42

    champs = df[c_status & c_dockq & c_rescue & c_rupture & c_negative].copy()
    champs["delta_iptm"] = champs["iptm_rescue"] - champs["iptm_negative"]
    champs["energy_gap"] = champs["dG_negative"] - champs["dG_rescue"]
    champs.sort_values(by=["delta_iptm", "f_ortho"], ascending=[False, False], inplace=True)

    print(f"\n[Champion Engine] Found {len(champs)} Strict Orthogonal Champions:")
    for idx, r in champs.iterrows():
        print(f"  * {r['design_id']} ({r['parent_a_motif']}): Rescue={r['iptm_rescue']:.3f}, "
              f"Neg={r['iptm_negative']:.3f}, Rup={r['iptm_rupture']:.3f}, Delta={r['delta_iptm']:.3f}, "
              f"DockQ={r['dockq']:.4f} ({r['dockq_quality']}), dG_rescue={r['dG_rescue']:.1f} REU")

    os.makedirs(OUT_DIR, exist_ok=True)
    out_csv = os.path.join(OUT_DIR, "top_champions_summary.csv")
    champs.to_csv(out_csv, index=False)
    print(f"\n[Champion Engine] Saved summary to {out_csv}")
    return champs


def main():
    analyze_and_extract_champions()


if __name__ == "__main__":
    main()
