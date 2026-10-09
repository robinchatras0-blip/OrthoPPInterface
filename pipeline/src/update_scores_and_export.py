#!/usr/bin/env python3
"""
update_scores_and_export.py - Synchronize corrected pass criteria and regenerate styled dashboards.
"""

import os
import shutil
import pandas as pd
from export_styled_reports import generate_styled_excel, generate_interactive_html

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN_DIR = os.path.join(BASE_DIR, "runs", "run_overnight_production")
EVAL_DIR = os.path.join(RUN_DIR, "05_final_eval")
CSV_PATH = os.path.join(EVAL_DIR, "orthogonality_scores.csv")


def main():
    if not os.path.exists(CSV_PATH):
        print(f"Error: {CSV_PATH} not found.")
        return

    df = pd.read_csv(CSV_PATH)

    # 1. Correct Pass Flags based on actual measurements
    df["pass_dockq"] = (df["dockq"] >= 0.23) & (df["lrms"] <= 15.0) & (df["status"] == "ok")
    
    # Ceiling-relative rescue threshold (0.80 * ceiling or 0.70)
    rescue_thr = df["iptm_ceiling"].apply(lambda c: min(0.70, 0.80 * c) if pd.notnull(c) and c > 0 else 0.70)
    df["pass_rescue"] = (df["iptm_rescue"] >= rescue_thr) & df["pass_dockq"]
    df["pass_rupture"] = df["iptm_rupture"] <= 0.35
    df["pass_negative"] = df["iptm_negative"] <= 0.55
    df["pass_energy"] = (df["b_rescue"] >= 0.5) & (df["f_energy"] >= -0.05)

    df["passes"] = (
        (df["status"] == "ok") &
        df["pass_rescue"] &
        df["pass_rupture"] &
        df["pass_negative"] &
        df["pass_energy"] &
        df["pass_dockq"] &
        (df["f_ortho"] >= 0.05)
    )

    # Sort: passes first, then dockq, then f_ortho
    df.sort_values(by=["passes", "dockq", "f_ortho"], ascending=[False, False, False], inplace=True)
    df.to_csv(CSV_PATH, index=False)
    print(f"[Update] Saved clean corrected rankings to {CSV_PATH}")
    print(f"  * pass_dockq: {df['pass_dockq'].sum()} / {len(df)}")
    print(f"  * pass_rescue: {df['pass_rescue'].sum()} / {len(df)}")
    print(f"  * pass_rupture (<= 0.35): {df['pass_rupture'].sum()} / {len(df)}")
    print(f"  * pass_negative (<= 0.55): {df['pass_negative'].sum()} / {len(df)}")
    print(f"  * pass_energy: {df['pass_energy'].sum()} / {len(df)}")
    print(f"  * TOTAL PASSES (all criteria): {df['passes'].sum()} / {len(df)}")

    # 2. Regenerate Styled Reports
    xlsx_eval = os.path.join(EVAL_DIR, "orthogonality_scores.xlsx")
    html_eval = os.path.join(EVAL_DIR, "orthogonality_scores.html")
    generate_styled_excel(CSV_PATH, xlsx_eval)
    generate_interactive_html(CSV_PATH, html_eval)

    xlsx_run = os.path.join(RUN_DIR, "SUMMARY.xlsx")
    html_run = os.path.join(RUN_DIR, "SUMMARY.html")
    shutil.copy2(xlsx_eval, xlsx_run)
    shutil.copy2(html_eval, html_run)
    print(f"[Update] Regenerated {html_run} and {xlsx_run} successfully.")


if __name__ == "__main__":
    main()
