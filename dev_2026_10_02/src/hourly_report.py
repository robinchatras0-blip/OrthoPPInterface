#!/usr/bin/env python3
"""Generates an hourly progress report for the overnight production run.
Reads real-time artifacts and outputs a structured Markdown summary.
"""
import os
import sys
import glob
import json
import pandas as pd
from datetime import datetime

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN_DIR = os.path.join(BASE_DIR, "runs", "run_overnight_production")
DESIGN_DIR = os.path.join(RUN_DIR, "04_rescue_design")
EVAL_DIR = os.path.join(RUN_DIR, "05_final_eval")
META_PATH = os.path.join(DESIGN_DIR, "diversity_pairs_metadata.json")
STAGE1_PATH = os.path.join(EVAL_DIR, "rescue_stage1_progress.csv")
FINAL_CSV_PATH = os.path.join(EVAL_DIR, "orthogonality_scores.csv")
HTML_PATH = os.path.join(RUN_DIR, "SUMMARY.html")


def main():
    now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"### Rapport Horaire de Production - {now_str}\n")

    # 1. Module 4 Inspection
    designed_pdbs = glob.glob(os.path.join(DESIGN_DIR, "*_B_cand_*.pdb"))
    meta = json.load(open(META_PATH)) if os.path.exists(META_PATH) else {}
    motifs = set(v.get("parent_a") for v in meta.values() if v.get("parent_a"))
    if len(motifs) >= 12:
        print(f"**Module 4 (Design B')** :  **Terminé** ({len(designed_pdbs)} candidats générés sur 12 motifs A')")
    elif len(motifs) > 0:
        print(f"**Module 4 (Design B')** : ⏳ **En cours** ({len(motifs)}/12 motifs terminés, **{len(designed_pdbs)}/600 candidats B' générés et relaxés**)")
    else:
        motifs_dirs = glob.glob(os.path.join(DESIGN_DIR, "motif_*"))
        print(f"**Module 4 (Design B')** : ⏳ **Initialisation / RFdiffusion** ({len(motifs_dirs)} dossiers motifs en cours...)")

    # 2. Module 5 Stage 1 Inspection
    if os.path.exists(STAGE1_PATH):
        df_s1 = pd.read_csv(STAGE1_PATH)
        n_tot = len(df_s1)
        survivors = df_s1[df_s1["status"] == "ok"]
        n_surv = len(survivors)
        med_dq = survivors["dockq"].median() if n_surv else 0.0
        best_dq = survivors["dockq"].max() if n_surv else 0.0
        best_iptm = survivors["iptm_rescue"].max() if n_surv else 0.0
        print(f"**Module 5 - Étape 1 (Sauvetage RF3)** :  **{n_tot} complexes évalués** | **{n_surv} survivants** à l'interface native ({n_surv/max(1,n_tot)*100:.1f}%)")
        print(f"  - Meilleur DockQ : **{best_dq:.4f}** | Médiane survivants : **{med_dq:.4f}** | Max iPTM sauvetage : **{best_iptm:.3f}**")

        # Top 5 Stage 1 preview
        top5_s1 = survivors.sort_values(by=["dockq", "iptm_rescue"], ascending=[False, False]).head(5)
        print("\n*Top 5 Pré-sélections Étape 1 (Amarrage natif vérifié)* :")
        print("| Design ID | Motif Parent | DockQ | L-RMSD (Å) | iPTM Sauvetage |")
        print("| :--- | :---: | :---: | :---: | :---: |")
        for _, r in top5_s1.iterrows():
            print(f"| `{r['design_id']}` | {r['parent_a']} | **{r['dockq']:.4f}** | {r['lrms']:.2f} Å | **{r['iptm_rescue']:.3f}** |")
    elif os.path.exists(EVAL_DIR):
        c_rescue = glob.glob(os.path.join(EVAL_DIR, "complex_rescue", "*_unrelaxed_rank_001.pdb"))
        print(f"**Module 5 - Étape 1** : ⏳ **Repliement RF3 en cours** ({len(c_rescue)} complexes résolus...)")

    # 3. Final Results Inspection
    if os.path.exists(FINAL_CSV_PATH):
        df = pd.read_csv(FINAL_CSV_PATH)
        passing = df[df["passes"] == True]
        print(f"\n**Module 5 - Évaluation Finale** :  **CALCUL TOTALEMENT TERMINÉ** !")
        print(f"  - **{len(passing)} / {len(df)} designs passent 100% des critères orthogonaux** !")
        print(f"  - Dashboard HTML : [{os.path.basename(HTML_PATH)}](file:///{HTML_PATH.replace(os.sep, '/')})")

        print("\n*Top Champions Orthogonaux Validés* :")
        print("| Rang | Design ID | Motif A' | DockQ | iPTM Rescue | Rupture (A'·B_WT) | Négatif (A_WT·B') | F_ortho |")
        print("| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
        for idx, (_, r) in enumerate(passing.head(5).iterrows()):
            print(f"| #{idx+1} | `{r['design_id']}` | {r['parent_a_motif']} | **{r['dockq']:.4f}** | **{r['iptm_rescue']:.3f}** | {r['iptm_rupture']:.3f} | {r['iptm_negative']:.3f} | **+{r['f_ortho']:.4f}** |")


if __name__ == "__main__":
    main()
