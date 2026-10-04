import os
import glob
import json
import pandas as pd
import numpy as np
from scipy import stats
from colabdesign import mk_af_model

# Force unbuffered output
import sys
sys.stdout.reconfigure(line_buffering=True)

params_dir = '/home/ommearo/BindCraft/params'
base_dir = '/mnt/c/Users/hommeacaro/Desktop/CodeDuR/OrthoPPInterface/DOSSIER_CHAMPIONS/structures_pdb'
summary_csv = '/mnt/c/Users/hommeacaro/Desktop/CodeDuR/OrthoPPInterface/DOSSIER_CHAMPIONS/tableaux_scores/champions_orthogonaux_summary.csv'
out_csv = '/mnt/c/Users/hommeacaro/Desktop/CodeDuR/OrthoPPInterface/DOSSIER_CHAMPIONS/tableaux_scores/af2_validation_19_champions.csv'
corr_csv = '/mnt/c/Users/hommeacaro/Desktop/CodeDuR/OrthoPPInterface/DOSSIER_CHAMPIONS/tableaux_scores/correlation_rf3_af2.csv'

df_champions = pd.read_csv(summary_csv)
print(f"Loaded {len(df_champions)} champions from summary CSV.")

print("Initializing AlphaFold2-Multimer v3 (JAX / CUDA RTX 5070 Ti)...")
af = mk_af_model(
    protocol='fixbb',
    use_multimer=True,
    use_templates=True,
    model_names=['model_1_multimer_v3'],
    data_dir=params_dir
)

results = []

for idx, r in df_champions.iterrows():
    rk = int(r['rank'])
    design_id = r['design_id']
    rank_str = f"{rk:02d}"
    print(f"\n=======================================================")
    print(f"[{rk}/19] Evaluating Champion: {design_id}")
    print(f"=======================================================")
    
    pdb_rescue = os.path.join(base_dir, f"{rank_str}_{design_id}_design.pdb")
    pdb_rupture = os.path.join(base_dir, f"{rank_str}_{design_id}_Ap_BWT_rupture.pdb")
    pdb_cross = os.path.join(base_dir, f"{rank_str}_{design_id}_AWT_Bp_cross.pdb")
    
    row = {
        'rank': rk,
        'design_id': design_id,
        'parent_a_motif': r['parent_a_motif'],
        'dockq_rf3': r['dockq'],
        'lrms_rf3': r['lrms'],
        'iptm_rescue_rf3': r['iptm_rescue'],
        'iptm_rupture_rf3': r['iptm_rupture'],
        'iptm_negative_rf3': r['iptm_negative'],
        'delta_iptm_rf3': r['delta_iptm'],
        'dG_rescue_pyrosetta': r['dG_rescue'],
        'f_ortho_pyrosetta': r['f_ortho'],
    }
    
    states = [
        ('rescue', pdb_rescue),
        ('rupture', pdb_rupture),
        ('cross', pdb_cross),
    ]
    
    for state_name, pdb_path in states:
        if not os.path.exists(pdb_path):
            print(f"Warning: {pdb_path} not found!")
            continue
        print(f"  Predicting state '{state_name}' with AF2...")
        af.prep_inputs(pdb_filename=pdb_path, chain='A,B', rm_template=False)
        af.set_seq(mode='wt')
        af.predict(num_recycles=3, verbose=False)
        m = af.aux['log']
        
        row[f'{state_name}_iptm_af2'] = round(float(m.get('i_ptm', 0.0)), 4)
        row[f'{state_name}_plddt_af2'] = round(float(m.get('plddt', 0.0)), 4)
        row[f'{state_name}_ptm_af2'] = round(float(m.get('ptm', 0.0)), 4)
        row[f'{state_name}_rmsd_af2'] = round(float(m.get('rmsd', 0.0)), 2)
        print(f"    -> iPTM: {row[f'{state_name}_iptm_af2']:.3f} | pLDDT: {row[f'{state_name}_plddt_af2']:.3f} | RMSD: {row[f'{state_name}_rmsd_af2']:.2f} A")
        
    row['delta_iptm_af2'] = round(row.get('rescue_iptm_af2', 0.0) - max(row.get('rupture_iptm_af2', 0.0), row.get('cross_iptm_af2', 0.0)), 4)
    print(f"  >>> Delta iPTM (AF2): {row['delta_iptm_af2']:+.3f} (vs RF3: {row['delta_iptm_rf3']:+.3f})")
    results.append(row)

df_res = pd.DataFrame(results)
df_res.to_csv(out_csv, index=False)
print(f"\n[DONE] All 19 AF2 validation results written to {out_csv}")

# ----------------- Statistical Correlations -----------------
pairs_to_correlate = [
    ("iPTM Rescue (RF3 vs AF2)", "iptm_rescue_rf3", "rescue_iptm_af2"),
    ("iPTM Rupture (RF3 vs AF2)", "iptm_rupture_rf3", "rupture_iptm_af2"),
    ("iPTM Croise (RF3 vs AF2)", "iptm_negative_rf3", "cross_iptm_af2"),
    ("Delta iPTM (RF3 vs AF2)", "delta_iptm_rf3", "delta_iptm_af2"),
    ("DockQ RF3 vs RMSD AF2", "dockq_rf3", "rescue_rmsd_af2"),
    ("PyRosetta dG vs iPTM Rescue AF2", "dG_rescue_pyrosetta", "rescue_iptm_af2"),
]

corr_rows = []
print("\n" + "="*70)
print("CORRELATION ANALYSIS: RoseTTAFold-3 vs AlphaFold2-Multimer v3")
print("="*70)

for label, col_x, col_y in pairs_to_correlate:
    x = df_res[col_x].values
    y = df_res[col_y].values
    
    r_pearson, p_pearson = stats.pearsonr(x, y)
    r_spearman, p_spearman = stats.spearmanr(x, y)
    
    print(f"\nMetric: {label}")
    print(f"  Pearson  r: {r_pearson:+.3f} (p-value: {p_pearson:.3e})")
    print(f"  Spearman rho: {r_spearman:+.3f} (p-value: {p_spearman:.3e})")
    
    corr_rows.append({
        'comparison': label,
        'metric_rf3': col_x,
        'metric_af2': col_y,
        'pearson_r': round(r_pearson, 4),
        'pearson_pvalue': round(p_pearson, 6),
        'spearman_rho': round(r_spearman, 4),
        'spearman_pvalue': round(p_spearman, 6),
    })

df_corr = pd.DataFrame(corr_rows)
df_corr.to_csv(corr_csv, index=False)
print(f"\n[DONE] Correlation matrix written to {corr_csv}")
