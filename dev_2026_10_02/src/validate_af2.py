import os
import glob
import pandas as pd
from colabdesign import mk_af_model

params_dir = '/home/ommearo/BindCraft/params'
base_dir = '/mnt/c/Users/hommeacaro/Desktop/CodeDuR/OrthoPPInterface/DOSSIER_CHAMPIONS/structures_pdb'
out_csv = '/mnt/c/Users/hommeacaro/Desktop/CodeDuR/OrthoPPInterface/DOSSIER_CHAMPIONS/tableaux_scores/af2_validation_top5.csv'

top5_prefixes = [
    ("01", "A_prime_candidate_01_B_cand_40"),
    ("02", "A_prime_candidate_17_B_cand_29"),
    ("03", "A_prime_candidate_01_B_cand_36"),
    ("04", "A_prime_candidate_13_B_cand_22"),
    ("05", "A_prime_candidate_11_B_cand_33"),
]

print("Initializing AlphaFold2-Multimer v3 (JAX / CUDA)...")
af = mk_af_model(
    protocol='fixbb',
    use_multimer=True,
    use_templates=True,
    model_names=['model_1_multimer_v3'],
    data_dir=params_dir
)

results = []

for rank_str, design_id in top5_prefixes:
    print(f"\n=======================================================")
    print(f"Evaluating Champion {rank_str}: {design_id}")
    print(f"=======================================================")
    
    pdb_rescue = os.path.join(base_dir, f"{rank_str}_{design_id}_design.pdb")
    pdb_rupture = os.path.join(base_dir, f"{rank_str}_{design_id}_Ap_BWT_rupture.pdb")
    pdb_cross = os.path.join(base_dir, f"{rank_str}_{design_id}_AWT_Bp_cross.pdb")
    
    row = {
        'rank': int(rank_str),
        'design_id': design_id
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
        
        row[f'{state_name}_iptm'] = round(float(m.get('i_ptm', 0.0)), 4)
        row[f'{state_name}_plddt'] = round(float(m.get('plddt', 0.0)), 4)
        row[f'{state_name}_ptm'] = round(float(m.get('ptm', 0.0)), 4)
        row[f'{state_name}_rmsd'] = round(float(m.get('rmsd', 0.0)), 2)
        print(f"    -> iPTM: {row[f'{state_name}_iptm']:.3f} | pLDDT: {row[f'{state_name}_plddt']:.3f} | RMSD: {row[f'{state_name}_rmsd']:.2f} A")
        
    row['delta_iptm_af2'] = round(row.get('rescue_iptm', 0.0) - max(row.get('rupture_iptm', 0.0), row.get('cross_iptm', 0.0)), 4)
    print(f"  >>> Delta iPTM (AF2): {row['delta_iptm_af2']:+.3f}")
    results.append(row)

df = pd.DataFrame(results)
df.to_csv(out_csv, index=False)
print(f"\n[DONE] AF2 validation results written to {out_csv}")
print(df[['rank', 'design_id', 'rescue_iptm', 'rupture_iptm', 'cross_iptm', 'delta_iptm_af2', 'rescue_plddt']])
