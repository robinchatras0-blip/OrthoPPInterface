import yaml
import os

configfile: "config.yaml"

py = config['pipeline'].get('python_bin', 'python')
run_name = config['pipeline'].get('run_name', 'run_01_baseline')
runs_dir = config['pipeline'].get('runs_dir', 'runs')
run_path = os.path.join(runs_dir, run_name).replace('\\', '/')

rule all:
    input:
        f"{run_path}/05_final_eval/orthogonality_scores.csv",
        f"{run_path}/05_final_eval/results.db"

rule analyze_interface:
    input:
        pdb = config['pipeline']['input_pdb']
    output:
        contigs = f"{run_path}/01_analysis/inputs.json",
        fixed_pos = f"{run_path}/01_analysis/mpnn_fixed_positions.json",
        fixed_pos_B = f"{run_path}/01_analysis/mpnn_fixed_positions_B.json",
        bias = f"{run_path}/01_analysis/mpnn_bias.json",
        mapping = f"{run_path}/01_analysis/index_mapping.json"
    threads: 2
    resources:
        mem_mb = 8000,
        runtime = 15
    shell:
        "{py} src/01_analyze_interface.py --config config.yaml --out_dir {run_path}/01_analysis"

rule generate_A_prime:
    input:
        contigs = f"{run_path}/01_analysis/inputs.json",
        fixed_pos = f"{run_path}/01_analysis/mpnn_fixed_positions.json"
    output:
        design_pdb = f"{run_path}/02_rupture_design/A_prime_candidate_00.pdb"
    threads: config['slurm_defaults']['cpus_per_task']
    resources:
        mem_mb = config['slurm_defaults']['mem_gb'] * 1024,
        gpu = 1,
        runtime = 60
    shell:
        "{py} src/02_generate_A_prime.py --config config.yaml --analysis_dir {run_path}/01_analysis --out_dir {run_path}/02_rupture_design"

rule filter_A_prime:
    input:
        design_pdb = f"{run_path}/02_rupture_design/A_prime_candidate_00.pdb"
    output:
        passed_list = f"{run_path}/03_fail_fast/passed_candidates.txt"
    threads: config['slurm_defaults']['cpus_per_task']
    resources:
        mem_mb = config['slurm_defaults']['mem_gb'] * 1024,
        gpu = 1,
        runtime = 60
    shell:
        "{py} src/03_filter_A_prime.py --config config.yaml --design_dir {run_path}/02_rupture_design --analysis_dir {run_path}/01_analysis --out_dir {run_path}/03_fail_fast"

rule generate_B_prime:
    input:
        passed_list = f"{run_path}/03_fail_fast/passed_candidates.txt",
        fixed_pos_B = f"{run_path}/01_analysis/mpnn_fixed_positions_B.json"
    output:
        rescue_pdb = f"{run_path}/04_rescue_design/B_prime_rfd_mpnn.pdb"
    threads: config['slurm_defaults']['cpus_per_task']
    resources:
        mem_mb = config['slurm_defaults']['mem_gb'] * 1024,
        gpu = 1,
        runtime = 60
    shell:
        "{py} src/04_generate_B_prime.py --config config.yaml --analysis_dir {run_path}/01_analysis --passed_candidates {run_path}/03_fail_fast/passed_candidates.txt --out_dir {run_path}/04_rescue_design"

rule eval_final:
    input:
        rescue_pdb = f"{run_path}/04_rescue_design/B_prime_rfd_mpnn.pdb"
    output:
        csv = f"{run_path}/05_final_eval/orthogonality_scores.csv",
        db = f"{run_path}/05_final_eval/results.db"
    threads: config['slurm_defaults']['cpus_per_task']
    resources:
        mem_mb = config['slurm_defaults']['mem_gb'] * 1024,
        gpu = 1,
        runtime = 60
    shell:
        "{py} src/05_eval_final.py --config config.yaml --analysis_dir {run_path}/01_analysis --design_dir {run_path}/04_rescue_design --out_dir {run_path}/05_final_eval"
