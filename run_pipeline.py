import argparse
import yaml
import subprocess
import sys
import os
import shutil
import datetime
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="OrthoInterface: End-to-End Experiment Pipeline Runner")
    parser.add_argument('--config', default='config.yaml', help="Path to config.yaml")
    parser.add_argument('--run_name', default=None, help="Name of the experiment run folder (overrides config)")
    parser.add_argument('--runs_dir', default=None, help="Root directory for runs (overrides config)")
    parser.add_argument('--steps', default='1-5', help="Pipeline steps to run (e.g. 1-5, all, 1,2,3, 4-5)")
    args = parser.parse_args()

    # Load configuration
    if not os.path.exists(args.config):
        print(f"Error: Configuration file {args.config} not found!")
        sys.exit(1)

    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    # Determine run directory
    runs_dir = args.runs_dir or config.get('pipeline', {}).get('runs_dir', 'runs')
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    run_name = args.run_name or config.get('pipeline', {}).get('run_name', f"run_{timestamp}")
    
    run_path = os.path.join(runs_dir, run_name).replace('\\', '/')
    os.makedirs(run_path, exist_ok=True)

    # Save a timestamped copy of the config used for full reproducibility
    config_copy_path = os.path.join(run_path, "config_used.yaml")
    with open(config_copy_path, 'w') as f:
        yaml.dump(config, f, sort_keys=False)

    print("=" * 70)
    print(f"  ORTHOINTERFACE PIPELINE - RUN: {run_name}")
    print(f"  Output Directory: {run_path}")
    print(f"  Configuration:    {args.config} (archived to {config_copy_path})")
    print(f"  Execution Mode:   {config.get('pipeline', {}).get('execution_mode', 'local')}")
    print("=" * 70)

    # Define module output subdirectories
    dir_01 = os.path.join(run_path, "01_analysis").replace('\\', '/')
    dir_02 = os.path.join(run_path, "02_rupture_design").replace('\\', '/')
    dir_03 = os.path.join(run_path, "03_fail_fast").replace('\\', '/')
    dir_04 = os.path.join(run_path, "04_rescue_design").replace('\\', '/')
    dir_05 = os.path.join(run_path, "05_final_eval").replace('\\', '/')

    # Parse steps
    steps_to_run = set()
    if args.steps.lower() in ['all', '1-5']:
        steps_to_run = {1, 2, 3, 4, 5}
    else:
        for part in args.steps.split(','):
            if '-' in part:
                start, end = part.split('-')
                steps_to_run.update(range(int(start), int(end) + 1))
            else:
                steps_to_run.add(int(part))

    python_exe = sys.executable

    # Module 1: Interface Analysis
    if 1 in steps_to_run:
        print("\n" + "#" * 60)
        print(f"  [STEP 1/5] MODULE 1: INTERFACE ANALYSIS")
        print("#" * 60)
        cmd_1 = [python_exe, "src/01_analyze_interface.py", "--config", config_copy_path, "--out_dir", dir_01]
        subprocess.run(cmd_1, check=True)

    # Module 2: A' Rupture Generation
    if 2 in steps_to_run:
        print("\n" + "#" * 60)
        print(f"  [STEP 2/5] MODULE 2: A' RUPTURE DESIGN (RFD3 + LigandMPNN)")
        print("#" * 60)
        cmd_2 = [python_exe, "src/02_generate_A_prime.py", "--config", config_copy_path, "--analysis_dir", dir_01, "--out_dir", dir_02]
        subprocess.run(cmd_2, check=True)

    # Module 3: Fail-Fast Filter
    if 3 in steps_to_run:
        print("\n" + "#" * 60)
        print(f"  [STEP 3/5] MODULE 3: FAIL-FAST SCREENING (AF2 / AF3)")
        print("#" * 60)
        cmd_3 = [python_exe, "src/03_filter_A_prime.py", "--config", config_copy_path, "--design_dir", dir_02, "--analysis_dir", dir_01, "--out_dir", dir_03]
        subprocess.run(cmd_3, check=True)

    # Module 4: B' Rescue Generation
    if 4 in steps_to_run:
        print("\n" + "#" * 60)
        print(f"  [STEP 4/5] MODULE 4: B' RESCUE DESIGN (RFD3 + LigandMPNN)")
        print("#" * 60)
        passed_txt = os.path.join(dir_03, "passed_candidates.txt")
        cmd_4 = [python_exe, "src/04_generate_B_prime.py", "--config", config_copy_path, "--analysis_dir", dir_01, "--passed_candidates", passed_txt, "--out_dir", dir_04]
        subprocess.run(cmd_4, check=True)

    # Module 5: Final Evaluation & Orthogonality Scoring
    if 5 in steps_to_run:
        print("\n" + "#" * 60)
        print(f"  [STEP 5/5] MODULE 5: FINAL EVALUATION & ORTHOGONALITY MATRIX")
        print("#" * 60)
        cmd_5 = [python_exe, "src/05_eval_final.py", "--config", config_copy_path, "--analysis_dir", dir_01, "--design_dir", dir_04, "--filter_dir", dir_03, "--out_dir", dir_05]
        subprocess.run(cmd_5, check=True)

    # Generate Run Summary Markdown
    csv_scores = os.path.join(dir_05, "orthogonality_scores.csv")
    summary_md_path = os.path.join(run_path, "SUMMARY.md")
    
    if os.path.exists(csv_scores):
        df = pd.read_csv(csv_scores)
        with open(summary_md_path, 'w') as f:
            f.write(f"# OrthoInterface Experiment Summary: `{run_name}`\n\n")
            f.write(f"- **Execution Date**: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"- **Execution Mode**: `{config.get('pipeline', {}).get('execution_mode', 'local')}`\n")
            f.write(f"- **Folding Engine**: `{config.get('folding', {}).get('engine', 'af2')}`\n")
            f.write(f"- **RFD3 Batches**: `{config.get('pipeline', {}).get('foundry_n_batches', 5)}`\n\n")
            f.write("## 📊 Orthogonality Scores Table\n\n")
            f.write(df.to_markdown(index=False) + "\n\n")
            f.write("## 🧬 Artifacts\n\n")
            f.write(f"- CSV Scores: [`orthogonality_scores.csv`](file:///{os.path.abspath(csv_scores)})\n")
            f.write(f"- SQLite DB: [`results.db`](file:///{os.path.abspath(os.path.join(dir_05, 'results.db'))})\n")
            f.write(f"- Configuration Archive: [`config_used.yaml`](file:///{os.path.abspath(config_copy_path)})\n")

    print("\n" + "=" * 70)
    print(f"  🎉 RUN COMPLETED SUCCESSFULLY: {run_name}")
    print(f"  All results saved in: {run_path}")
    if os.path.exists(summary_md_path):
        print(f"  Summary generated:    {summary_md_path}")
    print("=" * 70)

if __name__ == "__main__":
    main()
