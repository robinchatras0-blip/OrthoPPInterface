#!/usr/bin/env python3
"""Master Overnight Pipeline Runner for OrthoPPInterface.
Executes Module 4 (Rescue B' Design with GPU batching) then Module 5 (High-throughput RF3 evaluation + PyRosetta).
"""
import os
import sys
import shutil
import subprocess
import time

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN_DIR = os.path.join(BASE_DIR, "runs", "run_overnight_production")
CONFIG_PATH = os.path.join(RUN_DIR, "config_used.yaml")
REF_ANALYSIS_DIR = os.path.join(BASE_DIR, "runs", "run_negative_dockq_optimized", "01_analysis")
REF_FILTER_DIR = os.path.join(BASE_DIR, "runs", "run_negative_dockq_optimized", "03_fail_fast")
PASSED_CANDIDATES = os.path.join(REF_FILTER_DIR, "passed_candidates.txt")

RUN_ANALYSIS_DIR = os.path.join(RUN_DIR, "01_analysis")
RUN_DESIGN_DIR = os.path.join(RUN_DIR, "04_rescue_design")
RUN_EVAL_DIR = os.path.join(RUN_DIR, "05_final_eval")


def run_step(cmd, desc):
    print(f"\n{'='*70}\n[OVERNIGHT PIPELINE] STARTING: {desc}\n{'='*70}", flush=True)
    t0 = time.time()
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    env["TORCH_CUDA_MATMUL_PRECISION"] = "high"
    res = subprocess.run(cmd, env=env)
    elapsed = time.time() - t0
    if res.returncode != 0:
        print(f"\n[OVERNIGHT PIPELINE] ERROR: {desc} failed with code {res.returncode} after {elapsed:.1f}s", flush=True)
        sys.exit(res.returncode)
    print(f"\n[OVERNIGHT PIPELINE] COMPLETED: {desc} in {elapsed/60:.2f} minutes\n", flush=True)


def main():
    os.makedirs(RUN_DIR, exist_ok=True)
    os.makedirs(RUN_DESIGN_DIR, exist_ok=True)
    os.makedirs(RUN_EVAL_DIR, exist_ok=True)

    # 1. Setup analysis dir
    if not os.path.exists(RUN_ANALYSIS_DIR) and os.path.exists(REF_ANALYSIS_DIR):
        print(f"[Setup] Mirroring analysis directory to {RUN_ANALYSIS_DIR}...", flush=True)
        shutil.copytree(REF_ANALYSIS_DIR, RUN_ANALYSIS_DIR)

    # 2. Module 4: Rescue B' Design
    py_bin = sys.executable
    cmd_m4 = [
        py_bin, os.path.join(BASE_DIR, "src", "04_generate_B_prime.py"),
        "--config", CONFIG_PATH,
        "--analysis_dir", RUN_ANALYSIS_DIR,
        "--passed_candidates", PASSED_CANDIDATES,
        "--out_dir", RUN_DESIGN_DIR
    ]
    run_step(cmd_m4, "Module 4: Co-adaptation & Diverse Rescue B' Generation (RFD3 + LigandMPNN)")

    # 3. Module 5: High-throughput RF3 Evaluation & PyRosetta
    cmd_m5 = [
        py_bin, os.path.join(BASE_DIR, "src", "05_eval_final.py"),
        "--config", CONFIG_PATH,
        "--analysis_dir", RUN_ANALYSIS_DIR,
        "--design_dir", RUN_DESIGN_DIR,
        "--filter_dir", REF_FILTER_DIR,
        "--out_dir", RUN_EVAL_DIR
    ]
    run_step(cmd_m5, "Module 5: High-Throughput RF3 Orthogonality Evaluation & PyRosetta")

    print(f"\n{'='*70}\n[OVERNIGHT PIPELINE] ALL MODULES COMPLETED SUCCESSFULLY!\n{'='*70}\n", flush=True)


if __name__ == "__main__":
    main()
