"""Calibrates the RF3 MSA regime on in-silico controls BEFORE spending GPU time on designs.

Positive controls : WT chain A with a few CONSERVATIVE substitutions at interface positions (should still bind B_WT).
Negative controls : WT chain A with disruptive substitutions at interface positions (should not bind B_WT).
For every MSA regime (mask mode x scope) we fold all controls and report how well iPTM separates the
two groups. Pick the regime with the largest gap; if no regime separates them, iPTM is not a usable
filter for this system and structural / energy metrics must carry the decision.

Usage: python src/calibrate_predictor.py --config config.yaml --analysis_dir runs/<run>/01_analysis
"""
import argparse
import copy
import json
import os
import random
import sys

import pandas as pd

sys.path.append(os.path.dirname(__file__))
from design_utils import load_config  # noqa: E402
from folding_engine import get_chain_sequence, prepare_msa, predict_structure  # noqa: E402

CONSERVATIVE = {'I': 'V', 'V': 'I', 'L': 'I', 'M': 'L', 'K': 'R', 'R': 'K', 'D': 'E', 'E': 'D',
                'S': 'T', 'T': 'S', 'N': 'Q', 'Q': 'N', 'F': 'Y', 'Y': 'F', 'A': 'S', 'H': 'N', 'W': 'F', 'C': 'S', 'G': 'A', 'P': 'A'}
DISRUPTIVE = {'K': 'E', 'R': 'D', 'D': 'K', 'E': 'K', 'H': 'D', 'N': 'W', 'Q': 'W', 'S': 'W', 'T': 'W',
              'Y': 'D', 'F': 'D', 'W': 'D', 'L': 'D', 'I': 'K', 'V': 'K', 'M': 'E', 'A': 'W', 'C': 'W', 'G': 'W', 'P': 'W'}


def mutate(seq, positions, table):
    out = list(seq)
    for p in positions:
        out[p] = table.get(out[p], 'A')
    return "".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default='config.yaml')
    ap.add_argument('--analysis_dir', default='results/01_analysis')
    ap.add_argument('--out_dir', default='results/00_calibration')
    ap.add_argument('--n_controls', type=int, default=6)
    ap.add_argument('--n_pos_mut', type=int, default=5)
    ap.add_argument('--n_neg_mut', type=int, default=8)
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()

    config = load_config(args.config)
    p = config['pipeline']
    wt_pdb, chA, chB = p['input_pdb'], p.get('chain_A', 'A'), p.get('chain_B', 'B')
    a3m_A, a3m_B = p['input_msa_A'], p['input_msa_B']
    seq_A, seq_B = get_chain_sequence(wt_pdb, chA), get_chain_sequence(wt_pdb, chB)
    mapping = json.load(open(os.path.join(args.analysis_dir, 'index_mapping.json')))
    iface = [r['sequential_index'] - 1 for r in mapping['residues_A'] if r['is_interface']]
    neigh = {r['sequential_index'] - 1 for r in mapping['residues_A'] if r['is_neighborhood']}
    if len(iface) < args.n_neg_mut:
        sys.exit(f"Only {len(iface)} interface residues on A; lower --n_neg_mut.")

    rng = random.Random(args.seed)
    controls = []
    for i in range(args.n_controls):
        controls.append((f"pos_{i}", 1, mutate(seq_A, rng.sample(iface, args.n_pos_mut), CONSERVATIVE)))
        controls.append((f"neg_{i}", 0, mutate(seq_A, rng.sample(iface, args.n_neg_mut), DISRUPTIVE)))

    regimes = [("gap", "diff"), ("substitute", "diff"), ("keep", "diff"), ("gap", "mutable")]
    rows = []
    for mode, scope in regimes:
        cfg = copy.deepcopy(config)
        cfg['folding'].update(msa_mask_mode=mode, msa_mask_scope=scope)
        for name, label, seq in controls:
            d = os.path.join(args.out_dir, f"{mode}_{scope}", name)
            os.makedirs(d, exist_ok=True)
            msa_a = os.path.join(d, "A.a3m")
            prepare_msa(a3m_A, seq, msa_a, cfg, modified_indices=sorted(neigh))
            m = predict_structure([('A_var', seq, msa_a), ('B_wt', seq_B, a3m_B)], d, cfg)
            rows.append({"mode": mode, "scope": scope, "control": name, "expected_binder": label,
                         "iptm": m["iptm"], "ptm": m["ptm"], "plddt": m["plddt"]})
            print(f"  {mode:10s}/{scope:7s} {name:6s} iPTM={m['iptm']:.3f}")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.out_dir, "calibration_raw.csv"), index=False)
    summ = (df.groupby(["mode", "scope", "expected_binder"]).iptm.mean().unstack()
              .rename(columns={0: "iptm_neg_mean", 1: "iptm_pos_mean"}))
    summ["gap"] = summ["iptm_pos_mean"] - summ["iptm_neg_mean"]
    summ = summ.sort_values("gap", ascending=False)
    summ.to_csv(os.path.join(args.out_dir, "calibration_summary.csv"))
    print("\n" + summ.to_string())
    best = summ.index[0]
    print(f"\nBest regime: msa_mask_mode={best[0]}, msa_mask_scope={best[1]} (gap {summ.iloc[0]['gap']:.3f})")
    if summ.iloc[0]["gap"] < 0.15:
        print("WARNING: no regime separates conservative from disruptive variants (gap < 0.15). "
              "iPTM cannot be trusted as the decisive filter for this system.")


if __name__ == "__main__":
    main()
