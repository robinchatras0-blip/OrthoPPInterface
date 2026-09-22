import os
import re
import sys
import json
import glob
import hashlib
import shutil
import subprocess
from Bio.PDB import PDBParser

sys.path.append(os.path.dirname(__file__))
from design_utils import THREE_TO_ONE, cif_to_pdb, diff_positions  # noqa: E402


def _use_wsl(config):
    fcfg = (config or {}).get('folding', {})
    return bool(fcfg.get('use_wsl', os.name == 'nt'))


def to_wsl_path(path):
    """Converts a Windows absolute or relative path to a WSL compatible path."""
    abs_path = os.path.abspath(path).replace('\\', '/')
    if len(abs_path) > 1 and abs_path[1] == ':':
        drive = abs_path[0].lower()
        return f"/mnt/{drive}{abs_path[2:]}"
    return abs_path


def get_chain_sequence(pdb_path, chain_id='A'):
    """Extracts amino acid sequence (1-letter code) of a chain from a PDB file."""
    if not os.path.exists(pdb_path):
        return ""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_path)
    model = structure[0]
    if chain_id not in model:
        chains = list(model.get_chains())
        if not chains:
            return ""
        chain = chains[0]
    else:
        chain = model[chain_id]
    return "".join(THREE_TO_ONE.get(res.get_resname(), 'X') for res in chain if res.id[0] == ' ')


# --------------------------------------------------------------------------- #
# MSA handling
# --------------------------------------------------------------------------- #
def read_a3m(path):
    """Returns [(header, sequence)] (headers without '>'); tolerant to wrapped sequences,
    blank lines and ColabFold '#' comment lines."""
    records, header, chunks = [], None, []
    with open(path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header, chunks = line[1:], []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def build_hybrid_msa(wt_a3m_path, target_sequence, modified_indices, output_a3m_path,
                     mode='gap', strict=True):
    """Builds the MSA used to fold a *designed* sequence from the WT alignment.

    modified_indices: 0-based alignment columns to treat as re-designed. Pass None to derive them
        automatically as the positions where `target_sequence` differs from the a3m query (WT).
    mode:
        'gap'        - homologs get '-' at re-designed columns (no evolutionary info there),
        'substitute' - homologs get the designed residue at those columns (column looks conserved),
        'keep'       - homologs are untouched, only the query changes (WT-biased, over-optimistic).
    strict: raise if the alignment width does not match the target sequence (a silent width
        mismatch produces a corrupt MSA); otherwise fall back to a single-sequence MSA.
    Returns a small stats dict.
    """
    if mode not in ('gap', 'substitute', 'keep'):
        raise ValueError(f"Unknown MSA mask mode: {mode}")
    stats = {"n_rows": 1, "n_masked": 0, "frac_masked": 0.0, "mode": mode, "fallback": False}

    def single_sequence():
        with open(output_a3m_path, 'w') as f:
            f.write(f">query\n{target_sequence}\n")
        stats["fallback"] = True
        return stats

    if not wt_a3m_path or not os.path.exists(wt_a3m_path):
        return single_sequence()
    records = read_a3m(wt_a3m_path)
    if not records:
        return single_sequence()

    q_header, q_seq = records[0]
    q_wt = "".join(c for c in q_seq if not c.islower() and c != '-')
    if len(q_wt) != len(target_sequence):
        msg = (f"MSA width {len(q_wt)} != sequence length {len(target_sequence)} "
               f"({os.path.basename(wt_a3m_path)}). Chain cropped/renumbered or wrong a3m?")
        if strict:
            raise ValueError(msg)
        print(f"  [MSA WARNING] {msg} -> falling back to single-sequence input.")
        return single_sequence()

    if modified_indices is None:
        modified_indices = diff_positions(q_wt, target_sequence)
    modified = {i for i in modified_indices if 0 <= i < len(target_sequence)}

    kept = [i for i in range(len(target_sequence)) if i not in modified]
    if kept:
        identity = sum(q_wt[i] == target_sequence[i] for i in kept) / len(kept)
        if identity < 0.9:
            print(f"  [MSA WARNING] Only {identity:.0%} identity between a3m query and design outside "
                  f"masked columns -> residue numbering / column mapping is probably wrong.")

    out = [f">{q_header}\n{target_sequence}\n"]
    for header, seq in records[1:]:
        chars, col = [], 0
        for ch in seq:
            if ch.islower():                      # insertion relative to the query
                if mode != 'keep' and (col - 1) in modified:
                    continue
                chars.append(ch)
                continue
            if col in modified and mode != 'keep':
                if mode == 'gap' or ch == '-':
                    chars.append('-')
                else:
                    chars.append(target_sequence[col])
            else:
                chars.append(ch)
            col += 1
        if col != len(target_sequence):
            continue                              # malformed homolog row, drop it
        out.append(f">{header}\n{''.join(chars)}\n")

    with open(output_a3m_path, 'w') as f:
        f.writelines(out)
    stats.update(n_rows=len(out), n_masked=len(modified),
                 frac_masked=len(modified) / max(1, len(target_sequence)))
    return stats


def prepare_msa(wt_a3m, target_sequence, out_path, config, modified_indices=None,
                reference_sequence=None, interface_columns=None):
    """Regime-aware MSA preparation driven by config['folding'].

    msa_mask_scope: 'interface' -> mask the interface columns (`interface_columns`) plus every position that differs
                                   from WT: the alignment cannot vouch for the interface, whichever sequence is folded
                                   (default; WT and designs are compared in the same information regime),
                    'diff'      -> mask only positions that really differ from WT,
                    'mutable'   -> mask every position in `modified_indices` (legacy).
    msa_mask_mode : 'gap' | 'substitute' | 'keep'.
    `reference_sequence` lets the caller mask the columns where ANOTHER design differs (used to
    put a WT chain in the same information regime as its designed counterpart).
    """
    fcfg = config.get('folding', {})
    mode = fcfg.get('msa_mask_mode', 'gap')
    scope = fcfg.get('msa_mask_scope', 'interface')
    if scope not in ('interface', 'diff', 'mutable'):
        raise ValueError(f"Unknown folding.msa_mask_scope: {scope}")
    if not fcfg.get('use_msa', True):
        with open(out_path, 'w') as f:
            f.write(f">query\n{target_sequence}\n")
        return {"fallback": True, "mode": "none", "n_masked": 0, "frac_masked": 0.0, "n_rows": 1}
    wt_seq = _a3m_query(wt_a3m) or target_sequence
    ref = target_sequence if reference_sequence is None else reference_sequence
    cols = diff_positions(wt_seq, ref) if len(wt_seq) == len(ref) else []
    if scope == 'interface':
        if interface_columns is None:
            raise ValueError("folding.msa_mask_scope 'interface' requires the interface columns of the chain")
        cols = sorted(set(cols) | set(interface_columns))
    elif scope == 'mutable' and modified_indices is not None and reference_sequence is None:
        cols = modified_indices
    return build_hybrid_msa(wt_a3m, target_sequence, cols, out_path, mode=mode,
                            strict=fcfg.get('strict_msa', True))


def _a3m_query(path):
    if not path or not os.path.exists(path):
        return ""
    recs = read_a3m(path)
    return "".join(c for c in recs[0][1] if not c.islower() and c != '-') if recs else ""


# --------------------------------------------------------------------------- #
# RF3 wrapper
# --------------------------------------------------------------------------- #
def _sha1_file(path):
    h = hashlib.sha1()
    with open(path, 'rb') as f:
        h.update(f.read())
    return h.hexdigest()


def _signature(components_raw, rf3_params):
    """Cache key: sequences + MSA *contents* + inference parameters."""
    payload = []
    for c in components_raw:
        msa = c.get("msa_path")
        payload.append({"seq": c["seq"], "chain_id": c["chain_id"],
                        "msa": _sha1_file(msa) if msa and os.path.exists(msa) else None})
    blob = json.dumps({"c": payload, "p": rf3_params}, sort_keys=True)
    return hashlib.sha1(blob.encode()).hexdigest()


def _f(data, key, default=0.0):
    v = data.get(key, default)
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def _parse_summary(path):
    with open(path, 'r') as f:
        data = json.load(f)
    plddt = _f(data, "overall_plddt", _f(data, "plddt", 0.0))
    plddt = plddt * 100.0 if plddt <= 1.0 else plddt
    iptm = _f(data, "iptm")
    pae_min = data.get("chain_pair_pae_min")
    try:                                     # AF3-style nested matrix -> inter-chain entry
        pae_min = float(pae_min[0][1]) if pae_min else None
    except (TypeError, IndexError, ValueError):
        pae_min = None
    return {
        "plddt": plddt, "iptm": iptm, "ptm": _f(data, "ptm"),
        "ranking_score": _f(data, "ranking_score", iptm),
        "has_clash": bool(data.get("has_clash", False)),
        "chain_pair_pae_min": pae_min,
    }


def _collect_metrics(out_dir):
    """Aggregates every diffusion sample RF3 wrote. Reports the best-ranked sample (headline
    metrics) plus mean/std of iPTM over samples so that noise is visible."""
    files = sorted(glob.glob(os.path.join(out_dir, "**", "*_summary_confidences.json"), recursive=True))
    # RF3 also writes a top-level summary that duplicates the best sample: count samples only once
    files = [p for p in files if re.search(r"seed-\d+_sample-\d+", os.path.basename(p))] or files
    if not files:
        return None, None
    parsed = [(p, _parse_summary(p)) for p in files]
    best_path, best = max(parsed, key=lambda t: t[1]["ranking_score"])
    iptms = [m["iptm"] for _, m in parsed]
    mean = sum(iptms) / len(iptms)
    std = (sum((x - mean) ** 2 for x in iptms) / len(iptms)) ** 0.5
    best = dict(best, iptm_mean=mean, iptm_std=std, n_samples=len(parsed))
    return best_path, best


def _find_model(out_dir, tag, summary_path):
    prefix = summary_path[:-len("_summary_confidences.json")]
    for ext in ("_model.cif", "_model.cif.gz"):
        if os.path.exists(prefix + ext):
            return prefix + ext
    for pat in (f"{tag}_model.cif*", "*_model.cif*", "*.cif", "*.cif.gz"):
        hits = sorted(glob.glob(os.path.join(out_dir, "**", pat), recursive=True))
        if hits:
            return hits[0]
    return None


def _rf3_settings(config):
    fcfg = config.get('folding', {})
    wsl = _use_wsl(config)
    pth = to_wsl_path if wsl else (lambda p: os.path.abspath(p).replace('\\', '/'))
    params = {k: fcfg[k] for k in ('diffusion_batch_size', 'n_recycles', 'num_steps', 'seed') if k in fcfg}
    return fcfg, wsl, pth, params


def _prepare_request(fasta_sequences, out_dir, config):
    """Manifest, signature and cache state of one prediction (nothing is run here)."""
    os.makedirs(out_dir, exist_ok=True)
    tag = os.path.basename(out_dir.rstrip('/\\'))
    fcfg, wsl, pth, params = _rf3_settings(config)
    use_msa = fcfg.get('use_msa', True)
    raw = []
    for idx, (_, seq, msa_p) in enumerate(fasta_sequences):
        comp = {"seq": seq, "chain_id": chr(ord('A') + idx)}
        if use_msa and msa_p and os.path.exists(msa_p) and os.path.getsize(msa_p) > 0:
            comp["msa_path"] = msa_p
        raw.append(comp)
    entry = {"name": tag, "components": [{**c, **({"msa_path": pth(c["msa_path"])} if "msa_path" in c else {})} for c in raw]}
    manifest_path = os.path.join(out_dir, f"{tag}_rf3_input.json")
    with open(manifest_path, 'w') as f:
        json.dump([entry], f, indent=2)
    sig_path = os.path.join(out_dir, f"{tag}.signature")
    signature = _signature(raw, params)
    cached_sig = open(sig_path).read().strip() if os.path.exists(sig_path) else None
    has_output = _collect_metrics(out_dir)[0] is not None
    return {"out_dir": out_dir, "tag": tag, "entry": entry, "signature": signature, "sig_path": sig_path,
            "manifest_path": manifest_path, "cached": has_output and cached_sig == signature,
            "stale": has_output and cached_sig != signature}


def _run_rf3(pending, config):
    """ONE RF3 process for every prediction that is not cached: the model is loaded once, not once per fold."""
    fcfg, wsl, pth, params = _rf3_settings(config)
    for p in pending:
        if p["stale"]:
            print(f"  [Foundry RF3] Stale cache in {p['out_dir']} (inputs changed) -> recomputing.")
            for f in glob.glob(os.path.join(p["out_dir"], "**", "*"), recursive=True):
                if os.path.isfile(f) and not f.endswith(("_rf3_input.json", ".signature", ".a3m")):
                    os.remove(f)
    batch_dir = None
    if len(pending) == 1:
        manifest, out = pending[0]["manifest_path"], pending[0]["out_dir"]
    else:
        batch_dir = os.path.join(os.path.dirname(os.path.abspath(pending[0]["out_dir"])), f"_rf3_batch_{os.getpid()}")
        os.makedirs(batch_dir, exist_ok=True)
        entries = []
        for k, p in enumerate(pending):
            p["uname"] = f"{p['tag']}__{k:03d}"                     # tags repeat across folds (rescue / negative of a design)
            entries.append({**p["entry"], "name": p["uname"]})
        manifest, out = os.path.join(batch_dir, "batch_rf3_input.json"), batch_dir
        with open(manifest, 'w') as f:
            json.dump(entries, f, indent=2)
    rf3_bin = fcfg.get('rf3_bin', '/home/ommearo/miniforge3/envs/foundry_env/bin/rf3')
    rf3_ckpt = fcfg.get('rf3_ckpt', '/home/ommearo/.foundry/checkpoints/rf3_foundry_01_24_latest.ckpt')
    cmd = [rf3_bin, "fold", f"inputs={pth(manifest)}", f"out_dir={pth(out)}", f"ckpt_path={rf3_ckpt}"]
    cmd += [f"{k}={v}" for k, v in params.items()]
    if wsl:
        cmd = ["wsl", "-d", fcfg.get('wsl_distro', 'Ubuntu'), "--"] + cmd
    print(f"  [Foundry RF3 GPU Invoc] {len(pending)} prediction(s) in one process: {' '.join(cmd[:6])} ...", flush=True)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        tail = (res.stderr or res.stdout or "")[-600:]
        raise RuntimeError(f"RF3 execution failed with code {res.returncode}: {tail}")
    if batch_dir:
        for p in pending:
            src, dst = os.path.join(batch_dir, p["uname"]), os.path.join(p["out_dir"], p["uname"])
            if not os.path.isdir(src):
                raise RuntimeError(f"RF3 wrote no output for {p['tag']} (batch of {len(pending)}, see {batch_dir})")
            if os.path.exists(dst):
                shutil.rmtree(dst)
            shutil.move(src, dst)
        shutil.rmtree(batch_dir, ignore_errors=True)
    for p in pending:
        if _collect_metrics(p["out_dir"])[0] is None:
            raise RuntimeError(f"RF3 output *_summary_confidences.json not found in {p['out_dir']}")
        with open(p["sig_path"], 'w') as f:
            f.write(p["signature"])


def _finalize(p):
    summary_path, metrics = _collect_metrics(p["out_dir"])
    out_pdb = os.path.join(p["out_dir"], f"{p['tag']}_unrelaxed_rank_001.pdb")
    if not os.path.exists(out_pdb):
        cif = _find_model(p["out_dir"], p["tag"], summary_path)
        if cif:
            try:
                cif_to_pdb(cif, out_pdb)
                print(f"  [Foundry RF3] Converted {os.path.basename(cif)} -> {os.path.basename(out_pdb)}")
            except Exception as e:
                print(f"  Warning: could not convert {cif} to PDB: {e}")
    metrics["model_pdb"] = out_pdb if os.path.exists(out_pdb) else None
    return metrics


def predict_structures(requests, config):
    """Runs RoseTTAFold-3 on many monomers/complexes with ONE RF3 process.

    requests: [(fasta_sequences, out_dir)] with fasta_sequences = [(label, sequence, msa_path_or_None), ...] (chain ids
    A, B, ...). Each prediction is cached in its own out_dir together with a signature of sequences + MSA contents +
    parameters, so a change of the MSA strategy invalidates stale results and only missing predictions are run.
    Returns one metrics dict per request (best-ranked sample, iPTM mean/std over samples, `model_pdb`).
    """
    preps = [_prepare_request(fasta, out_dir, config) for fasta, out_dir in requests]
    for p in preps:
        if p["cached"]:
            print(f"  [Foundry RF3] Cache hit ({p['tag']}).")
    pending = [p for p in preps if not p["cached"]]
    if pending:
        _run_rf3(pending, config)
    return [_finalize(p) for p in preps]


def predict_structure(fasta_sequences, out_dir, config):
    """Single prediction (see predict_structures)."""
    return predict_structures([(fasta_sequences, out_dir)], config)[0]


def calculate_ca_rmsd(ref_pdb, pred_pdb_or_dir, chain_ref='A', chain_pred=None, subset_res_ids=None):
    """Computes C-alpha RMSD between a reference and a predicted structure (0.0 if not computable).

    Residues are paired by residue id when the numbering overlaps, otherwise by order. Optionally
    restricted to `subset_res_ids` (reference numbering).
    """
    from design_utils import match_residues, ca_rmsd_after_fit

    pred_pdb = pred_pdb_or_dir
    if os.path.isdir(pred_pdb_or_dir):
        ranked = (sorted(glob.glob(os.path.join(pred_pdb_or_dir, "*unrelaxed_rank_001*.pdb")))
                  + sorted(glob.glob(os.path.join(pred_pdb_or_dir, "*rank_001*.pdb")))
                  + sorted(glob.glob(os.path.join(pred_pdb_or_dir, "*.pdb"))))
        if not ranked:
            return 0.0
        pred_pdb = ranked[0]
    if not os.path.exists(ref_pdb) or not os.path.exists(pred_pdb):
        return 0.0

    parser = PDBParser(QUIET=True)
    model_ref = parser.get_structure("ref", ref_pdb)[0]
    model_pred = parser.get_structure("pred", pred_pdb)[0]
    ch_ref = model_ref[chain_ref] if chain_ref in model_ref else list(model_ref.get_chains())[0]
    ch_pred = model_pred[chain_pred] if chain_pred and chain_pred in model_pred else list(model_pred.get_chains())[0]

    pairs = match_residues(ch_ref, ch_pred)
    subset = set(subset_res_ids) if subset_res_ids is not None else None
    rmsd, _, _ = ca_rmsd_after_fit(pairs, fit_ids=subset)
    return float(rmsd) if rmsd is not None else 0.0
