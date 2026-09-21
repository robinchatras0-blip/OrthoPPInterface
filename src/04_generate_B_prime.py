import argparse
import yaml
import json
import os
import shutil
import glob
import subprocess
import gzip
import sys
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Superimposer, Structure, Model, Chain

sys.path.append(os.path.dirname(__file__))
from design_utils import (  # noqa: E402
    THREE_TO_ONE, select_diverse, std_residues, match_residues, scaffold_metrics,
    ca_rmsd_after_fit, build_rescue_bias, candidate_specificity_score, parse_mpnn_confidence,
    load_chain,
)


def save_structure(obj, path):
    io = PDBIO()
    io.set_structure(obj)
    io.save(path)


def strip_nan_lines(path):
    with open(path, 'r') as f:
        lines = [l for l in f if 'nan' not in l]
    with open(path, 'w') as f:
        f.writelines(lines)


def cif_to_pdb(cif_path):
    pdb_path = cif_path.replace(".cif.gz", ".pdb").replace(".cif", ".pdb")
    parser = MMCIFParser(QUIET=True)
    if cif_path.endswith('.gz'):
        with gzip.open(cif_path, 'rt') as f:
            struct = parser.get_structure("cif", f)
    else:
        struct = parser.get_structure("cif", cif_path)
    save_structure(struct, pdb_path)
    return pdb_path


def collect_pdbs(directory):
    """All PDB outputs in a directory tree, converting CIF(.gz) when no PDB was written."""
    pdbs = sorted(glob.glob(os.path.join(directory, "**", "*.pdb"), recursive=True))
    if pdbs:
        return pdbs
    cifs = sorted(glob.glob(os.path.join(directory, "**", "*.cif.gz"), recursive=True)
                  + glob.glob(os.path.join(directory, "**", "*.cif"), recursive=True))
    out = []
    for c in cifs:
        try:
            out.append(cif_to_pdb(c))
        except Exception as e:
            print(f"Warning: could not convert {c} to PDB: {e}")
    return out


def renumber_like(chain, ref_ids):
    """RFD3 may renumber residues. If the chain has as many residues as the reference id list but
    different numbering, re-apply the reference ids so that fixed-residue bookkeeping stays valid."""
    res = std_residues(chain)
    if len(res) != len(ref_ids) or {r.id[1] for r in res} == set(ref_ids):
        return chain.copy()
    new = Chain.Chain(chain.id)
    for r, rid in zip(res, ref_ids):
        nr = r.copy()
        nr.id = (' ', rid, ' ')
        new.add(nr)
    return new


def assemble_complex(scaffold_pdb, a_chain, wt_b_chain, crop_ids, fixed_b_ids, chain_A_id, chain_B_id, out_pdb):
    """A' (all-atom) + full-length B built from a diffused/native cropped scaffold.

    The WT head (< crop) and tail (> crop) of B are re-attached after fitting a *copy* of WT B onto
    the scaffold framework (the WT structure itself is never modified). Returns the B chain or None.
    """
    st = PDBParser(QUIET=True).get_structure("sc", scaffold_pdb)
    model = st[0]
    if chain_B_id not in model:
        return None
    sc_b = renumber_like(model[chain_B_id], crop_ids)
    wt_b = wt_b_chain.copy()
    pairs = [(r, m) for r, m in match_residues(sc_b, wt_b)
             if r.id[1] in fixed_b_ids and 'CA' in r and 'CA' in m]
    if len(pairs) < 3:
        print(f"    Skipping {os.path.basename(scaffold_pdb)}: only {len(pairs)} framework residues matched WT B.")
        return None
    sup = Superimposer()
    sup.set_atoms([r['CA'] for r, _ in pairs], [m['CA'] for _, m in pairs])
    sup.apply(list(wt_b.get_atoms()))

    lo, hi = min(crop_ids), max(crop_ids)
    new_b = Chain.Chain(chain_B_id)
    for r in std_residues(wt_b):
        if r.id[1] < lo:
            new_b.add(r.copy())
    for r in std_residues(sc_b):
        new_b.add(r.copy())
    for r in std_residues(wt_b):
        if r.id[1] > hi:
            new_b.add(r.copy())

    cx = Structure.Structure("assembled")
    cm = Model.Model(0)
    cx.add(cm)
    cm.add(a_chain.copy())
    cm.add(new_b)
    save_structure(cx, out_pdb)
    return new_b


def pick_scaffolds(rfd_pdbs, a_chain, wt_b_chain, crop_ids, fixed_b_ids, chain_B_id, cfg, native_pdb):
    """Filters, ranks and de-duplicates RFD3 scaffolds. Returns [(pdb, is_native, metrics)]."""
    max_sc = int(cfg.get('max_rescue_scaffolds', 4))
    include_native = bool(cfg.get('include_native_scaffold', False))
    max_clash = int(cfg.get('scaffold_max_clashes', 5))
    max_fw = float(cfg.get('scaffold_max_framework_rmsd', 1.5))
    min_div = float(cfg.get('scaffold_min_flex_rmsd', 0.5))

    scored = []
    for p in rfd_pdbs:
        try:
            ch = PDBParser(QUIET=True).get_structure("s", p)[0][chain_B_id]
        except Exception:
            continue
        ch = renumber_like(ch, crop_ids)
        m = scaffold_metrics(a_chain, ch, wt_b_chain, fixed_b_ids)
        m["_chain"] = ch
        if m["clashes"] > max_clash:
            continue
        if m["framework_rmsd"] is not None and m["framework_rmsd"] > max_fw:
            continue
        scored.append((m["contacts"] - 5 * m["clashes"], p, m))
    scored.sort(key=lambda t: -t[0])

    picked = []
    n_rfd_target = max_sc - (1 if include_native else 0)
    for _, p, m in scored:
        if len(picked) >= n_rfd_target:
            break
        # keep backbones that are actually different from those already selected
        is_dup = False
        for _, _, m2 in picked:
            pairs = match_residues(m2["_chain"], m["_chain"])
            flex = {r.id[1] for r, _ in pairs} - set(fixed_b_ids)
            _, ev, _ = ca_rmsd_after_fit(pairs, fit_ids=set(fixed_b_ids), eval_ids=flex or None)
            if ev is not None and ev < min_div:
                is_dup = True
                break
        if not is_dup:
            picked.append((p, False, m))

    result = [(p, False, {k: v for k, v in m.items() if k != "_chain"}) for p, _, m in picked]
    if include_native or not result:
        if not result:
            print("  WARNING: no RFD3 scaffold passed the filters -> falling back to the native WT backbone.")
        result.insert(0, (native_pdb, True, {"contacts": None, "clashes": None, "framework_rmsd": 0.0, "flex_rmsd": 0.0}))
    return result


def zscore(values):
    vals = [v for v in values if v is not None]
    if len(vals) < 2:
        return [0.0] * len(values)
    mean = sum(vals) / len(vals)
    sd = (sum((v - mean) ** 2 for v in vals) / len(vals)) ** 0.5 or 1.0
    return [0.0 if v is None else (v - mean) / sd for v in values]


def main():
    parser = argparse.ArgumentParser(description="Module 4: Diverse Rescue B' Generation")
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--analysis_dir', default='results/01_analysis')
    parser.add_argument('--passed_candidates', default='results/03_fail_fast/passed_candidates.txt')
    parser.add_argument('--out_dir', default='results/04_rescue_design')
    parser.add_argument('--force', action='store_true', help="Ignore cached RFD3/B' results")
    args = parser.parse_args()

    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    os.makedirs(args.out_dir, exist_ok=True)

    with open(args.passed_candidates, 'r') as f:
        raw_lines = [l.strip() for l in f if l.strip()]
    if not raw_lines:
        raise ValueError("No passed candidates found!")

    cand_dir = os.path.dirname(os.path.abspath(args.passed_candidates))
    parent_dir = os.path.dirname(cand_dir)
    resolved_candidates = []
    for p in raw_lines:
        for c in (os.path.join(parent_dir, "02_rupture_design", os.path.basename(p)), p,
                  os.path.join(cand_dir, p), os.path.join("results/02_rupture_design", os.path.basename(p))):
            if os.path.exists(c):
                resolved_candidates.append(c)
                break

    pcfg = config['pipeline']
    execution_mode = pcfg.get('execution_mode', 'mock')
    chain_A_id, chain_B_id = pcfg.get('chain_A', 'A'), pcfg.get('chain_B', 'B')
    wt_pdb = pcfg.get('input_pdb', 'data/inputs/complex_S1_S2.pdb')

    with open(os.path.join(args.analysis_dir, 'index_mapping.json')) as f:
        mapping = json.load(f)
    crop_min_b, crop_max_b = mapping.get('crop_min_B', 1), mapping.get('crop_max_B', 88)
    interface_indices_A = mapping.get('neighborhood_ids', [])
    fixed_ids_A = set(mapping.get('fixed_ids_A', []))
    with open(os.path.join(args.analysis_dir, 'mpnn_fixed_positions_B.json')) as f:
        mpnn_fixed = json.load(f)
    fixed_b_ids = set(mpnn_fixed.get(chain_B_id, []))

    # Diversity selection over A' motifs (interface sequence, correct Hamming this time)
    div_cfg = config.get('diversity_selection', {})
    max_k_motifs = div_cfg.get('max_a_prime_motifs', 3)
    if div_cfg.get('enabled', True) and len(resolved_candidates) > 1:
        seqs = []
        for p in resolved_candidates:
            ch = load_chain(p, chain_A_id)
            seqs.append({r.id[1]: THREE_TO_ONE.get(r.get_resname(), 'X') for r in std_residues(ch)
                         if r.id[1] in interface_indices_A})
        idx = select_diverse(seqs, max_k_motifs)
        selected_candidates = [resolved_candidates[i] for i in idx]
        print(f"  [Diversity Selector] Kept {len(idx)} A' motifs (indices {idx}) out of {len(resolved_candidates)}.")
    else:
        selected_candidates = resolved_candidates[:max_k_motifs]

    lmpnn_cfg = config.get('ligandmpnn', {})
    coadapt_enabled = lmpnn_cfg.get('rescue_coadaptation', True)
    samples_per_scaffold = int(lmpnn_cfg.get('samples_per_scaffold', lmpnn_cfg.get('rescue_n_batches', 4) * 5))
    top_k_per_scaffold = int(lmpnn_cfg.get('top_k_per_scaffold', 3))
    top_k_per_motif = int(lmpnn_cfg.get('top_k_per_motif', top_k_per_scaffold * 2))
    max_clash_a = int(lmpnn_cfg.get('max_clashes_with_a_prime', 3))
    temp_rescue = str(lmpnn_cfg.get('temperature_rescue', 0.20))
    model_type = lmpnn_cfg.get('model_type', "ligand_mpnn")
    is_legacy = str(lmpnn_cfg.get('is_legacy_weights', "True"))
    default_ckpt = ("OrthoIntRob/ligandmpnn/model_params/ligandmpnn_v_32_010_25.pt" if execution_mode == 'local'
                    else "/content/drive/MyDrive/OrthoPPInterface_Data/FoundryModels/LigandMPNN/model_params/ligandmpnn_v_32_010_25.pt")
    checkpoint_path = lmpnn_cfg.get('checkpoint_path', default_ckpt)

    metadata_pairs, all_saved_b_pdbs = {}, []
    print(f"\nModule 4: Co-adapting and designing rescue B' for {len(selected_candidates)} diverse A' motifs...")

    wt_struct = PDBParser(QUIET=True).get_structure("WT", wt_pdb)
    wt_A_chain = wt_struct[0][chain_A_id]
    wt_B_chain = wt_struct[0][chain_B_id]
    crop_ids = sorted(r.id[1] for r in std_residues(wt_B_chain) if crop_min_b <= r.id[1] <= crop_max_b)

    for a_idx, in_pdb in enumerate(selected_candidates):
        a_cand_name = os.path.splitext(os.path.basename(in_pdb))[0]
        print(f"\n{'=' * 60}\n  [A' Motif {a_idx + 1}/{len(selected_candidates)}]: {a_cand_name}\n{'=' * 60}")
        motif_work_dir = os.path.join(args.out_dir, f"motif_{a_cand_name}")
        os.makedirs(motif_work_dir, exist_ok=True)

        existing = sorted(glob.glob(os.path.join(args.out_dir, f"{a_cand_name}_B_cand_*.pdb")))
        meta_path = os.path.join(args.out_dir, "diversity_pairs_metadata.json")
        if existing and not args.force:
            print(f"  [Resume Cache] {len(existing)} existing B' candidates for {a_cand_name} (use --force to redo).")
            old_meta = json.load(open(meta_path)) if os.path.exists(meta_path) else {}
            for p in existing:
                pid = os.path.splitext(os.path.basename(p))[0]
                all_saved_b_pdbs.append(p)
                metadata_pairs[pid] = old_meta.get(pid, {"parent_a": a_cand_name, "pdb": p})
            continue

        # ---- A' in the WT frame (fit on the unchanged framework, not on the redesigned interface)
        struct_A = PDBParser(QUIET=True).get_structure("A_prime", in_pdb)
        a_chain = struct_A[0][chain_A_id]
        pairs = [(r, m) for r, m in match_residues(wt_A_chain, a_chain) if 'CA' in r and 'CA' in m]
        fw_pairs = [(r, m) for r, m in pairs if r.id[1] in fixed_ids_A] or pairs
        sup = Superimposer()
        sup.set_atoms([r['CA'] for r, _ in fw_pairs], [m['CA'] for _, m in fw_pairs])
        sup.apply(list(struct_A.get_atoms()))
        print(f"  Superimposed {a_cand_name} onto WT frame on {len(fw_pairs)} framework CA (RMSD = {sup.rms:.3f} A)")

        # ---- A' + cropped WT B as RFD3 input
        cx = Structure.Structure("complex_A_prime_B_wt")
        cm = Model.Model(0)
        cx.add(cm)
        cm.add(a_chain.copy())
        b_crop = Chain.Chain(chain_B_id)
        for r in std_residues(wt_B_chain):
            if crop_min_b <= r.id[1] <= crop_max_b:
                b_crop.add(r.copy())
        cm.add(b_crop)
        complex_pdb = os.path.join(motif_work_dir, "complex_cropped.pdb")
        save_structure(cx, complex_pdb)

        # ---- Step 1: RFD3 co-adaptation of B around A'
        rfd_pdbs = []
        if coadapt_enabled and execution_mode != 'mock':
            rfd_dir = os.path.join(motif_work_dir, "rfd3_coadapt_out")
            os.makedirs(rfd_dir, exist_ok=True)
            rfd_pdbs = [p for p in collect_pdbs(rfd_dir) if os.path.abspath(p) != os.path.abspath(complex_pdb)]
            if rfd_pdbs and not args.force:
                print(f"  [Resume Cache] Reusing {len(rfd_pdbs)} RFD3 backbones.")
            else:
                print(f"  Running RFD3 backbone co-adaptation on B around {a_cand_name}...")
                with open(os.path.join(args.analysis_dir, "inputs_rescue.json")) as f:
                    rescue_json = json.load(f)
                rescue_json['design']['input'] = os.path.relpath(complex_pdb, rfd_dir).replace('\\', '/')
                inputs_json = os.path.join(rfd_dir, "inputs.json").replace('\\', '/')
                with open(inputs_json, 'w') as f:
                    json.dump(rescue_json, f, indent=2)
                rfd_bin = pcfg.get('local_rfdiffusion', 'OrthoIntRob/bin/rfd3') if execution_mode == 'local' else pcfg['colab_rfdiffusion']
                rfd_cmd = rfd_bin.split() + [
                    f"inputs={inputs_json}", f"out_dir={rfd_dir.replace(chr(92), '/')}",
                    f"n_batches={pcfg.get('rescue_diffusion_n_batches', 10)}", "diffusion_batch_size=1"]
                env = os.environ.copy()
                env['PYTHONUNBUFFERED'] = '1'
                subprocess.run(rfd_cmd, check=True, env=env)
                rfd_pdbs = [p for p in collect_pdbs(rfd_dir) if os.path.abspath(p) != os.path.abspath(complex_pdb)]

        scaffolds = pick_scaffolds(rfd_pdbs, a_chain, wt_B_chain, crop_ids, fixed_b_ids, chain_B_id, pcfg, complex_pdb)
        print(f"  RFD3 produced {len(rfd_pdbs)} backbone(s); {len(scaffolds)} scaffold(s) retained:")
        for i, (p, native, m) in enumerate(scaffolds):
            print(f"    scaffold {i}: {'NATIVE' if native else os.path.basename(p)}  "
                  f"contacts={m['contacts']} clashes={m['clashes']} "
                  f"fw_rmsd={m['framework_rmsd']} flex_rmsd={m['flex_rmsd']}")
        with open(os.path.join(motif_work_dir, "scaffold_report.json"), 'w') as f:
            json.dump([{"pdb": p, "native": n, **m} for p, n, m in scaffolds], f, indent=2)

        # ---- Step 2: LigandMPNN on every scaffold, with specificity-by-difference bias
        pool = []
        per_scaffold = max(15, samples_per_scaffold // max(1, len(scaffolds)))
        for sc_idx, (sc_pdb, is_native, sc_m) in enumerate(scaffolds):
            assembled = os.path.join(motif_work_dir, f"complex_scaffold_{sc_idx}_allatom.pdb")
            b_scaffold = assemble_complex(sc_pdb, a_chain, wt_B_chain, crop_ids, fixed_b_ids,
                                          chain_A_id, chain_B_id, assembled)
            if b_scaffold is None:
                continue

            existing_ids = {f"{ch.id}{r.id[1]}" for ch in PDBParser(QUIET=True).get_structure("c", assembled)[0]
                            for r in std_residues(ch)}
            fixed_list = [f"{ch}{res}" for ch, lst in mpnn_fixed.items() for res in lst if f"{ch}{res}" in existing_ids]
            rescue_bias = build_rescue_bias(a_chain, b_scaffold, wt_A_chain, fixed_b_ids)
            n_spec = sum(1 for v in rescue_bias.values() if len(v) > 0)
            print(f"  Scaffold {sc_idx + 1}/{len(scaffolds)}: bias on {n_spec} B positions, "
                  f"{per_scaffold} LigandMPNN samples...")

            mpnn_out_dir = os.path.join(motif_work_dir, f"mpnn_out_scaffold_{sc_idx}")
            os.makedirs(mpnn_out_dir, exist_ok=True)
            if execution_mode == 'mock':
                continue
            mpnn_bin = pcfg.get('local_ligandmpnn', 'OrthoIntRob/bin/mpnn') if execution_mode == 'local' else pcfg['colab_ligandmpnn']
            mpnn_cmd = mpnn_bin.split() + [
                "--structure_path", assembled.replace('\\', '/'),
                "--out_directory", mpnn_out_dir.replace('\\', '/'),
                "--model_type", model_type,
                "--checkpoint_path", checkpoint_path.replace('\\', '/'),
                "--is_legacy_weights", is_legacy,
                "--batch_size", "1",
                "--number_of_batches", str(per_scaffold),
                "--temperature", temp_rescue,
                "--fixed_residues", ",".join(fixed_list),
                "--bias_per_residue", json.dumps(rescue_bias),
                "--write_structures", "True",
            ]
            subprocess.run(mpnn_cmd, check=True, env=os.environ.copy())

            outs = [f for f in collect_pdbs(mpnn_out_dir)
                    if os.path.abspath(f) != os.path.abspath(assembled)
                    and not f.endswith("B_prime_rfd_mpnn.pdb") and "candidate" not in os.path.basename(f)]
            packed = [f for f in outs if "packed" in f]
            outs = packed or outs

            mutable_ids = {r.id[1] for r in std_residues(b_scaffold) if r.id[1] not in fixed_b_ids}
            for p in outs:
                try:
                    st = PDBParser(QUIET=True).get_structure("o", p)[0]
                    ch_b = st[chain_B_id] if chain_B_id in st else list(st.get_chains())[1]
                except Exception as e:
                    print(f"    Warning: unreadable MPNN output {os.path.basename(p)}: {e}")
                    continue
                spec = candidate_specificity_score(a_chain, wt_A_chain, ch_b)
                seq_mut = "".join(THREE_TO_ONE.get(r.get_resname(), 'X') for r in std_residues(ch_b)
                                  if r.id[1] in mutable_ids)
                pool.append({"path": p, "scaffold": sc_idx, "native": is_native, "seq": seq_mut,
                             "conf": parse_mpnn_confidence(p), "spec": spec,
                             "scaffold_flex_rmsd": sc_m.get("flex_rmsd")})

        if execution_mode == 'mock' or not pool:
            if not pool:
                print(f"  WARNING: no LigandMPNN candidates produced for {a_cand_name}.")
            continue

        # ---- Step 3: rank (negative-design proxy + LigandMPNN confidence) then diversify
        ok = [c for c in pool if c["spec"]["clashes_a_prime"] <= max_clash_a] or pool
        z_spec = zscore([c["spec"]["spec_score"] for c in ok])
        z_conf = zscore([c["conf"] for c in ok])
        for c, zs, zc in zip(ok, z_spec, z_conf):
            c["score"] = zs + zc
        chosen = select_diverse([c["seq"] for c in ok], top_k_per_motif,
                                scores=[c["score"] for c in ok],
                                groups=[c["scaffold"] for c in ok], min_per_group=1)
        final = [ok[i] for i in chosen]
        by_sc = {}
        for c in final:
            by_sc[c["scaffold"]] = by_sc.get(c["scaffold"], 0) + 1
        print(f"  Selected {len(final)}/{len(pool)} B' designs; per-scaffold counts: {by_sc}")

        for b_idx, cand in enumerate(final):
            pair_id = f"{a_cand_name}_B_cand_{b_idx:02d}"
            dest = os.path.join(args.out_dir, f"{pair_id}.pdb")
            st_c = PDBParser(QUIET=True).get_structure("cand", cand["path"])[0]
            ch_b = st_c[chain_B_id] if chain_B_id in st_c else list(st_c.get_chains())[1]
            ch_b = ch_b.copy()
            ch_b.id = chain_B_id
            ps = Structure.Structure(pair_id)
            pm = Model.Model(0)
            ps.add(pm)
            pm.add(a_chain.copy())
            pm.add(ch_b)
            save_structure(ps, dest)
            strip_nan_lines(dest)
            all_saved_b_pdbs.append(dest)
            metadata_pairs[pair_id] = {
                "parent_a": a_cand_name, "pdb": dest, "scaffold_idx": cand["scaffold"],
                "scaffold_is_native": cand["native"], "scaffold_flex_rmsd": cand["scaffold_flex_rmsd"],
                "mpnn_confidence": cand["conf"], **cand["spec"]}
            if a_idx == 0:   # backward-compatible aliases for the primary motif
                shutil.copy(dest, os.path.join(args.out_dir, f"B_prime_candidate_{b_idx:02d}.pdb"))
                if b_idx == 0:
                    shutil.copy(dest, os.path.join(args.out_dir, "B_prime_rfd_mpnn.pdb"))

    with open(os.path.join(args.out_dir, "diversity_pairs_metadata.json"), 'w') as f:
        json.dump(metadata_pairs, f, indent=2)
    print(f"\nModule 4 Complete: {len(all_saved_b_pdbs)} A'/B' pairs across {len(selected_candidates)} A' motifs.")


if __name__ == "__main__":
    main()
