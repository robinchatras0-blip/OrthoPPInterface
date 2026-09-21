import json
import os

import pytest
from Bio.PDB import PDBParser

import folding_engine as fe
from conftest import load_module, make_chain, write_pdb
from dockq import calculate_dockq

WT = "ACDEFGHIKL"


def write_a3m(path, rows):
    with open(path, "w") as f:
        for i, r in enumerate(rows):
            f.write(f">seq{i}\n{r}\n")
    return str(path)


def read_rows(path):
    return [s for _, s in fe.read_a3m(path)]


@pytest.fixture
def a3m(tmp_path):
    return write_a3m(tmp_path / "wt.a3m", [WT, "ACDEFGHIKL", "AGDEFGHVKL", "ac-DEFGHIKL"])


def test_gap_mode_masks_only_requested_columns(a3m, tmp_path):
    out = str(tmp_path / "o.a3m")
    fe.build_hybrid_msa(a3m, "ACWEFGHIKL", [2], out, mode="gap")
    rows = read_rows(out)
    assert rows[0] == "ACWEFGHIKL"
    assert rows[1] == "AC-EFGHIKL" and rows[2] == "AG-EFGHVKL"
    assert all(len(r) == 10 for r in rows[:3])


def test_auto_scope_derives_mutations_from_query(a3m, tmp_path):
    out = str(tmp_path / "o.a3m")
    stats = fe.build_hybrid_msa(a3m, "ACWEFGHIKK", None, out)
    assert stats["n_masked"] == 2                      # only the 2 real mutations, not a whole neighbourhood
    assert read_rows(out)[1] == "AC-EFGHIK-"


def test_substitute_and_keep_modes(a3m, tmp_path):
    out = str(tmp_path / "o.a3m")
    fe.build_hybrid_msa(a3m, "ACWEFGHIKL", None, out, mode="substitute")
    assert read_rows(out)[1] == "ACWEFGHIKL"
    fe.build_hybrid_msa(a3m, "ACWEFGHIKL", None, out, mode="keep")
    assert read_rows(out)[1] == "ACDEFGHIKL"


def test_insertions_next_to_masked_columns_are_dropped(a3m, tmp_path):
    src = write_a3m(tmp_path / "ins.a3m", [WT, "ACDxxEFGHIKL", "ACDEFGHIKlL"])
    out = str(tmp_path / "o.a3m")
    fe.build_hybrid_msa(src, "ACWEFGHIKL", [2], out)
    assert read_rows(out)[1] == "AC-EFGHIKL"           # insertion after masked col 2 removed
    assert read_rows(out)[2] == "AC-EFGHIKlL"          # insertion elsewhere kept


def test_width_mismatch_raises_instead_of_corrupting_msa(a3m, tmp_path):
    with pytest.raises(ValueError, match="MSA width"):
        fe.build_hybrid_msa(a3m, "ACDEFGH", None, str(tmp_path / "o.a3m"))
    stats = fe.build_hybrid_msa(a3m, "ACDEFGH", None, str(tmp_path / "o2.a3m"), strict=False)
    assert stats["fallback"] and read_rows(str(tmp_path / "o2.a3m")) == ["ACDEFGH"]


def test_wrapped_and_commented_a3m_is_parsed(tmp_path):
    p = tmp_path / "w.a3m"
    p.write_text("#10\t1\n>q\nACDEF\nGHIKL\n>h\nACDEF\nGHIKL\n")
    assert read_a3m_rows(p) == [WT, WT]


def read_a3m_rows(p):
    return [s for _, s in fe.read_a3m(str(p))]


def test_prepare_msa_reference_masks_where_the_other_design_differs(a3m, tmp_path):
    cfg = {"folding": {"use_msa": True, "msa_mask_mode": "gap", "msa_mask_scope": "diff"}}
    out = str(tmp_path / "wtmatched.a3m")
    fe.prepare_msa(a3m, WT, out, cfg, reference_sequence="AAAEFGHIKL")   # design differs at cols 1,2
    rows = read_rows(out)
    assert rows[0] == WT and rows[1] == "A--EFGHIKL"


def test_prepare_msa_without_msa_is_single_sequence(a3m, tmp_path):
    out = str(tmp_path / "n.a3m")
    fe.prepare_msa(a3m, WT, out, {"folding": {"use_msa": False}})
    assert read_rows(out) == [WT]


def test_cache_signature_tracks_msa_content(tmp_path):
    m1 = write_a3m(tmp_path / "m.a3m", [WT])
    comps = [{"seq": WT, "chain_id": "A", "msa_path": m1}]
    s1 = fe._signature(comps, {})
    write_a3m(tmp_path / "m.a3m", [WT, "AAAAAAAAAA"])
    assert fe._signature(comps, {}) != s1
    assert fe._signature(comps, {"diffusion_batch_size": 5}) != fe._signature(comps, {})


def test_rf3_output_parsing_picks_best_ranked_sample(tmp_path):
    for name, iptm, rank in (("s0", 0.30, 0.31), ("s1", 0.70, 0.72), ("s2", 0.50, 0.52)):
        (tmp_path / f"x_{name}_summary_confidences.json").write_text(json.dumps(
            {"iptm": iptm, "ptm": 0.8, "ranking_score": rank, "overall_plddt": 0.85, "has_clash": False,
             "chain_pair_pae_min": [[1.0, 4.2], [4.2, 1.0]]}))
    path, m = fe._collect_metrics(str(tmp_path))
    assert "s1" in path and m["iptm"] == 0.70 and m["n_samples"] == 3
    assert m["plddt"] == pytest.approx(85.0) and m["chain_pair_pae_min"] == 4.2
    assert m["iptm_mean"] == pytest.approx(0.5) and m["iptm_std"] > 0


def test_stale_cache_is_recomputed_when_msa_changes(tmp_path, monkeypatch):
    calls = []

    class R:
        returncode, stdout, stderr = 0, "", ""

    def fake_run(cmd, **kw):
        calls.append(cmd)
        out = [c.split("=", 1)[1] for c in cmd if c.startswith("out_dir=")][0]
        with open(os.path.join(out, "t_summary_confidences.json"), "w") as f:
            json.dump({"iptm": 0.5, "ptm": 0.5, "ranking_score": 0.5, "overall_plddt": 80}, f)
        return R()

    monkeypatch.setattr(fe.subprocess, "run", fake_run)
    cfg = {"folding": {"use_wsl": False, "use_msa": True}}
    msa = write_a3m(tmp_path / "m.a3m", [WT])
    out = str(tmp_path / "t")
    fe.predict_structure([("A", WT, msa)], out, cfg)
    fe.predict_structure([("A", WT, msa)], out, cfg)
    assert len(calls) == 1                               # cache hit
    write_a3m(tmp_path / "m.a3m", [WT, "AAAAAAAAAA"])
    fe.predict_structure([("A", WT, msa)], out, cfg)
    assert len(calls) == 2                               # MSA changed -> recomputed, not the old score


def test_dockq_is_numbering_agnostic(helix_pair, tmp_path):
    path, seq_a, seq_b = helix_pair
    model = write_pdb(tmp_path / "renum.pdb", [make_chain("A", seq_a, first_id=1),
                                                make_chain("B", seq_b, first_id=1, origin=(5.0, 0, 0))])
    native = write_pdb(tmp_path / "shifted.pdb", [make_chain("A", seq_a, first_id=40),
                                                   make_chain("B", seq_b, first_id=77, origin=(5.0, 0, 0))])
    res = calculate_dockq(native, model)
    assert res["dockq"] > 0.99 and res["quality"] == "High"


def test_module4_assembles_full_length_b_with_head_and_tail(tmp_path):
    m4 = load_module("04_generate_B_prime.py")
    names = ["ALA", "LEU", "LYS", "GLU", "ASP", "VAL", "ILE", "SER", "THR", "ARG"] * 2      # 20 residues
    wt_b = make_chain("B", names, origin=(5, 0, 0))
    a = make_chain("A", names)
    crop_ids = list(range(5, 16))
    # scaffold = cropped B, renumbered from 1 the way some diffusion outputs are
    scaffold = write_pdb(tmp_path / "sc.pdb", [make_chain("A", names), make_chain("B", names[4:15], first_id=1, origin=(5, 0, 0))])
    out = str(tmp_path / "assembled.pdb")
    fixed = {5, 6, 7, 8, 12, 13, 14, 15}
    new_b = m4.assemble_complex(scaffold, a, wt_b, crop_ids, fixed, "A", "B", out)
    ids = [r.id[1] for r in new_b]
    assert ids == list(range(1, 21))                       # N-terminal head (1-4) and tail (16-20) restored
    st = PDBParser(QUIET=True).get_structure("x", out)[0]
    assert len(list(st["B"])) == 20 and len(list(st["A"])) == 20
    assert len({r.id[1] for r in wt_b}) == 20              # WT chain untouched (no in-place drift)


def test_module4_pick_scaffolds_drops_clashing_and_duplicate_backbones(tmp_path):
    m4 = load_module("04_generate_B_prime.py")
    names = ["ALA", "LEU", "LYS", "GLU", "ASP", "VAL", "ILE", "SER", "THR", "ARG"] * 2
    a = make_chain("A", names)
    wt_b = make_chain("B", names, origin=(5, 0, 0))
    native = write_pdb(tmp_path / "native.pdb", [a, wt_b])
    good = write_pdb(tmp_path / "good.pdb", [make_chain("A", names), make_chain("B", names, origin=(5, 0, 0))])
    dup = write_pdb(tmp_path / "dup.pdb", [make_chain("A", names), make_chain("B", names, origin=(5, 0, 0))])
    clash = write_pdb(tmp_path / "clash.pdb", [make_chain("A", names), make_chain("B", names, origin=(0.3, 0, 0))])
    cfg = {"max_rescue_scaffolds": 4, "include_native_scaffold": False, "scaffold_max_clashes": 5}
    fixed = set(range(1, 9))
    picked = m4.pick_scaffolds([good, dup, clash], a, wt_b, list(range(1, 21)), fixed, "B", cfg, native)
    paths = [os.path.basename(p) for p, native_flag, _ in picked]
    assert "clash.pdb" not in paths and len(paths) == 1        # clash rejected, duplicate collapsed
    cfg["include_native_scaffold"] = True
    picked = m4.pick_scaffolds([good], a, wt_b, list(range(1, 21)), fixed, "B", cfg, native)
    assert picked[0][1] is True                                 # native is opt-in and flagged


def test_module4_align_scaffolds_restores_a_prime_frame(tmp_path):
    import numpy as np
    from design_utils import interface_stats
    m4 = load_module("04_generate_B_prime.py")
    names = ["ALA", "LEU", "LYS", "GLU", "ASP", "VAL", "ILE", "SER", "THR", "ARG"] * 2
    a = make_chain("A", names)
    b = make_chain("B", names, origin=(9, 0, 0))
    ref_contacts, _ = interface_stats(a, b)
    a2, b2 = make_chain("A", names), make_chain("B", names, origin=(9, 0, 0))
    ang = 1.3
    rot = np.array([[np.cos(ang), 0, np.sin(ang)], [0, 1, 0], [-np.sin(ang), 0, np.cos(ang)]])
    for at in list(a2.get_atoms()) + list(b2.get_atoms()):      # RFD3-style re-centring
        at.coord = at.coord @ rot + np.array([-30.0, 70.0, 5.0])
    moved = write_pdb(tmp_path / "moved.pdb", [a2, b2])
    assert interface_stats(a, PDBParser(QUIET=True).get_structure("m", moved)[0]["B"])[0] == 0
    out = m4.align_scaffolds([moved], a, "A", str(tmp_path / "aligned"))
    assert len(out) == 1
    b_aligned = PDBParser(QUIET=True).get_structure("x", out[0])[0]["B"]
    assert interface_stats(a, b_aligned)[0] == ref_contacts > 0


def test_rf3_top_level_summary_does_not_double_count_best_sample(tmp_path):
    def write(path, iptm, rank):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"iptm": iptm, "ptm": 0.8, "ranking_score": rank, "overall_plddt": 0.8}))
    write(tmp_path / "x" / "x_summary_confidences.json", 0.9, 0.95)               # duplicate of the best sample
    for i, v in enumerate((0.9, 0.5, 0.1)):
        write(tmp_path / "x" / f"seed-0_sample-{i}" / f"x_seed-0_sample-{i}_summary_confidences.json", v, v + 0.05)
    _, m = fe._collect_metrics(str(tmp_path))
    assert m["n_samples"] == 3 and m["iptm"] == 0.9
    assert m["iptm_mean"] == pytest.approx(0.5)


def test_interface_scope_masks_interface_columns_for_wt_and_designs(a3m, tmp_path):
    cfg = {"folding": {"use_msa": True, "msa_mask_mode": "gap", "msa_mask_scope": "interface"}}
    out = str(tmp_path / "o.a3m")
    st = fe.prepare_msa(a3m, WT, out, cfg, interface_columns=[0, 5])          # WT chain: interface only
    assert st["n_masked"] == 2 and read_rows(out)[0] == WT and read_rows(out)[1] == "-CDEF-HIKL"
    st = fe.prepare_msa(a3m, "ACWEFGHIKL", out, cfg, interface_columns=[5])   # design: mutated column + interface
    assert st["n_masked"] == 2 and read_rows(out)[1] == "AC-EF-HIKL"
    st = fe.prepare_msa(a3m, WT, out, cfg, reference_sequence="AAAEFGHIKL", interface_columns=[5])
    assert st["n_masked"] == 3                                                # matched regime: columns 1, 2 and 5


def test_interface_scope_requires_columns_and_scope_is_validated(a3m, tmp_path):
    out = str(tmp_path / "o.a3m")
    with pytest.raises(ValueError, match="interface columns"):
        fe.prepare_msa(a3m, WT, out, {"folding": {"msa_mask_scope": "interface"}})
    with pytest.raises(ValueError, match="Unknown folding.msa_mask_scope"):
        fe.prepare_msa(a3m, WT, out, {"folding": {"msa_mask_scope": "everything"}})


def test_interface_columns_fall_back_for_mappings_without_the_new_keys(helix_pair):
    from design_utils import interface_columns
    path, _, _ = helix_pair
    old_mapping = {"residues_A": [{"resseq": i, "is_interface": i <= 3} for i in range(1, 21)]}
    cols_a, cols_b = interface_columns(path, "A", "B", old_mapping, fixed_b_ids=list(range(1, 17)))
    assert cols_a == [0, 1, 2] and cols_b == [16, 17, 18, 19]                 # B interface = residues left mutable
    new_mapping = {"residues_A": [], "interface_ids_A": [5], "interface_ids_B": [2, 3]}
    assert interface_columns(path, "A", "B", new_mapping) == ([4], [1, 2])
