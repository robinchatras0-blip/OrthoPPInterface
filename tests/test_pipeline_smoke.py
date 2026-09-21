"""End-to-end wiring test for modules 3 -> 4 -> 5 with RFD3 / LigandMPNN / RF3 replaced by fakes.

It cannot judge scientific quality, but it executes every line of the three main() functions on
synthetic geometry and asserts the properties that were broken before: scaffold diversity survives
selection, WT head/tail are restored, no fabricated metrics, matched MSA regimes.
"""
import json
import os
import random
import shutil
import sys

import numpy as np
import pandas as pd
import pytest
import yaml
from Bio.PDB import PDBParser

from conftest import load_module, make_chain, write_pdb
from design_utils import THREE_TO_ONE

NAMES_A = ["ALA", "LEU", "LYS", "GLU", "ASP", "VAL", "ILE", "SER", "THR", "ARG"] * 2
NAMES_B = ["GLU", "ASP", "LYS", "ILE", "VAL", "ARG", "LEU", "GLN", "ASN", "SER"] * 2
ONE = lambda names: "".join(THREE_TO_ONE[n] for n in names)  # noqa: E731
B_MUTABLE = {9, 10, 11, 13}


@pytest.fixture
def workspace(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    data = tmp_path / "data"
    data.mkdir()
    wt = write_pdb(data / "wt.pdb", [make_chain("A", NAMES_A), make_chain("B", NAMES_B, origin=(9, 0, 0))])
    for tag, names in (("A", NAMES_A), ("B", NAMES_B)):
        (data / f"{tag}.a3m").write_text(f">q\n{ONE(names)}\n>h1\n{ONE(names)}\n>h2\n{ONE(names)[:-1]}L\n")

    ana = tmp_path / "run" / "01_analysis"
    ana.mkdir(parents=True)
    json.dump({"crop_min_B": 5, "crop_max_B": 15, "neighborhood_ids": list(range(1, 11)),
               "fixed_ids_A": list(range(11, 21)), "interface_ids_A": list(range(1, 7)),
               "interface_ids_B": sorted(B_MUTABLE),
               "residues_A": [{"sequential_index": i + 1, "is_interface": i < 6, "is_neighborhood": i < 10}
                              for i in range(20)]}, open(ana / "index_mapping.json", "w"))
    json.dump({"A": list(range(1, 21)), "B": [i for i in range(1, 21) if i not in B_MUTABLE]},
              open(ana / "mpnn_fixed_positions_B.json", "w"))
    json.dump({"design": {"input": "x.pdb", "contig": "A1-20,/0,B5-15"}}, open(ana / "inputs_rescue.json", "w"))
    json.dump({"design": {"input": "x.pdb", "contig": "A1-20,/0,B5-15"}}, open(ana / "inputs.json", "w"))
    json.dump({"A": list(range(11, 21)), "B": list(range(1, 21))}, open(ana / "mpnn_fixed_positions.json", "w"))
    json.dump({"A": {"3": {"ASP": 2.0}}}, open(ana / "mpnn_bias.json", "w"))

    d2 = tmp_path / "run" / "02_rupture_design"
    d2.mkdir()
    cands = []
    for i, muts in enumerate(({2: "GLU", 3: "LYS"}, {4: "ARG", 6: "TRP"})):
        names = list(NAMES_A)
        for k, v in muts.items():
            names[k] = v
        cands.append(write_pdb(d2 / f"A_prime_candidate_{i:02d}.pdb", [make_chain("A", names)]))
    d3 = tmp_path / "run" / "03_fail_fast"
    d3.mkdir()
    (d3 / "passed_candidates.txt").write_text("\n".join(cands) + "\n")
    for i in range(2):
        json.dump({"plddt_monomer": 88.0, "rmsd_monomer": 0.9, "iptm_rupture": 0.2}, open(d3 / f"metrics_A_prime_candidate_{i:02d}.json", "w"))

    cfg = {
        "pipeline": {"input_pdb": str(wt), "input_msa_A": str(data / "A.a3m"),
                     "input_msa_B": str(data / "B.a3m"), "chain_A": "A", "chain_B": "B",
                     "local_rfdiffusion": "rfd3", "local_ligandmpnn": "mpnn", "rescue_diffusion_n_batches": 4,
                     "max_rescue_scaffolds": 3, "include_native_scaffold": False},
        "diversity_selection": {"enabled": True, "max_a_prime_motifs": 2},
        "ligandmpnn": {"samples_per_scaffold": 15, "top_k_per_motif": 5, "temperature_rescue": 0.2},
        "folding": {"use_msa": True, "msa_mask_scope": "interface", "msa_mask_mode": "gap",
                    "negative_msa_regime": "matched", "wt_ceiling_control": True},
        "structural_constraints": {"interface_distance_threshold_angstroms": 6.0, "neighborhood_radius_angstroms": 10.0,
                                   "crop_chain_B_distance_angstroms": 15.0, "max_hotspot_mutations": 3},
        "thresholds": {"monomer_stability": {"plddt_min": 80, "rmsd_max": 2.0},
                       "negative_design_rupture": {"iptm_rupture_max": 0.35},
                       "positive_design_rescue": {"iptm_rescue_min": 0.75, "relative_to_ceiling": 0.85},
                       "fail_fast": {"max_passing_candidates": 12}},
    }
    cfg_path = tmp_path / "config.yaml"
    yaml.safe_dump(cfg, open(cfg_path, "w"))
    return {"root": tmp_path, "cfg": str(cfg_path), "ana": str(ana), "wt": wt, "d3": str(d3)}


def fake_subprocess_run(cmd, **kw):
    rng = random.Random(len(cmd) + hash(cmd[-1]) % 1000)
    if any(str(c).startswith("inputs=") for c in cmd):                        # fake RFD3
        out = [c.split("=", 1)[1] for c in cmd if str(c).startswith("out_dir=")][0]
        inp = [c.split("=", 1)[1] for c in cmd if str(c).startswith("inputs=")][0]
        src = os.path.join(os.path.dirname(inp), json.load(open(inp))["design"]["input"])
        model = PDBParser(QUIET=True).get_structure("s", src)[0]
        for n in range(4):
            b = model["B"].copy()
            for r in b:
                if r.id[1] in (9, 10, 11):                                   # flexible loop actually moves
                    for at in r:
                        at.coord = at.coord + [0.9 * (n + 1), 0.4 * n, 0.0]
            a = model["A"].copy()
            ang = 0.9 + n
            rot = np.array([[np.cos(ang), -np.sin(ang), 0], [np.sin(ang), np.cos(ang), 0], [0, 0, 1]])
            for at in list(a.get_atoms()) + list(b.get_atoms()):     # real RFD3 re-centres its output (~70 A)
                at.coord = at.coord @ rot + np.array([-40.0, 55.0, 12.0])
            write_pdb(os.path.join(out, f"rfd_{n}.pdb"), [a, b])
    else:                                                                     # fake LigandMPNN
        get = lambda flag: cmd[cmd.index(flag) + 1]  # noqa: E731
        model = PDBParser(QUIET=True).get_structure("s", get("--structure_path"))[0]
        for i in range(int(get("--number_of_batches"))):
            b = model["B"].copy()
            for r in b:
                if r.id[1] in B_MUTABLE:
                    r.resname = rng.choice(list(THREE_TO_ONE))
            write_pdb(os.path.join(get("--out_directory"), f"design_{i}_packed.pdb"), [model["A"].copy(), b])

    class R:
        returncode, stdout, stderr = 0, "", ""
    return R()


def fake_predict(fasta_sequences, out_dir, config):
    labels = tuple(x[0] for x in fasta_sequences)
    for item in fasta_sequences:                                              # MSAs must be well-formed
        if len(item) > 2 and item[2] and os.path.exists(item[2]):
            rows = [l.strip() for l in open(item[2]) if not l.startswith(">") and l.strip()]
            assert all(len(r) == len(item[1]) for r in rows), f"corrupt MSA {item[2]}"
    iptm = {("A_prime", "B_prime"): 0.8, ("A_wt", "B_prime"): 0.2, ("A_prime", "B_wt"): 0.25,
            ("A_wt", "B_wt"): 0.85}.get(labels, 0.9)
    return {"iptm": iptm, "ptm": 0.8, "plddt": 88.0, "ranking_score": iptm, "iptm_mean": iptm,
            "iptm_std": 0.0, "n_samples": 1, "has_clash": False, "chain_pair_pae_min": None, "model_pdb": None}


def run_main(module, argv, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["x"] + argv)
    module.main()


def test_modules_4_and_5_end_to_end(workspace, monkeypatch):
    ws = workspace
    m4 = load_module("04_generate_B_prime.py")
    monkeypatch.setattr(m4.subprocess, "run", fake_subprocess_run)
    out4 = str(ws["root"] / "run" / "04_rescue_design")
    run_main(m4, ["--config", ws["cfg"], "--analysis_dir", ws["ana"], "--passed_candidates",
                  os.path.join(ws["d3"], "passed_candidates.txt"), "--out_dir", out4], monkeypatch)

    meta = json.load(open(os.path.join(out4, "diversity_pairs_metadata.json")))
    assert len(meta) >= 6                                                     # 2 motifs x top_k 5 (some dedup)
    per_motif = {}
    for pid, m in meta.items():
        per_motif.setdefault(m["parent_a"], set()).add(m["scaffold_idx"])
        b = PDBParser(QUIET=True).get_structure("x", m["pdb"])[0]["B"]
        assert [r.id[1] for r in b] == list(range(1, 21))                     # head + tail restored
        assert m["scaffold_is_native"] is False
        assert m["contacts_a_prime"] > 0                                      # B' is designed in A''s frame
    # the bug: previously every retained design came from scaffold 0
    assert all(len(s) >= 2 for s in per_motif.values()), per_motif

    m5 = load_module("05_eval_final.py")
    monkeypatch.setattr(m5, "predict_structure", fake_predict)
    out5 = str(ws["root"] / "run" / "05_final_eval")
    run_main(m5, ["--config", ws["cfg"], "--analysis_dir", ws["ana"], "--design_dir", out4,
                  "--filter_dir", ws["d3"], "--out_dir", out5], monkeypatch)
    df = pd.read_csv(os.path.join(out5, "orthogonality_scores.csv"))
    assert len(df) == len(meta)
    assert df[["iptm_rescue", "iptm_rupture", "iptm_negative", "iptm_ceiling"]].notna().all().all()
    assert (df["f_ortho_iptm"].round(3) == 0.55).all()                       # 0.8 - max(0.25, 0.2)
    assert (df["f_ortho"].round(3) == round(0.55 / 0.85, 3)).all()            # energy stage off: F_iptm / ceiling
    assert df["passes"].all() and df["scaffold_idx"].nunique() >= 2
    assert (df["plddt_a_prime"] == 88.0).all()                                # real Module-3 metric, not a default
    assert df["n_mut_B"].max() > 0 and df["n_mut_A"].min() > 0
    assert df["rmsd_b_prime"].max() > 0                                       # B' is no longer pinned to the WT backbone


def test_module5_never_invents_missing_metrics(workspace, monkeypatch):
    ws = workspace
    os.remove(os.path.join(ws["d3"], "metrics_A_prime_candidate_00.json"))
    d4 = ws["root"] / "run" / "04_rescue_design"
    d4.mkdir()
    write_pdb(d4 / "A_prime_candidate_00_B_cand_00.pdb",
              [make_chain("A", NAMES_A), make_chain("B", NAMES_B, origin=(9, 0, 0))])
    m5 = load_module("05_eval_final.py")
    monkeypatch.setattr(m5, "predict_structure", fake_predict)
    out5 = str(ws["root"] / "run" / "05_final_eval")
    run_main(m5, ["--config", ws["cfg"], "--analysis_dir", ws["ana"], "--design_dir", str(d4),
                  "--filter_dir", ws["d3"], "--out_dir", out5], monkeypatch)
    df = pd.read_csv(os.path.join(out5, "orthogonality_scores.csv"))
    assert pd.isna(df.loc[0, "plddt_a_prime"])                                # NaN, not the old 90.0 default


def test_module3_filters_and_writes_metrics(workspace, monkeypatch):
    ws = workspace
    m3 = load_module("03_filter_A_prime.py")
    monkeypatch.setattr(m3, "predict_structure", fake_predict)
    out3 = str(ws["root"] / "run" / "03_out")
    run_main(m3, ["--config", ws["cfg"], "--design_dir", str(ws["root"] / "run" / "02_rupture_design"),
                  "--analysis_dir", ws["ana"], "--out_dir", out3], monkeypatch)
    assert json.load(open(os.path.join(out3, "wt_control.json")))["iptm_full_msa"] == 0.85
    passed = open(os.path.join(out3, "passed_candidates.txt")).read().split()
    assert len(passed) == 2
    m = json.load(open(os.path.join(out3, "metrics_A_prime_candidate_00.json")))
    assert m["passed"] and m["n_masked_columns"] == 6            # 2 mutations + the 6 interface columns (2 and 3 overlap)
    assert m["msa_scope"] == "interface"                          # lets Module 5 know which regime it can reuse
    assert json.load(open(os.path.join(out3, "wt_control.json")))["msa_scope"] == "interface"


def test_calibration_script_runs(workspace, monkeypatch):
    ws = workspace
    cal = load_module("calibrate_predictor.py")
    monkeypatch.setattr(cal, "predict_structure", fake_predict)
    out = str(ws["root"] / "cal")
    run_main(cal, ["--config", ws["cfg"], "--analysis_dir", ws["ana"], "--out_dir", out,
                   "--n_controls", "2", "--n_pos_mut", "3", "--n_neg_mut", "5"], monkeypatch)
    df = pd.read_csv(os.path.join(out, "calibration_summary.csv"))
    assert len(df) == 6 and {"gap", "iptm_pos_mean", "iptm_neg_mean"} <= set(df.columns)   # 6 regimes


def fake_mpnn_cif_only(cmd, n_mutations=3, expect_chain_b=False):
    """LigandMPNN (Foundry) writes CIF files only; the designed sequence differs from the input backbone."""
    from Bio.PDB import MMCIFIO
    get = lambda flag: cmd[cmd.index(flag) + 1]  # noqa: E731
    st = PDBParser(QUIET=True).get_structure("s", get("--structure_path"))
    assert ("B" in st[0]) == expect_chain_b, "chain B must be absent from the A' design input unless rupture_apo_design is off"
    for r in [r for r in st[0]["A"] if r.id[1] in (2, 3, 4)][:n_mutations]:
        r.resname = "TRP"
    io = MMCIFIO()
    io.set_structure(st)
    io.save(os.path.join(get("--out_directory"), "design_b0_d0.cif"))


def test_module1_writes_all_analysis_files(workspace, monkeypatch):
    m1 = load_module("01_analyze_interface.py")
    out = str(workspace["root"] / "analysis_new")
    run_main(m1, ["--config", workspace["cfg"], "--out_dir", out], monkeypatch)
    for name in ("inputs.json", "inputs_rescue.json", "mpnn_fixed_positions.json", "mpnn_fixed_positions_B.json",
                 "mpnn_bias.json", "index_mapping.json"):
        assert os.path.exists(os.path.join(out, name)), name
    mapping = json.load(open(os.path.join(out, "index_mapping.json")))
    assert mapping["neighborhood_ids"] and len(mapping["residues_A"]) == 20
    assert mapping["interface_ids_A"] and mapping["interface_ids_B"]                   # used by the MSA regime


def test_module2_uses_mpnn_design_not_raw_rfd3_backbone(workspace, monkeypatch):
    """Regression: mpnn writes CIF only and the old code silently kept the raw RFD3 backbone."""
    m2 = load_module("02_generate_A_prime.py")
    root = workspace["root"]
    wt = workspace["wt"]

    def run(cmd, **kw):
        if any(str(c).startswith("inputs=") for c in cmd):                 # fake RFD3: backbone = WT A + B
            out = [c.split("=", 1)[1] for c in cmd if str(c).startswith("out_dir=")][0]
            shutil.copy(wt, os.path.join(out, "inputs_design_0_model_0.pdb"))
        else:
            fake_mpnn_cif_only(cmd)

        class R:
            returncode = 0
        return R()

    monkeypatch.setattr(m2.subprocess, "run", run)
    out2 = str(root / "run" / "02_new")
    run_main(m2, ["--config", workspace["cfg"], "--analysis_dir", workspace["ana"], "--out_dir", out2], monkeypatch)
    cand = PDBParser(QUIET=True).get_structure("c", os.path.join(out2, "A_prime_candidate_00.pdb"))[0]["A"]
    seq = "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in cand)
    assert seq.count("W") >= 3 and seq != ONE(NAMES_A)                      # the MPNN design, not the RFD3 input


def test_module2_fails_loudly_when_mpnn_writes_nothing(workspace, monkeypatch):
    m2 = load_module("02_generate_A_prime.py")
    wt = workspace["wt"]

    def run(cmd, **kw):
        if any(str(c).startswith("inputs=") for c in cmd):
            out = [c.split("=", 1)[1] for c in cmd if str(c).startswith("out_dir=")][0]
            shutil.copy(wt, os.path.join(out, "inputs_design_0_model_0.pdb"))

        class R:
            returncode = 0
        return R()

    monkeypatch.setattr(m2.subprocess, "run", run)
    with pytest.raises(RuntimeError, match="LigandMPNN wrote no structure"):
        run_main(m2, ["--config", workspace["cfg"], "--analysis_dir", workspace["ana"],
                      "--out_dir", str(workspace["root"] / "run" / "02_empty")], monkeypatch)


def test_module2_can_keep_chain_b_in_the_mpnn_context(workspace, monkeypatch):
    m2 = load_module("02_generate_A_prime.py")
    cfg = yaml.safe_load(open(workspace["cfg"]))
    cfg["ligandmpnn"]["rupture_apo_design"] = False
    yaml.safe_dump(cfg, open(workspace["cfg"], "w"))
    wt = workspace["wt"]

    def run(cmd, **kw):
        if any(str(c).startswith("inputs=") for c in cmd):
            out = [c.split("=", 1)[1] for c in cmd if str(c).startswith("out_dir=")][0]
            shutil.copy(wt, os.path.join(out, "inputs_design_0_model_0.pdb"))
        else:
            fake_mpnn_cif_only(cmd, expect_chain_b=True)

        class R:
            returncode = 0
        return R()

    monkeypatch.setattr(m2.subprocess, "run", run)
    run_main(m2, ["--config", workspace["cfg"], "--analysis_dir", workspace["ana"],
                  "--out_dir", str(workspace["root"] / "run" / "02_holo")], monkeypatch)


def with_energy(workspace, **extra):
    cfg = yaml.safe_load(open(workspace["cfg"]))
    cfg["energy"] = {"enabled": True, "pack_designs": True, "weight_iptm": 0.5, "min_binding_fraction": 0.5,
                     "f_energy_min": 0.0, **extra}
    yaml.safe_dump(cfg, open(workspace["cfg"], "w"))


def run_module2(workspace, monkeypatch, packed):
    m2 = load_module("02_generate_A_prime.py")
    wt = workspace["wt"]

    def run(cmd, **kw):
        if any(str(c).startswith("inputs=") for c in cmd):
            out = [c.split("=", 1)[1] for c in cmd if str(c).startswith("out_dir=")][0]
            shutil.copy(wt, os.path.join(out, "inputs_design_0_model_0.pdb"))
        else:
            fake_mpnn_cif_only(cmd)

        class R:
            returncode = 0
        return R()

    monkeypatch.setattr(m2.subprocess, "run", run)
    monkeypatch.setattr(m2, "pack_structures", lambda jobs, config, status: packed.append(jobs))
    run_main(m2, ["--config", workspace["cfg"], "--analysis_dir", workspace["ana"],
                  "--out_dir", str(workspace["root"] / "run" / "02_pack")], monkeypatch)


def test_module2_packs_a_prime_designs_in_place_when_energy_is_enabled(workspace, monkeypatch):
    with_energy(workspace)
    packed = []
    run_module2(workspace, monkeypatch, packed)
    assert len(packed) == 1 and [j["name"] for j in packed[0]] == ["A_prime_candidate_00"]
    assert packed[0][0]["pdb"] == packed[0][0]["out"]


def test_module2_does_not_pack_without_the_energy_stage(workspace, monkeypatch):
    packed = []
    run_module2(workspace, monkeypatch, packed)
    assert packed == []


def test_module4_packs_b_prime_against_the_fixed_a_prime(workspace, monkeypatch):
    with_energy(workspace)
    ws = workspace
    m4 = load_module("04_generate_B_prime.py")
    monkeypatch.setattr(m4.subprocess, "run", fake_subprocess_run)
    packed = []
    monkeypatch.setattr(m4, "pack_structures", lambda jobs, config, status: packed.append(jobs))
    out4 = str(ws["root"] / "run" / "04_pack")
    run_main(m4, ["--config", ws["cfg"], "--analysis_dir", ws["ana"], "--passed_candidates",
                  os.path.join(ws["d3"], "passed_candidates.txt"), "--out_dir", out4], monkeypatch)
    assert len(packed) == 2                                                    # one batch per A' motif
    for jobs in packed:
        assert len(jobs) >= 6 and all(j["fixed_chains"] == ["A"] and j["pdb"] == j["out"] for j in jobs)
        assert all(os.path.exists(j["pdb"]) for j in jobs)
    meta = json.load(open(os.path.join(out4, "diversity_pairs_metadata.json")))
    assert meta and all(os.path.exists(m["pdb"]) for m in meta.values())     # final designs = packed pool complexes


def test_module5_adds_energy_columns_and_the_combined_score(workspace, monkeypatch):
    with_energy(workspace)
    ws = workspace
    d4 = ws["root"] / "run" / "04_energy"
    d4.mkdir()
    write_pdb(d4 / "A_prime_candidate_00_B_cand_00.pdb", [make_chain("A", NAMES_A), make_chain("B", NAMES_B, origin=(9, 0, 0))])
    m5 = load_module("05_eval_final.py")
    monkeypatch.setattr(m5, "predict_structure", fake_predict)
    seen = []

    def fake_energy(jobs, out_json, config):
        seen.append([j["name"] for j in jobs])
        dg = {"native": -60.0, "Ap.Bp": -45.0, "AWT.Bp": -30.0, "Ap.BWT": -36.0}
        return {j["name"]: {"dG": dg["native" if j["name"] == "native" else j["name"].split("__")[1]], "dSASA": 2000.0,
                            "unsat_hb": 12} for j in jobs}

    monkeypatch.setattr(m5, "score_interfaces", fake_energy)
    out5 = str(ws["root"] / "run" / "05_energy")
    run_main(m5, ["--config", ws["cfg"], "--analysis_dir", ws["ana"], "--design_dir", str(d4),
                  "--filter_dir", ws["d3"], "--out_dir", out5], monkeypatch)
    assert seen[0] == ["native", "A_prime_candidate_00_B_cand_00__Ap.Bp", "A_prime_candidate_00_B_cand_00__AWT.Bp",
                       "A_prime_candidate_00_B_cand_00__Ap.BWT"]
    assert os.path.exists(os.path.join(out5, "energy", "A_prime_candidate_00_B_cand_00", "AWT.Bp.pdb"))
    row = pd.read_csv(os.path.join(out5, "orthogonality_scores.csv")).iloc[0]
    assert row["b_rescue"] == pytest.approx(0.75) and row["f_energy"] == pytest.approx(0.15)
    assert row["dG_native"] == -60.0 and row["dsasa_rescue"] == 2000.0
    f_i = 0.55 / 0.85
    assert row["f_iptm_rel"] == pytest.approx(f_i) and row["f_ortho"] == pytest.approx(0.5 * f_i + 0.5 * 0.15)
    assert bool(row["pass_energy"]) and bool(row["passes"])
    assert os.path.exists(os.path.join(out5, "coherence.json"))


def test_module5_energy_can_veto_a_design_that_rf3_likes(workspace, monkeypatch):
    with_energy(workspace)
    ws = workspace
    d4 = ws["root"] / "run" / "04_veto"
    d4.mkdir()
    write_pdb(d4 / "A_prime_candidate_00_B_cand_00.pdb", [make_chain("A", NAMES_A), make_chain("B", NAMES_B, origin=(9, 0, 0))])
    m5 = load_module("05_eval_final.py")
    monkeypatch.setattr(m5, "predict_structure", fake_predict)
    dg = {"native": -60.0, "Ap.Bp": -45.0, "AWT.Bp": -55.0, "Ap.BWT": -36.0}     # B' binds A_WT better than A'
    monkeypatch.setattr(m5, "score_interfaces", lambda jobs, out_json, config: {
        j["name"]: {"dG": dg["native" if j["name"] == "native" else j["name"].split("__")[1]]} for j in jobs})
    out5 = str(ws["root"] / "run" / "05_veto")
    run_main(m5, ["--config", ws["cfg"], "--analysis_dir", ws["ana"], "--design_dir", str(d4),
                  "--filter_dir", ws["d3"], "--out_dir", out5], monkeypatch)
    row = pd.read_csv(os.path.join(out5, "orthogonality_scores.csv")).iloc[0]
    assert row["pass_rescue"] and row["pass_rupture"] and row["pass_negative"]        # RF3 alone would accept it
    assert row["f_energy"] < 0 and not row["pass_energy"] and not row["passes"]       # ...the energy does not


# --------------------------------------------------------------------------- no repeated work in module 5
def recording_predict(calls):
    def predict(fasta_sequences, out_dir, config):
        calls.append(tuple(x[0] for x in fasta_sequences))
        return fake_predict(fasta_sequences, out_dir, config)
    return predict


def set_cfg(workspace, **sections):
    cfg = yaml.safe_load(open(workspace["cfg"]))
    for k, v in sections.items():
        cfg.setdefault(k, {}).update(v)
    yaml.safe_dump(cfg, open(workspace["cfg"], "w"))


def run_module5_on_two_designs(workspace, monkeypatch, name, calls):
    ws = workspace
    d4 = ws["root"] / "run" / name
    d4.mkdir()
    for i in range(2):
        write_pdb(d4 / f"A_prime_candidate_00_B_cand_{i:02d}.pdb", [make_chain("A", NAMES_A), make_chain("B", NAMES_B, origin=(9, 0, 0))])
    m5 = load_module("05_eval_final.py")
    monkeypatch.setattr(m5, "predict_structure", recording_predict(calls))
    out5 = str(ws["root"] / "run" / (name + "_out"))
    run_main(m5, ["--config", ws["cfg"], "--analysis_dir", ws["ana"], "--design_dir", str(d4),
                  "--filter_dir", ws["d3"], "--out_dir", out5], monkeypatch)
    return m5, pd.read_csv(os.path.join(out5, "orthogonality_scores.csv"))


def test_module5_reuses_module3_rupture_and_ceiling_instead_of_refolding(workspace, monkeypatch):
    set_cfg(workspace, folding={"wt_ceiling_control": "global", "reuse_module3_rupture": True})
    json.dump({"plddt_monomer": 88.0, "rmsd_monomer": 0.9, "iptm_rupture": 0.20, "msa_scope": "interface"},
              open(os.path.join(workspace["d3"], "metrics_A_prime_candidate_00.json"), "w"))
    json.dump({"iptm_full_msa": 0.86, "iptm_regime": 0.60, "msa_scope": "interface"},
              open(os.path.join(workspace["d3"], "wt_control.json"), "w"))
    calls = []
    _, df = run_module5_on_two_designs(workspace, monkeypatch, "reuse", calls)
    assert ("A_prime", "B_wt") not in calls and ("A_wt", "B_wt") not in calls              # neither refolded
    assert calls.count(("A_prime", "B_prime")) == 2 and calls.count(("A_wt", "B_prime")) == 2  # what depends on B' still is
    assert (df["iptm_rupture"] == 0.20).all() and (df["iptm_ceiling"] == 0.60).all()


def test_module5_computes_the_global_ceiling_once_and_only_reuses_a_matching_regime(workspace, monkeypatch):
    set_cfg(workspace, folding={"wt_ceiling_control": "global", "reuse_module3_rupture": True})
    json.dump({"plddt_monomer": 88.0, "rmsd_monomer": 0.9, "iptm_rupture": 0.20, "msa_scope": "diff"},      # other regime
              open(os.path.join(workspace["d3"], "metrics_A_prime_candidate_00.json"), "w"))
    calls = []
    _, df = run_module5_on_two_designs(workspace, monkeypatch, "once", calls)
    assert calls.count(("A_wt", "B_wt")) == 1                       # one ceiling fold for two designs (was one each)
    assert calls.count(("A_prime", "B_wt")) == 2                    # regime differs -> rupture is recomputed, not reused
    assert (df["iptm_ceiling"] == 0.85).all() and (df["iptm_rupture"] == 0.25).all()


def test_module5_early_exit_spends_neither_folds_nor_energy_on_hopeless_designs(workspace, monkeypatch):
    with_energy(workspace)
    set_cfg(workspace, thresholds={"early_exit": {"rescue_iptm_floor": 0.9}})           # fake rescue iPTM is 0.8
    energy_calls = []
    calls = []
    ws = workspace
    m5_holder = {}

    def fake_energy(jobs, out_json, config):
        energy_calls.append(jobs)
        return {}

    real_loader = load_module

    def load_and_patch(name):
        mod = real_loader(name)
        if name == "05_eval_final.py":
            mod.score_interfaces = fake_energy
            m5_holder["m"] = mod
        return mod

    monkeypatch.setitem(globals(), "load_module", load_and_patch)
    _, df = run_module5_on_two_designs(ws, monkeypatch, "early", calls)
    assert ("A_wt", "B_prime") not in calls and ("A_prime", "B_wt") not in calls        # no cross / rupture fold
    assert energy_calls == []                                                          # nothing to score, no call at all
    assert (df["status"] == "early_exit_low_rescue").all() and not df["passes"].any()
