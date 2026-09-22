import json
import math
import os
import re

import pandas as pd
import pytest

import energy_engine
import interface_energy
from conftest import make_chain, write_pdb
from scoring import coherence_report, combined_f_ortho, energy_fields, f_iptm_rel, iptm_margin

NAN = float("nan")


# ------------------------------------------------------------------ scores
def test_energy_fields_are_normalised_to_the_native_complex():
    e = energy_fields(dG_rescue=-45.0, dG_negative=-30.0, dG_rupture=-36.0, dG_native=-60.0)
    assert e["b_rescue"] == pytest.approx(0.75) and e["b_negative"] == pytest.approx(0.5) and e["b_rupture"] == pytest.approx(0.6)
    assert e["f_energy"] == pytest.approx(0.75 - 0.6)                       # margin over the BEST wild-type cross pair
    assert e["gap_negative"] == pytest.approx(-15.0) and e["gap_rupture"] == pytest.approx(-9.0)


def test_energy_margin_is_negative_when_a_wild_type_cross_pair_binds_better():
    assert energy_fields(-45.0, -55.0, -36.0, -60.0)["f_energy"] < 0


def test_energy_fields_are_nan_without_a_binding_native_reference():
    for native in (NAN, 0.0, 12.0):
        assert math.isnan(energy_fields(-45.0, -30.0, -36.0, native)["f_energy"])
    assert math.isnan(energy_fields(NAN, -30.0, -36.0, -60.0)["f_energy"])


def test_iptm_margins_and_combination():
    assert iptm_margin(0.8, 0.25, 0.2) == pytest.approx(0.55)
    assert f_iptm_rel(0.8, 0.25, 0.2, 0.55) == pytest.approx(1.0)
    assert math.isnan(f_iptm_rel(0.8, 0.25, 0.2, 0.0)) and math.isnan(f_iptm_rel(0.8, NAN, 0.2, 0.55))
    assert combined_f_ortho(0.6, 0.2, 0.5, energy_enabled=True) == pytest.approx(0.4)
    assert combined_f_ortho(0.6, 0.2, 0.75, energy_enabled=True) == pytest.approx(0.5)
    assert combined_f_ortho(0.6, 0.2, 0.5, energy_enabled=False) == 0.6           # iPTM only without the energy stage
    assert math.isnan(combined_f_ortho(0.6, NAN, 0.5, energy_enabled=True))       # never a half-score


def test_coherence_report_detects_agreement_and_disagreement():
    df = pd.DataFrame({"design_id": list("abcde"),
                       "iptm_rescue": [0.6, 0.5, 0.4, 0.3, 0.2], "b_rescue": [0.9, 0.8, 0.6, 0.5, 0.3],
                       "iptm_negative": [0.2, 0.3, 0.4, 0.5, 0.6], "b_negative": [0.3, 0.4, 0.5, 0.6, 0.7],
                       "iptm_rupture": [0.2, 0.2, 0.2, 0.2, 0.2], "b_rupture": [0.5, 0.4, 0.6, 0.5, 0.3],
                       "f_iptm_rel": [0.9, 0.5, 0.1, -0.2, -0.6], "f_energy": [0.8, 0.5, 0.05, -0.3, -0.5]})
    rep = coherence_report(df)
    assert rep["spearman_rescue_iptm_vs_binding"] == pytest.approx(1.0)
    assert rep["spearman_negative_iptm_vs_binding"] == pytest.approx(1.0)
    assert rep["margin_sign_agreement"] == 1.0 and rep["disagreements"] == []
    flipped = df.assign(f_energy=df["f_energy"].iloc[::-1].values)
    assert coherence_report(flipped)["spearman_margins"] == pytest.approx(-1.0)
    assert "a" in coherence_report(flipped)["disagreements"]
    assert coherence_report(df.head(3))["spearman_margins"] is None               # too few designs to correlate


# ------------------------------------------------------------------ energy cache
def test_pending_jobs_reuse_only_matching_successful_results(tmp_path):
    pdb = tmp_path / "c.pdb"
    pdb.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n")
    jobs = [{"name": "c", "pdb": str(pdb)}]
    key = interface_energy.job_key(str(pdb), 1, 1, "A_B")
    assert interface_energy.pending_jobs(jobs, {"c": {"key": key, "dG": -10.0}}, 1, 1) == []
    assert len(interface_energy.pending_jobs(jobs, {}, 1, 1)) == 1                                   # never scored
    assert len(interface_energy.pending_jobs(jobs, {"c": {"key": key, "dG": -10.0}}, 3, 1)) == 1      # other settings
    assert len(interface_energy.pending_jobs(jobs, {"c": {"key": key, "error": "boom"}}, 1, 1)) == 1  # failed before
    pdb.write_text("ATOM      1  CA  GLY A   1       1.000   0.000   0.000  1.00  0.00           C\n")
    assert len(interface_energy.pending_jobs(jobs, {"c": {"key": key, "dG": -10.0}}, 1, 1)) == 1      # file changed


# ------------------------------------------------------------------ bridge
class _Done:
    returncode, stdout, stderr = 0, "", ""


def _fake_run(payload, calls):
    def run(cmd, **kw):
        calls.append(cmd)
        out = cmd[cmd.index("--out") + 1]
        if not os.path.isdir(os.path.dirname(out)):                     # a WSL path seen from Windows
            out = re.sub(r"^/mnt/([a-z])/", lambda m: m.group(1).upper() + ":/", out)
        with open(out, "w") as f:
            json.dump(payload, f)
        return _Done()
    return run


def test_bridge_calls_pyrosetta_through_wsl_or_directly(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(energy_engine.subprocess, "run", _fake_run({"n": {"dG": -50.0}}, calls))
    base = {"energy": {"enabled": True, "python": "/env/bin/python", "n_procs": 4, "relax_repeats": 2, "seed": 7}}
    res = energy_engine.score_interfaces([{"name": "n", "pdb": str(tmp_path / "a.pdb")}], str(tmp_path / "e.json"),
                                         {**base, "folding": {"use_wsl": True, "wsl_distro": "Ubuntu"}})
    assert res["n"]["dG"] == -50.0
    assert calls[0][:5] == ["wsl", "-d", "Ubuntu", "--", "/env/bin/python"]
    assert calls[0][calls[0].index("--procs") + 1] == "4" and calls[0][calls[0].index("--relax_repeats") + 1] == "2"
    assert calls[0][calls[0].index("--mode") + 1] == "score"
    energy_engine.score_interfaces([{"name": "n", "pdb": str(tmp_path / "a.pdb")}], str(tmp_path / "e2.json"),
                                   {**base, "folding": {"use_wsl": False}})
    assert calls[1][0] == "/env/bin/python" and "wsl" not in calls[1]


def test_packing_failure_is_never_silent(tmp_path, monkeypatch):
    monkeypatch.setattr(energy_engine.subprocess, "run", _fake_run({"x": {"ok": False, "error": "boom"}}, []))
    with pytest.raises(RuntimeError, match="side-chain packing failed"):
        energy_engine.pack_structures([{"name": "x", "pdb": str(tmp_path / "x.pdb")}],
                                      {"energy": {"enabled": True}, "folding": {"use_wsl": False}}, str(tmp_path / "s.json"))
    monkeypatch.setattr(energy_engine.subprocess, "run", _fake_run({"x": {"ok": True}}, []))
    assert energy_engine.pack_structures([{"name": "x", "pdb": str(tmp_path / "x.pdb")}],
                                         {"energy": {"enabled": True}, "folding": {"use_wsl": False}}, str(tmp_path / "s.json"))["x"]["ok"]


def test_energy_is_off_unless_configured():
    assert not energy_engine.enabled({}) and not energy_engine.enabled({"energy": {"enabled": False}})
    assert energy_engine.enabled({"energy": {"enabled": True}})


# ------------------------------------------------------------------ real PyRosetta (only where it is installed)
def test_packing_rebuilds_missing_side_chains(tmp_path, monkeypatch):
    pytest.importorskip("pyrosetta")
    names = ["LEU", "LYS", "GLU", "VAL", "ILE", "SER", "PHE", "ARG"] * 3
    pdb = write_pdb(tmp_path / "bb_only.pdb", [make_chain("A", names), make_chain("B", names, origin=(9.0, 0.0, 0.0))])
    interface_energy._init_worker(1)
    name, res = interface_energy.pack_one({"name": "x", "pdb": pdb, "out": str(tmp_path / "packed.pdb"), "fixed_chains": ["A"]})
    assert res["ok"], res
    from Bio.PDB import PDBParser
    st = PDBParser(QUIET=True).get_structure("p", str(tmp_path / "packed.pdb"))[0]
    heavy = lambda r: sum(1 for a in r if a.element != "H")  # noqa: E731
    assert all(heavy(r) > 5 for r in st["B"] if r.get_resname() in ("LYS", "ARG", "PHE"))        # side chains exist now
    assert not any(a.element == "H" for a in st.get_atoms())                                     # heavy atoms only
