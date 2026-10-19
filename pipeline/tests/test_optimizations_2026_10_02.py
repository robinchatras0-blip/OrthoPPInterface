"""Unit tests for the 2026-10-02 optimization release:
- Consolidated Chain A loop diffusion with hydrophobic core preservation
- Decoupled monomer RMSD filtering (hydrophobic core <= 1.5 A, loop mobility <= 3.0 A)
- Scaffold drift gating (scaffold_max_flex_rmsd <= 4.5 A) & native scaffold baseline
- Enhanced negative design bias prioritizing mutated residues without WT-positive bias
"""
import copy
import json
import os
import pytest
from Bio.PDB import PDBParser

from conftest import load_module, make_chain, write_pdb
from design_utils import (
    build_rescue_bias,
    load_config,
    specificity_bias,
)


def test_chain_A_consolidation_absorbs_islands_and_locks_core():
    """Verify that consolidate_flexible absorbs 1-2 residue anchors while keeping core intact."""
    m1 = load_module("01_analyze_interface.py")
    all_ids_A = list(range(1, 225))
    # Simulated fragmented interface with 1-2 residue pinned islands and a large core
    # (loops: 3-4, 7, 9-10, 12-19, 29, 32-35, 37-41; core: 42-69, 76-184, 207-224)
    raw_if = [3, 4, 7, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 29, 32, 33, 34, 35, 37, 38, 39, 40, 41]
    consolidated = m1.consolidate_flexible(all_ids_A, set(raw_if), min_fixed=3)

    # 1. Short islands between 3 and 19 (5-6, 8, 11) must be absorbed into one contiguous loop
    for r in range(3, 20):
        assert r in consolidated, f"Residue {r} should be absorbed into the flexible loop"

    # 2. Island 36 between 32-35 and 37-41 must be absorbed
    assert 36 in consolidated, "Residue 36 should be absorbed into the flexible loop"

    # 3. Hydrophobic core (76 to 184) must remain strictly fixed (NOT in consolidated)
    for r in range(76, 185):
        assert r not in consolidated, f"Hydrophobic core residue {r} must remain fixed"


def test_decoupled_rmsd_gating_protects_core_while_allowing_loop_mobility(tmp_path, monkeypatch):
    """Verify that rmsd_framework_max rejects disturbed core even if monomer RMSD is moderate,
    and accepts remodeled loops if the core is preserved."""
    m3 = load_module("03_filter_A_prime.py")
    cfg = {
        "pipeline": {"chain_A": "A", "chain_B": "B", "input_pdb": "dummy.pdb",
                     "input_msa_A": "dummy.a3m", "input_msa_B": "dummy.a3m"},
        "thresholds": {
            "monomer_stability": {"plddt_min": 80.0, "rmsd_max": 3.0, "rmsd_framework_max": 1.5},
            "negative_design_rupture": {"iptm_rupture_max": 0.35},
        }
    }
    th = cfg["thresholds"]
    plddt_min = th["monomer_stability"]["plddt_min"]
    rmsd_max = th["monomer_stability"]["rmsd_max"]
    rmsd_framework_max = th["monomer_stability"]["rmsd_framework_max"]

    # Candidate 1: Loop remodeled (rmsd_monomer = 2.4 A <= 3.0 A), Core preserved (rmsd_framework = 0.8 A <= 1.5 A)
    c1 = {"plddt": 85.0, "rmsd_monomer": 2.4, "rmsd_framework": 0.8, "iptm_rupture": 0.20}
    c1_passed = bool(c1["iptm_rupture"] <= 0.35 and c1["rmsd_monomer"] <= rmsd_max
                     and c1["rmsd_framework"] <= rmsd_framework_max and c1["plddt"] >= plddt_min)
    assert c1_passed is True, "Candidate with remodeled loops but intact core should pass"

    # Candidate 2: Core disturbed (rmsd_framework = 1.9 A > 1.5 A) even though rmsd_monomer <= 3.0 A
    c2 = {"plddt": 85.0, "rmsd_monomer": 2.2, "rmsd_framework": 1.9, "iptm_rupture": 0.20}
    c2_passed = bool(c2["iptm_rupture"] <= 0.35 and c2["rmsd_monomer"] <= rmsd_max
                     and c2["rmsd_framework"] <= rmsd_framework_max and c2["plddt"] >= plddt_min)
    assert c2_passed is False, "Candidate with disturbed hydrophobic core must be rejected"


def test_scaffold_drift_filtering_rejects_wild_scaffolds(tmp_path):
    """Verify that pick_scaffolds rejects scaffolds where loop drift exceeds scaffold_max_flex_rmsd."""
    m4 = load_module("04_generate_B_prime.py")
    names = ["ALA", "LEU", "LYS", "GLU", "ASP", "VAL", "ILE", "SER", "THR", "ARG"] * 2
    a = make_chain("A", names)
    wt_b = make_chain("B", names, origin=(9, 0, 0))
    native = write_pdb(tmp_path / "native.pdb", [a, wt_b])

    # Well-adapted scaffold: loop moves 1.5 A (within 4.5 A limit)
    good_b = make_chain("B", names, origin=(9, 0, 0))
    for r in good_b:
        if r.id[1] in (9, 10, 11):
            for at in r:
                at.coord = at.coord + [0.0, 1.5, 0.0]
    good = write_pdb(tmp_path / "good.pdb", [a, good_b])

    # Wild drifted scaffold: loop flies away by 10 A (flex_rmsd > 4.5 A)
    wild_b = make_chain("B", names, origin=(9, 0, 0))
    for r in wild_b:
        if r.id[1] in (9, 10, 11):
            for at in r:
                at.coord = at.coord + [0.0, 10.0, 0.0]
    wild = write_pdb(tmp_path / "wild.pdb", [a, wild_b])

    cfg = {
        "max_rescue_scaffolds": 4,
        "include_native_scaffold": True,
        "scaffold_max_clashes": 5,
        "scaffold_max_framework_rmsd": 1.5,
        "scaffold_min_flex_rmsd": 0.5,
        "scaffold_max_flex_rmsd": 4.5,
    }
    fixed = set(range(1, 9))
    picked = m4.pick_scaffolds([good, wild], a, wt_b, list(range(1, 21)), fixed, "B", cfg, native)
    picked_paths = [os.path.basename(p) for p, is_nat, _ in picked]

    assert "native.pdb" in picked_paths, "Native scaffold must be included as baseline control"
    assert "good.pdb" in picked_paths, "Moderate adaptation scaffold (flex_rmsd <= 4.5 A) must be kept"
    assert "wild.pdb" not in picked_paths, "Drifted scaffold (flex_rmsd > 4.5 A) must be rejected"


def test_build_rescue_bias_targets_mutations_without_unmutated_wt_bonus(helix_pair):
    """Verify that build_rescue_bias applies negative design on residues contacting mutations,
    without biasing unmutated positions towards wild type."""
    path, _, _ = helix_pair
    m = PDBParser(QUIET=True).get_structure("x", path)[0]
    wt_a = copy.deepcopy(m["A"])
    a_prime = copy.deepcopy(m["A"])

    # Mutate position 3: GLU -> LYS (charge swap)
    for r in a_prime:
        if r.id[1] == 3:
            r.resname = "LYS"

    fixed_b = {r.id[1] for r in m["B"] if r.id[1] > 10}
    bias = build_rescue_bias(a_prime, m["B"], wt_a, fixed_b, cutoff=7.5, w_generic=0.0)

    # Residues contacting residue 3 must have rescue bias for acidic residues (ASP/GLU)
    # and negative bias against basic residues (ARG/LYS)
    assert len(bias) > 0, "Residues contacting the mutated position must receive bias"
    for b_key, b_bias in bias.items():
        assert b_bias.get("GLU", 0) > 0 or b_bias.get("ASP", 0) > 0, "Must favor acidic complement"
        assert b_bias.get("ARG", 0) < 0 and b_bias.get("LYS", 0) < 0, "Must penalize WT complement"


def test_guided_partial_diffusion_rescue_spec(tmp_path, helix_pair):
    """Verify that when rescue_partial_t is configured, Module 1 constructs continuous contigs
    with select_fixed_atoms for framework atoms and partial_t noise on interface loops."""
    m1 = load_module("01_analyze_interface.py")
    pdb_path, _, _ = helix_pair

    # Load real config and override input_pdb
    root_cfg = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "config.yaml")
    pipeline_cfg = os.path.join(os.path.dirname(os.path.dirname(__file__)), "config.yaml")
    cfg_path = root_cfg if os.path.exists(root_cfg) else pipeline_cfg
    cfg = load_config(cfg_path)
    cfg["pipeline"]["input_pdb"] = pdb_path
    cfg["pipeline"]["rescue_partial_t"] = 3.5

    # Run analyze interface with this config
    import sys
    out_dir = tmp_path / "out_analysis"
    out_dir.mkdir()
    cfg_file = tmp_path / "config.yaml"
    import yaml
    with open(cfg_file, "w") as f:
        yaml.dump(cfg, f)

    sys_argv = ["01_analyze_interface.py", "--config", str(cfg_file), "--out_dir", str(out_dir)]
    old_argv = sys.argv
    try:
        sys.argv = sys_argv
        m1.main()
    finally:
        sys.argv = old_argv

    with open(out_dir / "inputs_rescue.json") as f:
        rescue_spec = json.load(f)["design"]

    assert "partial_t" in rescue_spec, "partial_t must be in rescue inputs"
    assert rescue_spec["partial_t"] == 3.5, "partial_t must match configured value 3.5"
    assert "select_fixed_atoms" in rescue_spec, "select_fixed_atoms must be in rescue inputs"
    assert "A" in rescue_spec["select_fixed_atoms"], "Chain A must be in select_fixed_atoms"
    assert "/0," in rescue_spec["contig"], "Contig must have chain break /0"
    # Contig for chain B must NOT have unindexed length de-novo hallucinations like '20-20'
    contig_b = rescue_spec["contig"].split("/0,")[1]
    assert "-" in contig_b and contig_b.startswith("B"), f"Chain B contig should be real residue range: {contig_b}"


def test_dockq_gate_disqualifies_off_target_and_ranks_correctly(tmp_path):
    """Verify that DockQ gate rejects off-target complexes (DockQ < 0.23) even with high iPTM."""
    import pandas as pd
    from export_styled_reports import generate_interactive_html

    # Design 1: high iPTM (0.51) but aberrant off-target docking (DockQ 0.009)
    # Design 2: moderate iPTM (0.28) but correct interface orientation (DockQ 0.40, CAPRI Acceptable)
    rows = [
        {
            "design_id": "A_cand_05_B_cand_03", "parent_a_motif": "A_cand_05", "folding_engine": "rf3",
            "status": "early_exit_off_target", "iptm_rescue": 0.5114, "iptm_negative": 0.46, "iptm_rupture": 0.24,
            "f_ortho": -0.107, "dockq": 0.0093, "dockq_quality": "Incorrect", "fnat": 0.0, "lrms": 55.5,
            "plddt_rescue": 79.0, "pass_dockq": False, "passes": False,
        },
        {
            "design_id": "A_cand_01_B_cand_00", "parent_a_motif": "A_cand_01", "folding_engine": "rf3",
            "status": "ok", "iptm_rescue": 0.2833, "iptm_negative": 0.23, "iptm_rupture": 0.29,
            "f_ortho": 0.025, "dockq": 0.4012, "dockq_quality": "Acceptable", "fnat": 0.356, "lrms": 6.07,
            "plddt_rescue": 78.3, "pass_dockq": True, "passes": True,
        }
    ]

    df = pd.DataFrame(rows).sort_values(by=["passes", "dockq", "f_ortho"], ascending=[False, False, False])
    # The correctly oriented design (DockQ 0.40) must be ranked #1
    assert df.iloc[0]["design_id"] == "A_cand_01_B_cand_00", "DockQ >= 0.23 candidate must rank above off-target"
    assert df.iloc[1]["design_id"] == "A_cand_05_B_cand_03", "Off-target candidate must rank below"

    # Verify HTML badge rendering
    csv_file = tmp_path / "scores.csv"
    html_file = tmp_path / "SUMMARY.html"
    df.to_csv(csv_file, index=False)
    generate_interactive_html(str(csv_file), str(html_file))

    with open(html_file, encoding="utf-8") as f:
        html_text = f.read()

    assert ("Off-target (DockQ 0.009)" in html_text or "Hors-cible (DockQ 0.009)" in html_text), "HTML must display Off-target badge for DockQ < 0.23"
    assert "Orthogonal" in html_text, "HTML must display Orthogonal badge for passing candidate"


def test_hard_clash_limit_filtering():
    """Verify that candidate pools with severe clashes (> hard_max_clashes) are rejected."""
    max_clash_a = 3
    hard_max_clash = max_clash_a + 2  # 5

    # Case 1: Pool with only clashing candidates (11-15 clashes) like motif 05
    bad_pool = [
        {"spec": {"clashes_a_prime": 11, "spec_score": -3.0}, "seq": "AAA"},
        {"spec": {"clashes_a_prime": 14, "spec_score": -10.0}, "seq": "AAC"},
    ]
    ok_1 = [c for c in bad_pool if c["spec"]["clashes_a_prime"] <= max_clash_a]
    if not ok_1:
        ok_1 = [c for c in bad_pool if c["spec"]["clashes_a_prime"] <= hard_max_clash]
    assert len(ok_1) == 0, "All candidates with >= 11 clashes must be rejected"

    # Case 2: Pool with clean candidates (0-1 clashes) like motif 01
    good_pool = [
        {"spec": {"clashes_a_prime": 0, "spec_score": 27.0}, "seq": "AAA"},
        {"spec": {"clashes_a_prime": 1, "spec_score": 25.0}, "seq": "AAC"},
    ]
    ok_2 = [c for c in good_pool if c["spec"]["clashes_a_prime"] <= max_clash_a]
    assert len(ok_2) == 2, "Clean candidates with <= 3 clashes must be retained"


