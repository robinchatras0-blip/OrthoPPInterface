import copy

import numpy as np
import pytest
from Bio.PDB import PDBParser

from design_utils import (
    build_rescue_bias, candidate_specificity_score, ca_rmsd_after_fit, diff_positions, hamming,
    interface_stats, ligand_rmsd_after_receptor_fit, match_residues, select_diverse, specificity_bias,
)
from conftest import make_chain


def test_hamming_dict_uses_values_not_keys():
    a = {1: 'A', 2: 'C', 3: 'D'}
    b = {1: 'A', 2: 'W', 3: 'W'}
    assert hamming(a, b) == 2
    assert hamming("ACD", "AWW") == 2


def test_old_selection_bug_is_gone():
    """Regression: zip() over {resid: aa} dicts iterated the KEYS, so every distance was 0 and the
    'farthest-point' pick degenerated to 'first K in list order' (= all from scaffold 0)."""
    keys = list(range(5))
    seqs = [dict(zip(keys, "AAAAA")), dict(zip(keys, "AAAAC")),      # scaffold 0 near-duplicates
            dict(zip(keys, "WWWWW")), dict(zip(keys, "KKKKK"))]      # scaffold 1 / 2, very different
    old = [0]
    while len(old) < 3:                                              # the buggy algorithm, verbatim
        best, best_d = None, -1
        for i in range(len(seqs)):
            if i in old:
                continue
            d = min(sum(1 for a, b in zip(seqs[i], seqs[s]) if a != b) for s in old)
            if d > best_d:
                best, best_d = i, d
        old.append(best)
    assert old == [0, 1, 2]                                          # picked the near-duplicate #1
    assert select_diverse(seqs, 3) != [0, 1, 2]
    assert 1 not in select_diverse(seqs, 3)


def test_group_quota_guarantees_every_scaffold_is_represented():
    seqs = ["AAAA", "AAAC", "AACC", "WWWW", "KKKK"]
    groups = [0, 0, 0, 1, 2]
    scores = [3.0, 2.9, 2.8, 0.1, 0.2]        # scaffold 0 dominates by score
    sel = select_diverse(seqs, 3, scores=scores, groups=groups, min_per_group=1)
    assert {groups[i] for i in sel} == {0, 1, 2}


def test_select_diverse_avoids_duplicates_and_handles_small_pools():
    assert select_diverse(["AA", "AA", "AC"], 2, scores=[1, 2, 0]) in ([1, 2], [0, 2])
    assert select_diverse([], 3) == []
    assert sorted(select_diverse(["AA", "AC"], 5)) == [0, 1]


def test_diff_positions():
    assert diff_positions("ACDE", "ACWE") == [2]
    with pytest.raises(ValueError):
        diff_positions("ACD", "ACDE")


def test_specificity_bias_charge_swap_is_negative_designed():
    b = specificity_bias("GLU", "LYS")          # A: Lys -> Glu
    assert b["ARG"] > 0 and b["LYS"] > 0        # B' should be positive to pair with the new Glu ...
    assert b["ASP"] < 0                         # ... and NOT acidic (which would still fit WT Lys)
    assert specificity_bias("LEU", "LEU") == {}


def test_specificity_bias_knob_hole():
    assert specificity_bias("TRP", "ALA").get("GLY", 0) > 0     # A' knob -> B' hole
    assert specificity_bias("ALA", "TRP").get("PHE", 0) > 0     # A' hole -> B' fills it


def test_match_residues_falls_back_to_order_when_numbering_differs():
    ref = make_chain("A", ["ALA"] * 6, first_id=10)
    mob = make_chain("A", ["ALA"] * 6, first_id=1)              # RF3-style renumbering
    pairs = match_residues(ref, mob)
    assert [(r.id[1], m.id[1]) for r, m in pairs] == [(10 + i, 1 + i) for i in range(6)]


def test_ca_rmsd_is_zero_for_rigid_transform():
    ref = make_chain("A", ["ALA"] * 12)
    mob = copy.deepcopy(ref)
    ang = 0.7
    rot = np.array([[np.cos(ang), -np.sin(ang), 0], [np.sin(ang), np.cos(ang), 0], [0, 0, 1]])
    for at in mob.get_atoms():
        at.coord = at.coord @ rot + np.array([3.0, -2.0, 5.0])
    fit, _, _ = ca_rmsd_after_fit(match_residues(ref, mob))
    assert fit < 1e-4


def test_ligand_rmsd_zero_for_identical_complex(helix_pair):
    path, _, _ = helix_pair
    m = PDBParser(QUIET=True).get_structure("x", path)[0]
    assert ligand_rmsd_after_receptor_fit(match_residues(m["A"], m["A"]), match_residues(m["B"], m["B"])) < 1e-6


def test_interface_stats_and_specificity_proxy(helix_pair):
    path, _, _ = helix_pair
    m = PDBParser(QUIET=True).get_structure("x", path)[0]
    contacts, _ = interface_stats(m["A"], m["B"])
    assert contacts > 0
    far = make_chain("B", ["ALA"] * 5, origin=(80.0, 0, 0))
    assert interface_stats(m["A"], far) == (0, 0)
    s = candidate_specificity_score(m["A"], m["A"], m["B"])   # A' == A_WT -> nothing gained
    assert s["contacts_a_prime"] == s["contacts_a_wt"]


def test_build_rescue_bias_only_touches_mutable_positions(helix_pair):
    path, _, _ = helix_pair
    m = PDBParser(QUIET=True).get_structure("x", path)[0]
    wt_a = copy.deepcopy(m["A"])
    a_prime = copy.deepcopy(m["A"])
    for r in a_prime:
        if r.id[1] == 3:
            r.resname = "LYS"                               # was GLU -> charge swap
    fixed_b = {r.id[1] for r in m["B"] if r.id[1] > 10}
    bias = build_rescue_bias(a_prime, m["B"], wt_a, fixed_b)
    assert bias and all(int(k[1:]) not in fixed_b for k in bias)
