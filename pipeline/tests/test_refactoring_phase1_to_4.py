"""Targeted unit tests validating Phase 1 to Phase 4 refactoring:
1. Exact parity between residue_distance_matrix and min_residue_distance.
2. strip_nan_lines preserves HEADER/TITLE containing 'nan' while filtering ATOM/HETATM with NaN coordinates.
3. calculate_ca_rmsd returns float('nan') on missing, empty, or unalignable structures.
4. Consistency of ortho_count and status badges with the 'passes' column in export_styled_reports.py.
"""
import math
import os
import numpy as np
import pandas as pd
import pytest
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue

from conftest import make_chain, write_pdb
from design_utils import min_residue_distance, residue_distance_matrix, strip_nan_lines
from export_styled_reports import generate_interactive_html, generate_styled_excel
from folding_engine import calculate_ca_rmsd


def _make_residue(name, res_id, atoms_spec):
    """Helper to create a BioPython Residue with specified atom names, coords, and elements."""
    res = Residue((' ', res_id, ' '), name, ' ')
    for atom_name, coord, element in atoms_spec:
        atom = Atom(atom_name, coord, 20.0, 1.0, ' ', atom_name, res_id, element=element)
        res.add(atom)
    return res


def test_residue_distance_matrix_parity_with_min_residue_distance():
    """Verify that vectorised residue_distance_matrix yields identical results to min_residue_distance."""
    # Build list of residues for chain A
    res_A = [
        _make_residue("ALA", 1, [
            ("N", np.array([0.0, 0.0, 0.0], dtype=np.float32), "N"),
            ("CA", np.array([1.5, 0.0, 0.0], dtype=np.float32), "C"),
            ("CB", np.array([2.0, 1.2, 0.0], dtype=np.float32), "C"),
            ("H", np.array([-0.5, 0.0, 0.0], dtype=np.float32), "H"),  # Hydrogen to ignore
        ]),
        _make_residue("GLY", 2, [
            ("N", np.array([3.0, 0.0, 0.0], dtype=np.float32), "N"),
            ("CA", np.array([4.5, 0.0, 0.0], dtype=np.float32), "C"),
            ("HA1", np.array([4.5, 1.0, 0.0], dtype=np.float32), "H"),  # Hydrogen to ignore
        ]),
        _make_residue("VAL", 3, [
            ("CA", np.array([7.0, 0.0, 0.0], dtype=np.float32), "C"),
            ("CB", np.array([8.0, 1.0, 0.0], dtype=np.float32), "C"),
            ("CG1", np.array([8.5, 2.0, 0.5], dtype=np.float32), "C"),
        ]),
    ]

    # Build list of residues for chain B
    res_B = [
        _make_residue("LEU", 10, [
            ("CA", np.array([2.0, 4.0, 0.0], dtype=np.float32), "C"),
            ("CB", np.array([2.0, 5.0, 0.0], dtype=np.float32), "C"),
            ("H", np.array([2.0, 3.5, 0.0], dtype=np.float32), "H"),
        ]),
        _make_residue("ILE", 11, [
            ("CA", np.array([6.0, 3.0, 0.0], dtype=np.float32), "C"),
            ("CD1", np.array([6.5, 4.0, 1.0], dtype=np.float32), "C"),
        ]),
    ]

    # Compute vectorized matrix
    dist_mat = residue_distance_matrix(res_A, res_B)
    assert dist_mat.shape == (3, 2)

    # Compare pairwise with min_residue_distance
    for i, rA in enumerate(res_A):
        for j, rB in enumerate(res_B):
            expected = min_residue_distance(rA, rB)
            actual = dist_mat[i, j]
            assert np.isclose(actual, expected, atol=1e-5), (
                f"Mismatch at ({i}, {j}): actual={actual}, expected={expected}"
            )


def test_strip_nan_lines_preserves_nan_headers_while_filtering_atom_nans(tmp_path):
    """Verify strip_nan_lines removes only ATOM/HETATM lines with coordinates containing NaN."""
    pdb_file = tmp_path / "test_nan.pdb"
    content = (
        "HEADER    NANOBODY COMPLEX                                              \n"
        "TITLE     MAINTENANCE AND RESCUE STUDY                                  \n"
        "REMARK    A non-atom remark mentioning nano-technology                  \n"
        "ATOM      1  CA  ALA A   1      10.000  20.000  30.000  1.00 20.00           C\n"
        "ATOM      2  CA  ALA A   2         nan     nan     nan  1.00 20.00           C\n"
        "ATOM      3  CA  ALA A   3      12.000     NaN  32.000  1.00 20.00           C\n"
        "HETATM    4  ZN   ZN A   4         nan  15.000  25.000  1.00 20.00          ZN\n"
        "ATOM      5  CA  ALA A   5      15.000  25.000  35.000  1.00 20.00           C\n"
        "END                                                                     \n"
    )
    pdb_file.write_text(content, encoding="utf-8")

    strip_nan_lines(str(pdb_file))

    cleaned = pdb_file.read_text(encoding="utf-8").splitlines()
    # Headers must be preserved
    assert any("HEADER    NANOBODY COMPLEX" in line for line in cleaned)
    assert any("TITLE     MAINTENANCE AND RESCUE" in line for line in cleaned)
    assert any("REMARK    A non-atom remark" in line for line in cleaned)
    # Valid ATOM lines must be preserved
    assert any("ATOM      1  CA  ALA A   1" in line for line in cleaned)
    assert any("ATOM      5  CA  ALA A   5" in line for line in cleaned)
    assert any("END" in line for line in cleaned)
    # Corrupt ATOM/HETATM lines must be stripped
    assert not any("ATOM      2" in line for line in cleaned)
    assert not any("ATOM      3" in line for line in cleaned)
    assert not any("HETATM    4" in line for line in cleaned)


def test_calculate_ca_rmsd_returns_nan_on_missing_or_empty_structures(tmp_path):
    """Verify calculate_ca_rmsd returns float('nan') on absent or non-alignable files."""
    valid_pdb = tmp_path / "valid.pdb"
    ch = make_chain("A", [1, 2, 3])
    write_pdb(str(valid_pdb), [ch])

    nonexistent_pdb = str(tmp_path / "does_not_exist.pdb")
    empty_dir = tmp_path / "empty_dir"
    empty_dir.mkdir()

    # Missing reference file
    rmsd_no_ref = calculate_ca_rmsd(nonexistent_pdb, str(valid_pdb))
    assert math.isnan(rmsd_no_ref), f"Expected NaN for missing reference, got {rmsd_no_ref}"

    # Missing predicted file
    rmsd_no_pred = calculate_ca_rmsd(str(valid_pdb), nonexistent_pdb)
    assert math.isnan(rmsd_no_pred), f"Expected NaN for missing prediction, got {rmsd_no_pred}"

    # Empty prediction directory
    rmsd_empty_dir = calculate_ca_rmsd(str(valid_pdb), str(empty_dir))
    assert math.isnan(rmsd_empty_dir), f"Expected NaN for empty directory, got {rmsd_empty_dir}"

    # Non-overlapping residue IDs with different chain length (cannot pair by id or order)
    disjoint_pdb = tmp_path / "disjoint.pdb"
    ch_disjoint = make_chain("A", ["ALA", "GLY", "VAL", "LEU"], first_id=101)
    write_pdb(str(disjoint_pdb), [ch_disjoint])
    rmsd_no_overlap = calculate_ca_rmsd(str(valid_pdb), str(disjoint_pdb))
    assert math.isnan(rmsd_no_overlap), f"Expected NaN for unalignable residues, got {rmsd_no_overlap}"

    # Too few residues (< 3 atoms)
    short_pdb = tmp_path / "short.pdb"
    ch_short = make_chain("A", ["ALA", "GLY"], first_id=1)
    write_pdb(str(short_pdb), [ch_short])
    rmsd_short = calculate_ca_rmsd(str(short_pdb), str(short_pdb))
    assert math.isnan(rmsd_short), f"Expected NaN for < 3 residues, got {rmsd_short}"

    # Non-existent subset_res_ids (< 3 atoms)
    rmsd_no_subset = calculate_ca_rmsd(str(valid_pdb), str(valid_pdb), subset_res_ids=[999])
    assert math.isnan(rmsd_no_subset), f"Expected NaN for non-existent subset, got {rmsd_no_subset}"


def test_export_styled_reports_adheres_to_passes_verdict(tmp_path):
    """Verify export_styled_reports uses 'passes' column for ortho_count and badge determination."""
    csv_file = tmp_path / "orthogonality_scores.csv"
    xlsx_file = tmp_path / "orthogonality_scores.xlsx"
    html_file = tmp_path / "orthogonality_scores.html"

    # 4 Candidates:
    # 1. passes=True, f_ortho=0.35 -> Champion Orthogonal
    # 2. passes=True, f_ortho=0.15 -> Orthogonal
    # 3. passes=False, f_ortho=0.25 -> Non Orthogonal (fails rupture or negative despite positive f_ortho)
    # 4. passes=False, f_ortho=0.05, status='early_exit_off_target', dockq=0.10 -> Hors-cible
    df = pd.DataFrame([
        {
            "design_id": "pair_champ",
            "parent_a_motif": "motif_1",
            "iptm_rescue": 0.85,
            "iptm_negative": 0.20,
            "iptm_rupture": 0.25,
            "f_ortho": 0.35,
            "dockq": 0.75,
            "dockq_quality": "High",
            "fnat": 0.80,
            "plddt_rescue": 88.0,
            "pass_dockq": True,
            "passes": True,
            "status": "ok"
        },
        {
            "design_id": "pair_ortho",
            "parent_a_motif": "motif_2",
            "iptm_rescue": 0.70,
            "iptm_negative": 0.30,
            "iptm_rupture": 0.32,
            "f_ortho": 0.15,
            "dockq": 0.50,
            "dockq_quality": "Medium",
            "fnat": 0.60,
            "plddt_rescue": 82.0,
            "pass_dockq": True,
            "passes": True,
            "status": "ok"
        },
        {
            "design_id": "pair_failed_passes",
            "parent_a_motif": "motif_3",
            "iptm_rescue": 0.60,
            "iptm_negative": 0.55,  # Too high negative binding
            "iptm_rupture": 0.50,
            "f_ortho": 0.25,        # positive margin, but passes is False!
            "dockq": 0.60,
            "dockq_quality": "Medium",
            "fnat": 0.65,
            "plddt_rescue": 79.0,
            "pass_dockq": True,
            "passes": False,
            "status": "ok"
        },
        {
            "design_id": "pair_off_target",
            "parent_a_motif": "motif_4",
            "iptm_rescue": 0.40,
            "iptm_negative": 0.30,
            "iptm_rupture": 0.30,
            "f_ortho": 0.05,
            "dockq": 0.10,
            "dockq_quality": "Incorrect",
            "fnat": 0.05,
            "plddt_rescue": 70.0,
            "pass_dockq": False,
            "passes": False,
            "status": "early_exit_off_target"
        },
    ])
    df.to_csv(csv_file, index=False)

    # Run Excel and HTML export
    generate_styled_excel(str(csv_file), str(xlsx_file))
    generate_interactive_html(str(csv_file), str(html_file))

    assert xlsx_file.exists()
    assert html_file.exists()

    html_content = html_file.read_text(encoding="utf-8")

    # 1. ortho_count must be exactly 2 (out of 4)
    assert "{ortho_count}" not in html_content  # properly formatted f-string
    assert ">2 <small style=\"font-size: 14px; color: var(--text-muted);\">" in html_content

    # 2. Status badges must match verdict
    assert '<span class="badge badge-champion">Champion Orthogonal</span>' in html_content
    assert '<span class="badge badge-success">Orthogonal</span>' in html_content
    assert '<span class="badge badge-negative">Non Orthogonal</span>' in html_content
    assert 'Hors-cible' in html_content
