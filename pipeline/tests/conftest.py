import importlib.util
import math
import os
import sys

import numpy as np
import pytest
from Bio.PDB import PDBIO, Atom, Chain, Model, Residue, Structure

ROOT = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(ROOT, "src"))


def load_module(filename):
    """Import a src/NN_name.py module (file names start with a digit)."""
    spec = importlib.util.spec_from_file_location(filename[:-3].replace("-", "_"), os.path.join(ROOT, "src", filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def make_chain(chain_id, resnames, first_id=1, origin=(0.0, 0.0, 0.0), radius=2.3):
    """Idealised helix with N/CA/C/O/CB atoms so contacts/RMSD code has real geometry."""
    chain = Chain.Chain(chain_id)
    for i, name in enumerate(resnames):
        ang = math.radians(100 * i)
        ca = np.array([origin[0] + radius * math.cos(ang), origin[1] + radius * math.sin(ang), origin[2] + 1.5 * i])
        res = Residue.Residue((' ', first_id + i, ' '), name, ' ')
        for j, (an, off) in enumerate([("N", (-1.2, 0.3, -0.5)), ("CA", (0, 0, 0)), ("C", (1.2, 0.4, 0.5)),
                                       ("O", (1.5, 1.4, 0.6)), ("CB", (0.2, -1.2, 0.3))]):
            res.add(Atom.Atom(an, ca + np.array(off), 1.0, 1.0, ' ', f" {an:<3}", j + 1, an[0]))
        chain.add(res)
    return chain


def write_pdb(path, chains):
    st = Structure.Structure("s")
    model = Model.Model(0)
    st.add(model)
    for c in chains:
        model.add(c)
    io = PDBIO()
    io.set_structure(st)
    io.save(str(path))
    return str(path)


@pytest.fixture
def helix_pair(tmp_path):
    seq_a = ["ALA", "LYS", "GLU", "LEU", "ARG", "ASP", "VAL", "SER", "ILE", "THR"] * 2
    seq_b = ["GLU", "ASP", "LYS", "ILE", "VAL", "ARG", "LEU", "GLN", "ASN", "SER"] * 2
    a = make_chain("A", seq_a)
    b = make_chain("B", seq_b, origin=(5.0, 0.0, 0.0))
    return write_pdb(tmp_path / "wt.pdb", [a, b]), seq_a, seq_b
