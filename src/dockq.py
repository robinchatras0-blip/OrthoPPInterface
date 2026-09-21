import glob
import math
import os
import sys

import numpy as np
from Bio.PDB import PDBParser, Superimposer

sys.path.append(os.path.dirname(__file__))
from design_utils import match_residues, min_residue_distance, std_residues

_FAILED = {'dockq': 0.0, 'fnat': 0.0, 'irms': 99.9, 'lrms': 99.9, 'quality': 'Incorrect'}


def calculate_dockq(native_pdb, model_pdb_or_dir, chain_A='A', chain_B='B', contact_threshold=5.0, interface_threshold=10.0):
    """DockQ (Basu & Wallner) of a model against a native complex, on C-alpha for iRMS / lRMS.

    Residues are paired between native and model by residue id when the numbering overlaps, otherwise by
    order (RF3 outputs are renumbered from 1). Returns dockq, fnat, irms, lrms and the CAPRI quality class.
    """
    model_pdb = model_pdb_or_dir
    if os.path.isdir(model_pdb_or_dir):
        ranked = (sorted(glob.glob(os.path.join(model_pdb_or_dir, '*unrelaxed_rank_001*.pdb')))
                  + sorted(glob.glob(os.path.join(model_pdb_or_dir, '*rank_001*.pdb')))
                  + sorted(glob.glob(os.path.join(model_pdb_or_dir, '*.pdb'))))
        if not ranked:
            return dict(_FAILED)
        model_pdb = ranked[0]
    if not os.path.exists(native_pdb) or not os.path.exists(model_pdb):
        return dict(_FAILED)

    parser = PDBParser(QUIET=True)
    m_nat = parser.get_structure('native', native_pdb)[0]
    m_mod = parser.get_structure('model', model_pdb)[0]
    if any(c not in m for m in (m_nat, m_mod) for c in (chain_A, chain_B)):
        return dict(_FAILED)

    nat_A, nat_B = std_residues(m_nat[chain_A]), std_residues(m_nat[chain_B])
    # model residues keyed by the NATIVE residue id
    mod_A = {n.id[1]: m for n, m in match_residues(m_nat[chain_A], m_mod[chain_A])}
    mod_B = {n.id[1]: m for n, m in match_residues(m_nat[chain_B], m_mod[chain_B])}

    # 1. Native contacts and interface residues
    native_contacts, iface_A, iface_B = set(), set(), set()
    for rA in nat_A:
        for rB in nat_B:
            d = min_residue_distance(rA, rB)
            if d <= contact_threshold:
                native_contacts.add((rA.id[1], rB.id[1]))
            if d <= interface_threshold:
                iface_A.add(rA.id[1])
                iface_B.add(rB.id[1])

    # 2. Fraction of native contacts recovered in the model
    if native_contacts:
        recovered = sum(1 for a, b in native_contacts
                        if a in mod_A and b in mod_B and min_residue_distance(mod_A[a], mod_B[b]) <= contact_threshold)
        fnat = recovered / len(native_contacts)
    else:
        fnat = 1.0

    # 3. Interface RMSD (C-alpha of interface residues of both chains)
    nat_ca, mod_ca = [], []
    for chain_id, ids, mod in ((chain_A, sorted(iface_A), mod_A), (chain_B, sorted(iface_B), mod_B)):
        for i in ids:
            r_nat, r_mod = m_nat[chain_id][i], mod.get(i)
            if r_mod is not None and 'CA' in r_nat and 'CA' in r_mod:
                nat_ca.append(r_nat['CA'])
                mod_ca.append(r_mod['CA'])
    if len(nat_ca) >= 3:
        sup = Superimposer()
        sup.set_atoms(nat_ca, mod_ca)
        irms = float(sup.rms)
    else:
        irms = 99.9

    # 4. Ligand RMSD: fit model chain A on native chain A, RMSD of chain B C-alpha
    ids_A = sorted(i for i, m in mod_A.items() if 'CA' in m and 'CA' in m_nat[chain_A][i])
    ids_B = sorted(i for i, m in mod_B.items() if 'CA' in m and 'CA' in m_nat[chain_B][i])
    if len(ids_A) >= 3 and ids_B:
        sup = Superimposer()
        sup.set_atoms([m_nat[chain_A][i]['CA'] for i in ids_A], [mod_A[i]['CA'] for i in ids_A])
        sup.apply(list(m_mod.get_atoms()))
        sq = [np.sum((m_nat[chain_B][i]['CA'].coord - mod_B[i]['CA'].coord) ** 2) for i in ids_B]
        lrms = float(math.sqrt(np.mean(sq)))
    else:
        lrms = 99.9

    # 5. DockQ
    scale_irms = 1.0 / (1.0 + (irms / 1.5) ** 2) if irms < 90.0 else 0.0
    scale_lrms = 1.0 / (1.0 + (lrms / 8.5) ** 2) if lrms < 90.0 else 0.0
    dockq_val = max(0.0, min(1.0, (fnat + scale_irms + scale_lrms) / 3.0))
    quality = 'High' if dockq_val >= 0.80 else 'Medium' if dockq_val >= 0.49 else 'Acceptable' if dockq_val >= 0.23 else 'Incorrect'
    return {'dockq': round(dockq_val, 4), 'fnat': round(fnat, 4), 'irms': round(irms, 3),
            'lrms': round(lrms, 3), 'quality': quality}
