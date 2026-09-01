import os
import math
import glob
import numpy as np
from Bio.PDB import PDBParser, Superimposer

def get_heavy_atoms(residue):
    return [atom for atom in residue if atom.element != 'H' and not atom.get_name().startswith('H')]

def calculate_dockq(native_pdb, model_pdb_or_dir, chain_A='A', chain_B='B', contact_threshold=5.0, interface_threshold=10.0):
    model_pdb = model_pdb_or_dir
    if os.path.isdir(model_pdb_or_dir):
        ranked = sorted(glob.glob(os.path.join(model_pdb_or_dir, '*unrelaxed_rank_001*.pdb'))) +                  sorted(glob.glob(os.path.join(model_pdb_or_dir, '*rank_001*.pdb')))
        if not ranked:
            ranked = sorted(glob.glob(os.path.join(model_pdb_or_dir, '*.pdb')))
        if ranked:
            model_pdb = ranked[0]
        else:
            return {'dockq': 0.0, 'fnat': 0.0, 'irms': 99.9, 'lrms': 99.9, 'quality': 'Incorrect'}

    if not os.path.exists(native_pdb) or not os.path.exists(model_pdb):
        return {'dockq': 0.0, 'fnat': 0.0, 'irms': 99.9, 'lrms': 99.9, 'quality': 'Incorrect'}

    parser = PDBParser(QUIET=True)
    struct_nat = parser.get_structure('native', native_pdb)
    struct_mod = parser.get_structure('model', model_pdb)

    m_nat = struct_nat[0]
    m_mod = struct_mod[0]

    if chain_A not in m_nat or chain_B not in m_nat or chain_A not in m_mod or chain_B not in m_mod:
        return {'dockq': 0.0, 'fnat': 0.0, 'irms': 99.9, 'lrms': 99.9, 'quality': 'Incorrect'}

    nat_A = [r for r in m_nat[chain_A] if r.id[0] == ' ']
    nat_B = [r for r in m_nat[chain_B] if r.id[0] == ' ']
    mod_A = [r for r in m_mod[chain_A] if r.id[0] == ' ']
    mod_B = [r for r in m_mod[chain_B] if r.id[0] == ' ']

    mod_A_dict = {r.id[1]: r for r in mod_A}
    mod_B_dict = {r.id[1]: r for r in mod_B}

    # 1. Native Contacts and Interface Residues
    native_contacts = set()
    interface_A_ids = set()
    interface_B_ids = set()

    for rA in nat_A:
        atomsA = get_heavy_atoms(rA)
        if not atomsA:
            continue
        for rB in nat_B:
            atomsB = get_heavy_atoms(rB)
            if not atomsB:
                continue
            min_d = float('inf')
            for aA in atomsA:
                for aB in atomsB:
                    d = aA - aB
                    if d < min_d:
                        min_d = d
            if min_d <= contact_threshold:
                native_contacts.add((rA.id[1], rB.id[1]))
            if min_d <= interface_threshold:
                interface_A_ids.add(rA.id[1])
                interface_B_ids.add(rB.id[1])

    # 2. Compute F_nat in Model
    if not native_contacts:
        fnat = 1.0
    else:
        recovered_contacts = 0
        for (idA, idB) in native_contacts:
            if idA in mod_A_dict and idB in mod_B_dict:
                atomsA = get_heavy_atoms(mod_A_dict[idA])
                atomsB = get_heavy_atoms(mod_B_dict[idB])
                min_d = float('inf')
                for aA in atomsA:
                    for aB in atomsB:
                        d = aA - aB
                        if d < min_d:
                            min_d = d
                if min_d <= contact_threshold:
                    recovered_contacts += 1
        fnat = recovered_contacts / len(native_contacts)

    # 3. Compute Interface RMSD (I_RMS)
    nat_iface_ca = []
    mod_iface_ca = []
    for idA in sorted(interface_A_ids):
        r_nat = m_nat[chain_A][idA] if idA in m_nat[chain_A] else None
        r_mod = mod_A_dict.get(idA)
        if r_nat and r_mod and 'CA' in r_nat and 'CA' in r_mod:
            nat_iface_ca.append(r_nat['CA'])
            mod_iface_ca.append(r_mod['CA'])
    for idB in sorted(interface_B_ids):
        r_nat = m_nat[chain_B][idB] if idB in m_nat[chain_B] else None
        r_mod = mod_B_dict.get(idB)
        if r_nat and r_mod and 'CA' in r_nat and 'CA' in r_mod:
            nat_iface_ca.append(r_nat['CA'])
            mod_iface_ca.append(r_mod['CA'])

    if len(nat_iface_ca) >= 3 and len(mod_iface_ca) >= 3:
        sup_i = Superimposer()
        sup_i.set_atoms(nat_iface_ca, mod_iface_ca)
        irms = float(sup_i.rms)
    else:
        irms = 99.9

    # 4. Compute Ligand RMSD (L_RMS)
    common_A_ids = sorted(set([r.id[1] for r in nat_A if 'CA' in r]) & set([r.id[1] for r in mod_A if 'CA' in r]))
    common_B_ids = sorted(set([r.id[1] for r in nat_B if 'CA' in r]) & set([r.id[1] for r in mod_B if 'CA' in r]))

    if len(common_A_ids) >= 3 and len(common_B_ids) >= 1:
        nat_A_ca = [m_nat[chain_A][i]['CA'] for i in common_A_ids]
        mod_A_ca = [mod_A_dict[i]['CA'] for i in common_A_ids]
        sup_r = Superimposer()
        sup_r.set_atoms(nat_A_ca, mod_A_ca)
        sup_r.apply(m_mod.get_atoms())
        diff_sq = []
        for i in common_B_ids:
            nat_ca = m_nat[chain_B][i]['CA']
            mod_ca = mod_B_dict[i]['CA']
            diff = nat_ca.coord - mod_ca.coord
            diff_sq.append(np.sum(diff ** 2))
        lrms = float(math.sqrt(np.mean(diff_sq))) if diff_sq else 99.9
    else:
        lrms = 99.9

    # 5. DockQ Formula
    scale_irms = 1.0 / (1.0 + (irms / 1.5) ** 2) if irms < 90.0 else 0.0
    scale_lrms = 1.0 / (1.0 + (lrms / 8.5) ** 2) if lrms < 90.0 else 0.0
    dockq_val = (fnat + scale_irms + scale_lrms) / 3.0
    dockq_val = max(0.0, min(1.0, float(dockq_val)))

    if dockq_val >= 0.80:
        quality = 'High'
    elif dockq_val >= 0.49:
        quality = 'Medium'
    elif dockq_val >= 0.23:
        quality = 'Acceptable'
    else:
        quality = 'Incorrect'

    return {
        'dockq': round(dockq_val, 4),
        'fnat': round(fnat, 4),
        'irms': round(irms, 3),
        'lrms': round(lrms, 3),
        'quality': quality
    }
