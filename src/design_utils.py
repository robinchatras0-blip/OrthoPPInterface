"""Shared, dependency-light helpers for the OrthoPPInterface design modules.

Everything here is pure Python / Biopython so it can be unit-tested without a GPU.
"""
import math
import os

import numpy as np
from Bio.PDB import NeighborSearch, PDBParser

THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

# Approximate side-chain volumes (A^3, Zamyatnin) used for knob/hole complementarity.
RESIDUE_VOLUME = {
    'GLY': 60.1, 'ALA': 88.6, 'SER': 89.0, 'CYS': 108.5, 'ASP': 111.1, 'PRO': 112.7,
    'ASN': 114.1, 'THR': 116.1, 'GLU': 138.4, 'VAL': 140.0, 'GLN': 143.8, 'HIS': 153.2,
    'MET': 162.9, 'ILE': 166.7, 'LEU': 166.7, 'LYS': 168.6, 'ARG': 173.4, 'PHE': 189.9,
    'TYR': 193.6, 'TRP': 227.8,
}

_POS, _NEG = {'ARG', 'LYS'}, {'ASP', 'GLU'}
_HYD = {'LEU', 'ILE', 'VAL', 'PHE', 'MET', 'TRP', 'TYR', 'ALA'}
_POLAR = {'ASN', 'GLN', 'SER', 'THR', 'HIS'}

# Residues that make a favourable contact with a given partner residue.
COMPLEMENT = {
    **{r: {'ASP', 'GLU'} for r in ('ARG', 'LYS')},
    'HIS': {'ASP', 'GLU', 'ASN', 'GLN'},
    **{r: {'ARG', 'LYS'} for r in ('ASP', 'GLU')},
    **{r: {'LEU', 'ILE', 'VAL', 'PHE', 'MET'} for r in ('PHE', 'TYR', 'TRP', 'LEU', 'ILE', 'VAL', 'MET')},
    'ALA': {'LEU', 'ILE', 'VAL', 'PHE'},
    **{r: {'GLN', 'ASN', 'SER', 'THR', 'ARG', 'GLU'} for r in ('ASN', 'GLN', 'SER', 'THR')},
}
_SMALL = {'ALA', 'GLY', 'SER', 'THR', 'VAL'}
_LARGE = {'TRP', 'PHE', 'TYR', 'MET', 'ARG', 'LEU'}


# --------------------------------------------------------------------------- #
# Sequence helpers
# --------------------------------------------------------------------------- #
def seq_to_str(seq):
    """Normalise a sequence given as str / list / {resid: aa} dict to a comparable str."""
    if isinstance(seq, dict):
        return "".join(seq[k] for k in sorted(seq))
    return "".join(seq)


def hamming(a, b):
    """Hamming distance. Dicts ({resid: aa}) are compared over the union of keys;
    strings position by position (a length difference counts as mismatches)."""
    if isinstance(a, dict) or isinstance(b, dict):
        a = a if isinstance(a, dict) else dict(enumerate(a))
        b = b if isinstance(b, dict) else dict(enumerate(b))
        return sum(1 for k in set(a) | set(b) if a.get(k) != b.get(k))
    return sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))


def diff_positions(seq_wt, seq_design):
    """0-based indices where the designed sequence differs from WT (same length required)."""
    if len(seq_wt) != len(seq_design):
        raise ValueError(f"Length mismatch WT={len(seq_wt)} vs design={len(seq_design)}")
    return [i for i, (a, b) in enumerate(zip(seq_wt, seq_design)) if a != b]


def select_diverse(seqs, k, scores=None, groups=None, min_per_group=0):
    """Score-aware farthest-point selection. Returns the selected indices (ordered).

    - `scores` (higher is better) seeds the selection and breaks distance ties.
    - `groups` (e.g. RFD3 scaffold id per candidate) with `min_per_group` >= 1 guarantees that
      each group is represented (best-scoring member first) before the diversity phase.
    - Exact duplicate sequences are only picked when nothing else is left.
    """
    n = len(seqs)
    if n == 0 or k <= 0:
        return []
    k = min(k, n)
    scores = list(scores) if scores is not None else [0.0] * n
    order = sorted(range(n), key=lambda i: (-scores[i], i))
    selected = []

    if groups is not None and min_per_group > 0:
        by_group = {}
        for i in order:
            by_group.setdefault(groups[i], []).append(i)
        # best-scoring group representative first
        reps = sorted(by_group.values(), key=lambda idxs: (-scores[idxs[0]], idxs[0]))
        for rank in range(min_per_group):
            for idxs in reps:
                if len(selected) >= k:
                    break
                if rank < len(idxs):
                    selected.append(idxs[rank])

    if not selected:
        selected.append(order[0])

    while len(selected) < k:
        best, best_key = None, None
        for i in order:
            if i in selected:
                continue
            d = min(hamming(seqs[i], seqs[s]) for s in selected)
            key = (d > 0, d, scores[i])
            if best_key is None or key > best_key:
                best, best_key = i, key
        if best is None:
            break
        selected.append(best)
    return selected


# --------------------------------------------------------------------------- #
# Structure helpers
# --------------------------------------------------------------------------- #
def std_residues(chain):
    return [r for r in chain if r.id[0] == ' ']


def heavy_atoms(entity):
    return [a for a in entity.get_atoms() if a.element != 'H' and not a.get_name().startswith('H')]


def load_chain(pdb_path, chain_id):
    st = PDBParser(QUIET=True).get_structure("s", pdb_path)
    model = st[0]
    if chain_id in model:
        return model[chain_id]
    return list(model.get_chains())[0]


def match_residues(ref_chain, mob_chain):
    """Pairs residues of two chains: by residue id when they overlap well, otherwise by order
    (only if the chains have identical length). Returns [(ref_res, mob_res), ...]."""
    ref, mob = std_residues(ref_chain), std_residues(mob_chain)
    mob_by_id = {r.id[1]: r for r in mob}
    pairs = [(r, mob_by_id[r.id[1]]) for r in ref if r.id[1] in mob_by_id]
    if ref and len(pairs) >= 0.9 * len(ref):
        return pairs
    if len(ref) == len(mob):
        return list(zip(ref, mob))
    return pairs


def interface_stats(chain_x, chain_y, contact_cut=4.5, clash_cut=2.2):
    """(#residue-residue contacts, #heavy-atom clashes) between two chains."""
    atoms_y = heavy_atoms(chain_y)
    if not atoms_y:
        return 0, 0
    ns = NeighborSearch(atoms_y)
    contacts, clashes = set(), 0
    for a in heavy_atoms(chain_x):
        for b in ns.search(a.coord, contact_cut):
            contacts.add((a.get_parent().id[1], b.get_parent().id[1]))
            if (a - b) < clash_cut:
                clashes += 1
    return len(contacts), clashes


def ca_rmsd_after_fit(pairs, fit_ids=None, eval_ids=None):
    """Fit `mob` onto `ref` on the CA atoms with resid in fit_ids (all if None), and return
    (rmsd_fit, rmsd_eval, rotran). Does not modify the inputs."""
    def cas(sel):
        return [(r['CA'].coord, m['CA'].coord) for r, m in pairs
                if 'CA' in r and 'CA' in m and (sel is None or r.id[1] in sel)]
    fit = cas(fit_ids)
    if len(fit) < 3:
        return None, None, None
    ref_xyz = np.array([f[0] for f in fit])
    mob_xyz = np.array([f[1] for f in fit])
    rc, mc = ref_xyz.mean(0), mob_xyz.mean(0)
    u, _, vt = np.linalg.svd((mob_xyz - mc).T @ (ref_xyz - rc))
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1, 1, d]) @ vt
    def apply(x):
        return (x - mc) @ rot + rc
    rmsd_fit = float(np.sqrt(np.mean(np.sum((apply(mob_xyz) - ref_xyz) ** 2, axis=1))))
    ev = cas(eval_ids)
    if not ev:
        return rmsd_fit, None, (rot, mc, rc)
    e_ref = np.array([f[0] for f in ev])
    e_mob = apply(np.array([f[1] for f in ev]))
    return rmsd_fit, float(np.sqrt(np.mean(np.sum((e_mob - e_ref) ** 2, axis=1)))), (rot, mc, rc)


def ligand_rmsd_after_receptor_fit(pairs_rec, pairs_lig):
    """CAPRI-style ligand RMSD on CA atoms: fit on the receptor pairs, evaluate on the ligand pairs."""
    _, _, tr = ca_rmsd_after_fit(pairs_rec)
    if tr is None:
        return None
    rot, mc, rc = tr
    lig = [(r['CA'].coord, m['CA'].coord) for r, m in pairs_lig if 'CA' in r and 'CA' in m]
    if not lig:
        return None
    ref = np.array([x[0] for x in lig])
    mob = (np.array([x[1] for x in lig]) - mc) @ rot + rc
    return float(np.sqrt(np.mean(np.sum((mob - ref) ** 2, axis=1))))


def scaffold_metrics(a_chain, b_chain, wt_b_chain, fixed_b_ids):
    """Quality descriptors of a diffused B scaffold against the (fixed) A' chain."""
    contacts, clashes = interface_stats(a_chain, b_chain)
    pairs = match_residues(wt_b_chain, b_chain)
    fixed = set(fixed_b_ids)
    flex = {r.id[1] for r, _ in pairs} - fixed
    fw_rmsd, flex_rmsd, _ = ca_rmsd_after_fit(pairs, fit_ids=fixed, eval_ids=flex or None)
    return {"contacts": contacts, "clashes": clashes,
            "framework_rmsd": fw_rmsd, "flex_rmsd": flex_rmsd, "n_pairs": len(pairs)}


# --------------------------------------------------------------------------- #
# Negative-design ("specificity by difference") bias for B'
# --------------------------------------------------------------------------- #
def specificity_bias(a_new, a_wt, w_pos=2.5, w_neg=2.0, w_size=1.0):
    """LigandMPNN bias (3-letter keys) for a B position contacting A' residue `a_new`
    whose WT counterpart on A was `a_wt`.

    Rewards residues complementary to the *new* A' residue, penalises those complementary to
    the *WT* residue (the negative-design term), and adds a knob/hole term when A' changed the
    side-chain size. Positions where A' == WT carry no specificity information -> {}.
    """
    if a_new == a_wt:
        return {}
    bias = {}
    for b in COMPLEMENT.get(a_new, _HYD):
        bias[b] = bias.get(b, 0.0) + w_pos
    for b in COMPLEMENT.get(a_wt, ()):
        bias[b] = bias.get(b, 0.0) - w_neg
    dv = RESIDUE_VOLUME.get(a_new, 130.0) - RESIDUE_VOLUME.get(a_wt, 130.0)
    if dv > 40:      # A' added a knob -> B' should offer a hole
        for b in _SMALL:
            bias[b] = bias.get(b, 0.0) + w_size
    elif dv < -40:   # A' made a hole -> B' fills it (would clash with the WT knob)
        for b in _LARGE:
            bias[b] = bias.get(b, 0.0) + w_size
    return {b: round(v, 2) for b, v in bias.items() if v != 0.0}


def build_rescue_bias(a_chain, b_chain, wt_a_chain, fixed_b_ids, cutoff=7.5):
    """Per-residue bias dict {"B123": {"ASP": 2.5, ...}} for the mutable positions of B.

    A' and the WT chain A must share the same frame (Module 4 superimposes A' onto WT).
    """
    wt_by_id = {r.id[1]: r for r in std_residues(wt_a_chain)}
    a_res = std_residues(a_chain)
    a_atoms = [(r, heavy_atoms(r)) for r in a_res]
    ns = NeighborSearch([at for _, atoms in a_atoms for at in atoms]) if a_atoms else None
    fixed = set(fixed_b_ids)
    bias = {}
    for res_b in std_residues(b_chain):
        if res_b.id[1] in fixed or ns is None:
            continue
        best, best_d = None, float('inf')
        for atom in heavy_atoms(res_b):
            for near in ns.search(atom.coord, cutoff):
                d = atom - near
                if d < best_d:
                    best, best_d = near.get_parent(), d
        if best is None:
            continue
        wt_res = wt_by_id.get(best.id[1])
        wt_name = wt_res.get_resname() if wt_res is not None else best.get_resname()
        spec = specificity_bias(best.get_resname(), wt_name)
        if not spec:
            # A' kept the WT residue here: fall back to a mild generic complementarity bias
            spec = {b: 1.0 for b in COMPLEMENT.get(best.get_resname(), _HYD)}
        bias[f"{b_chain.id}{res_b.id[1]}"] = spec
    return bias


def candidate_specificity_score(a_prime_chain, a_wt_chain, b_chain):
    """Cheap rigid-body proxy for orthogonality of a designed B' (same frame for all chains).

    Returns dict with contacts/clashes vs A' and vs A_WT and a scalar score (higher = better):
    keeps contacts with A', loses contacts with A_WT, and rewards steric clashes with A_WT.
    """
    c_new, k_new = interface_stats(a_prime_chain, b_chain)
    c_wt, k_wt = interface_stats(a_wt_chain, b_chain)
    score = (c_new - c_wt) + 0.5 * min(k_wt, 10) - 2.0 * k_new
    return {"contacts_a_prime": c_new, "clashes_a_prime": k_new,
            "contacts_a_wt": c_wt, "clashes_a_wt": k_wt, "spec_score": float(score)}


def parse_mpnn_confidence(pdb_path):
    """Best-effort read of a LigandMPNN 'overall_confidence' for an output structure.

    Looks for a sibling FASTA (seqs/<name>.fa) with `overall_confidence=` headers; returns
    None when unavailable so callers can degrade gracefully.
    """
    stem = os.path.splitext(os.path.basename(pdb_path))[0]
    stem = stem.replace("_packed", "").split("_packed")[0]
    root = os.path.dirname(pdb_path)
    for base in {root, os.path.dirname(root)}:
        for cand in (os.path.join(base, "seqs", stem + ".fa"), os.path.join(base, stem + ".fa")):
            if os.path.exists(cand):
                try:
                    with open(cand) as f:
                        for line in f:
                            if line.startswith(">") and "overall_confidence=" in line:
                                val = line.split("overall_confidence=")[1].split(",")[0].strip()
                                return float(val)
                except (OSError, ValueError):
                    return None
    return None


def nan_to_none(x):
    return None if x is None or (isinstance(x, float) and math.isnan(x)) else x
