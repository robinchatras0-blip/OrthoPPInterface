import argparse
import json
import os
import sys
from Bio.PDB import PDBParser

import numpy as np

sys.path.append(os.path.dirname(__file__))
from design_utils import (
    THREE_TO_ONE, build_contig_string, heavy_atoms, load_config,
    min_residue_distance, residue_distance_matrix,
)


def consolidate_flexible(all_ids, flexible, min_fixed):
    """Absorbs fixed runs shorter than `min_fixed` residues into the surrounding flexible region.

    Selecting the redesignable residues of B purely by distance leaves isolated 1-2 residue fixed anchors
    between short flexible stretches. RFD3 then has to thread a backbone between arbitrary pinned points
    while avoiding the fixed partner chain, and returns interpenetrating scaffolds. Runs that touch either
    terminus are kept: they anchor the chain. `min_fixed <= 1` is a no-op.
    """
    cur = set(flexible)
    while True:
        runs, start = [], 0
        for i in range(1, len(all_ids) + 1):
            if i == len(all_ids) or (all_ids[i] in cur) != (all_ids[start] in cur):
                runs.append((all_ids[start] in cur, start, i - 1))
                start = i
        grew = False
        for is_flex, i0, i1 in runs:
            if not is_flex and (i1 - i0 + 1) < min_fixed and i0 > 0 and i1 < len(all_ids) - 1:
                cur |= set(all_ids[i0:i1 + 1])
                grew = True
        if not grew:
            return sorted(cur)


def main():
    parser = argparse.ArgumentParser(description="Module 1: Interface Parsing & Hotspot Targeting")
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--out_dir', default='results/01_analysis')
    args = parser.parse_args()

    config = load_config(args.config)

    os.makedirs(args.out_dir, exist_ok=True)

    pdb_file = config['pipeline']['input_pdb']
    chain_A_id = config['pipeline']['chain_A']
    chain_B_id = config['pipeline']['chain_B']
    dist_thresh = config['structural_constraints']['interface_distance_threshold_angstroms']
    neighborhood_radius = config['structural_constraints']['neighborhood_radius_angstroms']
    crop_dist = config['structural_constraints'].get('crop_chain_B_distance_angstroms', 15.0)

    if not os.path.exists(pdb_file):
        sys.exit(f"Error: Input PDB {pdb_file} not found!")

    parser_pdb = PDBParser(QUIET=True)
    structure = parser_pdb.get_structure("complex", pdb_file)
    model = structure[0]

    if chain_A_id not in model or chain_B_id not in model:
        print(f"Error: Chains {chain_A_id} or {chain_B_id} missing in PDB {pdb_file}.")
        sys.exit(1)

    chain_A = model[chain_A_id]
    chain_B = model[chain_B_id]

    residues_A = [res for res in chain_A if res.id[0] == ' ' and len(heavy_atoms(res)) > 0]
    residues_B = [res for res in chain_B if res.id[0] == ' ' and len(heavy_atoms(res)) > 0]

    # Pre-compute inter-chain distance matrix
    dist_AB = residue_distance_matrix(residues_A, residues_B)

    # Interface residues of chain A, with the number of B residues each one contacts (the count ranks the
    # PPI hotspots handed to RFD3 for the rescue diffusion below, so it is not just a boolean test).
    interface_A, contacts_A = [], {}
    contacts_per_A = (dist_AB <= dist_thresh).sum(axis=1) if dist_AB.size else []
    for idx_A, res_A in enumerate(residues_A):
        n_contacts = int(contacts_per_A[idx_A])
        if n_contacts:
            interface_A.append(res_A)
            contacts_A[res_A.id[1]] = n_contacts

    if not interface_A:
        print("Warning: No interface residues found with current distance threshold. Using first 3 residues.")
        interface_A = residues_A[:3]

    # Pre-compute intra-chain distance matrix from all Chain A residues to interface_A
    dist_AA = residue_distance_matrix(residues_A, interface_A)

    # Define mutable residues on Chain A (ALL interface residues + immediate spatial neighbors)
    is_neighbor_A = (dist_AA <= neighborhood_radius).any(axis=1) if dist_AA.size else []
    mutable_A = {res_A for idx_A, res_A in enumerate(residues_A) if is_neighbor_A[idx_A]}
    mutable_A.update(interface_A)

    neighborhood_ids = sorted(list(set([res.id[1] for res in mutable_A])))

    # Fixed positions for LigandMPNN (freeze everything EXCEPT mutable interface residues)
    fixed_ids_A = [res.id[1] for res in residues_A if res.id[1] not in neighborhood_ids]

    # Mutable residues on B (for B' rescue design: all residues of B contacting interface_A)
    iface_set_A = set(interface_A)
    iface_idx_A = [i for i, r in enumerate(residues_A) if r in iface_set_A]
    dist_ifaceA_B = dist_AB[iface_idx_A, :] if iface_idx_A and dist_AB.size else np.empty((0, len(residues_B)), dtype=np.float32)

    is_mutable_B = (dist_ifaceA_B <= dist_thresh + 1.0).any(axis=0) if dist_ifaceA_B.size else np.array([False] * len(residues_B))
    mutable_B = {res_B for idx_B, res_B in enumerate(residues_B) if is_mutable_B[idx_B]}
    neighborhood_ids_B = sorted(list(set([res.id[1] for res in mutable_B])))
    # Turn that scattered selection into real loops before it becomes an RFD3 contig (see consolidate_flexible).
    # Measured on the reference system: 5 of 8 fixed contig segments were 1-2 residues and two diffused segments
    # were a single residue; consolidating gives 3 anchors (8/26/7 residues) and 3 loops (20/19/8) for only
    # 6 extra redesignable residues (41 -> 47 of 88).
    min_fixed_B = int(config['pipeline'].get('rescue_min_fixed_segment', 3))
    neighborhood_ids_B = consolidate_flexible([r.id[1] for r in residues_B], set(neighborhood_ids_B), min_fixed_B)

    # Comprehensive Context-Aware Soluble Rupture Bias for LigandMPNN
    # Applies H-bond Donor/Acceptor inversion, electrostatic repulsion, and soluble steric clashing
    mpnn_bias = {chain_A_id: {}}

    DONORS_B = {'ARG', 'LYS', 'HIS', 'TRP'}
    ACCEPTORS_B = {'ASP', 'GLU'}
    AMPHOTERIC_POLAR_B = {'ASN', 'GLN', 'SER', 'THR', 'TYR'}
    allow_steric = config['pipeline'].get('rupture_allow_steric_hydrophobic', True)

    for res_A in interface_A:
        aid = str(res_A.id[1])
        wt_resname = res_A.get_resname()
        bias_dict = {}

        # 1. Find interacting partners on Chain B
        idx_A = residues_A.index(res_A)
        dists_to_B = dist_AB[idx_A, :]
        partners_B = [(float(dists_to_B[j]), res_B) for j, res_B in enumerate(residues_B) if dists_to_B[j] <= dist_thresh]

        partners_B.sort(key=lambda x: x[0])
        if not partners_B:
            # General soluble disruption bias
            bias_dict.update({"ARG": 2.0, "GLU": 2.0, "TYR": 1.5, "GLN": 1.5})
        else:
            closest_dist, closest_B = partners_B[0]
            b_name = closest_B.get_resname()

            # Rule A: B is a H-bond DONOR / Positively Charged -> Charge Inversion + Tyrosine Bulk
            if b_name in DONORS_B:
                bias_dict.update({
                    "ARG": 2.5, "LYS": 2.0, "TYR": 1.5,
                    "ASP": -2.5, "GLU": -2.5
                })
            # Rule B: B is a H-bond ACCEPTOR / Negatively Charged -> Charge Inversion + Tyrosine Bulk
            elif b_name in ACCEPTORS_B:
                bias_dict.update({
                    "ASP": 2.5, "GLU": 2.5, "TYR": 1.5,
                    "ARG": -2.5, "LYS": -2.5
                })
            # Rule C: B is Amphoteric / Polar (Asn, Gln, Ser, Thr, Tyr) -> Bulky Polar & Soluble Steric Clashing
            elif b_name in AMPHOTERIC_POLAR_B:
                bias_dict.update({
                    "ARG": 2.0, "GLN": 2.0, "TYR": 1.5, "GLU": 1.5
                })
            # Rule D: B is Hydrophobic -> Insert Soluble Charge into the pocket (forces severe desolvation penalty for B_WT)
            else:
                if allow_steric and contacts_A.get(res_A.id[1], 0) >= 4:
                    if wt_resname in ('ALA', 'VAL', 'SER', 'GLY'):
                        bias_dict.update({"TRP": 2.0, "PHE": 2.0, "TYR": 2.0})
                    elif wt_resname in ('TRP', 'PHE', 'LEU', 'ILE', 'MET'):
                        bias_dict.update({"ALA": 1.5, "VAL": 1.5, "GLU": 1.5, "ARG": 1.5})
                    else:
                        bias_dict.update({"GLU": 2.0, "ARG": 2.0, "ASN": 1.5, "TYR": 1.5})
                else:
                    bias_dict.update({
                        "GLU": 2.0, "ARG": 2.0, "ASN": 1.5, "TYR": 1.5
                    })

        # 2. Discourage the wild-type residue and glycine LAST, and never leave them positive. Applying the
        # penalty first let a rule that happens to recommend the wild-type residue overwrite it and actively
        # favour keeping the wild type - the opposite of the rupture objective (measured on the reference
        # system: 8 of 36 interface positions, all ARG/TYR/GLU/ASN, i.e. exactly the contact-making residues
        # the rules recommend). Adding the penalty is not enough either, since a +2.0 rule bonus survives a
        # -1.5 penalty; a rule that recommends the residue already there is self-contradictory for a rupture,
        # so the rule bonus is dropped for those two residues rather than merely reduced.
        for res_name, penalty in ((wt_resname, -1.5), ("GLY", -2.0)):
            bias_dict[res_name] = round(min(bias_dict.get(res_name, 0.0), 0.0) + penalty, 2)

        mpnn_bias[chain_A_id][aid] = bias_dict

    # Contigs specification for RFdiffusion
    res_ids_A = sorted(list(set(res.id[1] for res in residues_A)))
    interface_ids = sorted(list(set(res.id[1] for res in interface_A)))
    
    if not res_ids_A:
        print("Error: No residues found on chain A.")
        sys.exit(1)
        
    if not interface_ids:
        print("Error: No interface residues found, cannot design 0 residues.")
        sys.exit(1)

    # Consolidate flexible interface segments of Chain A (absorbs 1-2 residue fixed anchors within loops).
    # This prevents rigid pinning of interface loops (e.g. A8-8, A11-11, A36-36) while strictly keeping
    # the entire hydrophobic core and non-interface framework locked.
    min_fixed_A = int(config['pipeline'].get('rupture_min_fixed_segment', 3))
    flexible_ids_A = set(consolidate_flexible(res_ids_A, set(interface_ids), min_fixed_A))

    contig_str = build_contig_string(chain_A_id, res_ids_A, flexible_ids_A)
    
    # Cropping Chain B
    is_cropped_B = (dist_ifaceA_B <= crop_dist).any(axis=0) if dist_ifaceA_B.size else np.array([False] * len(residues_B))
    cropped_residues_B = [res_B for idx_B, res_B in enumerate(residues_B) if is_cropped_B[idx_B]]
            
    if cropped_residues_B:
        min_B = cropped_residues_B[0].id[1]
        max_B = cropped_residues_B[-1].id[1]
    else:
        min_B = residues_B[0].id[1]
        max_B = residues_B[-1].id[1]

    # Append cropped Chain B context
    contig_str += f",/0,{chain_B_id}{min_B}-{max_B}"
    
    # Update fixed_ids_B and fixed_ids_B_for_rescue for full complex
    fixed_ids_B = [res.id[1] for res in residues_B]
    fixed_ids_B_for_rescue = [res.id[1] for res in residues_B if res.id[1] not in neighborhood_ids_B]
        
    print(f"Generated RFD3 contigs: {contig_str}")
    
    # No PPI hotspots here on purpose: hotspots are a *positive* design signal ("engage these residues"),
    # which is the opposite of the rupture objective for A'.
    rfd3_inputs = {
        "design": {
            "input": os.path.relpath(pdb_file, args.out_dir).replace("\\", "/"),
            "contig": contig_str
        }
    }

    # Structural mapping index dictionary
    index_mapping = {
        "chain_A": chain_A_id,
        "chain_B": chain_B_id,
        "residues_A": [
            {
                "sequential_index": idx + 1,
                "resseq": res.id[1],
                "resname": res.get_resname(),
                "one_letter": THREE_TO_ONE.get(res.get_resname(), 'X'),
                "is_interface": res in interface_A,
                "is_neighborhood": res in mutable_A
            }
            for idx, res in enumerate(residues_A)
        ],
        "neighborhood_ids": neighborhood_ids,
        "interface_ids_A": interface_ids,
        "flexible_ids_A": sorted(list(flexible_ids_A)),
        "interface_ids_B": neighborhood_ids_B,
        "interface_contacts_A": contacts_A,
        "fixed_ids_A": fixed_ids_A,
        "crop_min_B": min_B,
        "crop_max_B": max_B
    }

    mpnn_fixed_positions_A_prime = {chain_A_id: fixed_ids_A, chain_B_id: fixed_ids_B}
    mpnn_fixed_positions_B_prime = {chain_A_id: [res.id[1] for res in residues_A], chain_B_id: fixed_ids_B_for_rescue}

    # Build RFD3 Rescue specification for B' Co-adaptation
    # Chain A is 100% fixed, Chain B framework is fixed, while Chain B interface loops within cropped domain are flexible.
    partial_t = config['pipeline'].get('rescue_partial_t')

    def _to_contiguous_ranges(ids):
        if not ids:
            return []
        sorted_ids = sorted(list(set(ids)))
        ranges = []
        s = e = sorted_ids[0]
        for x in sorted_ids[1:]:
            if x == e + 1:
                e = x
            else:
                ranges.append((s, e))
                s = e = x
        ranges.append((s, e))
        return ranges

    fixed_ranges_B = _to_contiguous_ranges(fixed_ids_B_for_rescue)
    fixed_b_str = ",".join(f"{chain_B_id}{s}-{e}" for s, e in fixed_ranges_B)
    select_fixed_atoms_str = f"{chain_A_id}{res_ids_A[0]}-{res_ids_A[-1]}"
    if fixed_b_str:
        select_fixed_atoms_str += f",{fixed_b_str}"

    n_hotspots = int(config['pipeline'].get('rescue_ppi_hotspots', 8))
    hotspot_ids_A = sorted(sorted(contacts_A, key=lambda rid: (-contacts_A[rid], rid))[:n_hotspots])
    index_mapping["rescue_hotspot_ids_A"] = hotspot_ids_A

    if partial_t is not None and float(partial_t) > 0.0:
        # Guided partial diffusion: Chain B specified as real residue range,
        # with framework atoms fixed and interface loop atoms noised by partial_t Angstroms.
        contig_rescue_str = f"{chain_A_id}{res_ids_A[0]}-{res_ids_A[-1]},/0,{chain_B_id}{min_B}-{max_B}"
        rfd3_rescue_inputs = {
            "design": {
                "input": os.path.relpath(pdb_file, args.out_dir).replace("\\", "/"),
                "contig": contig_rescue_str,
                "select_fixed_atoms": select_fixed_atoms_str,
                "partial_t": float(partial_t),
            }
        }
        print(f"Generated RFD3 Rescue contigs (Partial Diffusion t={partial_t} A): {contig_rescue_str}")
        print(f"  Fixed atoms: {select_fixed_atoms_str}")
    else:
        # Fallback de-novo loop hallucination (infer_ori_strategy: hotspots)
        res_ids_B_cropped = sorted(list(set(res.id[1] for res in cropped_residues_B))) if cropped_residues_B else sorted(list(set(res.id[1] for res in residues_B)))
        contig_parts_B_str = build_contig_string(chain_B_id, res_ids_B_cropped, neighborhood_ids_B)
        contig_rescue_str = f"{chain_A_id}{res_ids_A[0]}-{res_ids_A[-1]},/0,{contig_parts_B_str}"
        rfd3_rescue_inputs = {
            "design": {
                "input": os.path.relpath(pdb_file, args.out_dir).replace("\\", "/"),
                "contig": contig_rescue_str,
                "select_hotspots": ",".join(f"{chain_A_id}{rid}" for rid in hotspot_ids_A),
                "infer_ori_strategy": "hotspots",
            }
        }
        print(f"Generated RFD3 Rescue contigs (De-novo): {contig_rescue_str}")
        print(f"Rescue PPI hotspots on chain {chain_A_id} (contact-ranked): "
              f"{hotspot_ids_A} (contacts {[contacts_A[r] for r in hotspot_ids_A]})")

    # Save outputs
    with open(os.path.join(args.out_dir, "inputs.json"), 'w') as f:
        json.dump(rfd3_inputs, f, indent=2)

    with open(os.path.join(args.out_dir, "inputs_rescue.json"), 'w') as f:
        json.dump(rfd3_rescue_inputs, f, indent=2)

    with open(os.path.join(args.out_dir, "mpnn_fixed_positions.json"), 'w') as f:
        json.dump(mpnn_fixed_positions_A_prime, f, indent=2)
        
    with open(os.path.join(args.out_dir, "mpnn_fixed_positions_B.json"), 'w') as f:
        json.dump(mpnn_fixed_positions_B_prime, f, indent=2)

    with open(os.path.join(args.out_dir, "mpnn_bias.json"), 'w') as f:
        json.dump(mpnn_bias, f, indent=2)

    with open(os.path.join(args.out_dir, "index_mapping.json"), 'w') as f:
        json.dump(index_mapping, f, indent=2)

    print(f"Module 1 Complete: {len(interface_ids)} interface and {len(neighborhood_ids)} redesignable residues on chain "
          f"{chain_A_id}, {len(neighborhood_ids_B)} redesignable on chain {chain_B_id}.")

if __name__ == "__main__":
    main()
