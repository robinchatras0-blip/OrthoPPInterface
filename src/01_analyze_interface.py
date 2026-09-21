import argparse
import json
import os
import sys
from Bio.PDB import PDBParser

sys.path.append(os.path.dirname(__file__))
from design_utils import THREE_TO_ONE, heavy_atoms, load_config, min_residue_distance


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
    max_hotspots = config['structural_constraints']['max_hotspot_mutations']
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

    # Calculate interface residues on Chain A
    interface_A = []
    for res_A in residues_A:
        for res_B in residues_B:
            if min_residue_distance(res_A, res_B) <= dist_thresh:
                interface_A.append(res_A)
                break

    if not interface_A:
        print("Warning: No interface residues found with current distance threshold. Using first 3 residues.")
        interface_A = residues_A[:3]

    # Select hotspots (charged/polar contact residues)
    charged_polar = ['ARG', 'LYS', 'ASP', 'GLU', 'HIS', 'ASN', 'GLN', 'SER', 'THR', 'TYR']
    hotspots = [res for res in interface_A if res.get_resname() in charged_polar]
    if not hotspots:
        hotspots = interface_A
    hotspots = hotspots[:max_hotspots]
    hotspot_ids = [res.id[1] for res in hotspots]

    # Define mutable residues on Chain A (ALL interface residues + immediate spatial neighbors)
    mutable_A = set(interface_A)
    for res_A in residues_A:
        for if_res in interface_A:
            if min_residue_distance(res_A, if_res) <= neighborhood_radius:
                mutable_A.add(res_A)
                break

    neighborhood_ids = sorted(list(set([res.id[1] for res in mutable_A])))

    # Fixed positions for LigandMPNN (freeze everything EXCEPT mutable interface residues)
    fixed_ids_A = [res.id[1] for res in residues_A if res.id[1] not in neighborhood_ids]
    fixed_ids_B = [res.id[1] for res in residues_B]
    
    # Mutable residues on B (for B' rescue design: all residues of B contacting interface_A)
    mutable_B = set()
    for res_B in residues_B:
        for if_res in interface_A:
            if min_residue_distance(res_B, if_res) <= dist_thresh + 1.0:
                mutable_B.add(res_B)
                break
    neighborhood_ids_B = sorted(list(set([res.id[1] for res in mutable_B])))

    # Comprehensive Context-Aware Soluble Rupture Bias for LigandMPNN
    # Applies H-bond Donor/Acceptor inversion, electrostatic repulsion, and soluble steric clashing
    mpnn_bias = {chain_A_id: {}}

    DONORS_B = {'ARG', 'LYS', 'HIS', 'TRP'}
    ACCEPTORS_B = {'ASP', 'GLU'}
    AMPHOTERIC_POLAR_B = {'ASN', 'GLN', 'SER', 'THR', 'TYR'}

    for res_A in interface_A:
        aid = str(res_A.id[1])
        wt_resname = res_A.get_resname()
        bias_dict = {}

        # 1. Softly discourage the wild-type amino acid and glycine
        bias_dict[wt_resname] = -1.5
        bias_dict["GLY"] = -2.0

        # 2. Find interacting partners on Chain B
        partners_B = []
        for res_B in residues_B:
            dist = min_residue_distance(res_A, res_B)
            if dist <= dist_thresh:
                partners_B.append((dist, res_B))

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
                bias_dict.update({
                    "GLU": 2.0, "ARG": 2.0, "ASN": 1.5, "TYR": 1.5
                })

        mpnn_bias[chain_A_id][aid] = bias_dict

    # Contigs specification for RFdiffusion
    res_ids_A = sorted(list(set(res.id[1] for res in residues_A)))
    interface_ids = sorted(list(set(res.id[1] for res in interface_A)))
    
    contig_parts = []
    if not res_ids_A:
        print("Error: No residues found on chain A.")
        sys.exit(1)
        
    if not interface_ids:
        print("Error: No interface residues found, cannot design 0 residues.")
        sys.exit(1)

    current_mode = 'fixed' if res_ids_A[0] not in interface_ids else 'flexible'
    current_start = res_ids_A[0]
    current_len = 0

    for i, rid in enumerate(res_ids_A):
        mode = 'flexible' if rid in interface_ids else 'fixed'
        
        if mode != current_mode:
            if current_mode == 'fixed':
                # end of fixed block
                contig_parts.append(f"{chain_A_id}{current_start}-{res_ids_A[i-1]}")
            else:
                # end of flexible block
                contig_parts.append(f"{current_len}-{current_len}")
            
            current_mode = mode
            current_start = rid
            current_len = 1
        else:
            current_len += 1

    # append the last block
    if current_mode == 'fixed':
        contig_parts.append(f"{chain_A_id}{current_start}-{res_ids_A[-1]}")
    else:
        contig_parts.append(f"{current_len}-{current_len}")

    contig_str = ",".join(contig_parts)
    
    # Cropping Chain B
    cropped_residues_B = []
    for res_B in residues_B:
        is_close = False
        for res_A in interface_A:
            if min_residue_distance(res_B, res_A) <= crop_dist:
                is_close = True
                break
        if is_close:
            cropped_residues_B.append(res_B)
            
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
                "is_hotspot": res in hotspots,
                "is_neighborhood": res in mutable_A
            }
            for idx, res in enumerate(residues_A)
        ],
        "hotspot_ids": hotspot_ids,
        "neighborhood_ids": neighborhood_ids,
        "interface_ids_A": interface_ids,
        "interface_ids_B": neighborhood_ids_B,
        "fixed_ids_A": fixed_ids_A,
        "crop_min_B": min_B,
        "crop_max_B": max_B
    }

    mpnn_fixed_positions_A_prime = {chain_A_id: fixed_ids_A, chain_B_id: fixed_ids_B}
    mpnn_fixed_positions_B_prime = {chain_A_id: [res.id[1] for res in residues_A], chain_B_id: fixed_ids_B_for_rescue}

    # Build RFD3 Contig for B' Co-adaptation (Chain A is 100% fixed, Chain B interface loops within cropped domain are flexible)
    res_ids_B_cropped = sorted(list(set(res.id[1] for res in cropped_residues_B))) if cropped_residues_B else sorted(list(set(res.id[1] for res in residues_B)))
    contig_parts_B = []
    if res_ids_B_cropped:
        curr_mode_B = 'fixed' if res_ids_B_cropped[0] not in neighborhood_ids_B else 'flexible'
        curr_start_B = res_ids_B_cropped[0]
        curr_len_B = 0
        for i, rid in enumerate(res_ids_B_cropped):
            m_B = 'flexible' if rid in neighborhood_ids_B else 'fixed'
            if m_B != curr_mode_B:
                if curr_mode_B == 'fixed':
                    contig_parts_B.append(f"{chain_B_id}{curr_start_B}-{res_ids_B_cropped[i-1]}")
                else:
                    contig_parts_B.append(f"{curr_len_B}-{curr_len_B}")
                curr_mode_B = m_B
                curr_start_B = rid
                curr_len_B = 1
            else:
                curr_len_B += 1
        if curr_mode_B == 'fixed':
            contig_parts_B.append(f"{chain_B_id}{curr_start_B}-{res_ids_B_cropped[-1]}")
        else:
            contig_parts_B.append(f"{curr_len_B}-{curr_len_B}")

    contig_rescue_str = f"{chain_A_id}{res_ids_A[0]}-{res_ids_A[-1]},/0," + ",".join(contig_parts_B)
    print(f"Generated RFD3 Rescue contigs: {contig_rescue_str}")

    rfd3_rescue_inputs = {
        "design": {
            "input": os.path.relpath(pdb_file, args.out_dir).replace("\\", "/"),
            "contig": contig_rescue_str
        }
    }

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

    print(f"Module 1 Complete: Identified {len(hotspots)} hotspots and {len(neighborhood_ids)} neighborhood residues on Chain {chain_A_id}.")

if __name__ == "__main__":
    main()
