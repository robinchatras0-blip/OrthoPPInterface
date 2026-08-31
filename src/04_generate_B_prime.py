import argparse
import yaml
import subprocess
import sys
import os
import json
import glob
import shutil
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Superimposer

THREE_TO_ONE = {
    'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F',
    'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',
    'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R',
    'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'
}

def extract_interface_sequence(pdb_path, chain_id, interface_indices):
    parser = PDBParser(QUIET=True)
    st = parser.get_structure("st", pdb_path)
    seq_dict = {}
    for model in st:
        for chain in model:
            if chain.id == chain_id:
                for res in chain:
                    if res.id[0] == ' ' and res.id[1] in interface_indices:
                        seq_dict[res.id[1]] = THREE_TO_ONE.get(res.get_resname(), 'X')
        break
    return seq_dict

def compute_sequence_distance(seq1, seq2):
    all_keys = set(seq1.keys()) | set(seq2.keys())
    if not all_keys:
        return 0
    return sum(1 for k in all_keys if seq1.get(k) != seq2.get(k))

def select_diverse_candidates(candidate_pdbs, chain_id, interface_indices, max_k=3):
    if len(candidate_pdbs) <= max_k:
        return candidate_pdbs
    
    seqs = [extract_interface_sequence(p, chain_id, interface_indices) for p in candidate_pdbs]
    selected_indices = [0] # Always keep top-scoring candidate 0
    
    while len(selected_indices) < max_k:
        best_candidate = None
        max_min_dist = -1
        for i in range(len(candidate_pdbs)):
            if i in selected_indices:
                continue
            min_dist_to_selected = min(compute_sequence_distance(seqs[i], seqs[s]) for s in selected_indices)
            if min_dist_to_selected > max_min_dist:
                max_min_dist = min_dist_to_selected
                best_candidate = i
        
        if best_candidate is not None:
            selected_indices.append(best_candidate)
        else:
            for i in range(len(candidate_pdbs)):
                if i not in selected_indices:
                    selected_indices.append(i)
                    break
        if len(selected_indices) == len(candidate_pdbs):
            break
            
    print(f"  [Diversity Selector] Selected {len(selected_indices)} diverse A' motifs (indices: {selected_indices}) from {len(candidate_pdbs)} passing candidates.")
    return [candidate_pdbs[i] for i in selected_indices]

def main():
    parser = argparse.ArgumentParser(description="Module 4: Diverse Rescue B' Generation")
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--analysis_dir', default='results/01_analysis')
    parser.add_argument('--passed_candidates', default='results/03_fail_fast/passed_candidates.txt')
    parser.add_argument('--out_dir', default='results/04_rescue_design')
    args = parser.parse_args()

    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    os.makedirs(args.out_dir, exist_ok=True)

    with open(args.passed_candidates, 'r') as f:
        raw_lines = [l.strip() for l in f.readlines() if l.strip()]
    if not raw_lines:
        raise ValueError("No passed candidates found!")

    # Resolve paths
    resolved_candidates = []
    cand_dir = os.path.dirname(os.path.abspath(args.passed_candidates))
    parent_dir = os.path.dirname(cand_dir)

    for p in raw_lines:
        if os.path.exists(p):
            resolved_candidates.append(p)
        elif os.path.exists(os.path.join(cand_dir, p)):
            resolved_candidates.append(os.path.join(cand_dir, p))
        elif os.path.exists(os.path.join(parent_dir, "02_rupture_design", os.path.basename(p))):
            resolved_candidates.append(os.path.join(parent_dir, "02_rupture_design", os.path.basename(p)))
        elif os.path.exists(os.path.join("results/02_rupture_design", os.path.basename(p))):
            resolved_candidates.append(os.path.join("results/02_rupture_design", os.path.basename(p)))

    execution_mode = config['pipeline'].get('execution_mode', 'mock')
    chain_A_id = config['pipeline'].get('chain_A', 'A')
    chain_B_id = config['pipeline'].get('chain_B', 'B')
    wt_pdb = config['pipeline'].get('input_pdb', 'data/inputs/complex_S1_S2.pdb')

    # Load dynamic cropped residue range for Chain B from Module 1 index mapping
    mapping_file = os.path.join(args.analysis_dir, 'index_mapping.json')
    crop_min_b = 1
    crop_max_b = 88
    interface_indices_A = []
    if os.path.exists(mapping_file):
        with open(mapping_file, 'r') as f:
            mapping_data = json.load(f)
            crop_min_b = mapping_data.get('crop_min_B', 1)
            crop_max_b = mapping_data.get('crop_max_B', 88)
            interface_indices_A = mapping_data.get('neighborhood_ids', [])

    # Diversity Selection
    div_cfg = config.get('diversity_selection', {})
    div_enabled = div_cfg.get('enabled', True)
    max_k_motifs = div_cfg.get('max_a_prime_motifs', 3)

    if div_enabled and len(resolved_candidates) > 1:
        selected_candidates = select_diverse_candidates(resolved_candidates, chain_A_id, interface_indices_A, max_k=max_k_motifs)
    else:
        selected_candidates = resolved_candidates[:max_k_motifs]

    lmpnn_cfg = config.get('ligandmpnn', {})
    coadapt_enabled = lmpnn_cfg.get('rescue_coadaptation', True)
    rescue_n_batches = str(lmpnn_cfg.get('rescue_n_batches', 3))
    temp_rescue = str(lmpnn_cfg.get('temperature_rescue', 0.10))
    model_type = lmpnn_cfg.get('model_type', "ligand_mpnn")
    is_legacy = str(lmpnn_cfg.get('is_legacy_weights', "True"))
    default_ckpt = "OrthoIntRob/ligandmpnn/model_params/ligandmpnn_v_32_010_25.pt" if execution_mode == 'local' else "/content/drive/MyDrive/OrthoInterface_Data/FoundryModels/LigandMPNN/model_params/ligandmpnn_v_32_010_25.pt"
    checkpoint_path = lmpnn_cfg.get('checkpoint_path', default_ckpt)

    metadata_pairs = {}
    all_saved_b_pdbs = []

    print(f"\nModule 4: Co-adapting and designing rescue B' for {len(selected_candidates)} diverse A' motifs...")

    for a_idx, in_pdb in enumerate(selected_candidates):
        a_cand_name = os.path.splitext(os.path.basename(in_pdb))[0]
        print(f"\n============================================================")
        print(f"  [A' Motif {a_idx+1}/{len(selected_candidates)}]: {a_cand_name}")
        print(f"============================================================")

        motif_work_dir = os.path.join(args.out_dir, f"motif_{a_cand_name}")
        os.makedirs(motif_work_dir, exist_ok=True)

        parser = PDBParser(QUIET=True)
        struct_A = parser.get_structure("A_prime", in_pdb)
        struct_wt = parser.get_structure("WT", wt_pdb)

        # Superimpose designed Chain A back onto the WT Chain A spatial frame
        wt_a_ca = [r['CA'] for r in struct_wt[0][chain_A_id] if r.id[0] == ' ' and 'CA' in r]
        cand_a_ca = [r['CA'] for r in struct_A[0][chain_A_id] if r.id[0] == ' ' and 'CA' in r]
        if wt_a_ca and cand_a_ca:
            min_len = min(len(wt_a_ca), len(cand_a_ca))
            sup = Superimposer()
            sup.set_atoms(wt_a_ca[:min_len], cand_a_ca[:min_len])
            sup.apply(struct_A.get_atoms())
            print(f"  Superimposed {a_cand_name} onto WT coordinate frame (RMSD = {sup.rms:.3f} Å)")

        complex_struct = Structure.Structure("complex_A_prime_B_wt")
        complex_model = Model.Model(0)
        complex_struct.add(complex_model)

        # Add redesigned Chain A
        for model in struct_A:
            for chain in model:
                if chain.id == chain_A_id:
                    complex_model.add(chain.copy())
                    break
            break

        # Add cropped wild-type Chain B for RFD3
        b_chain_crop = Chain.Chain(chain_B_id)
        complex_model.add(b_chain_crop)
        wt_chain_B_res = []
        for model in struct_wt:
            for chain in model:
                if chain.id == chain_B_id:
                    wt_chain_B_res = [res for res in chain if res.id[0] == ' ']
                    for res in wt_chain_B_res:
                        if crop_min_b <= res.id[1] <= crop_max_b:
                            b_chain_crop.add(res.copy())
                    break
            break

        complex_pdb = os.path.join(motif_work_dir, "complex_cropped.pdb")
        io = PDBIO()
        io.set_structure(complex_struct)
        io.save(complex_pdb)

        rfd_input_pdb = complex_pdb

        # Step 1: RFD3 Co-adaptation
        if coadapt_enabled and execution_mode != 'mock':
            print(f"  Running RFD3 Backbone Co-adaptation on B around {a_cand_name}...")
            rfd_rescue_dir = os.path.join(motif_work_dir, "rfd3_coadapt_out")
            os.makedirs(rfd_rescue_dir, exist_ok=True)

            inputs_rescue_src = os.path.join(args.analysis_dir, "inputs_rescue.json")
            with open(inputs_rescue_src, 'r') as f:
                rescue_json = json.load(f)

            rescue_json['design']['input'] = os.path.relpath(complex_pdb, rfd_rescue_dir).replace('\\', '/')
            curr_inputs_json = os.path.join(rfd_rescue_dir, "inputs.json").replace('\\', '/')
            with open(curr_inputs_json, 'w') as f:
                json.dump(rescue_json, f, indent=2)

            rfd_bin = config['pipeline'].get('local_rfdiffusion', 'OrthoIntRob/bin/rfd3') if execution_mode == 'local' else config['pipeline']['colab_rfdiffusion']
            rescue_dir_clean = rfd_rescue_dir.replace('\\', '/')
            rfd_cmd = rfd_bin.split() + [
                f"inputs={curr_inputs_json}",
                f"out_dir={rescue_dir_clean}",
                "n_batches=1",
                "diffusion_batch_size=1"
            ]
            env = os.environ.copy()
            env['PYTHONUNBUFFERED'] = '1'
            subprocess.run(rfd_cmd, check=True, env=env)

            rfd_pdbs = sorted(glob.glob(os.path.join(rfd_rescue_dir, "*.pdb")))
            rfd_cifs = sorted(glob.glob(os.path.join(rfd_rescue_dir, "*.cif.gz")) + glob.glob(os.path.join(rfd_rescue_dir, "*.cif")))

            if not rfd_pdbs and rfd_cifs:
                import gzip
                from Bio.PDB import MMCIFParser
                for cif_path in rfd_cifs:
                    pdb_path = cif_path.replace(".cif.gz", ".pdb").replace(".cif", ".pdb")
                    cif_parser = MMCIFParser(QUIET=True)
                    if cif_path.endswith('.gz'):
                        with gzip.open(cif_path, 'rt') as f:
                            cif_struct = cif_parser.get_structure("coadapt", f)
                    else:
                        with open(cif_path, 'rt') as f:
                            cif_struct = cif_parser.get_structure("coadapt", f)
                    io_c = PDBIO()
                    io_c.set_structure(cif_struct)
                    io_c.save(pdb_path)
                    rfd_pdbs.append(pdb_path)

            if rfd_pdbs:
                rfd_input_pdb = rfd_pdbs[0]

        # Reconstruct full-length Chain B (grafting constant domain tail)
        if wt_chain_B_res:
            tail_res = [r for r in wt_chain_B_res if r.id[1] > crop_max_b]
            if tail_res:
                p_full = PDBParser(QUIET=True)
                st_full = p_full.get_structure("assembled", rfd_input_pdb)
                m_full = st_full[0]
                
                mpnn_fixed_pos = os.path.join(args.analysis_dir, 'mpnn_fixed_positions_B.json')
                with open(mpnn_fixed_pos) as f:
                    old_fixed = json.load(f)
                fixed_b_set = set(old_fixed.get(chain_B_id, []))
                
                fixed_crop_ids = set([res.id[1] for res in wt_chain_B_res if crop_min_b <= res.id[1] <= crop_max_b and res.id[1] in fixed_b_set])
                wt_b_ca = [r['CA'] for r in struct_wt[0][chain_B_id] if r.id[1] in fixed_crop_ids and 'CA' in r]
                rfd_b_ca = [r['CA'] for r in m_full[chain_B_id] if r.id[1] in fixed_crop_ids and 'CA' in r]
                
                if wt_b_ca and rfd_b_ca:
                    min_len_b = min(len(wt_b_ca), len(rfd_b_ca))
                    sup_b = Superimposer()
                    sup_b.set_atoms(rfd_b_ca[:min_len_b], wt_b_ca[:min_len_b])
                    sup_b.apply(struct_wt[0][chain_B_id].get_atoms())
                
                if chain_B_id in m_full:
                    for r in struct_wt[0][chain_B_id]:
                        if r.id[0] == ' ' and r.id[1] > crop_max_b:
                            m_full[chain_B_id].add(r.copy())
                assembled_pdb = os.path.join(motif_work_dir, "complex_assembled_full.pdb")
                io_full = PDBIO()
                io_full.set_structure(st_full)
                io_full.save(assembled_pdb)
                rfd_input_pdb = assembled_pdb

        # Step 2: LigandMPNN Sequence Design
        mpnn_fixed_pos = os.path.join(args.analysis_dir, 'mpnn_fixed_positions_B.json')
        with open(mpnn_fixed_pos) as f:
            old_fixed = json.load(f)
        fixed_list = []
        for chain, res_list in old_fixed.items():
            fixed_list.extend([f"{chain}{res}" for res in res_list])

        # Filter fixed residues to only those actually present in rfd_input_pdb
        parser_chk_b = PDBParser(QUIET=True)
        st_chk_b = parser_chk_b.get_structure("chk_b", rfd_input_pdb)
        existing_res_b = set(f"{ch.id}{r.id[1]}" for m in st_chk_b for ch in m for r in ch if r.id[0] == ' ')
        fixed_list = [f for f in fixed_list if f in existing_res_b]
        fixed_str = ",".join(fixed_list)

        # Compute Complementary Rescue Bias
        def get_heavy_atoms(res):
            return [atom for atom in res if not atom.get_name().startswith('H') and atom.element != 'H']

        def min_dist_res(r1, r2):
            d_min = float('inf')
            for a1 in get_heavy_atoms(r1):
                for a2 in get_heavy_atoms(r2):
                    d = a1 - a2
                    if d < d_min:
                        d_min = d
            return d_min

        model_c = st_chk_b[0]
        chain_A_res = [r for r in model_c[chain_A_id] if r.id[0] == ' ']
        chain_B_res = [r for r in model_c[chain_B_id] if r.id[0] == ' ']
        fixed_b_set = set(old_fixed.get(chain_B_id, []))
        rescue_bias = {}

        for res_B in chain_B_res:
            bid = res_B.id[1]
            if bid in fixed_b_set:
                continue
            closest_A = None
            min_d = float('inf')
            for res_A in chain_A_res:
                d = min_dist_res(res_B, res_A)
                if d < min_d:
                    min_d = d
                    closest_A = res_A

            if closest_A is not None and min_d <= 7.5:
                a_resname = closest_A.get_resname()
                b_key = f"{chain_B_id}{bid}"
                if a_resname in {'ARG', 'LYS', 'HIS'}:
                    b_bias = {"ASP": 5.0, "GLU": 5.0, "ARG": -4.0, "LYS": -4.0}
                elif a_resname in {'ASP', 'GLU'}:
                    b_bias = {"ARG": 5.0, "LYS": 5.0, "ASP": -4.0, "GLU": -4.0}
                elif a_resname in {'PHE', 'TYR', 'TRP'}:
                    b_bias = {"LEU": 3.0, "ILE": 3.0, "VAL": 3.0, "PHE": 2.0, "TYR": 2.0}
                elif a_resname in {'ASN', 'GLN', 'SER', 'THR'}:
                    b_bias = {"GLN": 3.0, "ASN": 3.0, "SER": 2.0, "THR": 2.0, "ARG": 2.0, "GLU": 2.0}
                else:
                    b_bias = {"LEU": 2.0, "VAL": 2.0, "ALA": 1.5, "ILE": 2.0}
                rescue_bias[b_key] = b_bias

        mpnn_out_dir = os.path.join(motif_work_dir, "mpnn_out")
        os.makedirs(mpnn_out_dir, exist_ok=True)

        if execution_mode == 'mock':
            for b_idx in range(int(rescue_n_batches)):
                dest_cand_pdb = os.path.join(args.out_dir, f"{a_cand_name}_B_prime_cand_{b_idx:02d}.pdb")
                with open(dest_cand_pdb, 'w') as f:
                    f.write("DUMMY COMPLEX PDB\n")
                all_saved_b_pdbs.append(dest_cand_pdb)
                metadata_pairs[os.path.splitext(os.path.basename(dest_cand_pdb))[0]] = {"parent_a": a_cand_name, "pdb": dest_cand_pdb}
        else:
            mpnn_bin = config['pipeline'].get('local_ligandmpnn', 'OrthoIntRob/bin/mpnn') if execution_mode == 'local' else config['pipeline']['colab_ligandmpnn']
            mpnn_cmd = mpnn_bin.split() + [
                "--structure_path", rfd_input_pdb.replace('\\', '/'),
                "--out_directory", mpnn_out_dir.replace('\\', '/'),
                "--model_type", model_type,
                "--checkpoint_path", checkpoint_path.replace('\\', '/'),
                "--is_legacy_weights", is_legacy,
                "--batch_size", "1",
                "--number_of_batches", rescue_n_batches,
                "--temperature", temp_rescue,
                "--fixed_residues", fixed_str,
                "--bias_per_residue", json.dumps(rescue_bias),
                "--write_structures", "True"
            ]
            print(f"  Running LigandMPNN ({rescue_n_batches} variants)...")
            env = os.environ.copy()
            subprocess.run(mpnn_cmd, check=True, env=env)

            mpnn_outputs = [f for f in glob.glob(os.path.join(mpnn_out_dir, "**", "*.pdb"), recursive=True) if not f.endswith("B_prime_rfd_mpnn.pdb") and "candidate" not in f]
            mpnn_cifs = glob.glob(os.path.join(mpnn_out_dir, "**", "*.cif.gz"), recursive=True) + glob.glob(os.path.join(mpnn_out_dir, "**", "*.cif"), recursive=True)
            mpnn_outputs = [f for f in mpnn_outputs if os.path.abspath(f) != os.path.abspath(rfd_input_pdb)]
            mpnn_cifs = [f for f in mpnn_cifs if os.path.abspath(f) != os.path.abspath(rfd_input_pdb)]

            packed_outputs = [f for f in mpnn_outputs if "packed" in f]
            if packed_outputs:
                mpnn_outputs = packed_outputs

            if mpnn_outputs:
                for b_idx, cand_pdb in enumerate(sorted(mpnn_outputs)):
                    pair_id = f"{a_cand_name}_B_cand_{b_idx:02d}"
                    dest_cand_pdb = os.path.join(args.out_dir, f"{pair_id}.pdb")
                    shutil.copy(cand_pdb, dest_cand_pdb)
                    with open(dest_cand_pdb, 'r') as f:
                        lines = f.readlines()
                    with open(dest_cand_pdb, 'w') as f:
                        for line in lines:
                            if 'nan' not in line:
                                f.write(line)
                    all_saved_b_pdbs.append(dest_cand_pdb)
                    metadata_pairs[pair_id] = {"parent_a": a_cand_name, "pdb": dest_cand_pdb}
                    
                    # Backward compatibility aliases for primary motif
                    if a_idx == 0:
                        shutil.copy(dest_cand_pdb, os.path.join(args.out_dir, f"B_prime_candidate_{b_idx:02d}.pdb"))
                        if b_idx == 0:
                            shutil.copy(dest_cand_pdb, os.path.join(args.out_dir, "B_prime_rfd_mpnn.pdb"))

    # Save metadata dictionary for Module 5
    metadata_json_path = os.path.join(args.out_dir, "diversity_pairs_metadata.json")
    with open(metadata_json_path, 'w') as f:
        json.dump(metadata_pairs, f, indent=2)

    print(f"\nModule 4 Complete: Successfully generated {len(all_saved_b_pdbs)} total orthogonal pairs across {len(selected_candidates)} diverse A' motifs.")

if __name__ == "__main__":
    main()
