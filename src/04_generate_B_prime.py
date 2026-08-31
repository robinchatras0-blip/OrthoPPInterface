import argparse
import yaml
import subprocess
import sys
import os
import json

def main():
    parser = argparse.ArgumentParser(description="Module 4: Rescue B' Generation")
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--analysis_dir', default='results/01_analysis')
    parser.add_argument('--passed_candidates', default='results/03_fail_fast/passed_candidates.txt')
    parser.add_argument('--out_dir', default='results/04_rescue_design')
    args = parser.parse_args()

    # Read the first passed candidate
    with open(args.passed_candidates, 'r') as f:
        lines = f.readlines()
    if not lines:
        raise ValueError("No passed candidates found!")
    
    in_pdb = lines[0].strip()
    if not os.path.exists(in_pdb):
        # Resolve relative to passed_candidates directory or parent run structure
        cand_dir = os.path.dirname(os.path.abspath(args.passed_candidates))
        parent_dir = os.path.dirname(cand_dir)
        if os.path.exists(os.path.join(cand_dir, in_pdb)):
            in_pdb = os.path.join(cand_dir, in_pdb)
        elif os.path.exists(os.path.join(parent_dir, "02_rupture_design", os.path.basename(in_pdb))):
            in_pdb = os.path.join(parent_dir, "02_rupture_design", os.path.basename(in_pdb))
        elif os.path.exists(os.path.join("results/02_rupture_design", os.path.basename(in_pdb))):
            in_pdb = os.path.join("results/02_rupture_design", os.path.basename(in_pdb))

    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    os.makedirs(args.out_dir, exist_ok=True)

    execution_mode = config['pipeline'].get('execution_mode', 'mock')
    apptainer_mpnn = config['pipeline'].get('apptainer_ligandmpnn', 'ligandmpnn.sif')
    
    mpnn_fixed_pos = os.path.join(args.analysis_dir, 'mpnn_fixed_positions_B.json')

    # Load dynamic cropped residue range for Chain B from Module 1 index mapping
    mapping_file = os.path.join(args.analysis_dir, 'index_mapping.json')
    crop_min_b = 1
    crop_max_b = 88
    if os.path.exists(mapping_file):
        with open(mapping_file, 'r') as f:
            mapping_data = json.load(f)
            crop_min_b = mapping_data.get('crop_min_B', 1)
            crop_max_b = mapping_data.get('crop_max_B', 88)

    print(f"Module 4: Simulating Rescue Design B' (Crop B range: {crop_min_b}-{crop_max_b})...")

    from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Superimposer

    # 1. Build cropped complex [A' + B_WT(cropped)] for RFD3
    parser = PDBParser(QUIET=True)
    struct_A = parser.get_structure("A_prime", in_pdb)
    wt_pdb = config['pipeline'].get('input_pdb', 'data/inputs/complex_S1_S2.pdb')
    struct_wt = parser.get_structure("WT", wt_pdb)

    chain_A_id = config['pipeline'].get('chain_A', 'A')
    chain_B_id = config['pipeline'].get('chain_B', 'B')

    # Superimpose designed Chain A back onto the WT Chain A spatial frame
    wt_a_ca = [r['CA'] for r in struct_wt[0][chain_A_id] if r.id[0] == ' ' and 'CA' in r]
    cand_a_ca = [r['CA'] for r in struct_A[0][chain_A_id] if r.id[0] == ' ' and 'CA' in r]
    if wt_a_ca and cand_a_ca:
        min_len = min(len(wt_a_ca), len(cand_a_ca))
        sup = Superimposer()
        sup.set_atoms(wt_a_ca[:min_len], cand_a_ca[:min_len])
        sup.apply(struct_A.get_atoms())
        print(f"Module 4: Superimposed designed A' onto WT coordinate frame (RMSD = {sup.rms:.3f} Å)")

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

    complex_pdb = os.path.join(args.out_dir, "complex_A_prime_B_wt.pdb")
    io = PDBIO()
    io.set_structure(complex_struct)
    io.save(complex_pdb)
    print(f"Module 4: Assembled cropped complex [A' + B_WT({crop_min_b}-{crop_max_b})] -> {complex_pdb}")

    # Check if RFD3 co-adaptation is enabled
    lmpnn_cfg = config.get('ligandmpnn', {})
    coadapt_enabled = lmpnn_cfg.get('rescue_coadaptation', True)
    rfd_input_pdb = complex_pdb

    if coadapt_enabled and execution_mode != 'mock':
        print("\nModule 4 (Step 1/2): Running RFD3 Backbone Co-adaptation on B around fixed A'...")
        rfd_rescue_dir = os.path.join(args.out_dir, "rfd3_coadapt_out")
        if os.path.exists(rfd_rescue_dir):
            import shutil
            shutil.rmtree(rfd_rescue_dir)
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
        print(f"  Running RFD3 Co-adaptation -> {' '.join(rfd_cmd)}")
        env = os.environ.copy()
        env['PYTHONUNBUFFERED'] = '1'
        subprocess.run(rfd_cmd, check=True, env=env)

        import glob
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
            print(f"  RFD3 Co-adaptation successful -> {rfd_input_pdb}")

    # Reconstruct full-length Chain B (align WT Chain B to the co-adapted cropped frame, then attach tail)
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
                print(f"  Superimposed WT Chain B tail onto co-adapted framework (RMSD = {sup_b.rms:.3f} Å)")
            
            if chain_B_id in m_full:
                for r in struct_wt[0][chain_B_id]:
                    if r.id[0] == ' ' and r.id[1] > crop_max_b:
                        m_full[chain_B_id].add(r.copy())
            assembled_pdb = os.path.join(args.out_dir, "complex_assembled_full.pdb")
            io_full = PDBIO()
            io_full.set_structure(st_full)
            io_full.save(assembled_pdb)
            rfd_input_pdb = assembled_pdb
            print(f"  Reconstructed seamless full-length complex -> {assembled_pdb}")

    # Step 2: LigandMPNN Sequence Design with Complementary Electrostatic & Steric Rescue Bias
    print("\nModule 4 (Step 2/2): Running LigandMPNN for complementary rescue design (T = 0.10)...")
    mpnn_fixed_pos = os.path.join(args.analysis_dir, 'mpnn_fixed_positions_B.json')
    with open(mpnn_fixed_pos) as f:
        old_fixed = json.load(f)
    fixed_list = []
    for chain, res_list in old_fixed.items():
        fixed_list.extend([f"{chain}{res}" for res in res_list])
    fixed_str = ",".join(fixed_list)

    # Compute Complementary Rescue Bias for B' mutable residues against mutated A'
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

    parser_res = PDBParser(QUIET=True)
    struct_complex = parser_res.get_structure("coadapted_complex", rfd_input_pdb)
    model_c = struct_complex[0]
    chain_A_res = [r for r in model_c[chain_A_id] if r.id[0] == ' ']
    chain_B_res = [r for r in model_c[chain_B_id] if r.id[0] == ' ']

    fixed_b_set = set(old_fixed.get(chain_B_id, []))
    rescue_bias = {}

    for res_B in chain_B_res:
        bid = res_B.id[1]
        if bid in fixed_b_set:
            continue
        
        # Find closest interacting partner on mutated Chain A
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
            b_bias = {}

            # 1. Partner on A' is Positively Charged (Arg, Lys, His) -> Complement with Negatives (Asp, Glu)
            if a_resname in {'ARG', 'LYS', 'HIS'}:
                b_bias = {"ASP": 5.0, "GLU": 5.0, "ARG": -4.0, "LYS": -4.0}
            # 2. Partner on A' is Negatively Charged (Asp, Glu) -> Complement with Positives (Arg, Lys)
            elif a_resname in {'ASP', 'GLU'}:
                b_bias = {"ARG": 5.0, "LYS": 5.0, "ASP": -4.0, "GLU": -4.0}
            # 3. Partner on A' is Aromatic (Tyr, Trp, Phe) -> Complement with Hydrophobic/Polar Pocket
            elif a_resname in {'PHE', 'TYR', 'TRP'}:
                b_bias = {"LEU": 3.0, "ILE": 3.0, "VAL": 3.0, "PHE": 2.0, "TYR": 2.0}
            # 4. Partner on A' is Polar (Asn, Gln, Ser, Thr) -> Complement with H-Bond Networks
            elif a_resname in {'ASN', 'GLN', 'SER', 'THR'}:
                b_bias = {"GLN": 3.0, "ASN": 3.0, "SER": 2.0, "THR": 2.0, "ARG": 2.0, "GLU": 2.0}
            else:
                b_bias = {"LEU": 2.0, "VAL": 2.0, "ALA": 1.5, "ILE": 2.0}

            rescue_bias[b_key] = b_bias

    print(f"  Calculated Complementary Rescue Bias for {len(rescue_bias)} interface residues on Chain {chain_B_id}.")

    if execution_mode == 'mock':
        mpnn_cmd = [
            sys.executable, "data/mock_tools/mock_ligandmpnn.py",
            "--pdb_path", rfd_input_pdb,
            "--out_folder", args.out_dir,
            "--fixed_residues", fixed_str.strip()
        ]
        print(f"  Running Mock LigandMPNN (B'): {' '.join(mpnn_cmd)}")
        subprocess.run(mpnn_cmd, check=True)
    else:
        mpnn_bin = config['pipeline'].get('local_ligandmpnn', 'OrthoIntRob/bin/mpnn') if execution_mode == 'local' else config['pipeline']['colab_ligandmpnn']
        mpnn_out_dir = os.path.join(args.out_dir, "mpnn_out")
        os.makedirs(mpnn_out_dir, exist_ok=True)
        
        default_ckpt = "OrthoIntRob/ligandmpnn/model_params/ligandmpnn_v_32_010_25.pt" if execution_mode == 'local' else "/content/drive/MyDrive/OrthoInterface_Data/FoundryModels/LigandMPNN/model_params/ligandmpnn_v_32_010_25.pt"
        checkpoint_path = lmpnn_cfg.get('checkpoint_path', default_ckpt)
        model_type = lmpnn_cfg.get('model_type', "ligand_mpnn")
        is_legacy = str(lmpnn_cfg.get('is_legacy_weights', "True"))
        temp_rescue = str(lmpnn_cfg.get('temperature_rescue', 0.10))

        mpnn_cmd = mpnn_bin.split() + [
            "--structure_path", rfd_input_pdb.replace('\\', '/'),
            "--out_directory", mpnn_out_dir.replace('\\', '/'),
            "--model_type", model_type,
            "--checkpoint_path", checkpoint_path.replace('\\', '/'),
            "--is_legacy_weights", is_legacy,
            "--batch_size", "1",
            "--number_of_batches", "1",
            "--temperature", temp_rescue,
            "--fixed_residues", fixed_str,
            "--bias_per_residue", json.dumps(rescue_bias),
            "--write_structures", "True"
        ]

        print(f"  Running LigandMPNN (B') -> {' '.join(mpnn_cmd)}")
        env = os.environ.copy()
        subprocess.run(mpnn_cmd, check=True, env=env)

    import glob
    import shutil
    
    out_search_dir = mpnn_out_dir if execution_mode in ['colab', 'local'] else args.out_dir

    mpnn_outputs = [f for f in glob.glob(os.path.join(out_search_dir, "**", "*.pdb"), recursive=True) if not f.endswith("B_prime_rfd_mpnn.pdb")]
    mpnn_cifs = glob.glob(os.path.join(out_search_dir, "**", "*.cif.gz"), recursive=True) + glob.glob(os.path.join(out_search_dir, "**", "*.cif"), recursive=True)
    
    mpnn_outputs = [f for f in mpnn_outputs if os.path.abspath(f) != os.path.abspath(rfd_input_pdb)]
    mpnn_cifs = [f for f in mpnn_cifs if os.path.abspath(f) != os.path.abspath(rfd_input_pdb)]

    packed_outputs = [f for f in mpnn_outputs if "packed" in f]
    if packed_outputs:
        mpnn_outputs = packed_outputs

    final_pdb_path = os.path.join(args.out_dir, "B_prime_rfd_mpnn.pdb")
    if mpnn_outputs:
        best_pdb = mpnn_outputs[0]
        shutil.copy(best_pdb, final_pdb_path)
    elif mpnn_cifs:
        best_cif = mpnn_cifs[0]
        import gzip
        from Bio.PDB import MMCIFParser
        
        parser = MMCIFParser(QUIET=True)
        if best_cif.endswith('.gz'):
            with gzip.open(best_cif, 'rt') as f:
                structure = parser.get_structure("struct", f)
        else:
            with open(best_cif, 'rt') as f:
                structure = parser.get_structure("struct", f)
                
        io = PDBIO()
        io.set_structure(structure)
        io.save(final_pdb_path)
    else:
        raise FileNotFoundError("LigandMPNN did not produce any backbones!")

    # Clean nan sidechains
    with open(final_pdb_path, 'r') as f:
        lines = f.readlines()
    with open(final_pdb_path, 'w') as f:
        for line in lines:
            if 'nan' not in line:
                f.write(line)

    print(f"\nModule 4 Complete: Rescue B' generated successfully -> {final_pdb_path}")

if __name__ == "__main__":
    main()
