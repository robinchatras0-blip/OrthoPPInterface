import argparse
import yaml
import subprocess
import sys
import os
import json

def main():
    parser = argparse.ArgumentParser(description="Module 2: A' Generation (Rupture Design)")
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--analysis_dir', default='results/01_analysis')
    parser.add_argument('--out_dir', default='results/02_rupture_design')
    args = parser.parse_args()

    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    os.makedirs(args.out_dir, exist_ok=True)

    execution_mode = config['pipeline'].get('execution_mode', 'mock')
    input_pdb = config['pipeline']['input_pdb']
    apptainer_rfd = config['pipeline'].get('apptainer_rfdiffusion', 'rfdiffusion.sif')
    apptainer_mpnn = config['pipeline'].get('apptainer_ligandmpnn', 'ligandmpnn.sif')

    rfd_inputs_json = os.path.join(args.analysis_dir, 'inputs.json').replace('\\', '/')
    mpnn_fixed_pos = os.path.join(args.analysis_dir, 'mpnn_fixed_positions.json').replace('\\', '/')
    mpnn_bias = os.path.join(args.analysis_dir, 'mpnn_bias.json').replace('\\', '/')

    # 1. Execute RFD3 (Foundry)
    rfd_output_dir = os.path.join(args.out_dir, 'A_prime_rfd_out').replace('\\', '/')
    os.makedirs(rfd_output_dir, exist_ok=True)

    if execution_mode == 'mock':
        rfd_cmd = [
            sys.executable, "data/mock_tools/mock_rfdiffusion.py",
            "--inference.input_pdb", input_pdb,
            "--inference.output_prefix", os.path.join(rfd_output_dir, "design").replace('\\', '/'),
            "--contigmap.contigs", rfd_inputs_json
        ]
    elif execution_mode in ['colab', 'local']:
        rfd_bin = config['pipeline'].get('local_rfdiffusion', 'OrthoIntRob/bin/rfd3') if execution_mode == 'local' else config['pipeline']['colab_rfdiffusion']
        n_batches = config['pipeline'].get('foundry_n_batches', 1)
        diff_batch = config['pipeline'].get('foundry_diffusion_batch_size', 1)
        rfd_cmd = rfd_bin.split() + [
            f"inputs={rfd_inputs_json}",
            f"out_dir={rfd_output_dir}",
            f"n_batches={n_batches}",
            f"diffusion_batch_size={diff_batch}"
        ]
    else:
        rfd_cmd = [
            "apptainer", "exec", "--nv", apptainer_rfd,
            "foundry", "run", "rfd3", "design",
            f"inputs={rfd_inputs_json}",
            f"out_dir={rfd_output_dir}"
        ]

    print(f"Module 2: Running RFD3 -> {' '.join(rfd_cmd)}")
    import glob
    import shutil

    if os.path.exists(rfd_output_dir):
        shutil.rmtree(rfd_output_dir)
    os.makedirs(rfd_output_dir, exist_ok=True)

    subprocess.run(rfd_cmd, check=True)
    rfd_pdbs = sorted(glob.glob(os.path.join(rfd_output_dir, "*.pdb")))
    rfd_cifs = sorted(glob.glob(os.path.join(rfd_output_dir, "*.cif.gz")) + glob.glob(os.path.join(rfd_output_dir, "*.cif")))

    if not rfd_pdbs and rfd_cifs:
        import gzip
        from Bio.PDB import MMCIFParser, PDBIO
        for cif_path in rfd_cifs:
            pdb_path = cif_path.replace(".cif.gz", ".pdb").replace(".cif", ".pdb")
            if cif_path.endswith('.gz'):
                with gzip.open(cif_path, 'rt') as f:
                    parser = MMCIFParser(QUIET=True)
                    structure = parser.get_structure("struct", f)
            else:
                with open(cif_path, 'rt') as f:
                    parser = MMCIFParser(QUIET=True)
                    structure = parser.get_structure("struct", f)
            io = PDBIO()
            io.set_structure(structure)
            io.save(pdb_path)
            rfd_pdbs.append(pdb_path)
    elif rfd_pdbs:
        print(f"Found {len(rfd_pdbs)} RFD3 backbones.")

    # 2. Execute LigandMPNN for each candidate backbone
    if not rfd_pdbs:
        raise FileNotFoundError("RFD3 did not produce any backbones!")

    with open(mpnn_fixed_pos) as f:
        old_fixed = json.load(f)
    fixed_str = ""
    fixed_list = []
    for chain, res_list in old_fixed.items():
        fixed_list.extend([f"{chain}{res}" for res in res_list])
    fixed_str = ",".join(fixed_list)

    with open(mpnn_bias) as f:
        old_bias = json.load(f)
    new_bias = {}
    for chain, res_dict in old_bias.items():
        for res, bias_dict in res_dict.items():
            new_bias[f"{chain}{res}"] = bias_dict
    
    mpnn_bias_path = os.path.join(args.out_dir, "mpnn_bias_ligandmpnn.json").replace('\\', '/')
    with open(mpnn_bias_path, "w") as f:
        json.dump(new_bias, f)

    from Bio.PDB import PDBParser

    for idx, rfd_pdb_out in enumerate(sorted(rfd_pdbs)):
        print(f"\n--- Processing Candidate {idx+1}/{len(rfd_pdbs)}: {os.path.basename(rfd_pdb_out)} ---")
        curr_mpnn_out_dir = os.path.join(args.out_dir, f"mpnn_out_cand_{idx}").replace('\\', '/')
        os.makedirs(curr_mpnn_out_dir, exist_ok=True)

        # Filter fixed residues to only those actually present in this PDB
        parser_chk = PDBParser(QUIET=True)
        st_chk = parser_chk.get_structure("chk", rfd_pdb_out)
        existing_res = set(f"{ch.id}{r.id[1]}" for m in st_chk for ch in m for r in ch if r.id[0] == ' ')
        curr_fixed_list = [f for f in fixed_list if f in existing_res]
        curr_fixed_str = ",".join(curr_fixed_list)

        if execution_mode == 'mock':
            mpnn_cmd = [
                sys.executable, "data/mock_tools/mock_ligandmpnn.py",
                "--pdb_path", rfd_pdb_out.replace('\\', '/'),
                "--out_folder", curr_mpnn_out_dir,
                "--fixed_residues", curr_fixed_str,
                "--bias_AA_per_residue", str(mpnn_bias_path)
            ]
        elif execution_mode in ['colab', 'local']:
            mpnn_bin = config['pipeline'].get('local_ligandmpnn', 'OrthoIntRob/bin/mpnn') if execution_mode == 'local' else config['pipeline']['colab_ligandmpnn']
            lmpnn_cfg = config.get('ligandmpnn', {})
            default_ckpt = "OrthoIntRob/ligandmpnn/model_params/ligandmpnn_v_32_010_25.pt" if execution_mode == 'local' else "/content/drive/MyDrive/OrthoPPInterface_Data/FoundryModels/LigandMPNN/model_params/ligandmpnn_v_32_010_25.pt"
            checkpoint_path = lmpnn_cfg.get('checkpoint_path', default_ckpt)
            model_type = lmpnn_cfg.get('model_type', "ligand_mpnn")
            is_legacy = str(lmpnn_cfg.get('is_legacy_weights', "True"))

            temp_rupture = str(lmpnn_cfg.get('temperature_rupture', 0.20))

            mpnn_cmd = mpnn_bin.split() + [
                "--structure_path", rfd_pdb_out.replace('\\', '/'),
                "--out_directory", curr_mpnn_out_dir,
                "--model_type", model_type,
                "--checkpoint_path", checkpoint_path.replace('\\', '/'),
                "--is_legacy_weights", is_legacy,
                "--batch_size", "1",
                "--number_of_batches", "1",
                "--temperature", temp_rupture,
                "--fixed_residues", curr_fixed_str,
                "--bias_per_residue", json.dumps(new_bias),
                "--write_structures", "True"
            ]
        else:
            mpnn_cmd = [
                "apptainer", "exec", "--nv", apptainer_mpnn,
                "python", "/app/run.py",
                "--pdb_path", rfd_pdb_out.replace('\\', '/'),
                "--out_folder", curr_mpnn_out_dir,
                "--fixed_residues", curr_fixed_str.strip(),
                "--bias_AA_per_residue", str(mpnn_bias_path)
            ]

        print(f"Module 2: Running LigandMPNN -> {' '.join(mpnn_cmd)}")
        env = os.environ.copy()
        subprocess.run(mpnn_cmd, check=True, env=env)

        mpnn_outputs = [f for f in glob.glob(os.path.join(curr_mpnn_out_dir, "**", "*.pdb"), recursive=True)]
        mpnn_outputs = [f for f in mpnn_outputs if os.path.abspath(f) != os.path.abspath(rfd_pdb_out)]
        packed_outputs = [f for f in mpnn_outputs if "packed" in f]
        if packed_outputs:
            mpnn_outputs = packed_outputs

        final_cand_pdb = os.path.join(args.out_dir, f"A_prime_candidate_{idx:02d}.pdb")
        if mpnn_outputs:
            shutil.copy(mpnn_outputs[0], final_cand_pdb)
        else:
            shutil.copy(rfd_pdb_out, final_cand_pdb)

        # Clean nan sidechains
        with open(final_cand_pdb, 'r') as f:
            lines = f.readlines()
        with open(final_cand_pdb, 'w') as f:
            for line in lines:
                if 'nan' not in line:
                    f.write(line)

        # Default compatibility
        if idx == 0:
            shutil.copy(final_cand_pdb, os.path.join(args.out_dir, "A_prime_rfd_mpnn.pdb"))

    print(f"Module 2 Complete: {len(rfd_pdbs)} A' Rupture Designs generated successfully.")

if __name__ == "__main__":
    main()
