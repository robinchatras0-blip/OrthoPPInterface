import argparse
import glob
import json
import os
import shutil
import subprocess
import sys

from Bio.PDB import PDBParser

sys.path.append(os.path.dirname(__file__))
from design_utils import collect_pdbs, load_config, save_structure, std_residues, strip_nan_lines
from energy_engine import enabled as energy_enabled, pack_structures


def main():
    parser = argparse.ArgumentParser(description="Module 2: A' Generation (Rupture Design)")
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--analysis_dir', default='results/01_analysis')
    parser.add_argument('--out_dir', default='results/02_rupture_design')
    args = parser.parse_args()

    config = load_config(args.config)
    pcfg, lcfg = config['pipeline'], config.get('ligandmpnn', {})
    chain_A_id = pcfg.get('chain_A', 'A')
    os.makedirs(args.out_dir, exist_ok=True)

    # 1. RFD3 (Foundry): sample backbones with the interface segments of A regenerated de novo
    rfd_output_dir = os.path.join(args.out_dir, 'A_prime_rfd_out').replace('\\', '/')
    if os.path.exists(rfd_output_dir):
        shutil.rmtree(rfd_output_dir)
    os.makedirs(rfd_output_dir)
    rfd_cmd = pcfg.get('local_rfdiffusion', 'OrthoIntRob/bin/rfd3').split() + [
        f"inputs={os.path.join(args.analysis_dir, 'inputs.json').replace(chr(92), '/')}",
        f"out_dir={rfd_output_dir}",
        f"n_batches={pcfg.get('foundry_n_batches', 1)}",
        f"diffusion_batch_size={pcfg.get('foundry_diffusion_batch_size', 1)}",
    ]
    print(f"Module 2: Running RFD3 -> {' '.join(rfd_cmd)}")
    subprocess.run(rfd_cmd, check=True)
    rfd_pdbs = collect_pdbs(rfd_output_dir)
    if not rfd_pdbs:
        raise FileNotFoundError("RFD3 did not produce any backbones!")
    print(f"Found {len(rfd_pdbs)} RFD3 backbones.")

    # 2. LigandMPNN on every backbone: chain B and the framework of A are fixed, rupture bias on the interface
    with open(os.path.join(args.analysis_dir, 'mpnn_fixed_positions.json')) as f:
        fixed_list = [f"{chain}{res}" for chain, res_list in json.load(f).items() for res in res_list]
    with open(os.path.join(args.analysis_dir, 'mpnn_bias.json')) as f:
        bias = {f"{chain}{res}": bias_dict for chain, res_dict in json.load(f).items() for res, bias_dict in res_dict.items()}

    mpnn_bin = pcfg.get('local_ligandmpnn', 'OrthoIntRob/bin/mpnn')
    checkpoint = lcfg.get('checkpoint_path', 'OrthoIntRob/ligandmpnn/model_params/ligandmpnn_v_32_010_25.pt')

    for idx, rfd_pdb in enumerate(rfd_pdbs):
        print(f"\n--- Processing candidate {idx + 1}/{len(rfd_pdbs)}: {os.path.basename(rfd_pdb)} ---")
        mpnn_out = os.path.join(args.out_dir, f"mpnn_out_cand_{idx}").replace('\\', '/')
        os.makedirs(mpnn_out, exist_ok=True)

        structure = PDBParser(QUIET=True).get_structure("rfd", rfd_pdb)
        mpnn_input = rfd_pdb
        if lcfg.get('rupture_apo_design', True):
            # Design A' alone. With B_WT in the structure LigandMPNN packs A' against B_WT (a positive design for the
            # wild-type partner), which defeats the rupture (measured: iPTM(A'.B_WT) 0.38 vs 0.23 without B).
            for ch in [c for c in structure[0] if c.id != chain_A_id]:
                structure[0].detach_child(ch.id)
            mpnn_input = os.path.join(args.out_dir, 'mpnn_in', f"cand_{idx}.pdb")
            os.makedirs(os.path.dirname(mpnn_input), exist_ok=True)
            save_structure(structure, mpnn_input)
        existing = {f"{ch.id}{r.id[1]}" for ch in structure[0] for r in std_residues(ch)}
        mpnn_cmd = mpnn_bin.split() + [
            "--structure_path", mpnn_input.replace('\\', '/'),
            "--out_directory", mpnn_out,
            "--model_type", lcfg.get('model_type', 'ligand_mpnn'),
            "--checkpoint_path", checkpoint.replace('\\', '/'),
            "--is_legacy_weights", str(lcfg.get('is_legacy_weights', 'True')),
            "--batch_size", "1",
            "--number_of_batches", "1",
            "--temperature", str(lcfg.get('temperature_rupture', 0.20)),
            "--fixed_residues", ",".join(f for f in fixed_list if f in existing),
            "--bias_per_residue", json.dumps(bias),
            "--write_structures", "True",
        ]
        print(f"Module 2: Running LigandMPNN -> {' '.join(mpnn_cmd[:6])} ...")
        subprocess.run(mpnn_cmd, check=True)

        outputs = [f for f in collect_pdbs(mpnn_out) if os.path.abspath(f) != os.path.abspath(rfd_pdb)]
        outputs = [f for f in outputs if "packed" in f] or outputs
        if not outputs:
            # never fall back to the raw RFD3 backbone: it would silently discard the MPNN design
            raise RuntimeError(f"LigandMPNN wrote no structure for {rfd_pdb} (see {mpnn_out})")
        final_pdb = os.path.join(args.out_dir, f"A_prime_candidate_{idx:02d}.pdb")
        shutil.copy(outputs[0], final_pdb)
        strip_nan_lines(final_pdb)

    # LigandMPNN writes no side chain for the residues it redesigns: rebuild them (PyRosetta) so that every
    # downstream geometric measure and every visualisation sees complete residues.
    if energy_enabled(config) and config['energy'].get('pack_designs', True):
        cands = sorted(glob.glob(os.path.join(args.out_dir, "A_prime_candidate_*.pdb")))
        print(f"Module 2: packing the side chains of {len(cands)} A' designs (PyRosetta)...")
        pack_structures([{"name": os.path.basename(c)[:-4], "pdb": c, "out": c} for c in cands], config,
                        os.path.join(args.out_dir, "pack_status.json"))
    else:
        print("Module 2: WARNING - side-chain packing disabled (energy.enabled): redesigned residues have no side chain.")

    print(f"Module 2 Complete: {len(rfd_pdbs)} A' rupture designs generated.")


if __name__ == "__main__":
    main()
