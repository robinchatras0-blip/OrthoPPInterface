import os
import sys
import json
import glob
import subprocess
from Bio.PDB import PDBParser, MMCIFParser, PDBIO

def to_wsl_path(path):
    """Converts a Windows absolute or relative path to a WSL compatible path."""
    abs_path = os.path.abspath(path).replace('\\', '/')
    if len(abs_path) > 1 and abs_path[1] == ':':
        drive = abs_path[0].lower()
        return f"/mnt/{drive}{abs_path[2:]}"
    return abs_path

def convert_cif_to_pdb(cif_path, pdb_path):
    """Converts a CIF format file to standard PDB format using BioPython."""
    try:
        cif_parser = MMCIFParser(QUIET=True)
        struct = cif_parser.get_structure("cif_model", cif_path)
        io = PDBIO()
        io.set_structure(struct)
        io.save(pdb_path)
        return True
    except Exception as e:
        print(f"Warning: Could not convert {cif_path} to PDB: {e}")
        return False


THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def get_chain_sequence(pdb_path, chain_id='A'):
    """Extracts amino acid sequence (1-letter code) of a chain from a PDB file."""
    if not os.path.exists(pdb_path):
        return ""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_path)
    model = structure[0]
    if chain_id not in model:
        chains = list(model.get_chains())
        if not chains:
            return ""
        chain = chains[0]
    else:
        chain = model[chain_id]
    seq = "".join(THREE_TO_ONE.get(res.get_resname(), 'X') for res in chain if res.id[0] == ' ')
    return seq

def extract_all_sequences(pdb_path):
    """Extracts all chain sequences from a PDB file as a dict {chain_id: sequence}."""
    if not os.path.exists(pdb_path):
        return {}
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_path)
    model = structure[0]
    sequences = {}
    for chain in model.get_chains():
        seq = "".join(THREE_TO_ONE.get(res.get_resname(), 'X') for res in chain if res.id[0] == ' ')
        if seq:
            sequences[chain.id] = seq
    return sequences

def build_hybrid_msa(wt_a3m_path, target_sequence, modified_indices, output_a3m_path):
    """
    Parses WT .a3m MSA file, updates query sequence (line 1), and replaces modified
    residue columns with gaps ('-') for all homologous sequences in the alignment.
    """
    if not wt_a3m_path or not os.path.exists(wt_a3m_path):
        with open(output_a3m_path, 'w') as f:
            f.write(f">query\n{target_sequence}\n")
        return

    with open(wt_a3m_path, 'r') as f:
        lines = f.readlines()

    modified_set = set(modified_indices)
    output_lines = []
    is_query = True

    for line in lines:
        line_str = line.strip()
        if not line_str:
            continue
        if line_str.startswith(">"):
            output_lines.append(line)
        else:
            if is_query:
                output_lines.append(target_sequence + "\n")
                is_query = False
            else:
                new_chars = []
                match_col_idx = 0
                for char in line_str:
                    if char.isupper() or char == '-':
                        if match_col_idx in modified_set:
                            new_chars.append('-')
                        else:
                            new_chars.append(char)
                        match_col_idx += 1
                    else:
                        # Lowercase insertion
                        if match_col_idx in modified_set:
                            continue # Omit insertions within modified/de novo positions
                        else:
                            new_chars.append(char)
                output_lines.append("".join(new_chars) + "\n")

    with open(output_a3m_path, 'w') as f:
        f.writelines(output_lines)

def run_rf3_prediction(input_pdb, fasta_sequences, out_dir, config, msa_path=None, execution_mode="local"):
    """
    Runs RoseTTAFold-3 / All-Atom (Foundry RF3) for de novo complex or monomer structure prediction.
    Uses individual monomer MSAs for internal chain folding, bypassing cross-chain covariance bias.
    """
    os.makedirs(out_dir, exist_ok=True)
    tag = os.path.basename(out_dir)

    # Check if prediction is already completed and cached
    summary_files = sorted(glob.glob(os.path.join(out_dir, "**/*_summary_confidences.json"), recursive=True))
    rank_pdb_files = sorted(glob.glob(os.path.join(out_dir, "*rank_001*.pdb"))) + sorted(glob.glob(os.path.join(out_dir, "*_predicted.pdb")))

    if summary_files and rank_pdb_files:
        print(f"  [Foundry RF3] Found cached prediction in {out_dir}, skipping re-computation.")
        with open(summary_files[0], 'r') as f:
            data = json.load(f)
            overall_plddt = float(data.get("overall_plddt", 0.0))
            plddt = overall_plddt * 100.0 if overall_plddt <= 1.0 else float(overall_plddt)
            iptm = float(data.get("iptm", 0.0))
            ptm = float(data.get("ptm", 0.0))
            ranking_score = float(data.get("ranking_score", iptm))
            return {
                "plddt": plddt,
                "iptm": iptm,
                "ptm": ptm,
                "ranking_score": ranking_score
            }

    use_msa = config.get('folding', {}).get('use_msa', True)
    default_msa_B = config.get('pipeline', {}).get('input_msa_B', None)

    # Build components list from fasta_sequences: [('A_prime', seq_A), ('B_prime', seq_B)]
    # Each item can be (name, seq) or (name, seq, msa_path)
    components = []
    chain_ids = ['A', 'B', 'C', 'D']
    for idx, item in enumerate(fasta_sequences):
        label = item[0]
        seq = item[1]
        msa_p = item[2] if len(item) > 2 else None

        if not msa_p and use_msa:
            if idx == 0 and msa_path and os.path.exists(msa_path):
                msa_p = msa_path
            elif idx == 1 and default_msa_B and os.path.exists(default_msa_B):
                msa_p = default_msa_B

        ch = chain_ids[idx] if idx < len(chain_ids) else chr(ord('A') + idx)
        comp = {"seq": seq, "chain_id": ch}
        if use_msa and msa_p and os.path.exists(msa_p) and os.path.getsize(msa_p) > 0:
            comp["msa_path"] = to_wsl_path(msa_p)
        components.append(comp)

    manifest_data = [{
        "name": tag,
        "components": components
    }]

    manifest_path = os.path.join(out_dir, f"{tag}_rf3_input.json")
    with open(manifest_path, 'w') as f:
        json.dump(manifest_data, f, indent=2)

    if execution_mode == 'mock':
        return {
            "plddt": 85.0,
            "iptm": 0.82,
            "ptm": 0.80,
            "ranking_score": 0.82
        }

    rf3_cfg = config.get('folding', {})
    rf3_bin = rf3_cfg.get('rf3_bin', '/home/ommearo/miniforge3/envs/foundry_env/bin/rf3')
    rf3_ckpt = rf3_cfg.get('rf3_ckpt', '/home/ommearo/.foundry/checkpoints/rf3_foundry_01_24_latest.ckpt')

    wsl_manifest = to_wsl_path(manifest_path)
    wsl_out = to_wsl_path(out_dir)

    summary_files = sorted(glob.glob(os.path.join(out_dir, "**/*_summary_confidences.json"), recursive=True))
    if summary_files:
        print(f"  [Foundry RF3 Resume Cache] Found existing RF3 prediction in {out_dir}, skipping invocation.")
    else:
        cmd = f"wsl -d Ubuntu -- {rf3_bin} fold inputs={wsl_manifest} out_dir={wsl_out} ckpt_path={rf3_ckpt}"
        print(f"  [Foundry RF3 GPU Invoc] {cmd}")

        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if res.returncode != 0:
            print(f"  [Foundry RF3 Error]: {res.stderr[:400] if res.stderr else res.stdout[:400]}")
            raise RuntimeError(f"RF3 execution failed with code {res.returncode}: {res.stderr[:300]}")

    # Parse output
    summary_files = sorted(glob.glob(os.path.join(out_dir, "**/*_summary_confidences.json"), recursive=True))
    if not summary_files:
        raise RuntimeError(f"RF3 output _summary_confidences.json not found in {out_dir}")

    with open(summary_files[0], 'r') as f:
        data = json.load(f)
        overall_plddt = float(data.get("overall_plddt", 0.0))
        plddt = overall_plddt * 100.0 if overall_plddt <= 1.0 else float(overall_plddt)
        iptm = float(data.get("iptm", 0.0))
        ptm = float(data.get("ptm", 0.0))
        ranking_score = float(data.get("ranking_score", iptm))

    # Find top predicted CIF model and convert to rank_001_predicted.pdb
    cif_files = sorted(glob.glob(os.path.join(out_dir, f"**/{tag}_model.cif"), recursive=True))
    if not cif_files:
        cif_files = sorted(glob.glob(os.path.join(out_dir, "**/*_model.cif"), recursive=True))
    if not cif_files:
        cif_files = sorted(glob.glob(os.path.join(out_dir, "**/*.cif"), recursive=True))

    if cif_files:
        top_cif = cif_files[0]
        out_pdb = os.path.join(out_dir, f"{tag}_unrelaxed_rank_001.pdb")
        convert_cif_to_pdb(top_cif, out_pdb)
        print(f"  [Foundry RF3] Converted {os.path.basename(top_cif)} -> {os.path.basename(out_pdb)}")

    return {
        "plddt": plddt,
        "iptm": iptm,
        "ptm": ptm,
        "ranking_score": ranking_score
    }

def predict_structure(input_pdb, fasta_sequences, out_dir, msa_path, config):
    """
    Unified folding prediction dispatching to RoseTTAFold-3 (RF3).
    fasta_sequences: list of tuples [('A_prime', 'SEQUENCE_A'), ('B', 'SEQUENCE_B')]
    """
    execution_mode = config.get('pipeline', {}).get('execution_mode', 'local')
    return run_rf3_prediction(input_pdb, fasta_sequences, out_dir, config, msa_path=msa_path, execution_mode=execution_mode)

def calculate_ca_rmsd(ref_pdb, pred_pdb_or_dir, chain_ref='A', chain_pred=None, subset_res_ids=None):
    """
    Computes C-alpha backbone RMSD between reference structure and predicted structure.
    pred_pdb_or_dir can be a direct PDB path or an output directory containing *rank_001*.pdb.
    Optionally filters by specific residue sequential IDs (e.g. fixed framework residues).
    """
    from Bio.PDB import PDBParser, Superimposer
    
    pred_pdb = pred_pdb_or_dir
    if os.path.isdir(pred_pdb_or_dir):
        ranked_pdbs = sorted(glob.glob(os.path.join(pred_pdb_or_dir, "*unrelaxed_rank_001*.pdb"))) + sorted(glob.glob(os.path.join(pred_pdb_or_dir, "*rank_001*.pdb")))
        if not ranked_pdbs:
            ranked_pdbs = sorted(glob.glob(os.path.join(pred_pdb_or_dir, "*.pdb")))
        if ranked_pdbs:
            pred_pdb = ranked_pdbs[0]
        else:
            return 0.0

    if not os.path.exists(ref_pdb) or not os.path.exists(pred_pdb):
        return 0.0

    parser = PDBParser(QUIET=True)
    struct_ref = parser.get_structure("ref", ref_pdb)
    struct_pred = parser.get_structure("pred", pred_pdb)

    model_ref = struct_ref[0]
    model_pred = struct_pred[0]

    # Reference chain
    if chain_ref in model_ref:
        ch_ref = model_ref[chain_ref]
    else:
        ch_ref = list(model_ref.get_chains())[0]

    # Pred chain
    if chain_pred and chain_pred in model_pred:
        ch_pred = model_pred[chain_pred]
    else:
        ch_pred = list(model_pred.get_chains())[0]

    res_ref = [r for r in ch_ref if r.id[0] == ' ']
    res_pred = [r for r in ch_pred if r.id[0] == ' ']

    subset_set = set(subset_res_ids) if subset_res_ids is not None else None

    ca_ref = []
    ca_pred = []

    for r_ref, r_pred in zip(res_ref, res_pred):
        if subset_set is not None and r_ref.id[1] not in subset_set:
            continue
        if 'CA' in r_ref and 'CA' in r_pred:
            ca_ref.append(r_ref['CA'])
            ca_pred.append(r_pred['CA'])

    if not ca_ref or len(ca_ref) < 3:
        return 0.0

    sup = Superimposer()
    sup.set_atoms(ca_ref, ca_pred)
    return float(sup.rms)
