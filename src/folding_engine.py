import os
import sys
import json
import glob
import subprocess
from Bio.PDB import PDBParser

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

def build_unpaired_complex_a3m(a3m_A_path, a3m_B_path, seq_A, seq_B, modified_indices_A, output_a3m_path, modified_indices_B=None, max_seqs=500):
    """
    Builds a multi-chain unpaired ColabFold .a3m file from local A and B MSAs.
    - Chains A and B have modified positions masked with gaps ('-').
    - Sequences are placed in block-diagonal (unpaired) format.
    - ColabFold uses this local MSA directly and completely skips MMseqs2 server lookup.
    """
    len_A = len(seq_A)
    len_B = len(seq_B)
    modified_set_A = set(modified_indices_A) if modified_indices_A else set()
    modified_set_B = set(modified_indices_B) if modified_indices_B else set()
    out_lines = []

    # ColabFold multi-chain header
    out_lines.append(f"#{len_A},{len_B}\t1,1\n")
    out_lines.append(">101\t102\n")
    out_lines.append(f"{seq_A}\t{seq_B}\n")

    # Process Chain A sequences (block 1: A aligned, B padded with gaps)
    if a3m_A_path and os.path.exists(a3m_A_path):
        with open(a3m_A_path, 'r') as f:
            lines_A = f.readlines()
        count_A = 0
        is_query = True
        curr_hdr = None
        for line in lines_A:
            l = line.strip()
            if not l:
                continue
            if l.startswith(">"):
                curr_hdr = l
            else:
                if is_query:
                    is_query = False
                    continue
                if count_A >= max_seqs:
                    break
                # Apply column masking for modified positions
                new_chars = []
                match_col_idx = 0
                for char in l:
                    if char.isupper() or char == '-':
                        if match_col_idx in modified_set_A:
                            new_chars.append('-')
                        else:
                            new_chars.append(char)
                        match_col_idx += 1
                    else:
                        if match_col_idx in modified_set_A:
                            continue
                        else:
                            new_chars.append(char)
                aln_A = "".join(new_chars)
                out_lines.append(f"{curr_hdr}\t102\n")
                out_lines.append(f"{aln_A}\t" + "-" * len_B + "\n")
                count_A += 1

    # Process Chain B sequences (block 2: A padded with gaps, B aligned)
    if a3m_B_path and os.path.exists(a3m_B_path):
        with open(a3m_B_path, 'r') as f:
            lines_B = f.readlines()
        count_B = 0
        is_query = True
        curr_hdr = None
        for line in lines_B:
            l = line.strip()
            if not l:
                continue
            if l.startswith(">"):
                curr_hdr = l
            else:
                if is_query:
                    is_query = False
                    continue
                if count_B >= max_seqs:
                    break
                new_chars_B = []
                match_col_idx = 0
                for char in l:
                    if char.isupper() or char == '-':
                        if match_col_idx in modified_set_B:
                            new_chars_B.append('-')
                        else:
                            new_chars_B.append(char)
                        match_col_idx += 1
                    else:
                        if match_col_idx in modified_set_B:
                            continue
                        else:
                            new_chars_B.append(char)
                aln_B = "".join(new_chars_B)
                out_lines.append(f">101\t{curr_hdr[1:]}\n")
                out_lines.append("-" * len_A + f"\t{aln_B}\n")
                count_B += 1

    with open(output_a3m_path, 'w') as f:
        f.writelines(out_lines)

def run_colabfold_prediction(fasta_path, out_dir, config, execution_mode):
    """Runs ColabFold (AlphaFold 2) and extracts pLDDT, iPTM, and ranking score."""
    os.makedirs(out_dir, exist_ok=True)
    fasta_path = fasta_path.replace('\\', '/')
    out_dir = out_dir.replace('\\', '/')
    colabfold_cfg = config.get('folding', {})
    colabfold_cmd = colabfold_cfg.get('colabfold_cmd', 'colabfold_batch')
    num_recycle = str(colabfold_cfg.get('colabfold_num_recycle', 3))
    model_type = str(colabfold_cfg.get('colabfold_model_type', 'alphafold2_multimer_v3'))
    pair_mode = str(colabfold_cfg.get('colabfold_pair_mode', 'paired'))

    if execution_mode == 'mock':
        mock_script = "data/mock_tools/mock_alphafold3.py"
        cmd = [sys.executable, mock_script, "--input", fasta_path, "--out_dir", out_dir]
    else:
        cmd = colabfold_cmd.split() + [
            fasta_path, out_dir,
            "--num-recycle", num_recycle,
            "--model-type", model_type,
            "--pair-mode", pair_mode
        ]

    print(f"  [AF2/ColabFold] Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    # Parse ColabFold output JSONs
    # In mock mode: *_metrics.json
    mock_metrics_files = glob.glob(os.path.join(out_dir, "*_metrics.json"))
    if mock_metrics_files and execution_mode == 'mock':
        with open(mock_metrics_files[0], 'r') as f:
            data = json.load(f)
            return {
                "plddt": data.get("plddt", 85.0),
                "iptm": data.get("iptm", 0.8),
                "ptm": data.get("ptm", 0.8),
                "ranking_score": data.get("ranking_score", 0.8)
            }

    # In real ColabFold: look for *scores_rank_001*.json or *scores*.json
    score_files = sorted(glob.glob(os.path.join(out_dir, "*scores*.json")))
    if score_files:
        with open(score_files[0], 'r') as f:
            data = json.load(f)
            plddt_list = data.get("plddt", [85.0])
            mean_plddt = sum(plddt_list) / len(plddt_list) if isinstance(plddt_list, list) and plddt_list else float(plddt_list)
            iptm = float(data.get("iptm", data.get("ptm", 0.0)))
            ptm = float(data.get("ptm", 0.0))
            ranking_score = float(data.get("ranking_confidence", data.get("multimer", iptm)))
            return {
                "plddt": mean_plddt,
                "iptm": iptm,
                "ptm": ptm,
                "ranking_score": ranking_score
            }

    return {"plddt": 85.0, "iptm": 0.5, "ptm": 0.5, "ranking_score": 0.5}

def run_af3_prediction(pdb_or_json, out_dir, msa_path, config, execution_mode):
    """Runs AlphaFold 3 and extracts pLDDT, iPTM, and ranking score."""
    os.makedirs(out_dir, exist_ok=True)
    pdb_or_json = pdb_or_json.replace('\\', '/')
    out_dir = out_dir.replace('\\', '/')
    if msa_path:
        msa_path = msa_path.replace('\\', '/')
    af3_cfg = config.get('folding', {})
    af3_cmd = af3_cfg.get('af3_cmd', config.get('pipeline', {}).get('af3_cmd', 'alphafold3'))

    if execution_mode == 'mock':
        mock_script = "data/mock_tools/mock_alphafold3.py"
        cmd = [sys.executable, mock_script, "--predict", pdb_or_json, "--out_dir", out_dir]
        if msa_path:
            cmd.extend(["--msa", msa_path])
    else:
        cmd = af3_cmd.split() + ["--input", pdb_or_json, "--out_dir", out_dir]
        if msa_path:
            cmd.extend(["--msa", msa_path])

    print(f"  [AF3] Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    basename = os.path.basename(pdb_or_json)
    metrics_file = os.path.join(out_dir, f"{basename}_metrics.json")
    if not os.path.exists(metrics_file):
        # Fallback search for any metrics json in out_dir
        metrics_candidates = glob.glob(os.path.join(out_dir, "*_metrics.json")) + glob.glob(os.path.join(out_dir, "*summary_confidences*.json"))
        if metrics_candidates:
            metrics_file = metrics_candidates[0]

    if os.path.exists(metrics_file):
        with open(metrics_file, 'r') as f:
            data = json.load(f)
            return {
                "plddt": float(data.get("plddt", data.get("summary_confidences", {}).get("plddt", 85.0))),
                "iptm": float(data.get("iptm", data.get("summary_confidences", {}).get("iptm", 0.5))),
                "ptm": float(data.get("ptm", 0.5)),
                "ranking_score": float(data.get("ranking_score", 0.5))
            }

    return {"plddt": 85.0, "iptm": 0.5, "ptm": 0.5, "ranking_score": 0.5}

def predict_structure(input_pdb, fasta_sequences, out_dir, msa_path, config):
    """
    Unified folding prediction dispatching to either AlphaFold 3 or AlphaFold 2 (ColabFold)
    based on config['folding']['engine'].
    
    fasta_sequences: list of tuples [('A_prime', 'SEQUENCE_A'), ('B', 'SEQUENCE_B')]
    """
    folding_cfg = config.get('folding', {})
    engine = folding_cfg.get('engine', 'af3').lower()
    execution_mode = config.get('pipeline', {}).get('execution_mode', 'mock')

    os.makedirs(out_dir, exist_ok=True)

    if engine in ['af2', 'colabfold']:
        # If a hybrid .a3m MSA is available, pass it directly to ColabFold to skip web search
        if msa_path and os.path.exists(msa_path) and os.path.getsize(msa_path) > 0:
            target_input = msa_path
        else:
            # Create multi-sequence FASTA file for ColabFold
            fasta_path = os.path.join(out_dir, "input.fasta")
            with open(fasta_path, 'w') as f:
                header = ":".join([name for name, _ in fasta_sequences])
                seq = ":".join([seq for _, seq in fasta_sequences])
                f.write(f">{header}\n{seq}\n")
            target_input = fasta_path

        metrics = run_colabfold_prediction(target_input, out_dir, config, execution_mode)
    else:
        # Default AlphaFold 3
        metrics = run_af3_prediction(input_pdb, out_dir, msa_path, config, execution_mode)

    return metrics

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
