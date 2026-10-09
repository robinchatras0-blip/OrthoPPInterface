"""PyRosetta helpers. Run inside the PyRosetta environment (see energy_engine.py for the pipeline-side bridge).

  --mode pack   Rebuilds the side chains of designed complexes. LigandMPNN writes no side chain for the residues it
                redesigns, so designs are backbone-only at the interface until they are packed here. Optionally keeps
                some chains fixed (`fixed_chains`), e.g. A' while B' is packed against it.
  --mode score  Interface energy: FastRelax constrained to the start coordinates, then InterfaceAnalyzer. Energies are
                only comparable between complexes scored by this script with identical settings. Results are cached by
                content hash + settings in the output JSON.

Usage: python interface_energy.py --mode {pack,score} --jobs jobs.json --out result.json
                                  [--procs 16] [--relax_repeats 1] [--seed 1]
jobs.json: [{"name": "...", "pdb": "/path/in.pdb", "out": "/path/packed.pdb" (pack), "fixed_chains": ["A"] (pack),
             "partners": "A_B" (score)}, ...]
"""
import argparse
import hashlib
import json
import os
from multiprocessing import Pool


def job_key(pdb_path, relax_repeats, seed, partners):
    """Cache key: file content + every setting that changes the result."""
    h = hashlib.sha1()
    with open(pdb_path, 'rb') as f:
        h.update(f.read())
    h.update(f"|{relax_repeats}|{seed}|{partners}".encode())
    return h.hexdigest()


def pending_jobs(jobs, cache, relax_repeats, seed):
    """Jobs whose cached result is missing, failed or computed with different inputs/settings."""
    todo = []
    for j in jobs:
        key = job_key(j["pdb"], relax_repeats, seed, j.get("partners", "A_B"))
        old = cache.get(j["name"])
        if old and old.get("key") == key and "error" not in old:
            continue
        todo.append({**j, "key": key})
    return todo


def strip_hydrogens(pdb_path):
    """Rosetta writes explicit hydrogens; the pipeline only needs heavy atoms."""
    with open(pdb_path) as f:
        lines = [l for l in f if not (l.startswith(("ATOM", "HETATM")) and l[76:78].strip() == "H")]
    with open(pdb_path, "w") as f:
        f.writelines(lines)


def _init_worker(seed):
    import pyrosetta
    pyrosetta.init(f"-mute all -ex1 -ex2aro -run:constant_seed -run:jran {seed}", silent=True)


def pack_one(job):
    """Repack all side chains (current rotamers included), then minimise the chi angles; backbone untouched."""
    try:
        import pyrosetta
        from pyrosetta.rosetta.core.pack.task import TaskFactory
        from pyrosetta.rosetta.core.pack.task.operation import (IncludeCurrent, InitializeFromCommandline, OperateOnResidueSubset,
                                                                PreventRepackingRLT, RestrictToRepacking)
        from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
        from pyrosetta.rosetta.core.kinematics import MoveMap
        from pyrosetta.rosetta.protocols.minimization_packing import MinMover, PackRotamersMover
        sfxn = pyrosetta.create_score_function("ref2015")
        pose = pyrosetta.pose_from_pdb(job["pdb"])
        tf = TaskFactory()
        tf.push_back(InitializeFromCommandline())
        tf.push_back(IncludeCurrent())
        tf.push_back(RestrictToRepacking())
        for chain in job.get("fixed_chains", []):
            tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), ChainSelector(chain)))
        packer = PackRotamersMover(sfxn)
        packer.task_factory(tf)
        packer.apply(pose)
        mm = MoveMap()
        mm.set_bb(False)
        mm.set_chi(True)
        MinMover(mm, sfxn, "lbfgs_armijo_nonmonotone", 0.01, True).apply(pose)
        out = job.get("out") or job["pdb"]
        pose.dump_pdb(out)
        strip_hydrogens(out)
        return job["name"], {"ok": True, "total": sfxn(pose)}
    except Exception as e:
        return job["name"], {"ok": False, "error": f"{type(e).__name__}: {e}"[:300]}


def score_one(args):
    job, relax_repeats = args
    try:
        import pyrosetta
        from pyrosetta.rosetta.core.pose import DockingPartners
        from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
        sfxn = pyrosetta.create_score_function("ref2015")
        pose = pyrosetta.pose_from_pdb(job["pdb"])
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(sfxn, int(relax_repeats))
        relax.constrain_relax_to_start_coords(True)
        relax.apply(pose)
        partners = DockingPartners.docking_partners_from_string(job.get("partners", "A_B"))
        iam = InterfaceAnalyzerMover(partners, False, sfxn, False, False, True, False)   # repack the separated chains
        iam.apply(pose)
        return job["name"], {"key": job["key"], "dG": iam.get_separated_interface_energy(),
                             "dSASA": iam.get_interface_delta_sasa(), "unsat_hb": iam.get_interface_delta_hbond_unsat(),
                             "n_interface_res": iam.get_num_interface_residues(), "total": sfxn(pose)}
    except Exception as e:                                      # one bad complex must not kill the batch
        return job["name"], {"key": job.get("key"), "error": f"{type(e).__name__}: {e}"[:300]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=['score', 'pack'], default='score')
    ap.add_argument('--jobs', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--relax_repeats', type=int, default=1)
    ap.add_argument('--seed', type=int, default=1)
    args = ap.parse_args()

    jobs = json.load(open(args.jobs))
    results = json.load(open(args.out)) if os.path.exists(args.out) and args.mode == 'score' else {}
    if args.mode == 'score':
        todo = pending_jobs(jobs, results, args.relax_repeats, args.seed)
        work, fn = [(j, args.relax_repeats) for j in todo], score_one
    else:
        todo, work, fn = jobs, jobs, pack_one
    print(f"[interface_energy:{args.mode}] {len(jobs)} complexes, {len(jobs) - len(todo)} cached, {len(todo)} to process "
          f"on {min(args.procs, max(1, len(todo)))} processes", flush=True)
    if todo:
        with Pool(min(args.procs, len(todo)), initializer=_init_worker, initargs=(args.seed,)) as pool:
            for name, res in pool.imap_unordered(fn, work):
                results[name] = res
                bad = res.get('error')
                print(f"[interface_energy:{args.mode}] {name}: " + (bad if bad else (f"dG={res['dG']:.1f} REU" if 'dG' in res else "ok")), flush=True)
    with open(args.out, 'w') as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
