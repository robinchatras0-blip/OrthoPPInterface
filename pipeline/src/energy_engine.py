"""Pipeline-side bridge to interface_energy.py, which needs PyRosetta and therefore runs in the PyRosetta environment
(through WSL on Windows, directly on Linux), the same way RF3 is called."""
import json
import os
import subprocess
import sys

sys.path.append(os.path.dirname(__file__))
from folding_engine import to_wsl_path, _use_wsl  # noqa: E402

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "interface_energy.py")


def enabled(config):
    return bool(config.get('energy', {}).get('enabled', False))


def _run(mode, jobs, result_json, config, path_keys):
    ecfg = config.get('energy', {})
    wsl = _use_wsl(config)
    pth = to_wsl_path if wsl else (lambda p: os.path.abspath(p).replace('\\', '/'))
    os.makedirs(os.path.dirname(os.path.abspath(result_json)), exist_ok=True)
    jobs_path = os.path.splitext(result_json)[0] + "_jobs.json"
    with open(jobs_path, 'w') as f:
        json.dump([{**j, **{k: pth(j[k]) for k in path_keys if k in j}} for j in jobs], f, indent=2)
    python = ecfg.get('python', 'python' if wsl else sys.executable)
    cmd = [python, pth(SCRIPT), "--mode", mode, "--jobs", pth(jobs_path), "--out", pth(result_json),
           "--procs", str(ecfg.get('n_procs', 8)), "--relax_repeats", str(ecfg.get('relax_repeats', 1)),
           "--seed", str(ecfg.get('seed', 1))]
    if wsl:
        cmd = ["wsl", "-d", config.get('folding', {}).get('wsl_distro', 'Ubuntu'), "--"] + cmd
    print(f"  [PyRosetta {mode}] {len(jobs)} complexes", flush=True)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"PyRosetta {mode} failed ({res.returncode}): {(res.stderr or res.stdout)[-500:]}")
    with open(result_json) as f:
        results = json.load(f)
    failed = [n for n, r in results.items() if r.get("error")]
    if failed:
        print(f"  [PyRosetta {mode}] WARNING: {len(failed)} complex(es) failed: {failed[:3]}")
    return results


def score_interfaces(jobs, out_json, config):
    """Interface energy of complexes. jobs: [{"name", "pdb", "partners"?}] -> {name: {"dG", "dSASA", ...}} (cached)."""
    return _run("score", jobs, out_json, config, ("pdb",))


def pack_structures(jobs, config, status_json):
    """Rebuilds side chains in place (or to job["out"]). jobs: [{"name", "pdb", "out"?, "fixed_chains"?}].
    Raises if any complex could not be packed: a half-packed design set must not go on silently."""
    results = _run("pack", jobs, status_json, config, ("pdb", "out"))
    failed = {n: r.get("error") for n, r in results.items() if not r.get("ok")}
    if failed:
        raise RuntimeError(f"side-chain packing failed for {len(failed)} structure(s), e.g. {next(iter(failed.items()))}")
    return results
