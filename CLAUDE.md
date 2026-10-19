# OrthoPPInterface — Project Memory

Pipeline for **orthogonal** protein-protein interface redesign: starting from a native complex A·B, engineer A' and B' such that A'·B' bind strongly while neither A'·B_WT nor A_WT·B' bind.
This is a problem of **specificity**, not affinity: any suggestion that only improves the affinity of A'·B' without widening the gap with cross-pairs misses the mark.

## Commands

```bash
conda activate orthoppinterface          # always activate, never invoke the env's python directly
python run_pipeline.py --run_name <name>              # launches the complete pipeline (modules 1-5)
python run_pipeline.py --run_name <name> --steps 4-5  # direct subset of modules
python pipeline/src/04_generate_B_prime.py --config config.yaml --force   # isolated module
python pipeline/src/calibrate_predictor.py --config config.yaml --analysis_dir runs/<run>/01_analysis
pytest pipeline/tests                     # fast, GPU-free: engines are mocked with fakes (78 tests)
```

## Repository Architecture

- `pipeline/`: Core algorithmic bioengineering (`pipeline/src/`, `pipeline/tests/`, `pipeline/bin/`).
- `config.yaml`: Unified pipeline configuration at repository root.
- `run_pipeline.py`: Master pipeline runner at repository root.
- `runs/`: Unified storage directory for all experiment run results (git-ignored).
- `wiki/`: Complete iGEM 2026 wiki package (`model.md`, `index.html`, `src_gitlab/model.tsx`, figures).
- `dossier_champions/`: Profiles and structures of the 19 all-atom elite champions.

| Module | File | Role |
| --- | --- | --- |
| 1 | `pipeline/src/01_analyze_interface.py` | Interface residues, rescue hotspots, loop consolidation for A (`rupture_min_fixed_segment`) & B (`rescue_min_fixed_segment`), loop drift ceiling (`scaffold_max_flex_rmsd`), LigandMPNN AA omissions (`ban_cysteines`, `omit_AAs`) |
| 2 | `pipeline/src/02_generate_A_prime.py` | RFD3 apo diffusion + LigandMPNN (cysteine-free) → A' candidates (rupture) |
| 3 | `pipeline/src/03_filter_A_prime.py` | Fail-fast RF3: monomer stability with dual-RMSD gating (global $\le 3.0$ Å and structural core `rmsd_framework_max` $\le 1.5$ Å), A'·B_WT rupture, WT ceiling |
| 4 | `pipeline/src/04_generate_B_prime.py` | RFD3 multi-scaffold co-adaptation (SVD realignment) + specificity-biased LigandMPNN + multi-objective composite ranking ($Z_{\text{spec}} + Z_{\text{conf}} + Z_{\text{energy}}$ with PyRosetta `dG_pack`) + clash filtering with fallback `hard_max_clashes_with_a_prime` → B' |
| 5 | `pipeline/src/05_eval_final.py` | Stage 1 (rescue A'·B' + CAPRI geometric early exit `dockq_min: 0.23`, `lrms_max: 15.0 Å` with real-time log `rescue_stage1_progress.csv`) + Stage 2 (cross-rupture + PyRosetta Ref2015 energy + CAPRI DockQ `pass_dockq` + $F_{\text{ortho}}$) + CSV/SQLite export + reports `export_styled_reports.py` (`SUMMARY.html`, `SUMMARY.xlsx`, `figures/` directory with 6 vector SVG/PNG 300 DPI figures) |

Transversal modules: `design_utils.py` (sequences/structures, `get_mpnn_omit_args`), `folding_engine.py` (RF3, MSA, cache), `energy_engine.py` + `interface_energy.py` (PyRosetta), `scoring.py` (F_ortho, ipSAE), `dockq.py`, `export_styled_reports.py`.

Every experiment lives in `runs/<run_name>/` (git-ignored) with `config_used.yaml` archived.
Key outputs: `runs/<run>/05_final_eval/orthogonality_scores.csv`, `results.db`, `coherence.json`, `runs/<run>/03_fail_fast/wt_control.json`, `rescue_stage1_progress.csv`, `SUMMARY.html`, `SUMMARY.xlsx`, `figures/`.

## Invariants Never to Break

- **Same information regime across all four folds.** A'·B', A_WT·B', A'·B_WT and the WT/WT ceiling must be predicted under the same MSA masks. Otherwise, F_ortho measures MSA asymmetry rather than binding. Any change to `msa_mask_scope`/`msa_mask_mode`/`negative_msa_regime` applies to all four.
- **Missing values are NaN.** Never substitute a default value for a failed measurement, neither in code nor in analyses.
- **RF3 cache is cryptographically signed** (sequences + MSA content + inference parameters). Do not bypass the signature; any regime shift must invalidate the cache.
- **A' is designed alone** (`rupture_apo_design: true`). With B_WT visible, LigandMPNN packs A' against B_WT and performs positive design for the wild-type partner.
- **RFD3 re-centres its output** (~70 Å). Every scaffold must be realigned onto A' via the fixed framework before any measurement.
- **LigandMPNN does not write side chains** for redesigned positions and outputs CIF files. PyRosetta packing is mandatory (`energy.pack_designs`) before any geometric or energetic evaluation.
- **Strict cysteine omission** (`ban_cysteines: true`, `omit_AAs: "C"`) to prevent reactive free cysteines at the synthetic interface.
- **CAPRI geometric early exit**: any rescue complex that binds off-target (`DockQ < dockq_min` or `L-RMSD > lrms_max`) is immediately disqualified (`early_exit_off_target`), saving 75% GPU compute on cross-folds and energy calculations.
- **No silent fallback**: if LigandMPNN fails to generate output, the module terminates with an error; it never falls back to the raw RFD3 backbone.

## Noise and Metric Interpretation

- iPTM varies by ±0.05–0.1 between runs. Never conclude on an iPTM difference under ~0.1 between two designs. `iptm_rescue_std` is logged; use it.
- F_energy has a standard deviation ≈ 0.05 across seeds for a spread of 0.22 between designs. Useful for ranking and repeated shortlisting, not for deciding marginal differences.
- Thresholds `min_binding_fraction: 0.5` and `f_energy_min: 0.0` are **provisional**, uncalibrated. Treat them as such in any analysis.
- `ipsae_*` and `f_ipsae` are **strictly diagnostic**: tracked in `coherence.json`, not yet used in pass/fail criteria. Do not promote them to thresholds without demonstrated agreement.
- RF3 and energy fail differently: a design must satisfy both. Strong disagreement between the two on a design is a biological signal, not a bug to hide.

## Code Conventions

- Unified Python 3.12 for pipeline and engines (Foundry RF3/RFD3/LigandMPNN + PyRosetta). Biopython, pandas, numpy, pyyaml.
- Modules in `src/` start with digits: they must be imported via `tests/conftest.py::load_module`, not standard `import`.
- Config is parsed in UTF-8 (BOM tolerated); missing `pipeline`/`thresholds` sections trigger an explicit error, never a silent fallback.
- Forward slashes `/` normalized for all paths (WSL/Windows compatibility).
- Tests cover plumbing (module wiring, parsing, cache, DockQ numbering) using mock engines. They do **not** validate scientific quality — never present a passing `pytest` suite as proof of methodological correctness.

## Working Rules

- Before modifying a threshold or parameter in `config.yaml`, state which run and measurement justifies it. An uncalibrated threshold change is blindness, not tuning.
- Never launch a GPU module (`--steps 2`, `4`, `5`) without explicit user request: these represent significant GPU compute time.
- For any modification to `src/`, run `pytest tests` and report the outcome.
- `runs/`, `data/outputs/`, `*.pdb`, `*.a3m`, `*.pt` are git-ignored: never propose committing them.
