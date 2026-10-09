# Official iGEM 2026 Wiki Package (Evry-Paris-Saclay Team)

This repository contains the cleaned, rigorously structured, and validated assets for the **Model** page of the iGEM 2026 wiki.

## Directory Structure

```
wiki/
├── figures/                                    # High-resolution PNG/GIF/MP4 figures generated directly from real pipeline data
│   ├── figure1_triad_2d_vs_3d.png             # 2D (RNA) vs 3D (PPI + 16S rRNA pocket) architectural coupling
│   ├── figure2_universal_pathways.png         # Universal downstream metabolic & synthetic pathway applications
│   ├── figure3_cross_reactivity_matrix.png    # 4-state thermodynamic affinity matrix & F_ortho definition
│   ├── fig1_iptm_distributions_and_orthogonality.png # Empirical selectivity landscape across 600 complexes & 19 champions
│   ├── fig2_correlations_rf3_af2_pyrosetta.png # Multi-engine correlation matrix: RF3 vs AF2 vs PyRosetta
│   ├── fig3_mutation_frequency_and_conservation.png # Positional mutation frequencies & evolutionary conservation
│   ├── fig4_plddt_profiles_per_residue.png    # Per-residue pLDDT confidence profiles across structural domains
│   ├── fig5_champion_radar_profiles.png       # Multi-parametric radar profiles for top candidate pairs
│   ├── pipeline_diaporama_animation.gif       # Looping animated GIF of the all-atom structural trajectory
│   └── pipeline_diaporama_video.mp4           # HD MP4 video (28s) with cross-dissolve transitions and parameter overlays
│
├── slides/                                     # HD 1600x900 PyMOL 3D all-atom renders per module (slide_1.png to slide_6.png)
├── diaporama_pipeline.html                    # Standalone interactive Web viewer (timeline, video toggle, technical HUD)
├── src_gitlab/                                 # Production GitLab wiki source files ready for direct deployment
│   ├── model.tsx                              # React 19 / Bootstrap 5 component (for src/contents/model.tsx)
│   └── model.bib                              # Complete BibTeX bibliography (for src/references/model.bib)
│
├── model.md                                   # Comprehensive, mathematically and biophysically formalized Markdown documentation
├── index.html                                 # Standalone local interactive browser preview (KaTeX, Mermaid)
└── README.md                                  # This organizational overview
```

## Integrated Data & Generative Modeling Pipeline

1. **5-Module Computational Architecture:**
   - **Module 1 (Topological Mapping):** Heavy-atom distance $\le 6.0\text{ \AA}$ defines interface residues; Chain B cropped to $15.0\text{ \AA}$ to preserve the 16S rRNA steric pocket under 24 GB GPU VRAM constraints. Loop consolidation (`rescue_min_fixed_segment: 3`) merges fragmented anchors into 3 continuous loops (20, 19, 8 residues). Outward contact normal shifts the diffusion origin $+10.0\text{ \AA}$ into solvent, preventing initial clashing backbones.
   - **Module 2 (A' Rupture Design):** Apo-state diffusion (`rupture_apo_design: true`) uncouples Chain A from any favorable packing against $B_{\text{WT}}$, preventing positive wild-type design leakage. LigandMPNN electrostatic inversion and steric truncation biases remodel the interface, relaxed via PyRosetta Ref2015.
   - **Module 3 (Fail-Fast Screening):** Fast rejection filters: autonomous monomer folding ($\text{pLDDT} \ge 80.0$, $\text{RMSD} \le 2.0\text{ \AA}$) and confirmed host disruption with $B_{\text{WT}}$ ($\text{iPTM} \le 0.35$). Farthest-point Hamming diversity clustering selects distinct motifs.
   - **Module 4 (B' Rescue Co-Adaptation):** SVD rigid-body realignment resets ~70 Å coordinate shifts. Loop drift filtering (`scaffold_max_flex_rmsd} \le 4.5\text{ \AA}`) ensures local adaptation while freezing the native framework ($\text{framework RMSD} \le 1.5\text{ \AA}$). Specificity-by-difference LigandMPNN designs reciprocal salt bridges and aromatic packing against fixed A'.
   - **Module 5 (Final Multi-Engine Evaluation):** Unified masked-MSA inference across all 4 cross-reactivity states (wild-type ceiling = 0.5446) removes ancestral co-evolutionary bias. PyRosetta FastRelax binding free energy, CAPRI DockQ scoring, and SQLite results database export.

2. **Full Production Campaign Results (600 Complexes):**
   - 600 all-atom complexes sampled across 12 distinct $A'$ lineage families.
   - 70.3% (422/600) preserve native binding geometry ($\text{DockQ} \ge 0.23$), including 20 CAPRI High-Quality complexes ($\text{DockQ} \ge 0.80$, peaking at 0.9052).
   - 19 Elite Pareto Champions strictly validated: complete wild-type disruption ($\text{iPTM} \le 0.35$), high cognate rescue ($\text{iPTM} \ge 0.60 - 0.78$), selectivity margin up to $\Delta\text{iPTM} = +0.4936$ (`cand_01_B_cand_40`), and binding energies down to $-74.69\text{ REU}$.

3. **Independent AlphaFold2-Multimer v3 Cross-Validation:**
   - 57 independent full-length multimeric predictions ($A' \cdot B'$, $A' \cdot B_{\text{WT}}$, $A_{\text{WT}} \cdot B'$).
   - Strong rank correlation with RoseTTAFold-3 All-Atom ($\rho = 0.754$, $p = 0.00019$).
   - Consensus on the `cand_01` lineage maintaining native docking geometry under 1.8 Å while expelling the wild-type partner to 13.26 Å.

4. **Biophysical Pareto Gating & Proposed Wet-Lab Cascade:**
   - 7-gate Pareto cascade guaranteeing selectivity ($\Delta\text{iPTM} \ge +0.19$), abolition of host cross-talk ($\text{iPTM} \le 0.35$), and native CAPRI geometry.
   - Proposed 3-tier experimental validation: 4x4 matrix BLI/SPR kinetics, 1:1 stoichiometry SEC-MALS, and in vivo dual-reporter sfGFP/mCherry cross-repression assays.
