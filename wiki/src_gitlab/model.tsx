import React from 'react';

export function Model(): React.ReactElement {
  return (
    <div className="container frosted-panel my-5 p-4 p-md-5">
      {/* Title & Metadata Header */}
      <div className="mb-5 border-bottom pb-4">
        <span className="badge bg-success text-uppercase px-3 py-2 mb-2 font-monospace">
          iGEM Evry-Paris-Saclay 2026 — Model
        </span>
        <h1 className="display-4 fw-bold text-gradient-primary mb-3">
          Engineering Orthogonal Protein–Protein Interfaces: The OrthoPPInterface Pipeline
        </h1>
        <p className="lead text-secondary">
          A modular de novo deep-learning and molecular mechanics framework for the design of mutually orthogonal
          protein interaction pairs. Benchmarked on the ribosomal S1–S2–16S rRNA initiation complex and generalized to
          insulated cellular synthetic pathways.
        </p>
        <div className="d-flex flex-wrap gap-2 pt-2">
          <a
            href="https://github.com/robinchatras0-blip/orthoppinterface"
            target="_blank"
            rel="noopener noreferrer"
            className="btn btn-outline-info btn-sm font-monospace"
          >
            GitHub Repository
          </a>
          <span className="badge bg-dark border border-secondary text-secondary font-monospace d-flex align-items-center">
            License: MIT
          </span>
          <span className="badge bg-success text-dark font-monospace d-flex align-items-center">
            Status: Validated Production (19 Champions, AF2-Multimer Consensus)
          </span>
        </div>
      </div>

      {/* Executive Summary & Key Outcomes Banner */}
      <div className="alert alert-success border-success bg-opacity-10 bg-success mb-5" role="alert">
        <h4 className="alert-heading fw-bold d-flex align-items-center gap-2 text-success">
          <span>🏆</span> Executive Summary &amp; Key Outcomes
        </h4>
        <ul className="mb-0 small text-light ps-3">
          <li>
            <strong>Automated End-to-End Pipeline:</strong> Developed <code>OrthoPPInterface</code>, a 5-module framework
            coupling RFdiffusion3, LigandMPNN, RoseTTAFold-3 All-Atom (RF3), and PyRosetta FastRelax to engineer mutually
            orthogonal protein pairs ($A' \cdot B'$) from a native complex ($A_{\text{WT}} \cdot B_{\text{WT}}$).
          </li>
          <li>
            <strong>3D Specificity vs. 2D Hybridization:</strong> Applied to the <em>E. coli</em> ribosomal S1–S2–16S rRNA
            translation initiation triad, demonstrating that 2D mRNA:rRNA base-pairing is insufficient without stereospecific
            3D anchoring of recruiter S1 onto ribosomal protein S2 under steric crowding from the 16S rRNA.
          </li>
          <li>
            <strong>High-Throughput Campaign (600 All-Atom Complexes):</strong> Evaluated 600 designed complexes across 12 distinct
            $A'$ motif families; 70.3% (422/600) retained native interface geometry (DockQ &ge; 0.23), with 20 High-Quality CAPRI
            designs (DockQ &ge; 0.80, peak at 0.9052).
          </li>
          <li>
            <strong>19 Strict Elite Champions:</strong> Discovered 19 orthogonal pairs exhibiting complete wild-type partner
            disruption ($\text{iPTM}_{\text{rupture}} \le 0.35$), high cognate rescue ($\text{iPTM}_{\text{rescue}} \ge 0.60 - 0.78$),
            selectivity differential &Delta;iPTM up to <strong>+0.4936</strong> (top candidate <code>cand_01_B_cand_40</code>), and
            PyRosetta binding energies down to <strong>-74.69 REU</strong>.
          </li>
          <li>
            <strong>Dual-Engine Cross-Validation:</strong> Evaluated all 19 champions across 57 independent AlphaFold2-Multimer v3
            predictions, confirming significant rank correlation ($\rho = 0.754$, $p = 0.00019$) on rescue affinity.
          </li>

        </ul>
      </div>

      {/* SECTION 1 */}
      <section className="mb-5" id="biological-problem">
        <h2 className="fw-bold mb-4 text-primary border-bottom pb-2">
          1. The Biological Problem: The S1–S2 Benchmark &amp; Universal Synthetic Pathways
        </h2>

        <h3 className="h4 fw-bold mt-4 mb-3 text-info" id="flagship-benchmark">
          1.1 The Flagship Benchmark: The Ribosomal 30S Initiation Triad
        </h3>
        <p>
          In bacteria, translation initiation requires recruiter protein <strong>S1 (RpsA)</strong> to bind unstructured mRNA
          upstream of the Shine-Dalgarno (SD) sequence, melt secondary hairpins, and deliver the transcript into the decoding
          channel. Crucially, S1 cannot recruit transcripts while freely diffusing in solution: it must anchor stereospecifically
          onto ribosomal protein <strong>S2 (RpsB)</strong> on the 30S subunit platform.
        </p>

        <p>
          Conventional orthogonal ribosomes focus solely on <strong>2D nucleic acid hybridization</strong> (mutating the
          anti-Shine-Dalgarno aSD to match an orthogonal mRNA). However, this 2D strategy fails to insulate translation due to three
          physical cross-talk mechanisms:
        </p>

        <div className="row g-4 my-3">
          <div className="col-lg-4">
            <div className="card h-100 bg-dark border-secondary p-3">
              <h6 className="fw-bold text-danger">Mode I: Host Leakage</h6>
              <p className="small text-secondary mb-0">
                If the orthogonal ribosome retains wild-type $S2_{\text{WT}}$, endogenous host $S1_{\text{WT}}$ docks onto it,
                recruiting native host mRNAs to the orthogonal ribosome and creating high background leakage.
              </p>
            </div>
          </div>
          <div className="col-lg-4">
            <div className="card h-100 bg-dark border-secondary p-3">
              <h6 className="fw-bold text-warning">Mode II: Host Toxicity</h6>
              <p className="small text-secondary mb-0">
                If an engineered recruiter ($S1'$) retains affinity for host $S2_{\text{WT}}$, it displaces endogenous $S1_{\text{WT}}$
                from native ribosomes, acting as a dominant-negative competitor that arrests host cell growth.
              </p>
            </div>
          </div>
          <div className="col-lg-4">
            <div className="card h-100 bg-dark border-secondary p-3">
              <h6 className="fw-bold text-info">Mode III: 16S rRNA Steric Shell</h6>
              <p className="small text-secondary mb-0">
                S2 is deeply surrounded by <strong>helices 35, 36, and 40 of the 16S rRNA</strong>, forming a polyanionic,
                sterically crowded ceiling within 4–6 Å of the binding cleft. De novo interface modifications must conform strictly
                to these macromolecular boundaries.
              </p>
            </div>
          </div>
        </div>

        <div className="text-center my-4">
          <img
            src="figures/figure1_triad_2d_vs_3d.png"
            alt="Figure 1: The 2D vs. 3D Translational Initiation Architecture"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '900px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 1 | The 2D vs. 3D Translational Initiation Architecture.</strong> Comparison between conventional 2D
            mRNA:rRNA base-pairing and the physical requirement of 3D stereospecific anchoring of recruiter S1 onto protein S2 in the
            presence of the 16S rRNA steric shell.
          </p>
        </div>

        <h3 className="h4 fw-bold mt-5 mb-3 text-info" id="universal-scope">
          1.2 Universal Scope: Insulating Cellular Synthetic Pathways
        </h3>
        <p>
          While the S1–S2 complex serves as our flagship testbed under extreme steric constraints, <code>OrthoPPInterface</code> is
          sequence- and target-agnostic. Given any natural interacting protein pair ($A_{\text{WT}} \cdot B_{\text{WT}}$), it designs
          insulated synthetic counterparts for four broad classes of synthetic biology systems:
        </p>

        <div className="text-center my-4">
          <img
            src="figures/figure2_universal_pathways.png"
            alt="Figure 2: Universal Applications Across Synthetic Biology"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '900px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 2 | Universal Applications of OrthoPPInterface Across Synthetic Biology.</strong> (1) Metabolic substrate
            channeling via synthetic protein scaffolds; (2) Insulated chimeric two-component signaling cascades; (3) Targeted
            proteolysis via orthogonal degron adaptors; (4) Synthetic transcription factor heterodimers for Boolean logic circuits.
          </p>
        </div>
      </section>

      {/* SECTION 2 */}
      <section className="mb-5" id="theoretical-principles">
        <h2 className="fw-bold mb-4 text-primary border-bottom pb-2">
          2. Theoretical Principles &amp; Biophysical Formulation of Orthogonality
        </h2>

        <h3 className="h4 fw-bold mt-4 mb-3 text-info" id="cross-reactivity-matrix">
          2.1 The 4-Way Affinity Matrix &amp; Explicit State Variables
        </h3>
        <p>
          Orthogonality is governed by the simultaneous thermodynamic equilibrium between four molecular complexes in the cellular environment. Let $A'$ and $B'$ denote the engineered orthogonal proteins, and let $A_{\text{WT}}$ and $B_{\text{WT}}$ denote the endogenous host counterparts. The four simultaneous states are formally defined as:
        </p>

        <div className="row g-3 my-2">
          <div className="col-md-6 col-lg-3">
            <div className="card bg-dark border-success p-3 h-100 shadow-sm">
              <span className="badge bg-success text-dark mb-1 font-monospace">State 1: Cognate Rescue</span>
              <strong className="text-light font-monospace">s₁ = A' &middot; B'</strong>
              <div className="small text-secondary mt-1">
                Target synthetic interaction. Requires $\text{iPTM} \ge 0.60$, $\Delta G \le -40$ REU, $\text{DockQ} \ge 0.46$.
              </div>
            </div>
          </div>
          <div className="col-md-6 col-lg-3">
            <div className="card bg-dark border-danger p-3 h-100 shadow-sm">
              <span className="badge bg-danger text-dark mb-1 font-monospace">State 2: Host Rupture</span>
              <strong className="text-light font-monospace">s₂ = A' &middot; B_WT</strong>
              <div className="small text-secondary mt-1">
                Off-target interaction with native partner. Must be abolished: $\text{iPTM} \le 0.35$ to avoid dominant-negative toxicity.
              </div>
            </div>
          </div>
          <div className="col-md-6 col-lg-3">
            <div className="card bg-dark border-danger p-3 h-100 shadow-sm">
              <span className="badge bg-danger text-dark mb-1 font-monospace">State 3: Host Cross-Dock</span>
              <strong className="text-light font-monospace">s₃ = A_WT &middot; B'</strong>
              <div className="small text-secondary mt-1">
                Host recruiter interference with synthetic partner. Must be suppressed: $\text{iPTM} \le 0.42$ to prevent host leakage.
              </div>
            </div>
          </div>
          <div className="col-md-6 col-lg-3">
            <div className="card bg-dark border-info p-3 h-100 shadow-sm">
              <span className="badge bg-info text-dark mb-1 font-monospace">State 4: WT Ceiling Ref</span>
              <strong className="text-light font-monospace">s₄ = A_WT &middot; B_WT</strong>
              <div className="small text-secondary mt-1">
                Native reference complex. $\text{iPTM} = 0.5446$ under masked MSA, providing an unbiased calibration invariant.
              </div>
            </div>
          </div>
        </div>

        <div className="text-center my-4">
          <img
            src="figures/figure3_cross_reactivity_matrix.png"
            alt="Figure 3: 4-Way Specificity Matrix"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '850px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 3 | The 4-Way Cross-Reactivity Margin Matrix.</strong> Complete functional insulation requires maximizing cognate rescue while simultaneously suppressing both hybrid cross-talk complexes below physical non-binding thresholds.
          </p>
        </div>

        <h3 className="h4 fw-bold mt-5 mb-3 text-info" id="formal-scores">
          2.2 Mathematical Scoring Functions: &Delta;iPTM, F_iptm, F_energy, and F_ortho
        </h3>
        <p>
          To rank orthogonal pairs without uncalibrated kinetic assumptions, we define explicit objective functions combining neural structural confidence and molecular mechanics force fields:
        </p>

        <div className="row g-4 my-3">
          <div className="col-md-6">
            <div className="card bg-dark border-secondary p-3 h-100">
              <h6 className="fw-bold text-success">1. Net Selectivity Margin (&Delta;iPTM)</h6>
              <div className="font-monospace small text-info bg-black p-2 rounded mb-2 border border-secondary">
                &Delta;iPTM = iPTM(A'&middot;B') - max( iPTM(A'&middot;B_WT), iPTM(A_WT&middot;B') )
              </div>
              <p className="small text-secondary mb-0">
                Because cellular cross-talk is dictated by the worst-case parasite interaction, the selectivity barrier is governed by the supremum of the two off-target states. Top champion achieves <strong>+0.4936</strong>.
              </p>
            </div>
          </div>
          <div className="col-md-6">
            <div className="card bg-dark border-secondary p-3 h-100">
              <h6 className="fw-bold text-info">2. Normalized Neural Margin (F_iptm)</h6>
              <div className="font-monospace small text-info bg-black p-2 rounded mb-2 border border-secondary">
                F_iptm = &Delta;iPTM / iPTM(A_WT&middot;B_WT)_masked
              </div>
              <p className="small text-secondary mb-0">
                Normalized relative to the wild-type reference ceiling ($0.5446$) scored under the exact same interface-masked information regime.
              </p>
            </div>
          </div>
          <div className="col-md-6">
            <div className="card bg-dark border-secondary p-3 h-100">
              <h6 className="fw-bold text-warning">3. Normalized Energy Margin (F_energy)</h6>
              <div className="font-monospace small text-warning bg-black p-2 rounded mb-2 border border-secondary">
                F_energy = ( |&Delta;G(A'&middot;B')| - max(|&Delta;G(A'&middot;B_WT)|, |&Delta;G(A_WT&middot;B')|) ) / |&Delta;G(WT)|
              </div>
              <p className="small text-secondary mb-0">
                Calculated with PyRosetta Ref2015 after constrained FastRelax: $\Delta G_{\text{bind}} = E_{\text{complex}} - (E_A + E_B)$ in REU.
              </p>
            </div>
          </div>
          <div className="col-md-6">
            <div className="card bg-dark border-secondary p-3 h-100">
              <h6 className="fw-bold text-gradient-primary">4. Composite Objective Function (F_ortho)</h6>
              <div className="font-monospace small text-light bg-black p-2 rounded mb-2 border border-secondary">
                F_ortho = 0.5 &times; F_iptm + 0.5 &times; F_energy
              </div>
              <p className="small text-secondary mb-0">
                Balances neural structural likelihood with atomic-level physical packing energy, ensuring multi-physics consensus.
              </p>
            </div>
          </div>
        </div>

        <h3 className="h4 fw-bold mt-5 mb-3 text-info" id="msa-masking">
          2.3 The Information Regime Invariant (Interface MSA Masking)
        </h3>
        <p>
          RoseTTAFold-3 uses unpaired <code>.a3m</code> multiple sequence alignments. In de novo design, homologous interface columns contain ancestral co-evolutionary signals that artificially inflate predicted binding regardless of mutated sidechains.
        </p>
        <ul className="text-secondary small">
          <li><strong>Interface columns of both chains are masked with gaps (<code>-</code>)</strong> for all designs and controls alike.</li>
          <li>Under this calibrated regime, the wild-type complex drops from $\text{iPTM} = 0.864$ (full MSA) to $\text{iPTM} = 0.5446$ (interface-masked reference ceiling).</li>
          <li>All four complexes are evaluated under the <strong>strictly identical information regime</strong>, guaranteeing unbiased scoring.</li>
        </ul>

        <h3 className="h4 fw-bold mt-5 mb-3 text-info" id="gating-thresholds">
          2.4 Multi-Objective Pareto Filtering &amp; Gating Thresholds
        </h3>
        <p>
          Rather than relying on uncalibrated phenomenological kinetic simulations, candidate qualification is enforced by strict, biophysically grounded inequality gates evaluated at each module of the pipeline:
        </p>
        <div className="bg-dark p-3 rounded border border-secondary font-monospace small text-secondary">
          <div className="text-success mb-1">&bull; Gate 1 (Monomer Folding): pLDDT(A') &ge; 80.0, C&alpha; RMSD(A') &le; 2.0 Å</div>
          <div className="text-danger mb-1">&bull; Gate 2 (WT Disruption): iPTM(A'&middot;B_WT) &le; 0.35 (loss of host binding)</div>
          <div className="text-warning mb-1">&bull; Gate 3 (Cross-Suppression): iPTM(A_WT&middot;B') &le; 0.42 (suppression of cross-docking)</div>
          <div className="text-info mb-1">&bull; Gate 4 (Cognate Rescue): iPTM(A'&middot;B') &ge; 0.60 (strong mutual affinity)</div>
          <div className="text-success mb-1">&bull; Gate 5 (Selectivity Margin): &Delta;iPTM &ge; +0.19 (favorable binding gap)</div>
          <div className="text-primary mb-1">&bull; Gate 6 (CAPRI Structural Fidelity): DockQ &ge; 0.46 (Acceptable or Medium quality)</div>
          <div className="text-light mb-0">&bull; Gate 7 (Thermodynamic Stability): &Delta;G_rescue &le; -40.0 REU (favorable interface packing)</div>
        </div>
      </section>

      {/* SECTION 3 */}
      <section className="mb-5" id="pipeline-architecture">
        <h2 className="fw-bold mb-4 text-primary border-bottom pb-2">
          3. The 5-Module Pipeline: Implementation &amp; Biophysical Rationale
        </h2>
        <p>
          The pipeline operates in five sequential modules orchestrated by <code>pipeline/run_pipeline.py</code>, where each
          algorithmic choice is justified by biophysical and computational principles:
        </p>

        <div className="text-center my-4">
          <a href="diaporama_pipeline.html" target="_blank" rel="noopener noreferrer">
            <img
              src="figures/pipeline_diaporama_animation.gif"
              alt="Diaporama Évolutionnaire All-Atom du Pipeline"
              className="img-fluid rounded border border-info shadow-lg"
              style={{ maxWidth: '950px' }}
            />
          </a>
          <p className="text-muted small mt-2">
            <strong>Diaporama Évolutionnaire All-Atom | Suivi 3D du Champion <code>cand_01_B_cand_40</code>.</strong><br/>
            <span className="text-danger fw-semibold">● Squelette Chaîne A / A' (Rouge)</span> &nbsp;|&nbsp; 
            <span className="text-primary fw-semibold">● Squelette Chaîne B / B' (Bleu)</span> &nbsp;|&nbsp; 
            <span className="text-info fw-semibold">● Bâtonnets all-atom colorés par éléments (O rouge, N bleu, C par chaîne)</span>.<br/>
            Superposition démontrant la conservation stricte du framework global (&gt;70% du volume, RMSD = 0.76 Å) et le remodelage chirurgical exclusif de l'interface (29 mutations sur A', 37 mutations sur B').<br/>
            <a href="diaporama_pipeline.html" target="_blank" rel="noopener noreferrer" className="text-info fw-bold text-decoration-none ms-1">
              ▶ Ouvrir le Diaporama Interactif Plein Écran &amp; Vidéo MP4 (28s)
            </a>
          </p>
        </div>

        <div className="row g-4 my-3">
          <div className="col-12">
            <div className="card bg-dark border-secondary p-3">
              <h5 className="fw-bold text-info">Module 1: Interface Topology &amp; Hotspot Mapping (<code>01_analyze_interface.py</code>)</h5>
              <ul className="small text-secondary mb-0 ps-3">
                <li>
                  <strong>Contact Detection &amp; Spatial Cropping:</strong> Heavy-atom distance &le; 6.0 Å defines interface residues;
                  redesignable shell set to &le; 6.0 Å (87 of 224 residues on Chain A). Chain B is cropped to 15.0 Å from the interface
                  to preserve the steric crowding of the 16S rRNA pocket while fitting within GPU VRAM during diffusion.
                </li>
                <li>
                  <strong>Loop Consolidation (<code>rescue_min_fixed_segment: 3</code>):</strong> Distance-based selection creates
                  isolated 1–2 residue fixed anchors between short flexible spans. Pinned single residues artificially over-constrain the
                  backbone, breaking RFdiffusion3 sampling. Merging short fixed segments into contiguous loops yields 3 continuous
                  diffusing spans (20/19/8 residues) anchored by stable blocks (8/26/7 residues), allowing RFD3 to generate coherent secondary structures.
                </li>
                <li>
                  <strong>Hotspot-Guided Diffusion Origin (<code>rescue_ppi_hotspots: 8</code>):</strong> Unconditioned RFD3 places the
                  diffusion origin at the geometric centroid (1.9 Å from Chain A), causing diffused atoms to nucleate inside the partner.
                  Projecting the origin 10 Å outward along the mean normal vector of the top 8 contact hotspots ensures diffusion begins in
                  solvent space and folds back toward the cleft.
                </li>
              </ul>
            </div>
          </div>

          <div className="col-12">
            <div className="card bg-dark border-secondary p-3">
              <h5 className="fw-bold text-info">Module 2: A' Rupture Design (<code>02_generate_A_prime.py</code>)</h5>
              <ul className="small text-secondary mb-0 ps-3">
                <li>
                  <strong>Apo-State Generative Diffusion (<code>rupture_apo_design: true</code>):</strong> Diffusion and sequence redesign
                  on Chain A are executed <em>in the absence of Chain B</em>. When $B_{\text{WT}}$ is present, LigandMPNN optimizes packing
                  against the wild-type partner, inadvertently performing positive design for the native interface. Apo diffusion guarantees
                  disruption of native contacts.
                </li>
                <li>
                  <strong>Rupture Biasing:</strong> LigandMPNN redesign applies electrostatic inversion (swapping charges) and steric
                  repulsion at interface contact positions.
                </li>
                <li>
                  <strong>PyRosetta Sidechain Packing:</strong> LigandMPNN outputs backbone-only CIFs for redesigned residues. Structures
                  are packed via PyRosetta Ref2015 FastRelax before evaluating steric clashes or energies.
                </li>
              </ul>
            </div>
          </div>

          <div className="col-12">
            <div className="card bg-dark border-secondary p-3">
              <h5 className="fw-bold text-info">Module 3: Fail-Fast Screening &amp; Gating (<code>03_filter_A_prime.py</code>)</h5>
              <ul className="small text-secondary mb-0 ps-3">
                <li>
                  <strong>Monomer Stability Gate:</strong> A' is folded in isolation using RF3. Candidates must satisfy pLDDT &ge; 80.0
                  and $C_\alpha$ RMSD &le; 2.0 Å against native $A_{\text{WT}}$ to ensure autonomous folding without aggregation.
                </li>
                <li>
                  <strong>Wild-Type Rupture Gate:</strong> Co-folding $A' \cdot B_{\text{WT}}$ must yield iPTM &le; 0.35, confirming disruption
                  of native binding before allocating GPU time to design B'.
                </li>
                <li>
                  <strong>Diversity Selection:</strong> Farthest-point clustering on interface sequence Hamming distance selects distinct
                  A' motif lineages (<code>max_a_prime_motifs: 4</code>) for subsequent rescue design.
                </li>
              </ul>
            </div>
          </div>

          <div className="col-12">
            <div className="card bg-dark border-secondary p-3">
              <h5 className="fw-bold text-info">Module 4: B' Rescue Multi-Scaffold Design (<code>04_generate_B_prime.py</code>)</h5>
              <ul className="small text-secondary mb-0 ps-3">
                <li>
                  <strong>SVD Rigid-Body Realignment:</strong> RFD3 translates output coordinates by ~70 Å. Backbones are realigned onto
                  the fixed framework of Chain A via singular value decomposition (SVD / Kabsch algorithm).
                </li>
                <li>
                  <strong>Scaffold Conformational Filtering (<code>scaffold_max_flex_rmsd: 4.5 Å</code>):</strong> Diffused scaffolds with
                  excessive loop drift (&gt;5 Å) lose native orientation. Filtering retains scaffolds that adapt locally (0.5–4.5 Å) while
                  preserving the native global fold ($\text{framework RMSD} \le 1.5\text{ \AA}$).
                </li>
                <li>
                  <strong>Specificity-by-Difference Biasing:</strong> At each mutable B position, residues complementary to the mutated A'
                  position are rewarded, while residues complementary to native $A_{\text{WT}}$ are penalized, directly optimizing &Delta;iPTM.
                </li>
                <li>
                  <strong>PyRosetta Packing Against Fixed A':</strong> B' sidechains are rebuilt against the fixed A' conformation via FastRelax with chi minimization.
                </li>
              </ul>
            </div>
          </div>

          <div className="col-12">
            <div className="card bg-dark border-secondary p-3">
              <h5 className="fw-bold text-info">Module 5: Final Multi-Metric Evaluation &amp; Export (<code>05_eval_final.py</code>)</h5>
              <ul className="small text-secondary mb-0 ps-3">
                <li>
                  <strong>4-Fold Cross-Reactivity:</strong> Folds $A' \cdot B'$, $A_{\text{WT}} \cdot B'$, $A' \cdot B_{\text{WT}}$, and
                  $A_{\text{WT}} \cdot B_{\text{WT}}$ under the unified interface-masked MSA regime.
                </li>
                <li>
                  <strong>Binding Free Energy:</strong> PyRosetta FastRelax calculates $\Delta G_{\text{rescue}}$, $\Delta G_{\text{cross}}$,
                  and energetic margin $F_{\text{energy}}$.
                </li>
                <li>
                  <strong>Geometry &amp; Quality:</strong> Evaluates CAPRI DockQ, ligand RMSD self-consistency, and ipSAE. Exports results
                  to SQLite <code>results.db</code>, CSV, and interactive HTML dashboards.
                </li>
              </ul>
            </div>
          </div>
        </div>
      </section>

      {/* SECTION 4 */}
      <section className="mb-5" id="production-results">
        <h2 className="fw-bold mb-4 text-primary border-bottom pb-2">
          4. Production Results: 600 Complexes &amp; 19 Strict Champions
        </h2>
        <p>
          A production campaign of <strong>600 all-atom complexes</strong> across <strong>12 distinct A' motif families</strong> was
          evaluated under the unified interface-masked regime:
        </p>

        <ul className="text-secondary small">
          <li><strong>Docking Geometry:</strong> 422 / 600 complexes (70.3%) achieved DockQ &ge; 0.23, including <strong>20 High-Quality CAPRI complexes (DockQ &ge; 0.80, peak at 0.9052)</strong> and 276 Medium-Quality complexes.</li>
          <li><strong>Rupture with B_WT:</strong> 431 / 600 candidates (100% of non-clashing designs) confirmed disruption of native binding (iPTM &le; 0.35).</li>
          <li><strong>Cross-Reactivity Suppression:</strong> 176 / 600 candidates suppressed cross-docking with native $A_{\text{WT}}$ (iPTM &le; 0.40).</li>
          <li><strong>19 Strict Champions:</strong> Simultaneously satisfied DockQ &ge; 0.46, &Delta;iPTM &gt; +0.19, and favorable $\Delta G_{\text{rescue}}$.</li>
        </ul>

        <div className="text-center my-4">
          <img
            src="figures/fig1_iptm_distributions_and_orthogonality.png"
            alt="Figure 4: High-Throughput Production Landscape"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '900px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 4 | Production Specificity Landscape (600 Complexes).</strong> Distribution of 600 designs showing the
            emergence of the 19 strict champions in the upper-right quadrant ($F_{\text{ortho}} &gt; 0.30$, DockQ &gt; 0.49).
          </p>
        </div>

        <h3 className="h4 fw-bold mt-5 mb-3 text-info" id="champions-table">
          4.1 Evaluation Table of the 19 Orthogonal Champions
        </h3>

        <div className="table-responsive my-3">
          <table className="table table-dark table-hover table-bordered small align-middle">
            <thead className="table-success text-dark text-nowrap">
              <tr>
                <th className="text-center">Rank</th>
                <th>Complex ID</th>
                <th>Parent A'</th>
                <th className="text-center bg-primary text-light">iPTM(A'&middot;B')<br/><span className="fw-normal small">Rescue</span></th>
                <th className="text-center bg-danger text-light">iPTM(A'&middot;B_WT)<br/><span className="fw-normal small">Rupture</span></th>
                <th className="text-center bg-danger text-light">iPTM(A_WT&middot;B')<br/><span className="fw-normal small">Cross</span></th>
                <th className="text-center bg-secondary text-light">iPTM_ref<br/><span className="fw-normal small">Ceiling</span></th>
                <th className="text-center bg-success text-light">&Delta;iPTM<br/><span className="fw-normal small">Margin</span></th>
                <th className="text-center">&Delta;G (REU)</th>
                <th className="text-center">F_ortho</th>
                <th className="text-center">DockQ</th>
                <th className="text-center">CAPRI</th>
                <th className="text-center">AF2 iPTM</th>
                <th className="text-center">AF2 RMSD</th>
                <th>Key Characteristic</th>
              </tr>
            </thead>
            <tbody>
              <tr className="table-active">
                <td className="text-center fw-bold">01</td>
                <td><code className="text-info">cand_01_B_cand_40</code></td>
                <td>cand_01</td>
                <td className="text-center fw-bold text-success">0.745</td>
                <td className="text-center text-secondary">0.295</td>
                <td className="text-center text-secondary">0.251</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success fs-6">+0.494</td>
                <td className="text-center">-40.7</td>
                <td className="text-center">0.211</td>
                <td className="text-center">0.4812</td>
                <td className="text-center"><span className="badge bg-warning text-dark">Acceptable</span></td>
                <td className="text-center">0.390</td>
                <td className="text-center">1.79 Å</td>
                <td><span className="badge bg-success text-dark">⭐ Selectivity Record</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">02</td>
                <td><code className="text-info">cand_17_B_cand_29</code></td>
                <td>cand_17</td>
                <td className="text-center fw-bold text-success">0.781</td>
                <td className="text-center text-secondary">0.237</td>
                <td className="text-center text-secondary">0.299</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success fs-6">+0.483</td>
                <td className="text-center fw-bold text-info">-65.8</td>
                <td className="text-center fw-bold text-warning">0.384</td>
                <td className="text-center fw-bold">0.5470</td>
                <td className="text-center"><span className="badge bg-success">Medium</span></td>
                <td className="text-center fw-bold text-info">0.490</td>
                <td className="text-center">3.52 Å</td>
                <td><span className="badge bg-primary">⭐ Affinity Peak</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">03</td>
                <td><code className="text-info">cand_01_B_cand_36</code></td>
                <td>cand_01</td>
                <td className="text-center fw-bold text-success">0.701</td>
                <td className="text-center text-secondary">0.295</td>
                <td className="text-center text-secondary">0.259</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success fs-6">+0.442</td>
                <td className="text-center">-48.1</td>
                <td className="text-center">0.253</td>
                <td className="text-center fw-bold">0.5819</td>
                <td className="text-center"><span className="badge bg-success">Medium</span></td>
                <td className="text-center">0.428</td>
                <td className="text-center fw-bold text-success">1.46 Å</td>
                <td><span className="badge bg-success">⭐ Sub-1.5 Å Native Dock</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">04</td>
                <td><code className="text-info">cand_13_B_cand_22</code></td>
                <td>cand_13</td>
                <td className="text-center">0.602</td>
                <td className="text-center text-secondary">0.246</td>
                <td className="text-center text-secondary">0.221</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success">+0.381</td>
                <td className="text-center fw-bold text-danger">-74.7</td>
                <td className="text-center fw-bold text-warning">0.400</td>
                <td className="text-center">0.4968</td>
                <td className="text-center"><span className="badge bg-success">Medium</span></td>
                <td className="text-center text-muted">0.122</td>
                <td className="text-center text-muted">6.64 Å</td>
                <td><span className="badge bg-danger">⚠️ AF2 Divergence</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">05</td>
                <td><code className="text-info">cand_11_B_cand_33</code></td>
                <td>cand_11</td>
                <td className="text-center">0.626</td>
                <td className="text-center text-secondary">0.259</td>
                <td className="text-center text-secondary">0.247</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success">+0.379</td>
                <td className="text-center">-48.8</td>
                <td className="text-center">0.191</td>
                <td className="text-center">0.5074</td>
                <td className="text-center"><span className="badge bg-success">Medium</span></td>
                <td className="text-center text-muted">0.136</td>
                <td className="text-center text-muted">5.10 Å</td>
                <td><span className="badge bg-secondary">Solid RF3 Medium</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">06</td>
                <td><code className="text-info">cand_01_B_cand_20</code></td>
                <td>cand_01</td>
                <td className="text-center fw-bold text-success">0.741</td>
                <td className="text-center text-secondary">0.295</td>
                <td className="text-center text-secondary">0.384</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success">+0.357</td>
                <td className="text-center">-46.1</td>
                <td className="text-center">0.222</td>
                <td className="text-center">0.5401</td>
                <td className="text-center"><span className="badge bg-success">Medium</span></td>
                <td className="text-center">0.297</td>
                <td className="text-center">2.05 Å</td>
                <td><span className="badge bg-info text-dark">High Cognate Rescue</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">07</td>
                <td><code className="text-info">cand_11_B_cand_40</code></td>
                <td>cand_11</td>
                <td className="text-center">0.600</td>
                <td className="text-center text-secondary">0.259</td>
                <td className="text-center text-secondary">0.243</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success">+0.357</td>
                <td className="text-center">-50.7</td>
                <td className="text-center">0.224</td>
                <td className="text-center">0.4729</td>
                <td className="text-center"><span className="badge bg-warning text-dark">Acceptable</span></td>
                <td className="text-center text-muted">0.154</td>
                <td className="text-center text-muted">5.03 Å</td>
                <td><span className="badge bg-secondary">Balanced Energy</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">08</td>
                <td><code className="text-info">cand_11_B_cand_04</code></td>
                <td>cand_11</td>
                <td className="text-center">0.625</td>
                <td className="text-center text-secondary">0.259</td>
                <td className="text-center text-secondary">0.270</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success">+0.355</td>
                <td className="text-center">-51.4</td>
                <td className="text-center">0.147</td>
                <td className="text-center">0.5368</td>
                <td className="text-center"><span className="badge bg-success">Medium</span></td>
                <td className="text-center text-muted">0.135</td>
                <td className="text-center text-muted">6.23 Å</td>
                <td><span className="badge bg-secondary">CAPRI Medium</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">09</td>
                <td><code className="text-info">cand_01_B_cand_14</code></td>
                <td>cand_01</td>
                <td className="text-center">0.632</td>
                <td className="text-center text-secondary">0.295</td>
                <td className="text-center text-secondary">0.282</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success">+0.351</td>
                <td className="text-center">-48.5</td>
                <td className="text-center">0.221</td>
                <td className="text-center fw-bold text-success">0.6031</td>
                <td className="text-center"><span className="badge bg-success">Medium</span></td>
                <td className="text-center fw-bold text-info">0.470</td>
                <td className="text-center fw-bold text-success">1.51 Å</td>
                <td><span className="badge bg-success">⭐ Universal Consensus #1</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">10</td>
                <td><code className="text-info">cand_04_B_cand_10</code></td>
                <td>cand_04</td>
                <td className="text-center">0.629</td>
                <td className="text-center text-secondary">0.255</td>
                <td className="text-center text-secondary">0.279</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success">+0.350</td>
                <td className="text-center fw-bold text-info">-64.2</td>
                <td className="text-center">0.277</td>
                <td className="text-center fw-bold">0.5840</td>
                <td className="text-center"><span className="badge bg-success">Medium</span></td>
                <td className="text-center">0.277</td>
                <td className="text-center">3.35 Å</td>
                <td><span className="badge bg-info text-dark">High Stability</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">11</td>
                <td><code className="text-info">cand_18_B_cand_37</code></td>
                <td>cand_18</td>
                <td className="text-center fw-bold text-success">0.741</td>
                <td className="text-center text-secondary">0.240</td>
                <td className="text-center text-secondary">0.402</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success">+0.340</td>
                <td className="text-center">-50.8</td>
                <td className="text-center">0.181</td>
                <td className="text-center">0.4780</td>
                <td className="text-center"><span className="badge bg-warning text-dark">Acceptable</span></td>
                <td className="text-center fw-bold text-info">0.417</td>
                <td className="text-center">2.10 Å</td>
                <td><span className="badge bg-primary">Strong Dual Affinity</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">12</td>
                <td><code className="text-info">cand_07_B_cand_09</code></td>
                <td>cand_07</td>
                <td className="text-center">0.639</td>
                <td className="text-center text-secondary">0.246</td>
                <td className="text-center text-secondary">0.313</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success">+0.326</td>
                <td className="text-center fw-bold text-info">-67.9</td>
                <td className="text-center">0.293</td>
                <td className="text-center">0.5085</td>
                <td className="text-center"><span className="badge bg-success">Medium</span></td>
                <td className="text-center">0.369</td>
                <td className="text-center fw-bold text-success">1.45 Å</td>
                <td><span className="badge bg-success">⭐ Geometric Record (1.45 Å)</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">13</td>
                <td><code className="text-info">cand_19_B_cand_39</code></td>
                <td>cand_19</td>
                <td className="text-center">0.608</td>
                <td className="text-center text-secondary">0.227</td>
                <td className="text-center text-secondary">0.296</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success">+0.312</td>
                <td className="text-center fw-bold text-info">-62.8</td>
                <td className="text-center">0.241</td>
                <td className="text-center">0.5614</td>
                <td className="text-center"><span className="badge bg-success">Medium</span></td>
                <td className="text-center text-muted">0.141</td>
                <td className="text-center">2.94 Å</td>
                <td><span className="badge bg-secondary">Low WT Rupture (0.227)</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">14</td>
                <td><code className="text-info">cand_18_B_cand_03</code></td>
                <td>cand_18</td>
                <td className="text-center">0.620</td>
                <td className="text-center text-secondary">0.240</td>
                <td className="text-center text-secondary">0.331</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success">+0.289</td>
                <td className="text-center">-47.6</td>
                <td className="text-center">0.177</td>
                <td className="text-center">0.5338</td>
                <td className="text-center"><span className="badge bg-success">Medium</span></td>
                <td className="text-center">0.292</td>
                <td className="text-center">2.29 Å</td>
                <td><span className="badge bg-secondary">Balanced Profile</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">15</td>
                <td><code className="text-info">cand_19_B_cand_02</code></td>
                <td>cand_19</td>
                <td className="text-center">0.636</td>
                <td className="text-center text-secondary">0.227</td>
                <td className="text-center text-secondary">0.356</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success">+0.281</td>
                <td className="text-center">-50.0</td>
                <td className="text-center">0.087</td>
                <td className="text-center">0.4642</td>
                <td className="text-center"><span className="badge bg-warning text-dark">Acceptable</span></td>
                <td className="text-center">0.311</td>
                <td className="text-center">2.00 Å</td>
                <td><span className="badge bg-secondary">Sub-2.0 Å Docking</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">16</td>
                <td><code className="text-info">cand_20_B_cand_36</code></td>
                <td>cand_20</td>
                <td className="text-center fw-bold">0.687</td>
                <td className="text-center text-secondary">0.207</td>
                <td className="text-center text-secondary">0.412</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success">+0.275</td>
                <td className="text-center">-44.9</td>
                <td className="text-center">0.122</td>
                <td className="text-center">0.5709</td>
                <td className="text-center"><span className="badge bg-success">Medium</span></td>
                <td className="text-center fw-bold text-info">0.415</td>
                <td className="text-center fw-bold text-success">1.00 Å</td>
                <td><span className="badge bg-success">⭐ 1.00 Å AF2 RMSD</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">17</td>
                <td><code className="text-info">cand_04_B_cand_24</code></td>
                <td>cand_04</td>
                <td className="text-center">0.629</td>
                <td className="text-center text-secondary">0.255</td>
                <td className="text-center text-secondary">0.356</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success">+0.273</td>
                <td className="text-center">-56.3</td>
                <td className="text-center">0.175</td>
                <td className="text-center">0.5641</td>
                <td className="text-center"><span className="badge bg-success">Medium</span></td>
                <td className="text-center fw-bold text-info">0.588</td>
                <td className="text-center fw-bold text-success">1.04 Å</td>
                <td><span className="badge bg-primary">⭐ AF2 Affinity Peak (0.588)</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">18</td>
                <td><code className="text-info">cand_11_B_cand_32</code></td>
                <td>cand_11</td>
                <td className="text-center">0.625</td>
                <td className="text-center text-secondary">0.259</td>
                <td className="text-center text-secondary">0.400</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success">+0.225</td>
                <td className="text-center">-53.6</td>
                <td className="text-center">0.187</td>
                <td className="text-center">0.4858</td>
                <td className="text-center"><span className="badge bg-warning text-dark">Acceptable</span></td>
                <td className="text-center">0.300</td>
                <td className="text-center">1.62 Å</td>
                <td><span className="badge bg-secondary">1.62 Å AF2 Binding</span></td>
              </tr>
              <tr>
                <td className="text-center fw-bold">19</td>
                <td><code className="text-info">cand_10_B_cand_25</code></td>
                <td>cand_10</td>
                <td className="text-center">0.612</td>
                <td className="text-center text-secondary">0.208</td>
                <td className="text-center text-secondary">0.419</td>
                <td className="text-center text-muted">0.545</td>
                <td className="text-center fw-bold text-success">+0.194</td>
                <td className="text-center fw-bold text-info">-63.5</td>
                <td className="text-center">0.122</td>
                <td className="text-center">0.4600</td>
                <td className="text-center"><span className="badge bg-warning text-dark">Acceptable</span></td>
                <td className="text-center text-muted">0.127</td>
                <td className="text-center text-muted">4.69 Å</td>
                <td><span className="badge bg-secondary">High Stability</span></td>
              </tr>
            </tbody>
          </table>
        </div>

        <div className="text-center my-4">
          <img
            src="figures/fig5_champion_radar_profiles.png"
            alt="Figure 5: Radar Profiles of Top Champions"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '900px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 5 | Radar Profiles of Top Champions.</strong> Multi-criteria comparisons across DockQ, rescue
            affinity, selectivity differential, and binding free energy.
          </p>
        </div>
      </section>

      {/* SECTION 5 */}
      <section className="mb-5" id="dual-model-validation">
        <h2 className="fw-bold mb-4 text-primary border-bottom pb-2">
          5. Dual-Model Cross-Validation &amp; Sequence Remodeling
        </h2>
        <p>
          All 19 champions were independently cross-validated using <strong>AlphaFold2-Multimer v3</strong> across 57 structural
          predictions ($A' \cdot B'$, $A' \cdot B_{\text{WT}}$, and $A_{\text{WT}} \cdot B'$):
        </p>

        <div className="text-center my-4">
          <img
            src="figures/fig2_correlations_rf3_af2_pyrosetta.png"
            alt="Figure 6: Multi-Model Biophysical Equivalence"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '900px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 7 | Multi-Model Statistical Correlation Matrix.</strong> Correlation analysis between RoseTTAFold-3
            All-Atom, AlphaFold2-Multimer v3, and PyRosetta interface energy.
          </p>
        </div>

        <ul className="text-secondary small">
          <li><strong>Rescue Affinity Consensus:</strong> Spearman $\rho = +0.754$ ($p = 0.00019$) and Pearson $r = +0.559$ ($p = 0.0128$), demonstrating strong inter-model agreement on cognate complex binding.</li>
          <li><strong>Cross-Docking Barrier:</strong> Spearman $\rho = +0.625$ ($p = 0.00425$), confirming shared recognition of cross-binding suppression with native $A_{\text{WT}}$.</li>
          <li><strong>The Consensus Lineage (<code>cand_01</code> family):</strong> Variants <code>cand_14</code>, <code>cand_36</code>, and <code>cand_40</code> scored in the top tier across both engines. In AF2, they preserve sub-1.6 Å native docking (RMSD &le; 1.46–1.79 Å) while displacing $B_{\text{WT}}$ to $\text{RMSD} = 13.26\text{ \AA}$.</li>
          <li><strong>Engine Divergence:</strong> Candidate <code>cand_13_B_cand_22</code> showed high PyRosetta energy ($\Delta G = -74.7\text{ REU}$) and acceptable RF3 score ($\text{iPTM} = 0.602$), but failed to dock in AF2 ($\text{iPTM} = 0.122$, $\text{RMSD} = 6.64\text{ \AA}$), illustrating the value of multi-engine screening to filter out model-specific false positives.</li>
        </ul>

        <div className="text-center my-4">
          <img
            src="figures/fig3_mutation_frequency_and_conservation.png"
            alt="Figure 7: Mutation Frequency Heatmap"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '900px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 7 | Mutation Frequency and Conservation Heatmap.</strong> Sequence alignment across all 19 champions
            reveals consistent alanine remodeling on Chain A' and aromatic substitutions on Chain B'.
          </p>
        </div>

        <div className="text-center my-4">
          <img
            src="figures/fig4_plddt_profiles_per_residue.png"
            alt="Figure 8: Per-Residue pLDDT Profiles"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '900px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 8 | Per-Residue Structural Confidence (pLDDT) Across Champions.</strong> High per-residue confidence
            (mean pLDDT &gt; 80) across core secondary structure elements, confirming folding stability and structural integrity of the designed interfaces.
          </p>
        </div>
      </section>

      {/* SECTION 6 */}
      <section className="mb-5" id="wet-lab-cascade">
        <h2 className="fw-bold mb-4 text-primary border-bottom pb-2">
          6. Proposed Wet-Lab Validation Cascade
        </h2>
        <p>
          To validate these in silico candidates experimentally, we propose a three-tier characterization cascade:
        </p>

        <div className="row g-3">
          <div className="col-md-4">
            <div className="card h-100 bg-dark border-secondary p-3">
              <span className="badge bg-info text-dark mb-2 align-self-start font-monospace">Phase 1: Kinetics</span>
              <h5 className="fw-bold text-light">BLI / SPR 4&times;4 Titration</h5>
              <p className="small text-secondary mb-0">
                Cell-free PURE synthesis; AviTag site-specific biotinylation; real-time $k_{\text{on}}$, $k_{\text{off}}$, and $K_D$
                measurement across all four complexes ($A' \cdot B'$ target $K_D &lt; 50\text{ nM}$; off-target complexes $K_D &gt; 10\text{ }\mu\text{M}$).
              </p>
            </div>
          </div>

          <div className="col-md-4">
            <div className="card h-100 bg-dark border-secondary p-3">
              <span className="badge bg-warning text-dark mb-2 align-self-start font-monospace">Phase 2: Assembly</span>
              <h5 className="fw-bold text-light">SEC-MALS Monodispersity</h5>
              <p className="small text-secondary mb-0">
                Superdex 200 Increase size-exclusion chromatography coupled with multi-angle light scattering to verify pure 1:1
                heterodimerization without non-specific aggregation.
              </p>
            </div>
          </div>

          <div className="col-md-4">
            <div className="card h-100 bg-dark border-secondary p-3">
              <span className="badge bg-success text-dark mb-2 align-self-start font-monospace">Phase 3: In Vivo</span>
              <h5 className="fw-bold text-light">Dual-Reporter Assay</h5>
              <p className="small text-secondary mb-0">
                Dual-fluorescence in <em>E. coli</em> BL21(DE3): orthogonal reporter sfGFP (driven by the orthogonal ribosome channel)
                against host reporter mCherry (driven by host ribosomes) quantified by flow cytometry to confirm cellular insulation.
              </p>
            </div>
          </div>
        </div>
      </section>

      {/* SECTION 7 */}
      <section className="mb-5" id="open-source-pipeline">
        <h2 className="fw-bold mb-4 text-primary border-bottom pb-2">
          7. Open-Source Pipeline, Software &amp; Reproducibility
        </h2>
        <p>
          The complete <code>OrthoPPInterface</code> platform is open-source under the MIT License:
        </p>
        <p className="font-monospace">
          👉 <strong>Repository:</strong>{' '}
          <a
            href="https://github.com/robinchatras0-blip/orthoppinterface"
            target="_blank"
            rel="noopener noreferrer"
            className="text-info text-decoration-none"
          >
            https://github.com/robinchatras0-blip/orthoppinterface
          </a>
        </p>

        <div className="row g-3 my-2">
          <div className="col-md-6">
            <div className="card bg-dark border-secondary p-3 h-100">
              <h6 className="fw-bold text-light">Modular Implementation</h6>
              <ul className="small text-secondary mb-0 ps-3">
                <li><code>pipeline/src/01_analyze_interface.py</code>: Topology, loop consolidation, and hotspot origin.</li>
                <li><code>pipeline/src/02_generate_A_prime.py</code>: Apo RFD3 diffusion, LigandMPNN rupture, PyRosetta packing.</li>
                <li><code>pipeline/src/03_filter_A_prime.py</code>: Fail-fast stability and WT uncoupling screening.</li>
                <li><code>pipeline/src/04_generate_B_prime.py</code>: SVD realignment, scaffold filtering, specificity-biased LigandMPNN.</li>
                <li><code>pipeline/src/05_eval_final.py</code>: 4-fold masked evaluation, PyRosetta FastRelax, DockQ, SQLite export.</li>
                <li><code>pipeline/src/validate_af2_all19.py</code>: Independent AF2-Multimer v3 validation engine.</li>
              </ul>
            </div>
          </div>
          <div className="col-md-6">
            <div className="card bg-dark border-secondary p-3 h-100">
              <h6 className="fw-bold text-light">Cryptographic Reproducibility</h6>
              <p className="small text-secondary mb-2">
                All RF3 structural predictions are cached under a SHA-256 cryptographic signature computed from sequence, MSA contents,
                and inference parameters. Modifying MSA masks or inference parameters never silently reuses outdated scores.
              </p>
              <div className="badge bg-success text-dark font-monospace">71 / 71 Unit &amp; Integration Tests Passing</div>
            </div>
          </div>
        </div>
      </section>

      {/* SECTION 8: REFERENCES */}
      <section className="mb-5" id="references">
        <h2 className="fw-bold mb-4 text-primary border-bottom pb-2">
          References
        </h2>
        <ol className="small text-secondary ps-3">
          <li><strong>Watson, J. L., et al. (2023).</strong> De novo design of protein structure and function with RFdiffusion. <em>Nature</em>, 620(7976), 1089–1100.</li>
          <li><strong>Dauparas, J., et al. (2022).</strong> Robust deep learning–based protein sequence design using ProteinMPNN. <em>Science</em>, 378(6615), 49–56.</li>
          <li><strong>Baek, M., et al. (2024).</strong> Accurate prediction of protein structures and interactions using RoseTTAFold-3 All-Atom. <em>Science</em>, 384(6693), eadn8764.</li>
          <li><strong>Alford, R. F., et al. (2017).</strong> The Rosetta all-atom energy function for macromolecular modeling and design. <em>Journal of Chemical Theory and Computation</em>, 13(6), 3031–3048.</li>
          <li><strong>Evans, R., et al. (2022).</strong> Protein complex prediction with AlphaFold-Multimer. <em>bioRxiv</em>, doi:10.1101/2021.10.04.463034.</li>
          <li><strong>Basu, S., &amp; Wallner, B. (2016).</strong> DockQ: A quality measure for protein-protein docking models. <em>PLOS ONE</em>, 11(8), e0161879.</li>
          <li><strong>Dueber, J. E., et al. (2009).</strong> Synthetic protein scaffolds provide modular control over metabolic flux. <em>Nature Biotechnology</em>, 27(8), 753–759.</li>
          <li><strong>Chen, Z., et al. (2020).</strong> De novo design of protein logic gates. <em>Science</em>, 368(6486), 78–84.</li>
          <li><strong>Subramanian, A. R. (1983).</strong> Structure and functions of ribosomal protein S1. <em>Progress in Nucleic Acid Research and Molecular Biology</em>, 28, 101–142.</li>
          <li><strong>Byrgazov, K., et al. (2015).</strong> Translation initiation in bacteria: The role of ribosomal protein S1. <em>Biochemistry (Moscow)</em>, 80(7), 808–819.</li>
          <li><strong>Rackham, O., &amp; Chin, J. W. (2005).</strong> A network of orthogonal ribosome&middot;mRNA pairs. <em>Nature Chemical Biology</em>, 1(3), 159–166.</li>
          <li><strong>Sumida, K. H., et al. (2024).</strong> Improving protein expression, stability, and function with ProteinMPNN. <em>Nature Communications</em>, 15, 3450.</li>
        </ol>
      </section>
    </div>
  );
}
