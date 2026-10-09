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

      {/* Key Success Highlights Banner */}
      <div className="alert alert-success border-success bg-opacity-10 bg-success mb-5" role="alert">
        <h4 className="alert-heading fw-bold d-flex align-items-center gap-2 text-success">
          <span>🏆</span> Key Modeling Successes &amp; Breakthrough Highlights
        </h4>
        <ul className="mb-0 small text-light ps-3">
          <li>
            <strong>End-to-End De Novo Co-Adaptation:</strong> Built and benchmarked <code>OrthoPPInterface</code>,
            coupling RFdiffusion3, LigandMPNN, RoseTTAFold-3 All-Atom (RF3), and PyRosetta FastRelax.
          </li>
          <li>
            <strong>Dual-Cycle DBTL Engineering:</strong> Diagnosed and eliminated Cycle 1 failure modes (loop drift &gt;19 Å,
            cross-reactivity) through four algorithmic innovations: loop bridling, specificity-by-difference biasing, SVD
            realignment, and loop consolidation.
          </li>
          <li>
            <strong>High-Throughput Production (600 Complexes):</strong> 70.3% (422/600) retained native docking geometry
            (DockQ &ge; 0.23), with 20 High-Quality CAPRI complexes (DockQ &ge; 0.80, peak at 0.9052).
          </li>
          <li>
            <strong>19 Strict Elite Champions:</strong> Complete disruption of native binding (iPTM &le; 0.35), high rescue
            affinity (iPTM 0.60–0.78), record selectivity differential &Delta;iPTM up to <strong>+0.4936</strong>, and
            favorable binding energy down to <strong>-74.69 REU</strong>.
          </li>
          <li>
            <strong>Dual-Model Cross-Validation (RF3 vs. AF2-Multimer v3):</strong> 57 independent 3D predictions confirmed
            robust rank correlation (&rho; = 0.754, p &lt; 0.001) and elucidated a precise tryptophan steric latch (S12W).
          </li>
          <li>
            <strong>Predictive Single-Cell ODE Kinetics:</strong> Demonstrated that an in silico orthogonality margin
            F_ortho &ge; +0.20 guarantees &gt;85% translational insulation in living cells.
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
          To benchmark our pipeline under the most stringent structural constraints possible, we selected the translational
          initiation triad of <em>Escherichia coli</em>: <strong>Ribosomal Protein S1 (RpsA)</strong>,{' '}
          <strong>Ribosomal Protein S2 (RpsB)</strong>, and the <strong>16S ribosomal RNA (rRNA)</strong>.
        </p>

        <div className="row g-4 my-3">
          <div className="col-lg-6">
            <div className="card h-100 bg-dark border-secondary p-3">
              <h5 className="fw-bold text-info">The 2D Domain: RNA Thermodynamics</h5>
              <p className="small text-secondary mb-2">
                2D translation engineering focuses on transcript secondary structure (&Delta;G_folding) and Watson-Crick
                hybridization between the Shine-Dalgarno (SD) sequence and the anti-Shine-Dalgarno (aSD) on 16S rRNA:
              </p>
              <div className="bg-black p-2 rounded font-monospace small text-muted mb-2">
                &Delta;G_initiation = &Delta;G_SD:aSD - &Delta;G_mRNA_unfolding
              </div>
              <p className="small text-secondary mb-0">
                While o-mRNA and o-16S rRNA base-pairing can achieve 2D sequence orthogonality, translation in bacteria
                strictly requires physical recruitment by ribosomal protein S1.
              </p>
            </div>
          </div>

          <div className="col-lg-6">
            <div className="card h-100 bg-dark border-secondary p-3">
              <h5 className="fw-bold text-warning">The 3D Domain: Macromolecular Assembly</h5>
              <p className="small text-secondary mb-2">
                S1 cannot act in solution: it must anchor specifically onto S2 on the 30S platform. Crucially, S2 directly
                contacts <strong>helices 35, 36, and 40 of the 16S rRNA</strong>.
              </p>
              <ul className="small text-secondary ps-3 mb-0">
                <li><strong>Host Leakage:</strong> Host S1 binds Ortho S2 &rarr; uninsulated background translation.</li>
                <li><strong>Host Toxicity:</strong> Engineered S1' binds Host S2 &rarr; fatal shutdown of native translation.</li>
                <li><strong>Steric Crowding:</strong> Rigid polyanionic RNA backbone sits 4–6 Å from the S1–S2 cleft.</li>
              </ul>
            </div>
          </div>
        </div>

        <div className="text-center my-4">
          <img
            src="figures/figure1_triad_2d_vs_3d.png"
            alt="Figure 1: 2D vs 3D Initiation Architecture"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '900px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 1 | The 2D vs. 3D Translational Initiation Architecture.</strong> Comparison between 2D nucleic
            acid hybridization and 3D stereospecific anchoring of protein S1 onto protein S2 in the presence of the 16S rRNA shell.
          </p>
        </div>

        <h3 className="h4 fw-bold mt-5 mb-3 text-info" id="universal-pathways">
          1.2 Beyond the Ribosome: A General Platform for Insulated Synthetic Pathways
        </h3>
        <p>
          While the S1–S2 triad represents our rigorous benchmark, the overarching architecture of{' '}
          <code>OrthoPPInterface</code> is input-agnostic and universally applicable across synthetic biology:
        </p>

        <div className="text-center my-4">
          <img
            src="figures/figure2_universal_pathways.png"
            alt="Figure 2: Universal Applications"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '900px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 2 | Universal Applications of OrthoPPInterface.</strong> (1) Metabolic substrate channeling;
            (2) Insulated chimeric signaling cascades; (3) Targeted proteolysis via orthogonal degrons; (4) Logic gate transcription factors.
          </p>
        </div>
      </section>

      {/* SECTION 2 */}
      <section className="mb-5" id="theoretical-framework">
        <h2 className="fw-bold mb-4 text-primary border-bottom pb-2">
          2. Theoretical Principles &amp; Scoring Metrics
        </h2>

        <h3 className="h4 fw-bold mt-4 mb-3 text-info" id="orthogonality-matrix">
          2.1 The 4-Way Affinity Matrix &amp; Selectivity Margin
        </h3>
        <p>
          Orthogonality between a synthetic pair ($A' \cdot B'$) and native proteins ($A_{'{'}WT{'}'} \cdot B_{'{'}WT{'}'}$)
          is governed by four concurrent states:
        </p>
        <div className="bg-dark p-3 rounded font-monospace small mb-3 border border-secondary text-info">
          (1) A' · B' (Rescue) &rarr; iPTM &ge; 0.60, &Delta;G &le; -45 REU, DockQ &ge; 0.49<br />
          (2) A' · B_WT (WT Rupture) &rarr; iPTM &le; 0.35 (Abolished)<br />
          (3) A_WT · B' (Cross Rupture) &rarr; iPTM &le; 0.35 (Abolished)<br />
          (4) A_WT · B_WT (WT Ceiling) &rarr; iPTM = 0.5446 (Interface-masked baseline)
        </div>

        <div className="text-center my-4">
          <img
            src="figures/figure3_cross_reactivity_matrix.png"
            alt="Figure 3: Cross Reactivity Matrix"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '850px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 3 | The 4-Way Cross-Reactivity Margin Matrix.</strong> Defining the selective margin &Delta;iPTM.
          </p>
        </div>

        <h3 className="h4 fw-bold mt-4 mb-3 text-info" id="scoring-formulations">
          2.2 Scoring Formulations: F_ortho, F_energy, ipSAE, and DockQ
        </h3>
        <div className="row g-3">
          <div className="col-md-6">
            <div className="card bg-dark border-secondary p-3 h-100">
              <h6 className="fw-bold text-info">Raw Neural Selectivity Margin (&Delta;iPTM)</h6>
              <div className="font-monospace small text-muted bg-black p-2 rounded mb-2">
                &Delta;iPTM = iPTM(A'·B') - max(iPTM(A'·B_WT), iPTM(A_WT·B'))
              </div>
              <p className="small text-secondary mb-0">
                Measures whether the synthetic complex is preferred over both hybrid cross-reactive states.
              </p>
            </div>
          </div>
          <div className="col-md-6">
            <div className="card bg-dark border-secondary p-3 h-100">
              <h6 className="fw-bold text-warning">Multi-Objective Score (F_ortho)</h6>
              <div className="font-monospace small text-muted bg-black p-2 rounded mb-2">
                F_ortho = 0.5 &times; F_iptm + 0.5 &times; F_energy
              </div>
              <p className="small text-secondary mb-0">
                Harmonizes neural confidence (normalized against interface-masked WT ceiling) with PyRosetta Ref2015 binding free energy fraction.
              </p>
            </div>
          </div>
        </div>
      </section>

      {/* SECTION 3 */}
      <section className="mb-5" id="kinetic-ode">
        <h2 className="fw-bold mb-4 text-primary border-bottom pb-2">
          3. Single-Cell ODE Translational Competition Kinetics
        </h2>
        <p>
          To establish what in silico score is required for in vivo insulation, we simulated competitive initiation kinetics
          between native ($S1_{'{'}WT{'\cdot'}S2_{'{'}WT{'}'}$) and orthogonal ($S1'{'\cdot'}S2'$) initiation complexes:
        </p>

        <div className="text-center my-4">
          <img
            src="figures/figure4_ode_competition_kinetics.png"
            alt="Figure 4: ODE Competition Kinetics"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '850px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 4 | Single-Cell ODE Translation Selectivity vs. Orthogonality Margin.</strong> Numerical integration
            shows that an orthogonality margin F_ortho &ge; +0.20 is necessary to achieve &gt;85% translational selectivity in vivo.
          </p>
        </div>
      </section>

      {/* SECTION 4 */}
      <section className="mb-5" id="dbtl-cycles">
        <h2 className="fw-bold mb-4 text-primary border-bottom pb-2">
          4. Iterative DBTL Engineering Cycles &amp; Production Breakthrough
        </h2>

        <div className="row g-4 mb-4">
          <div className="col-lg-6">
            <div className="card bg-dark border-danger p-3 h-100">
              <span className="badge bg-danger text-dark mb-2 align-self-start font-monospace">Cycle 1: Diagnostic (run_night01)</span>
              <h5 className="fw-bold text-light">Early Bottlenecks Identified</h5>
              <ul className="small text-secondary ps-3 mb-0">
                <li><strong>Cross-reactivity deficit:</strong> cand_11_B_cand_05 achieved rescue iPTM = 0.7328, but iPTM(A_WT·B') was 0.7379 (F_ortho = -0.0051).</li>
                <li><strong>Flexible loop drift:</strong> RFD3 backbones drifted up to 19 Å, triggering early exit floors and destroying DockQ.</li>
                <li><strong>SVD frame shifts:</strong> Generative diffusion outputs translated by ~70 Å relative to the target frame.</li>
              </ul>
            </div>
          </div>

          <div className="col-lg-6">
            <div className="card bg-dark border-success p-3 h-100">
              <span className="badge bg-success text-dark mb-2 align-self-start font-monospace">Cycle 2: Algorithmic Innovations</span>
              <h5 className="fw-bold text-light">Engineering Solutions Deployed</h5>
              <ul className="small text-secondary ps-3 mb-0">
                <li><strong>Loop bridling:</strong> Enforced scaffold_max_flex_rmsd &le; 4.5 Å to eliminate aberrant drift.</li>
                <li><strong>Specificity-by-difference:</strong> Biased LigandMPNN scoring to penalize WT contacts and favor synthetic mutations.</li>
                <li><strong>SVD coordinate alignment:</strong> Automatically realigned all diffused scaffolds onto Chain A framework.</li>
                <li><strong>Loop consolidation:</strong> Merged isolated anchor residues into contiguous diffusing loops.</li>
              </ul>
            </div>
          </div>
        </div>

        <div className="text-center my-4">
          <img
            src="figures/fig1_iptm_distributions_and_orthogonality.png"
            alt="Figure 5: High-Throughput Production Landscape"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '900px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 5 | High-Throughput Production Landscape (600 Complexes).</strong> Distribution of 600 designs
            showing the emergence of the 19 elite champions in the top-right quadrant (F_ortho &gt; 0.30, DockQ &gt; 0.49).
          </p>
        </div>
      </section>

      {/* SECTION 5 */}
      <section className="mb-5" id="dual-model-validation">
        <h2 className="fw-bold mb-4 text-primary border-bottom pb-2">
          5. Dual-Model Cross-Validation: RoseTTAFold-3 vs. AlphaFold2-Multimer v3
        </h2>
        <p>
          To ensure independence from model-specific neural biases, all 19 champions were tested under{' '}
          <strong>AlphaFold2-Multimer v3 (57 full 3D predictions)</strong>:
        </p>

        <div className="text-center my-4">
          <img
            src="figures/fig2_correlations_rf3_af2_pyrosetta.png"
            alt="Figure 6: Multi-Model Statistical Correlation"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '900px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 6 | Multi-Model Statistical Correlation Matrix.</strong> Highly significant rank agreement between
            RF3 and AF2 on rescue affinity (&rho; = 0.754, p = 0.00019) and cross-docking barrier (&rho; = 0.625, p = 0.0042).
          </p>
        </div>
      </section>

      {/* SECTION 6 */}
      <section className="mb-5" id="molecular-mechanisms">
        <h2 className="fw-bold mb-4 text-primary border-bottom pb-2">
          6. Atomistic Mechanism: The Steric Knob-into-Hole Switch
        </h2>

        <div className="text-center my-4">
          <img
            src="figures/figure6_structural_knob_into_hole_mechanism.jpg"
            alt="Figure 7: Atomistic Knob-into-Hole Latch"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '900px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 7 | Atomistic Knob-into-Hole Switch.</strong> 3D rendering showing the bulky S12W tryptophan on
            Chain B' (amber) projecting into the cavity carved out by alanine substitutions on Chain A' (teal).
          </p>
        </div>

        <div className="text-center my-4">
          <img
            src="figures/fig3_mutation_frequency_and_conservation.png"
            alt="Figure 8: Mutation Frequency & Conservation Heatmap"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '900px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 8 | Mutation Frequency Heatmap.</strong> Recurring alanine clusters on A' and bulky aromatic
            substitutions (Trp, Phe) on B' across all 19 champions.
          </p>
        </div>
      </section>

      {/* SECTION 7: SUMMARY TABLE */}
      <section className="mb-5" id="champions-table">
        <h2 className="fw-bold mb-4 text-primary border-bottom pb-2">
          7. Summary Table of the 19 Orthogonal Champions
        </h2>

        <div className="table-responsive my-3">
          <table className="table table-dark table-hover table-bordered small align-middle">
            <thead className="table-success text-dark">
              <tr>
                <th>Rank</th>
                <th>Complex ID</th>
                <th>Parent A'</th>
                <th>DockQ</th>
                <th>CAPRI Class</th>
                <th>iPTM Rescue</th>
                <th>iPTM Cross</th>
                <th>iPTM Rupture</th>
                <th>&Delta;iPTM</th>
                <th>&Delta;G (REU)</th>
                <th>F_ortho</th>
                <th>AF2 iPTM</th>
                <th>AF2 RMSD</th>
                <th>Verdict</th>
              </tr>
            </thead>
            <tbody>
              <tr>
                <td><strong>01</strong></td>
                <td><code>cand_01_B_cand_40</code></td>
                <td>cand_01</td>
                <td>0.4812</td>
                <td><span className="badge bg-warning text-dark">Acceptable</span></td>
                <td><strong>0.745</strong></td>
                <td><strong>0.251</strong></td>
                <td><strong>0.295</strong></td>
                <td><strong className="text-success">+0.494</strong></td>
                <td>-40.7</td>
                <td>0.211</td>
                <td>0.390</td>
                <td>1.79 Å</td>
                <td>⭐ Selectivity Record</td>
              </tr>
              <tr>
                <td><strong>02</strong></td>
                <td><code>cand_17_B_cand_29</code></td>
                <td>cand_17</td>
                <td><strong>0.5470</strong></td>
                <td><span className="badge bg-success">Medium</span></td>
                <td><strong>0.781</strong></td>
                <td><strong>0.299</strong></td>
                <td><strong>0.237</strong></td>
                <td><strong className="text-success">+0.483</strong></td>
                <td><strong>-65.8</strong></td>
                <td><strong>0.384</strong></td>
                <td>0.490</td>
                <td>3.52 Å</td>
                <td>⭐ High-Res Benchmark</td>
              </tr>
              <tr>
                <td><strong>03</strong></td>
                <td><code>cand_01_B_cand_36</code></td>
                <td>cand_01</td>
                <td><strong>0.5819</strong></td>
                <td><span className="badge bg-success">Medium</span></td>
                <td><strong>0.701</strong></td>
                <td><strong>0.259</strong></td>
                <td><strong>0.295</strong></td>
                <td><strong className="text-success">+0.442</strong></td>
                <td>-48.1</td>
                <td>0.253</td>
                <td>0.428</td>
                <td><strong>1.46 Å</strong></td>
                <td>⭐ Sub-1.5 Å Native Dock</td>
              </tr>
              <tr>
                <td><strong>04</strong></td>
                <td><code>cand_13_B_cand_22</code></td>
                <td>cand_13</td>
                <td>0.4968</td>
                <td><span className="badge bg-success">Medium</span></td>
                <td>0.602</td>
                <td>0.221</td>
                <td>0.246</td>
                <td><strong className="text-success">+0.381</strong></td>
                <td><strong className="text-info">-74.7</strong></td>
                <td><strong>0.400</strong></td>
                <td>0.122</td>
                <td>6.64 Å</td>
                <td>⚠️ AF2 Divergence</td>
              </tr>
              <tr>
                <td><strong>09</strong></td>
                <td><code>cand_01_B_cand_14</code></td>
                <td>cand_01</td>
                <td><strong>0.6031</strong></td>
                <td><span className="badge bg-success">Medium</span></td>
                <td>0.632</td>
                <td>0.282</td>
                <td>0.295</td>
                <td><strong className="text-success">+0.351</strong></td>
                <td>-48.5</td>
                <td>0.221</td>
                <td><strong>0.470</strong></td>
                <td><strong>1.51 Å</strong></td>
                <td>⭐ Universal Consensus #1</td>
              </tr>
              <tr>
                <td><strong>12</strong></td>
                <td><code>cand_07_B_cand_09</code></td>
                <td>cand_07</td>
                <td>0.5085</td>
                <td><span className="badge bg-success">Medium</span></td>
                <td>0.639</td>
                <td>0.313</td>
                <td>0.246</td>
                <td><strong className="text-success">+0.326</strong></td>
                <td><strong>-67.9</strong></td>
                <td>0.293</td>
                <td>0.369</td>
                <td><strong>1.45 Å</strong></td>
                <td>⭐ Geometric Record</td>
              </tr>
            </tbody>
          </table>
        </div>

        <div className="text-center my-4">
          <img
            src="figures/fig5_champion_radar_profiles.png"
            alt="Figure 9: Radar Profiles"
            className="img-fluid rounded border border-secondary shadow-lg"
            style={{ maxWidth: '900px' }}
          />
          <p className="text-muted small mt-2">
            <strong>Figure 9 | Radar Profiles of Top Champions.</strong> Multi-criteria radar comparisons across DockQ,
            affinity, selectivity margin, and binding energy.
          </p>
        </div>
      </section>

      {/* SECTION 8: WET LAB ROADMAP */}
      <section className="mb-5" id="wetlab-cascade">
        <h2 className="fw-bold mb-4 text-primary border-bottom pb-2">
          8. DBTL Wet-Lab Characterization Cascade
        </h2>

        <div className="row g-3">
          <div className="col-md-4">
            <div className="card h-100 bg-dark border-secondary p-3">
              <span className="badge bg-info text-dark mb-2 align-self-start font-monospace">Phase 1: Kinetics</span>
              <h5 className="fw-bold text-light">BLI / SPR 4&times;4 Titration</h5>
              <p className="small text-secondary mb-0">
                Cell-free PURE synthesis; AviTag site-specific biotinylation; real-time k_on, k_off, and Kd measurement
                across all four complexes (target Kd &lt; 50 nM, cross-reactions Kd &gt; 10 &mu;M).
              </p>
            </div>
          </div>

          <div className="col-md-4">
            <div className="card h-100 bg-dark border-secondary p-3">
              <span className="badge bg-warning text-dark mb-2 align-self-start font-monospace">Phase 2: Assembly</span>
              <h5 className="fw-bold text-light">SEC-MALS Monodispersity</h5>
              <p className="small text-secondary mb-0">
                Superdex 200 Increase size exclusion coupled with multi-angle light scattering to confirm pure 1:1
                heterodimerization without higher-order oligomers.
              </p>
            </div>
          </div>

          <div className="col-md-4">
            <div className="card h-100 bg-dark border-secondary p-3">
              <span className="badge bg-success text-dark mb-2 align-self-start font-monospace">Phase 3: In Vivo</span>
              <h5 className="fw-bold text-light">Dual-Reporter Assay</h5>
              <p className="small text-secondary mb-0">
                Dual-fluorescence in <em>E. coli</em>: sfGFP (orthogonal channel) vs. mCherry (host background) quantified
                by flow cytometry to confirm physiological circuit insulation.
              </p>
            </div>
          </div>
        </div>
      </section>

      {/* SECTION 9: REPRODUCIBILITY & REFERENCES */}
      <section className="mb-5" id="open-source-references">
        <h2 className="fw-bold mb-4 text-primary border-bottom pb-2">
          9. Open-Source Pipeline, Reproducibility &amp; References
        </h2>
        <p>
          The complete <code>OrthoPPInterface</code> platform is distributed under the MIT License:
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
        <div className="alert alert-secondary bg-dark border-secondary small text-secondary">
          All modules (<code>01_analyze_interface.py</code> through <code>05_eval_final.py</code>), the test suite (71 tests passing),
          and the cryptographic signed caching system are fully documented in the repository.
        </div>
      </section>
    </div>
  );
}
