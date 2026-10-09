# Package Officiel Wiki iGEM 2026 (Équipe Evry-Paris-Saclay)

Ce répertoire contient l'ensemble des fichiers nettoyés, structurés et validés pour l'intégration de la page **Model** du wiki iGEM 2026.

## Structure du Répertoire

```
wiki/
├── figures/                                    # Figures PNG/JPG haute résolution (300 DPI, dark/light mode)
│   ├── figure1_triad_2d_vs_3d.png             # Architecture 2D (ARN) vs 3D (PPI + ARNr 16S)
│   ├── figure2_universal_pathways.png         # Applications universelles aux voies synthétiques
│   ├── figure3_cross_reactivity_matrix.png    # Matrice de spécificité 4-états & définition F_ortho
│   ├── figure4_ode_competition_kinetics.png   # Cinétique ODE et sélectivité de traduction
│   ├── fig1_iptm_distributions_and_orthogonality.png # Paysage empirique des 600 complexes & 19 champions
│   ├── fig2_correlations_rf3_af2_pyrosetta.png # Matrice de corrélation statistique RF3 vs AF2 vs PyRosetta
│   ├── fig3_mutation_frequency_and_conservation.png # Fréquence mutationnelle & conservation
│   ├── fig4_plddt_profiles_per_residue.png    # Profils pLDDT résidu par résidu
│   ├── fig5_champion_radar_profiles.png       # Profils radar multicritères des super-champions
│   └── figure6_structural_knob_into_hole_mechanism.jpg # Mécanisme atomique 3D du verrou Trp12
│
├── src_gitlab/                                 # Fichiers sources prêts à copier vers le dépôt GitLab iGEM
│   ├── model.tsx                              # Composant React 19 / Bootstrap 5 (pour src/contents/model.tsx)
│   └── model.bib                              # Références BibTeX complètes (pour src/references/model.bib)
│
├── model.md                                   # Version Markdown complète mise à jour (Cycles DBTL 1 & 2)
├── index.html                                 # Prévisualisation interactive standalone locale (KaTeX, Mermaid)
└── README.md                                  # Ce guide d'organisation
```

## Données et Modélisation Intégrées

1. **Cycle DBTL 1 (Diagnostic) :**
   - Échecs initiaux documentés : dérive des boucles RFD3 (>10–19 Å), amarrage croisé de B' avec $A_{\text{WT}}$, plancher d'early exit.

2. **Innovations Algorithmiques :**
   - Bridage conformationnel des squelettes (`scaffold_max_flex_rmsd \le 4.5\text{ \AA}`).
   - Biais LigandMPNN différentiel de spécificité (pénalisation des résidus complémentaires au sauvage).
   - Recalage rigide SVD des motifs diffusés (correction du décalage de ~70 Å).
   - Consolidation des boucles flexibles et ancrage d'origine par hotspots PPI.

3. **Cycle DBTL 2 (Production Haute Résolution & Percée des Champions) :**
   - 600 complexes explorés sur 12 familles de motifs $A'$.
   - 70.3% (422/600) conservent la géométrie native ($\text{DockQ} \ge 0.23$), dont 20 complexes de Haute Qualité CAPRI ($\text{DockQ} \ge 0.80$, pic à 0.9052).
   - 19 Champions d'Élite Stricts validés : rupture totale du sauvage ($\text{iPTM} \le 0.35$), affinité de sauvetage élevée ($\text{iPTM} \ge 0.60 - 0.78$), record de sélectivité $\Delta\text{iPTM} = +0.4936$ (`cand_01_B_cand_40`), énergies PyRosetta jusqu'à $-74.69\text{ REU}$.

4. **Contre-Validation Indépendante AlphaFold2-Multimer v3 :**
   - 57 prédictions atomiques 3D complètes.
   - Très forte corrélation de rang avec RoseTTAFold-3 All-Atom ($\rho = 0.754$, $p < 0.001$).
   - Élucidation du mécanisme atomique de verrou stérique (insertion du tryptophane `S12W` dans la cavité alanine).
