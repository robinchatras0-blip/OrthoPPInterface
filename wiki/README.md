# Package Officiel Wiki iGEM 2026 (Équipe Evry-Paris-Saclay)

Ce répertoire contient l'ensemble des fichiers nettoyés, structurés et rigoureusement validés pour l'intégration de la page **Model** du wiki iGEM 2026.

## Structure du Répertoire

```
wiki/
├── figures/                                    # Figures PNG haute résolution générées à partir des données réelles
│   ├── figure1_triad_2d_vs_3d.png             # Architecture 2D (ARN) vs 3D (PPI + ARNr 16S)
│   ├── figure2_universal_pathways.png         # Applications universelles aux voies synthétiques
│   ├── figure3_cross_reactivity_matrix.png    # Matrice de spécificité 4-états & définition F_ortho
│   ├── figure4_ode_competition_kinetics.png   # Cinétique ODE et sélectivité de traduction
│   ├── fig1_iptm_distributions_and_orthogonality.png # Paysage empirique des 600 complexes & 19 champions
│   ├── fig2_correlations_rf3_af2_pyrosetta.png # Matrice de corrélation statistique RF3 vs AF2 vs PyRosetta
│   ├── fig3_mutation_frequency_and_conservation.png # Fréquence mutationnelle & conservation
│   ├── fig4_plddt_profiles_per_residue.png    # Profils pLDDT résidu par résidu
│   └── fig5_champion_radar_profiles.png       # Profils radar multicritères des super-champions
│
├── src_gitlab/                                 # Fichiers sources prêts à copier vers le dépôt GitLab iGEM
│   ├── model.tsx                              # Composant React 19 / Bootstrap 5 (pour src/contents/model.tsx)
│   └── model.bib                              # Références BibTeX complètes (pour src/references/model.bib)
│
├── model.md                                   # Version Markdown complète, épurée et scientifiquement justifiée
├── index.html                                 # Prévisualisation interactive standalone locale (KaTeX, Mermaid)
└── README.md                                  # Ce guide d'organisation
```

## Données et Modélisation Intégrées

1. **Architecture du Pipeline en 5 Modules :**
   - **Module 1 :** Cartographie topologique de l'interface (rayon de contact 6.0 Å, rognage de B à 15.0 Å pour préserver l'enveloppe stérique sous contrainte VRAM), consolidation des boucles (`rescue_min_fixed_segment: 3`) et guidage de l'origine de diffusion par les hotspots PPI (+10 Å le long de la normale extérieure).
   - **Module 2 :** Génération de la rupture $A'$ en état apo (`rupture_apo_design: true`) pour découpler la séquence de tout amarrage favorable avec $B_{\text{WT}}$, biais électrostatique et stérique LigandMPNN, compactage PyRosetta Ref2015 FastRelax.
   - **Module 3 :** Filtres éliminatoires rapides : repliement monomère de stabilité ($\text{pLDDT} \ge 80$, $\text{RMSD} \le 2.0$ Å) et rupture confirmée avec $B_{\text{WT}}$ ($\text{iPTM} \le 0.35$). Sélection de diversité par clustering Hamming le plus éloigné.
   - **Module 4 :** Co-adaptation de $B'$ : recalage rigide SVD des coordonnées diffusées, filtrage des squelettes (`scaffold_max_flex_rmsd \le 4.5\text{ \AA}`), biais de séquence LigandMPNN par spécificité différentielle (récompense des contacts mutés, pénalisation des contacts sauvages), reconstruction et relaxation PyRosetta contre $A'$ fixé.
   - **Module 5 :** Évaluation globale sous régime d'alignement unifié avec MSA masquées à l'interface (plafond sauvage = 0.5446) pour éliminer les biais co-évolutifs ancestraux, calcul d'énergie libre PyRosetta, score CAPRI DockQ et exportation SQLite.

2. **Résultats de Production (Campagne de 600 Complexes) :**
   - 600 complexes all-atom explorés sur 12 familles de motifs $A'$.
   - 70.3% (422/600) conservent la géométrie native ($\text{DockQ} \ge 0.23$), dont 20 complexes de Haute Qualité CAPRI ($\text{DockQ} \ge 0.80$, pic à 0.9052).
   - 19 Champions d'Élite Stricts validés : rupture totale du sauvage ($\text{iPTM} \le 0.35$), affinité de sauvetage élevée ($\text{iPTM} \ge 0.60 - 0.78$), record de sélectivité $\Delta\text{iPTM} = +0.4936$ (`cand_01_B_cand_40`), énergies PyRosetta jusqu'à $-74.69\text{ REU}$.

3. **Contre-Validation Indépendante AlphaFold2-Multimer v3 :**
   - 57 prédictions atomiques 3D complètes ($A' \cdot B'$, $A' \cdot B_{\text{WT}}$, $A_{\text{WT}} \cdot B'$).
   - Très forte corrélation de rang avec RoseTTAFold-3 All-Atom ($\rho = 0.754$, $p = 0.00019$).
   - Consensus fort sur la lignée `cand_01` maintenant un amarrage natif sous les 1.6 Å tout en expulsant le partenaire sauvage à 13.26 Å.

4. **Modèle Cinétique ODE & Cascade Expérimentale :**
   - Simulation numérique ODE démontrant qu'une marge $F_{\text{ortho}} \ge +0.20$ garantit $>85\%$ d'isolation traductionnelle in vivo.
   - Protocole de validation wet-lab à 3 niveaux : cinétique BLI/SPR 4x4, SEC-MALS 1:1, et essai double rapporteur in vivo sfGFP/mCherry.
