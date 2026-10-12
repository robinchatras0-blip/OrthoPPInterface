# OrthoPPInterface — mémoire projet

Pipeline de redesign d'interface protéine-protéine **orthogonale** : à partir d'un complexe natif
A·B, concevoir A' et B' tels que A'·B' se lie fortement alors que A'·B_WT et A_WT·B' ne se lient pas.
C'est un problème de **spécificité**, pas d'affinité : toute suggestion qui améliore seulement
l'affinité de A'·B' sans creuser l'écart avec les paires croisées rate la cible.

## Commandes

```bash
conda activate orthoppinterface          # toujours activer, ne jamais appeler le python de l'env directement
python run_pipeline.py --run_name <nom>              # lance le pipeline via pipeline/run_pipeline.py
python pipeline/run_pipeline.py --run_name <nom> --steps 4-5  # sous-ensemble de modules direct
python pipeline/src/04_generate_B_prime.py --config pipeline/config.yaml --force   # module isolé
python pipeline/src/calibrate_predictor.py --config pipeline/config.yaml --analysis_dir runs/<run>/01_analysis
pytest pipeline/tests                     # rapide, pas de GPU : les moteurs sont remplacés par des fakes (71 tests)
```

## Architecture du Dépôt

- `pipeline/` : Cœur algorithmique de bio-ingénierie (`pipeline/src/`, `pipeline/tests/`, `pipeline/bin/`, `pipeline/config.yaml`).
- `wiki/` : Package complet du wiki iGEM 2026 (`model.md`, `index.html`, `src_gitlab/model.tsx`, figures).
- `dossier_champions/` : Fiches et structures des 19 champions d'élite all-atom.
- `runs/` : Dossier unifié de stockage de tous les résultats d'expériences (git-ignoré).

| Module | Fichier | Rôle |
| --- | --- | --- |
| 1 | `pipeline/src/01_analyze_interface.py` | Résidus d'interface, hotspots, positions mutables, contigs RFD3, consolidation boucles |
| 2 | `pipeline/src/02_generate_A_prime.py` | Diffusion RFD3 apo + LigandMPNN → candidats A' (rupture) |
| 3 | `pipeline/src/03_filter_A_prime.py` | Fail-fast RF3 : stabilité monomère, rupture A'·B_WT, plafond WT |
| 4 | `pipeline/src/04_generate_B_prime.py` | Co-adaptation RFD3 multi-scaffolds (recalage SVD) + LigandMPNN biaisé spécificité → B' |
| 5 | `pipeline/src/05_eval_final.py` | 4 folds RF3 + énergie PyRosetta Ref2015 + DockQ CAPRI + F_ortho + export CSV/SQLite |

Modules transverses : `design_utils.py` (séquences/structures), `folding_engine.py` (RF3, MSA, cache),
`energy_engine.py` + `interface_energy.py` (PyRosetta), `scoring.py` (F_ortho, ipSAE), `dockq.py`,
`export_styled_reports.py`.

Chaque expérience vit dans `runs/<run_name>/` (git-ignoré) avec `config_used.yaml` archivé.
Sorties clés : `runs/<run>/05_final_eval/orthogonality_scores.csv`, `results.db`, `coherence.json`,
`runs/<run>/03_fail_fast/wt_control.json`.

## Invariants à ne jamais casser

- **Même régime d'information pour les quatre folds.** A'·B', A_WT·B', A'·B_WT et le plafond WT/WT
  doivent être prédits sous les mêmes masques MSA. Sinon F_ortho mesure une asymétrie de MSA, pas une
  liaison. Toute modification de `msa_mask_scope`/`msa_mask_mode`/`negative_msa_regime` s'applique aux
  quatre.
- **Les valeurs manquantes sont NaN.** Jamais de valeur par défaut substituée à une mesure ratée, ni
  dans le code ni dans les analyses.
- **Le cache RF3 est signé** (séquences + contenu MSA + paramètres d'inférence). Ne pas contourner la
  signature ; un changement de régime doit invalider le cache.
- **A' est designé seul** (`rupture_apo_design: true`). Avec B_WT visible, LigandMPNN empaquette A'
  contre B_WT et fait du design positif pour le partenaire sauvage.
- **RFD3 recentre sa sortie** (~70 Å). Tout scaffold doit être réaligné sur A' via le framework fixe
  avant la moindre mesure.
- **LigandMPNN n'écrit pas les chaînes latérales** des positions redesignées, et sort du CIF. Packing
  PyRosetta obligatoire (`energy.pack_designs`) avant toute mesure géométrique ou énergétique.
- Pas de fallback silencieux : si LigandMPNN ne produit rien, le module s'arrête ; on ne retombe pas
  sur le backbone RFD3 brut.

## Bruit et interprétation des métriques

- iPTM varie de ±0.05–0.1 entre runs. Ne jamais conclure sur une différence d'iPTM inférieure à ~0.1
  entre deux designs. `iptm_rescue_std` est loggé, s'en servir.
- F_energy a un écart-type ≈ 0.05 entre seeds pour un étalement de 0.22 entre designs. Utile pour
  classer et pour un shortlist répété, pas pour trancher une petite différence.
- Les seuils `min_binding_fraction: 0.5` et `f_energy_min: 0.0` sont **provisoires**, non calibrés.
  Les traiter comme tels dans toute analyse.
- `ipsae_*` et `f_ipsae` sont **diagnostiques uniquement** : suivis dans `coherence.json`, pas encore
  dans les critères pass/fail. Ne pas les promouvoir en seuil sans démonstration d'accord.
- RF3 et l'énergie échouent différemment : un design doit satisfaire les deux. Un désaccord fort entre
  les deux sur un design est un signal, pas un bug à masquer.

## Conventions de code

- Python 3.10 côté pipeline (3.12 dans l'env des moteurs). Biopython, pandas, numpy, pyyaml.
- Les modules de `src/` commencent par un chiffre : ils s'importent par `tests/conftest.py::load_module`,
  pas par `import`.
- Config lue en UTF-8 (BOM toléré) ; sections `pipeline`/`thresholds` manquantes → erreur explicite,
  jamais de fallback silencieux.
- Chemins normalisés en `/` (compat WSL/Windows).
- Les tests couvrent la plomberie (câblage des modules, parsing, cache, numérotation DockQ) avec des
  moteurs factices. Ils ne valident **pas** la qualité scientifique — ne jamais présenter un `pytest`
  vert comme une preuve qu'un changement méthodologique est correct.

## Règles de travail

- Avant de modifier un seuil ou un paramètre de `config.yaml`, dire quelle mesure de quel run le
  justifie. Un seuil changé sans donnée est un aveuglement, pas un réglage.
- Ne jamais lancer un module GPU (`--steps 2`, `4`, `5`) sans demande explicite : ce sont des dizaines
  de minutes de calcul.
- Pour toute modification de `src/`, faire tourner `pytest tests` et le dire.
- `runs/`, `data/outputs/`, `*.pdb`, `*.a3m`, `*.pt` sont git-ignorés : ne jamais proposer de les
  committer.
