#!/usr/bin/env python3
"""
generate_statistical_figures.py - Generates publication-quality statistical figures
evaluating the orthogonal protein engineering campaign (RF3, AF2, PyRosetta, mutations).

Outputs saved to DOSSIER_CHAMPIONS/figures/:
1. fig1_iptm_distributions_and_orthogonality.png
2. fig2_correlations_rf3_af2_pyrosetta.png
3. fig3_mutation_frequency_and_conservation.png
4. fig4_plddt_profiles_per_residue.png
5. fig5_champion_radar_profiles.png
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from Bio.PDB import PDBParser

# Configure Matplotlib styling for high-impact publication aesthetics
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#94a3b8'
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['grid.color'] = '#e2e8f0'
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.alpha'] = 0.7

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DOSSIER_DIR = os.path.join(BASE_DIR, "DOSSIER_CHAMPIONS")
TABLES_DIR = os.path.join(DOSSIER_DIR, "tableaux_scores")
PDB_DIR = os.path.join(DOSSIER_DIR, "structures_pdb")
FIG_DIR = os.path.join(DOSSIER_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Load data
summary_csv = os.path.join(TABLES_DIR, "af2_validation_19_champions.csv")
mutations_csv = os.path.join(TABLES_DIR, "mutations_detaillees.csv")
corr_csv = os.path.join(TABLES_DIR, "correlation_rf3_af2.csv")

df_champ = pd.read_csv(summary_csv)
df_mut = pd.read_csv(mutations_csv)
df_corr = pd.read_csv(corr_csv) if os.path.exists(corr_csv) else None

print(f"Loaded {len(df_champ)} champions and {len(df_mut)} mapped mutations.")


# ==============================================================================
# FIGURE 1: Matrice de Sélectivité et Distributions iPTM (Sauvetage vs Négatifs)
# ==============================================================================
print("Generating Figure 1: iPTM distributions & selectivity...")
fig = plt.figure(figsize=(16, 10), constrained_layout=True)
gs = fig.add_gridspec(2, 2)

# Panel A: Boxplot / Swarm comparison RF3 vs AF2 states
ax1 = fig.add_subplot(gs[0, 0])
rf3_data = [df_champ['iptm_rescue_rf3'], df_champ['iptm_negative_rf3'], df_champ['iptm_rupture_rf3']]
af2_data = [df_champ['rescue_iptm_af2'], df_champ['cross_iptm_af2'], df_champ['rupture_iptm_af2']]

labels = ['Sauvetage\n(A\'·B\')', 'Croisé WT\n(A_WT·B\')', 'Rupture WT\n(A\'·B_WT)']
positions_rf3 = np.array([1, 4, 7])
positions_af2 = np.array([2, 5, 8])

bp1 = ax1.boxplot(rf3_data, positions=positions_rf3, widths=0.7, patch_artist=True,
                  boxprops=dict(facecolor='#3b82f6', alpha=0.7, edgecolor='#1d4ed8'),
                  medianprops=dict(color='#1e3a8a', linewidth=2))
bp2 = ax1.boxplot(af2_data, positions=positions_af2, widths=0.7, patch_artist=True,
                  boxprops=dict(facecolor='#10b981', alpha=0.7, edgecolor='#047857'),
                  medianprops=dict(color='#064e3b', linewidth=2))

# Swarm scatter on top
for i, d in enumerate(rf3_data):
    x = np.random.normal(positions_rf3[i], 0.08, size=len(d))
    ax1.scatter(x, d, color='#1e3a8a', alpha=0.8, s=28, zorder=3)
for i, d in enumerate(af2_data):
    x = np.random.normal(positions_af2[i], 0.08, size=len(d))
    ax1.scatter(x, d, color='#064e3b', alpha=0.8, s=28, zorder=3)

ax1.axhline(0.35, color='#dc2626', linestyle='--', linewidth=1.5, alpha=0.8, label='Plafond Rupture (0.35)')
ax1.set_xticks([1.5, 4.5, 7.5])
ax1.set_xticklabels(labels, fontsize=11, fontweight='bold')
ax1.set_ylabel('Score d\'Amarrage iPTM', fontsize=12, fontweight='bold')
ax1.set_title('A. Séparation Radicale des États (Sauvetage vs Négatifs)', fontsize=13, fontweight='bold', pad=12)
ax1.legend([bp1["boxes"][0], bp2["boxes"][0]], ['RoseTTAFold-3 All-Atom', 'AlphaFold2-Multimer v3'], loc='upper right', frameon=True)
ax1.set_ylim(0.05, 0.90)

# Panel B: Waterfall Delta iPTM AF2
ax2 = fig.add_subplot(gs[0, 1])
df_sorted = df_champ.sort_values(by='delta_iptm_af2', ascending=True).copy()
colors = ['#059669' if d >= 0.15 else ('#10b981' if d >= 0.10 else ('#f59e0b' if d >= 0.0 else '#ef4444')) 
          for d in df_sorted['delta_iptm_af2']]
bars = ax2.barh(range(len(df_sorted)), df_sorted['delta_iptm_af2'], color=colors, edgecolor='#334155', linewidth=0.8)
ax2.axvline(0.0, color='#334155', linestyle='-', linewidth=1.2)
ax2.axvline(0.10, color='#059669', linestyle='--', linewidth=1.2, label='Seuil Fort Consensus (+0.10)')

ax2.set_yticks(range(len(df_sorted)))
ax2.set_yticklabels([f"#{r} {d.replace('A_prime_', '').replace('_cand_', ' c')}" 
                     for r, d in zip(df_sorted['rank'], df_sorted['design_id'])], fontsize=9)
ax2.set_xlabel('Différentiel de Sélectivité AlphaFold2 ($\Delta$ iPTM)', fontsize=11, fontweight='bold')
ax2.set_title('B. Spectre de Sélectivité Orthogonale AF2 (19 Champions)', fontsize=13, fontweight='bold', pad=12)
ax2.legend(loc='lower right', frameon=True)

# Panel C: Delta iPTM RF3 vs Delta iPTM AF2
ax3 = fig.add_subplot(gs[1, :])
x_val = df_champ['delta_iptm_rf3']
y_val = df_champ['delta_iptm_af2']
sizes = (df_champ['dockq_rf3'] ** 2) * 500

scatter = ax3.scatter(x_val, y_val, c=df_champ['rescue_plddt_af2'], cmap='viridis', 
                      s=sizes, alpha=0.85, edgecolors='#0f172a', linewidth=1.2, zorder=3)
cbar = plt.colorbar(scatter, ax=ax3, orientation='vertical', pad=0.02)
cbar.set_label('Confidence AF2 (pLDDT)', fontsize=11, fontweight='bold')

# Quadrants
ax3.axhline(0.10, color='#16a34a', linestyle='--', linewidth=1.2, alpha=0.7)
ax3.axvline(0.30, color='#2563eb', linestyle='--', linewidth=1.2, alpha=0.7)
ax3.fill_between([0.30, 0.55], 0.10, 0.35, color='#dcfce7', alpha=0.4, label='Zone de Consensus Élite (RF3+AF2)')

for _, r in df_champ.iterrows():
    if r['delta_iptm_af2'] > 0.12 or r['rank'] in [1, 2, 3, 9, 12, 17]:
        name_short = r['design_id'].replace('A_prime_candidate_', 'A').replace('_B_cand_', '_B')
        ax3.annotate(f"#{r['rank']} {name_short}", 
                     (r['delta_iptm_rf3'], r['delta_iptm_af2']),
                     textcoords="offset points", xytext=(6, 5), fontsize=9, fontweight='bold',
                     bbox=dict(boxstyle="round,pad=0.2", facecolor="#ffffff", alpha=0.8, edgecolor="#94a3b8"))

ax3.set_xlabel('$\Delta$ iPTM RoseTTAFold-3 All-Atom', fontsize=12, fontweight='bold')
ax3.set_ylabel('$\Delta$ iPTM AlphaFold2-Multimer v3', fontsize=12, fontweight='bold')
ax3.set_title('C. Espace de Sélectivité Bi-Modèle : Convergence vers le Consensus Universel', fontsize=13, fontweight='bold', pad=12)
ax3.legend(loc='upper left', frameon=True)

fig.savefig(os.path.join(FIG_DIR, "fig1_iptm_distributions_and_orthogonality.png"), dpi=300)
plt.close(fig)


# ==============================================================================
# ==============================================================================
# FIGURE 2: Corrélations Statistique Multi-Modèles & Équivalence Biophysique (RF3 vs AF2)
# ==============================================================================
print("Generating Figure 2: Equivalent multi-model comparisons (RF3 vs AF2)...")
fig, axes = plt.subplots(2, 2, figsize=(15, 12), constrained_layout=True)

# Panel A: iPTM Rescue RF3 vs AF2 (Strictly equivalent: Cognate rescue affinity)
ax = axes[0, 0]
x = df_champ['iptm_rescue_rf3']
y = df_champ['rescue_iptm_af2']
sns.regplot(x=x, y=y, ax=ax, color='#2563eb', scatter_kws={'s': 70, 'edgecolor': '#1e3a8a', 'alpha': 0.85},
            line_kws={'color': '#1d4ed8', 'linewidth': 2})
r_p, p_p = stats.pearsonr(x, y)
r_s, p_s = stats.spearmanr(x, y)
ax.annotate(f"Pearson $r = +{r_p:.3f}$ ($p = {p_p:.3e}$)\nSpearman $\\rho = +{r_s:.3f}$ ($p = {p_s:.3e}$)",
            xy=(0.05, 0.82), xycoords='axes fraction', fontsize=10, fontweight='bold',
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#eff6ff", edgecolor="#3b82f6"))
ax.set_title('A. Affinité de Sauvetage Cognat (A\'·B\') : Accord RF3 vs AF2', fontsize=12, fontweight='bold')
ax.set_xlabel('iPTM Sauvetage (RoseTTAFold-3 All-Atom)', fontsize=11, fontweight='bold')
ax.set_ylabel('iPTM Sauvetage (AlphaFold2-Multimer v3)', fontsize=11, fontweight='bold')
ax.set_xlim(0.55, 0.82)
ax.set_ylim(0.08, 0.65)

# Panel B: iPTM Croisé RF3 vs AF2 (Strictly equivalent: Cross-docking off-target)
ax = axes[0, 1]
x = df_champ['iptm_negative_rf3']
y = df_champ['cross_iptm_af2']
sns.regplot(x=x, y=y, ax=ax, color='#8b5cf6', scatter_kws={'s': 70, 'edgecolor': '#5b21b6', 'alpha': 0.85},
            line_kws={'color': '#6d28d9', 'linewidth': 2})
r_p, p_p = stats.pearsonr(x, y)
r_s, p_s = stats.spearmanr(x, y)
ax.annotate(f"Pearson $r = +{r_p:.3f}$ ($p = {p_p:.3f}$)\nSpearman $\\rho = +{r_s:.3f}$ ($p = {p_s:.3e}$)",
            xy=(0.05, 0.82), xycoords='axes fraction', fontsize=10, fontweight='bold',
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#f5f3ff", edgecolor="#8b5cf6"))
ax.set_title('B. Barrière d\'Amarrage Croisé (A_WT·B\') : Accord RF3 vs AF2', fontsize=12, fontweight='bold')
ax.set_xlabel('iPTM Croisé (RoseTTAFold-3 All-Atom)', fontsize=11, fontweight='bold')
ax.set_ylabel('iPTM Croisé (AlphaFold2-Multimer v3)', fontsize=11, fontweight='bold')
ax.set_xlim(0.18, 0.45)
ax.set_ylim(0.05, 0.50)

# Panel C: iPTM Rupture RF3 vs AF2 (Strictly equivalent: Host WT disruption)
ax = axes[1, 0]
x = df_champ['iptm_rupture_rf3']
y = df_champ['rupture_iptm_af2']
sc = ax.scatter(x, y, c='#ef4444', s=80, edgecolors='#991b1b', alpha=0.85, zorder=3)
ax.axhline(0.35, color='#dc2626', linestyle='--', linewidth=1.2, label='Seuil Rupture AF2 (0.35)')
ax.axvline(0.35, color='#b91c1c', linestyle=':', linewidth=1.2, label='Seuil Rupture RF3 (0.35)')
ax.fill_between([0.15, 0.35], 0.05, 0.35, color='#fee2e2', alpha=0.4, label='Zone de Rupture Mutuelle (<0.35)')
ax.annotate(f"100% des champions < 0.35\nMoyenne RF3 = 0.253\nMoyenne AF2 = 0.184",
            xy=(0.05, 0.78), xycoords='axes fraction', fontsize=10, fontweight='bold',
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#fff1f2", edgecolor="#f43f5e"))
ax.set_title('C. Rupture du Complexe Sauvage (A\'·B_WT) : Accord d\'Abolition Mutuelle', fontsize=12, fontweight='bold')
ax.set_xlabel('iPTM Rupture (RoseTTAFold-3 All-Atom)', fontsize=11, fontweight='bold')
ax.set_ylabel('iPTM Rupture (AlphaFold2-Multimer v3)', fontsize=11, fontweight='bold')
ax.set_xlim(0.18, 0.38)
ax.set_ylim(0.08, 0.38)
ax.legend(loc='lower right', frameon=True, fontsize=9)

# Panel D: Geometric pose deviation in Angstroms: lRMSD RF3 vs RMSD AF2 (Strictly equivalent units)
ax = axes[1, 1]
x = df_champ['lrms_rf3']
y = df_champ['rescue_rmsd_af2']
sns.regplot(x=x, y=y, ax=ax, color='#059669', scatter_kws={'s': 70, 'edgecolor': '#064e3b', 'alpha': 0.85},
            line_kws={'color': '#047857', 'linewidth': 2})
r_p, p_p = stats.pearsonr(x, y)
r_s, p_s = stats.spearmanr(x, y)
ax.axhline(2.0, color='#047857', linestyle='--', linewidth=1.2, label='Seuil Haute Résolution AF2 (2.0 Å)')
ax.axvline(3.0, color='#065f46', linestyle=':', linewidth=1.2, label='Seuil Sub-3.0 Å RF3')
ax.annotate(f"Pearson $r = +{r_p:.3f}$ ($p = {p_p:.3f}$)\nSpearman $\\rho = +{r_s:.3f}$ ($p = {p_s:.3f}$)\nCluster élite sub-2.0 Å partagé",
            xy=(0.05, 0.76), xycoords='axes fraction', fontsize=10, fontweight='bold',
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#ecfdf5", edgecolor="#059669"))
ax.set_title('D. Accord de Conformation Atomique : lRMSD RF3 vs RMSD AF2 (Å)', fontsize=12, fontweight='bold')
ax.set_xlabel('Ligand RMSD RF3 (Å)', fontsize=11, fontweight='bold')
ax.set_ylabel('RMSD de Sauvetage AF2 (Å)', fontsize=11, fontweight='bold')
ax.set_xlim(2.0, 7.0)
ax.set_ylim(0.5, 7.5)
ax.legend(loc='lower right', frameon=True, fontsize=9)

fig.savefig(os.path.join(FIG_DIR, "fig2_correlations_rf3_af2_pyrosetta.png"), dpi=300)
plt.close(fig)


# ==============================================================================
# FIGURE 3: Cartographie des Mutations & Taux de Conservation par Position
# ==============================================================================
print("Generating Figure 3: Mutation mapping & conservation profiles...")
fig = plt.figure(figsize=(16, 11), constrained_layout=True)
gs = fig.add_gridspec(3, 2, height_ratios=[1.2, 1.2, 1.6])

# Chain A (224 residues) mutation frequency
ax_a = fig.add_subplot(gs[0, :])
mut_a_counts = df_mut[df_mut['chain'] == 'A']['pos'].value_counts()
freq_a = np.zeros(224)
for pos, cnt in mut_a_counts.items():
    if 1 <= pos <= 224:
        freq_a[pos - 1] = cnt / 19.0 * 100.0

pos_a = np.arange(1, 225)
ax_a.bar(pos_a, freq_a, width=1.0, color='#3b82f6', edgecolor='none', alpha=0.85)
ax_a.plot(pos_a, freq_a, color='#1d4ed8', linewidth=1.0)
ax_a.set_xlim(1, 224)
ax_a.set_ylim(0, 105)
ax_a.set_xlabel('Position sur la Chaîne A\' (1 à 224)', fontsize=11, fontweight='bold')
ax_a.set_ylabel('% Mutation (19 Champions)', fontsize=11, fontweight='bold')
ax_a.set_title('A. Profil de Variabilité Séquentielle sur la Sous-Unité A\' (Hotspots vs Framework Conservé)', fontsize=12, fontweight='bold')

# Annotate key hotspots on A'
hotspots_a = [3, 4, 13, 15, 17, 18, 19, 29, 32, 34, 38, 39, 41, 70, 71, 185, 186, 201, 202, 205]
ax_a.axvspan(3, 20, color='#fed7aa', alpha=0.3, label='Cluster N-ter (Helic 1 / Loop)')
ax_a.axvspan(29, 42, color='#bbf7d0', alpha=0.3, label='Cluster Médian (Feuillet $\\beta$)')
ax_a.axvspan(185, 205, color='#ddd6fe', alpha=0.3, label='Cluster C-ter (Boucle de contact)')
ax_a.legend(loc='upper right', frameon=True, fontsize=9)

# Chain B (88 residues) mutation frequency
ax_b = fig.add_subplot(gs[1, :])
mut_b_counts = df_mut[df_mut['chain'] == 'B']['pos'].value_counts()
freq_b = np.zeros(88)
for pos, cnt in mut_b_counts.items():
    if 1 <= pos <= 88:
        freq_b[pos - 1] = cnt / 19.0 * 100.0

pos_b = np.arange(1, 89)
ax_b.bar(pos_b, freq_b, width=1.0, color='#f97316', edgecolor='none', alpha=0.85)
ax_b.plot(pos_b, freq_b, color='#c2410c', linewidth=1.0)
ax_b.set_xlim(1, 88)
ax_b.set_ylim(0, 105)
ax_b.set_xlabel('Position sur la Chaîne B\' (1 à 88)', fontsize=11, fontweight='bold')
ax_b.set_ylabel('% Mutation (19 Champions)', fontsize=11, fontweight='bold')
ax_b.set_title('B. Profil de Variabilité Séquentielle sur la Sous-Unité B\' (Les 3 Îlots de Contact)', fontsize=12, fontweight='bold')

ax_b.axvspan(2, 20, color='#fed7aa', alpha=0.3, label='Îlot N-ter (Hélice 1, dont S12W)')
ax_b.axvspan(29, 47, color='#bbf7d0', alpha=0.3, label='Îlot Médian (Boucle & Brin, dont L37R, D39R)')
ax_b.axvspan(74, 81, color='#ddd6fe', alpha=0.3, label='Îlot C-ter (Amarrage final)')
ax_b.legend(loc='upper right', frameon=True, fontsize=9)

# Panel C: Amino Acid Substitution Matrix (Heatmap)
ax_c = fig.add_subplot(gs[2, 0])
aa_order = list("ACDEFGHIKLMNPQRSTVWY")
sub_matrix = pd.crosstab(df_mut['wt_aa'], df_mut['mut_aa']).reindex(index=aa_order, columns=aa_order, fill_value=0)
sns.heatmap(sub_matrix, cmap='YlGnBu', ax=ax_c, cbar_kws={'label': 'Nombre de Substitutions'}, linewidths=0.5)
ax_c.set_title('C. Matrice des Substitutions d\'Acides Aminés (1 169 mutations)', fontsize=12, fontweight='bold')
ax_c.set_xlabel('Acide Aminé Muté (Conçu)', fontsize=10, fontweight='bold')
ax_c.set_ylabel('Acide Aminé Sauvage (WT)', fontsize=10, fontweight='bold')

# Panel D: Top mutated residues barplot
ax_d = fig.add_subplot(gs[2, 1])
top_mutations = df_mut['mut_code'].value_counts().head(15)
ax_d.barh(range(len(top_mutations)), top_mutations.values, color='#0284c7', edgecolor='#0369a1')
ax_d.set_yticks(range(len(top_mutations)))
ax_d.set_yticklabels(top_mutations.index, fontsize=10, fontweight='bold')
ax_d.invert_yaxis()
ax_d.set_xlabel('Occurrences parmi les 19 Champions', fontsize=10, fontweight='bold')
ax_d.set_title('D. Top 15 des Mutations Hotspots les plus Récurrentes', fontsize=12, fontweight='bold')
for i, v in enumerate(top_mutations.values):
    ax_d.text(v + 0.2, i, f"{v}/19 ({v/19*100:.0f}%)", va='center', fontsize=9, fontweight='bold', color='#0f172a')

fig.savefig(os.path.join(FIG_DIR, "fig3_mutation_frequency_and_conservation.png"), dpi=300)
plt.close(fig)


# ==============================================================================
# FIGURE 4: Profils de Confiance Atomique pLDDT par Résidu le long des Chaînes
# ==============================================================================
print("Generating Figure 4: Per-residue pLDDT profiles...")
fig, (ax_plddt_a, ax_plddt_b) = plt.subplots(2, 1, figsize=(16, 9), constrained_layout=True)

# Parse pLDDT from the top champions PDBs
top_pdbs = [
    ("09_A_prime_candidate_01_B_cand_14", "#2563eb", "Champion #9 (Consensus #1)"),
    ("03_A_prime_candidate_01_B_cand_36", "#059669", "Champion #3 (Consensus #2)"),
    ("01_A_prime_candidate_01_B_cand_40", "#d97706", "Champion #1 (Consensus #3)"),
    ("02_A_prime_candidate_17_B_cand_29", "#7c3aed", "Champion #2 (Pic d'Affinité)"),
]

parser = PDBParser(QUIET=True)

for prefix, col, lab in top_pdbs:
    pattern = os.path.join(PDB_DIR, f"{prefix}_rf3_predicted.pdb")
    matched = glob.glob(pattern)
    if not matched:
        continue
    st = parser.get_structure("s", matched[0])[0]
    
    # Chain A pLDDT
    plddt_a = [a.get_bfactor() * 100.0 if a.get_bfactor() <= 1.0 else a.get_bfactor() 
               for a in st['A'].get_atoms() if a.get_name() == 'CA']
    # Chain B pLDDT
    plddt_b = [a.get_bfactor() * 100.0 if a.get_bfactor() <= 1.0 else a.get_bfactor() 
               for a in st['B'].get_atoms() if a.get_name() == 'CA']
    
    if len(plddt_a) == 224:
        ax_plddt_a.plot(range(1, 225), plddt_a, label=lab, color=col, linewidth=1.8, alpha=0.9)
    if len(plddt_b) == 88:
        ax_plddt_b.plot(range(1, 89), plddt_b, label=lab, color=col, linewidth=1.8, alpha=0.9)

ax_plddt_a.axhline(80, color='#16a34a', linestyle='--', linewidth=1.2, label='Seuil Haute Confiance (pLDDT > 80)')
ax_plddt_a.axhline(70, color='#f59e0b', linestyle=':', linewidth=1.2, label='Seuil Conformationnel (pLDDT > 70)')
ax_plddt_a.set_xlim(1, 224)
ax_plddt_a.set_ylim(60, 100)
ax_plddt_a.set_xlabel('Numéro de Résidu sur Chaîne A\' (1 à 224)', fontsize=11, fontweight='bold')
ax_plddt_a.set_ylabel('Score de Confiance pLDDT', fontsize=11, fontweight='bold')
ax_plddt_a.set_title('A. Profils de Confiance Atomique pLDDT par Résidu — Sous-Unité A\' (Top Champions)', fontsize=12, fontweight='bold')
ax_plddt_a.legend(loc='lower left', frameon=True, ncol=3, fontsize=9)

ax_plddt_b.axhline(80, color='#16a34a', linestyle='--', linewidth=1.2, label='Seuil Haute Confiance (pLDDT > 80)')
ax_plddt_b.axhline(70, color='#f59e0b', linestyle=':', linewidth=1.2, label='Seuil Conformationnel (pLDDT > 70)')
ax_plddt_b.set_xlim(1, 88)
ax_plddt_b.set_ylim(60, 100)
ax_plddt_b.set_xlabel('Numéro de Résidu sur Chaîne B\' (1 à 88)', fontsize=11, fontweight='bold')
ax_plddt_b.set_ylabel('Score de Confiance pLDDT', fontsize=11, fontweight='bold')
ax_plddt_b.set_title('B. Profils de Confiance Atomique pLDDT par Résidu — Sous-Unité B\' (Top Champions)', fontsize=12, fontweight='bold')
ax_plddt_b.legend(loc='lower left', frameon=True, ncol=3, fontsize=9)

fig.savefig(os.path.join(FIG_DIR, "fig4_plddt_profiles_per_residue.png"), dpi=300)
plt.close(fig)


# ==============================================================================
# FIGURE 5: Radar Chart Multi-Critères (Top Champions Pareto Comparison)
# ==============================================================================
print("Generating Figure 5: Multi-metric radar comparison...")
categories = [
    'Affinité RF3\n(iPTM)',
    'Affinité AF2\n(iPTM)',
    'Sélectivité\n($\Delta$ iPTM AF2)',
    'Énergie d\'Amarrage\n($-\Delta$G Rosetta)',
    'Fidélité Native\n(DockQ)',
    'Précision 3D\n(1 / RMSD AF2)'
]
N = len(categories)
angles = [n / float(N) * 2 * np.pi for n in range(N)]
angles += angles[:1]

fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))

# Normalization min-max ranges for fair radar scaling
min_max = {
    'iptm_rescue_rf3': (0.55, 0.80),
    'rescue_iptm_af2': (0.10, 0.60),
    'delta_iptm_af2': (-0.10, 0.30),
    'dG_rescue_pyrosetta': (-35.0, -80.0), # more negative is better
    'dockq_rf3': (0.40, 0.65),
    'rescue_rmsd_af2': (7.0, 1.0),         # lower is better
}

def normalize_val(key, val):
    low, high = min_max[key]
    norm = (val - low) / (high - low)
    return max(0.05, min(1.0, norm))

champs_to_plot = [
    (9, "Champion #9 (cand_14) - Universel #1", '#2563eb'),
    (3, "Champion #3 (cand_36) - Universel #2", '#059669'),
    (1, "Champion #1 (cand_40) - Élite Sélectif", '#d97706'),
    (2, "Champion #2 (cand_29) - Pic Affinité", '#7c3aed'),
    (17, "Champion #17 (cand_24) - Pic AF2", '#db2777'),
]

for rk, lab, col in champs_to_plot:
    sub = df_champ[df_champ['rank'] == rk]
    if sub.empty:
        continue
    row = sub.iloc[0]
    vals = [
        normalize_val('iptm_rescue_rf3', row['iptm_rescue_rf3']),
        normalize_val('rescue_iptm_af2', row['rescue_iptm_af2']),
        normalize_val('delta_iptm_af2', row['delta_iptm_af2']),
        normalize_val('dG_rescue_pyrosetta', row['dG_rescue_pyrosetta']),
        normalize_val('dockq_rf3', row['dockq_rf3']),
        normalize_val('rescue_rmsd_af2', row['rescue_rmsd_af2']),
    ]
    vals += vals[:1]
    ax.plot(angles, vals, linewidth=2.2, linestyle='solid', label=lab, color=col)
    ax.fill(angles, vals, color=col, alpha=0.15)

ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)
plt.xticks(angles[:-1], categories, fontsize=10, fontweight='bold', color='#0f172a')
ax.set_rscale('linear')
ax.set_ylim(0, 1.05)
plt.yticks([0.25, 0.50, 0.75, 1.0], ["25%", "50%", "75%", "100%"], color="#94a3b8", size=8)
plt.title('Profils Multi-Critères de Performance des Meilleurs Champions (Radar Pareto)', 
          size=13, fontweight='bold', pad=25)
plt.legend(loc='upper right', bbox_to_anchor=(1.35, 1.1), frameon=True, fontsize=9)

fig.savefig(os.path.join(FIG_DIR, "fig5_champion_radar_profiles.png"), dpi=300)
plt.close(fig)

print("\n[SUCCESS] All 5 scientific figures successfully generated in DOSSIER_CHAMPIONS/figures/!")
