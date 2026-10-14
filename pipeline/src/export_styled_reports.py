"""Export Styled Reports: Scientific HTML Dashboard with Publication-Grade Vector Figures and Enriched Excel Report.

Features:
- Headless HPC / WSL / SLURM matplotlib Agg backend configuration.
- Strict alignment on Module 5 'passes' verdict as single source of truth.
- 6 Scientific figures (Nature / Science / Bioinformatics standards) exported inline (SVG)
  and saved to <out_dir>/figures/ as .svg and .png (300 DPI).
- French UTF-8 encoding restored, dead code removed.
- Interactive dashboard with key takeaways, Top 3 podium cards, parent motif summary,
  fast JS filtering and column sorting, plus detailed per-row rejection diagnostics.
"""

import io
import math
import os
import sys

import numpy as np
import openpyxl
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
import pandas as pd

# ---------------------------------------------------------------------------
# HPC Headless Matplotlib configuration (Guarantee no $DISPLAY / X11 errors)
# ---------------------------------------------------------------------------
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def _clean_svg_for_html(svg_str):
    """Strip xml declarations / doctypes so SVG can be cleanly embedded inline in HTML."""
    svg_start = svg_str.find("<svg")
    if svg_start != -1:
        return svg_str[svg_start:]
    return svg_str


def _save_and_close_fig(fig, fig_dir, base_name):
    """Save figure to SVG and 300 DPI PNG if fig_dir is provided, return cleaned SVG string and close figure."""
    sio = io.StringIO()
    fig.savefig(sio, format='svg', bbox_inches='tight')
    svg_data = sio.getvalue()

    if fig_dir:
        os.makedirs(fig_dir, exist_ok=True)
        svg_file = os.path.join(fig_dir, f"{base_name}.svg")
        png_file = os.path.join(fig_dir, f"{base_name}.png")
        try:
            fig.savefig(svg_file, format='svg', bbox_inches='tight')
            fig.savefig(png_file, format='png', dpi=300, bbox_inches='tight')
        except Exception as e:
            print(f"Warning: Could not save figure {base_name}: {e}")

    plt.close(fig)
    return _clean_svg_for_html(svg_data)


def _empty_figure(title, message, fig_dir, base_name):
    """Fallback placeholder when data is insufficient or empty."""
    fig, ax = plt.subplots(figsize=(8, 4.5), facecolor='#ffffff')
    ax.text(0.5, 0.5, message, ha='center', va='center', fontsize=12, color='#64748b', fontweight='medium')
    ax.set_title(title, fontsize=13, fontweight='bold', color='#1e293b', pad=14)
    ax.axis('off')
    return _save_and_close_fig(fig, fig_dir, base_name)


def _safe_col(df, col, default=0.0):
    """Safely return a column as Series with fallback default if missing."""
    if df is not None and col in df.columns:
        return df[col].fillna(default)
    if df is not None:
        return pd.Series(default, index=df.index)
    return pd.Series([default])


# ===========================================================================
# FIGURE 1: Matrices de Réactivité Croisée 2x2 du Meilleur Candidat
# ===========================================================================
def plot_cross_reactivity_matrix(df, fig_dir=None):
    """Figure 1: Orthogonality Cross-Reactivity 2x2 Heatmaps (Boyken et al., Science 2016 style).

    Panel 1A: RF3 iPTM confidence matrix (WT ceiling, Negative, Rupture, Rescue).
    Panel 1B: PyRosetta interface dG in REU (or ipSAE Dunbrack 2025 if PyRosetta is off).
    """
    if df is None or df.empty:
        return _empty_figure("Figure 1 : Matrices de Réactivité Croisée",
                            "Aucune donnée de conception disponible", fig_dir, "figure1_cross_reactivity")

    # Select champion design
    sort_cols = [c for c in ['passes', 'f_ortho', 'dockq'] if c in df.columns]
    if 'passes' in df.columns and df['passes'].any():
        champ = df[df['passes'].fillna(False).astype(bool)].sort_values(
            by=sort_cols, ascending=[False]*len(sort_cols)
        ).iloc[0]
    elif sort_cols:
        champ = df.sort_values(by=sort_cols, ascending=[False]*len(sort_cols)).iloc[0]
    else:
        champ = df.iloc[0]

    did = str(champ.get('design_id', 'Champion'))

    # Subfigure 1A : RF3 iPTM
    iptm_ceil = float(champ.get('iptm_ceiling', 0.85)) if pd.notnull(champ.get('iptm_ceiling')) else 0.85
    iptm_neg = float(champ.get('iptm_negative', 0.0)) if pd.notnull(champ.get('iptm_negative')) else 0.0
    iptm_rup = float(champ.get('iptm_rupture', 0.0)) if pd.notnull(champ.get('iptm_rupture')) else 0.0
    iptm_resc = float(champ.get('iptm_rescue', 0.0)) if pd.notnull(champ.get('iptm_rescue')) else 0.0

    mat_iptm = np.array([
        [iptm_ceil, iptm_neg],
        [iptm_rup,  iptm_resc]
    ])

    # Check for PyRosetta dG or ipSAE in Subfigure 1B
    has_energy = any(pd.notnull(champ.get(k)) for k in ['dG_native', 'dG_rescue', 'dG_rupture', 'dG_negative'])
    has_ipsae = any(pd.notnull(champ.get(k)) for k in ['ipsae_ceiling', 'ipsae_rescue', 'ipsae_rupture', 'ipsae_negative'])

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.8), facecolor='#ffffff')
    plt.subplots_adjust(wspace=0.35, top=0.76, bottom=0.14)

    # --- 1A: RF3 iPTM ---
    ax_a = axes[0]
    im_a = ax_a.imshow(mat_iptm, cmap='YlGnBu', vmin=0.0, vmax=1.0, aspect='auto')
    ax_a.set_xticks([0, 1])
    ax_a.set_yticks([0, 1])
    ax_a.set_xticklabels([r'$B_{WT}$', r"$B'$"], fontsize=12, fontweight='bold', color='#1e293b')
    ax_a.set_yticklabels([r'$A_{WT}$', r"$A'$"], fontsize=12, fontweight='bold', color='#1e293b')
    ax_a.set_title("1A : Confiance RF3 iPTM\n(Boyken et al. 2016)", fontsize=11.5, fontweight='bold', color='#0f172a', pad=8)

    roles_a = [
        ["Plafond Sauvage\n(WT · WT)", "Hors-Cible Négatif\n(A_WT · B')"],
        ["Rupture Spécifique\n(A' · B_WT)", "Sauvetage Orthogonal\n(A' · B')"]
    ]

    for i in range(2):
        for j in range(2):
            val = mat_iptm[i, j]
            txt_color = '#ffffff' if val > 0.55 else '#0f172a'
            ax_a.text(j, i - 0.10, f"{val:.3f}", ha='center', va='center', fontsize=14, fontweight='bold', color=txt_color)
            ax_a.text(j, i + 0.18, roles_a[i][j], ha='center', va='center', fontsize=8.5, color=txt_color, style='italic')

    cbar_a = fig.colorbar(im_a, ax=ax_a, fraction=0.046, pad=0.04)
    cbar_a.set_label('Score iPTM (0-1)', fontsize=10, color='#334155')

    # --- 1B: PyRosetta dG or ipSAE ---
    ax_b = axes[1]
    if has_energy:
        dg_nat = float(champ.get('dG_native', np.nan))
        dg_neg = float(champ.get('dG_negative', np.nan))
        dg_rup = float(champ.get('dG_rupture', np.nan))
        dg_resc = float(champ.get('dG_rescue', np.nan))

        mat_energy = np.array([
            [dg_nat if pd.notnull(dg_nat) else 0.0, dg_neg if pd.notnull(dg_neg) else 0.0],
            [dg_rup if pd.notnull(dg_rup) else 0.0, dg_resc if pd.notnull(dg_resc) else 0.0]
        ])
        vabs = max(15.0, np.nanmax(np.abs(mat_energy)))
        im_b = ax_b.imshow(mat_energy, cmap='RdBu_r', vmin=-vabs, vmax=vabs, aspect='auto')
        ax_b.set_title("1B : Énergie Libre d'Interface PyRosetta\n(ΔG en REU - Bleu = Liaison Forte)",
                       fontsize=11.5, fontweight='bold', color='#0f172a', pad=8)

        for i in range(2):
            for j in range(2):
                val = mat_energy[i, j]
                # Dark text on light backgrounds near 0, white text on strong blue/red (|val| >= 35 REU)
                txt_color = '#ffffff' if (val <= -35.0 or val >= 35.0) else '#0f172a'
                ax_b.text(j, i - 0.10, f"{val:.1f} REU", ha='center', va='center', fontsize=13, fontweight='bold', color=txt_color)
                ax_b.text(j, i + 0.18, roles_a[i][j], ha='center', va='center', fontsize=8.5, color=txt_color, style='italic')

        cbar_b = fig.colorbar(im_b, ax=ax_b, fraction=0.046, pad=0.04)
        cbar_b.set_label('ΔG Interface (REU)', fontsize=10, color='#334155')

    elif has_ipsae:
        ips_ceil = float(champ.get('ipsae_ceiling', 0.8)) if pd.notnull(champ.get('ipsae_ceiling')) else 0.8
        ips_neg = float(champ.get('ipsae_negative', 0.0)) if pd.notnull(champ.get('ipsae_negative')) else 0.0
        ips_rup = float(champ.get('ipsae_rupture', 0.0)) if pd.notnull(champ.get('ipsae_rupture')) else 0.0
        ips_resc = float(champ.get('ipsae_rescue', 0.0)) if pd.notnull(champ.get('ipsae_rescue')) else 0.0

        mat_ips = np.array([
            [ips_ceil, ips_neg],
            [ips_rup,  ips_resc]
        ])
        im_b = ax_b.imshow(mat_ips, cmap='Purples', vmin=0.0, vmax=1.0, aspect='auto')
        ax_b.set_title("1B : Métrique ipSAE (Dunbrack 2025)\n(Résolution Locale PAE)",
                       fontsize=11.5, fontweight='bold', color='#0f172a', pad=8)

        for i in range(2):
            for j in range(2):
                val = mat_ips[i, j]
                txt_color = '#ffffff' if val > 0.55 else '#0f172a'
                ax_b.text(j, i - 0.10, f"{val:.3f}", ha='center', va='center', fontsize=14, fontweight='bold', color=txt_color)
                ax_b.text(j, i + 0.18, roles_a[i][j], ha='center', va='center', fontsize=8.5, color=txt_color, style='italic')

        cbar_b = fig.colorbar(im_b, ax=ax_b, fraction=0.046, pad=0.04)
        cbar_b.set_label('Score ipSAE (0-1)', fontsize=10, color='#334155')
    else:
        ax_b.text(0.5, 0.5, "PyRosetta & ipSAE désactivés\nSeul RF3 iPTM est disponible",
                  ha='center', va='center', fontsize=11, color='#64748b')
        ax_b.set_title("1B : Physique Explicite (Non Exécutée)", fontsize=11.5, fontweight='bold', color='#0f172a', pad=8)

    ax_b.set_xticks([0, 1])
    ax_b.set_yticks([0, 1])
    ax_b.set_xticklabels([r'$B_{WT}$', r"$B'$"], fontsize=12, fontweight='bold', color='#1e293b')
    ax_b.set_yticklabels([r'$A_{WT}$', r"$A'$"], fontsize=12, fontweight='bold', color='#1e293b')

    fig.suptitle(f"Figure 1 : Matrices d'Orthogonalité & Réactivité Croisée — Candidat Champion : {did}\n"
                 "Démontre la diagonalisation d'interaction (forte affinité WT:WT & A':B', extinction hors-cible)",
                 fontsize=12.5, fontweight='bold', color='#0f172a', y=0.96)

    return _save_and_close_fig(fig, fig_dir, "figure1_cross_reactivity")


# ===========================================================================
# FIGURE 2: Paysage de Spécificité & Isolignes de Marge
# ===========================================================================
def plot_specificity_landscape(df, fig_dir=None):
    """Figure 2: Specificity Landscape: On-Target vs Worst Off-Target with margin isolines.

    X-axis: Worst wild-type off-target cross-reaction max(iPTM_rupture, iPTM_negative).
    Y-axis: Orthogonal rescue affinity iPTM_rescue.
    """
    if df is None or df.empty:
        return _empty_figure("Figure 2 : Paysage de Spécificité & Isolignes de Marge",
                            "Aucune donnée de conception disponible", fig_dir, "figure2_specificity_landscape")

    x_vals = np.maximum(_safe_col(df, 'iptm_rupture', 0.0), _safe_col(df, 'iptm_negative', 0.0))
    y_vals = _safe_col(df, 'iptm_rescue', 0.0)

    fig, ax = plt.subplots(figsize=(9.6, 7.2), facecolor='#ffffff')
    plt.subplots_adjust(top=0.88, bottom=0.12)
    ax.set_facecolor('#f8fafc')
    ax.grid(True, linestyle='--', color='#cbd5e1', alpha=0.6, zorder=1)

    # Neutrality diagonal y = x
    line_range = np.linspace(0.0, 1.0, 100)
    ax.plot(line_range, line_range, color='#dc2626', linestyle='--', linewidth=1.8,
            label='Diagonale de Neutralité : y = x (F_ortho,raw = 0)', zorder=2)
    ax.fill_between(line_range, 0, line_range, color='#fee2e2', alpha=0.45,
                    label='Zone Défavorable : Préférence pour le WT (y < x)', zorder=1)

    # Margin isolines y = x + delta
    for delta, color, ls in [(0.10, '#94a3b8', ':'), (0.20, '#64748b', ':'), (0.30, '#059669', '-.')]:
        y_iso = line_range + delta
        mask = y_iso <= 1.0
        ax.plot(line_range[mask], y_iso[mask], color=color, linestyle=ls, linewidth=1.4,
                label=f'Isoligne Marge F_ortho = +{delta:.2f}', zorder=2)

    # Threshold horizontal / vertical reference lines
    ax.axhline(0.65, color='#10b981', linestyle=':', linewidth=1.3, alpha=0.75, label='Seuil Sauvetage min (0.65)')
    ax.axvline(0.35, color='#f59e0b', linestyle=':', linewidth=1.3, alpha=0.75, label='Seuil Rupture max (0.35)')

    # Ceiling if available
    ceilings = df['iptm_ceiling'].dropna() if 'iptm_ceiling' in df.columns else pd.Series(dtype=float)
    if len(ceilings) > 0 and ceilings.median() > 0:
        med_ceil = ceilings.median()
        ax.axhline(med_ceil, color='#7c3aed', linestyle='--', linewidth=1.3, alpha=0.7,
                   label=f'Plafond Sauvage Médian ({med_ceil:.2f})')

    # Categorize designs
    pass_col = df['passes'].fillna(False).astype(bool) if 'passes' in df.columns else (
        (_safe_col(df, 'f_ortho', 0.0) > 0) & (_safe_col(df, 'dockq', 0.0) >= 0.23)
    )
    f_ortho_col = _safe_col(df, 'f_ortho', 0.0)
    dockq_col = _safe_col(df, 'dockq', 0.0)
    status_col = df['status'].fillna('ok').astype(str) if 'status' in df.columns else pd.Series(['ok']*len(df), index=df.index)

    is_champ = pass_col & (f_ortho_col >= 0.30)
    is_ortho = pass_col & (~is_champ)
    is_off_target = status_col.str.contains('off_target') | (dockq_col < 0.23)
    is_rejected_affinity = (~pass_col) & (~is_off_target)

    # Scatter points by category
    if is_off_target.any():
        ax.scatter(x_vals[is_off_target], y_vals[is_off_target], marker='x', s=65, color='#ef4444',
                   linewidths=1.8, label=f'Hors-cible CAPRI (DockQ < 0.23, n={is_off_target.sum()})', zorder=4)

    if is_rejected_affinity.any():
        ax.scatter(x_vals[is_rejected_affinity], y_vals[is_rejected_affinity], marker='s', s=70,
                   facecolors='#fed7aa', edgecolors='#d97706', linewidths=1.3,
                   label=f'Rejeté Spécificité / Énergie (n={is_rejected_affinity.sum()})', zorder=5)

    if is_ortho.any():
        ax.scatter(x_vals[is_ortho], y_vals[is_ortho], marker='o', s=85,
                   facecolors='#bae6fd', edgecolors='#0284c7', linewidths=1.5,
                   label=f'Orthogonal Validé (n={is_ortho.sum()})', zorder=6)

    if is_champ.any():
        ax.scatter(x_vals[is_champ], y_vals[is_champ], marker='*', s=240,
                   facecolors='#34d399', edgecolors='#065f46', linewidths=1.6,
                   label=f'Champion Orthogonal (F_ortho ≥ 0.30, n={is_champ.sum()})', zorder=7)

    # Annotate Top 3 designs with non-overlapping offsets
    sorted_df = df.copy()
    sorted_df['_x'] = x_vals
    sorted_df['_y'] = y_vals
    sorted_df['_pass'] = pass_col
    sort_keys = ['_pass']
    if 'f_ortho' in sorted_df.columns:
        sort_keys.append('f_ortho')
    if 'dockq' in sorted_df.columns:
        sort_keys.append('dockq')

    top_candidates = sorted_df.sort_values(by=sort_keys, ascending=[False]*len(sort_keys)).head(3)
    # Directional divergence: #1 pointing up-left, #2 pointing down-right, #3 pointing up-right
    offsets = [(-150, 30), (30, -32), (35, 25)]

    for idx_pt, (_, row) in enumerate(top_candidates.iterrows()):
        rank = idx_pt + 1
        did = str(row.get('design_id', f'Design_{rank}'))
        tx = row['_x']
        ty = row['_y']
        offset = offsets[idx_pt] if idx_pt < len(offsets) else (18, 16)
        ax.annotate(
            f"#{rank} {did}\n(F_ortho: {row.get('f_ortho', 0.0):+.3f})",
            xy=(tx, ty), xytext=offset, textcoords="offset points",
            fontsize=8.5, fontweight='bold', color='#0f172a',
            bbox=dict(boxstyle="round,pad=0.25", fc="#ffffff", ec="#0284c7", lw=1.2, alpha=0.95),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.1", color="#0284c7", lw=1.2),
            zorder=8
        )

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlabel("Pire Réaction Croisée Sauvage : max(iPTM rupture, iPTM négatif)\n"
                  "← Meilleure Sélectivité Contre le WT (Extinction Désirée)",
                  fontsize=11, fontweight='bold', color='#1e293b', labelpad=8)
    ax.set_ylabel("Affinité de Sauvetage : iPTM rescue (A' · B')\n"
                  "↑ Meilleure Affinité Orthogonale Réelle",
                  fontsize=11, fontweight='bold', color='#1e293b', labelpad=8)

    ax.set_title("Figure 2 : Paysage de Spécificité & Isolignes de Marge d'Orthogonalité\n"
                 "Isolement de la fenêtre orthogonale supérieure face aux réactions croisées sauvages",
                 fontsize=12, fontweight='bold', color='#0f172a', pad=14)

    ax.legend(loc='lower right', fontsize=8.5, framealpha=0.95, facecolor='#ffffff', edgecolor='#cbd5e1')

    return _save_and_close_fig(fig, fig_dir, "figure2_specificity_landscape")


# ===========================================================================
# FIGURE 3: Distributions de la Triade d'Interaction (Violin / Box + Jitter)
# ===========================================================================
def plot_interaction_triad_distributions(df, fig_dir=None):
    """Figure 3: Interaction triad distributions (Violin / Box + Swarm Strip Plot + Trajectories).

    Panels: 3A (RF3 iPTM distributions) and 3B (ipSAE Dunbrack 2025 distributions).
    """
    if df is None or df.empty:
        return _empty_figure("Figure 3 : Distributions de la Triade d'Interaction",
                            "Aucune donnée de conception disponible", fig_dir, "figure3_triad_distributions")

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.8), facecolor='#ffffff')
    plt.subplots_adjust(wspace=0.28, top=0.75, bottom=0.16)

    # --- 3A: RF3 iPTM ---
    ax_a = axes[0]
    ax_a.set_facecolor('#f8fafc')
    ax_a.grid(True, linestyle='--', color='#e2e8f0', alpha=0.7, zorder=1)

    labels_a = [r'$WT \cdot WT$' + '\n(Plafond)',
                r"$A' \cdot B_{WT}$" + '\n(Rupture)',
                r"$A_{WT} \cdot B'$" + '\n(Négatif)',
                r"$A' \cdot B'$" + '\n(Sauvetage)']

    data_a = [
        df['iptm_ceiling'].dropna().values if 'iptm_ceiling' in df.columns else np.array([]),
        df['iptm_rupture'].dropna().values if 'iptm_rupture' in df.columns else np.array([]),
        df['iptm_negative'].dropna().values if 'iptm_negative' in df.columns else np.array([]),
        df['iptm_rescue'].dropna().values if 'iptm_rescue' in df.columns else np.array([])
    ]

    colors_a = ['#8b5cf6', '#ef4444', '#f97316', '#10b981']

    valid_data_a = [d if len(d) > 0 else np.array([0.0]) for d in data_a]
    bp_a = ax_a.boxplot(valid_data_a, positions=[1, 2, 3, 4], patch_artist=True, widths=0.45, zorder=2,
                        medianprops=dict(color='#0f172a', linewidth=2.0),
                        boxprops=dict(linewidth=1.2, edgecolor='#475569'))

    for patch, col in zip(bp_a['boxes'], colors_a):
        patch.set_facecolor(col)
        patch.set_alpha(0.25)

    np.random.seed(42)
    for idx, (vals, col) in enumerate(zip(data_a, colors_a), start=1):
        if len(vals) > 0:
            jitter = np.random.normal(0, 0.04, size=len(vals))
            ax_a.scatter(idx + jitter, vals, color=col, edgecolors='#0f172a', s=35, alpha=0.75, zorder=3)

    pass_mask = df['passes'].fillna(False).astype(bool) if 'passes' in df.columns else (
        (_safe_col(df, 'f_ortho', 0.0) > 0) & (_safe_col(df, 'dockq', 0.0) >= 0.23)
    )
    passed_df = df[pass_mask]

    for _, row in passed_df.iterrows():
        y_pts = [row.get('iptm_rupture', np.nan), row.get('iptm_negative', np.nan), row.get('iptm_rescue', np.nan)]
        if not any(pd.isna(y_pts)):
            ax_a.plot([2, 3, 4], y_pts, color='#10b981', linewidth=1.5, alpha=0.7, zorder=4)

    ax_a.axhline(0.35, color='#ef4444', linestyle=':', linewidth=1.3, label='Seuil Rupture max (0.35)')
    ax_a.axhline(0.65, color='#10b981', linestyle=':', linewidth=1.3, label='Seuil Sauvetage min (0.65)')

    ax_a.set_xticks([1, 2, 3, 4])
    ax_a.set_xticklabels(labels_a, fontsize=10, fontweight='bold', color='#1e293b')
    ax_a.set_ylim(-0.05, 1.05)
    ax_a.set_ylabel('Score iPTM (0-1)', fontsize=11, fontweight='bold', color='#1e293b')
    ax_a.set_title("3A : Triade d'Interaction RF3\n(Lignes vertes = Candidats Validés)",
                   fontsize=11.5, fontweight='bold', color='#0f172a', pad=8)
    ax_a.legend(loc='lower left', fontsize=8.0, framealpha=0.9)

    # --- 3B: ipSAE Score (Dunbrack 2025) ---
    ax_b = axes[1]
    ax_b.set_facecolor('#f8fafc')
    ax_b.grid(True, linestyle='--', color='#e2e8f0', alpha=0.7, zorder=1)

    has_ipsae = any(('ipsae_' + k in df.columns and df['ipsae_' + k].dropna().any())
                    for k in ['ceiling', 'rupture', 'negative', 'rescue'])

    if has_ipsae:
        labels_b = [r'$WT \cdot WT$' + '\n(ipSAE Ceil)',
                    r"$A' \cdot B_{WT}$" + '\n(ipSAE Rup)',
                    r"$A_{WT} \cdot B'$" + '\n(ipSAE Neg)',
                    r"$A' \cdot B'$" + '\n(ipSAE Resc)']

        data_b = [
            df['ipsae_ceiling'].dropna().values if 'ipsae_ceiling' in df.columns else np.array([]),
            df['ipsae_rupture'].dropna().values if 'ipsae_rupture' in df.columns else np.array([]),
            df['ipsae_negative'].dropna().values if 'ipsae_negative' in df.columns else np.array([]),
            df['ipsae_rescue'].dropna().values if 'ipsae_rescue' in df.columns else np.array([])
        ]

        valid_data_b = [d if len(d) > 0 else np.array([0.0]) for d in data_b]
        bp_b = ax_b.boxplot(valid_data_b, positions=[1, 2, 3, 4], patch_artist=True, widths=0.45, zorder=2,
                            medianprops=dict(color='#0f172a', linewidth=2.0),
                            boxprops=dict(linewidth=1.2, edgecolor='#475569'))

        for patch, col in zip(bp_b['boxes'], colors_a):
            patch.set_facecolor(col)
            patch.set_alpha(0.25)

        for idx, (vals, col) in enumerate(zip(data_b, colors_a), start=1):
            if len(vals) > 0:
                jitter = np.random.normal(0, 0.04, size=len(vals))
                ax_b.scatter(idx + jitter, vals, color=col, edgecolors='#0f172a', s=35, alpha=0.75, zorder=3)

        for _, row in passed_df.iterrows():
            y_pts = [row.get('ipsae_rupture', np.nan), row.get('ipsae_negative', np.nan), row.get('ipsae_rescue', np.nan)]
            if not any(pd.isna(y_pts)):
                ax_b.plot([2, 3, 4], y_pts, color='#10b981', linewidth=1.5, alpha=0.7, zorder=4)

        ax_b.set_xticks([1, 2, 3, 4])
        ax_b.set_xticklabels(labels_b, fontsize=10, fontweight='bold', color='#1e293b')
        ax_b.set_ylim(-0.05, 1.05)
        ax_b.set_ylabel('Score ipSAE (Dunbrack 2025)', fontsize=11, fontweight='bold', color='#1e293b')
        ax_b.set_title("3B : ipSAE Dunbrack 2025\n(Validation Résiduelle Hors-Biais Alignement)",
                       fontsize=11.5, fontweight='bold', color='#0f172a', pad=8)
    else:
        ax_b.text(0.5, 0.5, "Données ipSAE non calculées sur ce run\n(Option activable via config.yaml)",
                  ha='center', va='center', fontsize=11, color='#64748b', fontweight='medium')
        ax_b.set_title("3B : ipSAE Dunbrack 2025 (Non Calculé)", fontsize=11.5, fontweight='bold', color='#0f172a', pad=8)
        ax_b.set_xticks([])
        ax_b.set_yticks([])

    fig.suptitle("Figure 3 : Distributions de la Triade d'Interaction & Séparation des États\n"
                 "Mise en évidence de la vallée d'extinction cross-reactive chez les candidats validés (lignes vertes)",
                 fontsize=12, fontweight='bold', color='#0f172a', y=0.96)

    return _save_and_close_fig(fig, fig_dir, "figure3_triad_distributions")


# ===========================================================================
# FIGURE 4: Entonnoir Structural CAPRI & Auto-Cohérence (Funnel + Spearman)
# ===========================================================================
def plot_structural_funnel(df, fig_dir=None, coherence_path=None):
    """Figure 4: Structural Folding Funnel: iPTM vs L-RMSD & DockQ, and IA vs Physics concordance.

    Panel 4A: selfcons_lrms (or lrms) vs iptm_rescue colored by dockq_quality.
    Panel 4B: Concordance Structural AI (F_iptm_rel) vs Explicit Physics (F_energy) with Spearman rho.
    """
    if df is None or df.empty:
        return _empty_figure("Figure 4 : Entonnoir Structural CAPRI & Auto-Cohérence",
                            "Aucune donnée de conception disponible", fig_dir, "figure4_structural_funnel")

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.8), facecolor='#ffffff')
    plt.subplots_adjust(wspace=0.30, top=0.75, bottom=0.16)

    # --- 4A: Folding Funnel CAPRI ---
    ax_a = axes[0]
    ax_a.set_facecolor('#f8fafc')
    ax_a.grid(True, linestyle='--', color='#e2e8f0', alpha=0.7, zorder=1)

    if 'selfcons_lrms' in df.columns and df['selfcons_lrms'].dropna().any():
        x_col = 'selfcons_lrms'
        x_label = "Auto-Cohérence L-RMSD vs Modèle Cible (Å)"
    elif 'lrms' in df.columns and df['lrms'].dropna().any():
        x_col = 'lrms'
        x_label = "L-RMSD CAPRI vs Complexe Natif (Å)"
    elif 'rmsd_b_prime' in df.columns and df['rmsd_b_prime'].dropna().any():
        x_col = 'rmsd_b_prime'
        x_label = "RMSD de Chaîne B' vs WT (Å)"
    else:
        x_col = 'dockq'
        x_label = "Score DockQ"

    x_funnel = _safe_col(df, x_col, 5.0)
    y_funnel = _safe_col(df, 'iptm_rescue', 0.0)

    quality_map = {
        'High': ('#10b981', 'o', 100, 'CAPRI High (DockQ ≥ 0.80)'),
        'Medium': ('#0284c7', 'o', 85, 'CAPRI Medium (0.49 ≤ DockQ < 0.80)'),
        'Acceptable': ('#f59e0b', 'o', 70, 'CAPRI Acceptable (0.23 ≤ DockQ < 0.49)'),
        'Incorrect': ('#ef4444', 'x', 60, 'CAPRI Incorrect (DockQ < 0.23)')
    }

    max_x = max(12.0, x_funnel.max() + 1.0)
    ax_a.axvspan(0.0, 5.0, color='#dcfce7', alpha=0.45, label='Puits Natif CAPRI (L-RMSD < 5 Å)', zorder=1)
    ax_a.axvline(5.0, color='#10b981', linestyle='--', linewidth=1.4, zorder=2)
    ax_a.axhline(0.65, color='#64748b', linestyle=':', linewidth=1.3, zorder=2, label='Seuil Sauvetage (0.65)')

    dq_qualities = df.get('dockq_quality', pd.Series(['Acceptable'] * len(df), index=df.index)).fillna('Incorrect')
    for q_label, (col, marker, size, legend_lbl) in quality_map.items():
        mask = (dq_qualities == q_label)
        if mask.any():
            if marker == 'x':
                ax_a.scatter(x_funnel[mask], y_funnel[mask], color=col, marker=marker, s=size,
                             linewidths=1.2, label=legend_lbl, zorder=4)
            else:
                ax_a.scatter(x_funnel[mask], y_funnel[mask], color=col, marker=marker, s=size,
                             edgecolors='#0f172a', linewidths=1.2, label=legend_lbl, zorder=4)

    ax_a.set_xlim(-0.2, max_x)
    ax_a.set_ylim(-0.02, 1.05)
    ax_a.set_xlabel(x_label, fontsize=10.5, fontweight='bold', color='#1e293b')
    ax_a.set_ylabel("Affinité Sauvetage iPTM (A' · B')", fontsize=10.5, fontweight='bold', color='#1e293b')
    ax_a.set_title("4A : Entonnoir Structural CAPRI\n(Convergence géométrique vers le puits natif)",
                   fontsize=11.5, fontweight='bold', color='#0f172a', pad=8)
    ax_a.legend(loc='upper right', fontsize=8.0, framealpha=0.92)

    # --- 4B: AI vs Physics Concordance ---
    ax_b = axes[1]
    ax_b.set_facecolor('#f8fafc')
    ax_b.grid(True, linestyle='--', color='#e2e8f0', alpha=0.7, zorder=1)

    has_f_energy = ('f_energy' in df.columns) and df['f_energy'].dropna().any()

    if has_f_energy:
        x_ai = _safe_col(df, 'f_iptm_rel', 0.0) if 'f_iptm_rel' in df.columns else _safe_col(df, 'f_ortho_iptm', 0.0)
        y_phys = _safe_col(df, 'f_energy', 0.0)

        ax_b.axvline(0.0, color='#94a3b8', linestyle='--', linewidth=1.2, zorder=2)
        ax_b.axhline(0.0, color='#94a3b8', linestyle='--', linewidth=1.2, zorder=2)
        ax_b.fill_between([0.0, 1.2], 0.0, 1.2, color='#d1fae5', alpha=0.5,
                          label='Accord Binaire Favorable (IA > 0 & Physique > 0)', zorder=1)

        pass_mask = df['passes'].fillna(False).astype(bool) if 'passes' in df.columns else (
            (_safe_col(df, 'f_ortho', 0.0) > 0) & (_safe_col(df, 'dockq', 0.0) >= 0.23)
        )
        ax_b.scatter(x_ai[pass_mask], y_phys[pass_mask], color='#10b981', s=80, edgecolors='#064e3b',
                     label=f'Orthogonal Validé (n={pass_mask.sum()})', zorder=4)
        ax_b.scatter(x_ai[~pass_mask], y_phys[~pass_mask], color='#ef4444', s=60, marker='x', linewidths=1.5,
                     label=f'Rejeté (n={(~pass_mask).sum()})', zorder=3)

        spearman_rho = None
        if coherence_path and os.path.exists(coherence_path):
            try:
                import json
                with open(coherence_path, 'r', encoding='utf-8') as f:
                    coh_data = json.load(f)
                    spearman_rho = coh_data.get('spearman_margins')
            except Exception:
                pass

        if spearman_rho is None and len(df) > 2:
            try:
                from scipy.stats import spearmanr
                res = spearmanr(x_ai, y_phys)
                spearman_rho = res.statistic if hasattr(res, 'statistic') else res[0]
            except Exception:
                pass

        rho_str = f"Spearman ρ = {spearman_rho:+.2f}" if spearman_rho is not None else "Spearman ρ: calculé"
        ax_b.text(0.05, 0.92, rho_str, transform=ax_b.transAxes, fontsize=10, fontweight='bold',
                  bbox=dict(boxstyle="round,pad=0.3", fc="#ffffff", ec="#10b981", lw=1.2), zorder=5)

        ax_b.set_xlabel("Score IA Structurale RF3 (F_iptm,rel)", fontsize=10.5, fontweight='bold', color='#1e293b')
        ax_b.set_ylabel("Score Physique Explicite PyRosetta (F_energy)", fontsize=10.5, fontweight='bold', color='#1e293b')
        ax_b.set_title("4B : Concordance IA (RF3) vs Physique (PyRosetta)\n(Validation orthogonale croisée)",
                       fontsize=11.5, fontweight='bold', color='#0f172a', pad=8)
        ax_b.legend(loc='lower right', fontsize=8.0, framealpha=0.92)

    else:
        # Fallback: Fnat vs DockQ
        fnat_vals = _safe_col(df, 'fnat', 0.0)
        dockq_vals = _safe_col(df, 'dockq', 0.0)
        pass_mask = df['passes'].fillna(False).astype(bool) if 'passes' in df.columns else (
            (_safe_col(df, 'f_ortho', 0.0) > 0) & (_safe_col(df, 'dockq', 0.0) >= 0.23)
        )

        ax_b.axhline(0.23, color='#f59e0b', linestyle='--', linewidth=1.3, label='Seuil CAPRI Acceptable (0.23)')
        ax_b.axhline(0.49, color='#0284c7', linestyle=':', linewidth=1.3, label='Seuil CAPRI Medium (0.49)')

        ax_b.scatter(fnat_vals[pass_mask], dockq_vals[pass_mask], color='#10b981', s=80, edgecolors='#064e3b',
                     label='Orthogonal Validé', zorder=4)
        ax_b.scatter(fnat_vals[~pass_mask], dockq_vals[~pass_mask], color='#ef4444', marker='x', s=60,
                     label='Rejeté', zorder=3)

        ax_b.set_xlabel("Fraction de Contacts Natifs Conservés (F_nat)", fontsize=10.5, fontweight='bold', color='#1e293b')
        ax_b.set_ylabel("Score DockQ Global", fontsize=10.5, fontweight='bold', color='#1e293b')
        ax_b.set_title("4B : Préservation des Contacts Natifs vs DockQ\n(PyRosetta désactivé)",
                       fontsize=11.5, fontweight='bold', color='#0f172a', pad=8)
        ax_b.legend(loc='lower right', fontsize=8.0, framealpha=0.92)

    fig.suptitle("Figure 4 : Entonnoir Structural CAPRI & Validation Multi-Physique\n"
                 "Vérification de la convergence géométrique haute résolution et de la cohérence énergétique",
                 fontsize=12, fontweight='bold', color='#0f172a', y=0.96)

    return _save_and_close_fig(fig, fig_dir, "figure4_structural_funnel")


# ===========================================================================
# FIGURE 5: Impact du Remodelage de Squelette RFD3 & Qualité Biophysique
# ===========================================================================
def plot_induced_fit_and_biophysics(df, fig_dir=None):
    """Figure 5: RFD3 Induced-Fit Analysis & Interface Biophysical Quality.

    Panel 5A: Performance vs scaffold flexibility RMSD (scaffold_flex_rmsd & scaffold_is_native).
    Panel 5B: Interface dSASA vs unsat_hb (PyRosetta) or Performance by Parent Motif A'.
    """
    if df is None or df.empty:
        return _empty_figure("Figure 5 : Analyse Biomécanique RFD3 & Profil Biophysique",
                            "Aucune donnée de conception disponible", fig_dir, "figure5_induced_fit_biophysics")

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.8), facecolor='#ffffff')
    plt.subplots_adjust(wspace=0.30, top=0.75, bottom=0.16)

    # --- 5A: RFD3 Induced-Fit Analysis ---
    ax_a = axes[0]
    ax_a.set_facecolor('#f8fafc')
    ax_a.grid(True, linestyle='--', color='#e2e8f0', alpha=0.7, zorder=1)

    has_scaffold = 'scaffold_flex_rmsd' in df.columns and df['scaffold_flex_rmsd'].dropna().any()

    if has_scaffold:
        x_flex = df['scaffold_flex_rmsd'].fillna(0.0)
        y_perf = _safe_col(df, 'f_ortho', 0.0)
        is_native = df['scaffold_is_native'].fillna(False).astype(bool) if 'scaffold_is_native' in df.columns else pd.Series([False]*len(df), index=df.index)

        if is_native.any():
            ax_a.scatter(x_flex[is_native], y_perf[is_native], marker='s', s=80, color='#3b82f6',
                         edgecolors='#1d4ed8', linewidths=1.2, label=f'Squelette WT Rigide (n={is_native.sum()})', zorder=4)

        if (~is_native).any():
            ax_a.scatter(x_flex[~is_native], y_perf[~is_native], marker='o', s=85, color='#10b981',
                         edgecolors='#064e3b', linewidths=1.2, label=f'Squelette Remodelé RFD3 (n={(~is_native).sum()})', zorder=4)

        ax_a.axhline(0.0, color='#94a3b8', linestyle='--', linewidth=1.2, zorder=2)
        ax_a.axhline(0.30, color='#10b981', linestyle=':', linewidth=1.2, label='Marge Champion (0.30)', zorder=2)

        ax_a.set_xlabel("Déviation RMSD des Boucles Remodelées RFD3 (Å)", fontsize=10.5, fontweight='bold', color='#1e293b')
        ax_a.set_ylabel("Score d'Orthogonalité F_ortho", fontsize=10.5, fontweight='bold', color='#1e293b')
        ax_a.set_title("5A : Apport de la Co-Adaptation RFD3\n(Induced-Fit vs Squelette Rigide)",
                       fontsize=11.5, fontweight='bold', color='#0f172a', pad=8)
        ax_a.legend(loc='lower right', fontsize=8.0, framealpha=0.92)

    else:
        x_b = _safe_col(df, 'rmsd_b_prime', 0.0)
        y_perf = _safe_col(df, 'f_ortho', 0.0)
        ax_a.scatter(x_b, y_perf, color='#0284c7', edgecolors='#0f172a', s=75, alpha=0.85, zorder=3)
        ax_a.axhline(0.0, color='#94a3b8', linestyle='--', linewidth=1.2)
        ax_a.set_xlabel("RMSD Chaîne B' vs WT (Å)", fontsize=10.5, fontweight='bold', color='#1e293b')
        ax_a.set_ylabel("Score d'Orthogonalité F_ortho", fontsize=10.5, fontweight='bold', color='#1e293b')
        ax_a.set_title("5A : Déviation Structurale B' vs F_ortho\n(Données RFD3 flex_rmsd non fournies)",
                       fontsize=11.5, fontweight='bold', color='#0f172a', pad=8)

    # --- 5B: Interface Biophysical Quality or Parent Motif Breakdown ---
    ax_b = axes[1]
    ax_b.set_facecolor('#f8fafc')
    ax_b.grid(True, linestyle='--', color='#e2e8f0', alpha=0.7, zorder=1)

    has_sasa = 'dsasa_rescue' in df.columns and df['dsasa_rescue'].dropna().any()
    has_unsat = 'unsat_hb_rescue' in df.columns and df['unsat_hb_rescue'].dropna().any()

    if has_sasa and has_unsat:
        x_sasa = df['dsasa_rescue'].fillna(0.0)
        y_unsat = df['unsat_hb_rescue'].fillna(0.0)
        pass_mask = df['passes'].fillna(False).astype(bool) if 'passes' in df.columns else (
            (_safe_col(df, 'f_ortho', 0.0) > 0) & (_safe_col(df, 'dockq', 0.0) >= 0.23)
        )

        max_sasa = max(2500.0, x_sasa.max() + 200.0)
        max_unsat = max(12.0, y_unsat.max() + 2.0)
        ax_b.fill_between([1200.0, max_sasa], 0, 5.0, color='#d1fae5', alpha=0.45,
                          label='Fenêtre Optimale (ΔSASA ≥ 1200 Å² & H-bonds insat. ≤ 5)', zorder=1)
        ax_b.axvline(1200.0, color='#10b981', linestyle='--', linewidth=1.3, zorder=2)
        ax_b.axhline(5.0, color='#10b981', linestyle='--', linewidth=1.3, zorder=2)

        ax_b.scatter(x_sasa[pass_mask], y_unsat[pass_mask], color='#10b981', s=80, edgecolors='#064e3b',
                     label='Orthogonal Validé', zorder=4)
        ax_b.scatter(x_sasa[~pass_mask], y_unsat[~pass_mask], color='#ef4444', marker='x', s=60,
                     label='Rejeté', zorder=3)

        ax_b.set_xlim(min(0.0, x_sasa.min() - 100), max_sasa)
        ax_b.set_ylim(-0.5, max_unsat)
        ax_b.set_xlabel("Surface d'Interface Enfouie ΔSASA (Å²)", fontsize=10.5, fontweight='bold', color='#1e293b')
        ax_b.set_ylabel("Liaisons Hydrogène Enfouies Non Satisfaites", fontsize=10.5, fontweight='bold', color='#1e293b')
        ax_b.set_title("5B : Qualité Biophysique d'Interface (PyRosetta)\n(Risque de dénaturation / agrégation)",
                       fontsize=11.5, fontweight='bold', color='#0f172a', pad=8)
        # Position legend at upper left to avoid collision with points located in upper right
        ax_b.legend(loc='upper left', fontsize=8.0, framealpha=0.92)

    else:
        # Performance by Parent Motif A'
        motifs = df['parent_a_motif'].fillna('Principal').astype(str) if 'parent_a_motif' in df.columns else pd.Series(['Principal']*len(df), index=df.index)
        f_vals = _safe_col(df, 'f_ortho', 0.0)
        motif_stats = pd.DataFrame({'motif': motifs, 'f_ortho': f_vals}).groupby('motif')['f_ortho'].agg(['count', 'mean']).reset_index()

        y_pos = np.arange(len(motif_stats))
        bars = ax_b.barh(y_pos, motif_stats['mean'], color='#38bdf8', edgecolor='#0284c7', height=0.55, zorder=3)
        ax_b.axvline(0.0, color='#94a3b8', linestyle='--', linewidth=1.2, zorder=2)

        for idx, bar in enumerate(bars):
            val = motif_stats['mean'].iloc[idx]
            cnt = motif_stats['count'].iloc[idx]
            ax_b.text(val + (0.02 if val >= 0 else -0.06), bar.get_y() + bar.get_height() / 2,
                      f"F: {val:+.3f} (n={cnt})", va='center', fontsize=8.5, fontweight='bold', color='#0f172a')

        ax_b.set_yticks(y_pos)
        ax_b.set_yticklabels(motif_stats['motif'], fontsize=9.5, fontweight='bold', color='#1e293b')
        ax_b.set_xlabel("Score F_ortho Moyen", fontsize=10.5, fontweight='bold', color='#1e293b')
        ax_b.set_title("5B : Rendement d'Orthogonalité par Motif A' Parent\n(PyRosetta SASA/H-bond non disponible)",
                       fontsize=11.5, fontweight='bold', color='#0f172a', pad=8)

    fig.suptitle("Figure 5 : Analyse Biomécanique RFD3 & Profil Biophysique d'Interface\n"
                 "Évaluation de l'accommodation structurale et de la solvabilité expérimentale",
                 fontsize=12, fontweight='bold', color='#0f172a', y=0.96)

    return _save_and_close_fig(fig, fig_dir, "figure5_induced_fit_biophysics")


# ===========================================================================
# FIGURE 6: Cascade d'Attrition du Pipeline (Multi-Stage Funnel)
# ===========================================================================
def plot_pipeline_attrition_funnel(df, fig_dir=None):
    """Figure 6: Multi-Stage Filtering Funnel (Cascade d'Attrition du Pipeline).

    Tracks attrition from total candidates through each physical & orthogonal filter.
    """
    if df is None or df.empty:
        return _empty_figure("Figure 6 : Cascade d'Attrition du Pipeline",
                            "Aucune donnée de conception disponible", fig_dir, "figure6_attrition_funnel")

    total_eval = len(df)
    pass_early = (df['status'].fillna('ok') == 'ok').sum() if 'status' in df.columns else total_eval

    pass_dockq = (df['dockq'].fillna(0.0) >= 0.23).sum() if 'pass_dockq' not in df.columns else (
        df['pass_dockq'].fillna(False).astype(bool).sum()
    )

    pass_rupture = (df['iptm_rupture'].fillna(1.0) <= 0.35).sum() if 'pass_rupture' not in df.columns else (
        df['pass_rupture'].fillna(False).astype(bool).sum()
    )

    pass_negative = (df['iptm_negative'].fillna(1.0) <= 0.40).sum() if 'pass_negative' not in df.columns else (
        df['pass_negative'].fillna(False).astype(bool).sum()
    )

    pass_rescue = (df['iptm_rescue'].fillna(0.0) >= 0.65).sum() if 'pass_rescue' not in df.columns else (
        df['pass_rescue'].fillna(False).astype(bool).sum()
    )

    has_energy = 'pass_energy' in df.columns and df['pass_energy'].dropna().any()
    pass_energy = df['pass_energy'].fillna(False).astype(bool).sum() if has_energy else None

    pass_final = df['passes'].fillna(False).astype(bool).sum() if 'passes' in df.columns else (
        (_safe_col(df, 'f_ortho', 0.0) > 0) & (_safe_col(df, 'dockq', 0.0) >= 0.23)
    ).sum()

    stages = [
        ("1. Total Designs Évalués", total_eval),
        ("2. Statut Valide (Pass Early-Exit)", pass_early),
        ("3. Pose CAPRI Valide (DockQ ≥ 0.23)", pass_dockq),
        ("4. Rupture A'·B_WT Validée (≤ 0.35)", pass_rupture),
        ("5. Négatif A_WT·B' Validé (≤ 0.40)", pass_negative),
        ("6. Sauvetage A'·B' Validé (≥ 0.65)", pass_rescue),
    ]

    if has_energy:
        stages.append(("7. Énergie PyRosetta (pass_energy)", pass_energy))

    stages.append(("8. Paires Orthogonales Validées (passes)", pass_final))

    labels = [s[0] for s in stages]
    counts = [s[1] for s in stages]

    labels_rev = labels[::-1]
    counts_rev = counts[::-1]

    colors_grad = [
        '#10b981',  # Step 8 / Final (Emerald)
        '#059669',
        '#0284c7',
        '#0369a1',
        '#2563eb',
        '#1d4ed8',
        '#334155',
        '#1e293b'   # Step 1 / Total (Dark Slate)
    ]
    if not has_energy:
        colors_grad = [colors_grad[0]] + colors_grad[2:]

    fig, ax = plt.subplots(figsize=(10.5, 5.8), facecolor='#ffffff')
    plt.subplots_adjust(top=0.88, bottom=0.12)
    ax.set_facecolor('#f8fafc')
    ax.grid(True, axis='x', linestyle='--', color='#e2e8f0', alpha=0.7, zorder=1)

    y_pos = np.arange(len(stages))
    bars = ax.barh(y_pos, counts_rev, color=colors_grad[:len(stages)], edgecolor='#0f172a',
                   height=0.62, zorder=3, alpha=0.92)

    max_c = max(counts) if counts else 1
    ax.set_xlim(0, max_c * 1.30)

    for idx, (bar, cnt) in enumerate(zip(bars, counts_rev)):
        pct_init = (cnt / total_eval * 100.0) if total_eval > 0 else 0.0
        orig_idx = len(stages) - 1 - idx
        prev_cnt = stages[orig_idx - 1][1] if orig_idx > 0 else total_eval
        retention = (cnt / prev_cnt * 100.0) if prev_cnt > 0 else 0.0

        if orig_idx == 0:
            ret_str = "[Init]"
        elif prev_cnt >= cnt and cnt > 0:
            ret_str = f"[-{100 - retention:.0f}%]"
        else:
            ret_str = ""

        label_txt = f"{cnt} / {total_eval} ({pct_init:.1f}%)"
        if ret_str:
            label_txt += f"  {ret_str}"

        ax.text(cnt + max_c * 0.015, bar.get_y() + bar.get_height() / 2,
                label_txt, va='center', fontsize=9.5, fontweight='bold', color='#0f172a')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels_rev, fontsize=10, fontweight='bold', color='#1e293b')
    ax.set_xlabel("Nombre de Candidats Conservés", fontsize=11, fontweight='bold', color='#1e293b', labelpad=8)
    ax.set_title("Figure 6 : Cascade d'Attrition & Entonnoir de Filtrage Multi-Étapes\n"
                 "Identification quantitative des goulots d'étranglement de la spécificité orthogonale",
                 fontsize=12, fontweight='bold', color='#0f172a', pad=14)

    return _save_and_close_fig(fig, fig_dir, "figure6_attrition_funnel")


# ===========================================================================
# STYLED EXCEL EXPORT (Strict 'passes' truth & UTF-8 clean accents)
# ===========================================================================
def generate_styled_excel(csv_path, xlsx_path):
    """Generate professional, multi-tab styled Excel workbook aligning strictly with Module 5 'passes'."""
    if not os.path.exists(csv_path):
        return

    try:
        df = pd.read_csv(csv_path)
    except Exception:
        df = pd.DataFrame()

    wb = openpyxl.Workbook()

    # -------------------------------------------------------------
    # Sheet 1: Executive Orthogonality Dashboard
    # -------------------------------------------------------------
    ws1 = wb.active
    ws1.title = "Executive Summary"
    ws1.views.sheetView[0].showGridLines = True

    # Title Banner
    ws1.merge_cells("A1:K1")
    title_cell = ws1["A1"]
    title_cell.value = "ORTHOPPINTERFACE - RAPPORT D'ORTHOGONALITÉ PROTÉINE-PROTÉINE"
    title_cell.font = Font(name="Segoe UI", size=15, bold=True, color="FFFFFF")
    title_cell.fill = PatternFill(start_color="1E293B", end_color="1E293B", fill_type="solid")
    title_cell.alignment = Alignment(horizontal="center", vertical="center")
    ws1.row_dimensions[1].height = 35

    # Subtitle Banner
    ws1.merge_cells("A2:K2")
    sub_cell = ws1["A2"]
    sub_cell.value = "Évaluation Haute Résolution RoseTTAFold-3 All-Atom (RF3) & Conformité DockQ CAPRI"
    sub_cell.font = Font(name="Segoe UI", size=10, italic=True, color="E2E8F0")
    sub_cell.fill = PatternFill(start_color="334155", end_color="334155", fill_type="solid")
    sub_cell.alignment = Alignment(horizontal="center", vertical="center")
    ws1.row_dimensions[2].height = 20

    # Summary table columns
    cols_summary = [
        ("design_id", "Paire Orthogonale (A'·B')", 40),
        ("parent_a_motif", "Motif A' Parent", 24),
        ("iptm_rescue", "iPTM Sauvetage (A'·B')", 20),
        ("iptm_negative", "iPTM Négatif (A_WT·B')", 20),
        ("iptm_rupture", "iPTM Rupture (A'·B_WT)", 20),
        ("f_ortho", "Score F_ortho", 16),
        ("dockq", "Score DockQ", 14),
        ("dockq_quality", "Qualité CAPRI", 16),
        ("fnat", "F_nat Contacts", 14),
        ("plddt_rescue", "pLDDT Complexe (%)", 18),
        ("verdict", "Verdict Orthogonal", 26),
    ]

    header_fill = PatternFill(start_color="2B4C7E", end_color="2B4C7E", fill_type="solid")
    header_font = Font(name="Segoe UI", size=11, bold=True, color="FFFFFF")
    thin_border = Border(
        left=Side(style='thin', color='D0D7DE'),
        right=Side(style='thin', color='D0D7DE'),
        top=Side(style='thin', color='D0D7DE'),
        bottom=Side(style='thin', color='D0D7DE')
    )

    ws1.row_dimensions[4].height = 28
    for col_idx, (_, label, width) in enumerate(cols_summary, start=1):
        cell = ws1.cell(row=4, column=col_idx, value=label)
        cell.font = header_font
        cell.fill = header_fill
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = thin_border
        col_letter = get_column_letter(col_idx)
        ws1.column_dimensions[col_letter].width = width

    fill_champ = PatternFill(start_color="D1FAE5", end_color="D1FAE5", fill_type="solid")
    fill_pass = PatternFill(start_color="F0FDF4", end_color="F0FDF4", fill_type="solid")
    fill_alt = PatternFill(start_color="F8FAFC", end_color="F8FAFC", fill_type="solid")
    fill_white = PatternFill(start_color="FFFFFF", end_color="FFFFFF", fill_type="solid")

    font_bold = Font(name="Segoe UI", size=10, bold=True)
    font_regular = Font(name="Segoe UI", size=10)
    font_code = Font(name="Consolas", size=9, bold=True)

    for row_idx, row in enumerate(df.to_dict('records')):
        r_num = row_idx + 5
        ws1.row_dimensions[r_num].height = 20

        f_ortho_val = float(row.get('f_ortho', 0.0)) if pd.notnull(row.get('f_ortho')) else 0.0
        dockq_val = row.get('dockq', 0.0)
        dockq_q = row.get('dockq_quality', 'N/A')
        fnat_val = row.get('fnat', 0.0)
        plddt_val = row.get('plddt_rescue', 0.0)
        status = str(row.get('status', 'ok'))

        # STRICT ALIGNMENT ON 'passes'
        if 'passes' in df.columns:
            passes = bool(row.get('passes', False))
        else:
            pass_dockq = bool(row.get('pass_dockq', (pd.notnull(dockq_val) and dockq_val >= 0.23)))
            passes = bool(f_ortho_val > 0.0 and pass_dockq)

        is_ortho = passes
        is_champ = passes and (f_ortho_val >= 0.30)
        is_off_target = "off_target" in status or (pd.notnull(dockq_val) and dockq_val < 0.23)

        row_fill = fill_champ if is_champ else (fill_pass if is_ortho else (fill_alt if row_idx % 2 == 1 else fill_white))

        if is_off_target:
            status_text = f"HORS-CIBLE (DockQ {dockq_val if pd.notnull(dockq_val) else 0.0:.3f})"
        elif is_champ:
            status_text = "CHAMPION ORTHOGONAL"
        elif is_ortho:
            status_text = "ORTHOGONAL"
        else:
            status_text = "NON ORTHOGONAL"

        row_values = [
            row.get('design_id', ''),
            row.get('parent_a_motif', ''),
            row.get('iptm_rescue', 0.0),
            row.get('iptm_negative', 0.0),
            row.get('iptm_rupture', 0.0),
            f_ortho_val,
            dockq_val if pd.notnull(dockq_val) else 0.0,
            dockq_q if pd.notnull(dockq_q) else 'N/A',
            fnat_val if pd.notnull(fnat_val) else 0.0,
            plddt_val if pd.notnull(plddt_val) else 0.0,
            status_text
        ]

        for col_idx, val in enumerate(row_values, start=1):
            cell = ws1.cell(row=r_num, column=col_idx, value=val)
            cell.fill = row_fill
            cell.border = thin_border
            cell.font = font_code if col_idx == 1 else (font_bold if col_idx in [3, 6, 11] else font_regular)

            if col_idx in [3, 4, 5, 6, 7, 9, 10]:
                cell.alignment = Alignment(horizontal="right", vertical="center")
                if col_idx in [3, 4, 5, 6, 7]:
                    cell.number_format = "0.0000"
                elif col_idx == 9:
                    cell.number_format = "0.000"
                elif col_idx == 10:
                    cell.number_format = "0.0"
            else:
                cell.alignment = Alignment(horizontal="center" if col_idx in [2, 8, 11] else "left", vertical="center")

    # -------------------------------------------------------------
    # Sheet 2: Raw High-Resolution All-Atom Metrics
    # -------------------------------------------------------------
    ws2 = wb.create_sheet(title="Raw All-Atom Metrics")
    ws2.views.sheetView[0].showGridLines = True

    for col_idx, col_name in enumerate(df.columns, start=1):
        cell = ws2.cell(row=1, column=col_idx, value=col_name)
        cell.font = header_font
        cell.fill = header_fill
        cell.alignment = Alignment(horizontal="center", vertical="center")
        cell.border = thin_border
        col_letter = get_column_letter(col_idx)
        ws2.column_dimensions[col_letter].width = max(len(col_name) + 4, 14)

    for row_idx, row in enumerate(df.to_dict('records')):
        for col_idx, col_name in enumerate(df.columns, start=1):
            val = row[col_name]
            cell = ws2.cell(row=row_idx + 2, column=col_idx, value=val)
            cell.font = font_regular
            cell.border = thin_border
            if isinstance(val, float):
                cell.number_format = "0.0000"

    wb.save(xlsx_path)
    print(f"[Styled Excel Export] Saved clean interactive spreadsheet -> {xlsx_path}")


# ===========================================================================
# INTERACTIVE SCIENTIFIC HTML DASHBOARD (With Inline SVGs & Rich Controls)
# ===========================================================================
def generate_interactive_html(csv_path, html_path):
    """Generate modern, publication-grade interactive HTML dashboard.

    Integrates inline SVG figures, executive synthesis takeaways, podium top 3,
    parent motif statistical aggregation, and complete interactive table.
    """
    if not os.path.exists(csv_path):
        return

    try:
        df = pd.read_csv(csv_path)
    except Exception:
        df = pd.DataFrame()

    fig_dir = os.path.join(os.path.dirname(os.path.abspath(html_path)), "figures")

    coherence_path = os.path.join(os.path.dirname(os.path.abspath(csv_path)), "coherence.json")

    # -----------------------------------------------------------------------
    # Generate 6 Scientific Figures
    # -----------------------------------------------------------------------
    svg_fig1 = plot_cross_reactivity_matrix(df, fig_dir=fig_dir)
    svg_fig2 = plot_specificity_landscape(df, fig_dir=fig_dir)
    svg_fig3 = plot_interaction_triad_distributions(df, fig_dir=fig_dir)
    svg_fig4 = plot_structural_funnel(df, fig_dir=fig_dir, coherence_path=coherence_path)
    svg_fig5 = plot_induced_fit_and_biophysics(df, fig_dir=fig_dir)
    svg_fig6 = plot_pipeline_attrition_funnel(df, fig_dir=fig_dir)

    # -----------------------------------------------------------------------
    # Statistics & Executive Key Takeaways
    # -----------------------------------------------------------------------
    total_designs = len(df)
    ortho_count = int(df['passes'].fillna(False).astype(bool).sum()) if 'passes' in df.columns else (
        int(((_safe_col(df, 'f_ortho', 0.0) > 0) & (_safe_col(df, 'dockq', 0.0) >= 0.23)).sum())
    )
    pass_rate = (ortho_count / total_designs * 100.0) if total_designs > 0 else 0.0

    max_iptm_r = df['iptm_rescue'].max() if ('iptm_rescue' in df.columns and not df.empty) else 0.0
    max_f_ortho = df['f_ortho'].max() if ('f_ortho' in df.columns and not df.empty) else 0.0
    max_dockq = df['dockq'].max() if ('dockq' in df.columns and not df.empty) else 0.0

    # Ceiling health diagnostic
    ceilings = df['iptm_ceiling'].dropna() if 'iptm_ceiling' in df.columns else pd.Series(dtype=float)
    med_ceiling = ceilings.median() if len(ceilings) > 0 else np.nan
    if pd.notnull(med_ceiling) and med_ceiling < 0.50:
        ceiling_alert = (f"Alerte Plafond MSA : iPTM médian du sauvage = {med_ceiling:.2f} (< 0.50). "
                         "Le complexe natif est difficilement reconnu sous ce régime MSA ; "
                         "privilégiez F_rel ou les critères structuraux CAPRI.")
    elif pd.notnull(med_ceiling):
        ceiling_alert = f"Plafond MSA Sain : iPTM sauvage médian = {med_ceiling:.2f} (reconnaissance robuste)."
    else:
        ceiling_alert = "Plafond MSA : contrôle sauvage non exécuté ou indisponible."

    # Primary cause of rejection
    rejection_causes = {}
    for row in df.to_dict('records'):
        p = row.get('passes', False) if 'passes' in df.columns else (row.get('f_ortho', 0) > 0)
        if not p:
            dq = row.get('dockq', 0.0)
            if pd.notnull(dq) and dq < 0.23:
                rejection_causes['Pose Hors-cible CAPRI (DockQ < 0.23)'] = rejection_causes.get('Pose Hors-cible CAPRI (DockQ < 0.23)', 0) + 1
            elif row.get('pass_rupture') is False or row.get('iptm_rupture', 0) > 0.35:
                rejection_causes["Échec Rupture A'·B_WT (affinité WT résiduelle > 0.35)"] = rejection_causes.get("Échec Rupture A'·B_WT (affinité WT résiduelle > 0.35)", 0) + 1
            elif row.get('pass_negative') is False or row.get('iptm_negative', 0) > 0.40:
                rejection_causes["Échec Négatif A_WT·B' (croisement hors-cible > 0.40)"] = rejection_causes.get("Échec Négatif A_WT·B' (croisement hors-cible > 0.40)", 0) + 1
            elif row.get('pass_rescue') is False or row.get('iptm_rescue', 0) < 0.65:
                rejection_causes["Échec Sauvetage A'·B' (affinité orthogonale < 0.65)"] = rejection_causes.get("Échec Sauvetage A'·B' (affinité orthogonale < 0.65)", 0) + 1
            elif row.get('pass_energy') is False:
                rejection_causes["Échec Énergie PyRosetta (ΔG répulsif)"] = rejection_causes.get("Échec Énergie PyRosetta (ΔG répulsif)", 0) + 1
            else:
                rejection_causes["Marge F_ortho insuffisante (≤ 0)"] = rejection_causes.get("Marge F_ortho insuffisante (≤ 0)", 0) + 1

    top_rejection = max(rejection_causes.items(), key=lambda x: x[1])[0] if rejection_causes else "Aucun rejet (100% de succès)"

    # Best parent motif
    motifs = df['parent_a_motif'].fillna('Principal').astype(str) if 'parent_a_motif' in df.columns else pd.Series(['Principal']*len(df), index=df.index)
    motif_summary = []
    for m_name, grp in df.groupby(motifs):
        m_total = len(grp)
        m_pass = int(grp['passes'].fillna(False).astype(bool).sum()) if 'passes' in df.columns else int((grp['f_ortho'] > 0).sum())
        m_rate = (m_pass / m_total * 100.0) if m_total > 0 else 0.0
        m_dockq = grp['dockq'].mean() if ('dockq' in grp.columns and pd.notnull(grp['dockq'].mean())) else 0.0
        m_fortho = grp['f_ortho'].mean() if ('f_ortho' in grp.columns and pd.notnull(grp['f_ortho'].mean())) else 0.0
        m_max_rescue = grp['iptm_rescue'].max() if ('iptm_rescue' in grp.columns and pd.notnull(grp['iptm_rescue'].max())) else 0.0
        motif_summary.append({
            'motif': m_name, 'total': m_total, 'passed': m_pass, 'rate': m_rate,
            'dockq': m_dockq, 'fortho': m_fortho, 'max_rescue': m_max_rescue
        })
    motif_summary.sort(key=lambda x: (x['passed'], x['fortho']), reverse=True)
    best_motif = motif_summary[0]['motif'] if motif_summary else "N/A"

    # Champion Candidate
    sort_keys = [c for c in ['passes', 'f_ortho', 'dockq'] if c in df.columns]
    if 'passes' in df.columns and df['passes'].any():
        champ_row = df[df['passes'].fillna(False).astype(bool)].sort_values(
            by=sort_keys, ascending=[False]*len(sort_keys)
        ).iloc[0]
    elif not df.empty and sort_keys:
        champ_row = df.sort_values(by=sort_keys, ascending=[False]*len(sort_keys)).iloc[0]
    elif not df.empty:
        champ_row = df.iloc[0]
    else:
        champ_row = {}

    champ_id = champ_row.get('design_id', 'N/A')

    # -----------------------------------------------------------------------
    # Top 3 Podium Cards HTML
    # -----------------------------------------------------------------------
    podium_cards_html = []
    top_sort_keys = [c for c in ['passes', 'f_ortho', 'dockq'] if c in df.columns]
    top_3_df = df.sort_values(
        by=top_sort_keys,
        ascending=[False]*len(top_sort_keys)
    ).head(3) if top_sort_keys else df.head(3)

    medals = [("🥇 1er - Champion", "gold"), ("🥈 2e - Dauphin", "silver"), ("🥉 3e - Sélectionné", "bronze")]
    for idx, (_, row) in enumerate(top_3_df.iterrows()):
        medal_title, medal_class = medals[idx] if idx < 3 else (f"#{idx+1}", "neutral")
        p_val = bool(row.get('passes', False)) if 'passes' in df.columns else (row.get('f_ortho', 0) > 0)
        p_badge = '<span class="badge badge-champion">Champion Orthogonal</span>' if (p_val and row.get('f_ortho', 0) >= 0.30) else (
            '<span class="badge badge-success">Orthogonal</span>' if p_val else '<span class="badge badge-negative">Non Orthogonal</span>'
        )

        card = f"""
        <div class="podium-card {medal_class}">
            <div class="podium-header">
                <span class="medal-tag">{medal_title}</span>
                {p_badge}
            </div>
            <div class="podium-body">
                <h3 class="podium-title"><code>{row.get('design_id', 'N/A')}</code></h3>
                <div class="podium-meta">Motif A' : <b>{row.get('parent_a_motif', 'Co-adaptation')}</b> | Mutations : <b>{row.get('n_mut_A', '?')} (A') / {row.get('n_mut_B', '?')} (B')</b></div>
                <div class="metrics-grid">
                    <div class="metric-box">
                        <span class="m-label">Score F_ortho</span>
                        <span class="m-val" style="color: var(--green); font-size: 16px;">{row.get('f_ortho', 0.0):+.4f}</span>
                    </div>
                    <div class="metric-box">
                        <span class="m-label">DockQ (CAPRI)</span>
                        <span class="m-val" style="color: var(--accent); font-size: 16px;">{row.get('dockq', 0.0):.4f} <small>({row.get('dockq_quality', 'N/A')})</small></span>
                    </div>
                    <div class="metric-box">
                        <span class="m-label">iPTM Sauvetage</span>
                        <span class="m-val" style="color: #10b981;">{row.get('iptm_rescue', 0.0):.4f}</span>
                    </div>
                    <div class="metric-box">
                        <span class="m-label">iPTM Rupture</span>
                        <span class="m-val" style="color: #ef4444;">{row.get('iptm_rupture', 0.0):.4f}</span>
                    </div>
                    <div class="metric-box">
                        <span class="m-label">iPTM Négatif</span>
                        <span class="m-val" style="color: #f97316;">{row.get('iptm_negative', 0.0):.4f}</span>
                    </div>
                    <div class="metric-box">
                        <span class="m-label">pLDDT Complexe</span>
                        <span class="m-val">{row.get('plddt_rescue', 0.0):.1f}%</span>
                    </div>
                </div>
            </div>
        </div>
        """
        podium_cards_html.append(card)

    # -----------------------------------------------------------------------
    # Parent Motif Breakdown Table HTML
    # -----------------------------------------------------------------------
    motif_rows_html = []
    for m in motif_summary:
        m_row = f"""
        <tr>
            <td><b>{m['motif']}</b></td>
            <td>{m['total']}</td>
            <td><b style="color: var(--green);">{m['passed']}</b></td>
            <td>
                <div style="display: flex; align-items: center; gap: 8px;">
                    <div class="progress-bar-bg" style="width: 80px; margin-top: 0;">
                        <div class="progress-bar-fill green" style="width: {m['rate']:.1f}%"></div>
                    </div>
                    <span>{m['rate']:.1f}%</span>
                </div>
            </td>
            <td>{m['dockq']:.3f}</td>
            <td style="font-weight: 700; color: {'#10b981' if m['fortho'] >= 0.20 else '#38bdf8'};">{m['fortho']:+.3f}</td>
            <td><b>{m['max_rescue']:.3f}</b></td>
        </tr>
        """
        motif_rows_html.append(m_row)

    # -----------------------------------------------------------------------
    # Interactive Table Rows HTML
    # -----------------------------------------------------------------------
    rows_html = []
    dockq_pass_count = 0
    offtarget_count = 0

    for idx, row in enumerate(df.to_dict('records')):
        iptm_r = float(row.get('iptm_rescue', 0.0))
        iptm_neg = float(row.get('iptm_negative', 0.0))
        iptm_rup = float(row.get('iptm_rupture', 0.0))
        f_ortho = float(row.get('f_ortho', 0.0))
        dockq_val = row.get('dockq', 0.0)
        dockq_q = row.get('dockq_quality', 'N/A')
        f_nat = float(row.get('fnat', 0.0)) if pd.notnull(row.get('fnat')) else 0.0
        plddt = float(row.get('plddt_rescue', 0.0))
        status = str(row.get('status', 'ok'))

        # Verdict
        if 'passes' in df.columns:
            passes = bool(row.get('passes', False))
        else:
            pass_dockq = bool(row.get('pass_dockq', (pd.notnull(dockq_val) and dockq_val >= 0.23)))
            passes = bool(f_ortho > 0.0 and pass_dockq)

        is_ortho = passes
        is_champ = passes and (f_ortho >= 0.30)
        is_off_target = "off_target" in status or (pd.notnull(dockq_val) and dockq_val < 0.23)

        if pd.notnull(dockq_val) and dockq_val >= 0.23:
            dockq_pass_count += 1
        if is_off_target:
            offtarget_count += 1

        # Diagnostic failure reason
        diag_reason = "Validé"
        if is_off_target:
            diag_reason = f"Hors-cible CAPRI (DockQ {dockq_val if pd.notnull(dockq_val) else 0.0:.3f} < 0.23)"
            status_badge = f'<span class="badge badge-negative" style="background: rgba(239, 68, 68, 0.15); color: #f87171; border: 1px solid rgba(239, 68, 68, 0.3);">Hors-cible (DockQ {dockq_val if pd.notnull(dockq_val) else 0.0:.3f})</span>'
        elif is_champ:
            status_badge = '<span class="badge badge-champion">Champion Orthogonal</span>'
        elif is_ortho:
            status_badge = '<span class="badge badge-success">Orthogonal</span>'
        else:
            status_badge = '<span class="badge badge-negative">Non Orthogonal</span>'
            if row.get('pass_rupture') is False or iptm_rup > 0.35:
                diag_reason = f"Échec Rupture A'·B_WT ({iptm_rup:.3f} > 0.35)"
            elif row.get('pass_negative') is False or iptm_neg > 0.40:
                diag_reason = f"Échec Négatif A_WT·B' ({iptm_neg:.3f} > 0.40)"
            elif row.get('pass_rescue') is False or iptm_r < 0.65:
                diag_reason = f"Échec Sauvetage ({iptm_r:.3f} < 0.65)"
            elif row.get('pass_energy') is False:
                diag_reason = "Échec Énergie PyRosetta (ΔG répulsif)"
            else:
                diag_reason = f"Marge F_ortho insuffisante ({f_ortho:+.3f})"

        row_filter_class = "row-passed" if passes else ("row-offtarget" if is_off_target else "row-failed")
        if pd.notnull(dockq_val) and dockq_val >= 0.23:
            row_filter_class += " row-dockq-ok"

        row_str = f"""
        <tr class="{row_filter_class}">
            <td style="font-weight: 700; color: var(--text-muted);">{idx+1}</td>
            <td><code>{row.get('design_id', '')}</code></td>
            <td><span class="strategy-tag">{row.get('parent_a_motif', 'Co-adaptation')}</span></td>
            <td>
                <div class="metric-val" style="color: #10b981;">{iptm_r:.4f}</div>
                <div class="progress-bar-bg"><div class="progress-bar-fill green" style="width: {min(iptm_r*100, 100):.1f}%"></div></div>
            </td>
            <td>
                <div class="metric-val" style="color: #38bdf8;">{iptm_neg:.4f}</div>
                <div class="progress-bar-bg"><div class="progress-bar-fill blue" style="width: {min(iptm_neg*100, 100):.1f}%"></div></div>
            </td>
            <td>
                <div class="metric-val" style="color: #f87171;">{iptm_rup:.4f}</div>
            </td>
            <td style="font-weight: 700; font-size: 14px; color: {'#10b981' if f_ortho >= 0.30 else ('#38bdf8' if is_ortho else '#ef4444')};">{f_ortho:+.4f}</td>
            <td><b style="color: {'#10b981' if (pd.notnull(dockq_val) and dockq_val >= 0.49) else ('#f59e0b' if (pd.notnull(dockq_val) and dockq_val >= 0.23) else '#ef4444')};">{dockq_val if pd.notnull(dockq_val) else 0.0:.4f}</b> <small style="color: var(--text-muted);">({dockq_q})</small></td>
            <td>{f_nat:.3f}</td>
            <td><b>{plddt:.1f}%</b></td>
            <td>
                {status_badge}
                <div style="font-size: 11px; color: var(--text-muted); margin-top: 4px;">{diag_reason}</div>
            </td>
        </tr>
        """
        rows_html.append(row_str)

    html_content = f"""<!DOCTYPE html>
<html lang="fr">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>OrthoPPInterface - Tableau de Bord Scientifique & Visualisations de Publication</title>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;600&display=swap" rel="stylesheet">
    <style>
        :root {{
            --bg: #0f172a;
            --card-bg: #1e293b;
            --card-sub: #162032;
            --border: #334155;
            --border-light: #475569;
            --text: #f8fafc;
            --text-muted: #94a3b8;
            --accent: #38bdf8;
            --green: #10b981;
            --blue: #3b82f6;
            --yellow: #f59e0b;
            --red: #ef4444;
            --purple: #8b5cf6;
        }}
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background-color: var(--bg);
            color: var(--text);
            padding: 30px;
            line-height: 1.5;
        }}
        .header {{
            background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%);
            border: 1px solid var(--border);
            border-radius: 16px;
            padding: 28px;
            margin-bottom: 24px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.35);
        }}
        .header h1 {{
            font-size: 25px;
            font-weight: 700;
            color: var(--text);
            display: flex;
            align-items: center;
            gap: 12px;
        }}
        .header p {{
            color: var(--text-muted);
            margin-top: 6px;
            font-size: 14px;
        }}

        /* Executive Takeaways Box */
        .takeaways-box {{
            background: rgba(30, 41, 59, 0.7);
            border: 1px solid var(--accent);
            border-left: 5px solid var(--accent);
            border-radius: 12px;
            padding: 20px 24px;
            margin-bottom: 24px;
        }}
        .takeaways-title {{
            font-size: 14px;
            text-transform: uppercase;
            letter-spacing: 0.06em;
            color: var(--accent);
            font-weight: 700;
            margin-bottom: 10px;
        }}
        .takeaways-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
            gap: 14px;
            font-size: 13.5px;
        }}
        .takeaway-item {{
            background: rgba(15, 23, 42, 0.6);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 10px 14px;
        }}
        .takeaway-item b {{ color: var(--text); }}

        /* Stats Grid */
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
            gap: 16px;
            margin-bottom: 28px;
        }}
        .stat-card {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 12px;
            padding: 20px;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.2);
        }}
        .stat-card .label {{
            font-size: 11.5px;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: var(--text-muted);
            font-weight: 600;
        }}
        .stat-card .val {{
            font-size: 26px;
            font-weight: 700;
            margin-top: 6px;
            font-family: 'JetBrains Mono', monospace;
        }}

        /* Podium Cards */
        .section-title {{
            font-size: 18px;
            font-weight: 700;
            color: var(--text);
            margin: 28px 0 16px 0;
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        .podium-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
            gap: 18px;
            margin-bottom: 30px;
        }}
        .podium-card {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 14px;
            padding: 20px;
            box-shadow: 0 4px 16px rgba(0, 0, 0, 0.25);
            position: relative;
        }}
        .podium-card.gold {{ border-color: #f59e0b; background: linear-gradient(145deg, #1e293b 70%, rgba(245, 158, 11, 0.08) 100%); }}
        .podium-card.silver {{ border-color: #94a3b8; background: linear-gradient(145deg, #1e293b 70%, rgba(148, 163, 184, 0.08) 100%); }}
        .podium-card.bronze {{ border-color: #d97706; background: linear-gradient(145deg, #1e293b 70%, rgba(217, 119, 6, 0.08) 100%); }}
        .podium-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 12px;
        }}
        .medal-tag {{
            font-size: 13px;
            font-weight: 700;
            color: var(--text);
        }}
        .podium-title {{
            font-size: 15px;
            font-weight: 700;
            color: var(--accent);
            margin-bottom: 6px;
        }}
        .podium-meta {{
            font-size: 12px;
            color: var(--text-muted);
            margin-bottom: 14px;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 10px;
        }}
        .metric-box {{
            background: var(--card-sub);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 8px 10px;
        }}
        .metric-box .m-label {{
            font-size: 10.5px;
            text-transform: uppercase;
            color: var(--text-muted);
            display: block;
        }}
        .metric-box .m-val {{
            font-family: 'JetBrains Mono', monospace;
            font-weight: 700;
            font-size: 13.5px;
            margin-top: 2px;
            display: block;
        }}

        /* Scientific Figures Grid */
        .figures-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(540px, 1fr));
            gap: 22px;
            margin-bottom: 32px;
        }}
        .figure-card {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 14px;
            padding: 20px;
            box-shadow: 0 4px 16px rgba(0, 0, 0, 0.25);
            display: flex;
            flex-direction: column;
        }}
        .figure-card-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 10px;
            border-bottom: 1px solid var(--border);
            padding-bottom: 10px;
        }}
        .figure-card-title {{
            font-size: 14px;
            font-weight: 700;
            color: var(--text);
        }}
        .figure-desc {{
            font-size: 12px;
            color: var(--text-muted);
            margin-bottom: 14px;
            line-height: 1.4;
        }}
        .figure-svg-wrapper {{
            background-color: #ffffff;
            border-radius: 10px;
            padding: 12px;
            display: flex;
            justify-content: center;
            align-items: center;
            overflow-x: auto;
            flex: 1;
        }}
        .figure-svg-wrapper svg {{
            max-width: 100%;
            height: auto;
            display: block;
        }}
        .figure-actions {{
            display: flex;
            justify-content: flex-end;
            gap: 10px;
            margin-top: 12px;
        }}
        .btn-download {{
            font-size: 11px;
            font-weight: 600;
            color: var(--accent);
            background: rgba(56, 189, 248, 0.1);
            border: 1px solid rgba(56, 189, 248, 0.25);
            border-radius: 6px;
            padding: 5px 10px;
            text-decoration: none;
            transition: all 0.2s ease;
        }}
        .btn-download:hover {{
            background: rgba(56, 189, 248, 0.2);
            color: #ffffff;
        }}

        /* Tables */
        .table-container {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 12px;
            overflow-x: auto;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.2);
            margin-bottom: 30px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            text-align: left;
            font-size: 12.5px;
        }}
        th {{
            background-color: #0f172a;
            color: var(--text-muted);
            font-weight: 600;
            padding: 13px 15px;
            border-bottom: 1px solid var(--border);
            text-transform: uppercase;
            font-size: 11px;
            letter-spacing: 0.05em;
            cursor: pointer;
            user-select: none;
            white-space: nowrap;
        }}
        th:hover {{ color: var(--accent); }}
        td {{
            padding: 12px 15px;
            border-bottom: 1px solid var(--border);
            color: var(--text);
            vertical-align: middle;
        }}
        tr:hover td {{ background-color: rgba(255, 255, 255, 0.025); }}
        code {{
            font-family: 'JetBrains Mono', monospace;
            background-color: #0f172a;
            padding: 3px 6px;
            border-radius: 6px;
            color: #38bdf8;
            font-size: 11.5px;
        }}
        .strategy-tag {{
            background-color: rgba(56, 189, 248, 0.1);
            color: #38bdf8;
            border: 1px solid rgba(56, 189, 248, 0.2);
            padding: 3px 8px;
            border-radius: 20px;
            font-size: 11px;
            font-weight: 500;
            white-space: nowrap;
        }}
        .metric-val {{ font-weight: 600; font-family: 'JetBrains Mono', monospace; font-size: 12.5px; }}
        .progress-bar-bg {{
            background-color: #334155;
            height: 5px;
            border-radius: 3px;
            width: 100%;
            margin-top: 5px;
            overflow: hidden;
        }}
        .progress-bar-fill.green {{ background-color: #10b981; height: 100%; }}
        .progress-bar-fill.blue {{ background-color: #3b82f6; height: 100%; }}
        .badge {{
            padding: 3px 8px;
            border-radius: 20px;
            font-size: 10.5px;
            font-weight: 600;
            display: inline-block;
            white-space: nowrap;
        }}
        .badge-champion {{ background-color: rgba(16, 185, 129, 0.2); color: #10b981; border: 1px solid rgba(16, 185, 129, 0.4); }}
        .badge-success {{ background-color: rgba(56, 189, 248, 0.15); color: #38bdf8; border: 1px solid rgba(56, 189, 248, 0.3); }}
        .badge-negative {{ background-color: rgba(239, 68, 68, 0.15); color: #f87171; border: 1px solid rgba(239, 68, 68, 0.3); }}

        /* Filter Controls */
        .controls-bar {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 12px;
            margin-bottom: 16px;
        }}
        .filter-btn-group {{
            display: flex;
            gap: 8px;
            flex-wrap: wrap;
        }}
        .filter-btn {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            color: var(--text-muted);
            padding: 8px 14px;
            border-radius: 8px;
            font-size: 12px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s ease;
        }}
        .filter-btn:hover {{ border-color: var(--accent); color: var(--text); }}
        .filter-btn.active {{
            background-color: var(--accent);
            color: #0f172a;
            border-color: var(--accent);
        }}
        .search-input {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 8px 16px;
            color: var(--text);
            font-size: 13px;
            outline: none;
            width: 280px;
        }}
        .search-input:focus {{ border-color: var(--accent); }}
    </style>
</head>
<body>
    <div class="header">
        <h1>OrthoPPInterface - Rapport d'Orthogonalité Protéine-Protéine</h1>
        <p>Conception Rationnelle & Évaluation RoseTTAFold-3 All-Atom (RF3) GPU | Conformité CAPRI DockQ | Visualisations Vectorielles</p>
    </div>

    <!-- Executive Synthesis Diagnostic Box -->
    <div class="takeaways-box">
        <div class="takeaways-title">Synthèse Diagnostique &amp; Clés d'Analyse (Executive Key Takeaways)</div>
        <div class="takeaways-grid">
            <div class="takeaway-item">
                <b>Taux de Succès Global :</b> {ortho_count} / {total_designs} designs ({pass_rate:.1f}%) satisfont l'ensemble des critères d'orthogonalité stricte.
            </div>
            <div class="takeaway-item">
                <b>Santé Contrôle Sauvage :</b> {ceiling_alert}
            </div>
            <div class="takeaway-item">
                <b>Cause n°1 d'Attrition :</b> {top_rejection}.
            </div>
            <div class="takeaway-item">
                <b>Motif A' Optimal :</b> Motif <code>{best_motif}</code> présente le meilleur rendement de sauvetage.
            </div>
            <div class="takeaway-item">
                <b>Candidat Champion :</b> Paire <code>{champ_id}</code> retenue pour l'analyse structurale haute résolution.
            </div>
        </div>
    </div>

    <!-- Stats Grid -->
    <div class="stats-grid">
        <div class="stat-card">
            <div class="label">Paires Évaluées</div>
            <div class="val" style="color: var(--accent);">{total_designs}</div>
        </div>
        <div class="stat-card">
            <div class="label">Paires Orthogonales Valides</div>
            <div class="val" style="color: var(--green);">{ortho_count} <small style="font-size: 14px; color: var(--text-muted);">({pass_rate:.1f}%)</small></div>
        </div>
        <div class="stat-card">
            <div class="label">Affinité Sauvetage Max (iPTM)</div>
            <div class="val" style="color: var(--green);">{max_iptm_r:.4f}</div>
        </div>
        <div class="stat-card">
            <div class="label">Score F_ortho Maximal</div>
            <div class="val" style="color: var(--green);">+{max_f_ortho:.4f}</div>
        </div>
        <div class="stat-card">
            <div class="label">DockQ Maximal (CAPRI)</div>
            <div class="val" style="color: var(--accent);">{max_dockq:.4f}</div>
        </div>
    </div>

    <!-- Podium Section -->
    <div class="section-title">🏆 Podium des 3 Meilleurs Candidats</div>
    <div class="podium-grid">
        {''.join(podium_cards_html)}
    </div>

    <!-- Scientific Figures Section -->
    <div class="section-title">📊 Figures Scientifiques de Publication (Matplotlib Agg Vectorielles)</div>
    <div class="figures-grid">
        <!-- Figure 1 -->
        <div class="figure-card">
            <div class="figure-card-header">
                <span class="figure-card-title">Figure 1 : Matrices de Réactivité Croisée 2×2 (Boyken et al. 2016)</span>
            </div>
            <div class="figure-desc">Démontre visuellement la diagonalisation de la matrice de liaison (sauvetage A'·B' et contrôle WT·WT forts, extinction des hybrides croisés A'·B_WT et A_WT·B').</div>
            <div class="figure-svg-wrapper">
                {svg_fig1}
            </div>
            <div class="figure-actions">
                <a href="figures/figure1_cross_reactivity.svg" download class="btn-download">Télécharger SVG</a>
                <a href="figures/figure1_cross_reactivity.png" download class="btn-download">Télécharger PNG (300 DPI)</a>
            </div>
        </div>

        <!-- Figure 2 -->
        <div class="figure-card">
            <div class="figure-card-header">
                <span class="figure-card-title">Figure 2 : Paysage de Spécificité &amp; Isolignes de Marge</span>
            </div>
            <div class="figure-desc">Nuage 2D on-target vs pire off-target sauvage avec diagonale de neutralité y = x, isolignes de marge F_ortho et plafond WT.</div>
            <div class="figure-svg-wrapper">
                {svg_fig2}
            </div>
            <div class="figure-actions">
                <a href="figures/figure2_specificity_landscape.svg" download class="btn-download">Télécharger SVG</a>
                <a href="figures/figure2_specificity_landscape.png" download class="btn-download">Télécharger PNG (300 DPI)</a>
            </div>
        </div>

        <!-- Figure 3 -->
        <div class="figure-card">
            <div class="figure-card-header">
                <span class="figure-card-title">Figure 3 : Distributions de la Triade d'Interaction</span>
            </div>
            <div class="figure-desc">Boxplots et trajectoires individuelles reliant chaque candidat validé le long de la vallée d'extinction cross-réactive (RF3 iPTM et ipSAE Dunbrack).</div>
            <div class="figure-svg-wrapper">
                {svg_fig3}
            </div>
            <div class="figure-actions">
                <a href="figures/figure3_triad_distributions.svg" download class="btn-download">Télécharger SVG</a>
                <a href="figures/figure3_triad_distributions.png" download class="btn-download">Télécharger PNG (300 DPI)</a>
            </div>
        </div>

        <!-- Figure 4 -->
        <div class="figure-card">
            <div class="figure-card-header">
                <span class="figure-card-title">Figure 4 : Entonnoir Structural CAPRI &amp; Auto-Cohérence</span>
            </div>
            <div class="figure-desc">Puits de repliement natif (&lt; 5 Å) face aux poses hors-cible et corrélation de concordance Spearman entre IA Structurale (RF3) et Physique Explicite (PyRosetta).</div>
            <div class="figure-svg-wrapper">
                {svg_fig4}
            </div>
            <div class="figure-actions">
                <a href="figures/figure4_structural_funnel.svg" download class="btn-download">Télécharger SVG</a>
                <a href="figures/figure4_structural_funnel.png" download class="btn-download">Télécharger PNG (300 DPI)</a>
            </div>
        </div>

        <!-- Figure 5 -->
        <div class="figure-card">
            <div class="figure-card-header">
                <span class="figure-card-title">Figure 5 : Induced-Fit RFD3 &amp; Qualité Biophysique d'Interface</span>
            </div>
            <div class="figure-desc">Apport de la plasticité du squelette RFD3 face au squelette natif rigide, et profil atomique de surface enfouie ΔSASA vs liaisons H insatisfaites.</div>
            <div class="figure-svg-wrapper">
                {svg_fig5}
            </div>
            <div class="figure-actions">
                <a href="figures/figure5_induced_fit_biophysics.svg" download class="btn-download">Télécharger SVG</a>
                <a href="figures/figure5_induced_fit_biophysics.png" download class="btn-download">Télécharger PNG (300 DPI)</a>
            </div>
        </div>

        <!-- Figure 6 -->
        <div class="figure-card">
            <div class="figure-card-header">
                <span class="figure-card-title">Figure 6 : Cascade d'Attrition du Pipeline</span>
            </div>
            <div class="figure-desc">Entonnoir quantitatif montrant le taux de rétention des designs étape par étape à travers les filtres géométriques, négatifs et énergétiques.</div>
            <div class="figure-svg-wrapper">
                {svg_fig6}
            </div>
            <div class="figure-actions">
                <a href="figures/figure6_attrition_funnel.svg" download class="btn-download">Télécharger SVG</a>
                <a href="figures/figure6_attrition_funnel.png" download class="btn-download">Télécharger PNG (300 DPI)</a>
            </div>
        </div>
    </div>

    <!-- Parent Motif Breakdown Table Section -->
    <div class="section-title">🧬 Rendement par Motif Parent A'</div>
    <div class="table-container">
        <table>
            <thead>
                <tr>
                    <th>Motif A' Parent</th>
                    <th>Designs Évalués</th>
                    <th>Validés (passes)</th>
                    <th>Taux de Succès</th>
                    <th>DockQ Moyen</th>
                    <th>F_ortho Moyen</th>
                    <th>iPTM Sauvetage Max</th>
                </tr>
            </thead>
            <tbody>
                {''.join(motif_rows_html)}
            </tbody>
        </table>
    </div>

    <!-- Complete Cohort Table Section -->
    <div class="section-title">📋 Cohorte Complète &amp; Diagnostics Détaillés</div>
    <div class="controls-bar">
        <div class="filter-btn-group">
            <button class="filter-btn active" onclick="setQuickFilter('all', this)">Tous ({total_designs})</button>
            <button class="filter-btn" onclick="setQuickFilter('passed', this)">Validés (passes) ({ortho_count})</button>
            <button class="filter-btn" onclick="setQuickFilter('dockq', this)">DockQ ≥ 0.23 ({dockq_pass_count})</button>
            <button class="filter-btn" onclick="setQuickFilter('offtarget', this)">Hors-cible ({offtarget_count})</button>
        </div>
        <input type="text" id="searchInput" class="search-input" placeholder="Filtrer par nom ou motif..." onkeyup="applyFilter()">
    </div>

    <div class="table-container">
        <table id="orthoTable">
            <thead>
                <tr>
                    <th onclick="sortTable(0)">#</th>
                    <th onclick="sortTable(1)">Paire Orthogonale (A'·B')</th>
                    <th onclick="sortTable(2)">Motif A' Parent</th>
                    <th onclick="sortTable(3)">iPTM Sauvetage (A'·B')</th>
                    <th onclick="sortTable(4)">iPTM Négatif (A_WT·B')</th>
                    <th onclick="sortTable(5)">iPTM Rupture (A'·B_WT)</th>
                    <th onclick="sortTable(6)">Score F_ortho</th>
                    <th onclick="sortTable(7)">Score DockQ</th>
                    <th onclick="sortTable(8)">F_nat</th>
                    <th onclick="sortTable(9)">pLDDT Complexe</th>
                    <th onclick="sortTable(10)">Verdict &amp; Diagnostic</th>
                </tr>
            </thead>
            <tbody>
                {''.join(rows_html)}
            </tbody>
        </table>
    </div>

    <script>
        let currentQuickFilter = 'all';

        function setQuickFilter(type, btn) {{
            currentQuickFilter = type;
            document.querySelectorAll('.filter-btn').forEach(b => b.classList.remove('active'));
            if (btn) btn.classList.add('active');
            applyFilter();
        }}

        function applyFilter() {{
            const searchVal = document.getElementById('searchInput').value.toUpperCase();
            const table = document.getElementById('orthoTable');
            const tr = table.getElementsByTagName('tr');

            for (let i = 1; i < tr.length; i++) {{
                const row = tr[i];
                let matchesQuick = true;
                if (currentQuickFilter === 'passed') {{
                    matchesQuick = row.classList.contains('row-passed');
                }} else if (currentQuickFilter === 'dockq') {{
                    matchesQuick = row.classList.contains('row-dockq-ok');
                }} else if (currentQuickFilter === 'offtarget') {{
                    matchesQuick = row.classList.contains('row-offtarget');
                }}

                const text = row.textContent || row.innerText;
                const matchesSearch = text.toUpperCase().indexOf(searchVal) > -1;

                row.style.display = (matchesQuick && matchesSearch) ? "" : "none";
            }}
        }}

        function sortTable(n) {{
            const table = document.getElementById("orthoTable");
            let rows, switching, i, x, y, shouldSwitch, dir, switchcount = 0;
            switching = true;
            dir = "desc";
            while (switching) {{
                switching = false;
                rows = table.rows;
                for (i = 1; i < (rows.length - 1); i++) {{
                    shouldSwitch = false;
                    x = rows[i].getElementsByTagName("TD")[n];
                    y = rows[i + 1].getElementsByTagName("TD")[n];
                    let xVal = parseFloat(x.innerText) || x.innerText.toLowerCase();
                    let yVal = parseFloat(y.innerText) || y.innerText.toLowerCase();
                    if (dir === "asc") {{
                        if (xVal > yVal) {{ shouldSwitch = true; break; }}
                    }} else if (dir === "desc") {{
                        if (xVal < yVal) {{ shouldSwitch = true; break; }}
                    }}
                }}
                if (shouldSwitch) {{
                    rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
                    switching = true;
                    switchcount++;
                }} else {{
                    if (switchcount === 0 && dir === "desc") {{
                        dir = "asc";
                        switching = true;
                    }}
                }}
            }}
        }}
    </script>
</body>
</html>
"""
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"[Interactive Scientific HTML Export] Saved clean visual dashboard -> {html_path}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Export Styled Orthogonality Reports & Scientific Publication Figures.")
    parser.add_argument("csv_path", help="Path to orthogonality_scores.csv")
    parser.add_argument("--out-dir", "-o", default=None, help="Output directory (defaults to dirname of csv_path)")
    args = parser.parse_args()

    csv_p = args.csv_path
    out_d = args.out_dir or os.path.dirname(os.path.abspath(csv_p))
    os.makedirs(out_d, exist_ok=True)

    xlsx_p = os.path.join(out_d, "SUMMARY.xlsx")
    html_p = os.path.join(out_d, "SUMMARY.html")

    print(f"Generating styled Excel -> {xlsx_p}")
    generate_styled_excel(csv_p, xlsx_p)

    print(f"Generating interactive HTML & figures -> {html_p}")
    generate_interactive_html(csv_p, html_p)
    print("Export complete!")
