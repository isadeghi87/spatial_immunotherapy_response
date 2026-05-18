#!/usr/bin/env python3
"""
imc_04_validation_strips.py

Statistical validation and strip/violin plots comparing Responders vs
Non-Responders on IMC-derived features.

Tests: Mann-Whitney U (two-sided), BH correction across features within
each feature group.

Output figures:
  out/imc_validation/figures/imc_strips_composition.{png,pdf}
  out/imc_validation/figures/imc_strips_functional.{png,pdf}
  out/imc_validation/figures/imc_strips_density.{png,pdf}

Output tables:
  out/imc_validation/tables/imc_mwu_results.csv
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')

# ── Paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)

sys.path.insert(0, PROJECT_DIR)
from config.publication_style import apply_style, PAL, save_fig

TABLE_DIR  = os.path.join(PROJECT_DIR, 'out', '03_validation', 'imc_validation', 'tables')
FIG_DIR    = os.path.join(PROJECT_DIR, 'out', '03_validation', 'imc_validation', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

FEAT_PATH = os.path.join(TABLE_DIR, 'imc_patient_features.csv')
MWU_PATH  = os.path.join(TABLE_DIR, 'imc_mwu_results.csv')

# ── Feature groups for separate panels ───────────────────────────────────
FEATURE_GROUPS = {
    'Composition': {
        'features': [
            'frac_CD8_T', 'frac_CD4_T', 'frac_Treg', 'frac_NK',
            'frac_B_cell', 'frac_Macrophage', 'frac_APC_DC', 'frac_Tumor',
            'Treg_CD4_ratio',
        ],
        'labels': {
            'frac_CD8_T':       'CD8 T (%)',
            'frac_CD4_T':       'CD4 T (%)',
            'frac_Treg':        'Treg (%)',
            'frac_NK':          'NK (%)',
            'frac_B_cell':      'B cell (%)',
            'frac_Macrophage':  'Macrophage (%)',
            'frac_APC_DC':      'APC / DC (%)',
            'frac_Tumor':       'Tumor (%)',
            'Treg_CD4_ratio':   'Treg / (CD4 + Treg)',
        },
        'scale_pct': {'frac_CD8_T', 'frac_CD4_T', 'frac_Treg', 'frac_NK',
                      'frac_B_cell', 'frac_Macrophage', 'frac_APC_DC', 'frac_Tumor'},
    },
    'Functional': {
        'features': [
            'CD8_pd1_frac', 'CD8_gzmb_frac', 'CD8_tcf1_frac',
            'CD8_tox_frac', 'CD8_tim3_frac', 'CD8_cd39_frac',
            'CD8_tex_like_frac', 'CD8_tbet_frac', 'CD8_cd103_frac',
            'CD4_CXCR5_frac', 'Tumor_PDL1_frac', 'Mac_PDL1_mean',
            'Mac_CD163_mean',
        ],
        'labels': {
            'CD8_pd1_frac':      'CD8 PD1+ (%)',
            'CD8_gzmb_frac':     'CD8 GranzymeB+ (%)',
            'CD8_tcf1_frac':     'CD8 TCF1+ (%)',
            'CD8_tox_frac':      'CD8 TOX+ (%)',
            'CD8_tim3_frac':     'CD8 Tim3+ (%)',
            'CD8_cd39_frac':     'CD8 CD39+ (%)',
            'CD8_tex_like_frac': 'CD8 TEX-like (%)',
            'CD8_tbet_frac':     'CD8 Tbet+ (%)',
            'CD8_cd103_frac':    'CD8 CD103+ (%)',
            'CD4_CXCR5_frac':    'CD4 CXCR5+ (%)',
            'Tumor_PDL1_frac':   'Tumor PDL1+ (%)',
            'Mac_PDL1_mean':     'Mac PDL1 (mean)',
            'Mac_CD163_mean':    'Mac CD163 (mean)',
        },
        'scale_pct': {'CD8_pd1_frac', 'CD8_gzmb_frac', 'CD8_tcf1_frac',
                      'CD8_tox_frac', 'CD8_tim3_frac', 'CD8_cd39_frac',
                      'CD8_tex_like_frac', 'CD8_tbet_frac', 'CD8_cd103_frac',
                      'CD4_CXCR5_frac', 'Tumor_PDL1_frac'},
    },
    'Density': {
        'features': [
            'CD8_density_mm2', 'immune_density_mm2', 'tumor_density_mm2',
            'CD8_parenchyma_mm2', 'CD8_stroma_mm2',
            'CD8_parenchyma_stroma_ratio',
        ],
        'labels': {
            'CD8_density_mm2':              'CD8 density (cells/mm²)',
            'immune_density_mm2':           'Immune density (cells/mm²)',
            'tumor_density_mm2':            'Tumor density (cells/mm²)',
            'CD8_parenchyma_mm2':           'CD8 parenchyma (cells/mm²)',
            'CD8_stroma_mm2':               'CD8 stroma (cells/mm²)',
            'CD8_parenchyma_stroma_ratio':  'CD8 paren./stroma ratio',
        },
        'scale_pct': set(),
    },
}

RESP_COLORS = {
    'Responder':     PAL['Responder'],
    'Non-Responder': PAL['Non-Responder'],
}
RESP_ORDER = ['Responder', 'Non-Responder']


def run_mwu_per_feature(df: pd.DataFrame,
                         features: list,
                         group_col: str = 'response_group') -> pd.DataFrame:
    """Mann-Whitney U test for each feature, BH-corrected within group."""
    r  = df[df[group_col] == 'Responder']
    nr = df[df[group_col] == 'Non-Responder']

    rows = []
    for feat in features:
        x_r  = r[feat].dropna().values
        x_nr = nr[feat].dropna().values

        n_r  = len(x_r)
        n_nr = len(x_nr)

        if n_r < 3 or n_nr < 3:
            rows.append({
                'feature': feat,
                'n_R': n_r, 'n_NR': n_nr,
                'median_R': np.nanmedian(x_r) if n_r > 0 else np.nan,
                'median_NR': np.nanmedian(x_nr) if n_nr > 0 else np.nan,
                'pval': np.nan, 'qval': np.nan,
                'direction': 'undetermined',
                'status': 'insufficient_n',
            })
            continue

        stat, pval = mannwhitneyu(x_r, x_nr, alternative='two-sided')

        med_r  = np.median(x_r)
        med_nr = np.median(x_nr)
        direction = 'higher_in_R' if med_r > med_nr else 'higher_in_NR'

        rows.append({
            'feature': feat,
            'n_R': n_r, 'n_NR': n_nr,
            'median_R':  med_r,
            'median_NR': med_nr,
            'pval': pval,
            'direction': direction,
            'status': 'tested',
        })

    result_df = pd.DataFrame(rows)
    tested = result_df[result_df['status'] == 'tested']

    if len(tested) > 0:
        reject, qvals, _, _ = multipletests(tested['pval'].values,
                                             method='fdr_bh')
        result_df.loc[result_df['status'] == 'tested', 'qval'] = qvals
    else:
        result_df['qval'] = np.nan

    return result_df


def sig_stars(pval: float, qval: float) -> str:
    """Return significance annotation string based on BH-corrected q-value."""
    if np.isnan(pval) or np.isnan(qval):
        return ''
    if qval < 0.001:
        return '***'
    if qval < 0.01:
        return '**'
    if qval < 0.05:
        return '*'
    if pval < 0.05:
        return '†'   # nominal only (uncorrected)
    return 'ns'


def make_strip_figure(df: pd.DataFrame,
                       features: list,
                       labels: dict,
                       scale_pct: set,
                       mwu_df: pd.DataFrame,
                       title: str,
                       ncols: int = 3) -> plt.Figure:
    """
    Create a multi-panel strip + violin plot, one panel per feature.
    """
    apply_style()
    n_feat = len(features)
    nrows  = int(np.ceil(n_feat / ncols))

    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols * 2.3, nrows * 3.1),
                              squeeze=False)
    axes_flat = axes.flatten()

    for ax_idx, feat in enumerate(features):
        ax = axes_flat[ax_idx]

        if feat not in df.columns:
            ax.set_visible(False)
            continue

        plot_df = df[['response_group', feat]].dropna()
        y_col   = feat + '_pct' if feat in scale_pct else feat

        if feat in scale_pct:
            plot_df = plot_df.copy()
            plot_df[y_col] = plot_df[feat] * 100

        # Violin
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            sns.violinplot(data=plot_df, x='response_group', y=y_col,
                           order=RESP_ORDER, palette=RESP_COLORS,
                           inner=None, linewidth=0.8, cut=0, ax=ax,
                           alpha=0.35)

        # Strip
        sns.stripplot(data=plot_df, x='response_group', y=y_col,
                      order=RESP_ORDER, palette=RESP_COLORS,
                      size=3.5, alpha=0.85, jitter=True, ax=ax, zorder=3)

        # Median bar
        for xi, grp in enumerate(RESP_ORDER):
            vals = plot_df[plot_df['response_group'] == grp][y_col].values
            if len(vals) > 0:
                med = np.median(vals)
                ax.plot([xi - 0.2, xi + 0.2], [med, med],
                        color='black', lw=1.5, zorder=4)

        # Significance annotation
        stat_row = mwu_df[mwu_df['feature'] == feat]
        pv = stat_row['pval'].values[0] if len(stat_row) else np.nan
        qv = stat_row['qval'].values[0] if len(stat_row) else np.nan
        # rank-biserial from MWU u-stat if available, else compute
        n_r_pts  = (plot_df['response_group'] == 'Responder').sum()
        n_nr_pts = (plot_df['response_group'] == 'Non-Responder').sum()

        is_bh_sig = (not np.isnan(qv)) and qv < 0.05
        is_nom    = (not np.isnan(pv)) and pv < 0.05 and not is_bh_sig

        # Highlight significant panels
        if is_bh_sig:
            ax.set_facecolor('#fff8e1')

        stars = sig_stars(pv, qv) if len(stat_row) else ''
        ymax   = plot_df[y_col].max()
        ymin   = plot_df[y_col].min()
        yrange = max(ymax - ymin, 1e-9)

        bar_y2 = ymax + yrange * 0.08
        bar_y  = ymax + yrange * 0.14
        star_y = bar_y + yrange * 0.02
        r_y    = bar_y + yrange * 0.15
        p_y    = bar_y + yrange * 0.28

        y_top = ymax + yrange * (0.70 if is_bh_sig else 0.52)
        ax.set_ylim(ymin - yrange * 0.05, y_top)

        # Bracket
        if stars:
            bc = '#c0392b' if is_bh_sig else ('#8a6d00' if is_nom else '0.45')
            ax.plot([0, 0, 1, 1], [bar_y2, bar_y, bar_y, bar_y2],
                    lw=0.9, color='0.30' if is_bh_sig else '0.45')

            # Stars label — bold/red for BH-sig, italic amber for nominal
            ax.text(0.5, star_y, stars,
                    ha='center', va='bottom',
                    fontsize=11 if is_bh_sig else 9,
                    fontweight='bold' if is_bh_sig else 'normal',
                    fontstyle='normal' if is_bh_sig else 'italic',
                    color='#c0392b' if is_bh_sig else ('#8a6d00' if is_nom else '0.5'))

            # Rank-biserial r (compute from medians direction only here)
            med_r  = np.median(plot_df[plot_df['response_group'] == 'Responder'][y_col].values)
            med_nr = np.median(plot_df[plot_df['response_group'] == 'Non-Responder'][y_col].values)
            r_rb   = stat_row['direction'].values[0] if len(stat_row) else ''
            # compute actual rank-biserial if u-stat available; approximate from direction
            r_vals  = plot_df[plot_df['response_group'] == 'Responder'][y_col].dropna().values
            nr_vals = plot_df[plot_df['response_group'] == 'Non-Responder'][y_col].dropna().values
            if len(r_vals) >= 2 and len(nr_vals) >= 2:
                u_stat, _ = mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
                rb_val = float(2 * u_stat / (len(r_vals) * len(nr_vals)) - 1)
            else:
                rb_val = 0.0
            ax.text(0.5, r_y, f'r = {rb_val:+.2f}',
                    ha='center', va='bottom', fontsize=7,
                    fontstyle='italic',
                    color='#c0392b' if is_bh_sig else '0.50')

            # BH p-value — only for BH-significant panels
            if is_bh_sig:
                p_str = 'p_BH < 0.001' if qv < 0.001 else f'p_BH = {qv:.3f}'
                ax.text(0.5, p_y, p_str,
                        ha='center', va='bottom', fontsize=6.5,
                        fontweight='semibold', color='#c0392b')

        # Panel styling
        ax.set_title(labels.get(feat, feat), fontsize=9.5, fontweight='bold', pad=4)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['R', 'NR'], fontsize=9, fontweight='semibold')
        sns.despine(ax=ax)

        # Orange spines for BH-significant panels
        if is_bh_sig:
            for sp in ['left', 'bottom']:
                ax.spines[sp].set_edgecolor('#e67e22')
                ax.spines[sp].set_linewidth(1.8)

        # Sample counts below axis
        for xi, grp in enumerate(RESP_ORDER):
            n = (plot_df['response_group'] == grp).sum()
            ax.text(xi, ax.get_ylim()[0] - 0.08 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
                    f'n={n}', ha='center', va='top', fontsize=6.5,
                    transform=ax.get_xaxis_transform())

    # Hide unused axes
    for ax_idx in range(len(features), len(axes_flat)):
        axes_flat[ax_idx].set_visible(False)

    # Legend
    handles = [mpatches.Patch(facecolor=RESP_COLORS[k], label=k)
               for k in RESP_ORDER]
    fig.legend(handles=handles, loc='lower right',
               bbox_to_anchor=(1.0, 0.01), ncol=1, fontsize=8,
               frameon=False)

    # Split title into two lines: main title + stats note
    parts = title.split('  (')
    if len(parts) == 2:
        line1 = parts[0].strip()
        line2 = '(' + parts[1].strip() + '   |   Mann–Whitney U + BH correction'
    else:
        line1 = title
        line2 = 'Mann–Whitney U + BH correction'
    fig.suptitle(f'{line1}\n{line2}',
                 fontsize=10, fontweight='bold', y=1.03, linespacing=1.6)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    return fig


def main():
    apply_style()
    print(f"Loading patient features from {FEAT_PATH}")
    df = pd.read_csv(FEAT_PATH)
    print(f"  {len(df)} patients  |  columns: {len(df.columns)}")

    # Only patients with known response
    df = df[df['response_binary'].notna()].copy()
    df['response_group'] = df['response_binary'].map(
        {1: 'Responder', 0: 'Non-Responder'})
    df = df[df['response_group'].notna()]

    n_R  = (df['response_group'] == 'Responder').sum()
    n_NR = (df['response_group'] == 'Non-Responder').sum()
    print(f"  Responders: {n_R}  |  Non-Responders: {n_NR}")

    # ── Run MWU for all features ──────────────────────────────────────
    all_mwu = []
    for group_name, group_def in FEATURE_GROUPS.items():
        feats = [f for f in group_def['features'] if f in df.columns]
        if not feats:
            continue
        mwu = run_mwu_per_feature(df, feats)
        mwu['feature_group'] = group_name
        all_mwu.append(mwu)

    mwu_df = pd.concat(all_mwu, ignore_index=True)
    mwu_df.to_csv(MWU_PATH, index=False)
    print(f"\nMWU results saved: {MWU_PATH}")

    # Print significant features
    sig = mwu_df[mwu_df['qval'] < 0.05].sort_values('qval')
    nom = mwu_df[(mwu_df['pval'] < 0.05) & (mwu_df['qval'] >= 0.05)]
    print(f"\nFDR-significant features (q<0.05):  {len(sig)}")
    if len(sig):
        print(sig[['feature', 'median_R', 'median_NR', 'pval', 'qval', 'direction']].to_string(index=False))
    print(f"\nNominal-only features (p<0.05, q≥0.05): {len(nom)}")
    if len(nom):
        print(nom[['feature', 'pval', 'qval', 'direction']].to_string(index=False))

    # ── Generate figures ──────────────────────────────────────────────
    for group_name, group_def in FEATURE_GROUPS.items():
        feats = [f for f in group_def['features'] if f in df.columns]
        if not feats:
            print(f"  Skipping {group_name}: no features present")
            continue

        group_mwu = mwu_df[mwu_df['feature_group'] == group_name]

        fig = make_strip_figure(
            df=df,
            features=feats,
            labels=group_def['labels'],
            scale_pct=group_def['scale_pct'],
            mwu_df=group_mwu,
            title=f'IMC Validation — {group_name}  (R n={n_R}, NR n={n_NR})',
            ncols=3,
        )
        out_stem = os.path.join(FIG_DIR, f'imc_strips_{group_name.lower()}')
        save_fig(fig, out_stem, formats=('png', 'pdf'))
        print(f"  Saved: {out_stem}.{{png,pdf}}")

    print("\nDone.")


if __name__ == '__main__':
    main()
