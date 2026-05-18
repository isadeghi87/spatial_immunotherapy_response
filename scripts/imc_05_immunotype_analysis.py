#!/usr/bin/env python3
"""
imc_05_immunotype_analysis.py

Three analyses to reproduce and validate the core findings of Salie et al. Gut 2025:

  A. Immunotype × response association
     - Contingency table: paper immunotype (E/C/D) × RECIST response (PR vs SD/PD)
     - Fisher's exact test (2×2: Enriched vs non-Enriched)
     - Chi-squared test across 3×2 table
     - Stacked bar chart

  B. IMC-derived features by paper immunotype
     - Key features: CD8_density_mm2, CD8_parenchyma_stroma_ratio, frac_CD8_T,
       frac_Macrophage, frac_Tumor, CD8_pd1_frac, CD8_gzmb_frac
     - Kruskal-Wallis test across E/C/D groups
     - Pairwise Dunn's test with BH correction
     - Strip + violin plots, colored by immunotype

  C. Survival analysis: KM curves for OS and PFS
     - Stratified by paper immunotype (E/C/D)
     - Stratified by response (R vs NR)
     - Restricted to survival_eligible patients (Child-Pugh < 7)
     - Log-rank tests with Bonferroni correction
     - Median OS/PFS with 95% CI from KM

Outputs:
  out/imc_validation/figures/imc_immunotype_response_bar.{png,pdf}
  out/imc_validation/figures/imc_features_by_immunotype.{png,pdf}
  out/imc_validation/figures/imc_km_os_immunotype.{png,pdf}
  out/imc_validation/figures/imc_km_pfs_immunotype.{png,pdf}
  out/imc_validation/figures/imc_km_os_response.{png,pdf}
  out/imc_validation/figures/imc_km_pfs_response.{png,pdf}
  out/imc_validation/tables/imc_contingency_results.csv
  out/imc_validation/tables/imc_kw_results.csv
  out/imc_validation/tables/imc_km_results.csv

Notes:
  - CP discrepancy: paper reports n=5 excluded (CP≥7); our CSV has n=7 with CP≥7.
    ICI_05(CP=12), ICI_07(CP=11), ICI_12(CP=10) may be PDF mis-extractions.
    We apply CP<7 strictly (n=7 excluded from survival, n=42 for response/composition).
  - n=5 Responders severely limits statistical power; all statistics are exploratory.
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
from scipy.stats import fisher_exact, chi2_contingency, kruskal
from statsmodels.stats.multitest import multipletests
from itertools import combinations

warnings.filterwarnings('ignore')

# ── Paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)

sys.path.insert(0, PROJECT_DIR)
from config.publication_style import apply_style, PAL, save_fig

TABLE_DIR = os.path.join(PROJECT_DIR, 'out', '03_validation', 'imc_validation', 'tables')
FIG_DIR   = os.path.join(PROJECT_DIR, 'out', '03_validation', 'imc_validation', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

FEAT_PATH = os.path.join(TABLE_DIR, 'imc_patient_features.csv')

IMMUNOTYPE_ORDER = ['Enriched', 'Compartmentalized', 'Depleted']
IMMUNOTYPE_COLORS = {
    'Enriched':          '#0077BB',
    'Compartmentalized': '#EE7733',
    'Depleted':          '#CC3311',
}

FEATURES_BY_IT = [
    ('CD8_density_mm2',             'CD8 density (cells/mm²)'),
    ('CD8_parenchyma_stroma_ratio', 'CD8 paren./stroma ratio'),
    ('CD8_parenchyma_mm2',          'CD8 parenchyma (cells/mm²)'),
    ('frac_CD8_T',                  'CD8 T (frac.)'),
    ('frac_Macrophage',             'Macrophage (frac.)'),
    ('frac_Tumor',                  'Tumor (frac.)'),
    ('CD8_pd1_frac',                'CD8 PD1+ (frac.)'),
    ('CD8_gzmb_frac',               'CD8 GranzymeB+ (frac.)'),
    ('CD8_tox_frac',                'CD8 TOX+ (frac.)'),
]


# ══════════════════════════════════════════════════════════════════════════
# Part A: Immunotype × response association
# ══════════════════════════════════════════════════════════════════════════

def run_contingency(df: pd.DataFrame) -> pd.DataFrame:
    """Fisher's exact and chi2 tests on immunotype × response."""
    df_valid = df[df['immunotype_paper'].isin(IMMUNOTYPE_ORDER)].copy()
    df_valid['is_PR'] = df_valid['response_binary'] == 1

    rows = []

    # 3×2 chi2
    ct = pd.crosstab(df_valid['immunotype_paper'], df_valid['is_PR'])
    ct = ct.reindex(index=IMMUNOTYPE_ORDER, columns=[False, True], fill_value=0)

    chi2, p_chi2, dof, _ = chi2_contingency(ct)
    rows.append({
        'test': 'chi2_3x2',
        'comparison': 'E vs C vs D',
        'n_R': ct[True].sum(), 'n_NR': ct[False].sum(),
        'statistic': chi2, 'df': dof, 'pval': p_chi2,
        'note': 'Chi2 across all 3 immunotypes'
    })

    # 2×2 Fisher: Enriched vs non-Enriched
    e_pr  = ct.loc['Enriched', True]
    e_nr  = ct.loc['Enriched', False]
    ne_pr = ct.loc[['Compartmentalized', 'Depleted'], True].sum()
    ne_nr = ct.loc[['Compartmentalized', 'Depleted'], False].sum()
    table_2x2 = [[e_pr, e_nr], [ne_pr, ne_nr]]
    _, p_fisher = fisher_exact(table_2x2, alternative='greater')
    rows.append({
        'test': 'fisher_Evs_CE+D',
        'comparison': 'Enriched vs (C+D)',
        'n_R': e_pr + ne_pr, 'n_NR': e_nr + ne_nr,
        'statistic': np.nan, 'df': 1, 'pval': p_fisher,
        'note': f'2x2: E=[{e_pr}PR,{e_nr}NR] vs C+D=[{ne_pr}PR,{ne_nr}NR]'
    })

    # Individual Fisher for each immunotype vs rest
    for it in IMMUNOTYPE_ORDER:
        rest = [x for x in IMMUNOTYPE_ORDER if x != it]
        it_pr  = ct.loc[it, True]
        it_nr  = ct.loc[it, False]
        rest_pr = ct.loc[rest, True].sum()
        rest_nr = ct.loc[rest, False].sum()
        _, pv = fisher_exact([[it_pr, it_nr], [rest_pr, rest_nr]],
                              alternative='greater')
        rows.append({
            'test': f'fisher_{it}_vs_rest',
            'comparison': f'{it} vs rest',
            'n_R': it_pr + rest_pr, 'n_NR': it_nr + rest_nr,
            'statistic': np.nan, 'df': 1, 'pval': pv,
            'note': f'{it}=[{it_pr}PR,{it_nr}NR] vs rest=[{rest_pr}PR,{rest_nr}NR]'
        })

    print("\n=== Part A: Immunotype × Response ===")
    print(f"Contingency table (paper immunotype × RECIST PR):")
    print(ct.to_string())
    print(f"\nPR rate by immunotype:")
    for it in IMMUNOTYPE_ORDER:
        pr = ct.loc[it, True]; nr = ct.loc[it, False]; total = pr + nr
        pct = 100 * pr / total if total > 0 else 0
        print(f"  {it:<20s}: {pr}/{total} = {pct:.0f}% PR")

    result_df = pd.DataFrame(rows)
    print(f"\nStatistics:")
    print(result_df[['test', 'pval', 'note']].to_string(index=False))
    return result_df, ct


def plot_immunotype_response_bar(ct: pd.DataFrame) -> plt.Figure:
    """Stacked bar: proportion PR/NR by immunotype."""
    apply_style()
    fig, ax = plt.subplots(figsize=(3.5, 3.2))

    totals = ct.sum(axis=1)
    pr_pct = 100 * ct[True] / totals
    nr_pct = 100 * ct[False] / totals

    x = np.arange(len(IMMUNOTYPE_ORDER))
    colors_pr = [IMMUNOTYPE_COLORS[it] for it in IMMUNOTYPE_ORDER]

    ax.bar(x, pr_pct, color=colors_pr, width=0.55, label='PR', alpha=0.9)
    ax.bar(x, nr_pct, bottom=pr_pct, color=colors_pr, width=0.55,
           alpha=0.25, label='SD/PD', edgecolor='none')

    # Annotate counts
    for i, it in enumerate(IMMUNOTYPE_ORDER):
        pr = ct.loc[it, True]; total = int(totals[it])
        ax.text(i, pr_pct[it] + 2, f'{pr}/{total}',
                ha='center', va='bottom', fontsize=7.5)

    ax.set_xticks(x)
    ax.set_xticklabels(IMMUNOTYPE_ORDER, rotation=30, ha='right', fontsize=8)
    ax.set_ylabel('Patients (%)', fontsize=8)
    ax.set_ylim(0, 115)
    ax.set_title('Response rate by immunotype\n(IMC ICI cohort, n=42)', fontsize=8)

    handles = [mpatches.Patch(facecolor='#555555', alpha=0.9, label='PR (Responder)'),
               mpatches.Patch(facecolor='#555555', alpha=0.25, label='SD/PD')]
    ax.legend(handles=handles, fontsize=7, loc='upper right', frameon=False)

    sns.despine(ax=ax)
    fig.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════════
# Part B: IMC features by paper immunotype
# ══════════════════════════════════════════════════════════════════════════

def dunn_test(df: pd.DataFrame, feature: str, group_col: str,
              groups: list) -> pd.DataFrame:
    """
    Dunn's pairwise test (z-statistic approximation via rank sums).
    Conservative version using Bonferroni correction.
    """
    from scipy.stats import rankdata
    all_vals = df[feature].dropna().values
    all_groups = df.loc[df[feature].notna(), group_col].values

    n_total = len(all_vals)
    ranks = rankdata(all_vals)

    rows = []
    for g1, g2 in combinations(groups, 2):
        mask1 = all_groups == g1
        mask2 = all_groups == g2
        n1 = mask1.sum(); n2 = mask2.sum()
        if n1 < 2 or n2 < 2:
            rows.append({'g1': g1, 'g2': g2, 'pval': np.nan})
            continue
        mean_rank1 = ranks[mask1].mean()
        mean_rank2 = ranks[mask2].mean()
        # Variance of difference in mean ranks (Dunn 1964)
        # Using tie-corrected formula
        var = (n_total * (n_total + 1) / 12) * (1/n1 + 1/n2)
        z = (mean_rank1 - mean_rank2) / np.sqrt(var)
        from scipy.stats import norm
        pval = 2 * norm.sf(abs(z))
        rows.append({'g1': g1, 'g2': g2, 'z': z, 'pval': pval})

    result = pd.DataFrame(rows)
    if len(result) > 0 and 'pval' in result.columns:
        valid = result['pval'].notna()
        if valid.sum() > 0:
            _, qvals, _, _ = multipletests(
                result.loc[valid, 'pval'].values, method='bonferroni')
            result.loc[valid, 'qval_bonf'] = qvals
    return result


def run_kruskal_by_immunotype(df: pd.DataFrame) -> pd.DataFrame:
    """Kruskal-Wallis + pairwise Dunn across paper immunotype groups."""
    df_valid = df[df['immunotype_paper'].isin(IMMUNOTYPE_ORDER)].copy()
    groups = IMMUNOTYPE_ORDER

    kw_rows = []
    dunn_all = []

    print("\n=== Part B: Features by Immunotype (Kruskal-Wallis) ===")
    for feat, label in FEATURES_BY_IT:
        if feat not in df_valid.columns:
            continue
        samples = [df_valid[df_valid['immunotype_paper'] == g][feat].dropna().values
                   for g in groups]
        ns = [len(s) for s in samples]
        medians = [np.median(s) if len(s) else np.nan for s in samples]

        if all(n >= 2 for n in ns):
            stat, pval = kruskal(*samples)
        else:
            stat, pval = np.nan, np.nan

        kw_rows.append({
            'feature': feat,
            'label': label,
            'n_E': ns[0], 'n_C': ns[1], 'n_D': ns[2],
            'median_E': medians[0], 'median_C': medians[1], 'median_D': medians[2],
            'KW_stat': stat, 'pval': pval,
        })

        # Dunn's test
        dunn = dunn_test(df_valid, feat, 'immunotype_paper', groups)
        dunn['feature'] = feat
        dunn_all.append(dunn)

        med_str = '  '.join([f'{g[:1]}={m:.3f}' for g, m in zip(groups, medians)])
        print(f"  {feat:<34s}  KW p={pval:.3f}  {med_str}")

    kw_df = pd.DataFrame(kw_rows)
    # BH correction across features
    valid = kw_df['pval'].notna()
    if valid.sum() > 0:
        _, qvals, _, _ = multipletests(kw_df.loc[valid, 'pval'].values, method='fdr_bh')
        kw_df.loc[valid, 'qval'] = qvals

    dunn_df = pd.concat(dunn_all, ignore_index=True) if dunn_all else pd.DataFrame()

    print(f"\n  Nominally significant (p<0.05): "
          f"{(kw_df['pval'] < 0.05).sum()}/{len(kw_df)}")
    return kw_df, dunn_df


def plot_features_by_immunotype(df: pd.DataFrame,
                                 kw_df: pd.DataFrame,
                                 dunn_df: pd.DataFrame) -> plt.Figure:
    """Strip + violin per feature, colored by immunotype."""
    apply_style()
    df_valid = df[df['immunotype_paper'].isin(IMMUNOTYPE_ORDER)].copy()

    n_feat = len(FEATURES_BY_IT)
    ncols = 3
    nrows = int(np.ceil(n_feat / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols * 2.5, nrows * 3.2),
                              squeeze=False)
    axes_flat = axes.flatten()

    for ax_idx, (feat, label) in enumerate(FEATURES_BY_IT):
        ax = axes_flat[ax_idx]
        if feat not in df_valid.columns:
            ax.set_visible(False)
            continue

        plot_df = df_valid[['immunotype_paper', feat]].dropna()

        # Determine significance before plotting
        kw_row = kw_df[kw_df['feature'] == feat]
        pv  = kw_row['pval'].values[0]  if len(kw_row) else np.nan
        qv  = kw_row['qval'].values[0]  if (len(kw_row) and 'qval' in kw_row.columns) else np.nan
        is_bh_sig  = (not np.isnan(pv) if not np.isnan(pv) else False) and (not np.isnan(qv) if not np.isnan(qv) else True) and (
            (qv < 0.05) if not np.isnan(qv) else (pv < 0.05))
        is_nom     = (not np.isnan(pv)) and pv < 0.05 and not is_bh_sig
        is_any_sig = is_bh_sig or is_nom

        # Amber background for significant panels
        if is_any_sig:
            ax.set_facecolor('#fff8e1')

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            sns.violinplot(data=plot_df, x='immunotype_paper', y=feat,
                           order=IMMUNOTYPE_ORDER,
                           palette=IMMUNOTYPE_COLORS,
                           inner=None, linewidth=0.7, cut=0, ax=ax, alpha=0.3)

        sns.stripplot(data=plot_df, x='immunotype_paper', y=feat,
                      order=IMMUNOTYPE_ORDER,
                      palette=IMMUNOTYPE_COLORS,
                      size=3.5, alpha=0.85, jitter=True, ax=ax, zorder=3)

        # Median bars
        for xi, it in enumerate(IMMUNOTYPE_ORDER):
            vals = plot_df[plot_df['immunotype_paper'] == it][feat].values
            if len(vals):
                ax.plot([xi - 0.2, xi + 0.2], [np.median(vals)]*2,
                        color='black', lw=1.5, zorder=4)

        # KW p-value label (top-right corner)
        if not np.isnan(pv):
            if pv < 0.001:
                pv_str = 'KW: p<0.001'
            else:
                pv_str = f'KW: p={pv:.3f}'
            sig_color = '#c0392b' if is_any_sig else '0.50'
            ax.text(0.98, 0.97, pv_str, transform=ax.transAxes,
                    ha='right', va='top', fontsize=7,
                    fontweight='bold' if is_any_sig else 'normal',
                    fontstyle='italic',
                    color=sig_color)

        # Pairwise Dunn brackets for significant panels
        if is_any_sig and len(dunn_df) > 0:
            feat_dunn = dunn_df[dunn_df['feature'] == feat].copy()
            feat_dunn = feat_dunn[feat_dunn['pval'].notna()].copy()

            y_max   = plot_df[feat].max()
            y_min   = plot_df[feat].min()
            y_range = max(y_max - y_min, 1e-9)

            # Build map: immunotype label → x position
            it_pos = {it: i for i, it in enumerate(IMMUNOTYPE_ORDER)}

            # Sort pairs by x-distance so brackets don't overlap
            sig_pairs = feat_dunn[feat_dunn['pval'] < 0.05].sort_values('pval')
            y_offset = y_max + y_range * 0.10
            step     = y_range * 0.16
            for _, row in sig_pairs.iterrows():
                g1, g2 = row['g1'], row['g2']
                if g1 not in it_pos or g2 not in it_pos:
                    continue
                x1, x2 = sorted([it_pos[g1], it_pos[g2]])
                pv_d = row['pval']
                stars_d = '***' if pv_d < 0.001 else ('**' if pv_d < 0.01 else '*')
                ax.plot([x1, x1, x2, x2],
                        [y_offset, y_offset + y_range*0.04,
                         y_offset + y_range*0.04, y_offset],
                        lw=0.8, color='0.30')
                ax.text((x1 + x2)/2, y_offset + y_range*0.05, stars_d,
                        ha='center', va='bottom', fontsize=9,
                        fontweight='bold', color='#c0392b')
                y_offset += step

            # Extend y-axis to fit brackets
            ax.set_ylim(y_min - y_range*0.05, y_offset + y_range*0.08)

        ax.set_title(label, fontsize=9.5, fontweight='bold', pad=4)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(['E', 'C', 'D'], fontsize=9, fontweight='semibold')
        sns.despine(ax=ax)

        # Orange spines for significant panels
        if is_any_sig:
            for sp in ['left', 'bottom']:
                ax.spines[sp].set_edgecolor('#e67e22')
                ax.spines[sp].set_linewidth(1.8)

    # Hide unused
    for ax_idx in range(len(FEATURES_BY_IT), len(axes_flat)):
        axes_flat[ax_idx].set_visible(False)

    # Legend
    handles = [mpatches.Patch(facecolor=IMMUNOTYPE_COLORS[it], label=it)
               for it in IMMUNOTYPE_ORDER]
    fig.legend(handles=handles, loc='lower right',
               bbox_to_anchor=(1.0, 0.01), ncol=1, fontsize=8, frameon=False)

    n_by_it = df_valid['immunotype_paper'].value_counts()
    ns_str = '  '.join([f'{it[0]}:n={n_by_it.get(it,0)}' for it in IMMUNOTYPE_ORDER])
    fig.suptitle(
        f'IMC features by immunotype\n{ns_str}   |   Kruskal\u2013Wallis + BH · pairwise Dunn',
        fontsize=10, fontweight='bold', y=1.03, linespacing=1.6,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return fig


# ══════════════════════════════════════════════════════════════════════════
# Part C: Survival analysis
# ══════════════════════════════════════════════════════════════════════════

def run_km(df: pd.DataFrame,
            duration_col: str,
            event_col: str,
            group_col: str,
            groups: list,
            colors: dict,
            title: str,
            xlabel: str,
            out_path: str):
    """
    Kaplan-Meier curves with log-rank test.
    event_col: 1=event occurred (death/progression), 0=censored.
    """
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test, multivariate_logrank_test

    apply_style()
    fig, ax = plt.subplots(figsize=(4.5, 3.5))

    df_valid = df[df[group_col].isin(groups)].copy()
    df_valid = df_valid[df_valid[duration_col].notna() &
                         df_valid[event_col].notna()].copy()

    km_results = {}
    for grp in groups:
        sub = df_valid[df_valid[group_col] == grp]
        if len(sub) < 2:
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(sub[duration_col], event_observed=sub[event_col], label=grp)
        kmf.plot_survival_function(ax=ax, ci_show=True,
                                    color=colors[grp], linewidth=1.5)
        med = kmf.median_survival_time_
        km_results[grp] = {'n': len(sub), 'median': med,
                           'n_events': int(sub[event_col].sum())}

    # Log-rank test
    lr_result = None
    if len(groups) == 2:
        g1, g2 = groups
        s1 = df_valid[df_valid[group_col] == g1]
        s2 = df_valid[df_valid[group_col] == g2]
        if len(s1) >= 2 and len(s2) >= 2:
            lr_result = logrank_test(
                s1[duration_col], s2[duration_col],
                s1[event_col], s2[event_col])
            pval = lr_result.p_value
    else:
        # Multivariate log-rank
        try:
            mv_result = multivariate_logrank_test(
                df_valid[duration_col], df_valid[group_col], df_valid[event_col])
            pval = mv_result.p_value
        except Exception:
            pval = np.nan

    # Annotation
    pv_str = f'Log-rank p={pval:.3f}' if (pval is not None and not np.isnan(pval)) else ''
    ax.text(0.98, 0.97, pv_str, transform=ax.transAxes,
            ha='right', va='top', fontsize=7.5, style='italic')

    # At-risk table below
    at_risk_y = -0.08
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel('Survival probability', fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=7.5, frameon=False, loc='lower left')
    sns.despine(ax=ax)
    fig.tight_layout()

    save_fig(fig, out_path)
    print(f"  Saved: {out_path}")

    return pval, km_results


def run_all_survival(df: pd.DataFrame) -> pd.DataFrame:
    """Run KM for all stratifications."""
    # Censoring convention: censored_0yes=0 means event occurred, 1 means censored
    # Event = NOT censored = censored_0yes == 0
    df = df.copy()
    df['event_os']  = (df['censored_0yes'] == 0).astype(int)
    df['event_pfs'] = (df['censored_0yes'] == 0).astype(int)

    surv = df[df['survival_eligible'] == 1].copy()
    n_surv = len(surv)
    n_excl = len(df) - n_surv
    print(f"\n=== Part C: Survival Analysis ===")
    print(f"  Survival-eligible (CP<7): n={n_surv}  |  Excluded (CP≥7): n={n_excl}")

    it_colors = IMMUNOTYPE_COLORS
    resp_colors = {'Responder': PAL['Responder'],
                   'Non-Responder': PAL['Non-Responder']}

    surv['response_group'] = surv['response_binary'].map(
        {1: 'Responder', 0: 'Non-Responder'})

    surv_it = surv[surv['immunotype_paper'].isin(IMMUNOTYPE_ORDER)].copy()

    km_records = []

    # OS by immunotype
    pval, results = run_km(
        surv_it, 'os_months', 'event_os', 'immunotype_paper',
        IMMUNOTYPE_ORDER, it_colors,
        title=f'Overall Survival by Immunotype\n(n={len(surv_it)}, CP<7)',
        xlabel='Months',
        out_path=os.path.join(FIG_DIR, 'imc_km_os_immunotype'))
    for g, r in results.items():
        km_records.append({'endpoint': 'OS', 'stratification': 'immunotype',
                           'group': g, 'n': r['n'], 'n_events': r['n_events'],
                           'median_months': r['median'], 'logrank_pval': pval})

    # PFS by immunotype
    pval, results = run_km(
        surv_it, 'pfs_months', 'event_pfs', 'immunotype_paper',
        IMMUNOTYPE_ORDER, it_colors,
        title=f'Progression-Free Survival by Immunotype\n(n={len(surv_it)}, CP<7)',
        xlabel='Months',
        out_path=os.path.join(FIG_DIR, 'imc_km_pfs_immunotype'))
    for g, r in results.items():
        km_records.append({'endpoint': 'PFS', 'stratification': 'immunotype',
                           'group': g, 'n': r['n'], 'n_events': r['n_events'],
                           'median_months': r['median'], 'logrank_pval': pval})

    # OS by response
    pval, results = run_km(
        surv, 'os_months', 'event_os', 'response_group',
        ['Responder', 'Non-Responder'], resp_colors,
        title=f'Overall Survival by Response\n(n={len(surv)}, CP<7)',
        xlabel='Months',
        out_path=os.path.join(FIG_DIR, 'imc_km_os_response'))
    for g, r in results.items():
        km_records.append({'endpoint': 'OS', 'stratification': 'response',
                           'group': g, 'n': r['n'], 'n_events': r['n_events'],
                           'median_months': r['median'], 'logrank_pval': pval})

    # PFS by response
    pval, results = run_km(
        surv, 'pfs_months', 'event_pfs', 'response_group',
        ['Responder', 'Non-Responder'], resp_colors,
        title=f'Progression-Free Survival by Response\n(n={len(surv)}, CP<7)',
        xlabel='Months',
        out_path=os.path.join(FIG_DIR, 'imc_km_pfs_response'))
    for g, r in results.items():
        km_records.append({'endpoint': 'PFS', 'stratification': 'response',
                           'group': g, 'n': r['n'], 'n_events': r['n_events'],
                           'median_months': r['median'], 'logrank_pval': pval})

    km_df = pd.DataFrame(km_records)
    print("\nSurvival summary:")
    print(km_df[['endpoint', 'stratification', 'group', 'n',
                  'n_events', 'median_months', 'logrank_pval']].to_string(index=False))
    return km_df


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    apply_style()
    print(f"Loading {FEAT_PATH}")
    df = pd.read_csv(FEAT_PATH)
    print(f"  {len(df)} patients  |  columns: {len(df.columns)}")

    # ── Add survival_eligible flag ─────────────────────────────────────
    if 'child_pugh' in df.columns:
        df['survival_eligible'] = (df['child_pugh'] < 7).astype(int)
    else:
        # Merge from clinical CSV
        clin_path = os.path.join(
            PROJECT_DIR, 'data', 'spatial_data',
            'salie_bengsch_imc', 'ici_cohort', 'ICICohort_clinical.csv')
        clin = pd.read_csv(clin_path)[['patient_id', 'child_pugh']]
        df = df.merge(clin, on='patient_id', how='left')
        df['survival_eligible'] = (df['child_pugh'] < 7).astype(int)

    n_elig  = df['survival_eligible'].sum()
    n_excl  = (~df['survival_eligible'].astype(bool)).sum()
    cp_excl = df[df['survival_eligible'] == 0][['patient_id', 'child_pugh', 'recist']]
    print(f"\nSurvival eligible (CP<7): n={n_elig}")
    print(f"Excluded (CP≥7):           n={n_excl}")
    print("Excluded patients:")
    print(cp_excl.to_string(index=False))
    print(f"\n  NOTE: Paper reports n=5 excluded; our data has n={n_excl}.")
    print("  ICI_05(CP=12), ICI_07(CP=11), ICI_12(CP=10) may be PDF mis-extractions.")

    # Re-save updated features CSV with new columns
    df.to_csv(FEAT_PATH, index=False)
    print(f"\nUpdated patient features saved: {FEAT_PATH}")

    # ── Part A ─────────────────────────────────────────────────────────
    contingency_df, ct = run_contingency(df)
    fig_a = plot_immunotype_response_bar(ct)
    save_fig(fig_a, os.path.join(FIG_DIR, 'imc_immunotype_response_bar'))
    print(f"  Saved: {os.path.join(FIG_DIR, 'imc_immunotype_response_bar.{{png,pdf}}')}")
    contingency_df.to_csv(
        os.path.join(TABLE_DIR, 'imc_contingency_results.csv'), index=False)

    # ── Part B ─────────────────────────────────────────────────────────
    kw_df, dunn_df = run_kruskal_by_immunotype(df)
    fig_b = plot_features_by_immunotype(df, kw_df, dunn_df)
    save_fig(fig_b, os.path.join(FIG_DIR, 'imc_features_by_immunotype'))
    print(f"  Saved: {os.path.join(FIG_DIR, 'imc_features_by_immunotype.{{png,pdf}}')}")
    kw_df.to_csv(os.path.join(TABLE_DIR, 'imc_kw_results.csv'), index=False)
    if len(dunn_df):
        dunn_df.to_csv(os.path.join(TABLE_DIR, 'imc_dunn_results.csv'), index=False)

    # ── Part C ─────────────────────────────────────────────────────────
    km_df = run_all_survival(df)
    km_df.to_csv(os.path.join(TABLE_DIR, 'imc_km_results.csv'), index=False)

    print("\n=== Done ===")
    print(f"Output directory: {FIG_DIR}")


if __name__ == '__main__':
    main()
