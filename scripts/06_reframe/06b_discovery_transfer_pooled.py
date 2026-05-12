"""06b: Immune-ecology reframe — discovery, external consistency, pooled sensitivity.

Locked module panel (platform-measurable on both cohorts):

  CD8 axis:
    - M_CD8_dysfunction       = mean(sig_T_cell_exhaustion, sig_Terminal_exhaustion)
    - M_CD8_activation        = mean(sig_T_cell_activation, sig_T_cell_cytotoxicity, sig_IFN_gamma)
    - M_CD8_progenitor        = sig_Progenitor_CD8
  Tfh / plasma / TLS axis:
    - M_Tfh                   = sig_Tfh_like
    - M_plasma_proxy          = NMF10 (plasma/Ig program)
    - M_TLS                   = mean(sig_Tfh_like, plasma proxy, sig_Antigen_presentation)
  Other immune context:
    - M_Treg                  = sig_Treg_signature
    - M_mregDC                = sig_mregDC
    - M_Macrophage_M1         = sig_M1_macrophage
    - M_Macrophage_M2         = sig_M2_macrophage
  Tumour / stromal:
    - M_Angiogenesis          = sig_Angiogenesis
    - M_TGFbeta               = sig_TGF_beta
    - M_WNT                   = sig_WNT_beta_catenin
    - M_Immune_checkpoint     = sig_Immune_checkpoint

Units are patient-level for both cohorts. Discovery = GSE206325 n=24;
external consistency = Cappuyns2023 n=32; Pooled = 56 with cohort as binary covariate.
All tests nonparametric + BH within the module family.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
import warnings
warnings.filterwarnings('ignore')

BASE = Path('/home/sadeghii/projects/spatial_immunotherapy_response')
OUT  = BASE / 'out/06_reframe/tables'
OUT.mkdir(parents=True, exist_ok=True)

feat = pd.read_csv(BASE / 'out/05_multilayer/tables/patient_multilayer_features.csv')
nmf  = pd.read_csv(BASE / 'out/05_multilayer/tables/nmf_programs_patient_activity.csv')[['patient_id','NMF10']]
df = feat.merge(nmf, on='patient_id', how='left')

# ----------------- Module construction ------------------
def m(cols):
    cols = [c for c in cols if c in df.columns]
    return df[cols].mean(axis=1)

df['M_CD8_dysfunction']   = m(['sig_T_cell_exhaustion','sig_Terminal_exhaustion'])
df['M_CD8_activation']    = m(['sig_T_cell_activation','sig_T_cell_cytotoxicity','sig_IFN_gamma'])
df['M_CD8_progenitor']    = df['sig_Progenitor_CD8']
df['M_Tfh']               = df['sig_Tfh_like']
df['M_plasma_proxy']      = df['NMF10']
df['M_TLS']               = m(['sig_Tfh_like','NMF10','sig_Antigen_presentation'])
df['M_Treg']              = df['sig_Treg_signature']
df['M_mregDC']            = df['sig_mregDC']
df['M_Macrophage_M1']     = df['sig_M1_macrophage']
df['M_Macrophage_M2']     = df['sig_M2_macrophage']
df['M_Angiogenesis']      = df['sig_Angiogenesis']
df['M_TGFbeta']           = df['sig_TGF_beta']
df['M_WNT']               = df['sig_WNT_beta_catenin']
df['M_Immune_checkpoint'] = df['sig_Immune_checkpoint']

MODULES = [c for c in df.columns if c.startswith('M_')]

# sanity print
print('Locked modules:', MODULES)
print('N total:', len(df))
for ds,g in df.groupby('dataset'):
    print(f"  {ds} n={len(g)} (R={(g['response']=='Responder').sum()}, NR={(g['response']=='Non-Responder').sum()})")

# ----------------- Helper ------------------
def mw_effect(x_r, x_nr):
    """Mann-Whitney U + cliff's-delta-style effect (common-language)."""
    x_r = np.asarray(x_r, dtype=float)
    x_nr = np.asarray(x_nr, dtype=float)
    if len(x_r) < 2 or len(x_nr) < 2:
        return np.nan, np.nan, np.nan, np.nan
    u, p = stats.mannwhitneyu(x_r, x_nr, alternative='two-sided')
    n1, n2 = len(x_r), len(x_nr)
    # rank-biserial correlation / common-language effect
    r_b = 2 * u / (n1 * n2) - 1
    md = float(np.median(x_r) - np.median(x_nr))
    return float(u), float(p), float(r_b), md

# ----------------- (1) Discovery GSE206325 ------------------
disc = df[df['dataset']=='GSE206325'].copy()
disc['y'] = (disc['response']=='Responder').astype(int)
rows = []
for mod in MODULES:
    u, p, rb, md = mw_effect(disc.loc[disc['y']==1, mod], disc.loc[disc['y']==0, mod])
    rows.append(dict(module=mod, n_R=int(disc['y'].sum()), n_NR=int((disc['y']==0).sum()),
                     median_R=float(np.median(disc.loc[disc['y']==1, mod])),
                     median_NR=float(np.median(disc.loc[disc['y']==0, mod])),
                     median_diff=md, rank_biserial=rb, U=u, p_raw=p))
disc_res = pd.DataFrame(rows)
disc_res['p_BH'] = multipletests(disc_res['p_raw'], method='fdr_bh')[1]
# NOTE: no auxiliary class column is written here. The canonical
# strength_class for each module lives in results_index.yaml under the
# reframe_2026_04_21 block and must not be duplicated in table-level
# labels (schema compliance: allowed values only {strong, supported_limited,
# exploratory, null, inconclusive, likely_artifact}).
disc_res = disc_res.sort_values('p_BH')
disc_res.to_csv(OUT / 'discovery_module_associations.csv', index=False)
print('\n[1] DISCOVERY (GSE206325, n=24) patient-level module associations:')
print(disc_res.to_string(index=False))

# Redundancy diagnostics
corr = disc[MODULES].corr(method='spearman')
corr.to_csv(OUT / 'discovery_module_correlation.csv')
# flag pairs with |rho| > 0.7
redun = []
for i,a in enumerate(MODULES):
    for b in MODULES[i+1:]:
        rho = corr.loc[a,b]
        if abs(rho) >= 0.7:
            redun.append(dict(module_a=a, module_b=b, spearman_rho=float(rho)))
redun_df = pd.DataFrame(redun).sort_values('spearman_rho', key=lambda s: s.abs(), ascending=False) if redun else pd.DataFrame(columns=['module_a','module_b','spearman_rho'])
redun_df.to_csv(OUT / 'discovery_module_redundancy.csv', index=False)
print('\n[1b] REDUNDANCY pairs |rho|>=0.7:')
print(redun_df.to_string(index=False))

# -------- Top-5 BH-strong axis: explicit per-pair rho and PC1 score --------
# The audit flagged that the summary "top-5 rho 0.83-0.99" overstates the
# collinearity because M_plasma_proxy is less tightly coupled. Here we
# compute (i) the exact rho matrix restricted to the 5 BH-strong modules,
# (ii) a "core-4" (CD8_dysfunction, Tfh, TLS, Immune_checkpoint) block,
# and (iii) a patient-level immune-axis score = PC1 of the standardized
# top-5 modules, so downstream figures can cite an audited output.
top5 = ['M_plasma_proxy', 'M_Tfh', 'M_TLS', 'M_CD8_dysfunction', 'M_Immune_checkpoint']
core4 = ['M_Tfh', 'M_TLS', 'M_CD8_dysfunction', 'M_Immune_checkpoint']
top5_rho = disc[top5].corr(method='spearman')
top5_rho.to_csv(OUT / 'discovery_top5_rho_matrix.csv')
pairs_top5 = []
for i, a in enumerate(top5):
    for b in top5[i+1:]:
        pairs_top5.append(dict(module_a=a, module_b=b,
                               spearman_rho=float(top5_rho.loc[a, b]),
                               block='core4' if (a in core4 and b in core4) else 'plasma_vs_core4'))
pairs_top5_df = pd.DataFrame(pairs_top5).sort_values('spearman_rho', ascending=False)
pairs_top5_df.to_csv(OUT / 'discovery_top5_axis_pairs.csv', index=False)

# Audited summary: report true rho ranges (no rounding narrative).
core_vals = pairs_top5_df.loc[pairs_top5_df['block'] == 'core4', 'spearman_rho']
plasma_vals = pairs_top5_df.loc[pairs_top5_df['block'] == 'plasma_vs_core4', 'spearman_rho']
all_vals = pairs_top5_df['spearman_rho']
axis_summary = pd.DataFrame([
    dict(scope='top5_all', n_pairs=len(all_vals), rho_min=float(all_vals.min()),
         rho_median=float(all_vals.median()), rho_max=float(all_vals.max())),
    dict(scope='core4_only', n_pairs=len(core_vals), rho_min=float(core_vals.min()),
         rho_median=float(core_vals.median()), rho_max=float(core_vals.max())),
    dict(scope='plasma_vs_core4', n_pairs=len(plasma_vals), rho_min=float(plasma_vals.min()),
         rho_median=float(plasma_vals.median()), rho_max=float(plasma_vals.max())),
])
axis_summary.to_csv(OUT / 'discovery_top5_axis_summary.csv', index=False)
print('\n[1c] TOP-5 BH-strong axis rho summary (true ranges):')
print(axis_summary.to_string(index=False))

# Immune-axis PC1 score (patient-level, z-scored inputs, sign-aligned to CD8_dysfunction)
from sklearn.decomposition import PCA
Z = (disc[top5] - disc[top5].mean()) / disc[top5].std(ddof=0)
pca = PCA(n_components=1, random_state=0).fit(Z.values)
pc1 = pca.transform(Z.values).ravel()
# sign-align so a higher PC1 points toward "more CD8_dysfunctional"
if np.corrcoef(pc1, Z['M_CD8_dysfunction'].values)[0, 1] < 0:
    pc1 = -pc1
    loadings = -pca.components_[0]
else:
    loadings = pca.components_[0]
axis_df = pd.DataFrame({
    'patient_id': disc['patient_id'].astype(str).values,
    'response': disc['response'].values,
    'immune_axis_PC1': pc1,
})
axis_df.to_csv(OUT / 'discovery_immune_axis_pc1.csv', index=False)
load_df = pd.DataFrame({
    'module': top5,
    'pc1_loading': loadings,
    'variance_explained_pc1': [float(pca.explained_variance_ratio_[0])] * len(top5),
})
load_df.to_csv(OUT / 'discovery_immune_axis_pc1_loadings.csv', index=False)
print(f"\n[1d] Immune-axis PC1 variance explained = {pca.explained_variance_ratio_[0]:.3f}")

# ----------------- (2) External consistency: Cappuyns2023 ------------------
# n=32 filter: exclude DeathBeforeImaging patients (pat5, pat28, pat29, pat38).
# These 4 patients are not treatment-response-evaluable and inflate NR 13→17.
_DEATH_BEFORE_IMAGING = {'pat5', 'pat28', 'pat29', 'pat38'}
tra = df[df['dataset']=='Cappuyns2023'].copy()
tra = tra[~tra['patient_id'].astype(str).isin(_DEATH_BEFORE_IMAGING)].copy()
tra['y'] = (tra['response']=='Responder').astype(int)
print(f'[2] Cappuyns external consistency: n={len(tra)} patients '
      f'(R={tra["y"].sum()}, NR={(tra["y"]==0).sum()})')
rows = []
for mod in MODULES:
    disc_range = disc[mod].max() - disc[mod].min()
    tra_range  = tra[mod].max()  - tra[mod].min()
    dr_frac    = tra_range / disc_range if disc_range>0 else np.nan
    # overlap via shifted min/max
    disc_q5, disc_q95 = np.percentile(disc[mod], [5,95])
    tra_q5,  tra_q95  = np.percentile(tra[mod],  [5,95])
    overlap_lo = max(disc_q5, tra_q5); overlap_hi = min(disc_q95, tra_q95)
    overlap_frac = max(0.0, overlap_hi - overlap_lo) / max(disc_q95-disc_q5, 1e-9)
    u, p, rb, md = mw_effect(tra.loc[tra['y']==1, mod], tra.loc[tra['y']==0, mod])
    rows.append(dict(module=mod,
                     n_R=int(tra['y'].sum()), n_NR=int((tra['y']==0).sum()),
                     discovery_range=float(disc_range),
                     transfer_range=float(tra_range),
                     transfer_range_fraction=float(dr_frac) if np.isfinite(dr_frac) else np.nan,
                     dynamic_range_adequate=bool(dr_frac>=0.4) if np.isfinite(dr_frac) else False,
                     q5_q95_overlap_frac=float(overlap_frac),
                     median_R=float(np.median(tra.loc[tra['y']==1, mod])),
                     median_NR=float(np.median(tra.loc[tra['y']==0, mod])),
                     median_diff=md, rank_biserial=rb, U=u, p_raw=p))
tra_res = pd.DataFrame(rows)
tra_res['p_BH'] = multipletests(tra_res['p_raw'], method='fdr_bh')[1]
# concordance with discovery sign
tra_res = tra_res.merge(disc_res[['module','rank_biserial']].rename(columns={'rank_biserial':'rb_discovery'}), on='module')
tra_res['sign_concordant'] = np.sign(tra_res['rank_biserial']) == np.sign(tra_res['rb_discovery'])
def transfer_class(r):
    if not r['dynamic_range_adequate']:
        return 'non_transferable_dynamic_range'
    if r['p_BH'] < 0.10 and r['sign_concordant']:
        return 'transferable_concordant'
    if abs(r['rank_biserial']) >= 0.3 and r['sign_concordant']:
        return 'direction_consistent_underpowered'
    if r['sign_concordant']:
        return 'null_concordant'
    return 'null_or_opposite'
tra_res['transfer_class'] = tra_res.apply(transfer_class, axis=1)
tra_res['external_consistency_class'] = tra_res['transfer_class'].replace({
    'transferable_concordant': 'concordant_attenuated',
    'direction_consistent_underpowered': 'direction_consistent_limited',
    'non_transferable_dynamic_range': 'limited_dynamic_range',
    'null_concordant': 'null_but_direction_concordant',
    'null_or_opposite': 'null_or_opposite_direction',
})
tra_res = tra_res.sort_values(['transfer_class','p_BH'])
tra_res.to_csv(OUT / 'cappuyns_locked_external_consistency.csv', index=False)
# Compatibility alias retained for old figure/table readers.
tra_res.to_csv(OUT / 'cappuyns_locked_transfer.csv', index=False)
print('\n[2] EXTERNAL CONSISTENCY (Cappuyns2023, n=32) locked-module results:')
print(tra_res[['module','rank_biserial','rb_discovery','sign_concordant',
               'transfer_range_fraction','p_raw','p_BH','external_consistency_class']].to_string(index=False))

# ----------------- (3) Pooled sensitivity with cohort term ------------------
# Apply n=32 filter to Cappuyns rows before pooling
pool = df.copy()
_cap_mask = pool['dataset'] == 'Cappuyns2023'
_dead_mask = pool['patient_id'].astype(str).isin(_DEATH_BEFORE_IMAGING)
pool = pool[~(_cap_mask & _dead_mask)].copy()
pool['y']       = (pool['response']=='Responder').astype(int)
pool['cohort']  = (pool['dataset']=='Cappuyns2023').astype(int)  # 1=Cappuyns, 0=GSE206325
rows = []
for mod in MODULES:
    # standardise within-module across pooled sample for numerical stability
    x = pool[mod].values
    z = (x - np.mean(x)) / (np.std(x) + 1e-12)
    X = sm.add_constant(np.column_stack([z, pool['cohort'].values]))
    try:
        fit = sm.Logit(pool['y'], X).fit(disp=False, maxiter=200)
        beta_mod, p_mod     = float(fit.params[1]), float(fit.pvalues[1])
        beta_cohort, p_coh  = float(fit.params[2]), float(fit.pvalues[2])
        # cohort-only model to quantify incremental gain
        X0 = sm.add_constant(pool['cohort'].values)
        fit0 = sm.Logit(pool['y'], X0).fit(disp=False, maxiter=200)
        aic_full, aic_null = float(fit.aic), float(fit0.aic)
        delta_aic = aic_full - aic_null
        rows.append(dict(module=mod, beta_module_std=beta_mod, p_module=p_mod,
                         beta_cohort=beta_cohort, p_cohort=p_coh,
                         aic_full=aic_full, aic_cohort_only=aic_null, delta_aic=delta_aic))
    except Exception as e:
        rows.append(dict(module=mod, beta_module_std=np.nan, p_module=np.nan,
                         beta_cohort=np.nan, p_cohort=np.nan,
                         aic_full=np.nan, aic_cohort_only=np.nan, delta_aic=np.nan))
pooled = pd.DataFrame(rows)
pooled['p_BH'] = multipletests(pooled['p_module'].fillna(1.0), method='fdr_bh')[1]
def pooled_class(r):
    if not np.isfinite(r['p_module']): return 'fit_failed'
    if r['p_BH']>=0.10 and r['p_cohort']<0.05:
        return 'cohort_dominated'
    if r['p_BH']<0.10 and r['p_cohort']<0.05 and abs(r['beta_cohort'])>abs(r['beta_module_std']):
        return 'cohort_larger_than_module'
    if r['p_BH']<0.10:
        return 'module_significant'
    return 'null'
pooled['pooled_class'] = pooled.apply(pooled_class, axis=1)
pooled = pooled.sort_values(['pooled_class','p_BH'])
pooled.to_csv(OUT / 'pooled_sensitivity_cohort_adjusted.csv', index=False)
print('\n[3] POOLED LOGIT with cohort term (n=56, R=26, NR=30):')
print(pooled.to_string(index=False))

# ----------------- (3b) Pooled interaction sensitivity (module × cohort) ------
# Additive pooling can hide effect-sign reversal between cohorts. The transfer
# table already flags sign flips for M_Tfh, M_CD8_dysfunction, M_Immune_checkpoint
# and sign-concordant-but-attenuated M_TLS, M_plasma_proxy. This block fits
#   y ~ module_z + cohort + module_z:cohort
# per module so the reframe has a FORMAL interaction test showing the additive
# pooled model is not a safe effect estimator.
rows_int = []
for mod in MODULES:
    x = pool[mod].values
    z = (x - np.mean(x)) / (np.std(x) + 1e-12)
    c = pool['cohort'].values.astype(float)
    zc = z * c
    X_int = sm.add_constant(np.column_stack([z, c, zc]))
    try:
        fit_int = sm.Logit(pool['y'], X_int).fit(disp=False, maxiter=300)
        # LR test: interaction vs additive
        X_add = sm.add_constant(np.column_stack([z, c]))
        fit_add = sm.Logit(pool['y'], X_add).fit(disp=False, maxiter=300)
        lr = 2 * (fit_int.llf - fit_add.llf)
        p_lr = float(stats.chi2.sf(lr, df=1))
        rows_int.append(dict(
            module=mod,
            beta_module_GSE206325=float(fit_int.params[1]),              # cohort=0 baseline (GSE206325)
            beta_module_Cappuyns_sum=float(fit_int.params[1] + fit_int.params[3]),
            beta_interaction=float(fit_int.params[3]),
            p_interaction_wald=float(fit_int.pvalues[3]),
            p_interaction_LRT=p_lr,
            sign_reversal_between_cohorts=bool(
                np.sign(fit_int.params[1]) * np.sign(fit_int.params[1] + fit_int.params[3]) < 0
            ),
            aic_interaction=float(fit_int.aic),
            aic_additive=float(fit_add.aic),
            delta_aic_interaction_vs_additive=float(fit_int.aic - fit_add.aic),
        ))
    except Exception:
        rows_int.append(dict(
            module=mod, beta_module_GSE206325=np.nan,
            beta_module_Cappuyns_sum=np.nan, beta_interaction=np.nan,
            p_interaction_wald=np.nan, p_interaction_LRT=np.nan,
            sign_reversal_between_cohorts=False,
            aic_interaction=np.nan, aic_additive=np.nan,
            delta_aic_interaction_vs_additive=np.nan,
        ))
pooled_int = pd.DataFrame(rows_int)
pooled_int['p_interaction_BH'] = multipletests(
    pooled_int['p_interaction_LRT'].fillna(1.0), method='fdr_bh')[1]
pooled_int = pooled_int.sort_values('p_interaction_LRT')
pooled_int.to_csv(OUT / 'pooled_sensitivity_module_x_cohort.csv', index=False)
print('\n[3b] POOLED MODULE x COHORT INTERACTION (LRT):')
print(pooled_int[['module', 'beta_module_GSE206325', 'beta_module_Cappuyns_sum',
                  'sign_reversal_between_cohorts', 'p_interaction_LRT',
                  'p_interaction_BH']].to_string(index=False))

# ----------------- (4) NMF10 → plasma/Ig re-framing ------------------
# already have nmf10_partial_correlation_pooled; re-examine explicitly after
# adjusting for **total plasma activity** (frac_B_Plasma where available + NMF10-shared genes)
# Here we use frac_Plasma (coarse harmonised) as the adjusting variable.
_plasma_col = 'frac_Plasma' if 'frac_Plasma' in disc.columns else 'frac_B_Plasma'
disc['plasma_frac'] = disc[_plasma_col]
# discovery partial correlation
from scipy.stats import spearmanr
rows = []
_pool_plasma = df[_plasma_col].fillna(0) if _plasma_col in df.columns else pd.Series(0, index=df.index)
for comparison, g in [('GSE206325', disc),
                      ('Pooled', df.assign(plasma_frac=_pool_plasma))]:
    y = (g['response']=='Responder').astype(int).values
    nmf10 = g['NMF10'].values
    plas  = g['plasma_frac'].fillna(0).values
    # naive assoc
    rho1, p1 = spearmanr(nmf10, y)
    # residualise
    try:
        b = np.linalg.lstsq(np.column_stack([np.ones_like(plas), plas]), nmf10, rcond=None)[0]
        resid = nmf10 - (b[0] + b[1]*plas)
        rho2, p2 = spearmanr(resid, y)
    except Exception:
        rho2, p2 = np.nan, np.nan
    # logit with plasma control
    try:
        X = sm.add_constant(np.column_stack([nmf10, plas]))
        fit = sm.Logit(y, X).fit(disp=False, maxiter=200)
        beta_nmf, p_nmf = float(fit.params[1]), float(fit.pvalues[1])
        beta_plas, p_plas = float(fit.params[2]), float(fit.pvalues[2])
    except Exception:
        beta_nmf=p_nmf=beta_plas=p_plas=np.nan
    rows.append(dict(cohort=comparison, n=len(g),
                     spearman_nmf10_vs_y=float(rho1), p_naive=float(p1),
                     spearman_resid_nmf10_vs_y=float(rho2) if np.isfinite(rho2) else np.nan,
                     p_resid=float(p2) if np.isfinite(p2) else np.nan,
                     logit_beta_NMF10=beta_nmf, p_logit_NMF10=p_nmf,
                     logit_beta_plasma=beta_plas, p_logit_plasma=p_plas))
nmf10_refr = pd.DataFrame(rows)
nmf10_refr.to_csv(OUT / 'nmf10_plasma_reframe.csv', index=False)
print('\n[4] NMF10 PLASMA-ADJUSTED re-examination:')
print(nmf10_refr.to_string(index=False))

# ----------------- (5) Audited module-level external-consistency matrix ------
consistency = (
    disc_res[['module', 'rank_biserial', 'median_diff', 'p_BH']]
    .rename(columns={
        'rank_biserial': 'discovery_rank_biserial',
        'median_diff': 'discovery_median_diff',
        'p_BH': 'discovery_p_BH',
    })
    .merge(
        tra_res[['module', 'rank_biserial', 'median_diff', 'p_BH',
                 'transfer_range_fraction', 'q5_q95_overlap_frac',
                 'sign_concordant', 'external_consistency_class']],
        on='module',
        how='left',
    )
    .rename(columns={
        'rank_biserial': 'cappuyns_rank_biserial',
        'median_diff': 'cappuyns_median_diff',
        'p_BH': 'cappuyns_p_BH',
    })
)
consistency['effect_delta_rank_biserial'] = (
    consistency['cappuyns_rank_biserial'] - consistency['discovery_rank_biserial']
)
consistency['manuscript_routing'] = np.where(
    consistency['module'].isin(['M_CD8_dysfunction', 'M_Tfh', 'M_TLS', 'M_Immune_checkpoint']),
    'main_discovery_axis',
    np.where(consistency['module'].eq('M_plasma_proxy'), 'main_with_plasma_caveat', 'supplement')
)
consistency['interpretation'] = np.where(
    consistency['cappuyns_p_BH'] < 0.10,
    'external_consistency_supported',
    np.where(
        consistency['sign_concordant'],
        'directional_or_null_external_context',
        'cohort_or_regimen_boundary'
    )
)
consistency.to_csv(OUT / 'external_consistency_module_summary.csv', index=False)
print('\n[5] Module-level external consistency summary written:')
print(consistency[['module', 'discovery_rank_biserial', 'cappuyns_rank_biserial',
                   'cappuyns_p_BH', 'external_consistency_class',
                   'manuscript_routing']].to_string(index=False))

print('\nAll reframe tables written to', OUT)
