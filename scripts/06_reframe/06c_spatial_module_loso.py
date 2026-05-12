"""06c: Spatial localization of locked immune-ecology modules in GSE238264.

Does NOT attempt responder-vs-non-responder classification from n=7 sections.
Instead asks: where in tissue do the CD8-dysfunction / Tfh / plasma / TLS
programs localize, and are section-level patterns robust to leave-one-section-out?

Also runs section-level (n=7) permutation for R-vs-NR with LOSO sensitivity
so that any residual R-vs-NR claim is explicitly section-level and honest.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

BASE = Path('/home/sadeghii/projects/spatial_immunotherapy_response')
DATA = BASE / 'data/spatial_data/GSE238264_HCC_IO_Visium'
OUT  = BASE / 'out/06_reframe/tables'
OUT.mkdir(parents=True, exist_ok=True)

# --- Gene lists (same as atlas signature scoring; subset to capture key programs) ---
MODULES = {
    'CD8_dysfunction':   ['PDCD1','LAG3','HAVCR2','TIGIT','CTLA4','TOX','ENTPD1','LAYN','CXCL13','GZMB','PRF1','IFNG'],
    'CD8_activation':    ['GZMA','GZMK','GZMH','NKG7','PRF1','IFNG','GNLY','CST7','CCL5'],
    'Tfh':               ['CXCR5','PDCD1','ICOS','BCL6','IL21','CXCL13','TOX2','MAF'],
    'Plasma':            ['MZB1','XBP1','JCHAIN','IGHG1','IGHG2','IGHG3','IGHA1','IGKC','CD79A','CD38'],
    'TLS':               ['CXCL13','CCL19','CCL21','CR2','MS4A1','CD19','LAMP3','FDCSP','CXCR5','LTB'],
    'Antigen_presentation':['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1','CD74','CIITA','B2M','HLA-A','HLA-B','HLA-C'],
    'Stromal_Hepatocyte':['ALB','HNF4A','CYP3A4','APOA1','APOB','TTR','TF','AFP'],
    'Immune_checkpoint': ['PDCD1','CD274','LAG3','TIGIT','CTLA4','HAVCR2','TIM3','IDO1'],
}

def score_section(ad, modules):
    # ensure log1p CPM-like baseline
    ad.var_names_make_unique()
    if 'log1p' not in ad.uns_keys() and ad.X.max() > 50:
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
    results = {}
    for name, genes in modules.items():
        present = [g for g in genes if g in ad.var_names]
        if len(present) < 3:
            results[name] = (np.nan, np.nan, 0)
            continue
        sc.tl.score_genes(ad, gene_list=present, score_name=f'M_{name}', random_state=0)
        s = ad.obs[f'M_{name}'].values
        results[name] = (float(np.mean(s)), float(np.std(s)), len(present))
    return results, ad.obs

sections = sorted([p.stem for p in DATA.glob('*.h5ad')])
print('sections:', sections)

# --- 1. Per-section module score summary ---
per_section_rows = []
spot_frames = []
for sec in sections:
    ad = sc.read_h5ad(DATA / f'{sec}.h5ad')
    res, obs = score_section(ad, MODULES)
    resp = 'Responder' if sec.endswith('R') and not sec.endswith('NR') else 'Non_Responder'
    for name, (mu, sd, ng) in res.items():
        per_section_rows.append(dict(section=sec, response=resp, module=name,
                                     mean=mu, sd=sd, n_genes_used=ng,
                                     n_spots=ad.n_obs))
    obs_sub = obs[[c for c in obs.columns if c.startswith('M_')]].copy()
    obs_sub['section'] = sec
    obs_sub['response'] = resp
    spot_frames.append(obs_sub)
per_section = pd.DataFrame(per_section_rows)
per_section.to_csv(OUT / 'spatial_modules_per_section.csv', index=False)
print('\nSection-level module means:')
print(per_section.pivot_table(index='section', columns='module', values='mean').round(3))

spots = pd.concat(spot_frames, axis=0)
spots.to_csv(OUT / 'spatial_modules_per_spot.csv.gz', index=False, compression='gzip')

# --- 2. Section-level R-vs-NR test (n=7) with LOSO sensitivity ---
ps = per_section.copy()
r_rows = []
for mod in MODULES:
    sub = ps[ps['module']==mod]
    r_vals  = sub.loc[sub['response']=='Responder',     'mean'].values
    nr_vals = sub.loc[sub['response']=='Non_Responder', 'mean'].values
    if len(r_vals) < 2 or len(nr_vals) < 2:
        r_rows.append(dict(module=mod, n_R=len(r_vals), n_NR=len(nr_vals),
                           mean_R=np.nan, mean_NR=np.nan, p_section=np.nan,
                           loso_p_min=np.nan, loso_p_max=np.nan,
                           loso_sign_stable=False, loso_fraction_p_lt_0p10=0.0))
        continue
    # exact permutation of section-level labels (n=7, C(7,4)=35 combos only)
    from itertools import combinations
    all_vals = np.concatenate([r_vals, nr_vals])
    obs_t = np.mean(r_vals) - np.mean(nr_vals)
    n_r = len(r_vals)
    perm_ts = []
    for idx in combinations(range(len(all_vals)), n_r):
        idx = list(idx)
        r = all_vals[idx]
        nr = np.delete(all_vals, idx)
        perm_ts.append(np.mean(r) - np.mean(nr))
    perm_ts = np.array(perm_ts)
    p_sec = float(np.mean(np.abs(perm_ts) >= abs(obs_t) - 1e-12))
    # LOSO — drop each section in turn, recompute p
    loso_ps = []; loso_signs = []
    for leave in sub['section'].unique():
        sub2 = sub[sub['section']!=leave]
        r2  = sub2.loc[sub2['response']=='Responder',    'mean'].values
        nr2 = sub2.loc[sub2['response']=='Non_Responder','mean'].values
        if len(r2) < 2 or len(nr2) < 2:
            continue
        all2 = np.concatenate([r2, nr2])
        obs2 = np.mean(r2) - np.mean(nr2)
        ps2 = []
        for idx in combinations(range(len(all2)), len(r2)):
            idx = list(idx)
            ps2.append(np.mean(all2[idx]) - np.mean(np.delete(all2, idx)))
        ps2 = np.array(ps2)
        loso_ps.append(np.mean(np.abs(ps2) >= abs(obs2) - 1e-12))
        loso_signs.append(np.sign(obs2))
    loso_ps = np.array(loso_ps)
    r_rows.append(dict(module=mod, n_R=len(r_vals), n_NR=len(nr_vals),
                       mean_R=float(np.mean(r_vals)), mean_NR=float(np.mean(nr_vals)),
                       p_section=p_sec,
                       loso_p_min=float(loso_ps.min()),
                       loso_p_max=float(loso_ps.max()),
                       loso_sign_stable=bool(len(set(loso_signs))==1),
                       loso_fraction_p_lt_0p10=float(np.mean(loso_ps<0.10))))
r_df = pd.DataFrame(r_rows)
# Claim class — deliberately conservative under the reframe.
# With n=7 sections (4R/3NR) and max LOSO p often >= 0.10, nothing is
# "supported". Borderline cases (p_section < 0.10 AND sign-stable under LOSO)
# are labelled exploratory_borderline; everything else is exploratory or null.
r_df['claim_class'] = r_df.apply(
    lambda r: 'exploratory_borderline' if np.isfinite(r['p_section']) and r['p_section']<0.10 and r['loso_sign_stable']
    else ('exploratory' if np.isfinite(r['p_section']) and r['p_section']<0.20
          else 'null'),
    axis=1)
r_df.to_csv(OUT / 'spatial_RvsNR_section_level_LOSO.csv', index=False)
print('\nSection-level R-vs-NR + LOSO (exact permutation, 4R/3NR):')
print(r_df.to_string(index=False))

# --- 3. Niche coherence: within-section spot-level correlations between modules ---
niche_rows = []
for sec, g in spots.groupby('section'):
    M = g[[f'M_{m}' for m in MODULES if f'M_{m}' in g.columns]].dropna(axis=1, how='all')
    if M.shape[1] < 2: continue
    C = M.corr(method='spearman')
    pairs = [('M_CD8_dysfunction','M_TLS'),
             ('M_CD8_dysfunction','M_Plasma'),
             ('M_Tfh','M_Plasma'),
             ('M_TLS','M_Plasma'),
             ('M_CD8_dysfunction','M_Antigen_presentation'),
             ('M_Stromal_Hepatocyte','M_CD8_dysfunction'),
             ('M_Stromal_Hepatocyte','M_Plasma')]
    for a,b in pairs:
        if a in C.index and b in C.columns:
            niche_rows.append(dict(section=sec, module_a=a, module_b=b,
                                   spearman_rho=float(C.loc[a,b])))
niche = pd.DataFrame(niche_rows)
niche.to_csv(OUT / 'spatial_niche_coherence.csv', index=False)
piv = niche.pivot_table(index='section', columns=['module_a','module_b'], values='spearman_rho').round(3)
print('\nNiche coherence (per-section spot-level Spearman rho):')
print(piv)
print('\n[summary] section-median niche rho:')
print(niche.groupby(['module_a','module_b'])['spearman_rho'].median().round(3))
