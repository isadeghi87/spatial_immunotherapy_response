#!/usr/bin/env python3
"""
05j_stratified_permutation_null.py — T2.1

Adds a design-stratified permutation null to the existing
revision_permutation_null.csv. The original null uses np.random.permutation(y),
which already preserves the marginal R/NR count. This script additionally
builds a *design-stratified* null that permutes labels only WITHIN strata
defined by n_cells-quartile. This tests whether depth imbalance inflates
the observed AUCs.

Inputs
------
out/05_multilayer/tables/patient_multilayer_features.csv
out/05_multilayer/tables/nmf_programs_patient_activity.csv
out/05_multilayer/tables/revision_layer_metrics.csv
out/05_multilayer/tables/revision_permutation_null.csv

Outputs
-------
out/05_multilayer/tables/revision_permutation_null_stratified.csv
out/05_multilayer/figures/05j_stratified_null.png
"""
from __future__ import annotations
import warnings; warnings.filterwarnings("ignore")
import sys
from pathlib import Path
PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score
from scipy.stats import mannwhitneyu

from config.publication_style import apply_style, PAL
apply_style()

RNG = np.random.default_rng(42)
N_PERM = 500
TAB = PROJECT / "out" / "05_multilayer" / "tables"
FIG = PROJECT / "out" / "05_multilayer" / "figures"

df = pd.read_csv(TAB / "patient_multilayer_features.csv")
nmf = pd.read_csv(TAB / "nmf_programs_patient_activity.csv")
nmf_cols = [c for c in nmf.columns if c.startswith("NMF")]
df = df.merge(nmf[["patient_id"] + nmf_cols].rename(
    columns={c: f"nmf_{c}" for c in nmf_cols}), on="patient_id", how="left")
df = df[df["dataset"] == "GSE206325"].reset_index(drop=True)
y = (df["response"] == "Responder").astype(int).values
n = len(y)
print(f"[05j] n={n} R={y.sum()} NR={n-y.sum()}")

# Strata: n_cells quartile  (technical depth)
q = pd.qcut(df["n_cells"], 3, labels=False, duplicates="drop").values
print(f"   strata (n_cells tertile) sizes: {np.bincount(q)}")
print(f"   R per stratum: {[int(y[q==k].sum()) for k in range(q.max()+1)]}")

def impute_scale(X):
    X = SimpleImputer(strategy="median").fit_transform(X)
    return StandardScaler().fit_transform(X)

def loocv_probs(X, y, C=1.0):
    p = np.zeros(len(y))
    for tr, te in LeaveOneOut().split(X):
        clf = LogisticRegression(C=C, penalty="l2", solver="liblinear",
                                 max_iter=2000, class_weight="balanced")
        clf.fit(X[tr], y[tr])
        p[te] = clf.predict_proba(X[te])[:, 1]
    return p

def pick_top(cols, y, k=3):
    scores = []
    for c in cols:
        x = df[c].values
        if np.all(np.isnan(x)) or np.nanstd(x) == 0: continue
        r = y.astype(float); r -= r.mean()
        xx = np.nan_to_num(x, nan=np.nanmedian(x)); xx -= xx.mean()
        denom = (np.linalg.norm(r) * np.linalg.norm(xx))
        scores.append((c, abs((r*xx).sum()/denom) if denom > 0 else 0))
    scores.sort(key=lambda t: -t[1])
    return [c for c, _ in scores[:k]]

LAYERS = {"L3_composition": ("frac_",),
          "Ensemble": ("cnv_", "regulon_", "frac_", "sig_", "nmf_NMF")}

def all_cols(prefs, exclude=("frac_Unknown", "frac_Doublet", "frac_Empty")):
    out = []
    for c in df.columns:
        if c in ("patient_id", "response", "dataset"): continue
        if any(c.startswith(p) for p in prefs) and not any(c.startswith(e) for e in exclude):
            if df[c].dtype.kind in "fi": out.append(c)
    return out

# Observed AUCs (re-compute to be self-consistent)
obs_metrics = {}
for name, prefs in LAYERS.items():
    feats = pick_top(all_cols(prefs), y, k=3 if name == "L3_composition" else 15)
    X = impute_scale(df[feats].values)
    obs_metrics[name] = roc_auc_score(y, loocv_probs(X, y, C=1.0 if name=="L3_composition" else 0.5))
print(f"   observed AUCs: {obs_metrics}")

# Stratified null
def stratified_permute(y, strata, rng):
    yperm = y.copy()
    for s in np.unique(strata):
        idx = np.where(strata == s)[0]
        yperm[idx] = rng.permutation(y[idx])
    return yperm

# Also keep an unstratified null for direct comparison
def unstratified_permute(y, rng):
    return rng.permutation(y)

print(f"[Null] {N_PERM} permutations x 2 schemes ...")
null_auc = {n: {"strat": [], "unstr": []} for n in LAYERS}
feat_cache = {n: all_cols(p) for n, p in LAYERS.items()}
for b in range(N_PERM):
    for scheme, fn in [("strat", lambda y: stratified_permute(y, q, RNG)),
                       ("unstr", lambda y: unstratified_permute(y, RNG))]:
        yp = fn(y)
        if yp.sum() == 0 or yp.sum() == n: continue
        for name in LAYERS:
            feats = pick_top(feat_cache[name], yp,
                             k=3 if name == "L3_composition" else 15)
            if len(feats) < 2: continue
            Xp = impute_scale(df[feats].values)
            C = 1.0 if name == "L3_composition" else 0.5
            null_auc[name][scheme].append(roc_auc_score(yp, loocv_probs(Xp, yp, C=C)))
    if (b+1) % 50 == 0:
        print(f"   perm {b+1}/{N_PERM}")

def summarize(null_arr, obs):
    a = np.asarray(null_arr)
    return dict(n_perm=len(a), median=float(np.median(a)),
                p95=float(np.percentile(a, 95)),
                p_one_sided=float((np.sum(a >= obs) + 1) / (len(a) + 1)),
                obs_auc=float(obs))

rows = []
for name in LAYERS:
    for scheme in ("strat", "unstr"):
        s = summarize(null_auc[name][scheme], obs_metrics[name])
        rows.append(dict(target=name, scheme=scheme, **s))
null_df = pd.DataFrame(rows)
null_df.to_csv(TAB / "revision_permutation_null_stratified.csv", index=False)
print("\nNull summary:")
print(null_df.to_string(index=False))

# ── Figure ─────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
for ax, name in zip(axes, LAYERS):
    a_un = np.asarray(null_auc[name]["unstr"])
    a_st = np.asarray(null_auc[name]["strat"])
    ax.hist(a_un, bins=30, alpha=0.55, color="#888888",
            edgecolor="black", linewidth=0.3, label=f"Unstratified (n={len(a_un)})")
    ax.hist(a_st, bins=30, alpha=0.55, color=PAL.get("Responder", "#0077BB"),
            edgecolor="black", linewidth=0.3,
            label=f"Stratified × n_cells tertile (n={len(a_st)})")
    ax.axvline(obs_metrics[name], color=PAL.get("Non_Responder", "#CC3311"),
               lw=2, label=f"Observed = {obs_metrics[name]:.3f}")
    s_st = null_df.query("target==@name and scheme=='strat'").iloc[0]
    s_un = null_df.query("target==@name and scheme=='unstr'").iloc[0]
    ax.set_title(f"{name}\np_unstr={s_un['p_one_sided']:.3f}  "
                 f"p_strat={s_st['p_one_sided']:.3f}",
                 fontsize=10, fontweight="bold")
    ax.set_xlabel("LOOCV AUC (null)", fontsize=9)
    ax.set_ylabel("Frequency", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
fig.suptitle("Design-stratified vs unstratified permutation null (GSE206325)",
             fontsize=11.5, fontweight="bold")
fig.savefig(FIG / "05j_stratified_null.png", dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(FIG / "05j_stratified_null.pdf", bbox_inches="tight", facecolor="white")
print(f"\nwrote {FIG/'05j_stratified_null.png'}")
