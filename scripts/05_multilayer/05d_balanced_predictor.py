#!/usr/bin/env python3
"""
05d_balanced_predictor.py — Fair cross-layer comparison with top-k features
per layer.

Rationale
---------
When layers have unequal feature counts (L1=3 vs L3=44), an L1-penalised
ensemble tends to zero out low-dimensional layers that are, in fact,
informative (e.g. CNV). A fair comparison requires:

  (a) reduce each layer to the same k features (selected by univariate
      |t-statistic| vs response, using the full cohort) ;
  (b) re-run LOOCV per layer and on a balanced ensemble ;
  (c) drop-one-layer variance partition on balanced features.

This avoids the over-parameterised layer dominating the LASSO and yields an
interpretable per-layer ΔAUC decomposition.

Inputs:  out/05_multilayer/tables/patient_multilayer_features.csv
         out/05_multilayer/tables/nmf_programs_patient_activity.csv
Outputs: out/05_multilayer/tables/balanced_layer_aucs.csv
         out/05_multilayer/tables/balanced_variance_partition.csv
         out/05_multilayer/tables/balanced_selected_features.csv
         out/05_multilayer/figures/05d_balanced_roc.png
         out/05_multilayer/figures/05d_balanced_variance_partition.png
         out/05_multilayer/figures/05d_balanced_feature_effects.png
"""
import warnings; warnings.filterwarnings("ignore")
import sys
from pathlib import Path
PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))

import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score, roc_curve

from config.publication_style import apply_style, PAL
apply_style()

OUT = PROJECT / "out" / "05_multilayer"
TAB = OUT / "tables"; FIG = OUT / "figures"

# ── Load ────────────────────────────────────────────────────────────────
df = pd.read_csv(TAB / "patient_multilayer_features.csv")
nmf_file = TAB / "nmf_programs_patient_activity.csv"
if nmf_file.exists():
    nmf = pd.read_csv(nmf_file)
    nmf_cols = [c for c in nmf.columns if c.startswith("NMF")]
    nmf = nmf.rename(columns={c: f"nmf_{c}" for c in nmf_cols})
    df = df.merge(nmf[["patient_id"] + [f"nmf_{c}" for c in nmf_cols]],
                   on="patient_id", how="left")
if "dataset" not in df.columns:
    df["dataset"] = np.where(df["patient_id"].astype(str).str.startswith("pat"),
                             "Cappuyns2023", "GSE206325")
df = df[df["dataset"].astype(str) != "Cappuyns2023"].copy()
y = (df["response"] == "Responder").astype(int).values

layer_prefixes = {
    "L1_tumour_intrinsic": ("cnv_",),
    "L2_regulon":          ("regulon_",),
    "L3_composition":      ("frac_",),
    "L4_io_signatures":    ("sig_",),
    "L5_nmf_programs":     ("nmf_",),
}
if any(c.startswith("gf_cls_") for c in df.columns):
    layer_prefixes["L6_geneformer"] = ("gf_cls_",)
K = 3  # per-layer cap

def valid(col):
    v = df[col]
    return v.notna().sum() >= len(df) * 0.5 and v.nunique(dropna=True) > 1

IS_TECHNICAL = ("Unknown", "Doublet", "Empty", "doublet", "empty",
                "Proliferating_Proliferating", "Naive_Naive")

def is_technical(col):
    return any(tok in col for tok in IS_TECHNICAL)

def pick_top(cols, k):
    stats = []
    for c in cols:
        a = df.loc[y == 1, c].dropna()
        b = df.loc[y == 0, c].dropna()
        if len(a) < 2 or len(b) < 2:
            continue
        t, p = ttest_ind(a, b, equal_var=False)
        stats.append((c, abs(t), t, p))
    stats.sort(key=lambda r: r[1], reverse=True)
    return stats[:k]

print("[05d] Balanced (top-3 per layer) predictor")
selected = {}
sel_rows = []
for name, prefs in layer_prefixes.items():
    cols = [c for c in df.columns
             if any(c.startswith(p) for p in prefs) and valid(c)]
    if name == "L3_composition":
        cols = [c for c in cols if not is_technical(c)]
    top = pick_top(cols, K)
    if top:
        selected[name] = [c for c, *_ in top]
        for c, abst, t, p in top:
            sel_rows.append(dict(layer=name, feature=c,
                                   t_stat=t, abs_t=abst, p_value=p))
sel_df = pd.DataFrame(sel_rows)
sel_df["analysis_cohort"] = "GSE206325_discovery"
pd.DataFrame(sel_df).to_csv(TAB / "balanced_selected_features.csv", index=False)
for L, feats in selected.items():
    print(f"   {L}: {feats}")

def impute_scale(X):
    X = pd.DataFrame(X).fillna(pd.DataFrame(X).median()).fillna(0)
    return StandardScaler().fit_transform(X.values)

def loocv(X, y, C=1.0):
    loo = LeaveOneOut(); p = np.zeros(len(y))
    for tr, te in loo.split(X):
        m = LogisticRegression(penalty="l2", C=C, max_iter=3000,
                                class_weight="balanced")
        m.fit(X[tr], y[tr])
        p[te] = m.predict_proba(X[te])[0, 1]
    return p

# ── Per-layer AUCs ─────────────────────────────────────────────────────
rows = []; roc = {}
for name, feats in selected.items():
    X = impute_scale(df[feats].values)
    p = loocv(X, y)
    auc = roc_auc_score(y, p)
    fpr, tpr, _ = roc_curve(y, p)
    roc[name] = (fpr, tpr, auc)
    rows.append(dict(layer=name, n_features=len(feats), auc=auc))
    print(f"   {name:25s}  AUC={auc:.3f}")

# Balanced ensemble
all_feats = sum(selected.values(), [])
X_all = impute_scale(df[all_feats].values)
p_all = loocv(X_all, y, C=0.5)
auc_all = roc_auc_score(y, p_all)
fpr, tpr, _ = roc_curve(y, p_all)
roc["Balanced ensemble"] = (fpr, tpr, auc_all)
rows.append(dict(layer="Balanced ensemble", n_features=len(all_feats),
                   auc=auc_all))
print(f"   Balanced ensemble AUC = {auc_all:.3f}")
pd.DataFrame(rows).to_csv(TAB / "balanced_layer_aucs.csv", index=False)

# Variance partition
vp_rows = []
for name in selected:
    keep = [f for f in all_feats if f not in selected[name]]
    Xk = impute_scale(df[keep].values)
    pk = loocv(Xk, y, C=0.5)
    ak = roc_auc_score(y, pk)
    vp_rows.append(dict(layer=name, auc_without=ak, delta=auc_all - ak))
    print(f"   without {name:25s}  AUC={ak:.3f}  Δ={auc_all - ak:+.3f}")
vp_df = pd.DataFrame(vp_rows).sort_values("delta", ascending=False)
vp_df.to_csv(TAB / "balanced_variance_partition.csv", index=False)

# ── Figures ─────────────────────────────────────────────────────────────
palette = {
    "L1_tumour_intrinsic": "#882255",
    "L2_regulon":          "#DDCC77",
    "L3_composition":      "#44AA99",
    "L4_io_signatures":    "#332288",
    "L5_nmf_programs":     "#AA4499",
    "L6_geneformer":       "#117733",
    "Balanced ensemble":   "#000000",
}

fig, ax = plt.subplots(figsize=(6.2, 5.8), constrained_layout=True)
for name, (fpr, tpr, auc) in roc.items():
    c = palette.get(name, "#666")
    lw = 2.5 if name == "Balanced ensemble" else 1.6
    ls = "-" if name == "Balanced ensemble" else "--"
    ax.plot(fpr, tpr, color=c, linewidth=lw, linestyle=ls,
             label=f"{name}  AUC={auc:.2f}")
ax.plot([0, 1], [0, 1], color="grey", linewidth=0.8, linestyle=":")
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title(f"Balanced multi-layer LOOCV ROC (top-{K}/layer)",
              fontsize=11, fontweight="bold")
ax.legend(loc="lower right", fontsize=8, frameon=False)
ax.set_aspect("equal")
for fmt in ("png", "pdf"):
    fig.savefig(FIG / f"05d_balanced_roc.{fmt}",
                 dpi=300, bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 3.5), constrained_layout=True)
colors = [palette.get(n, "#888") for n in vp_df["layer"]]
ax.barh(vp_df["layer"], vp_df["delta"], color=colors, edgecolor="black",
        linewidth=0.5, height=0.6)
ax.axvline(0, color="black", linewidth=0.8)
ax.set_xlabel("ΔAUC when layer removed")
ax.set_title("Balanced variance partition",
              fontsize=11, fontweight="bold")
for i, d in enumerate(vp_df["delta"]):
    ax.text(d + 0.003 * np.sign(d or 1), i, f"{d:+.3f}",
             va="center", ha="left" if d > 0 else "right", fontsize=8)
for fmt in ("png", "pdf"):
    fig.savefig(FIG / f"05d_balanced_variance_partition.{fmt}",
                 dpi=300, bbox_inches="tight")
plt.close(fig)

# Feature effects forest (all top-3 features, |t|)
sel = pd.DataFrame(sel_rows)
sel["label"] = sel["feature"] + "  (" + sel["layer"].str.replace("_.*", "", regex=True) + ")"
sel = sel.sort_values("t_stat")
fig, ax = plt.subplots(figsize=(8, max(4, len(sel) * 0.3)),
                          constrained_layout=True)
colors = [PAL["Responder"] if t > 0 else PAL["Non_Responder"]
          for t in sel["t_stat"]]
ax.barh(sel["label"], sel["t_stat"], color=colors, edgecolor="black",
        linewidth=0.5, height=0.6)
ax.axvline(0, color="black", linewidth=0.8)
ax.set_xlabel("t-statistic (positive = higher in Responders)")
ax.set_title(f"Top-{K} features per layer (Welch's t-test)",
              fontsize=11, fontweight="bold")
for i, (_, r) in enumerate(sel.iterrows()):
    ax.text(r["t_stat"] + 0.1 * np.sign(r["t_stat"] or 1), i,
             f"p={r['p_value']:.2g}", va="center",
             ha="left" if r["t_stat"] > 0 else "right", fontsize=7)
for fmt in ("png", "pdf"):
    fig.savefig(FIG / f"05d_balanced_feature_effects.{fmt}",
                 dpi=300, bbox_inches="tight")
plt.close(fig)

print("[DONE]")
