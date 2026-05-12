#!/usr/bin/env python3
"""
05c_multilayer_predictor.py — Layer-wise + ensemble LOOCV predictor with
SHAP attribution for the multi-layer feature matrix.

Inputs
------
  out/05_multilayer/tables/patient_multilayer_features.csv
  out/05_multilayer/tables/nmf_programs_patient_activity.csv (optional add-on)

Design
------
Layers:
  L1  cnv_*        tumour-intrinsic (CNV burden + heterogeneity + clonality)
  L2  regulon_*    IRF1 vs MYC axis
  L3  frac_*       cell-type composition
  L4  sig_*        IO signatures
  L5  nmf_*        NMF program activity (appended if present)
Pipeline per configuration:
  - StandardScaler → L1 LogReg (ridge)  with LOOCV AUC
  - Ensemble (all layers concatenated) LASSO LOOCV AUC + ROC
  - Variance partition: drop-one-layer ΔAUC vs full
  - SHAP LinearExplainer on trained ensemble model for attribution

Outputs
-------
  out/05_multilayer/tables/layer_aucs.csv
  out/05_multilayer/tables/shap_layer_attribution.csv
  out/05_multilayer/tables/ensemble_coefficients.csv
  out/05_multilayer/figures/05c_roc_by_layer.png
  out/05_multilayer/figures/05c_variance_partition.png
  out/05_multilayer/figures/05c_shap_layer_bar.png
  out/05_multilayer/figures/05c_shap_beeswarm.png
"""
import warnings; warnings.filterwarnings("ignore")
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score, roc_curve
try:
    import shap
except ImportError:
    shap = None

from config.publication_style import apply_style, PAL
apply_style()

OUT = PROJECT / "out" / "05_multilayer"
TAB = OUT / "tables"; FIG = OUT / "figures"
for d in (FIG, TAB): d.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("[05c] Multi-layer LOOCV predictor + SHAP attribution")
print("=" * 70)

# ── Load features ───────────────────────────────────────────────────────
df = pd.read_csv(TAB / "patient_multilayer_features.csv")
print(f"\n[LOAD] feature matrix: {df.shape}")

# Optionally merge NMF
nmf_file = TAB / "nmf_programs_patient_activity.csv"
if nmf_file.exists():
    nmf = pd.read_csv(nmf_file)
    nmf_cols = [c for c in nmf.columns if c.startswith("NMF")]
    nmf = nmf.rename(columns={c: f"nmf_{c}" for c in nmf_cols})
    df = df.merge(nmf[["patient_id"] + [f"nmf_{c}" for c in nmf_cols]],
                   on="patient_id", how="left")
    print(f"   Merged {len(nmf_cols)} NMF program columns → {df.shape}")

if "dataset" not in df.columns:
    df["dataset"] = np.where(df["patient_id"].astype(str).str.startswith("pat"),
                             "Cappuyns2023", "GSE206325")
df = df[df["dataset"].astype(str) != "Cappuyns2023"].copy()
print(f"   Discovery-only cohort → {df.shape}")

# Response → 1 = Responder
y = (df["response"] == "Responder").astype(int).values
print(f"   Class balance: R={y.sum()}  NR={(1 - y).sum()}")

# ── Define layers by prefix ─────────────────────────────────────────────
layer_prefixes = {
    "L1_tumour_intrinsic": ("cnv_",),
    "L2_regulon":          ("regulon_",),
    "L3_composition":      ("frac_",),
    "L4_io_signatures":    ("sig_",),
    "L5_nmf_programs":     ("nmf_",),
}
layers = {}
for name, prefs in layer_prefixes.items():
    cols = [c for c in df.columns if any(c.startswith(p) for p in prefs)]
    # drop constant / all-NaN
    cols = [c for c in cols
             if df[c].notna().sum() >= len(df) * 0.5
             and df[c].nunique(dropna=True) > 1]
    if cols:
        layers[name] = cols
        print(f"   {name}: {len(cols)} features")

# ── Helpers ─────────────────────────────────────────────────────────────
def impute_scale(X):
    X = pd.DataFrame(X).fillna(pd.DataFrame(X).median())
    X = X.fillna(0)
    return StandardScaler().fit_transform(X.values)

def loocv_proba(X, y, C=1.0, penalty="l2"):
    loo = LeaveOneOut(); probas = np.zeros(len(y))
    solver = "liblinear" if penalty == "l1" else "lbfgs"
    for tr, te in loo.split(X):
        mdl = LogisticRegression(penalty=penalty, C=C, solver=solver,
                                  max_iter=2000, class_weight="balanced")
        mdl.fit(X[tr], y[tr])
        probas[te] = mdl.predict_proba(X[te])[0, 1]
    return probas

# ── Layer-wise AUCs ─────────────────────────────────────────────────────
print("\n[LAYER AUCs] LOOCV per layer …")
rows = []; roc_data = {}
for name, cols in layers.items():
    X = impute_scale(df[cols].values)
    if X.shape[1] > X.shape[0] - 2:
        # more features than obs: use L1
        proba = loocv_proba(X, y, C=0.5, penalty="l1")
    else:
        proba = loocv_proba(X, y, C=1.0, penalty="l2")
    auc = roc_auc_score(y, proba)
    fpr, tpr, _ = roc_curve(y, proba)
    rows.append(dict(layer=name, n_features=len(cols), auc=auc))
    roc_data[name] = (fpr, tpr, auc)
    print(f"   {name:25s}  n={len(cols):3d}  AUC={auc:.3f}")

# ── Ensemble: all features ─────────────────────────────────────────────
print("\n[ENSEMBLE] LOOCV LASSO on all features …")
all_cols = sum(layers.values(), [])
X_all = impute_scale(df[all_cols].values)
proba_all = loocv_proba(X_all, y, C=0.3, penalty="l1")
auc_all = roc_auc_score(y, proba_all)
fpr, tpr, _ = roc_curve(y, proba_all)
roc_data["Ensemble"] = (fpr, tpr, auc_all)
rows.append(dict(layer="Ensemble", n_features=X_all.shape[1], auc=auc_all))
print(f"   Ensemble AUC = {auc_all:.3f}")

pd.DataFrame(rows).to_csv(TAB / "layer_aucs.csv", index=False)

# Fit one ensemble on all data to extract coefficients + SHAP
final_mdl = LogisticRegression(penalty="l1", C=0.3, solver="liblinear",
                                 max_iter=3000, class_weight="balanced")
final_mdl.fit(X_all, y)
coef = pd.DataFrame({"feature": all_cols,
                      "coefficient": final_mdl.coef_.flatten()})
coef["layer"] = [next(L for L, cs in layers.items() if f in cs)
                  for f in coef["feature"]]
coef["abs_coef"] = coef["coefficient"].abs()
coef = coef.sort_values("abs_coef", ascending=False)
coef.to_csv(TAB / "ensemble_coefficients.csv", index=False)

# ── Variance partition: drop-one layer ──────────────────────────────────
print("\n[VARIANCE PARTITION] ΔAUC when each layer removed …")
vp = []
for name in layers:
    keep = [c for c in all_cols if c not in layers[name]]
    Xk = impute_scale(df[keep].values)
    p = loocv_proba(Xk, y, C=0.3, penalty="l1")
    a = roc_auc_score(y, p)
    vp.append(dict(layer=name, auc_without=a, delta=auc_all - a))
    print(f"   without {name:25s} AUC={a:.3f}  ΔAUC={auc_all - a:+.3f}")
vp_df = pd.DataFrame(vp).sort_values("delta", ascending=False)
vp_df.to_csv(TAB / "variance_partition.csv", index=False)

# ── SHAP attribution on ensemble ────────────────────────────────────────
print("\n[ATTRIBUTION] ensemble feature importance …")
attribution_source = "shap"
shap_values = None
if shap is not None:
    try:
        expl = shap.LinearExplainer(final_mdl, X_all)
        shap_values = expl.shap_values(X_all)
        if isinstance(shap_values, list):
            shap_values = shap_values[-1]
        shap_values = np.asarray(shap_values)
        per_feat = pd.DataFrame({"feature": all_cols,
                                 "mean_abs_shap": np.abs(shap_values).mean(axis=0)})
    except Exception as e:
        print(f"   [WARN] SHAP failed, falling back to |coefficient| attribution: {e}")
        attribution_source = "coefficient"
else:
    print("   [WARN] shap not installed, falling back to |coefficient| attribution")
    attribution_source = "coefficient"

if attribution_source == "coefficient":
    per_feat = coef[["feature", "abs_coef"]].rename(columns={"abs_coef": "mean_abs_shap"}).copy()

per_feat["layer"] = [next(L for L, cs in layers.items() if f in cs)
                      for f in per_feat["feature"]]
per_feat["attribution_source"] = attribution_source
per_feat = per_feat.sort_values("mean_abs_shap", ascending=False)
per_feat.to_csv(TAB / "shap_feature_attribution.csv", index=False)

per_layer = (per_feat.groupby("layer", as_index=False)["mean_abs_shap"]
             .sum().sort_values("mean_abs_shap", ascending=False))
per_layer["attribution_source"] = attribution_source
per_layer.to_csv(TAB / "shap_layer_attribution.csv", index=False)
print(per_layer.to_string(index=False))

# ─────────────────────────── Figures ────────────────────────────────────
# 1) ROC overlay
print("\n[FIG] ROC overlay …")
fig, ax = plt.subplots(figsize=(6.2, 5.8), constrained_layout=True)
palette = {
    "L1_tumour_intrinsic": "#882255",   # wine
    "L2_regulon":          "#DDCC77",   # sand
    "L3_composition":      "#44AA99",   # teal
    "L4_io_signatures":    "#332288",   # indigo
    "L5_nmf_programs":     "#AA4499",   # magenta
    "Ensemble":            "#000000",
}
for name, (fpr, tpr, auc) in roc_data.items():
    c = palette.get(name, "#666666")
    lw = 2.5 if name == "Ensemble" else 1.6
    ls = "-" if name == "Ensemble" else "--"
    ax.plot(fpr, tpr, color=c, linewidth=lw, linestyle=ls,
            label=f"{name}  AUC={auc:.2f}")
ax.plot([0, 1], [0, 1], color="grey", linewidth=0.8, linestyle=":")
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title("LOOCV ROC by layer and ensemble",
             fontsize=11, fontweight="bold")
ax.legend(loc="lower right", fontsize=8, frameon=False)
ax.set_aspect("equal")
for fmt in ("png", "pdf"):
    fig.savefig(FIG / f"05c_roc_by_layer.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)

# 2) Variance partition bar
print("[FIG] Variance partition …")
fig, ax = plt.subplots(figsize=(7, 3.5), constrained_layout=True)
colors = [palette.get(n, "#888") for n in vp_df["layer"]]
ax.barh(vp_df["layer"], vp_df["delta"], color=colors, edgecolor="black",
        linewidth=0.5, height=0.6)
ax.axvline(0, color="black", linewidth=0.8)
ax.set_xlabel("ΔAUC when layer removed  (larger = more unique contribution)")
ax.set_title("Variance partition across layers",
             fontsize=11, fontweight="bold")
for i, d in enumerate(vp_df["delta"]):
    ax.text(d + 0.002 * np.sign(d or 1), i, f"{d:+.3f}",
             va="center", ha="left" if d > 0 else "right", fontsize=8)
for fmt in ("png", "pdf"):
    fig.savefig(FIG / f"05c_variance_partition.{fmt}",
                dpi=300, bbox_inches="tight")
plt.close(fig)

# 3) Attribution layer bar
print("[FIG] Attribution layer bar …")
fig, ax = plt.subplots(figsize=(6.5, 3.5), constrained_layout=True)
colors = [palette.get(n, "#888") for n in per_layer["layer"]]
ax.barh(per_layer["layer"], per_layer["mean_abs_shap"], color=colors,
         edgecolor="black", linewidth=0.5, height=0.6)
xlab = "Σ mean |SHAP value|  (total attribution)" if attribution_source == "shap" \
    else "Σ |standardized coefficient|  (fallback attribution)"
ttl = "SHAP attribution per layer" if attribution_source == "shap" \
    else "Coefficient attribution per layer"
ax.set_xlabel(xlab)
ax.set_title(ttl, fontsize=11, fontweight="bold")
for fmt in ("png", "pdf"):
    fig.savefig(FIG / f"05c_shap_layer_bar.{fmt}",
                 dpi=300, bbox_inches="tight")
plt.close(fig)

# 4) SHAP beeswarm or coefficient fallback — top 15 features
print("[FIG] Feature attribution summary (top 15) …")
top_feat = per_feat.head(15).copy().iloc[::-1]
if attribution_source == "shap" and shap_values is not None:
    idx = [all_cols.index(f) for f in top_feat["feature"]]
    fig = plt.figure(figsize=(8, 6))
    shap.summary_plot(shap_values[:, idx], X_all[:, idx],
                       feature_names=top_feat["feature"].tolist(), show=False,
                       plot_type="dot", color_bar_label="Scaled feature value")
    plt.title("Top 15 features by mean |SHAP value|", fontsize=11,
               fontweight="bold")
    for fmt in ("png", "pdf"):
        plt.savefig(FIG / f"05c_shap_beeswarm.{fmt}",
                     dpi=300, bbox_inches="tight")
    plt.close("all")
else:
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    colors = [palette.get(layer, "#888") for layer in top_feat["layer"]]
    ax.barh(top_feat["feature"], top_feat["mean_abs_shap"], color=colors,
            edgecolor="black", linewidth=0.5, height=0.7)
    ax.set_xlabel("|standardized coefficient|")
    ax.set_title("Top 15 features by |coefficient|  (SHAP fallback)",
                 fontsize=11, fontweight="bold")
    for fmt in ("png", "pdf"):
        fig.savefig(FIG / f"05c_shap_beeswarm.{fmt}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)

print("\n[DONE] Multi-layer predictor ready.")
