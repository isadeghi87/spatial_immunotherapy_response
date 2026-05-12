#!/usr/bin/env python3
"""
03c_cappuyns_validation.py — External validation of L3 composition predictor
on the Cappuyns et al. 2023 (Nat Commun) atezo+bev cohort.

Design
------
Discovery cohort : GSE206325 (Magen 2023) — 24 anti-PD-1 patients (7 R / 17 NR)
Validation cohort: Cappuyns2023 — 32 atezo+bev patients (19 R / 13 NR)

The L3 composition model is trained on all GSE206325 patients (no inner CV
needed — the LOOCV was already reported in 05d). The full training-set model
is then applied to predict response probability for every Cappuyns patient.

Metrics reported:
  - Validation ROC-AUC with 2,000-iteration patient-level bootstrap 95% CI
  - Validation PR-AUC (baseline = Cappuyns prevalence)
  - Validation balanced accuracy (optimal threshold from training set)
  - Per-layer validation AUCs (L1–L5)
  - L3 composition waterfall plot
  - L3 ROC overlay: discovery vs validation

Inputs
------
  out/05_multilayer/tables/patient_multilayer_features.csv
  out/05_multilayer/tables/nmf_programs_patient_activity.csv
  out/05_multilayer/tables/balanced_selected_features.csv

Outputs
-------
  out/03_validation/tables/cappuyns_validation_results.csv
  out/03_validation/tables/cappuyns_per_patient_predictions.csv
  out/03_validation/figures/03c_validation_roc.png
  out/03_validation/figures/03c_validation_waterfall.png
  out/03_validation/figures/03c_validation_perlayer.png
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
from scipy.stats import mannwhitneyu
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (roc_auc_score, roc_curve,
                             precision_recall_curve, average_precision_score,
                             balanced_accuracy_score)

from config.publication_style import apply_style, PAL
apply_style()

OUT = PROJECT / "out" / "03_validation"
FIG = OUT / "figures"; TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

ML_TAB = PROJECT / "out" / "05_multilayer" / "tables"

# ═══════════════════════════════════════════════════════════════════════
# 1.  LOAD FEATURES
# ═══════════════════════════════════════════════════════════════════════

print("=" * 70)
print("[03c] External validation: Cappuyns 2023 (atezo+bev)")
print("=" * 70)

df = pd.read_csv(ML_TAB / "patient_multilayer_features.csv")

# Merge NMF features if available
nmf_file = ML_TAB / "nmf_programs_patient_activity.csv"
if nmf_file.exists():
    nmf = pd.read_csv(nmf_file)
    nmf_cols = [c for c in nmf.columns if c.startswith("NMF")]
    nmf = nmf.rename(columns={c: f"nmf_{c}" for c in nmf_cols})
    df = df.merge(nmf[["patient_id"] + [f"nmf_{c}" for c in nmf_cols]],
                  on="patient_id", how="left")

if "dataset" not in df.columns:
    # Fallback: Cappuyns patient IDs start with "pat"
    df["dataset"] = np.where(df["patient_id"].str.startswith("pat"),
                             "Cappuyns2023", "GSE206325")

# ── Split discovery / validation ────────────────────────────────────────
disc = df[df["dataset"] != "Cappuyns2023"].copy()
val = df[df["dataset"] == "Cappuyns2023"].copy()

assert len(disc) > 0, "No GSE206325 training patients found"
assert len(val) > 0, "No Cappuyns2023 validation patients found"

# n=32 filter: exclude DeathBeforeImaging patients (pat5, pat28, pat29, pat38).
# These 4 patients died before imaging and are not evaluable for treatment
# response. They were mapped to "Non-Responder" during atlas loading, which
# inflates the NR group from 13 → 17 and biases the validation statistics.
_DEATH_BEFORE_IMAGING = {"pat5", "pat28", "pat29", "pat38"}
val = val[~val["patient_id"].astype(str).isin(_DEATH_BEFORE_IMAGING)].copy()
print(f"  Validation set after n=32 filter: {len(val)} patients "
      f"({(val['response']=='Responder').sum()} R / "
      f"{(val['response']=='Non-Responder').sum()} NR)")
assert len(val) == 32, f"Expected 32 Cappuyns patients (19R/13NR), got {len(val)}"

y_disc = (disc["response"] == "Responder").astype(int).values
y_val = (val["response"] == "Responder").astype(int).values

print(f"\n  Discovery : {len(disc)} patients "
      f"({y_disc.sum()} R / {(1-y_disc).sum()} NR) — GSE206325")
print(f"  Validation: {len(val)} patients "
      f"({y_val.sum()} R / {(1-y_val).sum()} NR) — Cappuyns2023")

# ═══════════════════════════════════════════════════════════════════════
# 2.  LAYER DEFINITIONS & FEATURE SELECTION
# ═══════════════════════════════════════════════════════════════════════

layer_prefixes = {
    "L1_tumour_intrinsic": ("cnv_",),
    "L2_regulon":          ("regulon_",),
    "L3_composition":      ("frac_",),
    "L4_io_signatures":    ("sig_",),
    "L5_nmf_programs":     ("nmf_",),
}
if any(c.startswith("gf_cls_") for c in df.columns):
    layer_prefixes["L6_geneformer"] = ("gf_cls_",)
K = 3
IS_TECHNICAL = ("Unknown", "Doublet", "Empty", "doublet", "empty",
                "Proliferating_Proliferating", "Naive_Naive")

def is_technical(col):
    return any(tok in col for tok in IS_TECHNICAL)

def valid(col, data):
    v = data[col]
    return v.notna().sum() >= len(data) * 0.5 and v.nunique(dropna=True) > 1

def pick_top(cols, k, data, yy):
    from scipy.stats import ttest_ind
    stats = []
    for c in cols:
        a = data.loc[yy == 1, c].dropna()
        b = data.loc[yy == 0, c].dropna()
        if len(a) < 2 or len(b) < 2:
            continue
        t, p = ttest_ind(a, b, equal_var=False)
        stats.append((c, abs(t)))
    stats.sort(key=lambda r: r[1], reverse=True)
    return [c for c, _ in stats[:k]]

# Load precomputed feature selection (from 05d) if available,
# otherwise re-select on discovery cohort
sel_file = ML_TAB / "balanced_selected_features.csv"
if sel_file.exists():
    sel_df = pd.read_csv(sel_file)
    selected = {}
    for layer in sel_df["layer"].unique():
        feats = sel_df.loc[sel_df["layer"] == layer, "feature"].tolist()
        # Verify features exist in both dataframes
        feats = [f for f in feats if f in disc.columns and f in val.columns]
        if feats:
            selected[layer] = feats
    print(f"\n  Loaded precomputed features from 05d ({sum(len(v) for v in selected.values())} total)")
else:
    print("\n  Selecting top-3 features per layer on discovery cohort...")
    selected = {}
    for name, prefs in layer_prefixes.items():
        cols = [c for c in disc.columns
                if any(c.startswith(p) for p in prefs) and valid(c, disc)
                and c in val.columns]
        if name == "L3_composition":
            cols = [c for c in cols if not is_technical(c)]
        top = pick_top(cols, K, disc, y_disc)
        if top:
            selected[name] = top

if "L3_composition" in selected:
    l3_clean = [f for f in selected["L3_composition"] if not is_technical(f)]
    if len(l3_clean) < K:
        cols = [c for c in disc.columns if c.startswith("frac_") and c in val.columns
                and valid(c, disc) and not is_technical(c)]
        l3_clean = pick_top(cols, K, disc, y_disc)
    selected["L3_composition"] = l3_clean

for L, feats in selected.items():
    print(f"   {L}: {feats}")


# ═══════════════════════════════════════════════════════════════════════
# 3.  TRAIN ON DISCOVERY, VALIDATE ON CAPPUYNS
# ═══════════════════════════════════════════════════════════════════════

def impute_scale_train_val(X_train, X_val):
    """Impute missing values using training-set medians, then standardize."""
    Xt = pd.DataFrame(X_train)
    medians = Xt.median()
    Xt = Xt.fillna(medians).fillna(0)
    Xv = pd.DataFrame(X_val).fillna(medians).fillna(0)
    scaler = StandardScaler().fit(Xt.values)
    return scaler.transform(Xt.values), scaler.transform(Xv.values), scaler


def bootstrap_auc(y_true, y_score, n_boot=2000, seed=42):
    """Patient-level bootstrap 95% CI for AUC."""
    rng = np.random.RandomState(seed)
    n = len(y_true)
    aucs = []
    for _ in range(n_boot):
        idx = rng.choice(n, n, replace=True)
        if len(np.unique(y_true[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y_true[idx], y_score[idx]))
    aucs = np.array(aucs)
    return np.percentile(aucs, 2.5), np.percentile(aucs, 97.5)


results = []
roc_curves = {}
prediction_table = val[["patient_id", "response", "dataset"]].copy()

# ── Per-layer validation ────────────────────────────────────────────────
for name, feats in selected.items():
    X_tr = disc[feats].values
    X_va = val[feats].values
    X_tr_s, X_va_s, _ = impute_scale_train_val(X_tr, X_va)

    model = LogisticRegression(penalty="l2", C=1.0, max_iter=3000,
                               class_weight="balanced")
    model.fit(X_tr_s, y_disc)

    # Discovery LOOCV (for reference ROC)
    from sklearn.model_selection import LeaveOneOut
    p_disc = np.zeros(len(y_disc))
    loo = LeaveOneOut()
    for tr, te in loo.split(X_tr_s):
        m = LogisticRegression(penalty="l2", C=1.0, max_iter=3000,
                               class_weight="balanced")
        m.fit(X_tr_s[tr], y_disc[tr])
        p_disc[te] = m.predict_proba(X_tr_s[te])[0, 1]

    disc_auc = roc_auc_score(y_disc, p_disc)

    # Validation predictions (full model trained on all discovery)
    p_val = model.predict_proba(X_va_s)[:, 1]
    prediction_table[f"{name}_pred_prob"] = p_val
    if len(np.unique(y_val)) >= 2:
        val_auc = roc_auc_score(y_val, p_val)
        val_prauc = average_precision_score(y_val, p_val)
        lo, hi = bootstrap_auc(y_val, p_val)
        fpr_d, tpr_d, _ = roc_curve(y_disc, p_disc)
        fpr_v, tpr_v, _ = roc_curve(y_val, p_val)
        roc_curves[name] = dict(
            disc_fpr=fpr_d, disc_tpr=tpr_d, disc_auc=disc_auc,
            val_fpr=fpr_v, val_tpr=tpr_v, val_auc=val_auc,
        )

        # Balanced accuracy using 0.5 threshold
        pred_val = (p_val >= 0.5).astype(int)
        bal_acc = balanced_accuracy_score(y_val, pred_val)
    else:
        val_auc = val_prauc = bal_acc = lo = hi = np.nan

    results.append(dict(
        layer=name, n_features=len(feats),
        disc_loocv_auc=disc_auc,
        val_auc=val_auc, val_auc_ci_lo=lo, val_auc_ci_hi=hi,
        val_prauc=val_prauc, val_balanced_acc=bal_acc,
    ))
    print(f"   {name:25s}  disc={disc_auc:.3f}  val={val_auc:.3f} [{lo:.2f}–{hi:.2f}]")

# ── Balanced ensemble ───────────────────────────────────────────────────
all_feats = sum(selected.values(), [])
X_tr_all = disc[all_feats].values
X_va_all = val[all_feats].values
X_tr_s, X_va_s, _ = impute_scale_train_val(X_tr_all, X_va_all)

model_ens = LogisticRegression(penalty="l2", C=0.5, max_iter=3000,
                               class_weight="balanced")
model_ens.fit(X_tr_s, y_disc)

p_disc_ens = np.zeros(len(y_disc))
loo = LeaveOneOut()
for tr, te in loo.split(X_tr_s):
    m = LogisticRegression(penalty="l2", C=0.5, max_iter=3000,
                           class_weight="balanced")
    m.fit(X_tr_s[tr], y_disc[tr])
    p_disc_ens[te] = m.predict_proba(X_tr_s[te])[0, 1]

disc_ens_auc = roc_auc_score(y_disc, p_disc_ens)
p_val_ens = model_ens.predict_proba(X_va_s)[:, 1]
prediction_table["Balanced_ensemble_pred_prob"] = p_val_ens

if len(np.unique(y_val)) >= 2:
    val_ens_auc = roc_auc_score(y_val, p_val_ens)
    val_ens_prauc = average_precision_score(y_val, p_val_ens)
    lo_e, hi_e = bootstrap_auc(y_val, p_val_ens)
    pred_ens = (p_val_ens >= 0.5).astype(int)
    bal_acc_e = balanced_accuracy_score(y_val, pred_ens)
    fpr_d, tpr_d, _ = roc_curve(y_disc, p_disc_ens)
    fpr_v, tpr_v, _ = roc_curve(y_val, p_val_ens)
    roc_curves["Ensemble"] = dict(
        disc_fpr=fpr_d, disc_tpr=tpr_d, disc_auc=disc_ens_auc,
        val_fpr=fpr_v, val_tpr=tpr_v, val_auc=val_ens_auc,
    )
else:
    val_ens_auc = val_ens_prauc = lo_e = hi_e = bal_acc_e = np.nan

results.append(dict(
    layer="Balanced ensemble", n_features=len(all_feats),
    disc_loocv_auc=disc_ens_auc,
    val_auc=val_ens_auc, val_auc_ci_lo=lo_e, val_auc_ci_hi=hi_e,
    val_prauc=val_ens_prauc, val_balanced_acc=bal_acc_e,
))
print(f"   {'Ensemble':25s}  disc={disc_ens_auc:.3f}  val={val_ens_auc:.3f} [{lo_e:.2f}–{hi_e:.2f}]")

# ── Save results ────────────────────────────────────────────────────────
res_df = pd.DataFrame(results)
res_df.to_csv(TAB / "cappuyns_validation_results.csv", index=False)
prediction_table.to_csv(TAB / "cappuyns_per_patient_predictions.csv", index=False)
if "L3_composition_pred_prob" in prediction_table.columns:
    val["L3_pred_prob"] = prediction_table["L3_composition_pred_prob"].values


# ═══════════════════════════════════════════════════════════════════════
# 4.  FIGURES
# ═══════════════════════════════════════════════════════════════════════

layer_palette = {
    "L1_tumour_intrinsic": "#882255",
    "L2_regulon":          "#DDCC77",
    "L3_composition":      "#44AA99",
    "L4_io_signatures":    "#332288",
    "L5_nmf_programs":     "#AA4499",
    "L6_geneformer":       "#117733",
    "Ensemble":            "#000000",
}

# ── Fig A: L3 ROC overlay (discovery LOOCV vs Cappuyns external) ────────
if "L3_composition" in roc_curves:
    rc = roc_curves["L3_composition"]
    fig, ax = plt.subplots(figsize=(5.5, 5.2), constrained_layout=True)
    ax.plot(rc["disc_fpr"], rc["disc_tpr"],
            color=PAL["Responder"], linewidth=2, linestyle="--",
            label=f"Discovery LOOCV  AUC={rc['disc_auc']:.2f}")
    ax.plot(rc["val_fpr"], rc["val_tpr"],
            color="#44AA99", linewidth=2.5,
            label=f"Validation (Cappuyns)  AUC={rc['val_auc']:.2f}")
    ax.plot([0, 1], [0, 1], ":", color="grey", linewidth=0.8)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("L3 composition: discovery → Cappuyns validation",
                 fontsize=10, fontweight="bold")
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    ax.set_aspect("equal")
    for fmt in ("png", "pdf"):
        fig.savefig(FIG / f"03c_validation_roc.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)

# ── Fig B: Per-layer validation AUC comparison ────────────────────────
fig, ax = plt.subplots(figsize=(7, 3.8), constrained_layout=True)
layers = [r["layer"] for r in results]
disc_aucs = [r["disc_loocv_auc"] for r in results]
val_aucs = [r["val_auc"] for r in results]
x = np.arange(len(layers))
w = 0.35
bars_d = ax.barh(x + w/2, disc_aucs, w, color="#4477AA", edgecolor="black",
                 linewidth=0.5, label="Discovery (LOOCV)")
bars_v = ax.barh(x - w/2, val_aucs, w, color="#44AA99", edgecolor="black",
                 linewidth=0.5, label="Validation (Cappuyns)")
ax.set_yticks(x)
ax.set_yticklabels([l.replace("_", " ") for l in layers], fontsize=8)
ax.set_xlabel("ROC-AUC")
ax.axvline(0.5, color="grey", linewidth=0.8, linestyle=":")
ax.set_title("Per-layer AUC: discovery vs external validation",
             fontsize=10, fontweight="bold")
ax.legend(loc="lower right", fontsize=8).set_zorder(10)
for i, (d, v) in enumerate(zip(disc_aucs, val_aucs)):
    ax.text(d + 0.01, i + w/2, f"{d:.2f}", va="center", fontsize=7)
    if not np.isnan(v):
        ax.text(v + 0.01, i - w/2, f"{v:.2f}", va="center", fontsize=7)
for fmt in ("png", "pdf"):
    fig.savefig(FIG / f"03c_validation_perlayer.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)

# ── Fig C: L3 prediction waterfall (Cappuyns patients) ──────────────
if "L3_pred_prob" in val.columns:
    wf = val[["patient_id", "response", "L3_pred_prob"]].sort_values("L3_pred_prob")
    fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
    colors = [PAL["Responder"] if r == "Responder" else PAL["Non_Responder"]
              for r in wf["response"]]
    ax.bar(range(len(wf)), wf["L3_pred_prob"], color=colors, edgecolor="black",
           linewidth=0.4, width=0.8)
    ax.axhline(0.5, color="black", linewidth=1, linestyle="--", alpha=0.6)
    ax.set_ylabel("L3 predicted probability (Responder)")
    ax.set_xlabel("Cappuyns patients (ranked)")
    ax.set_title("L3 composition predictor — Cappuyns external validation",
                 fontweight="bold", fontsize=10)
    ax.set_xticks(range(len(wf)))
    ax.set_xticklabels(wf["patient_id"], rotation=90, fontsize=6)
    from matplotlib.patches import Patch
    ax.legend([Patch(color=PAL["Responder"]), Patch(color=PAL["Non_Responder"])],
              ["Responder", "Non-Responder"], loc="upper left", fontsize=8)
    for fmt in ("png", "pdf"):
        fig.savefig(FIG / f"03c_validation_waterfall.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)

# ── Summary ─────────────────────────────────────────────────────────────
print(f"\n{'='*70}")
print(f"[DONE] Validation results:")
print(res_df.to_string(index=False))
print(f"\n  Tables: {TAB}")
print(f"  Figures: {FIG}")
print(f"{'='*70}")
