#!/usr/bin/env python3
"""
05m_lasso_C_sensitivity.py — T3.2

Sensitivity of the L3 and ensemble LOOCV-AUC to the LASSO/L2 regularisation
strength C. The headline manuscript result uses L3 C=1.0 (L2) and Ensemble
C=0.3 (L1). This script sweeps C over a log grid and records:
- LOOCV-AUC on discovery (GSE206325)
- Transfer AUC on Cappuyns2023
- Number of retained (non-zero) features under L1
- Balanced accuracy and PR-AUC

Outputs
-------
out/03_response_modeling/tables/lasso_C_sensitivity.csv
out/03_response_modeling/figures/05m_lasso_C_sensitivity.{png,pdf}
"""
from __future__ import annotations
import sys, warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import (roc_auc_score, average_precision_score,
                              balanced_accuracy_score)

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
from config.publication_style import apply_style, PAL
apply_style()

TAB_IN = PROJECT / "out" / "05_multilayer" / "tables"
TAB = PROJECT / "out" / "03_response_modeling" / "tables"
FIG = PROJECT / "out" / "03_response_modeling" / "figures"
TAB.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)

# ── Load features ────────────────────────────────────────────────────
df = pd.read_csv(TAB_IN / "patient_multilayer_features.csv")
nmf = pd.read_csv(TAB_IN / "nmf_programs_patient_activity.csv")
nmf_cols = [c for c in nmf.columns if c.startswith("NMF")]
df = df.merge(nmf[["patient_id"] + nmf_cols].rename(
    columns={c: f"nmf_{c}" for c in nmf_cols}), on="patient_id", how="left")
# positive = Responder (both cohorts use 'Responder' label)
df["y"] = (df["response"] == "Responder").astype(int)
disc = df[df["dataset"] == "GSE206325"].reset_index(drop=True)
cap = df[df["dataset"] == "Cappuyns2023"].reset_index(drop=True)
print(f"[05m] discovery n={len(disc)} (R={disc['y'].sum()})  "
      f"cappuyns n={len(cap)} (R={cap['y'].sum()})")

# ── Feature selection helpers ────────────────────────────────────────
TECH = ("frac_Unknown", "frac_Doublet", "frac_Empty")
def layer_cols(prefixes, exclude_technical=True):
    out = []
    for c in df.columns:
        if c in ("patient_id", "response", "dataset", "y"): continue
        if not any(c.startswith(p) for p in prefixes): continue
        if exclude_technical and any(c.startswith(t) for t in TECH): continue
        if df[c].dtype.kind not in "fi": continue
        if df[c].notna().sum() < len(df) * 0.5: continue
        if df[c].nunique(dropna=True) <= 1: continue
        out.append(c)
    return out

def impute_scale_fit(X):
    imp = SimpleImputer(strategy="median").fit(X)
    sc = StandardScaler().fit(imp.transform(X))
    return imp, sc

def impute_scale_apply(X, imp, sc):
    return sc.transform(imp.transform(X))

def loocv_auc(X, y, C, penalty):
    solver = "liblinear" if penalty == "l1" else "lbfgs"
    p = np.zeros(len(y))
    for tr, te in LeaveOneOut().split(X):
        clf = LogisticRegression(penalty=penalty, C=C, solver=solver,
                                 max_iter=3000, class_weight="balanced")
        clf.fit(X[tr], y[tr])
        p[te] = clf.predict_proba(X[te])[0, 1]
    return p

# ── Build feature sets ───────────────────────────────────────────────
L3_cols = layer_cols(("frac_",))
ensemble_cols = layer_cols(("cnv_", "regulon_", "frac_", "sig_", "nmf_"))
print(f"   L3 features: {len(L3_cols)}  Ensemble features: {len(ensemble_cols)}")

# Fit imputer/scaler on discovery; apply to Cappuyns for transfer
def prepare(cols):
    Xd = disc[cols].values
    imp, sc = impute_scale_fit(Xd)
    Xd_s = impute_scale_apply(Xd, imp, sc)
    cols_cap = [c for c in cols if c in cap.columns]
    if len(cols_cap) < len(cols) * 0.5:
        Xc_s, cap_mask = None, None
    else:
        # align by intersection, zero-fill missing
        Xc = cap[cols].copy()
        for c in cols:
            if c not in cap.columns:
                Xc[c] = np.nan
        Xc_s = impute_scale_apply(Xc.values, imp, sc)
        cap_mask = cap["y"].values
    return Xd_s, Xc_s

Xd_L3, Xc_L3 = prepare(L3_cols)
Xd_EN, Xc_EN = prepare(ensemble_cols)
y_disc = disc["y"].values
y_cap = cap["y"].values

# ── Sweep ────────────────────────────────────────────────────────────
C_GRID = np.logspace(-2, 1.5, 15)  # 0.01 … ~31
rows = []
for pen_label, penalty, reg_set in [
        ("L2 ridge (L3)", "l2", [("L3_composition", Xd_L3, Xc_L3, L3_cols)]),
        ("L1 LASSO (Ensemble)", "l1", [("Ensemble", Xd_EN, Xc_EN, ensemble_cols)]),
        ("L1 LASSO (L3)", "l1", [("L3_composition", Xd_L3, Xc_L3, L3_cols)]),
        ("L2 ridge (Ensemble)", "l2", [("Ensemble", Xd_EN, Xc_EN, ensemble_cols)])]:
    for model_name, Xd_s, Xc_s, cols_used in reg_set:
        for C in C_GRID:
            # LOOCV discovery
            p = loocv_auc(Xd_s, y_disc, C=C, penalty=penalty)
            auc_d = roc_auc_score(y_disc, p)
            pr_d  = average_precision_score(y_disc, p)
            ba_d  = balanced_accuracy_score(y_disc, (p > 0.5).astype(int))
            # Fit once on full discovery, apply to Cappuyns
            solver = "liblinear" if penalty == "l1" else "lbfgs"
            clf = LogisticRegression(penalty=penalty, C=C, solver=solver,
                                     max_iter=3000, class_weight="balanced")
            clf.fit(Xd_s, y_disc)
            nz = int((np.abs(clf.coef_) > 1e-8).sum())
            if Xc_s is not None and y_cap.sum() not in (0, len(y_cap)):
                p_cap = clf.predict_proba(Xc_s)[:, 1]
                auc_c = roc_auc_score(y_cap, p_cap)
                pr_c  = average_precision_score(y_cap, p_cap)
                ba_c  = balanced_accuracy_score(y_cap, (p_cap > 0.5).astype(int))
            else:
                auc_c = pr_c = ba_c = np.nan
            rows.append(dict(
                model=pen_label, target=model_name, penalty=penalty,
                C=float(C), n_features_input=len(cols_used),
                n_features_retained=nz,
                discovery_loocv_auc=auc_d, discovery_pr_auc=pr_d, discovery_balacc=ba_d,
                cappuyns_transfer_auc=auc_c, cappuyns_pr_auc=pr_c, cappuyns_balacc=ba_c,
            ))

out = pd.DataFrame(rows)
out.to_csv(TAB / "lasso_C_sensitivity.csv", index=False)
print(f"\nwrote {TAB/'lasso_C_sensitivity.csv'} ({len(out)} rows)")
print(out.groupby("model")[["discovery_loocv_auc", "cappuyns_transfer_auc", "n_features_retained"]]
        .agg(["min", "median", "max"]).round(3).to_string())

# ── Figure ──────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), constrained_layout=True)

style = {"L2 ridge (L3)":       dict(color="#4477AA", ls="-",  m="o"),
         "L1 LASSO (L3)":        dict(color="#4477AA", ls="--", m="s"),
         "L1 LASSO (Ensemble)": dict(color="#CC6677", ls="-",  m="o"),
         "L2 ridge (Ensemble)":  dict(color="#CC6677", ls="--", m="s")}
HEADLINE_L3 = 1.0
HEADLINE_EN = 0.3

# A — discovery LOOCV AUC vs C
axA = axes[0]
for name, sub in out.groupby("model"):
    s = sub.sort_values("C")
    st = style[name]
    axA.plot(s["C"], s["discovery_loocv_auc"], color=st["color"], ls=st["ls"],
             marker=st["m"], ms=4, lw=1.5, label=name)
axA.axvline(HEADLINE_L3, color="#4477AA", lw=0.6, ls=":")
axA.axvline(HEADLINE_EN, color="#CC6677", lw=0.6, ls=":")
axA.set_xscale("log")
axA.set_xlabel("Regularisation strength C", fontsize=9)
axA.set_ylabel("Discovery LOOCV AUC (GSE206325, n=24)", fontsize=9)
axA.set_title("A  Discovery AUC vs C",
              loc="left", fontsize=10.5, fontweight="bold")
axA.axhline(0.5, color="grey", lw=0.6, ls="--")
axA.legend(frameon=False, fontsize=7.5, loc="lower right")

# B — Cappuyns transfer AUC vs C
axB = axes[1]
for name, sub in out.groupby("model"):
    s = sub.sort_values("C")
    st = style[name]
    axB.plot(s["C"], s["cappuyns_transfer_auc"], color=st["color"], ls=st["ls"],
             marker=st["m"], ms=4, lw=1.5, label=name)
axB.axvline(HEADLINE_L3, color="#4477AA", lw=0.6, ls=":")
axB.axvline(HEADLINE_EN, color="#CC6677", lw=0.6, ls=":")
axB.set_xscale("log")
axB.set_xlabel("Regularisation strength C", fontsize=9)
axB.set_ylabel("Cappuyns transfer AUC (n=32)", fontsize=9)
axB.set_title("B  Cappuyns transfer vs C",
              loc="left", fontsize=10.5, fontweight="bold")
axB.axhline(0.5, color="grey", lw=0.6, ls="--")
axB.legend(frameon=False, fontsize=7.5, loc="lower right")

# C — retained features vs C (L1 only)
axC = axes[2]
for name, sub in out.groupby("model"):
    if "L1" not in name: continue
    s = sub.sort_values("C")
    st = style[name]
    axC.plot(s["C"], s["n_features_retained"], color=st["color"], ls=st["ls"],
             marker=st["m"], ms=4, lw=1.5, label=name)
axC.axvline(HEADLINE_EN, color="#CC6677", lw=0.6, ls=":")
axC.set_xscale("log"); axC.set_yscale("log")
axC.set_xlabel("Regularisation strength C", fontsize=9)
axC.set_ylabel("# non-zero features", fontsize=9)
axC.set_title("C  L1 sparsity vs C",
              loc="left", fontsize=10.5, fontweight="bold")
axC.legend(frameon=False, fontsize=7.5, loc="lower right")

fig.suptitle("Figure: Sensitivity of LOOCV/transfer AUC to regularisation strength C",
             fontsize=12, fontweight="bold")
fig.savefig(FIG / "05m_lasso_C_sensitivity.png", dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(FIG / "05m_lasso_C_sensitivity.pdf", bbox_inches="tight", facecolor="white")
print(f"\nwrote {FIG/'05m_lasso_C_sensitivity.png'}")

# ── Headline-C robustness window ────────────────────────────────────
# For each model, find range of C where discovery AUC stays within 5% of max
summary_rows = []
for name, sub in out.groupby("model"):
    s = sub.sort_values("C")
    mx = s["discovery_loocv_auc"].max()
    mask = s["discovery_loocv_auc"] >= mx - 0.05
    c_range = s.loc[mask, "C"]
    summary_rows.append(dict(
        model=name,
        C_at_max_disc=float(s.loc[s["discovery_loocv_auc"].idxmax(), "C"]),
        max_discovery_auc=float(mx),
        C_window_lo=float(c_range.min()),
        C_window_hi=float(c_range.max()),
        window_fraction_of_grid=float(mask.mean()),
    ))
summary = pd.DataFrame(summary_rows)
summary.to_csv(TAB / "lasso_C_sensitivity_summary.csv", index=False)
print("\nRobustness window (discovery AUC within 0.05 of max):")
print(summary.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
