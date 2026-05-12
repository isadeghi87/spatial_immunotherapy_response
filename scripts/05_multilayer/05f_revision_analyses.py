#!/usr/bin/env python3
"""
05f_revision_analyses.py — Second-revision statistical strengthening.

Executes the seven minimal actions requested by the second-round review:

1.  Label-permutation null (1000 shuffles) on the balanced L3 and 15-feature
    ensemble pipelines.
2.  L3 re-run after dropping every 'Unknown' / 'Doublet' / 'Empty' feature,
    with bootstrap 95 % CI, permutation p, PR-AUC and balanced accuracy.
3.  PR-AUC and balanced accuracy for every layer (and for the ensemble).
4.  Real IMbrave150 CD8 IHC x arm interaction analysis on the per-patient
    data extracted from MOESM4 (n = 177 biomarker-evaluable patients).
5.  BH-adjusted NMF programme p-values (10 programmes) and an updated
    NMF summary table.
6-7. (handled in manuscript rewrite) Marker-gene definition of
    Dysfunctional-effector CD8 and a one-sentence power statement under
    §3.7 — both read from the outputs of this script.

Outputs:
    out/05_multilayer/tables/revision_layer_metrics.csv
    out/05_multilayer/tables/revision_permutation_null.csv
    out/05_multilayer/tables/revision_L3_without_unknown.csv
    out/05_multilayer/tables/revision_nmf_bh_adjusted.csv
    out/04_complementary/tables/revision_imbrave150_cd8_arm.csv
    out/05_multilayer/figures/Fig6A_balanced_layer_metrics.png
    out/04_complementary/figures/ExtendedData_Fig4_imbrave150_real.png
"""
from __future__ import annotations
import warnings; warnings.filterwarnings("ignore")
import sys
from pathlib import Path
PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))

import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import (roc_auc_score, average_precision_score,
                              balanced_accuracy_score)

RNG = np.random.default_rng(0)
N_PERM = 1000
N_BOOT = 2000

OUT5 = PROJECT / "out" / "05_multilayer"
TAB5 = OUT5 / "tables"
FIG5 = OUT5 / "figures"
OUT4 = PROJECT / "out" / "04_complementary"
TAB4 = OUT4 / "tables"; TAB4.mkdir(parents=True, exist_ok=True)
FIG4 = OUT4 / "figures"; FIG4.mkdir(parents=True, exist_ok=True)

# ──────────────────────────────────────────────────────────────────────
# Load patient-level feature matrix (same as 05d)
# ──────────────────────────────────────────────────────────────────────
df = pd.read_csv(TAB5 / "patient_multilayer_features.csv")
nmf = pd.read_csv(TAB5 / "nmf_programs_patient_activity.csv")
nmf_cols = [c for c in nmf.columns if c.startswith("NMF")]
nmf = nmf.rename(columns={c: f"nmf_{c}" for c in nmf_cols})
df = df.merge(nmf[["patient_id"] + [f"nmf_{c}" for c in nmf_cols]],
               on="patient_id", how="left")
if "dataset" not in df.columns:
    df["dataset"] = np.where(df["patient_id"].astype(str).str.startswith("pat"),
                             "Cappuyns2023", "GSE206325")
df = df[df["dataset"].astype(str) != "Cappuyns2023"].copy()
y = (df["response"] == "Responder").astype(int).values
n = len(y)
print(f"[05f] cohort n={n}  (R={y.sum()} / NR={n - y.sum()})")

LAYERS = {
    "L1_tumour_intrinsic": ("cnv_",),
    "L2_regulon":          ("regulon_",),
    "L3_composition":      ("frac_",),
    "L4_io_signatures":    ("sig_",),
    "L5_nmf_programs":     ("nmf_",),
}
if any(c.startswith("gf_cls_") for c in df.columns):
    LAYERS["L6_geneformer"] = ("gf_cls_",)
K = 3

def valid(col):
    v = df[col]
    return v.notna().sum() >= n * 0.5 and v.nunique(dropna=True) > 1

# Mask to strip flagged technical features from L3 composition
IS_TECHNICAL = ("Unknown", "Doublet", "Empty", "doublet", "empty",
                "Proliferating_Proliferating", "Naive_Naive")

def is_technical(col):
    return any(tok in col for tok in IS_TECHNICAL)

def layer_cols(prefs, exclude_technical=False):
    cols = [c for c in df.columns
             if any(c.startswith(p) for p in prefs) and valid(c)]
    if exclude_technical:
        cols = [c for c in cols if not is_technical(c)]
    return cols

def pick_top(cols, yy, k=K):
    stats = []
    for c in cols:
        a = df.loc[yy == 1, c].dropna()
        b = df.loc[yy == 0, c].dropna()
        if len(a) < 2 or len(b) < 2:
            continue
        t, p = ttest_ind(a, b, equal_var=False)
        stats.append((c, abs(t)))
    stats.sort(key=lambda r: r[1], reverse=True)
    return [c for c, _ in stats[:k]]

def impute_scale(X):
    X = pd.DataFrame(X)
    X = X.fillna(X.median()).fillna(0)
    return StandardScaler().fit_transform(X.values)

def loocv_probs(X, yy, C=1.0):
    loo = LeaveOneOut(); p = np.zeros(len(yy))
    for tr, te in loo.split(X):
        m = LogisticRegression(penalty="l2", C=C, max_iter=3000,
                                class_weight="balanced")
        m.fit(X[tr], yy[tr])
        p[te] = m.predict_proba(X[te])[0, 1]
    return p

def metrics(yy, p):
    auc = roc_auc_score(yy, p)
    pr  = average_precision_score(yy, p)
    pred = (p >= 0.5).astype(int)
    ba  = balanced_accuracy_score(yy, pred)
    return auc, pr, ba

# ──────────────────────────────────────────────────────────────────────
# 1 + 2 + 3.  Per-layer metrics with bootstrap CIs and permutation null
# ──────────────────────────────────────────────────────────────────────
def evaluate(cols_by_layer, yy, C_layer=1.0, C_ens=0.5):
    """Compute AUC/PR-AUC/BA per layer + 15-feature ensemble."""
    out = {}
    for name, feats in cols_by_layer.items():
        X = impute_scale(df[feats].values)
        p = loocv_probs(X, yy, C=C_layer)
        out[name] = (feats, p)
    # ensemble
    all_feats = sum([f for f, _ in out.values()], [])
    Xall = impute_scale(df[all_feats].values)
    p_ens = loocv_probs(Xall, yy, C=C_ens)
    out["Balanced ensemble"] = (all_feats, p_ens)
    return out

print("[1] Selecting top-3 features per layer on observed labels ...")
sel_obs = {name: pick_top(layer_cols(prefs), y)
            for name, prefs in LAYERS.items()}
# L3-no-technical variant
sel_obs_noTech = dict(sel_obs)
sel_obs_noTech["L3_composition"] = pick_top(
    layer_cols(LAYERS["L3_composition"], exclude_technical=True), y)
print("   L3 (with tech): ", sel_obs["L3_composition"])
print("   L3 (no  tech):  ", sel_obs_noTech["L3_composition"])

# ── Bootstrap CI for the observed pipeline ────────────────────────────
print(f"[CI] {N_BOOT}-iteration patient-level bootstrap ...")
def bootstrap_metrics(feats_by_layer, yy, seed=0):
    rng = np.random.default_rng(seed)
    eval_probs = evaluate(feats_by_layer, yy)
    # collect OOB / full-LOOCV probs; bootstrap over indices
    boots = {k: {"auc": [], "pr": [], "ba": []}
              for k in list(feats_by_layer) + ["Balanced ensemble"]}
    for b in range(N_BOOT):
        idx = rng.integers(0, n, n)
        if yy[idx].sum() < 2 or (yy[idx] == 0).sum() < 2:
            continue
        for k, (feats, p) in eval_probs.items():
            try:
                auc = roc_auc_score(yy[idx], p[idx])
                pr  = average_precision_score(yy[idx], p[idx])
                ba  = balanced_accuracy_score(
                    yy[idx], (p[idx] >= 0.5).astype(int))
                boots[k]["auc"].append(auc)
                boots[k]["pr"].append(pr)
                boots[k]["ba"].append(ba)
            except Exception:
                pass
    rows = []
    for k, (feats, p) in eval_probs.items():
        auc, pr, ba = metrics(yy, p)
        rows.append(dict(
            layer=k, n_features=len(feats),
            auc=auc,
            auc_ci_lo=np.percentile(boots[k]["auc"], 2.5),
            auc_ci_hi=np.percentile(boots[k]["auc"], 97.5),
            pr_auc=pr,
            pr_ci_lo=np.percentile(boots[k]["pr"], 2.5),
            pr_ci_hi=np.percentile(boots[k]["pr"], 97.5),
            balanced_acc=ba,
            ba_ci_lo=np.percentile(boots[k]["ba"], 2.5),
            ba_ci_hi=np.percentile(boots[k]["ba"], 97.5),
            prevalence_baseline=yy.mean()))
    return pd.DataFrame(rows), eval_probs

layer_tbl_obs, probs_obs = bootstrap_metrics(sel_obs, y, seed=0)
layer_tbl_obs.to_csv(TAB5 / "revision_layer_metrics.csv", index=False)
print(layer_tbl_obs.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

# ── L3 without technical features ─────────────────────────────────────
print("[2] L3 after dropping Unknown/Doublets/Empty features ...")
l3_noTech_tbl, _ = bootstrap_metrics(
    {"L3_composition_noTech": sel_obs_noTech["L3_composition"]}, y, seed=1)
l3_noTech_tbl.to_csv(TAB5 / "revision_L3_without_unknown.csv", index=False)
print(l3_noTech_tbl.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

# ──────────────────────────────────────────────────────────────────────
# 1'. Permutation null — L3 and ensemble
# ──────────────────────────────────────────────────────────────────────
print(f"[Null] {N_PERM} label-shuffle permutations ...")
null_auc = {name: [] for name in LAYERS}
null_auc["Balanced ensemble"] = []
null_auc["L3_composition_noTech"] = []
cols_L3 = layer_cols(LAYERS["L3_composition"])
cols_L3_noTech = layer_cols(LAYERS["L3_composition"], exclude_technical=True)
all_cols_by_layer = {name: layer_cols(prefs) for name, prefs in LAYERS.items()}

for b in range(N_PERM):
    yperm = RNG.permutation(y)
    sel_perm = {}
    for name in LAYERS:
        feats = pick_top(all_cols_by_layer[name], yperm)
        sel_perm[name] = feats
        if len(feats) < 2:
            continue
        Xp = impute_scale(df[feats].values)
        null_auc[name].append(roc_auc_score(yperm, loocv_probs(Xp, yperm, C=1.0)))

    feats_nt = pick_top(cols_L3_noTech, yperm)
    if len(feats_nt) >= 2:
        Xnt = impute_scale(df[feats_nt].values)
        null_auc["L3_composition_noTech"].append(
            roc_auc_score(yperm, loocv_probs(Xnt, yperm, C=1.0)))

    all_feats = sum(sel_perm.values(), [])
    if len(all_feats) >= 2:
        Xa = impute_scale(df[all_feats].values)
        null_auc["Balanced ensemble"].append(
            roc_auc_score(yperm, loocv_probs(Xa, yperm, C=0.5)))
    if (b + 1) % 100 == 0:
        print(f"   perm {b + 1}/{N_PERM}")

def perm_p(null_arr, obs):
    null_arr = np.asarray(null_arr)
    # one-sided (right tail): fraction of null >= observed
    return (np.sum(null_arr >= obs) + 1) / (len(null_arr) + 1)

obs_auc = {row.layer: row.auc for row in layer_tbl_obs.itertuples()}
obs_auc["L3_composition_noTech"] = l3_noTech_tbl.loc[
    l3_noTech_tbl.layer == "L3_composition_noTech", "auc"
].values[0]

null_rows = []
for target, values in null_auc.items():
    if target not in obs_auc or len(values) == 0:
        continue
    null_rows.append(
        dict(
            target=target,
            obs=obs_auc[target],
            null_median=np.median(values),
            null_95=np.percentile(values, 95),
            perm_p=perm_p(values, obs_auc[target]),
        )
    )
null_df = pd.DataFrame(null_rows)
null_df.to_csv(TAB5 / "revision_permutation_null.csv", index=False)
print(null_df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

# ──────────────────────────────────────────────────────────────────────
# 5. BH-adjusted NMF programme p-values
# ──────────────────────────────────────────────────────────────────────
print("[5] BH-adjusting NMF programme p-values ...")
nmf_rows = []
for col in [c for c in df.columns if c.startswith("nmf_")]:
    a = df.loc[y == 1, col].dropna()
    b = df.loc[y == 0, col].dropna()
    if len(a) < 2 or len(b) < 2: continue
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    # Cliff's delta
    na, nb = len(a), len(b)
    gt = sum(x > yy_ for x in a for yy_ in b)
    lt = sum(x < yy_ for x in a for yy_ in b)
    delta = (gt - lt) / (na * nb)
    nmf_rows.append(dict(programme=col.replace("nmf_NMF", "NMF"),
                          n_R=na, n_NR=nb,
                          mean_R=a.mean(), mean_NR=b.mean(),
                          cliff_delta=delta,
                          p_raw=p))
nmf_tbl = pd.DataFrame(nmf_rows).sort_values("p_raw")
nmf_tbl["p_BH"] = multipletests(nmf_tbl["p_raw"], method="fdr_bh")[1]
nmf_tbl.to_csv(TAB5 / "revision_nmf_bh_adjusted.csv", index=False)
print(nmf_tbl.to_string(index=False, float_format=lambda x: f"{x:.3g}"))

# ──────────────────────────────────────────────────────────────────────
# 4.  Real IMbrave150 analysis (n = 177, per-patient IHC + arm + PFS/OS)
# ──────────────────────────────────────────────────────────────────────
print("[4] IMbrave150 real per-patient analysis (arm × CD8) ...")
imb_f = PROJECT / "data" / "clinical_trials" / "IMbrave150" / \
        "imbrave150_fig3b_patient_data.csv"
imb = pd.read_csv(imb_f)
# column convention in MOESM: 'censor=1, event=0' — convert to
# lifelines event_observed where 1 = event occurred
imb["event_pfs"] = 1 - imb["pfs_event"].astype(int)
imb["event_os"]  = 1 - imb["os_event"].astype(int)
imb["cd8_high"]  = (imb["cd8_hl"].astype(str).str.strip() == "High").astype(int)

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

arm_rows = []
for arm in ["Atezo+Bev", "Sorafenib"]:
    d = imb.loc[imb["arm"].str.strip() == arm]
    h = d[d["cd8_high"] == 1]; lo = d[d["cd8_high"] == 0]
    if len(h) < 5 or len(lo) < 5: continue
    for endpoint, t, e in [("PFS", "pfs_m", "event_pfs"),
                             ("OS",  "os_m",  "event_os")]:
        res = logrank_test(h[t], lo[t], event_observed_A=h[e],
                            event_observed_B=lo[e])
        arm_rows.append(dict(arm=arm, endpoint=endpoint,
                              n_high=len(h), n_low=len(lo),
                              median_high=h[t].median(),
                              median_low=lo[t].median(),
                              logrank_p=res.p_value,
                              test_statistic=res.test_statistic))

arm_tbl = pd.DataFrame(arm_rows)

# Cox with arm × CD8_high interaction
cox_rows = []
for endpoint, t, e in [("PFS", "pfs_m", "event_pfs"),
                         ("OS",  "os_m",  "event_os")]:
    c = imb[[t, e, "cd8_high"]].copy()
    c["arm_AtezoBev"] = (imb["arm"].str.strip() == "Atezo+Bev").astype(int)
    c["arm_x_cd8"]    = c["arm_AtezoBev"] * c["cd8_high"]
    c = c.dropna()
    cph = CoxPHFitter().fit(c, duration_col=t, event_col=e)
    for var in ["cd8_high", "arm_AtezoBev", "arm_x_cd8"]:
        cox_rows.append(dict(endpoint=endpoint, covariate=var,
                              hr=np.exp(cph.params_[var]),
                              hr_lo=np.exp(cph.confidence_intervals_.loc[var].iloc[0]),
                              hr_hi=np.exp(cph.confidence_intervals_.loc[var].iloc[1]),
                              p=cph.summary.loc[var, "p"]))

cox_tbl = pd.DataFrame(cox_rows)

print("\n  Within-arm log-rank (CD8 high vs low, IHC, median split):")
print(arm_tbl.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print("\n  Cox model (CD8_high + arm + arm × CD8 interaction):")
print(cox_tbl.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

arm_tbl.to_csv(TAB4 / "revision_imbrave150_cd8_arm.csv", index=False)
cox_tbl.to_csv(TAB4 / "revision_imbrave150_cox_interaction.csv", index=False)

# ── Plots ────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(13.3, 8.4), constrained_layout=True)
for j, (arm, ax_pfs, ax_os) in enumerate([
        ("Atezo+Bev", axes[0, 0], axes[0, 1]),
        ("Sorafenib", axes[1, 0], axes[1, 1]) ]):
    d = imb.loc[imb["arm"].str.strip() == arm]
    for endpoint, t, e, ax in [("PFS", "pfs_m", "event_pfs", ax_pfs),
                                 ("OS",  "os_m",  "event_os",  ax_os)]:
        h = d[d["cd8_high"] == 1]; lo = d[d["cd8_high"] == 0]
        kmf = KaplanMeierFitter()
        for strat, col, lbl in [(h, "#0077BB", "CD8-High"),
                                  (lo, "#CC3311", "CD8-Low")]:
            kmf.fit(strat[t], strat[e],
                     label=f"{lbl} (n={len(strat)})")
            kmf.plot_survival_function(ax=ax, ci_show=True, linewidth=2.4)
        res = logrank_test(h[t], lo[t], event_observed_A=h[e],
                            event_observed_B=lo[e])
        ax.set_title(f"{arm} — {endpoint}\nlog-rank p = {res.p_value:.3f}",
                      fontweight="bold", fontsize=11)
        ax.set_xlabel("Months"); ax.set_ylabel(f"{endpoint} probability")
        ax.set_ylim(-0.02, 1.02)
        ax.grid(alpha=0.3, linewidth=0.5)
fig.suptitle("IMbrave150 CD8 IHC stratification — per arm, real per-patient data (n = 177)",
              fontweight="bold", fontsize=13, y=1.03)
for fmt in ("png", "pdf"):
    fig.savefig(FIG4 / f"ExtendedData_Fig4_imbrave150_real.{fmt}",
                 dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved ExtendedData_Fig4_imbrave150_real.png")

# ──────────────────────────────────────────────────────────────────────
# Fig 6A rewrite — AUC + PR-AUC + balanced accuracy with 95 % CIs,
#                   plus permutation-null line
# ──────────────────────────────────────────────────────────────────────
met = layer_tbl_obs.copy()
# prepend L3-noTech row for panel
met_noTech_row = l3_noTech_tbl.copy()
met_noTech_row["layer"] = "L3_composition_noTech"
met = pd.concat([met, met_noTech_row], ignore_index=True)

# Ordering (display from worst to best)
order = ["L2_regulon", "L5_nmf_programs", "L1_tumour_intrinsic",
         "L4_io_signatures", "L6_geneformer", "L3_composition",
         "L3_composition_noTech", "Balanced ensemble"]
order = [layer for layer in order if layer in set(met["layer"])]
met = met.set_index("layer").loc[order].reset_index()
met["label"] = (met["layer"]
                 .str.replace("_", " ")
                 .str.replace("L1 tumour intrinsic", "L1 CNV")
                 .str.replace("L2 regulon", "L2 Regulon")
                 .str.replace("L6 geneformer", "L6 Geneformer")
                 .str.replace("L3 composition noTech", "L3 (no Unknown)")
                 .str.replace("L3 composition", "L3 Composition")
                 .str.replace("L4 io signatures", "L4 IO signatures")
                 .str.replace("L5 nmf programs", "L5 NMF"))
null_lookup = null_df.set_index("target")["null_95"].to_dict()

fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True,
                           sharey=False)
colors = ["#DDCC77", "#AA4499", "#882255", "#332288",
          "#44AA99", "#117733", "#000000"]

for ax, metric, lo, hi, title in [
        (axes[0], "auc",          "auc_ci_lo", "auc_ci_hi",
         "ROC-AUC (LOOCV)"),
        (axes[1], "pr_auc",       "pr_ci_lo",  "pr_ci_hi",
         "Precision-Recall AUC"),
        (axes[2], "balanced_acc", "ba_ci_lo",  "ba_ci_hi",
         "Balanced accuracy")]:
    vals = met[metric].values
    los  = vals - met[lo].values
    his  = met[hi].values - vals
    bars = ax.bar(range(len(met)), vals, yerr=[los, his],
                   color=colors[:len(met)], edgecolor="black",
                   linewidth=0.8, capsize=4, width=0.7)
    ax.set_xticks(range(len(met)))
    ax.set_xticklabels(met["label"], rotation=35, ha="right", fontsize=9)
    ax.set_ylabel(title)
    ax.set_title(title, fontweight="bold", fontsize=11)
    if metric == "auc":
        ax.axhline(0.5, color="grey", linestyle=":", linewidth=0.8,
                    label="Chance (0.5)")
        first_null = True
        for i, layer in enumerate(met["layer"]):
            if layer not in null_lookup:
                continue
            ax.plot(
                i,
                null_lookup[layer],
                ls="",
                marker="_",
                markersize=25,
                markeredgewidth=2,
                color="red",
                label="95% permutation null" if first_null else None,
            )
            first_null = False
        ax.set_ylim(0, 1.05)
        ax.legend(loc="lower left", fontsize=7, frameon=False)
    elif metric == "pr_auc":
        prev = float(met["prevalence_baseline"].iloc[0])
        ax.axhline(prev, color="grey", linestyle=":", linewidth=0.8,
                    label=f"Prevalence ({prev:.2f})")
        ax.set_ylim(0, 1.05)
        ax.legend(loc="lower left", fontsize=7, frameon=False)
    else:
        ax.axhline(0.5, color="grey", linestyle=":", linewidth=0.8,
                    label="Chance (0.5)")
        ax.set_ylim(0, 1.05)
        ax.legend(loc="lower left", fontsize=7, frameon=False)
    # annotate value above bar
    for i, v in enumerate(vals):
        ax.text(i, v + 0.02, f"{v:.2f}", ha="center", va="bottom",
                 fontsize=7)

fig.suptitle("Balanced multi-layer LOOCV — ROC-AUC, PR-AUC and balanced "
              "accuracy with 2,000-iteration bootstrap 95 % CIs",
              fontweight="bold", fontsize=12, y=1.05)

for fmt in ("png", "pdf"):
    fig.savefig(FIG5 / f"Fig6A_balanced_layer_metrics.{fmt}",
                 dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved Fig6A_balanced_layer_metrics.png")

print("\n[05f DONE]")
