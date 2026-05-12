#!/usr/bin/env python3
"""
05k_nmf10_cappuyns_transfer.py — T2.2 extension

Port the NMF10 partial-correlation diagnostic to the Cappuyns cohort and
pool both cohorts with a harmonised B/plasma-total control feature.

Harmonisation
-------------
Discovery GSE206325 features `frac_B_Plasma + frac_B_Memory + frac_B_Naive`.
Cappuyns 2023 features `frac_B_cell + frac_Plasma` (bulk deconvolution).
Common control: `total_B_plasma = sum of the available features per cohort`.

Outputs
-------
out/05_multilayer/tables/nmf10_partial_correlation_cappuyns.csv
out/05_multilayer/tables/nmf10_partial_correlation_pooled.csv
out/05_multilayer/figures/05k_nmf10_cappuyns_transfer.{png,pdf}
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
from config.publication_style import apply_style, PAL
apply_style()

TAB = PROJECT / "out" / "05_multilayer" / "tables"
FIG = PROJECT / "out" / "05_multilayer" / "figures"
R_COL = PAL.get("Responder", "#0077BB")
NR_COL = PAL.get("Non_Responder", "#CC3311")

feat = pd.read_csv(TAB / "patient_multilayer_features.csv")
act = pd.read_csv(TAB / "nmf_programs_patient_activity.csv")[
    ["patient_id", "NMF10"]]
df = feat.merge(act, on="patient_id", how="left")
df["y"] = df["response"].isin(["Responder"]).astype(int)

# Harmonised total B/plasma
disc_cols = ["frac_B_Plasma", "frac_B_Memory", "frac_B_Naive"]
cap_cols  = ["frac_B_cell", "frac_Plasma"]
df["total_B_plasma"] = np.where(
    df["dataset"] == "GSE206325",
    df[disc_cols].sum(axis=1),
    df[cap_cols].sum(axis=1))

def partial_corr(x, y, z):
    d = pd.DataFrame({"x": x, "y": y, "z": z}).dropna()
    if len(d) < 4 or d["z"].std() == 0:
        return np.nan, np.nan, len(d)
    rxz = stats.pearsonr(d["x"], d["z"]).statistic
    ryz = stats.pearsonr(d["y"], d["z"]).statistic
    rxy = stats.pearsonr(d["x"], d["y"]).statistic
    denom = np.sqrt(max(1e-12, (1 - rxz**2) * (1 - ryz**2)))
    r = (rxy - rxz * ryz) / denom
    n = len(d)
    if abs(r) >= 1: return r, 0.0, n
    t = r * np.sqrt((n - 3) / max(1e-12, 1 - r**2))
    p = 2 * stats.t.sf(abs(t), n - 3)
    return float(r), float(p), int(n)

rows = []
for coh_name, coh_df in [("GSE206325 (discovery)", df[df["dataset"] == "GSE206325"]),
                          ("Cappuyns2023 (validation)", df[df["dataset"] == "Cappuyns2023"]),
                          ("Pooled (both cohorts)", df)]:
    # marginal point-biserial
    r_m, p_m = stats.pointbiserialr(coh_df["y"], coh_df["NMF10"])
    rows.append(dict(cohort=coh_name, control="(none / marginal)",
                     r=float(r_m), p=float(p_m), n=len(coh_df)))
    # partial given total_B_plasma
    r_p, p_p, n_p = partial_corr(coh_df["NMF10"], coh_df["y"],
                                  coh_df["total_B_plasma"])
    rows.append(dict(cohort=coh_name, control="total_B_plasma",
                     r=r_p, p=p_p, n=n_p))
out = pd.DataFrame(rows)
out.to_csv(TAB / "nmf10_partial_correlation_pooled.csv", index=False)
print(out.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

# Cappuyns-only table
cap_only = out[out["cohort"].str.startswith("Cappuyns")].copy()
cap_only.to_csv(TAB / "nmf10_partial_correlation_cappuyns.csv", index=False)

# Logistic regression: y ~ NMF10 + total_B_plasma  (pooled, with dataset fixed effect)
pooled = df.dropna(subset=["NMF10", "total_B_plasma"]).copy()
pooled["dataset_cap"] = (pooled["dataset"] == "Cappuyns2023").astype(int)
logit_rows = []
for name, cols in [("NMF10 only", ["NMF10"]),
                   ("NMF10 + total_B_plasma", ["NMF10", "total_B_plasma"]),
                   ("NMF10 + total_B_plasma + cohort",
                    ["NMF10", "total_B_plasma", "dataset_cap"])]:
    X = sm.add_constant(pooled[cols], has_constant="add")
    try:
        fit = sm.Logit(pooled["y"], X).fit(disp=False, maxiter=200)
        for c in cols:
            logit_rows.append(dict(model=name, term=c,
                                   coef=float(fit.params[c]),
                                   p=float(fit.pvalues[c]),
                                   aic=float(fit.aic), n=len(pooled)))
    except Exception as e:
        logit_rows.append(dict(model=name, term="FAIL", coef=np.nan,
                               p=np.nan, aic=np.nan, n=len(pooled)))
logit_df = pd.DataFrame(logit_rows)
logit_df.to_csv(TAB / "nmf10_pooled_logit_models.csv", index=False)
print("\nPooled logistic models:")
print(logit_df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

# ── Figure ──────────────────────────────────────────────────────
fig = plt.figure(figsize=(11.5, 7.5), constrained_layout=True)
gs = fig.add_gridspec(2, 3)

# A — NMF10 vs total_B_plasma scatter, both cohorts
axA = fig.add_subplot(gs[0, 0])
markers = {"GSE206325": "o", "Cappuyns2023": "s"}
for ds, sub in df.groupby("dataset"):
    for resp, col in [("Responder", R_COL), ("Non-Responder", NR_COL),
                       ("Non_Responder", NR_COL)]:
        s = sub[sub["response"] == resp]
        if s.empty: continue
        axA.scatter(s["total_B_plasma"], s["NMF10"], s=42,
                    marker=markers.get(ds, "o"),
                    color=col, edgecolor="black", linewidth=0.4,
                    alpha=0.85,
                    label=f"{ds} · {resp}" if resp in ("Responder", "Non-Responder") else None)
axA.set_xlabel("Total B/plasma fraction (harmonised)", fontsize=9)
axA.set_ylabel("NMF10 activity", fontsize=9)
axA.set_title("A  NMF10 vs composition\n(all 56 patients)", loc="left",
              fontsize=10.5, fontweight="bold")
axA.legend(frameon=False, fontsize=7, loc="upper left")
rho, pv = stats.spearmanr(df["total_B_plasma"], df["NMF10"])
axA.text(0.98, 0.02, f"Spearman ρ={rho:.2f}  p={pv:.1e}",
         transform=axA.transAxes, ha="right", va="bottom", fontsize=8,
         fontweight="bold")

# B — Cappuyns NMF10 by response
axB = fig.add_subplot(gs[0, 1])
cap = df[df["dataset"] == "Cappuyns2023"].copy()
for i, resp in enumerate(["Responder", "Non-Responder"]):
    s = cap[cap["response"] == resp]
    if s.empty:
        s = cap[cap["response"] == resp.replace("-", "_")]
    axB.scatter(np.random.uniform(i-0.1, i+0.1, len(s)), s["NMF10"],
                s=45, color=R_COL if resp == "Responder" else NR_COL,
                edgecolor="black", linewidth=0.4, alpha=0.85)
r_vals = cap[cap["response"] == "Responder"]["NMF10"]
nr_vals = cap[cap["response"].isin(["Non-Responder", "Non_Responder"])]["NMF10"]
if len(r_vals) and len(nr_vals):
    u, p = stats.mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
    axB.text(0.5, cap["NMF10"].max() * 1.03, f"p={p:.3f}",
             ha="center", transform=axB.transData, fontsize=9, fontweight="bold")
axB.set_xticks([0, 1]); axB.set_xticklabels(["R (n=19)", "NR (n=13)"], fontsize=9)
axB.set_ylabel("NMF10 activity", fontsize=9)
axB.set_title("B  Cappuyns 2023", loc="left", fontsize=10.5, fontweight="bold")

# C — Partial-correlation bars across all three cohorts
axC = fig.add_subplot(gs[0, 2])
plot_df = out.copy()
plot_df["label"] = plot_df["cohort"] + "\n" + plot_df["control"]
plot_df = plot_df.sort_values("r")
y = np.arange(len(plot_df))
axC.barh(y, plot_df["r"].fillna(0),
         color=[R_COL if r > 0 else NR_COL for r in plot_df["r"].fillna(0)],
         edgecolor="black", linewidth=0.4)
for i, (_, r) in enumerate(plot_df.iterrows()):
    if pd.notna(r["p"]):
        axC.text(r["r"] + (0.01 if r["r"] >= 0 else -0.01), i,
                 f"p={r['p']:.3f}",
                 va="center",
                 ha="left" if r["r"] >= 0 else "right",
                 fontsize=7,
                 fontweight=("bold" if r["p"] < 0.05 else "normal"))
axC.set_yticks(y); axC.set_yticklabels(plot_df["label"], fontsize=7.5)
axC.axvline(0, color="black", lw=0.7)
axC.set_xlabel("r  (NMF10 vs response)", fontsize=9)
axC.set_title("C  Partial-correlation ladder", loc="left",
              fontsize=10.5, fontweight="bold")

# D — coefficients from pooled logit
axD = fig.add_subplot(gs[1, :2])
ld = logit_df.dropna(subset=["coef"]).copy()
ld["label"] = ld["model"] + " · " + ld["term"]
y = np.arange(len(ld))
axD.barh(y, ld["coef"],
         color=[R_COL if v > 0 else NR_COL for v in ld["coef"]],
         edgecolor="black", linewidth=0.4)
for i, (_, r) in enumerate(ld.iterrows()):
    axD.text(r["coef"] + (0.1 if r["coef"] >= 0 else -0.1), i,
             f"p={r['p']:.3f} AIC={r['aic']:.1f}",
             va="center",
             ha="left" if r["coef"] >= 0 else "right",
             fontsize=7.5,
             fontweight=("bold" if r["p"] < 0.05 else "normal"))
axD.set_yticks(y); axD.set_yticklabels(ld["label"], fontsize=8)
axD.axvline(0, color="black", lw=0.7)
axD.set_xlabel("Logistic coefficient (pooled cohorts, n=56)", fontsize=9)
axD.set_title("D  Pooled logistic models — response ~ NMF10 + composition + cohort",
              loc="left", fontsize=10.5, fontweight="bold")

# E — summary text
axE = fig.add_subplot(gs[1, 2]); axE.axis("off")
txt = "Transfer-diagnostic verdict\n" + "─"*32 + "\n\n"
txt += "NMF10 = plasma/IG programme\n"
txt += " ─ 18/19 canonical IG genes\n\n"
txt += "Discovery (n=24):\n"
rr = out.query("cohort=='GSE206325 (discovery)'")
txt += f"  marginal r={rr.iloc[0]['r']:.2f} p={rr.iloc[0]['p']:.1e}\n"
txt += f"  | total_B r={rr.iloc[1]['r']:.2f} p={rr.iloc[1]['p']:.2f}\n\n"
txt += "Cappuyns (n=32):\n"
rc = out.query("cohort=='Cappuyns2023 (validation)'")
txt += f"  marginal r={rc.iloc[0]['r']:.2f} p={rc.iloc[0]['p']:.2f}\n"
txt += f"  | total_B r={rc.iloc[1]['r']:.2f} p={rc.iloc[1]['p']:.2f}\n\n"
txt += "Pooled (n=56):\n"
rp = out.query("cohort=='Pooled (both cohorts)'")
txt += f"  marginal r={rp.iloc[0]['r']:.2f} p={rp.iloc[0]['p']:.2e}\n"
txt += f"  | total_B r={rp.iloc[1]['r']:.2f} p={rp.iloc[1]['p']:.2f}\n"
axE.text(0.0, 0.98, txt, ha="left", va="top", family="monospace", fontsize=8.5)
axE.set_title("E  Summary", loc="left", fontsize=10.5, fontweight="bold")

fig.suptitle("Figure: NMF10 transfer diagnostic — discovery + Cappuyns + pooled",
             fontsize=12, fontweight="bold")
fig.savefig(FIG / "05k_nmf10_cappuyns_transfer.png", dpi=600,
            bbox_inches="tight", facecolor="white")
fig.savefig(FIG / "05k_nmf10_cappuyns_transfer.pdf",
            bbox_inches="tight", facecolor="white")
print(f"\nwrote {FIG/'05k_nmf10_cappuyns_transfer.png'}")
