#!/usr/bin/env python3
"""
05l_cnv_exhaustion_spatial_mediation.py — T2.3

Chain: CNV clone-diversity  ->  T-cell exhaustion  ->  response
with a spatial-section consistency check.

Scope
-----
Per-patient mediation runs within discovery cohort (GSE206325, n=24) because
CNV and exhaustion features are both measured per-patient there.  The
spatial cohorts (GSE238264+GSE281759, 7 sections) provide an orthogonal
check at the section level: does mean NMF3 activity (already linked to
response in T1.3) also track exhaustion-programme loading?

Inputs
------
out/05_multilayer/tables/patient_multilayer_features.csv
out/05_multilayer/tables/nmf_programs_spatial_activity.csv
out/05_multilayer/tables/nmf_programs_gene_loadings.csv

Outputs
-------
out/05_multilayer/tables/cnv_exhaustion_mediation_summary.csv
out/05_multilayer/tables/cnv_exhaustion_mediation_bootstrap.csv
out/05_multilayer/tables/spatial_exhaustion_consistency.csv
out/05_multilayer/figures/05l_cnv_exhaustion_chain.{png,pdf}
"""
from __future__ import annotations
import sys, warnings; warnings.filterwarnings("ignore")
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
feat = feat[feat["dataset"] == "GSE206325"].reset_index(drop=True)
feat["y"] = (feat["response"] == "Responder").astype(int)
print(f"[05l] GSE206325 n={len(feat)}  R={feat['y'].sum()} NR={(1-feat['y']).sum()}")

# Exposure, mediator, outcome
X_raw = feat["cnv_clone_diversity"].values
M_raw = feat["sig_T_cell_exhaustion"].values
Y     = feat["y"].values

def _z(a):
    a = np.asarray(a, float); m = np.nanmean(a); s = np.nanstd(a)
    return (a - m) / (s if s > 0 else 1)

X = _z(X_raw); M = _z(M_raw)

# Step a: X -> M  (OLS)
aX = sm.add_constant(X)
fit_a = sm.OLS(M, aX).fit()
a = float(fit_a.params[1]); a_p = float(fit_a.pvalues[1])

# Step b: logistic Y ~ X + M  -> coefficient b is mediator given exposure
fitb = sm.Logit(Y, sm.add_constant(np.column_stack([X, M]))).fit(disp=False, maxiter=200)
c_prime = float(fitb.params[1])          # direct effect of X
b        = float(fitb.params[2])          # mediator effect
c_prime_p = float(fitb.pvalues[1])
b_p       = float(fitb.pvalues[2])

# Step c: total effect logistic Y ~ X
fitc = sm.Logit(Y, sm.add_constant(X)).fit(disp=False, maxiter=200)
c_total = float(fitc.params[isadeghi1]); c_total_p = float(fitc.pvalues[1])

indirect = a * b
prop_mediated = indirect / c_total if abs(c_total) > 1e-9 else np.nan

# Nonparametric bootstrap CI (1000 resamples)
B = 1000
rng = np.random.default_rng(0)
boots = []
for _ in range(B):
    idx = rng.integers(0, len(feat), len(feat))
    Xi = X[idx]; Mi = M[idx]; Yi = Y[idx]
    if Yi.sum() in (0, len(Yi)):
        continue
    try:
        fa = sm.OLS(Mi, sm.add_constant(Xi)).fit()
        ai = float(fa.params[1])
        fb = sm.Logit(Yi, sm.add_constant(np.column_stack([Xi, Mi]))).fit(
            disp=False, maxiter=200)
        bi = float(fb.params[2]); cpi = float(fb.params[1])
        fc = sm.Logit(Yi, sm.add_constant(Xi)).fit(disp=False, maxiter=200)
        ci = float(fc.params[1])
        boots.append((ai, bi, cpi, ci, ai * bi, ai * bi / ci if abs(ci) > 1e-9 else np.nan))
    except Exception:
        continue
boots = pd.DataFrame(boots, columns=["a", "b", "c_prime", "c_total",
                                      "indirect", "prop_mediated"])
boots.to_csv(TAB / "cnv_exhaustion_mediation_bootstrap.csv", index=False)

def ci(v): return (np.nanpercentile(v, 2.5), np.nanpercentile(v, 97.5))
summary = pd.DataFrame([
    dict(step="a (X -> M, OLS)", estimate=a, p=a_p,
         ci_lo=ci(boots["a"])[0], ci_hi=ci(boots["a"])[1]),
    dict(step="b (M -> Y | X, logit)", estimate=b, p=b_p,
         ci_lo=ci(boots["b"])[0], ci_hi=ci(boots["b"])[1]),
    dict(step="c' direct (X -> Y | M, logit)", estimate=c_prime, p=c_prime_p,
         ci_lo=ci(boots["c_prime"])[0], ci_hi=ci(boots["c_prime"])[1]),
    dict(step="c total (X -> Y, logit)", estimate=c_total, p=c_total_p,
         ci_lo=ci(boots["c_total"])[0], ci_hi=ci(boots["c_total"])[1]),
    dict(step="a*b indirect effect", estimate=indirect, p=np.nan,
         ci_lo=ci(boots["indirect"])[0], ci_hi=ci(boots["indirect"])[1]),
    dict(step="proportion mediated", estimate=prop_mediated, p=np.nan,
         ci_lo=ci(boots["prop_mediated"])[0], ci_hi=ci(boots["prop_mediated"])[1]),
])
summary.to_csv(TAB / "cnv_exhaustion_mediation_summary.csv", index=False)
print("\nMediation summary:")
print(summary.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

# ── Spatial consistency check ────────────────────────────────────────
sp = pd.read_csv(TAB / "nmf_programs_spatial_activity.csv")
sp["response"] = np.where(sp["section"].str.endswith("R") &
                           ~sp["section"].str.endswith("NR"),
                           "Responder", "Non-Responder")
pivot = sp.pivot_table(index="section", columns="program", values="mean")
pivot["response"] = pivot.index.to_series().apply(
    lambda s: "Responder" if (s.endswith("R") and not s.endswith("NR"))
    else "Non-Responder")

# Use NMF3 (linked to response) and, where present, a proxy exhaustion
# programme: NMF with top genes enriched in CD8-dysfunctional markers.
loads = pd.read_csv(TAB / "nmf_programs_gene_loadings.csv", index_col=0)
EXH_GENES = ["TOX", "PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "ENTPD1",
             "GZMB", "IFNG", "EOMES", "TNFRSF9", "ZEB2"]
exh_scores = {p: loads.loc[loads.index.isin(EXH_GENES), p].sum()
              for p in pivot.columns if p != "response"}
exh_prog = max(exh_scores, key=exh_scores.get)
print(f"\n   exhaustion-proxy programme = {exh_prog} "
      f"(sum of top exh-gene loadings = {exh_scores[exh_prog]:.3f})")

sp_rows = []
for prog in ["NMF3", exh_prog]:
    r_v = pivot.loc[pivot["response"] == "Responder", prog].astype(float)
    nr_v = pivot.loc[pivot["response"] == "Non-Responder", prog].astype(float)
    u, p = stats.mannwhitneyu(r_v, nr_v, alternative="two-sided")
    sp_rows.append(dict(programme=prog, r_mean=r_v.mean(),
                        nr_mean=nr_v.mean(), u=u, p_value=p,
                        n_R=len(r_v), n_NR=len(nr_v)))
# Also the between-programme correlation across sections
rho, prho = stats.spearmanr(pivot["NMF3"].astype(float),
                             pivot[exh_prog].astype(float))
sp_rows.append(dict(programme=f"NMF3↔{exh_prog} (section-level Spearman)",
                    r_mean=rho, nr_mean=np.nan, u=np.nan, p_value=prho,
                    n_R=len(pivot), n_NR=np.nan))
sp_df = pd.DataFrame(sp_rows)
sp_df.to_csv(TAB / "spatial_exhaustion_consistency.csv", index=False)
print("\nSpatial consistency:")
print(sp_df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

# ── Figure ──────────────────────────────────────────────────────────
fig = plt.figure(figsize=(12.5, 7.2), constrained_layout=True)
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1])

# A — scatter X vs M, coloured by response
axA = fig.add_subplot(gs[0, 0])
for lab, col in [("Responder", R_COL), ("Non-Responder", NR_COL)]:
    sub = feat[feat["response"] == lab]
    axA.scatter(sub["cnv_clone_diversity"], sub["sig_T_cell_exhaustion"],
                s=50, color=col, edgecolor="black", linewidth=0.4,
                alpha=0.85, label=f"{lab} (n={len(sub)})")
axA.set_xlabel("CNV clone diversity", fontsize=9)
axA.set_ylabel("T-cell exhaustion signature", fontsize=9)
rho_ab, p_ab = stats.pearsonr(feat["cnv_clone_diversity"],
                               feat["sig_T_cell_exhaustion"])
axA.set_title(f"A  Step a:  CNV → exhaustion\nPearson r={rho_ab:.2f}  p={p_ab:.3f}",
              loc="left", fontsize=10, fontweight="bold")
axA.legend(frameon=False, fontsize=8)

# B — violin exhaustion vs response
axB = fig.add_subplot(gs[0, 1])
data_v = [feat[feat["y"] == 1]["sig_T_cell_exhaustion"].values,
           feat[feat["y"] == 0]["sig_T_cell_exhaustion"].values]
parts = axB.violinplot(data_v, positions=[0, 1], widths=0.7, showmeans=False,
                        showmedians=True, showextrema=False)
for body, c in zip(parts["bodies"], [R_COL, NR_COL]):
    body.set_facecolor(c); body.set_edgecolor("black"); body.set_alpha(0.55)
for i, (d, c) in enumerate(zip(data_v, [R_COL, NR_COL])):
    axB.scatter(np.random.uniform(i-0.05, i+0.05, len(d)), d,
                s=30, color=c, edgecolor="black", linewidth=0.3, alpha=0.85)
u, p_by = stats.mannwhitneyu(data_v[0], data_v[1])
axB.set_xticks([0, 1]); axB.set_xticklabels(["R (n=7)", "NR (n=17)"], fontsize=9)
axB.set_ylabel("T-cell exhaustion signature", fontsize=9)
axB.set_title(f"B  Step b:  exhaustion → response\np={p_by:.3f}",
              loc="left", fontsize=10, fontweight="bold")

# C — mediation coefficient forest
axC = fig.add_subplot(gs[0, 2])
rows_c = summary[summary["step"].str.startswith(("a ", "b ", "c'", "c total",
                                                  "a*b"))]
y = np.arange(len(rows_c))
axC.errorbar(rows_c["estimate"], y,
             xerr=[rows_c["estimate"] - rows_c["ci_lo"],
                   rows_c["ci_hi"] - rows_c["estimate"]],
             fmt="o", color="black", ecolor="black", capsize=3)
axC.axvline(0, color="grey", lw=0.7, ls="--")
axC.set_yticks(y); axC.set_yticklabels(rows_c["step"], fontsize=8)
for i, (_, r) in enumerate(rows_c.iterrows()):
    if pd.notna(r["p"]):
        axC.text(r["ci_hi"] + 0.05, i, f"p={r['p']:.3f}",
                 va="center", fontsize=7.5,
                 fontweight="bold" if r["p"] < 0.05 else "normal")
axC.set_xlabel("Coefficient (standardised)", fontsize=9)
axC.set_title("C  Mediation path coefficients\n(95% bootstrap CI, B=1000)",
              loc="left", fontsize=10, fontweight="bold")

# D — bootstrap indirect effect distribution
axD = fig.add_subplot(gs[1, 0])
vals = boots["indirect"].dropna().values
axD.hist(vals, bins=40, color="#DDCC77", edgecolor="black", linewidth=0.3)
axD.axvline(indirect, color=NR_COL, lw=2, label=f"Observed a·b = {indirect:.2f}")
lo, hi = ci(vals)
axD.axvline(lo, color="black", lw=0.7, ls="--")
axD.axvline(hi, color="black", lw=0.7, ls="--")
axD.set_xlabel("Indirect effect a·b", fontsize=9)
axD.set_ylabel("Bootstrap frequency", fontsize=9)
axD.set_title(f"D  Indirect effect distribution\n95% CI [{lo:.2f}, {hi:.2f}]",
              loc="left", fontsize=10, fontweight="bold")
axD.legend(frameon=False, fontsize=8)

# E — spatial consistency bars
axE = fig.add_subplot(gs[1, 1])
bar_df = sp_df[~sp_df["programme"].str.contains("Spearman")].copy()
y = np.arange(len(bar_df))
width = 0.35
axE.barh(y - width/2, bar_df["r_mean"], width, color=R_COL,
         edgecolor="black", linewidth=0.4, label="Responder sections")
axE.barh(y + width/2, bar_df["nr_mean"], width, color=NR_COL,
         edgecolor="black", linewidth=0.4, label="Non-R sections")
for i, (_, r) in enumerate(bar_df.iterrows()):
    axE.text(max(r["r_mean"], r["nr_mean"]) * 1.02, i,
             f"p={r['p_value']:.3f}",
             va="center", fontsize=8,
             fontweight="bold" if r["p_value"] < 0.05 else "normal")
axE.set_yticks(y); axE.set_yticklabels(bar_df["programme"], fontsize=8.5)
axE.set_xlabel("Mean spatial programme activity", fontsize=9)
axE.set_title("E  Spatial consistency\n(NMF3 vs exhaustion-proxy programme)",
              loc="left", fontsize=10, fontweight="bold")
axE.legend(frameon=False, fontsize=8, loc="lower right")

# F — DAG-like summary
axF = fig.add_subplot(gs[1, 2]); axF.axis("off")
txt = "Mediation chain\n" + "─"*28 + "\n\n"
txt += f"X: CNV clone diversity\n"
txt += f"M: T-cell exhaustion sig\n"
txt += f"Y: Responder (logit)\n\n"
txt += f"a (X→M)      = {a:+.2f}  p={a_p:.3f}\n"
txt += f"b (M→Y|X)    = {b:+.2f}  p={b_p:.3f}\n"
txt += f"c total      = {c_total:+.2f}  p={c_total_p:.3f}\n"
txt += f"c' direct    = {c_prime:+.2f}  p={c_prime_p:.3f}\n"
txt += f"indirect a·b = {indirect:+.2f}\n"
if pd.notna(prop_mediated):
    txt += f"% mediated   = {prop_mediated*100:.0f}%\n"
axF.text(0.0, 0.95, txt, ha="left", va="top", family="monospace", fontsize=9)
axF.set_title("F  Summary", loc="left", fontsize=10, fontweight="bold")

fig.suptitle("Figure: CNV clone diversity → T-cell exhaustion → response mediation",
             fontsize=12, fontweight="bold")
fig.savefig(FIG / "05l_cnv_exhaustion_chain.png", dpi=600,
            bbox_inches="tight", facecolor="white")
fig.savefig(FIG / "05l_cnv_exhaustion_chain.pdf",
            bbox_inches="tight", facecolor="white")
print(f"\nwrote {FIG/'05l_cnv_exhaustion_chain.png'}")
