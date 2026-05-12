#!/usr/bin/env python3
"""
05i_nmf10_transfer_diagnostic.py
================================

Tier-2 (T2.2): diagnostic panels for the NMF10 plasma/TLS programme.

This script intentionally treats NMF10 as a composition-linked plasma/IG
programme. The key diagnostic is attenuation of the apparent response
association after controlling for B/plasma composition, not persistence of an
independent responder programme.

Inputs
------
- out/05_multilayer/tables/nmf_programs_patient_activity.csv
- out/05_multilayer/tables/nmf_programs_gene_loadings.csv
- out/05_multilayer/tables/patient_multilayer_features.csv (L3 composition)

Outputs
-------
- out/05_multilayer/tables/nmf10_gene_overlap.csv
- out/05_multilayer/tables/nmf10_partial_correlation.csv
- out/05_multilayer/figures/05i_nmf10_plasma_diagnostic.{png,pdf}
- out/05_multilayer/figures/05i_nmf10_transfer.{png,pdf}  # compatibility alias
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
from config.publication_style import apply_style, PAL
apply_style()
np.random.seed(0)

TAB = PROJECT / "out" / "05_multilayer" / "tables"
FIG = PROJECT / "out" / "05_multilayer" / "figures"
R_COL = PAL.get("Responder", "#0077BB")
NR_COL = PAL.get("Non_Responder", "#CC3311")

act = pd.read_csv(TAB / "nmf_programs_patient_activity.csv")
loads = pd.read_csv(TAB / "nmf_programs_gene_loadings.csv", index_col=0)

# ------------------------------------------------------------------
# 1. Gene-overlap of NMF10 top genes with canonical gene sets
# ------------------------------------------------------------------
top_n = 50
nmf10_top = loads["NMF10"].sort_values(ascending=False).head(top_n).index.tolist()

# Canonical panels
PLASMA_PANEL = ["IGKC", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2",
                "IGHM", "IGJ", "JCHAIN", "MZB1", "CD79A", "CD79B",
                "XBP1", "DERL3", "FKBP11", "SSR4", "PRDX4", "SEC11C",
                "ITM2C", "TNFRSF17"]
TLS_PANEL = ["CXCL13", "CCL19", "CCL21", "CXCR5", "CR2", "CD79A",
             "MS4A1", "LAMP3", "CD83", "CCR7", "IGHG1", "IGHG2",
             "CXCL12", "SELL", "TNFRSF13B"]
BCELL_PANEL = ["CD19", "MS4A1", "CD79A", "CD79B", "BANK1", "HLA-DRA",
               "CD22", "CD37", "PAX5"]
CD8_PANEL = ["CD8A", "CD8B", "GZMB", "GZMA", "GZMK", "PRF1", "NKG7",
             "IFNG", "CXCL13", "TOX", "EOMES"]
MREG_PANEL = ["LAMP3", "CCR7", "CCL19", "CCL22", "FSCN1", "BIRC3"]

def _overlap_test(top, panel, universe_n=loads.shape[0]):
    top = set(top)
    panel = set(panel) & set(loads.index)
    k = len(top & panel)
    K = len(panel)
    N = universe_n
    n = len(top)
    if K == 0 or k == 0:
        p = 1.0
        enr = 0.0
    else:
        # hypergeometric
        rv = stats.hypergeom(N, K, n)
        p = rv.sf(k - 1)
        expected = K * n / N
        enr = k / expected if expected > 0 else np.nan
    return k, K, p, enr

rows = []
for name, panel in [("Plasma/IG", PLASMA_PANEL),
                    ("TLS (Cabrita/Helmink)", TLS_PANEL),
                    ("Naive B-cell", BCELL_PANEL),
                    ("CD8 cytotoxic", CD8_PANEL),
                    ("mregDC", MREG_PANEL)]:
    k, K, p, enr = _overlap_test(nmf10_top, panel)
    rows.append(dict(panel=name, panel_size=K, top_hits=k,
                     fold_enrichment=enr, p_hypergeom=p,
                     overlap_genes=", ".join(sorted(set(nmf10_top) & set(panel)))))

overlap = pd.DataFrame(rows)
overlap.to_csv(TAB / "nmf10_gene_overlap.csv", index=False)
print("Gene-overlap:")
print(overlap.to_string(index=False))

# ------------------------------------------------------------------
# 2. Partial correlation: NMF10 vs response controlling for B/plasma features
# ------------------------------------------------------------------
feat = pd.read_csv(TAB / "patient_multilayer_features.csv")
feat = feat.merge(act[["patient_id", "NMF10"]], on="patient_id", how="left",
                  suffixes=("", "_act"))
# Find a B-plasma composition feature
bcell_candidates = [c for c in feat.columns
                    if "B_" in c or "Plasma" in c or "plasma" in c.lower()]
print(f"B/Plasma composition candidates: {bcell_candidates}")

def _partial_corr(x, y, z):
    """Partial correlation of x with y controlling for z (all 1D)."""
    d = pd.DataFrame({"x": x, "y": y, "z": z}).dropna()
    if len(d) < 4:
        return np.nan, np.nan
    rxz = stats.pearsonr(d["x"], d["z"]).statistic
    ryz = stats.pearsonr(d["y"], d["z"]).statistic
    rxy = stats.pearsonr(d["x"], d["y"]).statistic
    denom = np.sqrt(max(1e-12, (1 - rxz**2) * (1 - ryz**2)))
    r = (rxy - rxz * ryz) / denom
    # t-stat p approximation
    n = len(d)
    t = r * np.sqrt((n - 3) / max(1e-12, 1 - r**2))
    p = 2 * stats.t.sf(abs(t), n - 3)
    return float(r), float(p)

disc = feat[feat["dataset"] == "GSE206325"].copy()
disc["resp_bin"] = (disc["response"] == "Responder").astype(int)
pc_rows = []
# Unconditional
r_uncond, p_uncond = stats.pointbiserialr(disc["resp_bin"], disc["NMF10"])
pc_rows.append(dict(control_feature="(none)", r=float(r_uncond),
                    p=float(p_uncond), n=len(disc)))
for cand in bcell_candidates:
    if cand not in disc.columns: continue
    r, p = _partial_corr(disc["NMF10"], disc["resp_bin"], disc[cand])
    pc_rows.append(dict(control_feature=cand, r=r, p=p, n=disc.shape[0]))
pc = pd.DataFrame(pc_rows)
pc.to_csv(TAB / "nmf10_partial_correlation.csv", index=False)
print("\nPartial correlations (NMF10 ↔ response | control):")
print(pc.to_string(index=False))

# ------------------------------------------------------------------
# 3. Figure
# ------------------------------------------------------------------
fig = plt.figure(figsize=(12.5, 8.8), constrained_layout=True)
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0])

# A — NMF10 activity by response, split by cohort
axA = fig.add_subplot(gs[0, 0])
cohort_map = {"GSE206325": "Discovery", "Cappuyns2023": "Cappuyns"}
act["cohort"] = act["dataset"].map(cohort_map)
rows_scatter = []
for i, (coh_name, coh_df) in enumerate(act.groupby("cohort")):
    for resp, col in [("Responder", R_COL), ("Non-Responder", NR_COL),
                       ("Non_Responder", NR_COL)]:
        sub = coh_df[coh_df["response"] == resp]
        if sub.empty: continue
        x = i + (-0.2 if "Non" in resp else 0.2)
        axA.scatter(np.random.uniform(x - 0.08, x + 0.08, len(sub)),
                    sub["NMF10"], s=34, color=col, edgecolor="black",
                    linewidth=0.4, alpha=0.85,
                    label=resp if coh_name == "Discovery" else None)
axA.set_xticks([0, 1])
axA.set_xticklabels(["Discovery (GSE206325)\n7R/17NR",
                     "Cappuyns 2023\n19R/13NR"], fontsize=8.5)
axA.set_ylabel("NMF10 activity", fontsize=9)
axA.set_title("A  NMF10 activity by response × cohort",
              loc="left", fontsize=10.5, fontweight="bold")
# dedupe legend
handles, labels = axA.get_legend_handles_labels()
seen = {}
for h, l in zip(handles, labels): seen.setdefault(l, h)
axA.legend(seen.values(), seen.keys(), frameon=False, fontsize=8,
           loc="upper right")
# R vs NR p-values by cohort
for i, (coh_name, coh_df) in enumerate(act.groupby("cohort")):
    r_vals = coh_df[coh_df["response"] == "Responder"]["NMF10"]
    nr_vals = coh_df[coh_df["response"].isin(["Non-Responder", "Non_Responder"])]["NMF10"]
    if len(r_vals) and len(nr_vals):
        u, pv = stats.mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
        axA.text(i, act["NMF10"].max() * 1.05, f"p={pv:.3f}",
                 ha="center", fontsize=8, fontweight="bold")

# B — gene overlap bars (-log10 p)
axB = fig.add_subplot(gs[0, 1])
ov = overlap.sort_values("fold_enrichment")
y = np.arange(len(ov))
axB.barh(y, ov["fold_enrichment"].fillna(0),
         color=["#117733" if f > 1 else "#B0B0B0" for f in ov["fold_enrichment"].fillna(0)],
         edgecolor="black", linewidth=0.4)
for i, (_, r) in enumerate(ov.iterrows()):
    axB.text(r["fold_enrichment"] + 0.2, i,
             f"{r['top_hits']}/{r['panel_size']}"
             + (f" p={r['p_hypergeom']:.1e}" if r['p_hypergeom'] < 0.05 else ""),
             va="center", fontsize=8)
axB.set_yticks(y); axB.set_yticklabels(ov["panel"], fontsize=8.5)
axB.set_xlabel("Fold-enrichment (observed / expected)", fontsize=9)
axB.set_title("B  NMF10 top-50 genes — panel overlap",
              loc="left", fontsize=10.5, fontweight="bold")
axB.axvline(1, ls="--", color="black", lw=0.6)

# C — NMF10 loadings top 20
axC = fig.add_subplot(gs[0, 2])
top20 = loads["NMF10"].sort_values(ascending=False).head(20)[::-1]
colors_c = [R_COL if g in PLASMA_PANEL else "#777777" for g in top20.index]
axC.barh(np.arange(len(top20)), top20.values, color=colors_c,
         edgecolor="black", linewidth=0.4)
axC.set_yticks(np.arange(len(top20))); axC.set_yticklabels(top20.index, fontsize=8)
axC.set_xlabel("NMF10 gene loading", fontsize=9)
axC.set_title("C  NMF10 top genes (blue=plasma/IG panel)",
              loc="left", fontsize=10.5, fontweight="bold")

# D — partial correlations / attenuation diagnostic
axD = fig.add_subplot(gs[1, :2])
pcp = pc.copy().sort_values("r")
y = np.arange(len(pcp))
axD.barh(y, pcp["r"],
         color=[R_COL if r > 0 else NR_COL for r in pcp["r"]],
         edgecolor="black", linewidth=0.4)
for i, (_, r) in enumerate(pcp.iterrows()):
    if pd.notna(r["p"]):
        axD.text(r["r"] + (0.01 if r["r"] >= 0 else -0.01), i,
                 f"p={r['p']:.3f}",
                 va="center",
                 ha="left" if r["r"] >= 0 else "right",
                 fontsize=8,
                 fontweight=("bold" if r["p"] < 0.05 else "normal"))
axD.set_yticks(y)
axD.set_yticklabels([f"control = {c}" for c in pcp["control_feature"]],
                    fontsize=8.5)
axD.axvline(0, color="black", lw=0.7)
axD.set_xlabel("Point-biserial / partial correlation of NMF10 with response",
               fontsize=9)
axD.set_title("D  NMF10 response association attenuates after "
              "B/plasma control",
              loc="left", fontsize=10.5, fontweight="bold")

# E — NMF10 discovery vs Cappuyns summary
axE = fig.add_subplot(gs[1, 2])
tbl = []
for (coh, sub) in act.groupby("cohort"):
    r_vals = sub[sub["response"] == "Responder"]["NMF10"]
    nr_vals = sub[sub["response"].isin(["Non-Responder", "Non_Responder"])]["NMF10"]
    if len(r_vals) and len(nr_vals):
        u, p = stats.mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
        # Cliff's delta
        n1, n2 = len(r_vals), len(nr_vals)
        d = (2 * u / (n1 * n2)) - 1
        tbl.append(dict(cohort=coh, r_mean=r_vals.mean(),
                        nr_mean=nr_vals.mean(),
                        cliff_delta=d, p=p))
tbldf = pd.DataFrame(tbl)
axE.axis("off")
text = "Cohort comparison\n" + "─" * 34 + "\n"
for _, r in tbldf.iterrows():
    text += (f"{r['cohort']}\n"
             f"  R mean = {r['r_mean']:.3f}\n"
             f"  NR mean = {r['nr_mean']:.3f}\n"
             f"  Cliff δ = {r['cliff_delta']:+.2f}\n"
             f"  p = {r['p']:.3f}\n\n")
axE.text(0.0, 0.95, text, ha="left", va="top", family="monospace",
         fontsize=9)
axE.set_title("E  Cohort context summary", loc="left",
              fontsize=10.5, fontweight="bold")

fig.suptitle("Figure: NMF10 plasma/TLS programme diagnostic",
             fontsize=12, fontweight="bold")

primary_png = FIG / "05i_nmf10_plasma_diagnostic.png"
primary_pdf = FIG / "05i_nmf10_plasma_diagnostic.pdf"
legacy_png = FIG / "05i_nmf10_transfer.png"
legacy_pdf = FIG / "05i_nmf10_transfer.pdf"
for out_path in [primary_png, legacy_png]:
    fig.savefig(out_path, dpi=600, bbox_inches="tight", facecolor="white")
for out_path in [primary_pdf, legacy_pdf]:
    fig.savefig(out_path, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"\nwrote {primary_png}")
print(f"wrote {primary_pdf}")
print(f"compatibility aliases: {legacy_png}, {legacy_pdf}")
