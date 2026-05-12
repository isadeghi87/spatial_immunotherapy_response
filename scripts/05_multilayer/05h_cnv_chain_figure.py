#!/usr/bin/env python3
"""
05h_cnv_chain_figure.py
=======================

Tier-1 deliverable (T1.3): consolidate the CNV → pathway → NMF / response
mechanistic chain into a single composite figure and a summary table.

Inputs
------
- out/05_multilayer/tables/cnv_chromosome_arm_stats.csv
- out/05_multilayer/tables/cnv_bin_response_stats.csv
- out/05_multilayer/tables/cnv_gene_response_ranked.csv
- out/05_multilayer/tables/cnv_pathway_enrichment.csv
- out/05_multilayer/tables/cnv_patient_trajectory.csv
- out/05_multilayer/tables/cnv_response_stats.csv
- out/05_multilayer/tables/nmf_programs_response_stats.csv

Outputs
-------
- out/05_multilayer/tables/cnv_chain_summary.csv
- out/05_multilayer/figures/05h_cnv_chain.{png,pdf}

Figure layout (hero = chromosome forest panel, information-density scaled):
    [ A  wide chromosome forest (responder-gained vs lost) ]
    [ B  bin-level volcano   |  C  pathway enrichment      ]
    [ D  CNV metric R vs NR  |  E  clone-div vs NMF scatter ]
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
from config.publication_style import apply_style, PAL
apply_style()

TAB = PROJECT / "out" / "05_multilayer" / "tables"
FIG = PROJECT / "out" / "05_multilayer" / "figures"

R_COL = PAL.get("Responder", "#0077BB")
NR_COL = PAL.get("Non_Responder", "#CC3311")

arm = pd.read_csv(TAB / "cnv_chromosome_arm_stats.csv")
bins = pd.read_csv(TAB / "cnv_bin_response_stats.csv")
genes = pd.read_csv(TAB / "cnv_gene_response_ranked.csv")
path = pd.read_csv(TAB / "cnv_pathway_enrichment.csv")
traj = pd.read_csv(TAB / "cnv_patient_trajectory.csv")
cnv_stats = pd.read_csv(TAB / "cnv_response_stats.csv")
nmf_stats = pd.read_csv(TAB / "nmf_programs_response_stats.csv")

# ---- Summary table ---------------------------------------------
chrom_order = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
               "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
               "chr18","chr19","chr20","chr21","chr22"]
arm["chromosome"] = pd.Categorical(arm["chromosome"], categories=chrom_order,
                                   ordered=True)
arm = arm.sort_values("chromosome")

summary = pd.DataFrame({
    "step": ["CNV clone diversity (R vs NR)",
             "CNV score SD (R vs NR)",
             "Top responder-gained chromosomes (≥2% bins BH<0.10)",
             "Hallmark IFNγ-response enrichment in R-gained genes",
             "Hallmark TNFα/NFκB directional split",
             "CNV clone diversity ↔ NMF3 (Spearman)",
             "CNV clone diversity ↔ fraction stromal-hepatocytes"],
    "evidence": [
        f"δ={cnv_stats.loc[cnv_stats.feature=='cnv_clone_diversity','effect_size'].iloc[0]:.2f}, "
        f"p={cnv_stats.loc[cnv_stats.feature=='cnv_clone_diversity','p_value'].iloc[0]:.3f}",
        f"δ={cnv_stats.loc[cnv_stats.feature=='cnv_score_sd','effect_size'].iloc[0]:.2f}, "
        f"p={cnv_stats.loc[cnv_stats.feature=='cnv_score_sd','p_value'].iloc[0]:.3f}",
        ", ".join(arm.loc[arm["frac_bins_BH_lt_0.10"]>=0.02,"chromosome"].astype(str).tolist()),
        f"GBP1, STAT1 (OR=18.4, BH={path.loc[path.Term=='HALLMARK_INTERFERON_GAMMA_RESPONSE','Adjusted P-value'].iloc[0]:.3f})",
        "NFKBIA gained_in_R; NFKB1 lost_in_R (both BH=0.116)",
        f"rho={traj.loc[traj.feature=='NMF3','rho_vs_clone_diversity'].iloc[0]:.2f}, "
        f"p={traj.loc[traj.feature=='NMF3','p_spearman'].iloc[0]:.3f}",
        f"rho={traj.loc[traj.feature=='frac_Stromal_Hepatocytes','rho_vs_clone_diversity'].iloc[0]:.2f}, "
        f"p={traj.loc[traj.feature=='frac_Stromal_Hepatocytes','p_spearman'].iloc[0]:.4f}",
    ]
})
summary.to_csv(TAB / "cnv_chain_summary.csv", index=False)
print(f"wrote {TAB/'cnv_chain_summary.csv'}")

# ---- Figure ----------------------------------------------------
fig = plt.figure(figsize=(13.5, 9.8), constrained_layout=True)
gs = fig.add_gridspec(3, 3, height_ratios=[1.3, 1.0, 1.0],
                       width_ratios=[1.15, 1.0, 1.0])

# Panel A — chromosome forest (mean delta R-NR + BH fraction)
axA = fig.add_subplot(gs[0, :])
y = np.arange(len(arm))
colors = [R_COL if d > 0 else NR_COL for d in arm["mean_delta_R_minus_NR"]]
axA.barh(y, arm["mean_delta_R_minus_NR"], color=colors,
         edgecolor="black", linewidth=0.4, height=0.78)
# annotate BH-significant fraction
for i, (_, r) in enumerate(arm.iterrows()):
    if r["frac_bins_BH_lt_0.10"] > 0:
        axA.text(r["mean_delta_R_minus_NR"]
                 + (0.00015 if r["mean_delta_R_minus_NR"] >= 0 else -0.00015),
                 i, f"{r['frac_bins_BH_lt_0.10']*100:.1f}%",
                 ha="left" if r["mean_delta_R_minus_NR"] >= 0 else "right",
                 va="center", fontsize=7.5, fontweight="bold")
axA.set_yticks(y)
axA.set_yticklabels(arm["chromosome"].astype(str), fontsize=8.5)
axA.axvline(0, color="black", lw=0.7)
axA.set_xlabel("Mean CNV score Δ (R − NR)  •  % bins BH<0.10 annotated",
               fontsize=9)
axA.set_title("A  Chromosome-level CNV shift by response",
              loc="left", fontsize=11, fontweight="bold")
axA.invert_yaxis()
# legend
from matplotlib.patches import Patch
axA.legend(handles=[Patch(color=R_COL, label="Gained in R"),
                     Patch(color=NR_COL, label="Gained in NR")],
           frameon=False, loc="lower right", fontsize=8)

# Panel B — bin-level volcano
axB = fig.add_subplot(gs[1, 0])
bins["neglog_p"] = -np.log10(bins["p_BH"].clip(lower=1e-10))
sig = bins["p_BH"] < 0.10
axB.scatter(bins.loc[~sig, "delta_R_minus_NR"], bins.loc[~sig, "neglog_p"],
            s=6, c="lightgrey", edgecolor="none", alpha=0.6)
axB.scatter(bins.loc[sig & (bins["delta_R_minus_NR"] > 0), "delta_R_minus_NR"],
            bins.loc[sig & (bins["delta_R_minus_NR"] > 0), "neglog_p"],
            s=22, c=R_COL, edgecolor="black", linewidth=0.3,
            label="R-gained (BH<0.10)")
axB.scatter(bins.loc[sig & (bins["delta_R_minus_NR"] < 0), "delta_R_minus_NR"],
            bins.loc[sig & (bins["delta_R_minus_NR"] < 0), "neglog_p"],
            s=22, c=NR_COL, edgecolor="black", linewidth=0.3,
            label="NR-gained (BH<0.10)")
axB.axhline(-np.log10(0.10), ls="--", color="black", lw=0.6)
axB.set_xlabel("Δ CNV score (R − NR)", fontsize=9)
axB.set_ylabel("-log10(BH-adjusted p)", fontsize=9)
axB.set_title("B  Bin-level volcano",
              loc="left", fontsize=10.5, fontweight="bold")
axB.legend(frameon=False, fontsize=7.5, loc="upper left")

# Panel C — pathway enrichment
axC = fig.add_subplot(gs[1, 1])
path_plot = path.copy()
path_plot["label"] = path_plot["Term"].str.replace("HALLMARK_", "").str.replace("_", " ")
path_plot["dir_col"] = np.where(path_plot["gene_set_label"] == "gained_in_R", R_COL, NR_COL)
path_plot["neglog"] = -np.log10(path_plot["Adjusted P-value"].clip(lower=1e-10))
ord_idx = np.argsort(path_plot["neglog"].values)
axC.barh(np.arange(len(path_plot)),
         path_plot["neglog"].values[ord_idx],
         color=path_plot["dir_col"].values[ord_idx],
         edgecolor="black", linewidth=0.4)
axC.set_yticks(np.arange(len(path_plot)))
dir_short = pd.Series(path_plot["gene_set_label"].values[ord_idx]).replace(
    {"gained_in_R": "R↑", "lost_in_R": "R↓"}).values
labs = [f"{lab} ({d})" for lab, d in zip(path_plot["label"].values[ord_idx], dir_short)]
axC.set_yticklabels(labs, fontsize=8)
axC.axvline(-np.log10(0.10), ls="--", color="black", lw=0.6)
axC.set_xlabel("-log10(BH p)", fontsize=9)
axC.set_title("C  Hallmark pathway enrichment",
              loc="left", fontsize=10.5, fontweight="bold")

# Panel D — CNV metric R vs NR
axD = fig.add_subplot(gs[1, 2])
metrics_show = ["cnv_score_mean", "cnv_score_sd", "cnv_clone_diversity"]
labels = ["CNV score\nmean", "CNV score\nSD", "CNV clone\ndiversity"]
x = np.arange(len(metrics_show))
mean_R = [cnv_stats.loc[cnv_stats.feature == m, "mean_R"].iloc[0] for m in metrics_show]
mean_NR = [cnv_stats.loc[cnv_stats.feature == m, "mean_NR"].iloc[0] for m in metrics_show]
# Normalize to 0-1 for visual parity (each metric scaled by its own max)
mx = np.maximum(mean_R, mean_NR)
rn = np.array(mean_R) / mx
nrn = np.array(mean_NR) / mx
axD.bar(x - 0.2, rn, 0.4, color=R_COL, edgecolor="black", linewidth=0.4,
        label="R (n=7)")
axD.bar(x + 0.2, nrn, 0.4, color=NR_COL, edgecolor="black", linewidth=0.4,
        label="NR (n=17)")
for i, m in enumerate(metrics_show):
    p = cnv_stats.loc[cnv_stats.feature == m, "p_value"].iloc[0]
    mark = "***" if p < 0.01 else ("**" if p < 0.05 else ("†" if p < 0.1 else "ns"))
    axD.text(i, max(rn[i], nrn[i]) + 0.05, f"p={p:.3f}\n{mark}",
             ha="center", fontsize=7.5, fontweight="bold")
axD.set_xticks(x)
axD.set_xticklabels(labels, fontsize=8.5)
axD.set_ylabel("Normalised to per-metric max", fontsize=9)
axD.set_ylim(0, 1.35)
axD.set_title("D  CNV metric summary", loc="left",
              fontsize=10.5, fontweight="bold")
axD.legend(frameon=False, fontsize=7.5, loc="upper left")

# Panel E — clone-diversity vs top-correlated features
axE = fig.add_subplot(gs[2, :])
tt = traj.sort_values("rho_vs_clone_diversity", ascending=False)
tt = pd.concat([tt.head(8), tt.tail(4)]).drop_duplicates()
y = np.arange(len(tt))
colors_e = [R_COL if r > 0 else NR_COL for r in tt["rho_vs_clone_diversity"]]
axE.barh(y, tt["rho_vs_clone_diversity"], color=colors_e,
         edgecolor="black", linewidth=0.4)
axE.set_yticks(y)
axE.set_yticklabels(tt["feature"], fontsize=8.5)
axE.invert_yaxis()
axE.axvline(0, color="black", lw=0.7)
axE.set_xlabel("Spearman ρ vs CNV clone diversity (across 24 patients)",
               fontsize=9)
axE.set_title("E  CNV clone diversity couples to hepatocyte/stromal "
              "programmes, inversely to immune composition",
              loc="left", fontsize=10.5, fontweight="bold")
for i, (_, r) in enumerate(tt.iterrows()):
    if r["p_spearman"] < 0.05:
        axE.text(r["rho_vs_clone_diversity"]
                 + (0.02 if r["rho_vs_clone_diversity"] >= 0 else -0.02),
                 i, f"p={r['p_spearman']:.3f}",
                 va="center", fontsize=7.5,
                 ha="left" if r["rho_vs_clone_diversity"] >= 0 else "right",
                 fontweight="bold")

fig.suptitle(
    "Figure: CNV → pathway → TME-program mechanistic chain  "
    "(Tier 1 composite)",
    fontsize=12, fontweight="bold")

out_png = FIG / "05h_cnv_chain.png"
out_pdf = FIG / "05h_cnv_chain.pdf"
fig.savefig(out_png, dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"wrote {out_png}")
print(f"wrote {out_pdf}")
print("\nSummary:")
print(summary.to_string(index=False))
