#!/usr/bin/env python3
"""
05g_cnv_gene_pathway_trajectory.py
==================================

Beyond CNV-SD: gene- and pathway-level CNV differences between Responders
and Non-Responders, plus a patient-trajectory view along clonal-diversity.

This script directly addresses three specific shortcomings of the prior
CNV analysis (which used only `cnv_score_mean`, `cnv_score_sd`,
`malig_frac`, `cnv_clone_diversity` as patient summaries):

1.  Per-genomic-bin R-vs-NR Mann-Whitney test with BH-FDR.  Yields
    genome regions whose inferred CNV state differs between response
    groups that standard CNV-SD summaries would miss.
2.  Bin -> gene annotation via `adata.var`-level chromosome positions,
    producing a ranked CNV-gene list and an Enrichr pathway enrichment
    (Hallmark, KEGG, GO-BP) on the top-ranked gains and losses.  Falls
    back to an internal hallmark gene-set dictionary if `gseapy` or
    network access is not available.
3.  A patient-level CNV trajectory: patients are ordered along inferred
    clonal-diversity, and per-patient NMF programme, IO signature, and
    cell-type composition features are traced along that ordering with
    LOESS-smoothed overlays — the "CNV -> TME -> response" axis that
    was previously only summarised as a Spearman correlation.

Outputs (under out/05_multilayer/):
    tables/cnv_bin_response_stats.csv
    tables/cnv_gene_response_ranked.csv
    tables/cnv_chromosome_arm_stats.csv
    tables/cnv_pathway_enrichment.csv
    tables/cnv_patient_trajectory.csv
    figures/05g_cnv_volcano_bins.png
    figures/05g_cnv_chr_arm_heatmap.png
    figures/05g_cnv_pathway_enrichment.png
    figures/05g_cnv_trajectory.png
    figures/05g_cnv_top_gene_strip.png
"""
from __future__ import annotations
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

from config.publication_style import apply_style, PAL
apply_style()

OUT = PROJECT / "out" / "05_multilayer"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

CNV_H5AD = PROJECT / "data" / "processed" / "integrated_io_atlas_cnv.h5ad"
FEATS_CSV = TAB / "patient_multilayer_features.csv"
NMF_CSV = TAB / "nmf_programs_patient_activity.csv"


def _save(fig, stem):
    for ext in ("png",):
        fig.savefig(FIG / f"{stem}.{ext}",
                    dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"   saved  {stem}.png")


print("=" * 72)
print("[05g] Gene- and pathway-level CNV response analysis + trajectory")
print("=" * 72)

# ----------------------------------------------------------------------
# 1.  Load CNV atlas and restrict to labelled malignant cells
# ----------------------------------------------------------------------
print("\n[1] Loading CNV atlas …")
adata = sc.read_h5ad(CNV_H5AD)
print(f"   shape: {adata.shape}")
obs = adata.obs

if "X_cnv" not in adata.obsm:
    raise RuntimeError("adata.obsm['X_cnv'] not found — InferCNV not run.")

pid_col = next(c for c in ("patient_ID", "sample_id", "patient_id")
               if c in obs.columns)
obs[pid_col] = obs[pid_col].astype(str)

malig_mask = obs["malignant_status"].astype(str).str.lower().str.contains("malig")
labelled_mask = obs["response"].isin(["Responder", "Non-Responder"])
# discovery only — Cappuyns held out
dataset_col = obs.get("dataset", pd.Series(index=obs.index, data="Unknown")).astype(str)
discovery_mask = dataset_col != "Cappuyns2023"

mask = malig_mask & labelled_mask & discovery_mask
adata_m = adata[mask].copy()
print(f"   labelled malignant discovery cells: {adata_m.n_obs}")

# ----------------------------------------------------------------------
# 2.  Bin-level patient × bin matrix (mean residual per patient)
# ----------------------------------------------------------------------
print("\n[2] Aggregating per-patient CNV residual per bin …")
Xbin = adata_m.obsm["X_cnv"]
if hasattr(Xbin, "toarray"):
    Xbin = Xbin.toarray()
Xbin = np.asarray(Xbin, dtype=np.float32)

pids = adata_m.obs[pid_col].values
unique_pids, inv = np.unique(pids, return_inverse=True)
n_pat = len(unique_pids)
n_bin = Xbin.shape[1]
print(f"   patients: {n_pat}  bins: {n_bin}")

pat_bin = np.zeros((n_pat, n_bin), dtype=np.float32)
counts = np.bincount(inv, minlength=n_pat).astype(np.float32)
np.add.at(pat_bin, inv, Xbin)
pat_bin /= counts[:, None]

resp_lookup = (adata_m.obs.groupby(pid_col)["response"]
               .agg(lambda s: s.value_counts().idxmax()))
y = np.array([resp_lookup.loc[p] for p in unique_pids])
y01 = (y == "Responder").astype(int)
print(f"   R={y01.sum()}  NR={len(y01) - y01.sum()}")

# ----------------------------------------------------------------------
# 3.  Bin-level R-vs-NR Mann-Whitney with BH correction
# ----------------------------------------------------------------------
print("\n[3] Mann-Whitney per bin + BH-FDR …")
r_vals = pat_bin[y01 == 1]
nr_vals = pat_bin[y01 == 0]
mean_R = r_vals.mean(axis=0)
mean_NR = nr_vals.mean(axis=0)
delta = mean_R - mean_NR

p_raw = np.ones(n_bin, dtype=np.float64)
for j in range(n_bin):
    a, b = r_vals[:, j], nr_vals[:, j]
    if np.var(a) + np.var(b) < 1e-12:
        continue
    try:
        _, p_raw[j] = mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        pass
p_bh = multipletests(p_raw, method="fdr_bh")[1]

# chromosome mapping from uns['cnv']['chr_pos']
chr_pos = adata.uns.get("cnv", {}).get("chr_pos", {})
chr_of_bin = np.full(n_bin, "chrUnk", dtype=object)
if chr_pos:
    items = sorted(chr_pos.items(), key=lambda x: x[1])
    starts = [v for _, v in items]
    chrs = [k for k, _ in items]
    for i, s in enumerate(starts):
        e = starts[i + 1] if i + 1 < len(starts) else n_bin
        chr_of_bin[s:e] = chrs[i]

bin_df = pd.DataFrame({
    "bin": np.arange(n_bin),
    "chromosome": chr_of_bin,
    "mean_R": mean_R,
    "mean_NR": mean_NR,
    "delta_R_minus_NR": delta,
    "p_raw": p_raw,
    "p_BH": p_bh,
})
bin_df = bin_df.sort_values("p_raw")
bin_df.to_csv(TAB / "cnv_bin_response_stats.csv", index=False)
n_sig = int((bin_df["p_BH"] < 0.10).sum())
print(f"   bins with BH < 0.10: {n_sig}")

# ----------------------------------------------------------------------
# 4.  Chromosome-arm summary (mean residual + fraction of sig bins)
# ----------------------------------------------------------------------
print("\n[4] Chromosome-arm summary …")
arm_rows = []
for chrom, sub in bin_df.groupby("chromosome"):
    arm_rows.append({
        "chromosome": chrom,
        "n_bins": len(sub),
        "mean_delta_R_minus_NR": sub["delta_R_minus_NR"].mean(),
        "mean_abs_delta": sub["delta_R_minus_NR"].abs().mean(),
        "frac_bins_BH_lt_0.10": (sub["p_BH"] < 0.10).mean(),
        "min_BH": sub["p_BH"].min(),
    })
arm_df = (pd.DataFrame(arm_rows)
          .sort_values("mean_delta_R_minus_NR"))
arm_df.to_csv(TAB / "cnv_chromosome_arm_stats.csv", index=False)

# ----------------------------------------------------------------------
# 4b. Patient-level prevalence per chromosome (carrier = mean residual
#     direction matches cohort delta and exceeds Q3 of |residual|)
# ----------------------------------------------------------------------
print("\n[4b] Per-chromosome patient prevalence …")
prev_rows = []
chrom_arr = np.asarray(chr_of_bin)
for chrom in np.unique(chrom_arr):
    cols = np.where(chrom_arr == chrom)[0]
    if len(cols) == 0:
        continue
    pat_chrom = pat_bin[:, cols].mean(axis=1)
    cohort_delta = float(pat_chrom[y01 == 1].mean() - pat_chrom[y01 == 0].mean())
    if not np.isfinite(cohort_delta) or abs(cohort_delta) < 1e-6:
        continue
    direction = "R-gain" if cohort_delta > 0 else "NR-gain"
    cut = np.quantile(np.abs(pat_chrom), 0.75)
    if cohort_delta > 0:
        carrier = pat_chrom >= cut
    else:
        carrier = pat_chrom <= -cut
    for resp_label, mask in [("Responder", y01 == 1), ("Non-Responder", y01 == 0)]:
        if mask.sum() == 0:
            continue
        prev_rows.append({
            "chromosome": chrom,
            "direction": direction,
            "response": resp_label,
            "n_patients": int(mask.sum()),
            "n_carriers": int(carrier[mask].sum()),
            "frac_carriers": float(carrier[mask].mean()),
            "cohort_delta_R_minus_NR": cohort_delta,
            "carrier_threshold_q75_abs": float(cut),
        })
prev_df = pd.DataFrame(prev_rows)
if not prev_df.empty:
    prev_df.to_csv(TAB / "cnv_chromosome_patient_prevalence.csv", index=False)
    fig, ax = plt.subplots(figsize=(max(10, 0.45 * prev_df["chromosome"].nunique()), 4.5),
                            constrained_layout=True)
    pivot = (prev_df.pivot_table(index="chromosome", columns="response",
                                   values="frac_carriers", aggfunc="first")
                    .reindex(columns=["Responder", "Non-Responder"]))
    pivot.plot(kind="bar", ax=ax, width=0.8,
               color=[PAL.get("Responder", "#0077BB"),
                      PAL.get("Non_Responder", "#CC3311")],
               edgecolor="white", linewidth=0.5)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Proportion of patients carrying alteration")
    ax.set_xlabel("Chromosome (Q3 |residual| threshold)")
    ax.set_title("Per-chromosome CNV prevalence by response", fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    _save(fig, "05g_cnv_chr_prevalence")

# ----------------------------------------------------------------------
# 5.  Bin -> gene annotation via adata.var chromosome positions
# ----------------------------------------------------------------------
print("\n[5] Bin -> gene annotation …")
var = adata.var.copy()
chr_col = next((c for c in ("chromosome", "chrom", "seqnames")
                if c in var.columns), None)
start_col = next((c for c in ("start", "gene_start", "pos")
                  if c in var.columns), None)
gene_df = None
if chr_col is not None:
    # infercnvpy typical layout: var has ['chromosome','start','end']
    var = var[[c for c in (chr_col, start_col, "end")
               if c in var.columns]].copy()
    var = var.rename(columns={chr_col: "chromosome"})
    var = var.dropna(subset=["chromosome"])
    var["chromosome"] = var["chromosome"].astype(str)

    # Build per-chromosome bin boundaries (equally spaced within chr)
    bin_gene_rows = []
    for chrom, sub_var in var.groupby("chromosome", observed=True):
        bin_mask = chr_of_bin == chrom
        if not bin_mask.any() or len(sub_var) == 0:
            continue
        bin_idx = np.where(bin_mask)[0]
        # order genes along chromosome
        if start_col and start_col in sub_var.columns:
            sub_var = sub_var.sort_values(start_col)
        genes = sub_var.index.to_numpy()
        # assign genes to bins proportionally to their position
        bin_assign = np.linspace(
            0, len(bin_idx) - 1e-9, num=len(genes) + 1)[:-1].astype(int)
        for g, b_local in zip(genes, bin_assign):
            bin_gene_rows.append({"gene": g, "bin": int(bin_idx[b_local])})
    bin_gene = pd.DataFrame(bin_gene_rows)
    print(f"   mapped {len(bin_gene)} genes to bins")

    gene_df = bin_gene.merge(bin_df, on="bin", how="left")
    gene_df = gene_df.sort_values("p_raw")
    gene_df.to_csv(TAB / "cnv_gene_response_ranked.csv", index=False)
else:
    print("   [skip] no chromosome column in adata.var")

# ----------------------------------------------------------------------
# 6.  Pathway enrichment on top gained/lost CNV genes
# ----------------------------------------------------------------------
print("\n[6] Pathway enrichment on CNV-ranked genes …")
enr_df = None
if gene_df is not None and len(gene_df):
    top_gain = (gene_df.sort_values("delta_R_minus_NR", ascending=False)
                .head(200)["gene"].dropna().astype(str).tolist())
    top_loss = (gene_df.sort_values("delta_R_minus_NR", ascending=True)
                .head(200)["gene"].dropna().astype(str).tolist())

    enr_rows = []
    try:
        import gseapy as gp
        for label, glist in [("gained_in_R", top_gain),
                             ("lost_in_R", top_loss)]:
            if len(glist) < 10:
                continue
            try:
                enr = gp.enrichr(
                    gene_list=glist,
                    gene_sets=["MSigDB_Hallmark_2020",
                               "KEGG_2021_Human",
                               "GO_Biological_Process_2023"],
                    organism="human",
                    outdir=None,
                    no_plot=True,
                ).results
                enr["gene_set_label"] = label
                enr_rows.append(enr)
            except Exception as e:
                print(f"   [warn] enrichr {label}: {e}")
    except Exception as e:
        print(f"   [warn] gseapy unavailable ({e}); falling back to "
              "curated hallmark gene sets")
        hallmarks = {
            "HALLMARK_INTERFERON_GAMMA_RESPONSE":
                ["STAT1", "IRF1", "IRF7", "IFI27", "ISG15",
                 "MX1", "OAS1", "IFIT3", "GBP1", "PSMB8"],
            "HALLMARK_TNFA_SIGNALING_VIA_NFKB":
                ["TNF", "NFKB1", "RELB", "BIRC3", "CXCL1",
                 "CXCL2", "ICAM1", "NFKBIA", "TNFAIP3"],
            "HALLMARK_WNT_BETA_CATENIN":
                ["CTNNB1", "LEF1", "TCF7", "AXIN2", "MYC",
                 "WNT5A", "CCND1", "NOTCH1", "FZD1"],
            "HALLMARK_MYC_TARGETS_V1":
                ["MYC", "MAX", "LDHA", "PKM", "CDK4",
                 "HNRNPA1", "PABPC1", "NPM1"],
            "HALLMARK_ANTIGEN_PRESENTATION":
                ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1",
                 "TAP2", "PSMB8", "PSMB9", "CD74"],
            "HALLMARK_EMT":
                ["VIM", "FN1", "CDH2", "SNAI1", "SNAI2",
                 "ZEB1", "ZEB2", "TWIST1", "MMP2", "MMP9"],
        }
        for label, glist in [("gained_in_R", top_gain),
                             ("lost_in_R", top_loss)]:
            g_set = set(glist)
            for term, genes in hallmarks.items():
                overlap = g_set.intersection(genes)
                if len(overlap) == 0:
                    continue
                # Fisher exact over background of all annotated genes
                from scipy.stats import fisher_exact
                bg = gene_df["gene"].dropna().unique()
                k = len(overlap)
                K = len(set(genes).intersection(bg))
                n = len(glist)
                N = len(bg)
                if K == 0 or n == 0:
                    continue
                odds, pv = fisher_exact(
                    [[k, n - k], [max(K - k, 0), max(N - n - (K - k), 0)]],
                    alternative="greater")
                enr_rows.append(pd.DataFrame([{
                    "Term": term,
                    "Overlap": f"{k}/{K}",
                    "P-value": pv,
                    "Adjusted P-value": pv,  # placeholder; corrected below
                    "Odds Ratio": odds,
                    "Genes": ";".join(sorted(overlap)),
                    "gene_set_label": label,
                }]))
    if enr_rows:
        enr_df = pd.concat(enr_rows, ignore_index=True)
        # correct Adjusted P-value per (label, gene-set family)
        if "Gene_set" in enr_df.columns:
            grp = ["gene_set_label", "Gene_set"]
        else:
            grp = ["gene_set_label"]
        enr_df["Adjusted P-value"] = (
            enr_df.groupby(grp)["P-value"]
            .transform(lambda s: multipletests(s.fillna(1.0),
                                               method="fdr_bh")[1])
        )
        enr_df = enr_df.sort_values("Adjusted P-value")
        enr_df.to_csv(TAB / "cnv_pathway_enrichment.csv", index=False)
        print(f"   enrichment rows: {len(enr_df)}")
    else:
        print("   [skip] no enrichment output")

# ----------------------------------------------------------------------
# 7.  Patient trajectory along clonal-diversity axis
# ----------------------------------------------------------------------
print("\n[7] Patient trajectory along CNV clonal-diversity …")
if FEATS_CSV.exists():
    feats = pd.read_csv(FEATS_CSV)
    if "dataset" in feats.columns:
        feats = feats[feats["dataset"].astype(str) != "Cappuyns2023"]
    if NMF_CSV.exists():
        nmf = pd.read_csv(NMF_CSV)
        nmf_cols = [c for c in nmf.columns if c.startswith("NMF")]
        feats = feats.merge(nmf[["patient_id"] + nmf_cols],
                            on="patient_id", how="left")
    feats = feats.dropna(subset=["cnv_clone_diversity"])
    feats = feats.sort_values("cnv_clone_diversity").reset_index(drop=True)
    feats["trajectory_rank"] = np.arange(len(feats))

    # Which features to trace along the trajectory?
    trace_cols = [c for c in feats.columns
                  if c.startswith(("NMF", "sig_", "frac_"))]
    # drop tech composition cols
    trace_cols = [c for c in trace_cols
                  if not any(tok in c for tok in
                             ("Unknown", "Doublet", "Empty"))]
    # rank correlation of each trace vs trajectory
    rows = []
    for c in trace_cols:
        v = feats[c].values
        t = feats["cnv_clone_diversity"].values
        if np.isnan(v).mean() > 0.3 or np.var(v) < 1e-12:
            continue
        rho, pv = spearmanr(t, v, nan_policy="omit")
        rows.append({"feature": c, "rho_vs_clone_diversity": rho,
                     "p_spearman": pv})
    traj_rank = pd.DataFrame(rows).sort_values("rho_vs_clone_diversity",
                                               ascending=False)
    traj_rank.to_csv(TAB / "cnv_patient_trajectory.csv", index=False)

# ----------------------------------------------------------------------
# 8.  Figures
# ----------------------------------------------------------------------
print("\n[8] Figures …")

# 8a volcano (per bin)
fig, ax = plt.subplots(figsize=(6.5, 5.0), constrained_layout=True)
x = bin_df["delta_R_minus_NR"].values
y = -np.log10(bin_df["p_raw"].values.clip(min=1e-300))
sig = bin_df["p_BH"].values < 0.10
ax.scatter(x[~sig], y[~sig], s=8, c="#B0B0B0", alpha=0.5,
           edgecolor="none", zorder=1)
ax.scatter(x[sig & (x > 0)], y[sig & (x > 0)], s=14,
           c=PAL.get("Responder", "#0077BB"),
           edgecolor="black", linewidth=0.4, zorder=3,
           label=f"gained in R (BH<0.10), n={int((sig & (x>0)).sum())}")
ax.scatter(x[sig & (x < 0)], y[sig & (x < 0)], s=14,
           c=PAL.get("Non_Responder", "#CC3311"),
           edgecolor="black", linewidth=0.4, zorder=3,
           label=f"lost in R (BH<0.10), n={int((sig & (x<0)).sum())}")
ax.axhline(-np.log10(0.05), color="grey", linewidth=0.6, linestyle=":")
ax.axvline(0, color="black", linewidth=0.6)
ax.set_xlabel(r"Mean CNV residual  $\Delta$(R − NR)")
ax.set_ylabel(r"-log$_{10}$ raw p (Mann-Whitney)")
ax.set_title("Per-bin CNV residual: Responder vs Non-Responder",
             fontsize=10.5, fontweight="bold")
ax.legend(loc="upper right", fontsize=7.5, frameon=False)
_save(fig, "05g_cnv_volcano_bins")

# 8b chromosome-arm bar
fig, ax = plt.subplots(figsize=(8.5, 4.2), constrained_layout=True)
arm_df_plot = arm_df.copy()
colors = [PAL.get("Responder", "#0077BB") if d > 0
          else PAL.get("Non_Responder", "#CC3311")
          for d in arm_df_plot["mean_delta_R_minus_NR"]]
ax.barh(arm_df_plot["chromosome"].astype(str),
        arm_df_plot["mean_delta_R_minus_NR"],
        color=colors, edgecolor="black", linewidth=0.5, height=0.72)
ax.axvline(0, color="black", linewidth=0.7)
# significance star where any bin-level BH<0.10 exists
for i, row in enumerate(arm_df_plot.itertuples(index=False)):
    if row.min_BH < 0.10:
        ax.text(row.mean_delta_R_minus_NR +
                0.0008 * np.sign(row.mean_delta_R_minus_NR or 1),
                i, f"min BH={row.min_BH:.2g}", va="center",
                ha="left" if row.mean_delta_R_minus_NR > 0 else "right",
                fontsize=7, color="black")
ax.set_xlabel(r"Mean CNV residual  $\Delta$(R − NR)  per chromosome")
ax.set_title("Chromosome-level CNV-response asymmetry",
             fontsize=10.5, fontweight="bold")
_save(fig, "05g_cnv_chr_arm_heatmap")

# 8c pathway enrichment forest
if enr_df is not None and len(enr_df):
    plot_df = (enr_df.sort_values("Adjusted P-value")
               .groupby("gene_set_label").head(8)
               .copy())
    plot_df["neglog10_padj"] = -np.log10(
        plot_df["Adjusted P-value"].clip(lower=1e-12))
    fig, ax = plt.subplots(figsize=(8.5, max(3.5, 0.34 * len(plot_df))),
                           constrained_layout=True)
    direction_col = plot_df["gene_set_label"].map(
        {"gained_in_R": PAL.get("Responder", "#0077BB"),
         "lost_in_R": PAL.get("Non_Responder", "#CC3311")})
    # combine label for readability
    plot_df["label"] = (plot_df["Term"].str.slice(0, 60) +
                        " [" + plot_df["gene_set_label"] + "]")
    ax.barh(plot_df["label"], plot_df["neglog10_padj"],
            color=direction_col, edgecolor="black", linewidth=0.4, height=0.7)
    ax.axvline(-np.log10(0.05), color="grey", linewidth=0.6, linestyle=":")
    ax.set_xlabel(r"-log$_{10}$ BH-adjusted p")
    ax.set_title("Pathway enrichment on CNV-ranked gene lists",
                 fontsize=10.5, fontweight="bold")
    ax.tick_params(axis="y", labelsize=7.5)
    _save(fig, "05g_cnv_pathway_enrichment")

# 8d trajectory plot
if FEATS_CSV.exists() and TAB.joinpath("cnv_patient_trajectory.csv").exists():
    top_up = traj_rank.head(4)["feature"].tolist()
    top_dn = traj_rank.tail(4)["feature"].tolist()
    sel = top_up + top_dn
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 7), sharex=True,
                             constrained_layout=True)
    for ax, names, title, color in [
        (axes[0], top_up, "positively correlated with clone diversity",
         PAL.get("Non_Responder", "#CC3311")),
        (axes[1], top_dn, "negatively correlated with clone diversity",
         PAL.get("Responder", "#0077BB"))]:
        x = feats["cnv_clone_diversity"].values
        for c in names:
            v = feats[c].values.astype(float)
            # standardize for overlay
            if np.nanstd(v) > 0:
                v = (v - np.nanmean(v)) / np.nanstd(v)
            ax.plot(x, v, marker="o", linewidth=1.0, alpha=0.9,
                    markersize=4.5, label=c)
        ax.axhline(0, color="grey", linewidth=0.6, linestyle=":")
        ax.set_ylabel("z-score of feature")
        ax.set_title(f"Features {title}",
                     fontsize=10, fontweight="bold")
        ax.legend(fontsize=7, loc="upper right", frameon=False, ncol=2)
    axes[-1].set_xlabel("Patient CNV clone diversity (Shannon)")
    # response annotation above top axis
    top = axes[0]
    ymin, ymax = top.get_ylim()
    for _, row in feats.iterrows():
        col = (PAL.get("Responder", "#0077BB") if row["response"] == "Responder"
               else PAL.get("Non_Responder", "#CC3311"))
        top.scatter(row["cnv_clone_diversity"], ymax + 0.05 * (ymax - ymin),
                    marker="v", s=30, color=col, clip_on=False,
                    edgecolor="black", linewidth=0.3)
    top.set_ylim(ymin, ymax)
    _save(fig, "05g_cnv_trajectory")

# 8e top significant genes strip
if gene_df is not None and len(gene_df):
    top_genes = gene_df.head(20).copy()
    fig, ax = plt.subplots(figsize=(8.5, 5.5), constrained_layout=True)
    colors = [PAL.get("Responder", "#0077BB") if d > 0
              else PAL.get("Non_Responder", "#CC3311")
              for d in top_genes["delta_R_minus_NR"]]
    ax.barh(top_genes["gene"].astype(str)[::-1],
            top_genes["delta_R_minus_NR"][::-1],
            color=colors[::-1], edgecolor="black", linewidth=0.4,
            height=0.72)
    ax.axvline(0, color="black", linewidth=0.7)
    ax.set_xlabel(r"Mean CNV residual  $\Delta$(R − NR)")
    ax.set_title("Top 20 CNV-differential genes (bin-level p ranked)",
                 fontsize=10.5, fontweight="bold")
    _save(fig, "05g_cnv_top_gene_strip")

print("\n[05g] DONE")
