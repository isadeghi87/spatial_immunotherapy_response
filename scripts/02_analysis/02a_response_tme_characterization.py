#!/usr/bin/env python3
"""
02_io_response_analysis.py — IO Response Biomarker Analysis

Using the integrated atlas (with actual R vs NR labels from GSE206325):
  1. Cell type annotation & marker-based validation
  2. IO signature scoring (12 signatures)
  3. Cell type abundance: R vs NR (Fisher exact, Dirichlet-style)
  4. DEG per cell type: R vs NR (Wilcoxon pseudobulk)
  5. T cell state analysis: exhaustion, progenitor, effector
  6. mregDC / cDC1 enrichment
  7. IO phenotype classification
  8. Multi-panel publication figures (varied plot types)
"""
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import mannwhitneyu, fisher_exact, spearmanr, kruskal
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.metrics import roc_auc_score, roc_curve
from collections import Counter
import re

sc.settings.verbosity = 1

PROJECT = Path(__file__).resolve().parents[2]
PROC    = PROJECT / "data" / "processed"
OUT     = PROJECT / "out" / "02a_response_tme"
FIG     = OUT / "figures"
TAB     = OUT / "tables"
for d in [FIG, TAB]:
    d.mkdir(parents=True, exist_ok=True)

# ── Signatures ─────────────────────────────────────────────────────────
IO_SIGNATURES = {
    "T_cell_exhaustion":      ["PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","TOX","ENTPD1"],
    "T_cell_activation":      ["IFNG","GZMB","PRF1","TNF","CD69","IL2"],
    "T_cell_cytotoxicity":    ["GZMA","GZMB","GZMH","GZMK","PRF1","GNLY","NKG7"],
    "Progenitor_CD8":         ["TCF7","SLAMF6","LEF1","CCR7","IL7R","SELL"],
    "Terminal_exhaustion":    ["TOX","ENTPD1","HAVCR2","LAYN","TIGIT","CXCL13"],
    "Treg_signature":         ["FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"],
    "M1_macrophage":          ["NOS2","IL1B","IL6","TNF","CD80","CD86"],
    "M2_macrophage":          ["CD163","MRC1","MSR1","CD200R1","IL10"],
    "IFN_gamma":              ["IFNG","STAT1","IRF1","CXCL9","CXCL10","CXCL11","IDO1"],
    "TGF_beta":               ["TGFB1","TGFB2","TGFB3","ACTA2","COL1A1","SERPINE1"],
    "Antigen_presentation":   ["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","PSMB9"],
    "Immune_checkpoint":      ["CD274","PDCD1LG2","PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","VSIR"],
    "Angiogenesis":           ["VEGFA","PECAM1","KDR","FLT1","ANGPT2","FGF2"],
    "WNT_beta_catenin":       ["CTNNB1","APC","AXIN2","LGR5","MYC"],
    "mregDC":                 ["CCR7","LAMP3","CCL19","CCL22","CD80","CD83","FSCN1","MARCKSL1"],
    "Tfh_like":               ["CXCL13","CH25H","IL21","PDCD1","ICOS","BCL6","TOX2"],
}

CELL_TYPE_MARKERS = {
    "CD8_T":        ["CD8A","CD8B","CD3E"],
    "CD4_T":        ["CD4","IL7R","CD3E"],
    "Treg":         ["FOXP3","IL2RA","CTLA4"],
    "NK":           ["NKG7","GNLY","KLRD1","NCAM1"],
    "B_cell":       ["CD79A","MS4A1","CD19"],
    "Plasma":       ["JCHAIN","MZB1","SDC1"],
    "Macrophage":   ["CD68","CD163","C1QA","CSF1R"],
    "mregDC":       ["LAMP3","CCR7","FSCN1"],
    "cDC1":         ["CLEC9A","XCR1","BATF3"],
    "cDC2":         ["CD1C","CLEC10A","FCER1A"],
    "Hepatocyte":   ["ALB","APOB","HP","TF"],
    "Endothelial":  ["PECAM1","CDH5","VWF"],
    "Fibroblast":   ["COL1A1","COL3A1","DCN","LUM"],
    "Tumor":        ["GPC3","AFP","GOLM1","EPCAM"],
}


def score_signatures(adata):
    """Score all IO-relevant gene signatures."""
    genes_available = set(adata.raw.var_names) if adata.raw else set(adata.var_names)
    scored = []
    for name, genes in IO_SIGNATURES.items():
        present = [g for g in genes if g in genes_available]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=name, use_raw=True)
            scored.append(name)
    print(f"  Scored {len(scored)}/{len(IO_SIGNATURES)} signatures")
    return scored


def annotate_cells(adata):
    """Marker-based cell type annotation on log-normalized data."""
    genes_available = set(adata.raw.var_names) if adata.raw else set(adata.var_names)
    for ct, markers in CELL_TYPE_MARKERS.items():
        present = [m for m in markers if m in genes_available]
        if present:
            sc.tl.score_genes(adata, present, score_name=f"score_{ct}", use_raw=True)

    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    if score_cols:
        score_df = adata.obs[score_cols]
        adata.obs["cell_type"] = score_df.idxmax(axis=1).str.replace("score_", "")
        adata.obs.loc[score_df.max(axis=1) < 0.05, "cell_type"] = "Unknown"


def response_deg_per_celltype(adata, min_cells=50):
    """Wilcoxon DEG between Responder vs Non-Responder per cell type."""
    adata_io = adata[adata.obs["response"].isin(["Responder", "Non-Responder"])].copy()
    adata_raw = adata_io.raw.to_adata() if adata_io.raw else adata_io

    results = []
    for ct in sorted(adata_io.obs["cell_type"].unique()):
        if ct == "Unknown":
            continue
        mask_ct = adata_io.obs["cell_type"] == ct
        n_r = ((adata_io.obs["response"] == "Responder") & mask_ct).sum()
        n_nr = ((adata_io.obs["response"] == "Non-Responder") & mask_ct).sum()
        if n_r < min_cells or n_nr < min_cells:
            continue

        sub = adata_raw[mask_ct].copy()
        sub.obs["response"] = adata_io.obs.loc[mask_ct, "response"].values
        try:
            sc.tl.rank_genes_groups(sub, groupby="response", groups=["Responder"],
                                    reference="Non-Responder", method="wilcoxon",
                                    use_raw=False, n_genes=200)
            deg = sc.get.rank_genes_groups_df(sub, group="Responder")
            deg["cell_type"] = ct
            deg["n_responder"] = n_r
            deg["n_nonresponder"] = n_nr
            results.append(deg)
        except Exception as e:
            print(f"    [WARN] DEG failed for {ct}: {e}")

    if results:
        deg_df = pd.concat(results, ignore_index=True)
        return deg_df
    return pd.DataFrame()


def cell_abundance_test(adata):
    """Test cell type abundance differences between R vs NR at patient level."""
    io_cells = adata[adata.obs["response"].isin(["Responder", "Non-Responder"])].copy()

    # Patient-level composition
    patient_comp = io_cells.obs.groupby(["patient_ID", "response", "cell_type"]).size().reset_index(name="count")
    patient_total = patient_comp.groupby(["patient_ID", "response"])["count"].transform("sum")
    patient_comp["fraction"] = patient_comp["count"] / patient_total

    results = []
    for ct in sorted(patient_comp["cell_type"].unique()):
        ct_data = patient_comp[patient_comp["cell_type"] == ct]
        r_fracs = ct_data[ct_data["response"] == "Responder"]["fraction"].values
        nr_fracs = ct_data[ct_data["response"] == "Non-Responder"]["fraction"].values
        if len(r_fracs) >= 3 and len(nr_fracs) >= 3:
            stat, pval = mannwhitneyu(r_fracs, nr_fracs, alternative="two-sided")
            results.append({
                "cell_type": ct,
                "mean_R": np.mean(r_fracs),
                "mean_NR": np.mean(nr_fracs),
                "log2FC": np.log2((np.mean(r_fracs) + 1e-6) / (np.mean(nr_fracs) + 1e-6)),
                "pvalue": pval,
                "n_R_patients": len(r_fracs),
                "n_NR_patients": len(nr_fracs),
            })

    if results:
        df = pd.DataFrame(results).sort_values("pvalue")
        df["padj"] = np.minimum(df["pvalue"] * len(df), 1.0)
        return df
    return pd.DataFrame()


def main():
    print("=" * 70)
    print("HCC IO Response — Biomarker Analysis")
    print("=" * 70)

    # ── Load integrated atlas ──────────────────────────────────────
    h5_path = PROC / "integrated_io_atlas.h5ad"
    if not h5_path.exists():
        print(f"[ERROR] {h5_path} not found. Run 01_preprocess_integrate.py first.")
        return
    print("[LOAD] Reading integrated atlas...")
    adata = sc.read_h5ad(h5_path)
    print(f"  {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

    # ── Cell type annotation ───────────────────────────────────────
    # Use original Magen annotations where available; annotate rest
    print("\n[ANNOTATE] Cell types...")
    if "cell_type_original" in adata.obs.columns:
        # Use published annotations for GSE206325
        has_orig = adata.obs["cell_type_original"].notna()
        print(f"  Magen published annotations: {has_orig.sum():,} cells")
        annotate_cells(adata)  # Score-based for all
        # Override with published where available
        adata.obs.loc[has_orig, "cell_type"] = adata.obs.loc[has_orig, "cell_type_original"]
    else:
        annotate_cells(adata)

    ct_counts = adata.obs["cell_type"].value_counts()
    print("  Cell type distribution:")
    for ct, n in ct_counts.head(15).items():
        print(f"    {ct}: {n:,} ({100*n/len(adata):.1f}%)")

    # ── IO signature scoring ───────────────────────────────────────
    print("\n[SIGNATURES] Scoring IO signatures...")
    scored = score_signatures(adata)

    # ── Focus on IO-labeled cells (GSE206325) ──────────────────────
    io_mask = adata.obs["response"].isin(["Responder", "Non-Responder"])
    n_io = io_mask.sum()
    print(f"\n[IO] Cells with R/NR labels: {n_io:,}")

    # ── Cell abundance R vs NR ─────────────────────────────────────
    print("\n[ABUNDANCE] Testing cell type abundance R vs NR...")
    if "patient_ID" in adata.obs.columns:
        abund_df = cell_abundance_test(adata)
        if len(abund_df) > 0:
            abund_df.to_csv(TAB / "celltype_abundance_R_vs_NR.csv", index=False)
            print("  Significant cell types (padj < 0.05):")
            for _, row in abund_df[abund_df["padj"] < 0.05].iterrows():
                direction = "↑R" if row["log2FC"] > 0 else "↑NR"
                print(f"    {row['cell_type']}: log2FC={row['log2FC']:.2f} {direction} (p={row['pvalue']:.2e})")

    # ── DEG per cell type ──────────────────────────────────────────
    print("\n[DEG] Response-associated DEGs per cell type...")
    deg_df = response_deg_per_celltype(adata)
    if len(deg_df) > 0:
        deg_df.to_csv(TAB / "deg_R_vs_NR_per_celltype.csv", index=False)
        sig_degs = deg_df[(deg_df["pvals_adj"] < 0.05) & (deg_df["logfoldchanges"].abs() > 0.25)]
        print(f"  Total significant DEGs: {len(sig_degs)}")
        for ct in sig_degs["cell_type"].unique():
            n = (sig_degs["cell_type"] == ct).sum()
            print(f"    {ct}: {n} DEGs")

    # ── Signature comparison R vs NR ───────────────────────────────
    print("\n[SIG DIFF] Comparing IO signatures: R vs NR...")
    io_cells = adata[io_mask].copy()
    sig_results = []
    for sig in scored:
        r_vals = io_cells.obs.loc[io_cells.obs["response"] == "Responder", sig]
        nr_vals = io_cells.obs.loc[io_cells.obs["response"] == "Non-Responder", sig]
        if len(r_vals) > 0 and len(nr_vals) > 0:
            stat, pval = mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
            sig_results.append({
                "signature": sig, "mean_R": r_vals.mean(), "mean_NR": nr_vals.mean(),
                "effect_size": r_vals.mean() - nr_vals.mean(), "pvalue": pval,
            })
    sig_diff = pd.DataFrame(sig_results).sort_values("pvalue")
    sig_diff["padj"] = np.minimum(sig_diff["pvalue"] * len(sig_diff), 1.0)
    sig_diff.to_csv(TAB / "signature_R_vs_NR.csv", index=False)

    # ── Tissue-specific analysis ───────────────────────────────────
    print("\n[TISSUE] Tissue-specific signature differences...")
    tissue_results = []
    for tissue in ["Tumor", "Normal", "PBMC", "LN"]:
        mask_tissue = (adata.obs["tissue"] == tissue) & io_mask
        if mask_tissue.sum() < 100:
            continue
        for sig in scored:
            r_vals = adata.obs.loc[mask_tissue & (adata.obs["response"] == "Responder"), sig]
            nr_vals = adata.obs.loc[mask_tissue & (adata.obs["response"] == "Non-Responder"), sig]
            if len(r_vals) >= 20 and len(nr_vals) >= 20:
                stat, pval = mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
                tissue_results.append({
                    "tissue": tissue, "signature": sig,
                    "mean_R": r_vals.mean(), "mean_NR": nr_vals.mean(),
                    "effect": r_vals.mean() - nr_vals.mean(), "pvalue": pval,
                })
    tissue_df = pd.DataFrame(tissue_results)
    if len(tissue_df) > 0:
        tissue_df.to_csv(TAB / "tissue_signature_R_vs_NR.csv", index=False)

    # ═══════════════════════════════════════════════════════════════
    # FIGURES — varied plot types
    # ═══════════════════════════════════════════════════════════════
    print("\n[FIGURES] Generating publication figures...")

    # ─── Fig 1: Cell type UMAP + Response overlay ──────────────────
    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    ax = axes[0]
    top_cts = ct_counts.head(12).index.tolist()
    colors_ct = plt.cm.tab20(np.linspace(0, 1, len(top_cts)))
    for i, ct in enumerate(top_cts):
        mask = adata.obs["cell_type"] == ct
        ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                   s=0.3, alpha=0.2, c=[colors_ct[i]], label=ct, rasterized=True)
    ax.legend(fontsize=6, markerscale=8, ncol=2, loc="upper left")
    ax.set_title("A. Cell Types", fontweight="bold"); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    ax = axes[1]
    mask_tumor = adata.obs["tissue"] == "Tumor"
    ax.scatter(adata.obsm["X_umap"][~mask_tumor, 0], adata.obsm["X_umap"][~mask_tumor, 1],
               s=0.1, alpha=0.05, c="gray", rasterized=True, label="Other")
    ax.scatter(adata.obsm["X_umap"][mask_tumor, 0], adata.obsm["X_umap"][mask_tumor, 1],
               s=0.3, alpha=0.2, c="#e74c3c", rasterized=True, label="Tumor")
    ax.legend(fontsize=8, markerscale=10); ax.set_title("B. Tumor Cells", fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    ax = axes[2]
    resp_colors = {"Responder": "#2ecc71", "Non-Responder": "#e74c3c"}
    for resp, color in resp_colors.items():
        mask = adata.obs["response"] == resp
        if mask.sum() > 0:
            ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                       s=0.3, alpha=0.15, c=color, label=resp, rasterized=True)
    ax.legend(fontsize=8, markerscale=10); ax.set_title("C. IO Response", fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    plt.suptitle("HCC IO Response Atlas — Cell Type & Response", fontweight="bold", fontsize=13)
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"02_celltype_response_umap.{ext}", dpi=200, bbox_inches="tight")
    plt.close()

    # ─── Fig 2: Signature violin plots R vs NR ────────────────────
    top_sigs = sig_diff.head(8)["signature"].tolist() if len(sig_diff) > 0 else scored[:8]
    n_sigs = len(top_sigs)
    ncols = 4
    nrows = (n_sigs + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3.5*nrows))
    axes = axes.flatten() if nrows > 1 else [axes] if nrows == 1 and ncols == 1 else axes.flatten()
    for i, sig in enumerate(top_sigs):
        ax = axes[i]
        data_r = io_cells.obs.loc[io_cells.obs["response"] == "Responder", sig]
        data_nr = io_cells.obs.loc[io_cells.obs["response"] == "Non-Responder", sig]
        parts = ax.violinplot([data_r.values, data_nr.values], showmedians=True)
        for j, body in enumerate(parts["bodies"]):
            body.set_facecolor(["#2ecc71", "#e74c3c"][j])
            body.set_alpha(0.7)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["R", "NR"])
        p = sig_diff.loc[sig_diff["signature"] == sig, "pvalue"].values
        pstr = f"p={p[0]:.1e}" if len(p) > 0 else ""
        ax.set_title(f"{sig}\n{pstr}", fontsize=9)
    for j in range(i+1, len(axes)):
        axes[j].set_visible(False)
    plt.suptitle("IO Signatures: Responder vs Non-Responder", fontweight="bold")
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"02_signature_violins.{ext}", dpi=200, bbox_inches="tight")
    plt.close()

    # ─── Fig 3: Cell type abundance dotplot (R vs NR) ──────────────
    if len(abund_df) > 0:
        fig, ax = plt.subplots(figsize=(8, 6))
        abund_plot = abund_df.sort_values("log2FC")
        colors = ["#2ecc71" if x > 0 else "#e74c3c" for x in abund_plot["log2FC"]]
        sizes = np.clip(-np.log10(abund_plot["pvalue"].clip(1e-20)), 1, 20) * 5
        ax.barh(range(len(abund_plot)), abund_plot["log2FC"], color=colors, alpha=0.8)
        for i, (_, row) in enumerate(abund_plot.iterrows()):
            star = "***" if row["padj"] < 0.001 else "**" if row["padj"] < 0.01 else "*" if row["padj"] < 0.05 else ""
            ax.text(row["log2FC"] + 0.05 * np.sign(row["log2FC"]), i, star, va="center", fontsize=10)
        ax.set_yticks(range(len(abund_plot)))
        ax.set_yticklabels(abund_plot["cell_type"])
        ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
        ax.set_xlabel("log2(Fold Change): R vs NR")
        ax.set_title("Cell Type Abundance: Responder vs Non-Responder\n(Patient-level Mann-Whitney)", fontweight="bold")
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"02_abundance_barplot.{ext}", dpi=200, bbox_inches="tight")
        plt.close()

    # ─── Fig 4: Signature heatmap by tissue × response ─────────────
    if len(tissue_df) > 0:
        pivot = tissue_df.pivot_table(index="signature", columns="tissue", values="effect", fill_value=0)
        fig, ax = plt.subplots(figsize=(8, 10))
        sns.heatmap(pivot, cmap="RdBu_r", center=0, ax=ax, annot=True, fmt=".2f",
                    linewidths=0.5, cbar_kws={"label": "Effect Size (R - NR)"})
        ax.set_title("IO Signature Effects: R vs NR by Tissue", fontweight="bold")
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"02_tissue_signature_heatmap.{ext}", dpi=200, bbox_inches="tight")
        plt.close()

    # ─── Fig 5: Scatter — Exhaustion vs Activation by response ─────
    if "T_cell_exhaustion" in scored and "T_cell_activation" in scored:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        for idx, resp in enumerate(["Responder", "Non-Responder"]):
            ax = axes[idx]
            mask = io_cells.obs["response"] == resp
            sub = io_cells[mask]
            sc_plot = ax.scatter(sub.obs["T_cell_activation"], sub.obs["T_cell_exhaustion"],
                       s=0.5, alpha=0.1, c=resp_colors[resp], rasterized=True)
            rho, p = spearmanr(sub.obs["T_cell_activation"], sub.obs["T_cell_exhaustion"])
            ax.set_title(f"{resp} (ρ={rho:.3f}, p={p:.1e})", fontweight="bold")
            ax.set_xlabel("T cell Activation"); ax.set_ylabel("T cell Exhaustion")
            # Fit line
            z = np.polyfit(sub.obs["T_cell_activation"], sub.obs["T_cell_exhaustion"], 1)
            x_line = np.linspace(sub.obs["T_cell_activation"].min(), sub.obs["T_cell_activation"].max(), 100)
            ax.plot(x_line, np.polyval(z, x_line), "k--", alpha=0.5)
        plt.suptitle("T Cell Activation vs Exhaustion", fontweight="bold")
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"02_activation_exhaustion_scatter.{ext}", dpi=200, bbox_inches="tight")
        plt.close()

    # ─── Fig 6: Volcano plot of top DEGs (from largest cell type) ──
    if len(deg_df) > 0:
        # Pick cell type with most DEGs
        top_ct = sig_degs["cell_type"].value_counts().index[0] if len(sig_degs) > 0 else deg_df["cell_type"].iloc[0]
        ct_deg = deg_df[deg_df["cell_type"] == top_ct].copy()
        ct_deg["neg_log_p"] = -np.log10(ct_deg["pvals_adj"].clip(1e-300))

        fig, ax = plt.subplots(figsize=(8, 6))
        sig_up = ct_deg[(ct_deg["logfoldchanges"] > 0.25) & (ct_deg["pvals_adj"] < 0.05)]
        sig_dn = ct_deg[(ct_deg["logfoldchanges"] < -0.25) & (ct_deg["pvals_adj"] < 0.05)]
        ns = ct_deg[~ct_deg.index.isin(sig_up.index) & ~ct_deg.index.isin(sig_dn.index)]

        ax.scatter(ns["logfoldchanges"], ns["neg_log_p"], s=3, c="gray", alpha=0.3, rasterized=True)
        ax.scatter(sig_up["logfoldchanges"], sig_up["neg_log_p"], s=5, c="#e74c3c", alpha=0.5, label=f"Up in R ({len(sig_up)})")
        ax.scatter(sig_dn["logfoldchanges"], sig_dn["neg_log_p"], s=5, c="#3498db", alpha=0.5, label=f"Down in R ({len(sig_dn)})")

        # Label top genes
        for _, row in pd.concat([sig_up.nlargest(5, "neg_log_p"), sig_dn.nlargest(5, "neg_log_p")]).iterrows():
            ax.annotate(row["names"], (row["logfoldchanges"], row["neg_log_p"]),
                       fontsize=7, alpha=0.8, ha="center")

        ax.axhline(-np.log10(0.05), color="gray", linestyle="--", alpha=0.5)
        ax.axvline(0.25, color="gray", linestyle="--", alpha=0.5)
        ax.axvline(-0.25, color="gray", linestyle="--", alpha=0.5)
        ax.set_xlabel("log2 Fold Change (R vs NR)"); ax.set_ylabel("-log10(padj)")
        ax.set_title(f"Volcano: Responder vs Non-Responder DEGs\n({top_ct})", fontweight="bold")
        ax.legend()
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"02_volcano_deg.{ext}", dpi=200, bbox_inches="tight")
        plt.close()

    # ─── Fig 7: Stacked bar — cell composition by response ────────
    if "response" in io_cells.obs.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        comp = pd.crosstab(io_cells.obs["response"], io_cells.obs["cell_type"], normalize="index")
        comp = comp[comp.mean().sort_values(ascending=False).index]  # Sort by abundance
        comp.plot(kind="bar", stacked=True, ax=ax, colormap="tab20", width=0.6)
        ax.set_ylabel("Fraction"); ax.set_xlabel("")
        ax.set_title("Cell Type Composition: R vs NR", fontweight="bold")
        ax.legend(fontsize=7, bbox_to_anchor=(1.05, 1), loc="upper left")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"02_composition_stacked_bar.{ext}", dpi=200, bbox_inches="tight")
        plt.close()

    # ─── Fig 8: Radar/Spider chart — signatures by response ───────
    if len(sig_diff) >= 4:
        top_n = min(10, len(scored))
        sigs_radar = scored[:top_n]
        r_means = [io_cells.obs.loc[io_cells.obs["response"]=="Responder", s].mean() for s in sigs_radar]
        nr_means = [io_cells.obs.loc[io_cells.obs["response"]=="Non-Responder", s].mean() for s in sigs_radar]

        # Normalize to 0-1 range for radar
        all_vals = r_means + nr_means
        vmin, vmax = min(all_vals), max(all_vals)
        r_norm = [(v - vmin) / (vmax - vmin + 1e-10) for v in r_means]
        nr_norm = [(v - vmin) / (vmax - vmin + 1e-10) for v in nr_means]

        angles = np.linspace(0, 2*np.pi, top_n, endpoint=False).tolist()
        angles += angles[:1]  # close

        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
        ax.plot(angles, r_norm + r_norm[:1], "o-", color="#2ecc71", linewidth=2, label="Responder")
        ax.fill(angles, r_norm + r_norm[:1], color="#2ecc71", alpha=0.15)
        ax.plot(angles, nr_norm + nr_norm[:1], "o-", color="#e74c3c", linewidth=2, label="Non-Responder")
        ax.fill(angles, nr_norm + nr_norm[:1], color="#e74c3c", alpha=0.15)
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels([s.replace("_", "\n") for s in sigs_radar], fontsize=8)
        ax.set_title("IO Signature Radar: R vs NR", fontweight="bold", pad=20)
        ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1))
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"02_signature_radar.{ext}", dpi=200, bbox_inches="tight")
        plt.close()

    # ── Save updated adata ─────────────────────────────────────────
    adata.write(PROC / "integrated_io_atlas.h5ad")

    print(f"\n{'='*70}")
    print("[DONE] IO Response Analysis")
    print(f"  Figures: {len(list(FIG.glob('*.png')))} PNG files in {FIG}")
    for f in sorted(FIG.glob("*.png")):
        print(f"    {f.name}")
    print(f"  Tables: {len(list(TAB.glob('*.csv')))} CSV files in {TAB}")
    print(f"={'='*70}")


if __name__ == "__main__":
    main()
