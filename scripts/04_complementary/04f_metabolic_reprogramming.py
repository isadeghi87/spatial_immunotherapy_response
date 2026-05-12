#!/usr/bin/env python3
"""04f – Metabolic reprogramming profiles from real atlas data.

Scores glycolysis and OXPHOS gene signatures per cell, then plots
violin distributions stratified by cell type and response.
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
import sys; sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
from config.publication_style import apply_style, PAL, CMAP, label_clean
apply_style()
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
ATLAS   = PROJECT / "data" / "processed" / "integrated_io_atlas.h5ad"
OUT     = PROJECT / "out" / "04_complementary" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# MSigDB Hallmark gene sets (core members present in most scRNA-seq)
GLYCOLYSIS_GENES = [
    "HK1", "HK2", "GPI", "PFKL", "PFKM", "ALDOA", "ALDOB", "TPI1",
    "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "LDHB",
    "SLC2A1", "SLC2A3", "PFKFB3",
]
OXPHOS_GENES = [
    "NDUFA1", "NDUFA2", "NDUFB1", "NDUFB2", "NDUFS1", "NDUFS2",
    "SDHA", "SDHB", "SDHC", "UQCRB", "UQCRH", "UQCRC1", "UQCRC2",
    "COX5A", "COX5B", "COX6A1", "COX7A2", "COX7B", "ATP5F1A", "ATP5F1B",
    "ATP5MC1", "ATP5PO",
]


def _resolve_celltype_col(adata):
    for col in ("cell_type", "celltype_sub", "Cell.type", "celltype_global", "celltype"):
        if col in adata.obs.columns:
            return col
    raise KeyError("No cell-type column found")


def plot_metabolism(adata_path=None):
    if adata_path is None:
        adata_path = str(ATLAS)

    print("[04f] Loading atlas …")
    adata = sc.read_h5ad(adata_path)
    ct_col = _resolve_celltype_col(adata)
    print(f"  Shape: {adata.shape}, cell-type col: {ct_col}")

    # Score gene sets
    glyc_present = [g for g in GLYCOLYSIS_GENES if g in adata.var_names]
    oxphos_present = [g for g in OXPHOS_GENES if g in adata.var_names]
    print(f"  Glycolysis genes found: {len(glyc_present)}/{len(GLYCOLYSIS_GENES)}")
    print(f"  OXPHOS genes found: {len(oxphos_present)}/{len(OXPHOS_GENES)}")

    sc.tl.score_genes(adata, gene_list=glyc_present, score_name="glycolysis_score")
    sc.tl.score_genes(adata, gene_list=oxphos_present, score_name="oxphos_score")

    # Focus on key immune cell types with response info
    has_response = "response" in adata.obs.columns
    top_cts = adata.obs[ct_col].value_counts().nlargest(6).index.tolist()
    mask = adata.obs[ct_col].isin(top_cts)
    if has_response:
        mask = mask & adata.obs["response"].notna()

    df = adata.obs.loc[mask, [ct_col, "glycolysis_score", "oxphos_score"]].copy()
    if has_response:
        df["response"] = adata.obs.loc[mask, "response"]

    # Subsample for plotting speed (max 50k)
    if len(df) > 50000:
        df = df.sample(50000, random_state=42)

    if has_response:
        # ── Separate violin per cell type, coloured by response ──────
        from scipy.stats import mannwhitneyu

        fig, axes = plt.subplots(2, len(top_cts), figsize=(3.2 * len(top_cts), 9),
                                 squeeze=False, sharey="row")

        palette = {"Responder": "#0077BB", "Non-Responder": "#CC3311"}

        for row_idx, (score_col, score_label) in enumerate([
            ("glycolysis_score", "Glycolysis"),
            ("oxphos_score", "OXPHOS"),
        ]):
            for col_idx, ct_name in enumerate(top_cts):
                ax = axes[row_idx][col_idx]
                ct_df = df[df[ct_col] == ct_name]

                sns.violinplot(data=ct_df, x="response", y=score_col,
                               palette=palette, ax=ax, inner="box", cut=0,
                               linewidth=0.8, order=["Responder", "Non-Responder"],
                               width=0.7)

                # Mann-Whitney p-value
                r_vals = ct_df.loc[ct_df["response"] == "Responder", score_col].dropna()
                nr_vals = ct_df.loc[ct_df["response"] == "Non-Responder", score_col].dropna()
                if len(r_vals) >= 10 and len(nr_vals) >= 10:
                    _, p = mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
                    stars = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                    ax.set_title(f"{ct_name}\n{stars} (p={p:.2g})", fontsize=8, fontweight="bold")
                else:
                    ax.set_title(ct_name, fontsize=8, fontweight="bold")

                ax.set_xlabel("")
                ax.set_ylabel(f"{score_label} Score" if col_idx == 0 else "")
                ax.tick_params(axis="x", rotation=30, labelsize=7)
                ax.tick_params(axis="y", labelsize=8)

        plt.suptitle("Metabolic Reprogramming: Glycolysis & OXPHOS by Cell Type × Response",
                      fontweight="bold", fontsize=12, y=1.01)
    else:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        df["Category"] = df[ct_col].astype(str)

        sns.violinplot(data=df, x="Category", y="glycolysis_score", ax=axes[0],
                       palette="Set2", inner="box", cut=0, linewidth=0.8)
        axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=45, ha="right", fontsize=9)
        axes[0].set_title("A. Glycolysis Signature Score", fontweight="bold")
        axes[0].set_ylabel("Score")
        axes[0].set_xlabel("")

        sns.violinplot(data=df, x="Category", y="oxphos_score", ax=axes[1],
                       palette="Set2", inner="box", cut=0, linewidth=0.8)
        axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=45, ha="right", fontsize=9)
        axes[1].set_title("B. OXPHOS Signature Score", fontweight="bold")
        axes[1].set_ylabel("Score")
        axes[1].set_xlabel("")

    plt.tight_layout()
    plt.savefig(OUT / "04f_metabolic_reprogramming.png", dpi=300, bbox_inches="tight")
    plt.savefig(OUT / "04f_metabolic_reprogramming.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {OUT / '04f_metabolic_reprogramming.png'}")


if __name__ == "__main__":
    plot_metabolism()
