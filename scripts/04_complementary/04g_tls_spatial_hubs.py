#!/usr/bin/env python3
"""04g – Tertiary lymphoid structure (TLS) signature analysis from real data.

Scores a TLS 12-gene signature per cell and visualises spatial enrichment
on the UMAP manifold, plus TLS-score distribution by response.
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
from scipy.stats import mannwhitneyu

PROJECT = Path(__file__).resolve().parents[2]
ATLAS   = PROJECT / "data" / "processed" / "integrated_io_atlas.h5ad"
OUT     = PROJECT / "out" / "04_complementary" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# 12-gene TLS signature (Cabrita et al. 2020 / Helmink et al. 2020)
TLS_GENES = [
    "CCL19", "CCL21", "CXCL13", "CCR7", "SELL", "LAMP3",
    "CD79A", "MS4A1", "CR2", "CD3D", "CD3E", "CXCR5",
]


def plot_tls(adata_path=None):
    if adata_path is None:
        adata_path = str(ATLAS)

    print("[04g] Loading atlas …")
    adata = sc.read_h5ad(adata_path)
    print(f"  Shape: {adata.shape}")

    tls_present = [g for g in TLS_GENES if g in adata.var_names]
    print(f"  TLS genes found: {len(tls_present)}/{len(TLS_GENES)}")

    sc.tl.score_genes(adata, gene_list=tls_present, score_name="tls_score")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # ── Panel A: UMAP coloured by TLS score ──────────────────────────
    ax = axes[0]
    umap = adata.obsm["X_umap"]
    tls = adata.obs["tls_score"].values

    # Subsample for plotting speed
    n = len(umap)
    step = max(1, n // 50000)
    idx = np.arange(0, n, step)

    # Clip extreme outliers for colormap contrast
    tls_plot = tls[idx]
    vmin, vmax = np.nanpercentile(tls_plot, [2, 98])

    # Sort by score for proper layering (low in back, high in front)
    order = np.argsort(tls_plot)
    sc_plot = ax.scatter(umap[idx[order], 0], umap[idx[order], 1],
                         c=tls_plot[order], cmap="viridis",
                         vmin=vmin, vmax=vmax,
                         s=3, alpha=0.7, edgecolor="none", rasterized=True)
    ax.set_title("A. TLS 12-Gene Signature on UMAP", fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(sc_plot, ax=ax, label="TLS Enrichment Score", shrink=0.7)

    # ── Panel B: TLS score distribution by response ──────────────────
    ax = axes[1]
    has_response = "response" in adata.obs.columns
    if has_response:
        df = adata.obs[["response", "tls_score"]].dropna()
        if len(df) > 50000:
            df = df.sample(50000, random_state=42)
        groups = sorted(df["response"].unique())
        palette = {"Responder": "#9b59b6", "Non-Responder": "#34495e"}
        sns.violinplot(data=df, x="response", y="tls_score", palette=palette,
                       ax=ax, inner="box", cut=0, linewidth=0.8)

        # Statistical test
        if len(groups) == 2:
            g1 = df.loc[df["response"] == groups[0], "tls_score"]
            g2 = df.loc[df["response"] == groups[1], "tls_score"]
            stat, pval = mannwhitneyu(g1, g2, alternative="two-sided")
            # Format p-value properly (avoid 0.00e+00)
            if pval < 1e-300:
                p_str = "p < 1e-300"
            elif pval < 0.001:
                p_str = f"p = {pval:.2e}"
            else:
                p_str = f"p = {pval:.4f}"
            ax.text(0.5, 0.95, f"Mann-Whitney {p_str}", transform=ax.transAxes,
                    ha="center", fontsize=11, fontweight="bold",
                    bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.3"))
    else:
        ax.text(0.5, 0.5, "No response column", ha="center", va="center",
                transform=ax.transAxes, fontsize=14)

    ax.set_title("B. TLS Score by IO Response", fontweight="bold")
    ax.set_ylabel("TLS Enrichment Score")
    ax.set_xlabel("")

    plt.tight_layout()
    plt.savefig(OUT / "04g_tls_spatial_hubs.png", dpi=300, bbox_inches="tight")
    plt.savefig(OUT / "04g_tls_spatial_hubs.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {OUT / '04g_tls_spatial_hubs.png'}")


if __name__ == "__main__":
    plot_tls()
