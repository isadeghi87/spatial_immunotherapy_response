#!/usr/bin/env python3
"""04c – T-cell pseudotime trajectory analysis using real atlas data.

Subsets CD8 T cells, computes diffusion pseudotime (DPT), and plots:
  A. UMAP coloured by pseudotime with stream arrows
  B. Smoothed gene dynamics along pseudotime for key markers
"""
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
from config.publication_style import apply_style, PAL, CMAP, label_clean
apply_style()
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter1d
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
ATLAS   = PROJECT / "data" / "processed" / "integrated_io_atlas.h5ad"
OUT     = PROJECT / "out" / "04_complementary" / "figures"
OUT.mkdir(parents=True, exist_ok=True)


def _resolve_celltype_col(adata):
    for col in ("cell_type", "celltype_sub", "Cell.type", "celltype_global", "celltype"):
        if col in adata.obs.columns:
            return col
    raise KeyError("No cell-type column found")


def set_style():
    plt.style.use("seaborn-v0_8-whitegrid")
    matplotlib.rcParams.update({
        "axes.facecolor": "#fafafa", "figure.facecolor": "white",
        "axes.edgecolor": "#dddddd", "grid.color": "#ebebeb",
        "font.family": "sans-serif", "axes.titlesize": 16,
        "axes.labelsize": 13, "xtick.labelsize": 11,
        "ytick.labelsize": 11, "legend.fontsize": 11,
    })


def plot_tcell_trajectory(adata_path=None):
    if adata_path is None:
        adata_path = str(ATLAS)

    print("[04c] Loading atlas …")
    adata = sc.read_h5ad(adata_path)
    print(f"  Shape: {adata.shape}")

    ct_col = _resolve_celltype_col(adata)

    # ── subset to CD8 T cells ────────────────────────────────────────
    cd8_mask = adata.obs[ct_col].astype(str).str.contains("CD8|T cell", case=False, na=False)
    adata_cd8 = adata[cd8_mask].copy()
    print(f"  CD8/T subset: {adata_cd8.n_obs} cells")

    if adata_cd8.n_obs < 100:
        raise ValueError(f"Too few CD8 T cells ({adata_cd8.n_obs}) for trajectory analysis")

    # Subsample for computational efficiency
    if adata_cd8.n_obs > 50_000:
        sc.pp.subsample(adata_cd8, n_obs=50_000, random_state=42)
        print(f"  Subsampled to {adata_cd8.n_obs} cells")

    # ── recompute neighbors + diffmap + DPT ──────────────────────────
    print("  Computing neighbors on CD8 subset …")
    if "X_scVI" in adata_cd8.obsm:
        sc.pp.neighbors(adata_cd8, use_rep="X_scVI", n_neighbors=15)
    elif "X_pca_harmony" in adata_cd8.obsm:
        sc.pp.neighbors(adata_cd8, use_rep="X_pca_harmony", n_neighbors=15)
    else:
        sc.pp.neighbors(adata_cd8, n_neighbors=15)
    sc.tl.umap(adata_cd8)
    sc.tl.diffmap(adata_cd8)

    # Pick root cell: the one with highest expression of naive marker (CCR7 or TCF7)
    naive_genes = [g for g in ["CCR7", "TCF7", "SELL", "IL7R"] if g in adata_cd8.var_names]
    if naive_genes:
        expr = np.array(adata_cd8[:, naive_genes].X.mean(axis=1)).flatten()
        root_idx = int(np.argmax(expr))
    else:
        root_idx = 0

    adata_cd8.uns["iroot"] = root_idx
    sc.tl.dpt(adata_cd8)
    print("  Diffusion pseudotime computed.")

    # ── Plot ─────────────────────────────────────────────────────────
    set_style()
    fig, axes = plt.subplots(1, 2, figsize=(18, 7),
                              gridspec_kw={"width_ratios": [1, 1.3], "wspace": 0.28})

    cmap_pt = LinearSegmentedColormap.from_list(
        "pseudotime", ["#3b4cc0", "#8a41c6", "#d2408c", "#f37a59", "#fac05e"])

    # ── Panel A: UMAP with pseudotime + smooth trajectory arrows ─────
    ax = axes[0]
    umap = adata_cd8.obsm["X_umap"]
    dpt = adata_cd8.obs["dpt_pseudotime"].values

    # Clip infinite pseudotime values for cleaner plotting
    dpt_clipped = dpt.copy()
    finite_mask_plot = np.isfinite(dpt_clipped)
    if finite_mask_plot.sum() < len(dpt_clipped):
        max_finite = np.nanmax(dpt_clipped[finite_mask_plot])
        dpt_clipped[~finite_mask_plot] = max_finite

    # Sort by pseudotime for proper layering
    order = np.argsort(dpt_clipped)
    sc_plot = ax.scatter(umap[order, 0], umap[order, 1], c=dpt_clipped[order],
                         cmap=cmap_pt, s=6, alpha=0.7, edgecolor="none", rasterized=True)

    # Annotate cell states along the trajectory
    ct_col = _resolve_celltype_col(adata_cd8)
    ct_vals = adata_cd8.obs[ct_col].astype(str).values
    unique_cts = pd.Series(ct_vals).value_counts()
    # Label the top cell subtypes at their UMAP centroids
    for ct_name in unique_cts.index[:6]:  # top 6 cell states
        ct_mask = ct_vals == ct_name
        if ct_mask.sum() < 20:
            continue
        cx, cy = umap[ct_mask, 0].mean(), umap[ct_mask, 1].mean()
        short_name = ct_name.split("_")[-1] if "_" in ct_name else ct_name
        short_name = short_name.replace("-", "\n", 1) if len(short_name) > 15 else short_name
        ax.annotate(short_name, (cx, cy), fontsize=7, fontweight="bold",
                    ha="center", va="center", color="black",
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                              edgecolor="gray", alpha=0.75))

    # Smooth trajectory arrows via pseudotime-binned UMAP centroids
    finite = np.isfinite(dpt)
    dpt_fin, umap_fin = dpt[finite], umap[finite]
    n_arrow_bins = 25
    edges = np.linspace(dpt_fin.min(), dpt_fin.max(), n_arrow_bins + 1)
    centroids = []
    for i in range(n_arrow_bins):
        m = (dpt_fin >= edges[i]) & (dpt_fin < edges[i + 1])
        if m.sum() > 5:
            centroids.append(umap_fin[m].mean(axis=0))
    if len(centroids) > 3:
        centroids = np.array(centroids)
        centroids[:, 0] = gaussian_filter1d(centroids[:, 0], sigma=1.5)
        centroids[:, 1] = gaussian_filter1d(centroids[:, 1], sigma=1.5)
        for i in range(len(centroids) - 1):
            ax.annotate("", xy=centroids[i + 1], xytext=centroids[i],
                         arrowprops=dict(arrowstyle="-|>", color="black",
                                         lw=2.0, alpha=0.7,
                                         connectionstyle="arc3,rad=0.08"))

    ax.set_title("A. CD8 T Cell Trajectory Manifold", fontweight="bold", fontsize=16, pad=15)
    ax.set_xlabel("UMAP 1", fontweight="bold")
    ax.set_ylabel("UMAP 2", fontweight="bold")
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

    cb = plt.colorbar(sc_plot, ax=ax, shrink=0.8, pad=0.03)
    cb.set_label("Diffusion Pseudotime", fontweight="bold")
    cb.outline.set_visible(False)

    # ── Panel B: Gene dynamics (Gaussian-smoothed) ───────────────────
    ax = axes[1]
    marker_genes = {
        "IL7R":   ("#3498db", "IL7R (Naive/Memory)"),
        "GZMB":   ("#0077BB", "GZMB (Effector)"),
        "PRF1":   ("#27ae60", "PRF1 (Cytotoxic)"),
        "TOX":    ("#CC3311", "TOX (Exhaustion)"),
        "TCF7":   ("#9b59b6", "TCF7 (Progenitor)"),
        "HAVCR2": ("#e67e22", "HAVCR2/TIM3 (Terminal Exh)"),
        "PDCD1":  ("#c0392b", "PDCD1/PD-1 (Checkpoint)"),
    }

    n_bins = 60
    sigma_smooth = 4
    finite_mask = np.isfinite(dpt)
    dpt_f = dpt[finite_mask]
    bins = np.linspace(dpt_f.min(), dpt_f.max(), n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_idx = np.clip(np.digitize(dpt_f, bins) - 1, 0, n_bins - 1)

    for gene, (color, label) in marker_genes.items():
        if gene not in adata_cd8.var_names:
            continue
        gene_idx = list(adata_cd8.var_names).index(gene)
        X_sub = adata_cd8.X[finite_mask]
        expr = np.array(X_sub[:, gene_idx].toarray()).flatten() \
               if hasattr(X_sub, "toarray") else np.array(X_sub[:, gene_idx]).flatten()

        bin_means = np.zeros(n_bins)
        bin_sems  = np.zeros(n_bins)
        for b in range(n_bins):
            m = bin_idx == b
            if m.sum() > 0:
                v = expr[m]
                bin_means[b] = v.mean()
                bin_sems[b]  = v.std() / np.sqrt(max(m.sum(), 1))

        smoothed     = gaussian_filter1d(bin_means, sigma=sigma_smooth)
        smoothed_sem = gaussian_filter1d(bin_sems,  sigma=sigma_smooth)

        ax.plot(bin_centers, smoothed, label=label, color=color, lw=2.5, alpha=0.95)
        ax.fill_between(bin_centers, smoothed - smoothed_sem, smoothed + smoothed_sem,
                        color=color, alpha=0.15, lw=0)

    ax.set_title("B. Gene Dynamics Along Pseudotime", fontweight="bold", fontsize=16, pad=15)
    ax.set_xlabel("Diffusion Pseudotime", fontweight="bold")
    ax.set_ylabel("Mean Expression (smoothed)", fontweight="bold")
    ax.legend(frameon=True, fancybox=True, shadow=True, facecolor="white",
              edgecolor="#eee", loc="upper right", fontsize=10)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="y", ls=":", alpha=0.4)

    plt.tight_layout()
    outpath = OUT / "04c_tcell_pseudotime_trajectories.png"
    plt.savefig(outpath, dpi=300, facecolor=fig.get_facecolor())
    plt.close()
    print(f"  Saved → {outpath}")


if __name__ == "__main__":
    plot_tcell_trajectory()
