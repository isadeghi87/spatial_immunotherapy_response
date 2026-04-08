#!/usr/bin/env python3
"""02g – InferCNV analysis on the integrated IO atlas.

Identifies malignant clones via inferred copy-number variation using
immune reference cells, scores malignancy, and produces:
  1. Chromosome heatmap grouped by malignant status / response / cell type
  2. CNV score UMAP
  3. Multi-panel attractive dashboard (heatmap + violin)
"""
import warnings
warnings.filterwarnings("ignore")

import scanpy as sc
try:
    import infercnvpy as cnv
except ImportError:
    cnv = None
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mc

# ── paths ────────────────────────────────────────────────────────────
PROJECT  = Path(__file__).resolve().parents[2]
ATLAS    = PROJECT / "data" / "processed" / "integrated_io_atlas.h5ad"
GTF      = PROJECT / "data" / "references" / "gencode.v44.basic.annotation.gtf.gz"
OUT      = PROJECT / "out" / "02_analysis" / "figures"
OUT.mkdir(parents=True, exist_ok=True)


def _resolve_celltype_col(adata):
    for col in ("cell_type", "celltype_sub", "Cell.type", "celltype_global", "celltype", "leiden"):
        if col in adata.obs.columns:
            return col
    raise KeyError("No cell-type annotation column found in adata.obs")


def _pick_reference_cats(adata, ct_col):
    """Select ALL immune / stromal cell types as reference for InferCNV.

    In HCC, only Tumor / Hepatocyte / epithelial cells should be in the
    'test' group.  Everything else (T cells, NK, B cells, macrophages,
    DCs, monocytes, fibroblasts, mast cells …) is normal reference.
    """
    available = set(str(c) for c in adata.obs[ct_col].dropna().unique())

    # Keywords that mark likely MALIGNANT / epithelial cells (test group)
    tumor_kw = ("Tumor", "Hepatocyte", "Epithelial", "Cancer", "Malig")

    # Everything that does NOT match a tumor keyword is reference
    ref_cats = sorted(
        c for c in available
        if not any(k.lower() in c.lower() for k in tumor_kw)
        and c not in ("Unknown", "Unknown_Doublets", "Unknown_Empty", "nan")
    )

    if len(ref_cats) < 2:
        raise ValueError(
            f"Could not find >= 2 immune/stromal reference categories. "
            f"Available: {available}"
        )

    print(f"  Reference categories ({len(ref_cats)}): {ref_cats[:10]} …")
    test_cats = sorted(available - set(ref_cats))
    print(f"  Test (potential malignant) categories ({len(test_cats)}): {test_cats}")
    return ref_cats


def run_infercnv(adata_path=None):
    if adata_path is None:
        adata_path = str(ATLAS)

    print("Loading atlas …")
    adata = sc.read_h5ad(adata_path)
    print(f"  Shape: {adata.shape}")

    ct_col = _resolve_celltype_col(adata)
    print(f"  Cell-type column: {ct_col}")

    ref_cats = _pick_reference_cats(adata, ct_col)
    print(f"  Reference categories: {ref_cats}")

    if cnv is None:
        raise ImportError("infercnvpy is required – install with: pip install infercnvpy")

    # ── Add genomic positions from GENCODE GTF ──────────────────────
    if "chromosome" not in adata.var.columns or "start" not in adata.var.columns:
        print("Adding genomic positions from GENCODE GTF …")
        cnv.io.genomic_position_from_gtf(str(GTF), adata, gtf_gene_id="gene_name")
        n_mapped = adata.var["chromosome"].notna().sum()
        print(f"  Mapped {n_mapped}/{adata.shape[1]} genes to genomic positions")
        # Drop genes without positions (infercnvpy requires them)
        mask = adata.var["chromosome"].notna()
        adata = adata[:, mask].copy()
        print(f"  After filtering: {adata.shape}")

    # ── Subsample for InferCNV (972k cells is too large) ─────────────
    MAX_CELLS = 50_000
    if adata.n_obs > MAX_CELLS:
        print(f"  Stratified subsampling {MAX_CELLS} cells (from {adata.n_obs}) …")
        rng = np.random.default_rng(42)
        # Stratify by cell_type so rare types (Tumor, Hepatocyte) are represented
        ct = adata.obs[ct_col].astype(str)
        per_type = max(50, MAX_CELLS // ct.nunique())
        keep = []
        for label in ct.unique():
            mask = np.where(ct.values == label)[0]
            n = min(len(mask), per_type)
            keep.append(rng.choice(mask, n, replace=False))
        idx = np.concatenate(keep)
        # If still above MAX_CELLS, trim; if below, top up randomly
        if len(idx) > MAX_CELLS:
            idx = rng.choice(idx, MAX_CELLS, replace=False)
        idx = np.sort(np.unique(idx))
        adata = adata[idx].copy()
        print(f"  Subsampled shape: {adata.shape}")
        ct_counts = adata.obs[ct_col].value_counts()
        print(f"  Top types: {ct_counts.head(10).to_dict()}")

    print("Running InferCNV (window=100, step=10) …")
    cnv.tl.infercnv(adata, reference_key=ct_col, reference_cat=ref_cats,
                     window_size=100, step=10)

    print("  PCA + neighbors + leiden on CNV space …")
    cnv.tl.pca(adata)
    cnv.pp.neighbors(adata)
    cnv.tl.leiden(adata, key_added="cnv_leiden")
    cnv.tl.cnv_score(adata)

    threshold = np.percentile(adata.obs["cnv_score"].dropna(), 85)
    adata.obs["malignant_status"] = np.where(
        adata.obs["cnv_score"] > threshold, "Malignant", "Normal"
    )
    print(f"  Malignant threshold (85th pct): {threshold:.4f}")
    print(adata.obs["malignant_status"].value_counts().to_string())

    groupby = ct_col

    # Fill NaN in groupby column
    adata.obs[ct_col] = adata.obs[ct_col].astype(str).replace("nan", "Unknown")

    # ── Plot 1: publication-quality chromosome heatmap ────────────────
    print("Plotting chromosome heatmap (publication style) …")
    _plot_publication_heatmap(adata, ct_col)

    # ── Plot 2: CNV score UMAP ───────────────────────────────────────
    print("Plotting CNV-score UMAP …")
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color="cnv_score", ax=ax, show=False, cmap="magma",
               title="CNV Score (Malignancy Proxy)")
    plt.tight_layout()
    plt.savefig(OUT / "02g_cnv_score_umap.png", dpi=300, bbox_inches="tight")
    plt.close()

    # ── Plot 3: attractive multi-panel dashboard ─────────────────────
    print("Generating attractive dashboard …")
    _plot_dashboard(adata, ct_col)

    # ── Save annotated atlas with CNV columns for downstream (02h) ──
    cnv_out = PROJECT / "data" / "processed" / "integrated_io_atlas_cnv.h5ad"
    print(f"Saving CNV-annotated atlas → {cnv_out}")
    adata.write_h5ad(str(cnv_out))

    print("Done – all InferCNV figures saved to", OUT)


def _major_group(ct):
    """Map fine cell types to broad biological groups for display."""
    ct_l = str(ct).lower()
    if any(k in ct_l for k in ("tumor", "cancer", "malig")):
        return "Tumor"
    if "hepatocyte" in ct_l or "stromal_hepatocyte" in ct_l:
        return "Hepatocyte"
    if any(k in ct_l for k in ("cd8", "cd4", "treg", "mait", "naive_naive",
                                "inkt", "gdt", "proliferating_prolif")):
        return "T / NK-T cells"
    if "nk" in ct_l or "ilc" in ct_l:
        return "NK / ILC"
    if any(k in ct_l for k in ("b_", "b cell", "plasma")):
        return "B / Plasma"
    if any(k in ct_l for k in ("macro", "mono", "kupffer")):
        return "Myeloid"
    if any(k in ct_l for k in ("dc", "mreg", "cdc")):
        return "Dendritic cells"
    if any(k in ct_l for k in ("endothelial", "fibroblast", "stromal", "hsc")):
        return "Stromal"
    if any(k in ct_l for k in ("mast",)):
        return "Mast"
    return "Other"


# Ordered groups for display (tumor/hepatocyte first, then immune, then stromal)
_GROUP_ORDER = [
    "Tumor", "Hepatocyte", "T / NK-T cells", "NK / ILC",
    "B / Plasma", "Myeloid", "Dendritic cells", "Mast", "Stromal", "Other",
]

_GROUP_COLORS = {
    "Tumor":            "#8B1A1A",
    "Hepatocyte":       "#CD853F",
    "T / NK-T cells":   "#2E8B57",
    "NK / ILC":         "#4682B4",
    "B / Plasma":       "#6A5ACD",
    "Myeloid":          "#D2691E",
    "Dendritic cells":  "#DAA520",
    "Mast":             "#DB7093",
    "Stromal":          "#708090",
    "Other":            "#A9A9A9",
}


def _get_chr_boundaries(adata):
    """Compute chromosome boundaries from adata.var for the CNV bins."""
    var = adata.var[["chromosome", "start"]].copy()
    var["chromosome"] = var["chromosome"].astype(str)
    # Filter to standard autosomes
    auto = [f"chr{i}" for i in range(1, 23)]
    var = var[var["chromosome"].isin(auto)]
    # infercnvpy orders genes by chromosome-then-position and bins them
    # We approximate bin boundaries from the var ordering
    chr_order = {f"chr{i}": i for i in range(1, 23)}
    var["chr_num"] = var["chromosome"].map(chr_order)
    var = var.sort_values(["chr_num", "start"])
    # Count genes per chromosome → proportional bin count
    gene_counts = var.groupby("chr_num").size()
    total_genes = gene_counts.sum()
    n_bins = adata.obsm["X_cnv"].shape[1] if "X_cnv" in adata.obsm else 1265
    boundaries = {}
    pos = 0
    for ch in range(1, 23):
        n = int(gene_counts.get(ch, 0) / total_genes * n_bins)
        boundaries[ch] = (pos, pos + n)
        pos += n
    return boundaries, n_bins


def _plot_publication_heatmap(adata, ct_col):
    """Publication-quality CNV heatmap grouped by cell type with chromosome
    annotations, inspired by Tirosh / Patel / Filbin style CNV figures."""
    from scipy.stats import mannwhitneyu

    # ── Assign major groups ──────────────────────────────────────────
    adata.obs["_cnv_group"] = adata.obs[ct_col].astype(str).map(_major_group)

    # Keep only groups that actually exist
    present = [g for g in _GROUP_ORDER if g in adata.obs["_cnv_group"].values]

    # ── Build ordered index by group → cnv_score within group ────────
    obs = adata.obs.copy()
    obs["_cnv_group"] = pd.Categorical(obs["_cnv_group"], categories=present, ordered=True)
    obs = obs.sort_values(["_cnv_group", "cnv_score"], ascending=[True, False])
    order = obs.index

    # ── Get CNV matrix ───────────────────────────────────────────────
    cnv_mat = adata[order].obsm["X_cnv"]
    if hasattr(cnv_mat, "toarray"):
        cnv_mat = cnv_mat.toarray()
    cnv_mat = np.asarray(cnv_mat, dtype=float)

    # Downsample rows for plotting if too many
    n_cells = cnv_mat.shape[0]
    max_rows = 4000
    if n_cells > max_rows:
        step = n_cells // max_rows
        row_idx = np.arange(0, n_cells, step)
        cnv_mat = cnv_mat[row_idx]
        obs = obs.iloc[row_idx]
    else:
        row_idx = np.arange(n_cells)

    n_bins = cnv_mat.shape[1]

    # ── Chromosome boundaries ────────────────────────────────────────
    chr_bounds, _ = _get_chr_boundaries(adata)

    # ── Figure layout ────────────────────────────────────────────────
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(
        2, 3,
        height_ratios=[3, 1],
        width_ratios=[0.03, 1, 0.03],
        hspace=0.25, wspace=0.04,
    )

    # ── A: Chromosome position bar (top) ─────────────────────────────
    ax_chr = fig.add_subplot(gs[0, 1])

    # Brown/teal CNV colormap (like the reference image)
    cnv_cmap = LinearSegmentedColormap.from_list(
        "cnv_pub", ["#1B4F72", "#5DADE2", "#F0F0F0", "#E59866", "#784212"], N=256
    )
    vmax = np.percentile(np.abs(cnv_mat), 99)
    im = ax_chr.imshow(
        cnv_mat, aspect="auto", cmap=cnv_cmap, vmin=-vmax, vmax=vmax,
        interpolation="none",
    )

    # Chromosome boundary lines + labels at top
    chr_centers = []
    for ch in range(1, 23):
        if ch not in chr_bounds:
            continue
        lo, hi = chr_bounds[ch]
        if ch > 1:
            ax_chr.axvline(x=lo, color="black", linewidth=0.6, linestyle="--", alpha=0.5)
        chr_centers.append((ch, (lo + hi) / 2))

    ax2 = ax_chr.twiny()
    ax2.set_xlim(ax_chr.get_xlim())
    ax2.set_xticks([c for _, c in chr_centers])
    ax2.set_xticklabels([str(c) for c, _ in chr_centers], fontsize=8, fontweight="bold")
    ax2.tick_params(length=0)
    ax2.set_xlabel("Chromosome position", fontsize=10, fontweight="bold", labelpad=8)
    ax_chr.set_xticks([])
    ax_chr.set_yticks([])

    # ── Left sidebar: cell type group color bar ──────────────────────
    ax_sidebar = fig.add_subplot(gs[0, 0])
    group_vals = obs["_cnv_group"].values
    group_numeric = np.array([present.index(g) if g in present else -1 for g in group_vals])
    sidebar_cmap = mc.ListedColormap([_GROUP_COLORS.get(g, "#999") for g in present])
    ax_sidebar.imshow(
        group_numeric.reshape(-1, 1), aspect="auto", cmap=sidebar_cmap,
        interpolation="none", vmin=0, vmax=len(present) - 1,
    )
    ax_sidebar.set_xticks([])
    ax_sidebar.set_yticks([])

    # Group name labels at the midpoint of each group band
    for g in present:
        rows = np.where(group_vals == g)[0]
        if len(rows) == 0:
            continue
        mid = (rows[0] + rows[-1]) / 2
        ax_sidebar.text(
            -0.3, mid, g, ha="right", va="center", fontsize=8,
            fontweight="bold", color=_GROUP_COLORS.get(g, "k"),
            transform=ax_sidebar.get_yaxis_transform(),
        )
        # Horizontal separator between groups
        if rows[-1] < len(group_vals) - 1:
            ax_chr.axhline(y=rows[-1] + 0.5, color="white", linewidth=1.5)

    # ── Right sidebar: colorbar ──────────────────────────────────────
    ax_cb = fig.add_subplot(gs[0, 2])
    cb = plt.colorbar(im, cax=ax_cb)
    cb.set_label("Inferred CNV", fontweight="bold", fontsize=9)

    # ── B: Violin of CNA burden per group ────────────────────────────
    ax_vio = fig.add_subplot(gs[1, 1])
    group_data = []
    group_labels = []
    group_colors = []
    for g in present:
        mask = adata.obs["_cnv_group"] == g
        vals = adata.obs.loc[mask, "cnv_score"].dropna().values
        if len(vals) > 0:
            group_data.append(vals)
            group_labels.append(g)
            group_colors.append(_GROUP_COLORS.get(g, "#999"))

    if group_data:
        parts = ax_vio.violinplot(
            group_data, showmeans=False, showmedians=False, showextrema=False
        )
        for pc, c in zip(parts["bodies"], group_colors):
            pc.set_facecolor(c)
            pc.set_edgecolor("black")
            pc.set_alpha(0.85)
            pc.set_linewidth(0.5)
        for i, d in enumerate(group_data):
            q1, med, q3 = np.percentile(d, [25, 50, 75])
            ax_vio.scatter(i + 1, med, marker="o", color="white", s=20, zorder=3,
                           edgecolors="black", linewidths=0.5)
            ax_vio.vlines(i + 1, q1, q3, color="k", lw=3)
            ax_vio.vlines(i + 1, d.min(), d.max(), color="k", lw=0.7)

        ax_vio.set_xticks(range(1, len(group_labels) + 1))
        ax_vio.set_xticklabels(group_labels, rotation=35, ha="right", fontsize=8,
                                fontweight="bold")

        # Statistical tests against first group (Tumor if present)
        ref_idx = 0
        for i in range(1, len(group_data)):
            if len(group_data[ref_idx]) > 3 and len(group_data[i]) > 3:
                _, p = mannwhitneyu(group_data[ref_idx], group_data[i], alternative="two-sided")
                if p < 0.001:
                    stars = f"p={p:.1e}"
                elif p < 0.05:
                    stars = f"p={p:.3f}"
                else:
                    stars = "ns"
                y_max = max(np.percentile(group_data[ref_idx], 95),
                            np.percentile(group_data[i], 95))
                ax_vio.text(i + 1, y_max * 1.02, stars, ha="center", fontsize=6,
                            color=group_colors[i])

        ax_vio.set_ylabel("CNV Score\n(CNA Burden)", fontweight="bold", fontsize=9)
        ax_vio.set_title("Relative CNA Burden by Cell Type Group", fontweight="bold",
                         fontsize=11, pad=10)
        for sp in ("top", "right"):
            ax_vio.spines[sp].set_visible(False)

    fig.savefig(OUT / "02g_infercnv_chromosome_heatmap.png", dpi=250, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {OUT / '02g_infercnv_chromosome_heatmap.png'}")


def _plot_dashboard(adata, ct_col):
    """Compact 2-panel dashboard: heatmap + cell-type violin."""
    plt.style.use("seaborn-v0_8-white")
    fig = plt.figure(figsize=(20, 7))
    gs = gridspec.GridSpec(1, 3, width_ratios=[0.06, 2, 1], wspace=0.12)

    # Assign groups
    adata.obs["_cnv_group"] = adata.obs[ct_col].astype(str).map(_major_group)
    present = [g for g in _GROUP_ORDER if g in adata.obs["_cnv_group"].values]

    order = adata.obs.sort_values(["_cnv_group", "cnv_score"], ascending=[True, False]).index
    obs_sorted = adata.obs.loc[order]

    n = len(obs_sorted)
    step = max(1, n // 3000)
    idx = obs_sorted.index[::step]
    obs_sub = obs_sorted.loc[idx]

    # ── sidebar: cell type group ─────────────────────────────────────
    ax_annot = fig.add_subplot(gs[0])
    group_vals = obs_sub["_cnv_group"].values
    group_numeric = np.array([present.index(g) if g in present else -1 for g in group_vals])
    sidebar_cmap = mc.ListedColormap([_GROUP_COLORS.get(g, "#999") for g in present])
    ax_annot.imshow(
        group_numeric.reshape(-1, 1), aspect="auto", cmap=sidebar_cmap,
        interpolation="none", vmin=0, vmax=len(present) - 1,
    )
    ax_annot.set_xticks([0])
    ax_annot.set_xticklabels(["Cell Group"], rotation=90, fontsize=9)
    ax_annot.set_yticks([])
    # Label groups
    for g in present:
        rows = np.where(group_vals == g)[0]
        if len(rows):
            mid = (rows[0] + rows[-1]) / 2
            ax_annot.text(-0.5, mid, g, ha="right", va="center", fontsize=7,
                          fontweight="bold", color=_GROUP_COLORS.get(g, "k"),
                          transform=ax_annot.get_yaxis_transform())

    # ── heatmap ──────────────────────────────────────────────────────
    ax_main = fig.add_subplot(gs[1])
    cnv_key = "X_cnv" if "X_cnv" in adata.obsm else None
    if cnv_key is None:
        ax_main.text(0.5, 0.5, "X_cnv not found", ha="center", va="center",
                     transform=ax_main.transAxes, fontsize=14)
    else:
        cnv_mat = adata[idx].obsm[cnv_key]
        if hasattr(cnv_mat, 'toarray'):
            cnv_mat = cnv_mat.toarray()
        cnv_mat = np.asarray(cnv_mat, dtype=float)
        n_bins = cnv_mat.shape[1]
        cmap_cnv = LinearSegmentedColormap.from_list(
            "cnv", ["#1B4F72", "#5DADE2", "#F0F0F0", "#E59866", "#784212"], N=256
        )
        vmax = np.percentile(np.abs(cnv_mat), 99)
        im = ax_main.imshow(cnv_mat, aspect="auto", cmap=cmap_cnv, vmin=-vmax, vmax=vmax,
                            interpolation="none")

        chr_bounds, _ = _get_chr_boundaries(adata)
        for ch in range(1, 23):
            if ch in chr_bounds and ch > 1:
                ax_main.axvline(x=chr_bounds[ch][0], color="black", lw=0.5,
                                ls="--", alpha=0.4)

        chr_centers = [(ch, (lo + hi) / 2) for ch, (lo, hi) in chr_bounds.items()]
        ax_main.set_xticks([c for _, c in chr_centers])
        ax_main.set_xticklabels([str(c) for c, _ in chr_centers], fontsize=8)
        ax_main.set_xlim(0, n_bins)
        cb = plt.colorbar(im, ax=ax_main, pad=0.02, shrink=0.7)
        cb.set_label("Inferred CNV", fontweight="bold")

        # White separators between groups
        for g in present:
            gv = np.where(group_vals == g)[0]
            if len(gv) and gv[-1] < len(group_vals) - 1:
                ax_main.axhline(y=gv[-1] + 0.5, color="white", lw=1.5)

    ax_main.set_yticks([])
    ax_main.set_title("A. InferCNV Chromosomal Heatmap", fontweight="bold", pad=12)
    ax_main.set_xlabel("Chromosome", fontweight="bold")

    # ── violin: CNA burden per group ─────────────────────────────────
    ax_v = fig.add_subplot(gs[2])
    group_data = []
    group_labels = []
    group_cols = []
    for g in present:
        mask = adata.obs["_cnv_group"] == g
        vals = adata.obs.loc[mask, "cnv_score"].dropna().values
        if len(vals) > 0:
            group_data.append(vals)
            group_labels.append(g)
            group_cols.append(_GROUP_COLORS.get(g, "#999"))
    if group_data:
        parts = ax_v.violinplot(group_data, showmeans=False, showmedians=False,
                                showextrema=False)
        for pc, c in zip(parts["bodies"], group_cols):
            pc.set_facecolor(c); pc.set_edgecolor("black"); pc.set_alpha(0.85)
        for i, d in enumerate(group_data):
            q1, med, q3 = np.percentile(d, [25, 50, 75])
            ax_v.scatter(i + 1, med, marker="o", color="white", s=25, zorder=3,
                         edgecolors="black", linewidths=0.5)
            ax_v.vlines(i + 1, q1, q3, color="#333", lw=4)
            ax_v.vlines(i + 1, d.min(), d.max(), color="#333", lw=0.7)
        ax_v.set_xticks(range(1, len(group_labels) + 1))
        ax_v.set_xticklabels(group_labels, rotation=35, ha="right", fontsize=7,
                              fontweight="bold")
    ax_v.set_ylabel("CNV Score", fontweight="bold")
    ax_v.set_title("B. CNA Burden by Cell Group", fontweight="bold", pad=12)
    for sp in ("top", "right"):
        ax_v.spines[sp].set_visible(False)

    fig.savefig(OUT / "02g_infercnv_attractive_dashboard.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved dashboard → {OUT / '02g_infercnv_attractive_dashboard.png'}")


def replot_only():
    """Reload saved CNV atlas and regenerate all plots (skip InferCNV)."""
    cnv_path = PROJECT / "data" / "processed" / "integrated_io_atlas_cnv.h5ad"
    print(f"Loading saved CNV atlas → {cnv_path}")
    adata = sc.read_h5ad(str(cnv_path))
    print(f"  Shape: {adata.shape}")

    ct_col = _resolve_celltype_col(adata)
    adata.obs[ct_col] = adata.obs[ct_col].astype(str).replace("nan", "Unknown")

    print("Plotting publication heatmap …")
    _plot_publication_heatmap(adata, ct_col)

    print("Plotting CNV-score UMAP …")
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color="cnv_score", ax=ax, show=False, cmap="magma",
               title="CNV Score (Malignancy Proxy)")
    plt.tight_layout()
    plt.savefig(OUT / "02g_cnv_score_umap.png", dpi=300, bbox_inches="tight")
    plt.close()

    print("Plotting dashboard …")
    _plot_dashboard(adata, ct_col)
    print("Done – all plots regenerated.")


if __name__ == "__main__":
    import sys
    if "--replot" in sys.argv:
        replot_only()
    else:
        run_infercnv()
