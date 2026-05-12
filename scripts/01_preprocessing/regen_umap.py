#!/usr/bin/env python
"""Regenerate the 4-panel integration UMAP from the saved atlas,
plus cell-type annotation UMAPs and per-marker feature plots."""
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from pathlib import Path
from scipy.sparse import issparse

FIG = Path("out/01_preprocessing/figures")
FIG.mkdir(parents=True, exist_ok=True)

# ── Cell-type markers (shared definition) ──────────────────────────
CELL_TYPE_MARKERS = {
    "CD8_T":        ["CD8A", "CD8B", "CD3E"],
    "CD4_T":        ["CD4", "IL7R", "CD3E"],
    "Treg":         ["FOXP3", "IL2RA", "CTLA4"],
    "NK":           ["NKG7", "GNLY", "KLRD1", "NCAM1"],
    "B_cell":       ["CD79A", "MS4A1", "CD19"],
    "Plasma":       ["JCHAIN", "MZB1", "SDC1"],
    "Macrophage":   ["CD68", "CD163", "C1QA", "CSF1R"],
    "mregDC":       ["LAMP3", "CCR7", "FSCN1"],
    "cDC1":         ["CLEC9A", "XCR1", "BATF3"],
    "cDC2":         ["CD1C", "CLEC10A", "FCER1A"],
    "Hepatocyte":   ["ALB", "APOB", "HP", "TF"],
    "Endothelial":  ["PECAM1", "CDH5", "VWF"],
    "Fibroblast":   ["COL1A1", "COL3A1", "DCN", "LUM"],
    "Tumor":        ["GPC3", "AFP", "GOLM1", "EPCAM"],
}

# ── Custom palettes — full subtype coverage ────────────────────────
CELLTYPE_PALETTE = {
    # CD8 T cells (red family)
    "CD8_T":                          "#e6194b",
    "CD8_Dysfunctional-terminal":   "#c62828",
    "CD8_Dysfunctional-progenitor": "#d32f2f",
    "CD8_Dysfunctional-Proliferating": "#e53935",
    "CD8_Dysfunctional-effector":   "#ef5350",
    "CD8_Effector":                 "#f44336",
    "CD8_Memory":                   "#ff5252",
    "CD8_Cytotoxic":                "#ff1744",
    # CD4 T cells (green family)
    "CD4_T":                          "#2e7d32",
    "CD4_Tfh-like":                 "#43a047",
    "CD4_Th1 GZMK":                "#66bb6a",
    "CD4_Th1/Th17":                 "#81c784",
    "CD4_Effector-Memory":          "#4caf50",
    "CD4_Memory":                   "#388e3c",
    # Tregs (yellow-amber)
    "Treg":                           "#f9a825",
    "Treg_Treg":                      "#fbc02d",
    # NK cells (blue family)
    "NK":                             "#1565c0",
                # Unconventional T (teal-cyan)
    "MAIT":                           "#00838f",
    "iNKT":                           "#00acc1",
    "gdT":                            "#26c6da",
    "ILC":                            "#4dd0e1",
    "Naive_Naive":                    "#80deea",
    # B cells (orange family)
    "B_cell":                         "#e65100",
    "B_Memory":                       "#ef6c00",
    "B_Naive":                        "#f57c00",
    "B_Plasma":                       "#fb8c00",
    "Plasma":                         "#ff9800",
    # Macrophages/Monocytes (purple family)
    "Macrophage":                     "#6a1b9a",
    "Macrophage_I":                   "#7b1fa2",
    "Macrophage_II":                  "#8e24aa",
    "Macrophage_FCN1/CXCL9":          "#9c27b0",
    "Macrophage_FOLR2":               "#ab47bc",
    "Macrophage_TREM2":               "#ba68c8",
    "Macrophage_Kupffer":             "#ce93d8",
    "Macrophage_SPP1":                "#7c4dff",
    "Monocyte_CD14":                  "#5e35b1",
    "Monocyte_CD14 SPP1":             "#673ab7",
    "Monocyte_CD16":                  "#7e57c2",
    "Monocyte_Intermediate":          "#9575cd",
    # DCs (pink-magenta family)
    "mregDC":                         "#c2185b",
    "DC_mregDC":                      "#d81b60",
    "cDC1":                           "#e91e63",
    "DC_cDC1":                        "#ec407a",
    "cDC2":                           "#f06292",
    "DC_cDC2":                        "#f48fb1",
    "DC_pDC":                         "#f8bbd0",
    # Stromal / Structural (earth tones)
    "Hepatocyte":                     "#00695c",
    "Stromal_Hepatocytes":            "#00897b",
    "Endothelial":                    "#7986cb",
    "Stromal_Endothelial cells":      "#9fa8da",
    "Fibroblast":                     "#8d6e63",
    "Stromal_HSC":                    "#a1887f",
    "Tumor":                          "#b71c1c",
    "Mast":                           "#ff6f00",
    # Proliferating / Other
    "Proliferating_Proliferating":    "#ffd600",
    "Unknown":                        "#bdbdbd",
    "Unknown_Doublets":               "#9e9e9e",
    "Unknown_Empty":                  "#e0e0e0",
}
MAGMA_DARK = LinearSegmentedColormap.from_list(
    "magma_dark", ["#000004", "#3b0f70", "#8c2981", "#de4968", "#fca50a", "#fcffa4"]
)

print("Loading atlas (backed)...")
adata = ad.read_h5ad("data/processed/integrated_io_atlas.h5ad", backed="r")
umap = np.array(adata.obsm["X_umap"])
obs = adata.obs.copy()
print(f"  {len(obs):,} cells, UMAP shape: {umap.shape}")

fig, axes = plt.subplots(2, 2, figsize=(14, 13))

# Panel A: Dataset
ax = axes[0, 0]
ds_list = sorted(obs["dataset"].dropna().unique())
cmap_ds = plt.cm.Set2(np.linspace(0, 1, len(ds_list)))
for i, ds in enumerate(ds_list):
    idx = np.where(obs["dataset"].values == ds)[0]
    if len(idx) > 0:
        ax.scatter(umap[idx, 0], umap[idx, 1],
                   s=0.5, alpha=0.2, c=[cmap_ds[i]], label=ds, rasterized=True)
        print(f"  A: {ds} -> {len(idx):,} cells")
ax.legend(fontsize=7, markerscale=10)
ax.set_title("A. Dataset", fontweight="bold")
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

# Panel B: Study
ax = axes[0, 1]
studies = sorted(obs["study"].dropna().unique())
cmap_st = plt.cm.tab10(np.linspace(0, 1, len(studies)))
for i, st in enumerate(studies):
    idx = np.where(obs["study"].values == st)[0]
    if len(idx) > 0:
        ax.scatter(umap[idx, 0], umap[idx, 1],
                   s=0.5, alpha=0.2, c=[cmap_st[i]], label=st, rasterized=True)
        print(f"  B: {st} -> {len(idx):,} cells")
ax.legend(fontsize=8, markerscale=10)
ax.set_title("B. Study", fontweight="bold")
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

# Panel C: Tissue
ax = axes[1, 0]
tissues = sorted(obs["tissue"].dropna().unique())
cmap_ti = plt.cm.Paired(np.linspace(0, 1, max(len(tissues), 2)))
for i, ti in enumerate(tissues):
    idx = np.where(obs["tissue"].values == ti)[0]
    if len(idx) > 0:
        ax.scatter(umap[idx, 0], umap[idx, 1],
                   s=0.5, alpha=0.2, c=[cmap_ti[i]], label=ti, rasterized=True)
        print(f"  C: {ti} -> {len(idx):,} cells")
ax.legend(fontsize=8, markerscale=10)
ax.set_title("C. Tissue", fontweight="bold")
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

# Panel D: IO Response
ax = axes[1, 1]
resp_colors = {"Responder": "#0077BB", "Non-Responder": "#CC3311"}
for resp, color in resp_colors.items():
    idx = np.where(obs["response"].values == resp)[0]
    if len(idx) > 0:
        ax.scatter(umap[idx, 0], umap[idx, 1],
                   s=0.5, alpha=0.2, c=color, label=f"{resp} ({len(idx):,})",
                   rasterized=True)
        print(f"  D: {resp} -> {len(idx):,} cells")
idx_na = np.where(obs["response"].isna().values)[0]
if len(idx_na) > 0:
    ax.scatter(umap[idx_na, 0], umap[idx_na, 1],
               s=0.3, alpha=0.05, c="gray", label=f"No label ({len(idx_na):,})",
               rasterized=True)
    print(f"  D: No label -> {len(idx_na):,} cells")
ax.legend(fontsize=8, markerscale=10)
ax.set_title("D. IO Response", fontweight="bold")
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

plt.suptitle(f"Integrated HCC IO Atlas \u2014 {len(obs):,} cells\n"
             f"Harmony batch correction, {len(ds_list)} datasets",
             fontweight="bold", fontsize=13)
plt.tight_layout()
for ext in ["png", "pdf"]:
    plt.savefig(FIG / f"01_integration_umap.{ext}", dpi=140, bbox_inches="tight")
plt.close()
print(f"\nSaved to {FIG / '01_integration_umap.png'}")

# ═══════════════════════════════════════════════════════════════════════
# FIGURE 2 — Cell-type annotation UMAP with density contours
# ═══════════════════════════════════════════════════════════════════════
print("\n[CELLTYPE UMAP] Generating cell-type annotation figure...")

ct_col = None
for candidate in ["cell_type", "cell_type_original", "leiden"]:
    if candidate in obs.columns:
        ct_col = candidate
        break
if ct_col is None:
    print("  [WARN] No cell type column found — skipping annotation UMAP")
else:
    ct_series = obs[ct_col].astype(str).fillna("Unknown")
    top_cts = ct_series.value_counts().head(14).index.tolist()

    fig = plt.figure(figsize=(20, 10))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.4, 1], wspace=0.25)

    # Left: Full cell-type scatter with KDE contours per major type
    ax_main = fig.add_subplot(gs[0])
    ax_main.scatter(umap[:, 0], umap[:, 1], s=0.08, c="#f0f0f0", alpha=0.15, rasterized=True)
    for ct in top_cts:
        idx = np.where(ct_series.values == ct)[0]
        color = CELLTYPE_PALETTE.get(ct, plt.cm.Spectral(hash(ct) % 256 / 256.0))
        ax_main.scatter(umap[idx, 0], umap[idx, 1],
                        s=0.6, alpha=0.35, c=[color] if isinstance(color, str) else [color],
                        label=f"{ct} ({len(idx):,})", rasterized=True)
        if len(idx) > 300:
            try:
                from scipy.stats import gaussian_kde
                xy = umap[idx].T
                kde = gaussian_kde(xy, bw_method=0.3)
                xg = np.linspace(xy[0].min(), xy[0].max(), 60)
                yg = np.linspace(xy[1].min(), xy[1].max(), 60)
                xx, yy = np.meshgrid(xg, yg)
                zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
                ax_main.contour(xx, yy, zz, levels=3, colors=[color] if isinstance(color, str) else ["#555"],
                                linewidths=0.6, alpha=0.5)
            except Exception:
                pass
    ax_main.legend(fontsize=7, markerscale=12, ncol=2,
                   loc="upper left", framealpha=0.85, edgecolor="#ccc")
    ax_main.set_title("Cell-Type Annotations", fontweight="bold", fontsize=13)
    ax_main.set_xlabel("UMAP 1"); ax_main.set_ylabel("UMAP 2")
    for spine in ax_main.spines.values():
        spine.set_visible(False)
    ax_main.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    # Right: Stacked horizontal bar showing composition
    ax_bar = fig.add_subplot(gs[1])
    ct_counts = ct_series.value_counts()
    ct_top = ct_counts.head(14)
    bars = ax_bar.barh(range(len(ct_top)), ct_top.values,
                       color=[CELLTYPE_PALETTE.get(c, "#999") for c in ct_top.index],
                       edgecolor="white", linewidth=0.5, height=0.7)
    ax_bar.set_yticks(range(len(ct_top)))
    ax_bar.set_yticklabels(ct_top.index, fontsize=9)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("Number of Cells", fontsize=10)
    ax_bar.set_title("Cell-Type Abundance", fontweight="bold", fontsize=12)
    for i, (ct, n) in enumerate(ct_top.items()):
        ax_bar.text(n + ct_top.max() * 0.01, i, f"{n:,}", va="center", fontsize=8)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    plt.suptitle(f"HCC IO Atlas Cell-Type Map — {len(obs):,} cells",
                 fontweight="bold", fontsize=14, y=1.01)
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"01_celltype_annotation_umap.{ext}", dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Saved {FIG / '01_celltype_annotation_umap.png'}")

# ═══════════════════════════════════════════════════════════════════════
# FIGURE 3 — Per-marker feature plots (dark-background magma style)
# ═══════════════════════════════════════════════════════════════════════
print("\n[MARKER FEATURES] Generating per-cell-type marker feature plots...")

# Determine available gene names
var_names = set(adata.var_names)
raw_available = adata.raw is not None
if raw_available:
    var_names_raw = set(adata.raw.var_names)
    _raw_idx_map = {g: i for i, g in enumerate(adata.raw.var_names)}
else:
    var_names_raw = var_names
    _raw_idx_map = {}
_var_idx_map = {g: i for i, g in enumerate(adata.var_names)}


def _get_gene_expr(gene):
    """Extract a 1-D expression array for *gene*."""
    if raw_available and gene in _raw_idx_map:
        idx = _raw_idx_map[gene]
        x = adata.raw.X[:, idx]
    elif gene in _var_idx_map:
        idx = _var_idx_map[gene]
        x = adata.X[:, idx]
    else:
        return None
    if issparse(x):
        x = np.asarray(x.todense()).ravel()
    else:
        x = np.asarray(x).ravel()
    return x


# Build list of panels: one representative marker per cell type
rep_markers = {}
for ct, markers in CELL_TYPE_MARKERS.items():
    for m in markers:
        if m in var_names or m in var_names_raw:
            rep_markers[ct] = m
            break

n_panels = len(rep_markers)
if n_panels == 0:
    print("  [WARN] No marker genes found in atlas — skipping feature plots")
else:
    ncols = 4
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows),
                             facecolor="#0e0e1a")
    axes = np.atleast_2d(axes)
    axf = axes.flatten()

    for i, (ct, gene) in enumerate(rep_markers.items()):
        ax = axf[i]
        ax.set_facecolor("#0e0e1a")
        expr = _get_gene_expr(gene)
        if expr is None:
            ax.set_visible(False)
            continue
        # Clip top 1% for better contrast
        vmax = np.percentile(expr[expr > 0], 99) if (expr > 0).any() else 1.0
        expr_clip = np.clip(expr, 0, vmax)
        # Sort so high-expression cells are drawn on top
        order = np.argsort(expr_clip)
        sc_pt = ax.scatter(umap[order, 0], umap[order, 1],
                           c=expr_clip[order], cmap=MAGMA_DARK,
                           s=0.4, alpha=0.6, vmin=0, vmax=vmax, rasterized=True)
        ax.set_title(f"{ct}\n{gene}", color="white", fontweight="bold", fontsize=10)
        ax.tick_params(colors="#555", labelsize=0)
        for spine in ax.spines.values():
            spine.set_edgecolor("#333")
        cb = plt.colorbar(sc_pt, ax=ax, fraction=0.04, pad=0.02)
        cb.ax.tick_params(colors="#aaa", labelsize=7)
        cb.outline.set_edgecolor("#555")

    # Hide empty axes
    for j in range(i + 1, len(axf)):
        axf[j].set_visible(False)

    fig.suptitle("Cell-Type Marker Expression — Feature Plots",
                 color="white", fontweight="bold", fontsize=15, y=1.01)
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(FIG / f"01_marker_feature_plots.{ext}", dpi=180,
                    bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    print(f"  Saved {FIG / '01_marker_feature_plots.png'}")

print("\n[DONE] All UMAP figures generated.")

# ═══════════════════════════════════════════════════════════════════════
# FIGURE 4 — Immune Major and Subtypes UMAPs
# ═══════════════════════════════════════════════════════════════════════
print("\n[IMMUNE UMAPs] Generating Major and Subtype figures...")

import pandas as pd
clust_ann = pd.read_csv("data/processed/GSE206325_cluster_annotations.csv")
clust_ann["cluster_ID"] = clust_ann["cluster_ID"].astype(str)

major_map = {}
sub_map = {}
for _, row in clust_ann.iterrows():
    cid = row["cluster_ID"]
    comp = row["compartment"]
    group = row["group"]
    type_ = row["type"] if pd.notna(row["type"]) else ""
    subg = row["subgroup"] if pd.notna(row["subgroup"]) else ""
    
    if comp in ["Lymphoid", "Myeloid"]:
        major_map[cid] = group
        st = f"{type_}_{subg}".strip("_")
        if not st: st = type_
        if not st: st = subg
        if not st: st = group
        sub_map[cid] = st
    else:
        major_map[cid] = "Other"
        sub_map[cid] = "Other"

obs["immune_major"] = obs["cluster_id"].map(major_map).fillna("Unknown")
obs["immune_subtype"] = obs["cluster_id"].map(sub_map).fillna("Unknown")

import scanpy as sc
adata_pl = ad.AnnData(obs=obs, obsm={"X_umap": umap})

fig, axes = plt.subplots(1, 2, figsize=(16, 7))


adata_pl.obs["immune_major"] = adata_pl.obs["immune_major"].astype("category")
adata_pl.obs["immune_subtype"] = adata_pl.obs["immune_subtype"].astype("category")

if "Other" in adata_pl.obs["immune_major"].cat.categories:
    import matplotlib.pyplot as plt
    colors = plt.cm.tab20.colors
    c_list = [colors[i % len(colors)] for i in range(len(adata_pl.obs["immune_major"].cat.categories))]
    c_list[list(adata_pl.obs["immune_major"].cat.categories).index("Other")] = "#E0E0E0"
    adata_pl.uns["immune_major_colors"] = c_list

if "Other" in adata_pl.obs["immune_subtype"].cat.categories:
    import matplotlib.pyplot as plt
    colors = plt.cm.tab20.colors
    c_list = [colors[i % len(colors)] for i in range(len(adata_pl.obs["immune_subtype"].cat.categories))]
    c_list[list(adata_pl.obs["immune_subtype"].cat.categories).index("Other")] = "#E0E0E0"
    adata_pl.uns["immune_subtype_colors"] = c_list

sc.pl.umap(adata_pl, color="immune_major", ax=axes[0], show=False, title="A. Major Immune Cell Types")
sc.pl.umap(adata_pl, color="immune_subtype", ax=axes[1], show=False, title="B. Immune Subtypes")


plt.tight_layout()
fig.savefig(FIG / "01_immune_umaps.png", dpi=200, bbox_inches="tight")
plt.close()
print("  Saved", FIG / "01_immune_umaps.png")

