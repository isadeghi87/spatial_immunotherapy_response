#!/usr/bin/env python3
"""
02b_visium_celltype_umap.py — Clustering, Cell-Type Annotation & UMAP for Visium ST

1. Merge all GSE238264 sections → PCA → harmony batch-correction → UMAP
2. Leiden clustering at multiple resolutions
3. Cell-type annotation via marker module scores (max-score assignment)
4. Figures:
   a) UMAP coloured by annotated cell type (all spots)
   b) UMAP coloured by sample / response
   c) Heatmap of cell-type marker genes (dot plot style)
   d) UMAP of immune-subtype spots only
   e) Spatial maps coloured by cell-type annotation
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from pathlib import Path
from collections import OrderedDict
from scipy.sparse import issparse

sc.settings.verbosity = 1

# ── Paths ──────────────────────────────────────────────────────────────
PROJECT = Path(__file__).resolve().parents[2]
DATA    = PROJECT / "data" / "spatial_data" / "GSE238264_HCC_IO_Visium"
OUT     = PROJECT / "out" / "02b_spatial_biomarkers"
FIG     = OUT / "figures"
TAB     = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

# ── High-contrast cell-type palette ───────────────────────────────────
CELLTYPE_PAL = {
    "CD8+ T cell":   "#E31A1C",   # strong red
    "CD4+ T cell":   "#FF7F00",   # orange
    "Treg":          "#B15928",   # brown
    "NK cell":       "#E7298A",   # magenta
    "B cell":        "#6A3D9A",   # deep purple
    "Macrophage":    "#1F78B4",   # strong blue
    "DC":            "#33A02C",   # green
    "Hepatocyte":    "#A6CEE3",   # light blue
    "Endothelial":   "#FDBF6F",   # light orange
    "Fibroblast":    "#CAB2D6",   # light purple
    "Tumor":         "#666666",   # dark grey
    "Unassigned":    "#CCCCCC",   # light grey
}

# Immune subtypes for subsetting
IMMUNE_TYPES = {"CD8+ T cell", "CD4+ T cell", "Treg", "NK cell",
                "B cell", "Macrophage", "DC"}

IMMUNE_SUBTYPE_PAL = {k: v for k, v in CELLTYPE_PAL.items() if k in IMMUNE_TYPES}

# Response palette
RESPONSE_PAL = {"Responder": "#2ecc71", "Non_Responder": "#e74c3c"}

# Sample palette (high-contrast for 7)
SAMPLE_PAL = {
    "HCC1R":  "#E31A1C",
    "HCC2R":  "#FF7F00",
    "HCC3R":  "#33A02C",
    "HCC4R":  "#1F78B4",
    "HCC5NR": "#6A3D9A",
    "HCC6NR": "#B15928",
    "HCC7NR": "#E7298A",
}

# ── Cell markers ──────────────────────────────────────────────────────
CELL_MARKERS = OrderedDict({
    "CD8+ T cell":  ["CD8A", "CD8B", "GZMB", "PRF1", "GZMA", "GZMK", "NKG7"],
    "CD4+ T cell":  ["CD4", "IL7R", "CCR7", "LEF1", "TCF7"],
    "Treg":         ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"],
    "NK cell":      ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1"],
    "B cell":       ["CD79A", "MS4A1", "CD19", "CD79B", "IGHG1", "IGKC"],
    "Macrophage":   ["CD68", "CD163", "C1QA", "CSF1R", "MSR1", "MARCO"],
    "DC":           ["CLEC9A", "XCR1", "BATF3", "IRF8", "LAMP3", "FSCN1"],
    "Hepatocyte":   ["ALB", "APOB", "HP", "TF", "APOA1", "SERPINA1"],
    "Endothelial":  ["PECAM1", "CDH5", "VWF", "ENG", "PLVAP"],
    "Fibroblast":   ["COL1A1", "COL3A1", "DCN", "LUM", "ACTA2"],
    "Tumor":        ["GPC3", "AFP", "EPCAM", "KRT19", "SPINK1"],
})

# Immune-specific IO signatures for the immune UMAP
IMMUNE_SIGS = OrderedDict({
    "T_cell_exhaustion":   ["PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","TOX","ENTPD1"],
    "T_cell_cytotoxicity": ["GZMA","GZMB","GZMH","GZMK","PRF1","GNLY","NKG7"],
    "Progenitor_CD8":      ["TCF7","SLAMF6","LEF1","CCR7","IL7R","SELL"],
    "Terminal_exhaustion":  ["TOX","ENTPD1","HAVCR2","LAYN","TIGIT","CXCL13"],
    "Treg_signature":      ["FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"],
    "M1_macrophage":       ["NOS2","IL1B","TNF","CXCL9","CXCL10"],
    "M2_macrophage":       ["CD163","MRC1","MSR1","CD200R1","IL10"],
    "TLS_signature":       ["CXCL13","CXCR5","CR2","CD79A","MS4A1","CD19","IGHG1"],
    "IFN_gamma":           ["IFNG","STAT1","IRF1","CXCL9","CXCL10","CXCL11","IDO1"],
    "Antigen_presentation":["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","PSMB9"],
})


# ═══════════════════════════════════════════════════════════════════════
#  DATA LOADING & MERGING
# ═══════════════════════════════════════════════════════════════════════

SAMPLES = OrderedDict({
    "HCC1R": "Responder", "HCC2R": "Responder",
    "HCC3R": "Responder", "HCC4R": "Responder",
    "HCC5NR": "Non_Responder", "HCC6NR": "Non_Responder",
    "HCC7NR": "Non_Responder",
})


def load_and_merge():
    """Load each section, merge into one AnnData."""
    adatas = []
    for name, resp in SAMPLES.items():
        path = DATA / f"{name}.h5ad"
        if not path.exists():
            print(f"  [SKIP] {path} not found")
            continue
        print(f"  Loading {name} ({resp})...")
        ad = sc.read_h5ad(str(path))
        ad.var_names_make_unique()
        ad.obs_names = [f"{name}_{b}" for b in ad.obs_names]
        ad.obsm["spatial"] = ad.obs[["pxl_col", "pxl_row"]].values.astype(float)
        ad.obs["sample"] = name
        ad.obs["response"] = resp
        # Keep raw counts in layer
        ad.layers["counts"] = ad.X.copy()
        adatas.append(ad)
        print(f"    {ad.n_obs} spots")

    adata = adatas[0].concatenate(adatas[1:], batch_key="batch", join="inner")
    adata.obs["sample"] = adata.obs["sample"].astype(str)
    adata.obs["response"] = adata.obs["response"].astype(str)
    print(f"\n  Merged: {adata.n_obs} spots × {adata.n_vars} genes")
    return adata


# ═══════════════════════════════════════════════════════════════════════
#  PREPROCESSING → PCA → UMAP → LEIDEN
# ═══════════════════════════════════════════════════════════════════════

def preprocess_and_embed(adata):
    """Normalize, HVG, PCA, neighbors, UMAP, Leiden."""
    print("\n[1] Preprocessing & embedding...")

    # Normalise
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # store for later gene expression queries

    # HVG
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="sample")
    print(f"  HVGs: {adata.var['highly_variable'].sum()}")

    # Scale and PCA
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=50)
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.uns["pca"] = adata_hvg.uns["pca"]

    # Harmony batch correction across samples
    try:
        import harmonypy
        print("  Running Harmony batch correction...")
        ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, "sample",
                                   max_iter_harmony=30)
        Z = ho.Z_corr
        # harmonypy (PyTorch) may return a tensor; ensure numpy (n_cells × n_pcs)
        if hasattr(Z, "cpu"):
            Z = Z.cpu().numpy()
        Z = np.asarray(Z, dtype=np.float32)
        if Z.shape[0] == adata.n_obs:
            adata.obsm["X_pca_harmony"] = Z
        else:
            adata.obsm["X_pca_harmony"] = Z.T
        print(f"  Harmony shape: {adata.obsm['X_pca_harmony'].shape}")
        use_rep = "X_pca_harmony"
    except ImportError:
        print("  [INFO] harmonypy not installed, using uncorrected PCA")
        use_rep = "X_pca"

    # Neighbors & UMAP
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30, use_rep=use_rep)
    sc.tl.umap(adata, min_dist=0.3, spread=1.2)
    print(f"  UMAP computed ({adata.obsm['X_umap'].shape})")

    # Leiden at multiple resolutions
    for res in [0.3, 0.5, 0.8, 1.0]:
        key = f"leiden_{str(res).replace('.','_')}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        n = adata.obs[key].nunique()
        print(f"  Leiden res={res}: {n} clusters")

    return adata


# ═══════════════════════════════════════════════════════════════════════
#  CELL-TYPE ANNOTATION
# ═══════════════════════════════════════════════════════════════════════

def annotate_celltypes(adata):
    """Score cell-type markers and assign dominant type per spot."""
    print("\n[2] Cell-type annotation...")
    avail = set(adata.var_names)

    score_cols = []
    for ct_name, genes in CELL_MARKERS.items():
        present = [g for g in genes if g in avail]
        col = f"score_{ct_name.replace(' ', '_').replace('+', 'p')}"
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=col, use_raw=False)
            score_cols.append((ct_name, col))
            print(f"  {ct_name}: {len(present)}/{len(genes)} markers")
        else:
            print(f"  [SKIP] {ct_name}: only {len(present)} markers found")

    # Assign cell type = max score (with minimum threshold)
    score_matrix = np.column_stack([adata.obs[col].values for _, col in score_cols])
    ct_names = [ct for ct, _ in score_cols]

    max_idx = np.argmax(score_matrix, axis=1)
    max_val = np.max(score_matrix, axis=1)

    # Threshold: if max score < 0.05 → Unassigned
    labels = np.array([ct_names[i] if max_val[j] > 0.01 else "Unassigned"
                        for j, i in enumerate(max_idx)])
    adata.obs["celltype"] = pd.Categorical(labels,
                                           categories=list(CELLTYPE_PAL.keys()))

    # Print composition
    counts = adata.obs["celltype"].value_counts()
    print("\n  Cell-type composition:")
    for ct, n in counts.items():
        pct = 100 * n / adata.n_obs
        print(f"    {ct:20s}: {n:6d} ({pct:5.1f}%)")

    # Also score immune signatures for later
    for sig_name, genes in IMMUNE_SIGS.items():
        present = [g for g in genes if g in avail]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=sig_name, use_raw=False)

    return adata, score_cols


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 1: UMAP — Annotated Cell Types (all spots)
# ═══════════════════════════════════════════════════════════════════════

def plot_umap_celltypes(adata):
    """UMAP coloured by annotated cell type — high-contrast."""
    print("\n[3a] UMAP — Annotated cell types...")

    fig, axes = plt.subplots(1, 2, figsize=(20, 8))

    # Panel A: all cell types
    ax = axes[0]
    umap = adata.obsm["X_umap"]
    # Plot Unassigned first (background), then each type on top
    order = ["Unassigned", "Hepatocyte", "Endothelial", "Fibroblast", "Tumor",
             "B cell", "Macrophage", "DC", "NK cell", "CD4+ T cell", "Treg",
             "CD8+ T cell"]
    for ct in order:
        mask = adata.obs["celltype"] == ct
        if mask.sum() == 0:
            continue
        color = CELLTYPE_PAL.get(ct, "#999999")
        size = 2 if ct in ("Unassigned", "Hepatocyte") else 4
        ax.scatter(umap[mask, 0], umap[mask, 1], c=color, s=size,
                   label=f"{ct} ({mask.sum():,})", alpha=0.7,
                   edgecolors="none", rasterized=True)
    ax.set_xlabel("UMAP 1", fontsize=11)
    ax.set_ylabel("UMAP 2", fontsize=11)
    ax.set_title("Cell-Type Annotation", fontsize=13, fontweight="bold")
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9, markerscale=3,
              handletextpad=0.3, borderpad=0.4)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Panel B: coloured by response
    ax = axes[1]
    for resp, color in [("Non_Responder", RESPONSE_PAL["Non_Responder"]),
                         ("Responder", RESPONSE_PAL["Responder"])]:
        mask = adata.obs["response"] == resp
        ax.scatter(umap[mask, 0], umap[mask, 1], c=color, s=2,
                   label=f"{resp.replace('_',' ')} ({mask.sum():,})",
                   alpha=0.5, edgecolors="none", rasterized=True)
    ax.set_xlabel("UMAP 1", fontsize=11)
    ax.set_ylabel("UMAP 2", fontsize=11)
    ax.set_title("IO Response", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right", framealpha=0.9, markerscale=4)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.suptitle("HCC Visium ST — UMAP of 17,292 Spots (7 Sections, GSE238264)",
                 fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(FIG / "02b_ct_umap_celltype_response.png", dpi=250, bbox_inches="tight")
    plt.close()
    print("  Saved.")


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 2: UMAP — By Sample
# ═══════════════════════════════════════════════════════════════════════

def plot_umap_sample(adata):
    """UMAP coloured by sample."""
    print("\n[3b] UMAP — By sample...")

    fig, ax = plt.subplots(figsize=(10, 8))
    umap = adata.obsm["X_umap"]

    for sname, color in SAMPLE_PAL.items():
        mask = adata.obs["sample"] == sname
        if mask.sum() == 0:
            continue
        resp = SAMPLES[sname]
        marker = "o" if resp == "Responder" else "s"
        ax.scatter(umap[mask, 0], umap[mask, 1], c=color, s=3, marker=marker,
                   label=f"{sname} ({resp[:1]})", alpha=0.6,
                   edgecolors="none", rasterized=True)
    ax.set_xlabel("UMAP 1", fontsize=11)
    ax.set_ylabel("UMAP 2", fontsize=11)
    ax.set_title("HCC Visium — UMAP Coloured by Sample\n"
                 "(● = Responder, ■ = Non-Responder)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9, markerscale=3)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    plt.tight_layout()
    plt.savefig(FIG / "02b_ct_umap_sample.png", dpi=250, bbox_inches="tight")
    plt.close()
    print("  Saved.")


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 3: MARKER HEATMAP (dot plot)
# ═══════════════════════════════════════════════════════════════════════

def plot_marker_heatmap(adata):
    """Dot plot of cell-type markers grouped by annotation."""
    print("\n[3c] Marker heatmap / dot plot...")

    # Collect top markers per cell type (only those present in data)
    avail = set(adata.var_names)
    marker_list = []
    type_groups = {}
    for ct_name, genes in CELL_MARKERS.items():
        present = [g for g in genes if g in avail][:5]  # top 5 per type
        marker_list.extend(present)
        for g in present:
            type_groups[g] = ct_name
    marker_list = list(dict.fromkeys(marker_list))  # deduplicate preserving order

    # Filter to cell types with spots
    valid_cts = [ct for ct in CELL_MARKERS.keys()
                 if (adata.obs["celltype"] == ct).sum() > 0]

    # Use scanpy dot plot
    fig, ax = plt.subplots(figsize=(18, 7))
    plt.close(fig)  # close the dummy

    sc.settings.figdir = str(FIG)

    # Custom dot plot via manual computation for better control
    # Compute mean expression and fraction expressing per cell type
    rows = []
    adata_raw = adata.raw if adata.raw is not None else adata
    for ct in valid_cts:
        mask = adata.obs["celltype"] == ct
        if mask.sum() < 5:
            continue
        sub = adata_raw[mask]
        for gene in marker_list:
            if gene not in adata_raw.var_names:
                continue
            if issparse(sub[:, gene].X):
                vals = np.asarray(sub[:, gene].X.todense()).flatten()
            else:
                vals = sub[:, gene].X.flatten()
            rows.append({
                "celltype": ct,
                "gene": gene,
                "mean_expr": np.mean(vals),
                "frac_expr": np.mean(vals > 0),
            })
    df = pd.DataFrame(rows)

    if df.empty:
        print("  [WARN] No data for dot plot")
        return

    # Pivot
    mean_piv = df.pivot_table(index="celltype", columns="gene", values="mean_expr",
                               aggfunc="first").reindex(index=valid_cts, columns=marker_list)
    frac_piv = df.pivot_table(index="celltype", columns="gene", values="frac_expr",
                               aggfunc="first").reindex(index=valid_cts, columns=marker_list)
    mean_piv = mean_piv.fillna(0)
    frac_piv = frac_piv.fillna(0)

    # Z-score across cell types for colour
    zscore_piv = (mean_piv - mean_piv.mean(axis=0)) / (mean_piv.std(axis=0) + 1e-8)
    zscore_piv = zscore_piv.clip(-2.5, 2.5)

    fig, ax = plt.subplots(figsize=(max(14, 0.45 * len(marker_list)),
                                     max(5, 0.5 * len(valid_cts))))

    cmap = plt.cm.RdYlBu_r
    norm = plt.Normalize(-2.5, 2.5)

    for i, ct in enumerate(valid_cts):
        for j, gene in enumerate(marker_list):
            z = zscore_piv.loc[ct, gene] if ct in zscore_piv.index and gene in zscore_piv.columns else 0
            frac = frac_piv.loc[ct, gene] if ct in frac_piv.index and gene in frac_piv.columns else 0
            size = max(frac * 300, 3)
            ax.scatter(j, i, s=size, c=[cmap(norm(z))],
                       edgecolors="#333333", linewidth=0.3, zorder=3)

    ax.set_xticks(range(len(marker_list)))
    ax.set_xticklabels(marker_list, rotation=90, fontsize=7, fontstyle="italic")
    ax.set_yticks(range(len(valid_cts)))
    ax.set_yticklabels(valid_cts, fontsize=9)

    # Colour the y-tick labels
    for tick_l in ax.get_yticklabels():
        ct_name = tick_l.get_text()
        tick_l.set_color(CELLTYPE_PAL.get(ct_name, "black"))
        tick_l.set_fontweight("bold")

    # Gene grouping lines
    prev_group = None
    for j, gene in enumerate(marker_list):
        grp = type_groups.get(gene, "")
        if grp != prev_group and prev_group is not None:
            ax.axvline(j - 0.5, color="#BBBBBB", linewidth=0.5, linestyle="--")
        prev_group = grp

    ax.set_xlim(-0.5, len(marker_list) - 0.5)
    ax.set_ylim(-0.5, len(valid_cts) - 0.5)
    ax.invert_yaxis()
    ax.set_title("Cell-Type Marker Expression — Dot Plot\n"
                 "Size = fraction expressing · Colour = z-scored mean expression",
                 fontsize=12, fontweight="bold")

    # Colour bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label("Z-scored expression", fontsize=9)

    # Size legend
    for frac_val, label in [(0.2, "20%"), (0.5, "50%"), (0.8, "80%")]:
        ax.scatter([], [], s=frac_val * 300, c="gray", edgecolors="black",
                   linewidth=0.3, label=f"{label} expressing")
    ax.legend(fontsize=7, loc="lower right", framealpha=0.9, title="Fraction",
              title_fontsize=8)

    plt.tight_layout()
    plt.savefig(FIG / "02b_ct_marker_dotplot.png", dpi=250, bbox_inches="tight")
    plt.close()
    print("  Saved.")


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 3b: MARKER HEATMAP (classic clustered heatmap)
# ═══════════════════════════════════════════════════════════════════════

def plot_marker_clustered_heatmap(adata):
    """Clustered heatmap: mean marker expression per cell type."""
    print("\n[3c-bis] Clustered heatmap of markers...")

    avail = set(adata.var_names)
    marker_list = []
    for ct_name, genes in CELL_MARKERS.items():
        present = [g for g in genes if g in avail][:4]
        marker_list.extend(present)
    marker_list = list(dict.fromkeys(marker_list))

    valid_cts = [ct for ct in CELL_MARKERS.keys()
                 if (adata.obs["celltype"] == ct).sum() >= 5]

    adata_raw = adata.raw if adata.raw is not None else adata
    mat = np.zeros((len(valid_cts), len(marker_list)))
    for i, ct in enumerate(valid_cts):
        mask = adata.obs["celltype"] == ct
        sub = adata_raw[mask]
        for j, gene in enumerate(marker_list):
            if gene not in adata_raw.var_names:
                continue
            if issparse(sub[:, gene].X):
                vals = np.asarray(sub[:, gene].X.todense()).flatten()
            else:
                vals = sub[:, gene].X.flatten()
            mat[i, j] = np.mean(vals)

    # Z-score per gene
    mat_z = (mat - mat.mean(axis=0)) / (mat.std(axis=0) + 1e-8)
    mat_z = np.clip(mat_z, -2.5, 2.5)

    df_z = pd.DataFrame(mat_z, index=valid_cts, columns=marker_list)

    # Row colours
    row_colors = [CELLTYPE_PAL.get(ct, "#999999") for ct in valid_cts]

    g = sns.clustermap(
        df_z, cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5,
        row_cluster=True, col_cluster=False,
        row_colors=row_colors,
        linewidths=0.3, linecolor="white",
        figsize=(max(12, 0.35 * len(marker_list)), max(4, 0.5 * len(valid_cts))),
        dendrogram_ratio=(0.12, 0),
        cbar_kws={"label": "Z-score", "shrink": 0.6},
        xticklabels=True, yticklabels=True,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(),
                                  rotation=90, fontsize=7, fontstyle="italic")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=9)
    g.fig.suptitle("Cell-Type Marker Expression — Clustered Heatmap\n"
                   "Z-scored mean expression per annotated type",
                   fontsize=12, fontweight="bold", y=1.04)

    g.savefig(FIG / "02b_ct_marker_heatmap.png", dpi=250, bbox_inches="tight")
    plt.close()
    print("  Saved.")


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 4: UMAP — Immune Subtypes Only
# ═══════════════════════════════════════════════════════════════════════

def plot_umap_immune(adata):
    """UMAP of immune spots only, coloured by subtype and key signatures."""
    print("\n[3d] UMAP — Immune subtypes...")

    immune_mask = adata.obs["celltype"].isin(IMMUNE_TYPES)
    n_immune = immune_mask.sum()
    print(f"  {n_immune} immune spots")

    if n_immune < 100:
        print("  [SKIP] Too few immune spots for separate UMAP")
        return

    adata_imm = adata[immune_mask].copy()

    # Re-compute neighbours and UMAP for immune subset
    use_rep = "X_pca_harmony" if "X_pca_harmony" in adata_imm.obsm else "X_pca"
    sc.pp.neighbors(adata_imm, n_neighbors=15, n_pcs=20, use_rep=use_rep)
    sc.tl.umap(adata_imm, min_dist=0.25, spread=1.0)

    umap = adata_imm.obsm["X_umap"]

    # --- Multi-panel figure ---
    fig, axes = plt.subplots(2, 3, figsize=(21, 13))

    # Panel A: immune subtype
    ax = axes[0, 0]
    for ct in ["B cell", "DC", "Macrophage", "NK cell",
               "CD4+ T cell", "Treg", "CD8+ T cell"]:
        mask = adata_imm.obs["celltype"] == ct
        if mask.sum() == 0:
            continue
        ax.scatter(umap[mask, 0], umap[mask, 1],
                   c=CELLTYPE_PAL[ct], s=6, label=f"{ct} ({mask.sum():,})",
                   alpha=0.75, edgecolors="none", rasterized=True)
    ax.set_title("Immune Subtype", fontsize=11, fontweight="bold")
    ax.legend(fontsize=7, loc="best", framealpha=0.9, markerscale=2.5)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Panel B: response
    ax = axes[0, 1]
    for resp, color in [("Non_Responder", RESPONSE_PAL["Non_Responder"]),
                         ("Responder", RESPONSE_PAL["Responder"])]:
        mask = adata_imm.obs["response"] == resp
        ax.scatter(umap[mask, 0], umap[mask, 1], c=color, s=4,
                   label=resp.replace("_", " "), alpha=0.6,
                   edgecolors="none", rasterized=True)
    ax.set_title("IO Response", fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, loc="best", framealpha=0.9, markerscale=3)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Panels C-F: feature plots for key immune signatures
    sig_keys = ["T_cell_exhaustion", "T_cell_cytotoxicity",
                "TLS_signature", "Antigen_presentation"]
    sig_keys = [s for s in sig_keys if s in adata_imm.obs.columns]
    remaining_axes = [axes[0, 2], axes[1, 0], axes[1, 1], axes[1, 2]]

    for ax, sig in zip(remaining_axes[:len(sig_keys)], sig_keys):
        vals = adata_imm.obs[sig].values
        vmin, vmax = np.percentile(vals, [2, 98])
        sp = ax.scatter(umap[:, 0], umap[:, 1], c=vals,
                        cmap="magma", s=4, vmin=vmin, vmax=vmax,
                        alpha=0.75, edgecolors="none", rasterized=True)
        plt.colorbar(sp, ax=ax, shrink=0.6, pad=0.01)
        ax.set_title(sig.replace("_", " "), fontsize=10, fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    # Turn off unused axes
    for ax in remaining_axes[len(sig_keys):]:
        ax.axis("off")

    plt.suptitle(f"Immune Cell UMAP — {n_immune:,} Spots\n"
                 f"GSE238264 HCC IO Visium (7 sections)",
                 fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(FIG / "02b_ct_umap_immune_subtypes.png", dpi=250, bbox_inches="tight")
    plt.close()
    print("  Saved.")


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 5: SPATIAL MAPS — Cell-Type Annotation per Section
# ═══════════════════════════════════════════════════════════════════════

def plot_spatial_celltypes(adata):
    """Spatial maps coloured by annotated cell type, one per section."""
    print("\n[3e] Spatial cell-type maps...")

    samples = list(SAMPLES.keys())
    n_samples = len(samples)
    fig, axes = plt.subplots(2, 4, figsize=(24, 11))
    axes_flat = axes.flatten()

    for idx, sname in enumerate(samples):
        ax = axes_flat[idx]
        mask = adata.obs["sample"] == sname
        sub = adata[mask]
        coords = sub.obsm["spatial"]
        resp = sub.obs["response"].iloc[0]

        # Plot each cell type
        for ct in CELLTYPE_PAL.keys():
            ct_mask = sub.obs["celltype"] == ct
            if ct_mask.sum() == 0:
                continue
            color = CELLTYPE_PAL[ct]
            size = 3 if ct in ("Unassigned", "Hepatocyte") else 6
            ax.scatter(coords[ct_mask.values, 0], coords[ct_mask.values, 1],
                       c=color, s=size, edgecolors="none",
                       rasterized=True, alpha=0.8)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        resp_short = "R" if resp == "Responder" else "NR"
        ax.set_title(f"{sname} ({resp_short})\n{mask.sum()} spots",
                     fontsize=10, fontweight="bold",
                     color=RESPONSE_PAL[resp])
        for spine in ax.spines.values():
            spine.set_visible(False)

    # Last subplot: legend
    ax = axes_flat[-1]
    ax.axis("off")
    for ct, color in CELLTYPE_PAL.items():
        ax.scatter([], [], c=color, s=40, label=ct, edgecolors="none")
    ax.legend(fontsize=9, loc="center", framealpha=0.95, markerscale=2,
              title="Cell Type", title_fontsize=10)

    plt.suptitle("Spatial Cell-Type Annotation — GSE238264 HCC Visium",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG / "02b_ct_spatial_celltype_maps.png", dpi=200, bbox_inches="tight")
    plt.close()
    print("  Saved.")


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 6: UMAP — Leiden Clusters with Annotation Overlay
# ═══════════════════════════════════════════════════════════════════════

def plot_umap_leiden(adata):
    """UMAP coloured by Leiden clusters."""
    print("\n[3f] UMAP — Leiden clusters...")

    key = "leiden_0_8"
    n_clust = adata.obs[key].nunique()
    cmap = plt.cm.get_cmap("tab20", n_clust)
    palette = {str(i): cmap(i) for i in range(n_clust)}

    fig, axes = plt.subplots(1, 2, figsize=(20, 8))

    # Panel A: Leiden
    ax = axes[0]
    umap = adata.obsm["X_umap"]
    for cl in sorted(adata.obs[key].unique(), key=int):
        mask = adata.obs[key] == cl
        ax.scatter(umap[mask, 0], umap[mask, 1], c=[palette[cl]], s=2,
                   label=f"C{cl}", alpha=0.6, edgecolors="none", rasterized=True)
    ax.set_title(f"Leiden Clustering (res=0.8, {n_clust} clusters)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=6, loc="upper right", ncol=2, framealpha=0.9, markerscale=3)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Panel B: cell type (same as before but side by side)
    ax = axes[1]
    order = ["Unassigned", "Hepatocyte", "Endothelial", "Fibroblast", "Tumor",
             "B cell", "Macrophage", "DC", "NK cell", "CD4+ T cell", "Treg",
             "CD8+ T cell"]
    for ct in order:
        mask = adata.obs["celltype"] == ct
        if mask.sum() == 0:
            continue
        ax.scatter(umap[mask, 0], umap[mask, 1], c=CELLTYPE_PAL[ct], s=2,
                   label=ct, alpha=0.6, edgecolors="none", rasterized=True)
    ax.set_title("Cell-Type Annotation", fontsize=12, fontweight="bold")
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9, markerscale=3)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.suptitle("HCC Visium ST — Leiden Clusters vs Cell-Type Labels",
                 fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(FIG / "02b_ct_umap_leiden.png", dpi=250, bbox_inches="tight")
    plt.close()
    print("  Saved.")


# ═══════════════════════════════════════════════════════════════════════
#  SAVE ANNOTATION TABLE
# ═══════════════════════════════════════════════════════════════════════

def save_annotations(adata):
    """Save per-spot annotation table."""
    print("\n[4] Saving annotation table...")
    cols = ["sample", "response", "celltype", "leiden_0_8"]
    cols += [c for c in adata.obs.columns if c.startswith("score_")]
    df = adata.obs[cols].copy()
    df.to_csv(TAB / "gse238264_celltype_annotations.csv")
    print(f"  Saved {df.shape[0]} spots × {df.shape[1]} columns")

    # Per-sample composition
    comp = adata.obs.groupby(["sample", "response", "celltype"]).size().reset_index(name="n_spots")
    comp.to_csv(TAB / "gse238264_celltype_composition.csv", index=False)
    print(f"  Composition table: {comp.shape[0]} rows")


# ═══════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("02b — Visium Cell-Type Annotation & UMAP")
    print("=" * 70)

    adata = load_and_merge()
    adata = preprocess_and_embed(adata)
    adata, score_cols = annotate_celltypes(adata)

    # Figures
    plot_umap_celltypes(adata)
    plot_umap_sample(adata)
    plot_marker_heatmap(adata)
    plot_marker_clustered_heatmap(adata)
    plot_umap_immune(adata)
    plot_spatial_celltypes(adata)
    plot_umap_leiden(adata)
    save_annotations(adata)

    print("\n" + "=" * 70)
    print("[DONE] All cell-type annotation & UMAP figures generated")
    print(f"  Figures: {FIG}")
    print(f"  Tables:  {TAB}")
    print("=" * 70)


if __name__ == "__main__":
    main()
