#!/usr/bin/env python3
"""
02b_improve_spatial_plots.py — High-contrast publication-quality remake
of Visium spatial & UMAP figures.

Regenerates the following figures with improved colour contrast,
dimensions, colour-maps, and labelling:

  1. UMAP cell-type + response       (02b_ct_umap_celltype_response.png)
  2. UMAP by sample                  (02b_ct_umap_sample.png)
  3. UMAP Leiden vs cell-type        (02b_ct_umap_leiden.png)
  4. UMAP immune subtypes            (02b_ct_umap_immune_subtypes.png)
  5. Spatial cell-type maps           (02b_ct_spatial_celltype_maps.png)
  6. Spatial IO signature maps        (02b_visium_spatial_signature_maps.png)
  7. Cell-type marker spatial         (02b_visium_celltype_spatial_maps.png)
  8. Neighbourhood enrichment         (02b_visium_nhood_enrichment.png)
  9. Marker dot plot                  (02b_ct_marker_dotplot.png)
 10. Marker heatmap                   (02b_ct_marker_heatmap.png)
"""

import warnings, gc
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
import seaborn as sns
import sys; sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
from config.publication_style import apply_style, PAL, CMAP, label_clean
apply_style()
from pathlib import Path
from collections import OrderedDict
from scipy.sparse import issparse

sc.settings.verbosity = 1

# ── Paths ──────────────────────────────────────────────────────────────
PROJECT = Path(__file__).resolve().parents[2]
DATA238 = PROJECT / "data" / "spatial_data" / "GSE238264_HCC_IO_Visium"
OUT     = PROJECT / "out" / "02b_spatial_biomarkers"
FIG     = OUT / "figures"
TAB     = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

# ── Global settings ────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.linewidth": 0.8,
    "axes.labelsize": 11,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
    "savefig.dpi": 300,
})

# ── Palettes ───────────────────────────────────────────────────────────
# Cell-type: muted hepatocyte/unassigned so immune pops
CELLTYPE_PAL = OrderedDict({
    "CD8+ T cell":   "#D62728",   # bold red
    "CD4+ T cell":   "#FF7F0E",   # vivid orange
    "Treg":          "#8C564B",   # brown
    "NK cell":       "#E377C2",   # pink
    "B cell":        "#7B2D8E",   # deep violet
    "Macrophage":    "#1F77B4",   # strong blue
    "DC":            "#2CA02C",   # green
    "Hepatocyte":    "#B8B8B8",   # ← visible muted grey (not too pale)
    "Endothelial":   "#BCBD22",   # olive
    "Fibroblast":    "#17BECF",   # teal
    "Tumor":         "#505050",   # dark grey
    "Unassigned":    "#D5D5D5",   # light grey
})

IMMUNE_TYPES = {"CD8+ T cell", "CD4+ T cell", "Treg", "NK cell",
                "B cell", "Macrophage", "DC"}

RESPONSE_PAL = {"Responder": "#0077BB", "Non_Responder": "#CC3311"}

SAMPLE_PAL = OrderedDict({
    "HCC1R":  "#E31A1C", "HCC2R":  "#FF7F00",
    "HCC3R":  "#33A02C", "HCC4R":  "#1F78B4",
    "HCC5NR": "#6A3D9A", "HCC6NR": "#B15928", "HCC7NR": "#E7298A",
})

SAMPLES = OrderedDict({
    "HCC1R": "Responder", "HCC2R": "Responder",
    "HCC3R": "Responder", "HCC4R": "Responder",
    "HCC5NR": "Non_Responder", "HCC6NR": "Non_Responder",
    "HCC7NR": "Non_Responder",
})

# ── Markers & signatures ──────────────────────────────────────────────
CELL_MARKERS = OrderedDict({
    "CD8+ T cell":  ["CD8A","CD8B","GZMB","PRF1","GZMA","GZMK","NKG7"],
    "CD4+ T cell":  ["CD4","IL7R","CCR7","LEF1","TCF7"],
    "Treg":         ["FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"],
    "NK cell":      ["NKG7","GNLY","KLRD1","KLRF1","NCAM1"],
    "B cell":       ["CD79A","MS4A1","CD19","CD79B","IGHG1","IGKC"],
    "Macrophage":   ["CD68","CD163","C1QA","CSF1R","MSR1","MARCO"],
    "DC":           ["CLEC9A","XCR1","BATF3","IRF8","LAMP3","FSCN1"],
    "Hepatocyte":   ["ALB","APOB","HP","TF","APOA1","SERPINA1"],
    "Endothelial":  ["PECAM1","CDH5","VWF","ENG","PLVAP"],
    "Fibroblast":   ["COL1A1","COL3A1","DCN","LUM","ACTA2"],
    "Tumor":        ["GPC3","AFP","EPCAM","KRT19","SPINK1"],
})

IO_SIGS = OrderedDict({
    "T_cell_cytotoxicity":  ["GZMA","GZMB","GZMH","GZMK","PRF1","GNLY","NKG7"],
    "T_cell_exhaustion":    ["PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","TOX","ENTPD1"],
    "TLS_signature":        ["CXCL13","CXCR5","CR2","CD79A","MS4A1","CD19","IGHG1"],
    "TGF_beta":             ["TGFB1","TGFB2","TGFB3","ACTA2","COL1A1","SERPINE1"],
    "IFN_gamma":            ["IFNG","STAT1","IRF1","CXCL9","CXCL10","CXCL11","IDO1"],
    "Antigen_presentation": ["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","PSMB9"],
    "M2_macrophage":        ["CD163","MRC1","MSR1","CD200R1","IL10"],
    "Immune_checkpoint":    ["CD274","PDCD1LG2","PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","VSIR"],
    "Angiogenesis":         ["VEGFA","PECAM1","KDR","FLT1","ANGPT2","FGF2"],
    "WNT_beta_catenin":     ["CTNNB1","APC","AXIN2","LGR5","MYC"],
})

IMMUNE_SIGS = OrderedDict({
    "T_cell_exhaustion":    ["PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","TOX","ENTPD1"],
    "T_cell_cytotoxicity":  ["GZMA","GZMB","GZMH","GZMK","PRF1","GNLY","NKG7"],
    "TLS_signature":        ["CXCL13","CXCR5","CR2","CD79A","MS4A1","CD19","IGHG1"],
    "Antigen_presentation": ["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","PSMB9"],
})

# Cell-type markers for old-style spatial maps
CT_MARKERS_OLD = {
    "CD8_T":       ["CD8A","CD8B","GZMB","PRF1"],
    "Macrophage":  ["CD68","CD163","C1QA","CSF1R"],
    "Tumor":       ["GPC3","AFP","EPCAM"],
    "B_cell":      ["CD79A","MS4A1","CD19"],
}


# ═══════════════════════════════════════════════════════════════════════
#  DATA LOADING
# ═══════════════════════════════════════════════════════════════════════

def load_section(name: str) -> sc.AnnData:
    """Load one Visium section, normalise, score."""
    path = DATA238 / f"{name}.h5ad"
    ad = sc.read_h5ad(str(path))
    ad.var_names_make_unique()
    ad.obs_names = [f"{name}_{b}" for b in ad.obs_names]
    ad.obsm["spatial"] = ad.obs[["pxl_col", "pxl_row"]].values.astype(float)
    ad.obs["sample"] = name
    ad.obs["response"] = SAMPLES[name]
    ad.layers["counts"] = ad.X.copy()
    return ad


def load_all_sections():
    """Load, normalise, merge all 7 sections."""
    adatas = {}
    for name in SAMPLES:
        ad = load_section(name)
        adatas[name] = ad
        print(f"  {name}: {ad.n_obs} spots")
    return adatas


def merge_and_embed(adatas: dict) -> sc.AnnData:
    """Merge → normalise → HVG → PCA → Harmony → UMAP → Leiden."""
    ads = list(adatas.values())
    adata = ads[0].concatenate(ads[1:], batch_key="batch", join="inner")
    adata.obs["sample"] = adata.obs["sample"].astype(str)
    adata.obs["response"] = adata.obs["response"].astype(str)
    print(f"  Merged: {adata.n_obs} spots × {adata.n_vars} genes")

    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="sample")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=50)
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

    # Harmony
    try:
        import harmonypy
        ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, "sample",
                                   max_iter_harmony=30)
        Z = ho.Z_corr
        if hasattr(Z, "cpu"):
            Z = Z.cpu().numpy()
        Z = np.asarray(Z, dtype=np.float32)
        adata.obsm["X_pca_harmony"] = Z if Z.shape[0] == adata.n_obs else Z.T
        use_rep = "X_pca_harmony"
        print(f"  Harmony done: {adata.obsm['X_pca_harmony'].shape}")
    except ImportError:
        use_rep = "X_pca"

    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30, use_rep=use_rep)
    sc.tl.umap(adata, min_dist=0.3, spread=1.2)

    for res in [0.5, 0.8]:
        key = f"leiden_{str(res).replace('.','_')}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        print(f"  Leiden res={res}: {adata.obs[key].nunique()} clusters")

    return adata


def annotate(adata):
    """Score markers, assign cell types, score IO signatures.

    Uses z-score normalisation across cell-type scores so that
    ubiquitous liver genes (hepatocyte markers) do not dominate every
    spot.  Each score is centred and scaled before the max-score
    assignment, giving non-hepatocyte types a fair chance.
    """
    avail = set(adata.var_names)

    # Cell-type scoring
    score_cols = []
    for ct_name, genes in CELL_MARKERS.items():
        present = [g for g in genes if g in avail]
        col = f"score_{ct_name.replace(' ','_').replace('+','p')}"
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=col, use_raw=False)
            score_cols.append((ct_name, col))

    # Build raw score matrix and z-score per cell type (column)
    raw_scores = np.column_stack([adata.obs[col].values for _, col in score_cols])
    col_mean = raw_scores.mean(axis=0, keepdims=True)
    col_std  = raw_scores.std(axis=0, keepdims=True)
    col_std[col_std < 1e-8] = 1.0
    z_scores = (raw_scores - col_mean) / col_std

    ct_names = [ct for ct, _ in score_cols]
    max_idx = np.argmax(z_scores, axis=1)
    max_val = np.max(z_scores, axis=1)
    # Require at least z > 0.3 to call a type; otherwise Unassigned
    labels = np.array([ct_names[i] if max_val[j] > 0.3 else "Unassigned"
                       for j, i in enumerate(max_idx)])
    adata.obs["celltype"] = pd.Categorical(labels, categories=list(CELLTYPE_PAL.keys()))

    # Immune signatures
    for sig_name, genes in IMMUNE_SIGS.items():
        present = [g for g in genes if g in avail]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=sig_name, use_raw=False)

    # IO signatures
    for sig_name, genes in IO_SIGS.items():
        present = [g for g in genes if g in avail]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=sig_name, use_raw=False)

    # Old-style cell-type markers for spatial
    for ct_name, genes in CT_MARKERS_OLD.items():
        present = [g for g in genes if g in avail]
        col = f"ct_{ct_name}"
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=col, use_raw=False)

    # Immunotype classification
    immune_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    immune_sum = adata.obs[immune_cols].sum(axis=1) if immune_cols else pd.Series(0, index=adata.obs_names)
    q33 = np.percentile(immune_sum, 33)
    q66 = np.percentile(immune_sum, 66)
    adata.obs["spatial_immunotype"] = pd.Categorical(
        np.where(immune_sum > q66, "Immune_Hot",
                 np.where(immune_sum < q33, "Immune_Desert", "Immune_Excluded")),
        categories=["Immune_Hot", "Immune_Excluded", "Immune_Desert", "Immune_Mixed"]
    )

    return adata, score_cols


# ═══════════════════════════════════════════════════════════════════════
#  HELPER: clean UMAP axes
# ═══════════════════════════════════════════════════════════════════════

def clean_umap_ax(ax, xlabel="UMAP 1", ylabel="UMAP 2"):
    ax.set_xlabel(xlabel, fontsize=10, labelpad=2)
    ax.set_ylabel(ylabel, fontsize=10, labelpad=2)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


# ═══════════════════════════════════════════════════════════════════════
#  FIG 1: UMAP — Cell Type + Response  (improved)
# ═══════════════════════════════════════════════════════════════════════

def fig_umap_celltype_response(adata):
    print("\n[Fig 1] UMAP — Cell types + IO response...")
    umap = adata.obsm["X_umap"]

    fig, axes = plt.subplots(1, 2, figsize=(20, 10))

    # ---- Panel A: Cell types (layered: background → rare on top) ----
    ax = axes[0]
    ax.set_facecolor("#F5F5F5")

    # Background types first (small, muted)
    bg_types = ["Unassigned", "Hepatocyte"]
    mid_types = ["Tumor", "Fibroblast", "Endothelial"]
    fg_types = ["B cell", "Macrophage", "DC", "NK cell",
                "CD4+ T cell", "Treg", "CD8+ T cell"]

    for ct in bg_types:
        m = adata.obs["celltype"] == ct
        if m.sum() == 0: continue
        ax.scatter(umap[m, 0], umap[m, 1], c=CELLTYPE_PAL[ct], s=1.5,
                   alpha=0.2, edgecolors="none", rasterized=True, label=f"{ct} ({m.sum():,})")
    for ct in mid_types:
        m = adata.obs["celltype"] == ct
        if m.sum() == 0: continue
        ax.scatter(umap[m, 0], umap[m, 1], c=CELLTYPE_PAL[ct], s=10,
                   alpha=0.8, edgecolors="white", linewidth=0.2, rasterized=True,
                   label=f"{ct} ({m.sum():,})")
    for ct in fg_types:
        m = adata.obs["celltype"] == ct
        if m.sum() == 0: continue
        ax.scatter(umap[m, 0], umap[m, 1], c=CELLTYPE_PAL[ct], s=22,
                   alpha=0.95, edgecolors="white", linewidth=0.4, rasterized=True,
                   zorder=5, label=f"{ct} ({m.sum():,})")

    ax.set_title("Cell-Type Annotation", fontsize=14, fontweight="bold")
    leg = ax.legend(fontsize=8, loc="upper right", framealpha=0.95,
                    markerscale=2.5, handletextpad=0.3, borderpad=0.5,
                    edgecolor="#CCCCCC")
    clean_umap_ax(ax)

    # ---- Panel B: Response ----
    ax = axes[1]
    ax.set_facecolor("#F5F5F5")
    for resp, color in [("Non_Responder", RESPONSE_PAL["Non_Responder"]),
                        ("Responder", RESPONSE_PAL["Responder"])]:
        m = adata.obs["response"] == resp
        ax.scatter(umap[m, 0], umap[m, 1], c=color, s=4,
                   alpha=0.6, edgecolors="none", rasterized=True,
                   label=f"{resp.replace('_',' ')} ({m.sum():,})")
    ax.set_title("IO Response", fontsize=14, fontweight="bold")
    ax.legend(fontsize=11, loc="upper right", framealpha=0.95, markerscale=4,
              edgecolor="#CCCCCC")
    clean_umap_ax(ax)

    fig.suptitle("HCC Visium ST — UMAP of 17,292 Spots (7 Sections, GSE238264)",
                 fontsize=15, fontweight="bold", y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG / "02b_ct_umap_celltype_response.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved.")


# ═══════════════════════════════════════════════════════════════════════
#  FIG 2: UMAP — By Sample
# ═══════════════════════════════════════════════════════════════════════

def fig_umap_sample(adata):
    print("\n[Fig 2] UMAP — By sample...")
    umap = adata.obsm["X_umap"]

    fig, ax = plt.subplots(figsize=(11, 9))
    ax.set_facecolor("#F8F8F8")

    for sname, color in SAMPLE_PAL.items():
        m = adata.obs["sample"] == sname
        if m.sum() == 0: continue
        resp = SAMPLES[sname]
        marker = "o" if resp == "Responder" else "D"
        ax.scatter(umap[m, 0], umap[m, 1], c=color, s=8, marker=marker,
                   alpha=0.65, edgecolors="none", rasterized=True,
                   label=f"{sname} ({'R' if resp=='Responder' else 'NR'})")

    ax.set_title("HCC Visium — UMAP by Sample\n(● Responder · ◆ Non-Responder)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right", framealpha=0.95, markerscale=3,
              edgecolor="#CCCCCC")
    clean_umap_ax(ax)
    plt.tight_layout()
    fig.savefig(FIG / "02b_ct_umap_sample.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved.")


# ═══════════════════════════════════════════════════════════════════════
#  FIG 3: UMAP — Leiden vs Cell Type
# ═══════════════════════════════════════════════════════════════════════

def fig_umap_leiden(adata):
    print("\n[Fig 3] UMAP — Leiden vs cell-type...")
    key = "leiden_0_8"
    n_clust = adata.obs[key].nunique()

    # High-contrast Leiden palette (tab20 + Set1)
    base_colors = (list(plt.cm.tab10.colors) + list(plt.cm.Set3.colors))[:n_clust]
    leiden_pal = {str(i): base_colors[i] for i in range(n_clust)}

    umap = adata.obsm["X_umap"]
    fig, axes = plt.subplots(1, 2, figsize=(20, 10))

    # ---- Panel A: Leiden ----
    ax = axes[0]
    ax.set_facecolor("#F5F5F5")
    for cl in sorted(adata.obs[key].unique(), key=int):
        m = adata.obs[key] == cl
        ax.scatter(umap[m, 0], umap[m, 1], c=[leiden_pal[cl]], s=4,
                   alpha=0.7, edgecolors="none", rasterized=True, label=f"C{cl}")
    ax.set_title(f"Leiden Clustering (res=0.8, {n_clust} clusters)",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=8, loc="upper right", ncol=2, framealpha=0.95,
              markerscale=4, edgecolor="#CCCCCC")
    clean_umap_ax(ax)

    # ---- Panel B: Cell type (layered) ----
    ax = axes[1]
    ax.set_facecolor("#F5F5F5")
    bg = ["Unassigned", "Hepatocyte"]
    mid = ["Tumor", "Fibroblast", "Endothelial"]
    fg = ["B cell", "Macrophage", "DC", "NK cell", "CD4+ T cell", "Treg", "CD8+ T cell"]
    for ct in bg:
        m = adata.obs["celltype"] == ct
        if m.sum() == 0: continue
        ax.scatter(umap[m, 0], umap[m, 1], c=CELLTYPE_PAL[ct], s=1.5,
                   alpha=0.2, edgecolors="none", rasterized=True, label=ct)
    for ct in mid:
        m = adata.obs["celltype"] == ct
        if m.sum() == 0: continue
        ax.scatter(umap[m, 0], umap[m, 1], c=CELLTYPE_PAL[ct], s=10,
                   alpha=0.8, edgecolors="white", linewidth=0.2,
                   rasterized=True, label=ct)
    for ct in fg:
        m = adata.obs["celltype"] == ct
        if m.sum() == 0: continue
        ax.scatter(umap[m, 0], umap[m, 1], c=CELLTYPE_PAL[ct], s=24,
                   alpha=0.95, edgecolors="white", linewidth=0.4,
                   rasterized=True, zorder=5, label=ct)
    ax.set_title("Cell-Type Annotation", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.95, markerscale=3,
              edgecolor="#CCCCCC")
    clean_umap_ax(ax)

    fig.suptitle("HCC Visium ST — Leiden Clusters vs Cell-Type Labels",
                 fontsize=14, fontweight="bold", y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIG / "02b_ct_umap_leiden.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved.")


# ═══════════════════════════════════════════════════════════════════════
#  FIG 4: UMAP — Immune Subtypes (improved)
# ═══════════════════════════════════════════════════════════════════════

def fig_umap_immune(adata):
    print("\n[Fig 4] UMAP — Immune subtypes...")
    mask = adata.obs["celltype"].isin(IMMUNE_TYPES)
    n_imm = mask.sum()
    print(f"  {n_imm} immune spots")
    if n_imm < 50:
        print("  [SKIP]")
        return

    adata_imm = adata[mask].copy()
    use_rep = "X_pca_harmony" if "X_pca_harmony" in adata_imm.obsm else "X_pca"
    sc.pp.neighbors(adata_imm, n_neighbors=15, n_pcs=20, use_rep=use_rep)
    sc.tl.umap(adata_imm, min_dist=0.25, spread=1.0)
    umap = adata_imm.obsm["X_umap"]

    # Adaptive point size: larger when fewer spots
    base_s = max(20, min(60, 12000 / max(n_imm, 1)))
    feat_s = max(14, min(40, 8000 / max(n_imm, 1)))

    fig, axes = plt.subplots(2, 3, figsize=(20, 13))
    fig.patch.set_facecolor("white")

    # ---- Panel A: Immune subtype ----
    ax = axes[0, 0]
    ax.set_facecolor("#F0F0F0")
    draw_order = ["B cell", "DC", "Macrophage", "NK cell",
                  "CD4+ T cell", "Treg", "CD8+ T cell"]
    for ct in draw_order:
        m = adata_imm.obs["celltype"] == ct
        if m.sum() == 0: continue
        ax.scatter(umap[m, 0], umap[m, 1], c=CELLTYPE_PAL[ct], s=base_s,
                   alpha=0.9, edgecolors="white", linewidth=0.3,
                   rasterized=True, zorder=3, label=f"{ct} ({m.sum():,})")
    ax.set_title("Immune Subtype", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="best", framealpha=0.95, markerscale=2,
              edgecolor="#CCCCCC", handletextpad=0.4)
    clean_umap_ax(ax)

    # ---- Panel B: Response ----
    ax = axes[0, 1]
    ax.set_facecolor("#F0F0F0")
    for resp, color in [("Non_Responder", RESPONSE_PAL["Non_Responder"]),
                        ("Responder", RESPONSE_PAL["Responder"])]:
        m = adata_imm.obs["response"] == resp
        ax.scatter(umap[m, 0], umap[m, 1], c=color, s=base_s * 0.7,
                   alpha=0.8, edgecolors="white", linewidth=0.2,
                   rasterized=True, label=resp.replace("_", " "))
    ax.set_title("IO Response", fontsize=13, fontweight="bold")
    ax.legend(fontsize=11, loc="best", framealpha=0.95, markerscale=2.5,
              edgecolor="#CCCCCC")
    clean_umap_ax(ax)

    # ---- Panels C-F: feature plots (consistent white background) ----
    sig_keys = [s for s in ["T_cell_exhaustion", "T_cell_cytotoxicity",
                            "TLS_signature", "Antigen_presentation"]
                if s in adata_imm.obs.columns]
    remaining_axes = [axes[0, 2], axes[1, 0], axes[1, 1], axes[1, 2]]
    feat_cmaps = ["YlOrRd", "YlGnBu", "PuBuGn", "YlOrBr"]

    for idx, (ax, sig) in enumerate(zip(remaining_axes[:len(sig_keys)], sig_keys)):
        ax.set_facecolor("#F0F0F0")
        vals = adata_imm.obs[sig].values
        vmin, vmax = np.percentile(vals, [2, 98])
        if vmax - vmin < 0.01:
            vmin, vmax = vals.min(), vals.max()
        cmap = feat_cmaps[idx % len(feat_cmaps)]
        sp = ax.scatter(umap[:, 0], umap[:, 1], c=vals,
                        cmap=cmap, s=feat_s, vmin=vmin, vmax=vmax,
                        alpha=0.9, edgecolors="none", rasterized=True)
        cbar = plt.colorbar(sp, ax=ax, shrink=0.6, pad=0.02, aspect=20)
        cbar.ax.tick_params(labelsize=9)
        title = sig.replace("_", " ")
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([])
        for sp_ax in ax.spines.values():
            sp_ax.set_visible(False)

    for ax in remaining_axes[len(sig_keys):]:
        ax.axis("off")

    fig.suptitle(f"Immune Cell UMAP — {n_imm:,} Spots · GSE238264 (7 sections)",
                 fontsize=15, fontweight="bold", y=0.99)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(FIG / "02b_ct_umap_immune_subtypes.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved.")


# ═══════════════════════════════════════════════════════════════════════
#  FIG 5: Spatial Cell-Type Maps  (improved)
# ═══════════════════════════════════════════════════════════════════════

def fig_spatial_celltypes(adata):
    print("\n[Fig 5] Spatial cell-type maps...")
    samples = list(SAMPLES.keys())
    fig, axes = plt.subplots(2, 4, figsize=(26, 12))
    axes_flat = axes.flatten()

    for idx, sname in enumerate(samples):
        ax = axes_flat[idx]
        ax.set_facecolor("white")
        m = adata.obs["sample"] == sname
        sub = adata[m]
        coords = sub.obsm["spatial"]
        resp = sub.obs["response"].iloc[0]

        # Adaptive spot size: scale inversely with coordinate range
        coord_span = max(coords[:, 0].max() - coords[:, 0].min(),
                         coords[:, 1].max() - coords[:, 1].min(), 1)
        base_s = max(6, min(30, 18000 / coord_span))

        # Draw order: background → mid → immune foreground
        bg_ct = ["Unassigned", "Hepatocyte"]
        mid_ct = ["Tumor", "Fibroblast", "Endothelial"]
        fg_ct = [ct for ct in CELLTYPE_PAL if ct not in bg_ct and ct not in mid_ct]

        for ct in bg_ct:
            ct_m = sub.obs["celltype"] == ct
            if ct_m.sum() == 0: continue
            ax.scatter(coords[ct_m.values, 0], coords[ct_m.values, 1],
                       c=CELLTYPE_PAL[ct], s=base_s * 0.9, alpha=0.7,
                       edgecolors="none", rasterized=True)
        for ct in mid_ct:
            ct_m = sub.obs["celltype"] == ct
            if ct_m.sum() == 0: continue
            ax.scatter(coords[ct_m.values, 0], coords[ct_m.values, 1],
                       c=CELLTYPE_PAL[ct], s=base_s * 1.6, alpha=0.85,
                       edgecolors="white", linewidth=0.3, rasterized=True)
        for ct in fg_ct:
            ct_m = sub.obs["celltype"] == ct
            if ct_m.sum() == 0: continue
            ax.scatter(coords[ct_m.values, 0], coords[ct_m.values, 1],
                       c=CELLTYPE_PAL[ct], s=base_s * 2.4, alpha=0.95,
                       edgecolors="white", linewidth=0.4, rasterized=True, zorder=5)

        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        resp_short = "R" if resp == "Responder" else "NR"
        title_col = RESPONSE_PAL[resp]
        ax.set_title(f"{sname} ({resp_short}) · {m.sum()} spots",
                     fontsize=10, fontweight="bold", color=title_col)
        for sp in ax.spines.values():
            sp.set_linewidth(0.5); sp.set_color("#CCCCCC")

    # Legend in last slot
    ax = axes_flat[-1]
    ax.axis("off")
    for ct, color in CELLTYPE_PAL.items():
        ax.scatter([], [], c=color, s=60, label=ct, edgecolors="white", linewidth=0.3)
    ax.legend(fontsize=9.5, loc="center", framealpha=0.95, markerscale=1.5,
              title="Cell Type", title_fontsize=11, edgecolor="#CCCCCC",
              handletextpad=0.4, borderpad=0.6)

    fig.suptitle("Spatial Cell-Type Annotation — GSE238264 HCC Visium",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    fig.savefig(FIG / "02b_ct_spatial_celltype_maps.png", dpi=250, bbox_inches="tight")
    plt.close()
    print("  Saved.")


# ═══════════════════════════════════════════════════════════════════════
#  FIG 6: Spatial IO Signature Maps  (SPLIT: R page + NR page)
# ═══════════════════════════════════════════════════════════════════════

def fig_io_signature_maps(adata, adatas):
    """Regenerate spatial IO maps — split into Responder and Non-Responder
    with quantile-based colour scaling and inferno/plasma colormaps."""
    print("\n[Fig 6] Spatial IO signature maps...")

    key_sigs = ["T_cell_cytotoxicity", "T_cell_exhaustion", "TLS_signature",
                "TGF_beta", "IFN_gamma", "Antigen_presentation"]
    key_sigs = [s for s in key_sigs if s in adata.obs.columns]
    n_sigs = len(key_sigs)

    # Score each individual section
    avail = set(adata.var_names)
    for name, ad in adatas.items():
        sc.pp.filter_genes(ad, min_cells=1)
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
        section_avail = set(ad.var_names)
        for sig_name, genes in IO_SIGS.items():
            present = [g for g in genes if g in section_avail]
            if len(present) >= 2 and sig_name not in ad.obs.columns:
                sc.tl.score_genes(ad, present, score_name=sig_name, use_raw=False)

    # Use sequential high-contrast cmaps (perceptually uniform, light bg)
    sig_cmaps = ["YlOrRd", "BuPu", "YlGnBu", "OrRd", "PuBuGn", "YlOrBr"]

    for tag, sample_names in [("Responder", ["HCC1R","HCC2R","HCC3R","HCC4R"]),
                               ("Non_Responder", ["HCC5NR","HCC6NR","HCC7NR"])]:
        n_samples = len(sample_names)
        fig, axes = plt.subplots(n_samples, n_sigs,
                                 figsize=(4.5 * n_sigs, 4.5 * n_samples))
        if n_samples == 1:
            axes = axes[np.newaxis, :]

        for i, sname in enumerate(sample_names):
            ad = adatas[sname]
            coords = ad.obsm["spatial"]
            for j, sig in enumerate(key_sigs):
                ax = axes[i, j]
                ax.set_facecolor("white")

                if sig in ad.obs.columns:
                    vals = ad.obs[sig].values
                    # Use tighter percentile clipping for better dynamic range
                    vmin, vmax = np.percentile(vals, [5, 95])
                    if vmax - vmin < 0.05:
                        vmin, vmax = vals.min(), vals.max()
                    cmap = sig_cmaps[j % len(sig_cmaps)]
                    sc_p = ax.scatter(coords[:, 0], coords[:, 1],
                                      c=vals, cmap=cmap, s=18,
                                      vmin=vmin, vmax=vmax,
                                      edgecolors="none", rasterized=True, alpha=0.9)
                    cbar = plt.colorbar(sc_p, ax=ax, shrink=0.7, pad=0.02, aspect=18)
                    cbar.ax.tick_params(labelsize=8)

                ax.set_aspect("equal"); ax.invert_yaxis()
                ax.set_xticks([]); ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_visible(False)
                if i == 0:
                    ax.set_title(sig.replace("_", " "), fontsize=10,
                                 fontweight="bold")
                if j == 0:
                    resp_short = "R" if "R" in sname and "NR" not in sname else "NR"
                    ax.set_ylabel(f"{sname} ({resp_short})", fontsize=10,
                                  fontweight="bold")

        tag_label = "Responder" if tag == "Responder" else "Non-Responder"
        fig.suptitle(f"GSE238264: Spatial IO Signatures — {tag_label} Sections",
                     fontsize=14, fontweight="bold", color="black", y=1.01)
        fig.patch.set_facecolor("white")
        plt.tight_layout()
        fname = f"02b_visium_spatial_signature_maps_{tag.lower()}.png"
        fig.savefig(FIG / fname, dpi=250, bbox_inches="tight")
        plt.close()
        print(f"  Saved {fname}")


# ═══════════════════════════════════════════════════════════════════════
#  FIG 7: Cell-Type Marker Spatial  (improved)
# ═══════════════════════════════════════════════════════════════════════

def fig_celltype_spatial_scores(adata, adatas):
    """Improved cell-type marker spatial maps for R vs NR representative."""
    print("\n[Fig 7] Cell-type marker spatial maps...")

    ct_cols = ["ct_CD8_T", "ct_Macrophage", "ct_Tumor", "ct_B_cell"]
    ct_labels = ["CD8+ T cells", "Macrophage", "Tumor", "B cells"]
    r_name, nr_name = "HCC1R", "HCC7NR"

    for name, ad in adatas.items():
        sc.pp.filter_genes(ad, min_cells=1)
        avail = set(ad.var_names)
        for ct_name, genes in CT_MARKERS_OLD.items():
            present = [g for g in genes if g in avail]
            col = f"ct_{ct_name}"
            if len(present) >= 2 and col not in ad.obs.columns:
                sc.pp.normalize_total(ad, target_sum=1e4)
                sc.pp.log1p(ad)
                sc.tl.score_genes(ad, present, score_name=col, use_raw=False)

    n_ct = len(ct_cols)
    cmaps = ["YlOrRd", "BuGn", "RdPu", "BuPu"]

    fig, axes = plt.subplots(2, n_ct, figsize=(5.5 * n_ct, 10))

    for row_idx, (sname, label) in enumerate([(r_name, "Responder"),
                                               (nr_name, "Non-Responder")]):
        ad = adatas[sname]
        coords = ad.obsm["spatial"]
        for col_idx, (ct_col, ct_label) in enumerate(zip(ct_cols, ct_labels)):
            ax = axes[row_idx, col_idx]
            ax.set_facecolor("white")

            if ct_col in ad.obs.columns:
                vals = ad.obs[ct_col].values
                vmin, vmax = np.percentile(vals, [2, 98])
                if vmax - vmin < 0.05:
                    vmin, vmax = vals.min(), vals.max()
                sc_p = ax.scatter(coords[:, 0], coords[:, 1],
                                  c=vals, cmap=cmaps[col_idx], s=15,
                                  vmin=vmin, vmax=vmax,
                                  edgecolors="none", rasterized=True, alpha=0.9)
                cbar = plt.colorbar(sc_p, ax=ax, shrink=0.65, pad=0.02)
                cbar.ax.tick_params(labelsize=8)
            ax.set_aspect("equal"); ax.invert_yaxis()
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_linewidth(0.5); sp.set_color("#CCCCCC")
            if row_idx == 0:
                ax.set_title(ct_label, fontsize=12, fontweight="bold")
            if col_idx == 0:
                title_color = RESPONSE_PAL[SAMPLES[sname]]
                ax.set_ylabel(f"{sname}\n({label})", fontsize=11,
                              fontweight="bold", color=title_color)

    fig.suptitle("GSE238264: Cell-Type Marker Scores — Responder vs Non-Responder",
                 fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    fig.savefig(FIG / "02b_visium_celltype_spatial_maps.png", dpi=250, bbox_inches="tight")
    plt.close()
    print("  Saved.")


# ═══════════════════════════════════════════════════════════════════════
#  FIG 8: Neighbourhood Enrichment  (fixed labels)
# ═══════════════════════════════════════════════════════════════════════

def fig_nhood_enrichment(adata):
    """Neighbourhood enrichment with proper labels and no overlap."""
    print("\n[Fig 8] Neighbourhood enrichment...")
    if "spatial_immunotype" not in adata.obs.columns:
        print("  [SKIP] no spatial_immunotype")
        return

    # Build spatial graph on combined
    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)
    sq.gr.nhood_enrichment(adata, cluster_key="spatial_immunotype")

    # Extract enrichment matrix from adata.uns
    nhood_key = "spatial_immunotype_nhood_enrichment"
    if nhood_key not in adata.uns:
        print("  [WARN] nhood enrichment key not found")
        return

    zscore = adata.uns[nhood_key]["zscore"]
    cats = adata.obs["spatial_immunotype"].cat.categories.tolist()
    # Remove categories with 0 observations
    cat_counts = adata.obs["spatial_immunotype"].value_counts()
    valid_cats = [c for c in cats if cat_counts.get(c, 0) > 0]
    valid_idx = [cats.index(c) for c in valid_cats]
    zscore_sub = zscore[np.ix_(valid_idx, valid_idx)]

    IMTYPE_PAL = {"Immune_Hot": "#D62728", "Immune_Excluded": "#FF7F0E",
                  "Immune_Desert": "#1F77B4", "Immune_Mixed": "#2CA02C"}

    fig, ax = plt.subplots(figsize=(8, 7))
    vabs = max(abs(zscore_sub.min()), abs(zscore_sub.max()))
    if vabs < 0.1:
        vabs = 1.0  # avoid degenerate colorbar
    im = ax.imshow(zscore_sub, cmap="RdBu_r", vmin=-vabs, vmax=vabs, aspect="equal")
    cbar = plt.colorbar(im, ax=ax, shrink=0.75, pad=0.04)
    cbar.set_label("Enrichment Z-score", fontsize=11)
    cbar.ax.tick_params(labelsize=10)

    # Annotate cells with z-scores
    for i in range(len(valid_cats)):
        for j in range(len(valid_cats)):
            val = zscore_sub[i, j]
            color = "white" if abs(val) > vabs * 0.5 else "black"
            ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                    fontsize=10, fontweight="bold", color=color)

    ax.set_xticks(range(len(valid_cats)))
    ax.set_yticks(range(len(valid_cats)))
    labels_clean = [c.replace("Immune_", "") for c in valid_cats]
    ax.set_xticklabels(labels_clean, fontsize=11, fontweight="bold",
                       rotation=35, ha="right")
    ax.set_yticklabels(labels_clean, fontsize=11, fontweight="bold")

    # Colour the tick labels
    for tick, cat in zip(ax.get_xticklabels(), valid_cats):
        tick.set_color(IMTYPE_PAL.get(cat, "black"))
    for tick, cat in zip(ax.get_yticklabels(), valid_cats):
        tick.set_color(IMTYPE_PAL.get(cat, "black"))

    ax.set_xlabel("spatial_immunotype", fontsize=11, labelpad=8)
    ax.set_ylabel("spatial_immunotype", fontsize=11, labelpad=8)
    ax.set_title("Neighbourhood Enrichment: Spatial Immunotypes",
                 fontsize=13, fontweight="bold", pad=14)
    plt.tight_layout()
    fig.savefig(FIG / "02b_visium_nhood_enrichment.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved.")


# ═══════════════════════════════════════════════════════════════════════
#  FIG 9: Marker Dot Plot  (improved dimensions and contrast)
# ═══════════════════════════════════════════════════════════════════════

def fig_marker_dotplot(adata):
    print("\n[Fig 9] Marker dot plot (improved)...")
    avail = set(adata.var_names)
    marker_list = []
    type_groups = {}
    for ct_name, genes in CELL_MARKERS.items():
        present = [g for g in genes if g in avail][:5]
        marker_list.extend(present)
        for g in present:
            type_groups[g] = ct_name
    marker_list = list(dict.fromkeys(marker_list))

    valid_cts = [ct for ct in CELL_MARKERS.keys()
                 if (adata.obs["celltype"] == ct).sum() >= 5]

    adata_raw = adata.raw if adata.raw is not None else adata
    rows = []
    for ct in valid_cts:
        m = adata.obs["celltype"] == ct
        sub = adata_raw[m]
        for gene in marker_list:
            if gene not in adata_raw.var_names:
                continue
            if issparse(sub[:, gene].X):
                vals = np.asarray(sub[:, gene].X.todense()).flatten()
            else:
                vals = sub[:, gene].X.flatten()
            rows.append({"celltype": ct, "gene": gene,
                         "mean_expr": np.mean(vals),
                         "frac_expr": np.mean(vals > 0)})
    df = pd.DataFrame(rows)
    if df.empty:
        print("  [WARN] no data"); return

    mean_piv = df.pivot_table(index="celltype", columns="gene", values="mean_expr",
                              aggfunc="first").reindex(index=valid_cts, columns=marker_list).fillna(0)
    frac_piv = df.pivot_table(index="celltype", columns="gene", values="frac_expr",
                              aggfunc="first").reindex(index=valid_cts, columns=marker_list).fillna(0)
    zscore_piv = (mean_piv - mean_piv.mean(axis=0)) / (mean_piv.std(axis=0) + 1e-8)
    zscore_piv = zscore_piv.clip(-2.5, 2.5)

    fig, ax = plt.subplots(figsize=(max(16, 0.5 * len(marker_list)),
                                    max(6, 0.65 * len(valid_cts))))
    ax.set_facecolor("#FAFAFA")

    cmap = plt.cm.RdYlBu_r
    norm = plt.Normalize(-2.5, 2.5)

    for i, ct in enumerate(valid_cts):
        for j, gene in enumerate(marker_list):
            z = zscore_piv.loc[ct, gene] if gene in zscore_piv.columns else 0
            frac = frac_piv.loc[ct, gene] if gene in frac_piv.columns else 0
            size = max(frac * 350, 5)
            ax.scatter(j, i, s=size, c=[cmap(norm(z))],
                       edgecolors="#333333", linewidth=0.4, zorder=3)

    ax.set_xticks(range(len(marker_list)))
    ax.set_xticklabels(marker_list, rotation=90, fontsize=8, fontstyle="italic")
    ax.set_yticks(range(len(valid_cts)))
    ax.set_yticklabels(valid_cts, fontsize=10)

    for tick_l in ax.get_yticklabels():
        ct_name = tick_l.get_text()
        tick_l.set_color(CELLTYPE_PAL.get(ct_name, "black"))
        tick_l.set_fontweight("bold")

    # Gene group separators
    prev_group = None
    for j, gene in enumerate(marker_list):
        grp = type_groups.get(gene, "")
        if grp != prev_group and prev_group is not None:
            ax.axvline(j - 0.5, color="#AAAAAA", linewidth=0.6, linestyle="--")
        prev_group = grp

    ax.set_xlim(-0.5, len(marker_list) - 0.5)
    ax.set_ylim(-0.5, len(valid_cts) - 0.5)
    ax.invert_yaxis()
    ax.set_title("Cell-Type Marker Expression — Dot Plot\n"
                 "Size = fraction expressing · Colour = z-scored mean",
                 fontsize=12, fontweight="bold")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label("Z-scored expression", fontsize=10)

    for frac_val, label in [(0.2, "20%"), (0.5, "50%"), (0.8, "80%")]:
        ax.scatter([], [], s=frac_val * 350, c="gray", edgecolors="black",
                   linewidth=0.4, label=f"{label}")
    ax.legend(fontsize=8, loc="lower right", framealpha=0.95, title="Fraction expressing",
              title_fontsize=9, edgecolor="#CCCCCC")

    plt.tight_layout()
    fig.savefig(FIG / "02b_ct_marker_dotplot.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved.")


# ═══════════════════════════════════════════════════════════════════════
#  FIG 10: Marker Clustered Heatmap  (improved)
# ═══════════════════════════════════════════════════════════════════════

def fig_marker_heatmap(adata):
    print("\n[Fig 10] Marker clustered heatmap (improved)...")
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
        m = adata.obs["celltype"] == ct
        sub = adata_raw[m]
        for j, gene in enumerate(marker_list):
            if gene not in adata_raw.var_names: continue
            if issparse(sub[:, gene].X):
                vals = np.asarray(sub[:, gene].X.todense()).flatten()
            else:
                vals = sub[:, gene].X.flatten()
            mat[i, j] = np.mean(vals)

    mat_z = (mat - mat.mean(axis=0)) / (mat.std(axis=0) + 1e-8)
    mat_z = np.clip(mat_z, -2.5, 2.5)
    df_z = pd.DataFrame(mat_z, index=valid_cts, columns=marker_list)

    row_colors = [CELLTYPE_PAL.get(ct, "#999999") for ct in valid_cts]

    g = sns.clustermap(
        df_z, cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5,
        row_cluster=True, col_cluster=False,
        row_colors=row_colors,
        linewidths=0.5, linecolor="white",
        figsize=(max(14, 0.4 * len(marker_list)), max(5, 0.6 * len(valid_cts))),
        dendrogram_ratio=(0.12, 0),
        cbar_kws={"label": "Z-score", "shrink": 0.6},
        xticklabels=True, yticklabels=True,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(),
                                  rotation=90, fontsize=8, fontstyle="italic")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=10)
    g.fig.suptitle("Cell-Type Marker Expression — Clustered Heatmap\n"
                   "(Z-scored mean per annotated type)",
                   fontsize=13, fontweight="bold", y=1.05)
    g.savefig(FIG / "02b_ct_marker_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved.")


# ═══════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("02b — Improve All Spatial Plots & UMAPs")
    print("=" * 70)

    # Load sections
    print("\n[0] Loading data...")
    adatas = load_all_sections()

    # Merge + embed
    print("\n[0b] Merging + embedding...")
    adata = merge_and_embed(adatas)

    # Annotate
    print("\n[0c] Annotating cell types + scoring...")
    adata, score_cols = annotate(adata)

    counts = adata.obs["celltype"].value_counts()
    print("\n  Cell-type composition:")
    for ct, n in counts.items():
        print(f"    {ct:20s}: {n:6d} ({100*n/adata.n_obs:5.1f}%)")

    # Generate all improved figures
    fig_umap_celltype_response(adata)
    fig_umap_sample(adata)
    fig_umap_leiden(adata)
    fig_umap_immune(adata)
    fig_spatial_celltypes(adata)
    fig_io_signature_maps(adata, adatas)
    fig_celltype_spatial_scores(adata, adatas)
    fig_nhood_enrichment(adata)
    fig_marker_dotplot(adata)
    fig_marker_heatmap(adata)

    print("\n" + "=" * 70)
    print("[DONE] All figures regenerated with improved contrast & dimensions")
    print(f"  Output: {FIG}")
    print("=" * 70)


if __name__ == "__main__":
    main()
