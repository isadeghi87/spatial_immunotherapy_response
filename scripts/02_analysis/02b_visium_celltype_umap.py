#!/usr/bin/env python3
"""
02b_visium_celltype_umap.py — Clustering, Cell-Type Annotation & UMAP for Visium ST

1. Merge all GSE238264 sections → PCA → UMAP
2. Leiden clustering at multiple resolutions
3. Cell-type annotation via cell2location with atlas-derived reference signatures
4. Figures:
   a) UMAP coloured by annotated cell type (all spots)
   b) UMAP coloured by sample / response
   c) Heatmap of cell-type marker genes (dot plot style)
   d) UMAP of immune-subtype spots only
   e) Spatial maps coloured by cell-type annotation
"""

import warnings
warnings.filterwarnings("ignore")

import os
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
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
DATA    = PROJECT / "data" / "spatial_data" / "GSE238264_HCC_IO_Visium"
ATLAS   = PROJECT / "data" / "processed" / "integrated_io_atlas.h5ad"
OUT     = PROJECT / "out" / "02b_spatial_biomarkers"
FIG     = OUT / "figures"
TAB     = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

C2L_MAX_EPOCHS = int(os.environ.get("C2L_MAX_EPOCHS", "30000"))
# Use C2L_MAX_EPOCHS=3000 in the environment for a quick CPU smoke test.
# Production runs default to 30,000 with early stopping/checkpoint reuse.
C2L_REFERENCE_MAX_CELLS_PER_TYPE = 4_000

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
RESPONSE_PAL = {"Responder": "#0077BB", "Non_Responder": "#CC3311"}

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

SYMLINK_ROOT_FALLBACKS = {
    "/mnt/projects": "/data/projects",
}


def _resolve_visium_section_path(path):
    if path.exists():
        return path

    if not path.is_symlink():
        raise FileNotFoundError(f"Required Visium section not found: {path}")

    raw_target = Path(os.readlink(path))
    candidates = []
    if raw_target.is_absolute():
        candidates.append(raw_target)
    else:
        candidates.append((path.parent / raw_target).resolve())

    target_str = str(raw_target)
    for src_root, dst_root in SYMLINK_ROOT_FALLBACKS.items():
        if target_str == src_root or target_str.startswith(f"{src_root}/"):
            candidates.append(Path(target_str.replace(src_root, dst_root, 1)))

    for candidate in candidates:
        if candidate.exists():
            if candidate != path:
                print(f"  Resolved broken symlink via fallback: {path} -> {candidate}")
            return candidate

    candidate_text = ", ".join(str(candidate) for candidate in candidates)
    raise FileNotFoundError(
        f"Broken Visium symlink: {path} -> {raw_target}. Tried: {candidate_text}"
    )


def load_and_merge():
    """Load each section, merge into one AnnData."""
    adatas = []
    for name, resp in SAMPLES.items():
        path = DATA / f"{name}.h5ad"
        resolved = _resolve_visium_section_path(path)
        print(f"  Loading {name} ({resp})...")
        ad = sc.read_h5ad(str(resolved))
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

    # Neighbors & UMAP
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30, use_rep="X_pca")
    sc.tl.umap(adata, min_dist=0.3, spread=1.2)
    adata.obsm["X_umap_expr"] = adata.obsm["X_umap"].copy()
    adata.uns["umap_basis_notes"] = {
        "X_umap_expr": "Expression/PCA UMAP computed before cell2location.",
    }
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

def _sanitize_celltype_name(name):
    return (str(name)
            .replace("+", "p")
            .replace("/", "_")
            .replace(" ", "_")
            .replace("-", "_"))


def _build_cell2location_reference():
    """Build a compact atlas-derived reference signature matrix for cell2location."""
    print("  Building atlas reference signatures for cell2location...")
    ref_csv = TAB / "gse238264_cell2location_reference_signatures.csv"
    if ref_csv.exists():
        print(f"  [CHECKPOINT] Loading pre-computed reference: {ref_csv}")
        reference = pd.read_csv(ref_csv, index_col=0)
        print(f"  Reference matrix: {reference.shape[0]} genes × {reference.shape[1]} cell types")
        return reference
    if not ATLAS.exists():
        raise FileNotFoundError(f"Atlas not found: {ATLAS}")

    atlas = sc.read_h5ad(ATLAS, backed="r")
    try:
        if "cell_type" not in atlas.obs.columns:
            raise KeyError("integrated atlas is missing 'cell_type' annotations")

        obs = atlas.obs[["cell_type"]].copy()
        obs["cell_type"] = obs["cell_type"].astype(str)
        obs["pos"] = np.arange(len(obs))  # integer row indices for backed-mode subsetting
        obs = obs[~obs["cell_type"].isin(["nan", "None", "Unknown", "Unassigned"])]

        rng = np.random.default_rng(0)
        selected_pos = []
        for cell_type, group in obs.groupby("cell_type", observed=True):
            pos = group["pos"].to_numpy()
            if len(pos) == 0:
                continue
            n_keep = min(len(pos), C2L_REFERENCE_MAX_CELLS_PER_TYPE)
            selected_pos.append(rng.choice(pos, n_keep, replace=False))

        if not selected_pos:
            raise RuntimeError("No reference cells available for cell2location")

        selected_pos = np.sort(np.unique(np.concatenate(selected_pos)))
        subset = atlas[selected_pos].to_memory()
    finally:
        atlas.file.close()

    layer = subset.layers["counts"] if "counts" in subset.layers else subset.X
    ref_rows = {}
    for cell_type in sorted(subset.obs["cell_type"].astype(str).unique()):
        mask = (subset.obs["cell_type"].astype(str) == cell_type).to_numpy()
        if mask.sum() < 20:
            continue
        matrix = layer[mask]
        mean_expr = np.asarray(matrix.mean(axis=0)).ravel()
        ref_rows[cell_type] = mean_expr

    reference = pd.DataFrame(ref_rows, index=subset.var_names)
    reference = reference.loc[reference.sum(axis=1) > 0]
    reference.to_csv(TAB / "gse238264_cell2location_reference_signatures.csv")
    print(f"  Reference matrix: {reference.shape[0]} genes × {reference.shape[1]} cell types")
    return reference


def annotate_celltypes(adata):
    """Annotate Visium spots with cell2location using atlas-derived signatures."""
    print("\n[2] Cell-type annotation with cell2location...")
    try:
        from cell2location.models import Cell2location
    except ImportError as exc:
        raise ImportError(
            "cell2location is required for Visium annotation. Install it in the analysis environment first."
        ) from exc

    reference = _build_cell2location_reference()
    shared_genes = adata.var_names.intersection(reference.index)
    if len(shared_genes) < 500:
        raise RuntimeError(f"Too few shared genes for cell2location: {len(shared_genes)}")

    adata = adata[:, shared_genes].copy()
    reference = reference.loc[shared_genes]

    abundance_csv = TAB / "gse238264_cell2location_abundances.csv"
    if abundance_csv.exists():
        print(f"  [CHECKPOINT] Loading pre-computed cell2location abundance: {abundance_csv}")
        abundance = pd.read_csv(abundance_csv, index_col=0)
        abundance = abundance.reindex(adata.obs_names).fillna(0.0)
        abundance_key = "checkpoint"
    else:
        if "counts" not in adata.layers:
            adata.layers["counts"] = adata.X.copy()

        Cell2location.setup_anndata(adata, layer="counts", batch_key="sample")
        # Force CPU training — GPU may be occupied by other jobs (e.g. Geneformer)
        import os as _os
        _os.environ["CUDA_VISIBLE_DEVICES"] = ""
        model = Cell2location(
            adata,
            cell_state_df=reference,
            N_cells_per_location=30,
            detection_alpha=20,
        )
        print(f"  Training cell2location for up to {C2L_MAX_EPOCHS:,} epochs (CPU, early stopping)...")
        model.train(
            max_epochs=C2L_MAX_EPOCHS,
            batch_size=2048,
            train_size=1,
            accelerator="cpu",
        )
        adata = model.export_posterior(
            adata,
            sample_kwargs={"num_samples": 1000, "batch_size": min(2048, adata.n_obs)},
        )

        abundance_key = None
        for candidate in ("q05_cell_abundance_w_sf", "means_cell_abundance_w_sf"):
            if candidate in adata.obsm:
                abundance_key = candidate
                break
        if abundance_key is None:
            raise RuntimeError("cell2location posterior abundance matrix not found in adata.obsm")

        abundance = adata.obsm[abundance_key]
        if not isinstance(abundance, pd.DataFrame):
            abundance = pd.DataFrame(abundance, index=adata.obs_names, columns=reference.columns)
        abundance = abundance.loc[adata.obs_names]
        abundance.columns = [str(c) for c in abundance.columns]
        abundance.to_csv(abundance_csv)

    # Normalize column names to match CELLTYPE_PAL keys.
    # Atlas cell_type labels may differ from display names (e.g. B_cell → B cell).
    # cell2location Export adds a prefix like "q05cell_abundance_w_sf_" — strip it first.
    import re as _re
    _prefix_pat = _re.compile(
        r'^q05cell_abundance_w_sf_|^q05_cell_abundance_w_sf_'
        r'|^means_cell_abundance_w_sf_|^means_cell_abundance_'
    )
    abundance.columns = [_prefix_pat.sub('', c) for c in abundance.columns]
    _CT_RENAME = {
        "B_cell": "B cell",
        "CD4_T":  "CD4+ T cell",
        "CD8_T":  "CD8+ T cell",
        "NK":     "NK cell",
    }
    # Merge DC subtypes to a single "DC" column
    _dc_cols = [c for c in abundance.columns if c in ("cDC1", "cDC2", "mregDC")]
    if _dc_cols:
        if "DC" not in abundance.columns:
            abundance.insert(0, "DC", abundance[_dc_cols].sum(axis=1))
        else:
            abundance["DC"] = abundance["DC"] + abundance[_dc_cols].sum(axis=1)
        abundance = abundance.drop(columns=_dc_cols)
    abundance = abundance.rename(columns={k: v for k, v in _CT_RENAME.items() if k in abundance.columns})

    for cell_type in abundance.columns:
        adata.obs[f"c2l_{_sanitize_celltype_name(cell_type)}"] = abundance[cell_type].values

    if "X_umap" in adata.obsm and "X_umap_expr" not in adata.obsm:
        adata.obsm["X_umap_expr"] = adata.obsm["X_umap"].copy()
    adata.obsm["X_cell2location"] = abundance.values.astype(np.float32)
    sc.pp.neighbors(adata, use_rep="X_cell2location", n_neighbors=20)
    sc.tl.umap(adata, min_dist=0.3, spread=1.15)
    adata.obsm["X_umap_cell2location"] = adata.obsm["X_umap"].copy()
    adata.uns.setdefault("umap_basis_notes", {})
    adata.uns["umap_basis_notes"].update({
        "X_umap_cell2location": "UMAP of continuous cell2location abundance vectors.",
        "X_umap": "Alias of X_umap_cell2location for legacy plotting functions.",
    })

    # Use only columns present in CELLTYPE_PAL for label assignment
    celltype_cols = [c for c in CELLTYPE_PAL if c != "Unassigned" and c in abundance.columns]
    values = abundance[celltype_cols].to_numpy(dtype=float)
    top_idx = np.argmax(values, axis=1)
    top_val = values[np.arange(values.shape[0]), top_idx]
    second_val = np.partition(values, -2, axis=1)[:, -2] if values.shape[1] > 1 else np.zeros(values.shape[0])
    margin = top_val - second_val
    labels = np.array(celltype_cols)[top_idx].astype(object)
    labels[(top_val < 1.0) | (margin < 0.5)] = "Unassigned"

    adata.obs["celltype"] = pd.Categorical(labels, categories=list(CELLTYPE_PAL.keys()))
    adata.obs["celltype_visual_label"] = adata.obs["celltype"].astype(str)
    adata.obs["celltype_confident"] = top_val >= 1.0
    adata.obs["celltype_top_abundance"] = top_val
    adata.obs["celltype_margin"] = margin
    adata.uns["celltype_assignment_method"] = "cell2location_dominant_visual_label"
    adata.uns["celltype_label_warning"] = (
        "celltype is the dominant cell2location abundance label for visualization only; "
        "quantitative analyses must use continuous c2l_* abundance columns."
    )
    adata.uns["quantitative_spatial_features"] = [
        f"c2l_{_sanitize_celltype_name(cell_type)}" for cell_type in abundance.columns
    ]
    adata.uns["cell2location_params"] = {
        "max_epochs": C2L_MAX_EPOCHS,
        "reference_max_cells_per_type": C2L_REFERENCE_MAX_CELLS_PER_TYPE,
        "n_cells_per_location": 30,
        "detection_alpha": 20,
        "abundance_key": abundance_key,
    }

    counts = adata.obs["celltype"].value_counts()
    print("\n  Cell-type composition:")
    for ct, n in counts.items():
        pct = 100 * n / adata.n_obs
        print(f"    {ct:20s}: {n:6d} ({pct:5.1f}%)")

    avail = set(adata.var_names)
    for sig_name, genes in IMMUNE_SIGS.items():
        present = [g for g in genes if g in avail]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=sig_name, use_raw=False)

    return adata, []


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

    if not valid_cts:
        print("  [SKIP] No cell types with ≥5 spots — clustered heatmap skipped")
        return

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
    if "X_cell2location" in adata_imm.obsm:
        use_rep = "X_cell2location"
    elif "X_pca" in adata_imm.obsm:
        use_rep = "X_pca"
    else:
        use_rep = None
    if use_rep is not None:
        sc.pp.neighbors(adata_imm, n_neighbors=15, use_rep=use_rep)
    else:
        sc.pp.neighbors(adata_imm, n_neighbors=15)
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
    """Spatial maps coloured by annotated cell type, one per section.

    HCC_omics reference style: dark background, bold white titles,
    rasterized point overlay, shared legend panel.
    """
    print("\n[3e] Spatial cell-type maps...")

    samples = list(SAMPLES.keys())
    fig, axes = plt.subplots(2, 4, figsize=(24, 12),
                              facecolor="#111111")
    axes_flat = axes.flatten()

    for idx, sname in enumerate(samples):
        ax = axes_flat[idx]
        ax.set_facecolor("#1a1a1a")
        mask = adata.obs["sample"] == sname
        sub = adata[mask]
        coords = sub.obsm["spatial"]
        resp = sub.obs["response"].iloc[0]

        # Plot non-immune first (background), immune on top
        base_order = ["Unassigned", "Hepatocyte", "Endothelial", "Fibroblast", "Tumor"]
        top_order = ["B cell", "Macrophage", "DC", "NK cell",
                     "CD4+ T cell", "Treg", "CD8+ T cell"]
        for ct in base_order + top_order:
            ct_mask = sub.obs["celltype"] == ct
            if ct_mask.sum() == 0:
                continue
            color = CELLTYPE_PAL[ct]
            size = 3 if ct in ("Unassigned", "Hepatocyte") else 7
            alpha = 0.55 if ct in ("Unassigned", "Hepatocyte") else 0.9
            ax.scatter(coords[ct_mask.values, 0], coords[ct_mask.values, 1],
                       c=color, s=size, edgecolors="none",
                       rasterized=True, alpha=alpha)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        resp_short = "R" if resp == "Responder" else "NR"
        ax.set_title(f"{sname} ({resp_short})\n{mask.sum()} spots",
                     fontsize=11, fontweight="bold",
                     color="white")
        for spine in ax.spines.values():
            spine.set_visible(False)

    # Last subplot: legend on dark background
    ax = axes_flat[-1]
    ax.set_facecolor("#1a1a1a")
    ax.axis("off")
    for ct, color in CELLTYPE_PAL.items():
        ax.scatter([], [], c=color, s=80, label=ct, edgecolors="none")
    leg = ax.legend(fontsize=10, loc="center", framealpha=0.0, markerscale=1.5,
                    title="Cell Type", title_fontsize=11,
                    labelcolor="white")
    plt.setp(leg.get_title(), color="white", fontweight="bold")

    plt.suptitle("Spatial Cell-Type Annotation — GSE238264 HCC Visium",
                 fontsize=14, fontweight="bold", color="white", y=0.995)
    plt.tight_layout()
    plt.savefig(FIG / "02b_ct_spatial_celltype_maps.png", dpi=250,
                bbox_inches="tight", facecolor="#111111")
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
    cols += [c for c in adata.obs.columns if c.startswith("c2l_")]
    cols += [c for c in ["celltype_confident", "celltype_top_abundance", "celltype_margin"] if c in adata.obs.columns]
    df = adata.obs[cols].copy()
    df.to_csv(TAB / "gse238264_celltype_annotations.csv")
    print(f"  Saved {df.shape[0]} spots × {df.shape[1]} columns")

    # Per-sample composition
    comp = adata.obs.groupby(["sample", "response", "celltype"]).size().reset_index(name="n_spots")
    comp.to_csv(TAB / "gse238264_celltype_composition.csv", index=False)
    print(f"  Composition table: {comp.shape[0]} rows")

    support = (adata.obs.assign(confident=lambda d: d["celltype"] != "Unassigned")
               .groupby(["sample", "response"], observed=True)
               .agg(
                   n_spots=("celltype", "size"),
                   confident_fraction=("confident", "mean"),
                   dominant_celltype=("celltype", lambda x: x.value_counts().index[0]),
                   dominant_fraction=("celltype", lambda x: x.value_counts(normalize=True).iloc[0]),
               )
               .reset_index())
    support.to_csv(TAB / "gse238264_sample_annotation_support.csv", index=False)
    print(f"  Sample support table: {support.shape[0]} rows")


def plot_sample_annotation_support(adata):
    """Summarize per-sample annotation support and dominant cell types."""
    print("\n[3g] Sample-level annotation support...")
    comp = (adata.obs.groupby(["sample", "response", "celltype"], observed=True)
            .size()
            .rename("n_spots")
            .reset_index())
    totals = comp.groupby("sample", observed=True)["n_spots"].sum().rename("total_spots")
    comp = comp.merge(totals, on="sample", how="left")
    comp["fraction"] = comp["n_spots"] / comp["total_spots"]

    top_celltypes = (comp.groupby("celltype", observed=True)["n_spots"]
                     .sum()
                     .sort_values(ascending=False)
                     .head(8)
                     .index.tolist())
    plot_df = comp[comp["celltype"].isin(top_celltypes)].copy()
    sample_order = list(SAMPLES.keys())

    support = pd.read_csv(TAB / "gse238264_sample_annotation_support.csv")

    fig, axes = plt.subplots(1, 2, figsize=(18, 6), gridspec_kw={"width_ratios": [1.4, 1]})

    ax = axes[0]
    pivot = (plot_df.pivot_table(index="sample", columns="celltype", values="fraction", fill_value=0)
             .reindex(sample_order)
             .fillna(0))
    bottom = np.zeros(len(pivot))
    for celltype in pivot.columns:
        ax.bar(pivot.index, pivot[celltype].values, bottom=bottom,
               color=CELLTYPE_PAL.get(celltype, "#999999"), label=celltype,
               edgecolor="white", linewidth=0.7)
        bottom += pivot[celltype].values
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of spots")
    ax.set_title("A. Per-sample cell2location annotation composition", fontweight="bold")
    ax.tick_params(axis="x", rotation=35)
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)

    ax = axes[1]
    support = support.sort_values("sample")
    colors = [RESPONSE_PAL.get(str(resp).replace("-", "_"), "#666666") for resp in support["response"]]
    ax.bar(support["sample"], support["confident_fraction"], color=colors,
           edgecolor="white", linewidth=0.8)
    for idx, row in support.iterrows():
        ax.text(idx, min(row["confident_fraction"] + 0.03, 1.02),
                f"{row['dominant_celltype']}\n{row['dominant_fraction']:.0%}",
                ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Confidently assigned fraction")
    ax.set_title("B. Sample-level support and dominant cell type", fontweight="bold")
    ax.tick_params(axis="x", rotation=35)

    plt.suptitle("GSE238264 cell2location annotation support", fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(FIG / "02b_ct_sample_annotation_support.png", dpi=250, bbox_inches="tight")
    plt.close()
    print("  Saved.")


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
    plot_sample_annotation_support(adata)

    print("\n" + "=" * 70)
    print("[DONE] All cell-type annotation & UMAP figures generated")
    print(f"  Figures: {FIG}")
    print(f"  Tables:  {TAB}")
    print("=" * 70)


if __name__ == "__main__":
    main()
