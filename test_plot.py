#!/usr/bin/env python3
"""
01_preprocess_integrate.py — Preprocess & Integrate All HCC IO Datasets

Loads, QCs, normalizes, and integrates:
  - GSE206325 (Magen 2023 Nat Med): 682k cells, anti-PD1 R vs NR, scRNA+TCR
  - GSE125449 (Ma 2019): 10k cells, HCC landscape
  - GSE140228 (Zhang 2019 Cell): CD45+ immune cells, HCC immune landscape
  - GSE149614 (Lu 2022 Nat Commun): 70k cells, HCC primary/metastatic
  - E-MTAB-10553 (Pfister 2021): ~60k cells, NASH-HCC + anti-PD1
  - GSE202069 (Zhu 2023): bulk RNA-seq (processed separately)

Pipeline:
  1. Load each dataset → standardize gene names & metadata
  2. QC filtering (genes/cell, mito%, doublets)
  3. Concatenate on shared genes
  4. Normalize, HVG, PCA
  5. Harmony batch integration
  6. UMAP, Leiden clustering
  7. Transfer Magen cluster annotations (cell types) to full atlas
  8. Save integrated h5ad
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
from scipy import sparse
from scipy.io import mmread

sc.settings.verbosity = 1
sc.settings.n_jobs = 8

PROJECT = Path(__file__).resolve().parents[2]
DATA    = PROJECT / "data"
PROC    = DATA / "processed"
OUT     = PROJECT / "out" / "01_preprocessing"
FIG     = OUT / "figures"
TAB     = OUT / "tables"
for d in [FIG, TAB, PROC]:
    d.mkdir(parents=True, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════
# 1.  LOAD DATASETS
# ═══════════════════════════════════════════════════════════════════════

def load_gse206325():
    """Load Magen et al. 2023 — 682k cells, anti-PD1 R vs NR."""
    print("[GSE206325] Loading Magen 2023 (Nat Med)...")
    d = PROC

    # Barcodes & genes
    barcodes = pd.read_csv(d / "GSE206325_barcodes.tsv", header=None, sep="\t")[0].values
    genes = pd.read_csv(d / "GSE206325_genes.tsv", header=None, sep="\t")
    gene_names = genes[1].values if genes.shape[1] >= 2 else genes[0].values

    # Matrix
    mat = mmread(d / "GSE206325_matrix.mtx.gz").T.tocsr()
    adata = sc.AnnData(X=mat, obs=pd.DataFrame(index=barcodes),
                       var=pd.DataFrame(index=gene_names))
    adata.var_names_make_unique()

    # Cell metadata (sample_ID, cluster)
    cell_meta = pd.read_csv(d / "GSE206325_metadata.csv", index_col=0)
    adata.obs["sample_id"] = cell_meta.loc[adata.obs_names, "cell_to_sample_ID"].values.astype(str)
    adata.obs["cluster_id"] = cell_meta.loc[adata.obs_names, "cell_to_cluster"].values.astype(str)

    # Sample annotations → map to cells
    sample_ann = pd.read_csv(d / "GSE206325_sample_annotations.csv", index_col=0)
    sample_ann.index = sample_ann.index.astype(str)
    for col in ["patient_ID", "tissue", "treatment_Resp", "Immune_Infiltration"]:
        mapping = sample_ann[col].to_dict()
        adata.obs[col] = adata.obs["sample_id"].map(mapping)

    # Derive response label
    resp_map = {"anti-PD1_R": "Responder", "anti-PD1_NR": "Non-Responder"}
    adata.obs["response"] = adata.obs["treatment_Resp"].map(resp_map)

    # Cluster annotations → cell type
    clust_ann = pd.read_csv(d / "GSE206325_cluster_annotations.csv", index_col=0)
    clust_ann["cluster_ID"] = clust_ann["cluster_ID"].astype(str)
    ct_map = {}
    for _, row in clust_ann.iterrows():
        cid = row["cluster_ID"]
        subgroup = row["subgroup"] if pd.notna(row["subgroup"]) else ""
        ct = row["type"] if pd.notna(row["type"]) else ""
        group = row["group"] if pd.notna(row["group"]) else ""
        # Build hierarchical label
        if subgroup:
            ct_map[cid] = f"{ct}_{subgroup}" if ct else subgroup
        elif ct:
            ct_map[cid] = ct
        else:
            ct_map[cid] = group
    adata.obs["cell_type_original"] = adata.obs["cluster_id"].map(ct_map)

    adata.obs["dataset"] = "GSE206325"
    adata.obs["study"] = "Magen2023"
    print(f"  {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    print(f"  Response: {adata.obs['response'].value_counts().to_dict()}")
    print(f"  Tissues: {adata.obs['tissue'].value_counts().to_dict()}")
    return adata


def load_gse125449():
    """Load Ma et al. 2019 — HCC scRNA-seq landscape."""
    print("[GSE125449] Loading Ma 2019...")
    d = DATA / "pre_treatment" / "GSE125449_Ma2019_antiPD1"
    datasets = []
    for set_name in ["Set1", "Set2"]:
        bc = d / f"GSE125449_{set_name}_barcodes.tsv.gz"
        gn = d / f"GSE125449_{set_name}_genes.tsv.gz"
        mt = d / f"GSE125449_{set_name}_matrix.mtx.gz"
        if not all(f.exists() for f in [bc, gn, mt]):
            continue
        barcodes = pd.read_csv(bc, header=None, sep="\t")[0].values
        genes = pd.read_csv(gn, header=None, sep="\t")
        gene_names = genes[1].values if genes.shape[1] >= 2 else genes[0].values
        mat = mmread(mt).T.tocsr()
        ad = sc.AnnData(X=mat, obs=pd.DataFrame(index=barcodes),
                        var=pd.DataFrame(index=gene_names))
        ad.var_names_make_unique()
        ad.obs["dataset"] = f"GSE125449_{set_name}"
        ad.obs["study"] = "Ma2019"
        datasets.append(ad)
        print(f"  {set_name}: {ad.shape[0]:,} cells × {ad.shape[1]:,} genes")

    if len(datasets) == 0:
        return None
    adata = sc.concat(datasets, join="outer", fill_value=0)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    return adata


def load_gse140228():
    """Load Zhang et al. 2019 — HCC CD45+ immune landscape."""
    print("[GSE140228] Loading Zhang 2019 (Cell)...")
    d = DATA / "pre_treatment" / "GSE140228_Zhang2019_immune"
    mt = d / "GSE140228_UMI_counts_Droplet.mtx.gz"
    bc = d / "GSE140228_UMI_counts_Droplet_barcodes.tsv.gz"
    gn = d / "GSE140228_UMI_counts_Droplet_genes.tsv.gz"
    ci = d / "GSE140228_UMI_counts_Droplet_cellinfo.tsv.gz"

    if not mt.exists():
        print("  [SKIP] Not downloaded yet")
        return None

    barcodes = pd.read_csv(bc, header=None, sep="\t")[0].values
    # Genes file has header: ENSEMBL SYMBOL CHR START END STRAND BIOTYPE
    genes = pd.read_csv(gn, sep="\t")
    gene_names = genes["SYMBOL"].values
    mat = mmread(mt).T.tocsr()

    adata = sc.AnnData(X=mat, obs=pd.DataFrame(index=barcodes),
                       var=pd.DataFrame(index=gene_names))
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # Cell info (index = Barcode column matching barcodes)
    if ci.exists():
        cell_info = pd.read_csv(ci, sep="\t", index_col=0)
        common = adata.obs_names.intersection(cell_info.index)
        print(f"  Cellinfo overlap: {len(common):,} / {len(adata):,} cells")
        if len(common) > 0:
            for col in ["Donor", "Tissue", "celltype_sub", "celltype_global", "Histology"]:
                if col in cell_info.columns:
                    adata.obs[col] = cell_info.reindex(adata.obs_names)[col].values

    adata.obs["dataset"] = "GSE140228"
    adata.obs["study"] = "Zhang2019"
    print(f"  {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    return adata


def load_gse149614():
    """Load Lu et al. 2022 — HCC multicellular ecosystem."""
    print("[GSE149614] Loading Lu 2022 (Nat Commun)...")
    d = DATA / "pre_treatment" / "GSE149614_Lu2022_landscape"
    count_f = d / "GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
    meta_f  = d / "GSE149614_HCC.metadata.updated.txt.gz"

    if not count_f.exists():
        print("  [SKIP] Not downloaded yet")
        return None

    # This is a dense txt matrix (gene × cell), genes in index auto-detected
    print("  Loading full count matrix (may take a while)...")
    counts = pd.read_csv(count_f, sep="\t")
    gene_names = counts.index.values
    cell_barcodes = counts.columns.values
    
    adata = sc.AnnData(X=sparse.csr_matrix(counts.values.astype(np.float32).T),
                       obs=pd.DataFrame(index=cell_barcodes),
                       var=pd.DataFrame(index=gene_names))
    adata.var_names_make_unique()

    # Metadata
    if meta_f.exists():
        meta = pd.read_csv(meta_f, sep="\t")
        meta = meta.set_index("Cell")
        common = adata.obs_names.intersection(meta.index)
        print(f"  Metadata overlap: {len(common):,} / {len(adata):,} cells")
        if len(common) > 0:
            meta = meta.reindex(adata.obs_names)
            for col in ["site", "patient", "celltype", "stage", "virus"]:
                if col in meta.columns:
                    adata.obs[col] = meta[col].values

    adata.obs["dataset"] = "GSE149614"
    adata.obs["study"] = "Lu2022"
    print(f"  {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    return adata


# ═══════════════════════════════════════════════════════════════════════
# 1e. LOAD E-MTAB-10553
# ═══════════════════════════════════════════════════════════════════════

def load_emtab10553():
    """Load Pfister et al. 2021 — NASH-HCC + anti-PD1 scRNA-seq."""
    print("[E-MTAB-10553] Loading Pfister 2021...")
    d = PROC / "E-MTAB-10553"
    mt = d / "matrix.mtx.gz"
    bc = d / "barcodes.tsv.gz"
    gn = d / "genes.tsv.gz"
    meta_f = d / "metadata.csv"

    if not mt.exists():
        print("  [SKIP] Not converted yet — run 00b_convert_emtab10553.R first")
        return None

    barcodes = pd.read_csv(bc, header=None, sep="\t")[0].values
    genes = pd.read_csv(gn, header=None, sep="\t")
    gene_names = genes[1].values if genes.shape[1] >= 2 else genes[0].values
    mat = mmread(mt).T.tocsr()

    adata = sc.AnnData(X=mat, obs=pd.DataFrame(index=barcodes),
                       var=pd.DataFrame(index=gene_names))
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    if meta_f.exists():
        meta = pd.read_csv(meta_f, index_col=0)
        common = adata.obs_names.intersection(meta.index)
        if len(common) > 0:
            meta = meta.reindex(adata.obs_names)
            for col in meta.columns:
                adata.obs[col] = meta[col].values

    adata.obs["dataset"] = "E-MTAB-10553"
    adata.obs["study"] = "Pfister2021"
    print(f"  {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    return adata


# ═══════════════════════════════════════════════════════════════════════
# 2.  QC FILTERING
# ═══════════════════════════════════════════════════════════════════════

def qc_filter(adata, label=""):
    """Standard scRNA-seq QC."""
    n_before = adata.shape[0]
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)

    # Adaptive thresholds
    min_genes = 200
    max_genes = np.percentile(adata.obs["n_genes_by_counts"], 99)
    max_mt = 25

    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs["n_genes_by_counts"] < max_genes].copy()
    adata = adata[adata.obs["pct_counts_mt"] < max_mt].copy()

    print(f"  [{label}] QC: {n_before:,} → {adata.shape[0]:,} cells "
          f"(removed {n_before - adata.shape[0]:,})")
    return adata


# ═══════════════════════════════════════════════════════════════════════
# 3.  MAIN PIPELINE
# ═══════════════════════════════════════════════════════════════════════

def main_orig():
    print("=" * 70)
    print("HCC IO Response — Preprocessing & Integration")
    print("=" * 70)

    # ── Load all datasets ──────────────────────────────────────────
    datasets = {}

    # Primary IO dataset (must have this)
    adata_magen = load_gse206325()
    adata_magen = qc_filter(adata_magen, "GSE206325")
    datasets["GSE206325"] = adata_magen

    # Reference datasets
    adata_ma = load_gse125449()
    if adata_ma is not None:
        adata_ma = qc_filter(adata_ma, "GSE125449")
        datasets["GSE125449"] = adata_ma

    adata_zhang = load_gse140228()
    if adata_zhang is not None:
        adata_zhang = qc_filter(adata_zhang, "GSE140228")
        datasets["GSE140228"] = adata_zhang

    adata_lu = load_gse149614()
    if adata_lu is not None:
        adata_lu = qc_filter(adata_lu, "GSE149614")
        datasets["GSE149614"] = adata_lu

    adata_pfister = load_emtab10553()
    if adata_pfister is not None:
        adata_pfister = qc_filter(adata_pfister, "E-MTAB-10553")
        datasets["E-MTAB-10553"] = adata_pfister

    # ── Find shared genes & concatenate ────────────────────────────
    print("\n[MERGE] Finding shared genes...")
    gene_sets = [set(ad.var_names) for ad in datasets.values()]
    shared_genes = gene_sets[0]
    for gs in gene_sets[1:]:
        shared_genes = shared_genes.intersection(gs)
    shared_genes = sorted(shared_genes)
    print(f"  Shared genes: {len(shared_genes)}")

    subsets = []
    for name, ad in datasets.items():
        sub = ad[:, shared_genes].copy()
        subsets.append(sub)
        print(f"  {name}: {sub.shape[0]:,} cells")

    adata = sc.concat(subsets, join="outer")
    adata.obs_names_make_unique()
    print(f"\n  Combined: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

    # ── Normalize ──────────────────────────────────────────────────
    print("\n[NORMALIZE] ...")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    # ── HVGs ───────────────────────────────────────────────────────
    print("[HVG] Selecting highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="dataset")
    n_hvg = adata.var["highly_variable"].sum()
    print(f"  HVGs: {n_hvg}")

    # ── PCA ────────────────────────────────────────────────────────
    print("[PCA] ...")
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True)

    # ── Harmony integration ────────────────────────────────────────
    print("[HARMONY] Batch integration by dataset...")
    import harmonypy
    ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, "dataset",
                               max_iter_harmony=30)
    # harmonypy >= 0.2: Z_corr is already (cells, PCs)
    Z = ho.Z_corr if ho.Z_corr.shape[0] == adata.n_obs else ho.Z_corr.T
    adata.obsm["X_pca_harmony"] = Z

    # ── Neighbors & UMAP ──────────────────────────────────────────
    print("[UMAP] ...")
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_pcs=30)
    sc.tl.umap(adata, min_dist=0.3)

    # ── Leiden clustering ──────────────────────────────────────────
    print("[LEIDEN] ...")
    sc.tl.leiden(adata, resolution=1.0, key_added="leiden")
    print(f"  Clusters: {adata.obs['leiden'].nunique()}")

    # ── QC summary table ──────────────────────────────────────────
    qc_summary = adata.obs.groupby("dataset").agg(
        n_cells=("n_genes_by_counts", "count"),
        median_genes=("n_genes_by_counts", "median"),
        median_counts=("total_counts", "median"),
        median_mito=("pct_counts_mt", "median"),
    ).round(1)
    qc_summary.to_csv(TAB / "qc_summary.csv")
    print("\n  QC Summary:")
    print(qc_summary.to_string())

    # ── Composition table ─────────────────────────────────────────
    if "response" in adata.obs.columns:
        resp_counts = adata.obs.groupby(["dataset", "response"]).size().unstack(fill_value=0)
        resp_counts.to_csv(TAB / "response_counts_by_dataset.csv")
        print("\n  Response counts:")
        print(resp_counts.to_string())

    # ═══════════════════════════════════════════════════════════════
    # FIGURES
    # ═══════════════════════════════════════════════════════════════
    print("\n[FIGURES] Generating QC & integration plots...")

    # ── Recover GSE206325 metadata if missing from concat ──────────
    if "response" not in adata.obs.columns or adata.obs["response"].isna().all():
        print("[META] Recovering GSE206325 metadata for integration plots...")
        cell_meta = pd.read_csv(PROC / "GSE206325_metadata.csv", index_col=0)
        sample_ann = pd.read_csv(PROC / "GSE206325_sample_annotations.csv", index_col=0)
        sample_ann.index = sample_ann.index.astype(str)

        shared = adata.obs_names.intersection(cell_meta.index)
        print(f"  Matched {len(shared):,} GSE206325 barcodes in atlas")

        if "sample_id" not in adata.obs.columns:
            adata.obs["sample_id"] = np.nan
        adata.obs.loc[shared, "sample_id"] = cell_meta.loc[shared, "cell_to_sample_ID"].values.astype(str)

        for col in ["patient_ID", "tissue", "treatment_Resp"]:
            mapping = sample_ann[col].to_dict()
            adata.obs[col] = adata.obs["sample_id"].map(mapping)

        resp_map = {"anti-PD1_R": "Responder", "anti-PD1_NR": "Non-Responder"}
        adata.obs["response"] = adata.obs["treatment_Resp"].map(resp_map)
        print(f"  Response: {adata.obs['response'].value_counts(dropna=False).to_dict()}")
        print(f"  Tissue: {adata.obs['tissue'].value_counts(dropna=False).to_dict()}")

    # ─── Figure 1: QC violin plots ────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, metric, title in zip(axes,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            ["Genes / Cell", "UMI Counts / Cell", "% Mitochondrial"]):
        data_list = []
        labels = []
        for ds in sorted(adata.obs["dataset"].unique()):
            vals = adata.obs.loc[adata.obs["dataset"] == ds, metric].values
            data_list.append(vals)
            labels.append(ds)
        vp = ax.violinplot(data_list, showmedians=True)
        for body in vp["bodies"]:
            body.set_alpha(0.7)
        ax.set_xticks(range(1, len(labels)+1))
        ax.set_xticklabels([l.replace("GSE", "").replace("_", "\n") for l in labels], fontsize=8)
        ax.set_title(title, fontweight="bold")
        ax.set_ylabel(title)
    plt.suptitle(f"QC Metrics — {adata.shape[0]:,} cells, {len(datasets)} datasets", fontweight="bold")
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"01_qc_violins.{ext}", dpi=200, bbox_inches="tight")
    plt.close()

    # ─── Figure 2: Integration UMAP (4-panel) ─────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(14, 13))

    # Panel A: UMAP by dataset
    ax = axes[0, 0]
    ds_list = sorted(adata.obs["dataset"].unique())
    cmap_ds = plt.cm.Set2(np.linspace(0, 1, len(ds_list)))
    for i, ds in enumerate(ds_list):
        mask = adata.obs["dataset"] == ds
        ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                   s=0.5, alpha=0.2, c=[cmap_ds[i]], label=ds, rasterized=True)
    ax.legend(fontsize=7, markerscale=10)
    ax.set_title("A. Dataset", fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    # Panel B: UMAP by study
    ax = axes[0, 1]
    studies = sorted(adata.obs["study"].dropna().unique())
    cmap_st = plt.cm.tab10(np.linspace(0, 1, len(studies)))
    for i, st in enumerate(studies):
        mask = adata.obs["study"] == st
        ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                   s=0.5, alpha=0.2, c=[cmap_st[i]], label=st, rasterized=True)
    ax.legend(fontsize=8, markerscale=10)
    ax.set_title("B. Study", fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    # Panel C: UMAP by tissue (GSE206325 has tumor/normal/PBMC/LN)
    ax = axes[1, 0]
    if "tissue" in adata.obs.columns:
        tissues = sorted(adata.obs["tissue"].dropna().unique())
        cmap_ti = plt.cm.Paired(np.linspace(0, 1, max(len(tissues), 2)))
        for i, ti in enumerate(tissues):
            mask = adata.obs["tissue"] == ti
            ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                       s=0.5, alpha=0.2, c=[cmap_ti[i]], label=ti, rasterized=True)
        ax.legend(fontsize=8, markerscale=10)
    ax.set_title("C. Tissue", fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    # Panel D: UMAP by response (R vs NR)
    ax = axes[1, 1]
    if "response" in adata.obs.columns:
        resp_colors = {"Responder": "#2ecc71", "Non-Responder": "#e74c3c"}
        for resp, color in resp_colors.items():
            mask = adata.obs["response"] == resp
            if mask.sum() > 0:
                ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                           s=0.5, alpha=0.2, c=color, label=f"{resp} ({mask.sum():,})",
                           rasterized=True)
        # Also plot cells without response labels (other datasets)
        mask_na = adata.obs["response"].isna()
        if mask_na.sum() > 0:
            ax.scatter(adata.obsm["X_umap"][mask_na, 0], adata.obsm["X_umap"][mask_na, 1],
                       s=0.3, alpha=0.05, c="gray", label=f"No label ({mask_na.sum():,})",
                       rasterized=True)
        ax.legend(fontsize=8, markerscale=10)
    ax.set_title("D. IO Response (GSE206325)", fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    plt.suptitle(f"Integrated HCC IO Atlas — {adata.shape[0]:,} cells\n"
                 f"Harmony batch correction, {len(datasets)} datasets",
                 fontweight="bold", fontsize=13)
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"01_integration_umap.{ext}", dpi=200, bbox_inches="tight")
    plt.close()

    # ─── Figure 3: Leiden clusters ────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 8))
    n_clusters = adata.obs["leiden"].nunique()
    cmap_cl = plt.cm.tab20(np.linspace(0, 1, min(n_clusters, 20)))
    for i, cl in enumerate(sorted(adata.obs["leiden"].unique(), key=int)):
        mask = adata.obs["leiden"] == cl
        ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                   s=0.5, alpha=0.3, c=[cmap_cl[i % 20]], rasterized=True)
        # Label centroid
        cx, cy = adata.obsm["X_umap"][mask].mean(axis=0)
        ax.text(cx, cy, cl, fontsize=6, ha="center", va="center",
                fontweight="bold", bbox=dict(boxstyle="round,pad=0.1",
                fc="white", ec="gray", alpha=0.7))
    ax.set_title(f"Leiden Clustering (n={n_clusters})", fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"01_leiden_clusters.{ext}", dpi=200, bbox_inches="tight")
    plt.close()

    # ── Save ──────────────────────────────────────────────────────
    out_path = PROC / "integrated_io_atlas.h5ad"
    print(f"\n[SAVE] Writing {out_path}...")
    adata.write(out_path)

    print(f"\n{'='*70}")
    print(f"[DONE] Integrated IO Atlas:")
    print(f"  Cells:    {adata.shape[0]:,}")
    print(f"  Genes:    {adata.shape[1]:,}")
    print(f"  Datasets: {adata.obs['dataset'].nunique()}")
    print(f"  Clusters: {adata.obs['leiden'].nunique()}")
    if "response" in adata.obs.columns:
        print(f"  Response: {adata.obs['response'].value_counts().to_dict()}")
    print(f"  Output:   {out_path}")
    print(f"  Figures:  {FIG}")
    print(f"={'='*70}")


if __name__ == "__main__":
    main()

def test_plot():
    import anndata as ad
    import matplotlib.pyplot as plt
    import numpy as np
    from pathlib import Path

    print("Loading H5AD...")
    adata = ad.read_h5ad("data/processed/integrated_io_atlas.h5ad", backed="r")
    print(f"Response: {adata.obs['response'].value_counts(dropna=False).to_dict()}")

    fig, axes = plt.subplots(2, 2, figsize=(14, 13))

    ax = axes[1, 1]
    if "response" in adata.obs.columns:
        print("Plotting response...")
        resp_colors = {"Responder": "#2ecc71", "Non-Responder": "#e74c3c"}
        for resp, color in resp_colors.items():
            mask = adata.obs["response"] == resp
            print(f"Mask {resp}: {mask.sum()}")
            if mask.sum() > 0:
                ax.scatter(adata.obsm["X_umap"][mask.values, 0], adata.obsm["X_umap"][mask.values, 1],
                           s=0.5, alpha=0.2, c=color, label=f"{resp} ({mask.sum():,})",
                           rasterized=True)
        ax.legend(fontsize=8, markerscale=10)
    ax.set_title("D. IO Response (GSE206325)", fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    
    print(f"Xlim: {ax.get_xlim()}, Ylim: {ax.get_ylim()}")
    plt.savefig("test_out.png", dpi=200, bbox_inches="tight")
    print("Done!")

if __name__ == "__main__":
    test_plot()
