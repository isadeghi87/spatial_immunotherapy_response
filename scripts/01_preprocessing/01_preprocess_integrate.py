#!/usr/bin/env python3
"""
01_preprocess_integrate.py — Preprocess & Integrate All HCC IO Datasets

Loads, QCs, normalizes, and integrates:
  - GSE206325  (Magen 2023 Nat Med):  682k cells, anti-PD1 R vs NR, scRNA+TCR
  - GSE125449  (Ma 2019):             10k cells, HCC landscape
  - GSE140228  (Zhang 2019 Cell):     CD45+ immune cells, HCC immune landscape
  - GSE149614  (Lu 2022 Nat Commun):  70k cells, HCC primary/metastatic
  - E-MTAB-10553 (Pfister 2021):      ~150k cells, NASH-HCC + anti-PD1
  - Cappuyns2023 (Cappuyns 2023):     ~87k cells, atezo+bev R vs NR

Pipeline:
  1. Load each dataset → standardize gene names & metadata
  2. QC filtering (genes/cell, mito%)
  3. Per-dataset label harmonization → unified cell_type_raw column
     (GSE206325: full fine-grained labels; GSE140228 / GSE149614: mapped;
      others: "Unknown" → inferred by scANVI)
  4. Optional clinical covariates: sex/age from
     data/clinical_metadata/sex_age_supplement.csv if present, and
     disease_etiology from dataset-specific annotation columns.
  5. Concatenate on shared genes
  6. Normalize, HVG, PCA (context for scVI init)
  7. scVI batch integration (n_latent=30, NB, batch_key=dataset, optional
     categorical_covariate_keys=[disease_etiology] / continuous=[age]).
     Model checkpoint saved to data/processed/scvi_model/.
  8. scANVI semi-supervised label transfer: fine-tunes scVI with cell_type_raw
     labels; predicts cell_type and confidence for ALL 1.06 M cells.
     Model checkpoint saved to data/processed/scanvi_model/.
     Outputs obsm["X_scANVI"], obs["cell_type"], obs["cell_type_confidence"].
  9. UMAP + Leiden on X_scANVI.
 10. Figures + save.

Re-run flags (env vars):
  RETRAIN=1   Force scVI + scANVI retraining even if checkpoints exist.
"""
import warnings
warnings.filterwarnings("ignore")
import argparse

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

# ── scANVI unlabeled sentinel ──────────────────────────────────────────
UNLABELED = "Unknown"

# ─────────────────────────────────────────────────────────────────────
#  LABEL HARMONIZATION MAP
#  GSE140228 celltype_sub → GSE206325-compatible vocabulary
# ─────────────────────────────────────────────────────────────────────
_ZHANG_MAP = {
    # CD4 T
    "CD4-C3-ANXA1":     "CD4_Effector-Memory",
    "CD4-C4-IL7R":      "CD4_Memory",
    "CD4-C5-TCF7":      "CD4_Memory",
    "CD4-C6-CXCL13":    "CD4_Tfh-like",
    "CD4-C7-FOXP3":     "Treg_Treg",
    # Dual / cycling
    "CD4/CD8-C1-CCR7":  "Naive_Naive",
    "CD4/CD8-C2-MKI67": "CD8_Dysfunctional-Proliferating",
    # CD8 T
    "CD8-C3-IL7R":      "CD8_Memory",
    "CD8-C4-CX3CR1":    "CD8_Cytotoxic",
    "CD8-C5-SELL":      "Naive_Naive",
    "CD8-C6-GZMK":      "CD8_Effector",
    "CD8-C7-KLRD1":     "CD8_Cytotoxic",
    "CD8-C8-PDCD1":     "CD8_Dysfunctional-effector",
    "CD8-C9-SLC4A10":   "MAIT",
    # DC
    "DC-C1-CD1C":       "DC_cDC2",
    "DC-C2-FCER1A":     "DC_cDC2",
    "DC-C3-CLEC9A":     "DC_cDC1",
    "DC-C4-LAMP3":      "DC_mregDC",
    "ILCs":             "ILC",
    # B
    "Lymphoid-B":       "B_Naive",
    "Lymphoid-B-Plasma": "B_Plasma",
    # Mast
    "Mast-C1-IL7R":     "Mast",
    "Mast-C2-CPA3":     "Mast",
    # Mono/Mac
    "Mono-C1-CD14":     "Monocyte_CD14",
    "Mono-C2-FCGR3A":   "Monocyte_CD16",
    "Mφ-C1-THBS1":      "Macrophage_SPP1",
    "Mφ-C2-C1QA":       "Macrophage_Kupffer",
    "Mφ-C3-APOE":       "Macrophage_FOLR2",
    "Mφ-C4-GPX3":       "Macrophage_FOLR2",
    "Mφ-C5-VCAN":       "Macrophage_FCN1/CXCL9",
    "Mφ-C6-MARCO":      "Macrophage_Kupffer",
    # NK
    "NK-C1-FCGR3A":     "NK",
    "NK-C2-SELL":       "NK",
    "NK-C3-IFNG":       "NK",
    "NK-C4-IFNG-HSPA1A": "NK",
    "NK-C5-CD69":       "NK",
    "NK-C6-IL7R":       "NK",
    "NK-C7-CD160":      "NK",
    "NK-C8-CD160-HSPA1A": "NK",
    "NK-C9-MKI67":      "NK",
}

# GSE149614 coarse labels → GSE206325 group-level vocabulary
_LU_MAP = {
    "T/NK":        "CD8_Effector",   # majority T/NK → CD8; scANVI refines within T/NK
    "Myeloid":     "Macrophage_SPP1",
    "B":           "B_Naive",
    "Endothelial": "Stromal_Endothelial cells",
    "Fibroblast":  "Stromal_HSC",
    "Hepatocyte":  "Stromal_Hepatocytes",
}


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
    # Unified label for scANVI: drop Doublets/Unknown/Empty clusters
    _exclude = {"Unknown", "Doublets", "Empty", "nan", ""}
    adata.obs["cell_type_raw"] = adata.obs["cell_type_original"].apply(
        lambda x: UNLABELED if (pd.isna(x) or str(x) in _exclude) else str(x)
    )

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

    # Harmonize labels using Zhang→GSE206325 vocabulary map
    if "celltype_sub" in adata.obs.columns:
        adata.obs["cell_type_raw"] = adata.obs["celltype_sub"].map(
            _ZHANG_MAP).fillna(UNLABELED).astype(str)
    else:
        adata.obs["cell_type_raw"] = UNLABELED
    # Keep disease_etiology covariate
    if "Histology" in adata.obs.columns:
        adata.obs["disease_etiology"] = adata.obs["Histology"].fillna("Unknown").astype(str)

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

    # Coarse harmonized labels (fine-grained resolution not available in Lu2022)
    if "celltype" in adata.obs.columns:
        adata.obs["cell_type_raw"] = adata.obs["celltype"].map(
            _LU_MAP).fillna(UNLABELED).astype(str)
    else:
        adata.obs["cell_type_raw"] = UNLABELED
    # disease_etiology from virus column
    if "virus" in adata.obs.columns:
        adata.obs["disease_etiology"] = adata.obs["virus"].fillna("Unknown").astype(str)

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

    adata.obs["cell_type_raw"] = UNLABELED   # no published cell-type labels
    adata.obs["dataset"] = "E-MTAB-10553"
    adata.obs["study"] = "Pfister2021"
    print(f"  {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    return adata


# ═══════════════════════════════════════════════════════════════════════
# 1f. LOAD CAPPUYNS 2023
# ═══════════════════════════════════════════════════════════════════════

def load_cappuyns():
    """Load Cappuyns et al. 2023 — atezo+bev R vs NR (32 patients)."""
    print("[Cappuyns2023] Loading Cappuyns 2023 (Nat Commun)...")
    d = DATA / "pre_treatment" / "Cappuyns2023_NatCommun_atezobev" / "HCC_data_VIB_CCB_portal"
    mtx_dir = d / "processed_feature_bc_matrix"
    meta_f  = d / "meta_data.txt"

    if not mtx_dir.exists():
        print("  [SKIP] Not downloaded yet")
        return None

    # Load 10x-formatted MTX directory
    mtx_files = list(mtx_dir.glob("*.mtx*"))
    if not mtx_files:
        print("  [SKIP] No MTX files found")
        return None

    # features.tsv.gz is single-column (gene names only, no Ensembl IDs).
    # sc.read_10x_mtx always tries genes[1] for a "gene_symbols" column → KeyError.
    # Build AnnData manually instead.
    import gzip, scipy.io
    with gzip.open(mtx_dir / "features.tsv.gz", "rt") as fh:
        gene_names = pd.read_csv(fh, sep="\t", header=None)[0].tolist()
    with gzip.open(mtx_dir / "barcodes.tsv.gz", "rt") as fh:
        barcodes = pd.read_csv(fh, sep="\t", header=None)[0].tolist()
    with gzip.open(mtx_dir / "matrix.mtx.gz", "rt") as fh:
        X_mat = scipy.io.mmread(fh).T.tocsr()   # 10x MTX is genes×cells; transpose to cells×genes
    import anndata as _ad
    adata = _ad.AnnData(
        X=X_mat,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=gene_names),
    )
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    print(f"    {mtx_dir.name}: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

    # Metadata
    if meta_f.exists():
        meta = pd.read_csv(meta_f, sep="\t", index_col=0)
        common = adata.obs_names.intersection(meta.index)
        print(f"  Metadata overlap: {len(common):,} / {len(adata):,} cells")
        if len(common) > 0:
            meta = meta.reindex(adata.obs_names)
            adata.obs["patient_ID"]    = meta["Patient_ID"].values.astype(str)
            adata.obs["Tumour_ID"]     = meta["Tumour_ID"].values.astype(str)
            # Normalize response to shared vocabulary (Responder / Non-Responder)
            _resp_map = {
                "Responder":           "Responder",
                "NonResponder":        "Non-Responder",
                "Non-Responder":       "Non-Responder",
                "DeathBeforeImaging":  "Non-Responder",   # treated as non-response for IO analysis
            }
            adata.obs["response"] = meta["Response"].map(_resp_map).fillna("Unknown").astype(str)
            adata.obs["treatment"]     = meta["Treatment"].values.astype(str)
            print(f"  Response: {adata.obs['response'].value_counts().to_dict()}")

    adata.obs["cell_type_raw"]    = UNLABELED   # no published cell-type labels
    adata.obs["dataset"]          = "Cappuyns2023"
    adata.obs["study"]            = "Cappuyns2023"
    adata.obs["disease_etiology"] = "HCC"
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

def main():
    import os, torch
    import scvi

    retrain = os.environ.get("RETRAIN", "0") == "1"
    scvi.settings.seed = 0
    use_gpu = torch.cuda.is_available()
    accelerator = "gpu" if use_gpu else "cpu"
    print(f"  accelerator={accelerator}")

    print("=" * 70)
    print("HCC IO Response — Preprocessing & Integration")
    print("=" * 70)

    # ── Load all datasets ──────────────────────────────────────────
    datasets = {}

    adata_magen = load_gse206325()
    adata_magen = qc_filter(adata_magen, "GSE206325")
    datasets["GSE206325"] = adata_magen

    adata_ma = load_gse125449()
    if adata_ma is not None:
        adata_ma.obs["cell_type_raw"] = UNLABELED
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

    adata_capp = load_cappuyns()
    if adata_capp is not None:
        adata_capp = qc_filter(adata_capp, "Cappuyns2023")
        datasets["Cappuyns2023"] = adata_capp

    # ── Clinical covariates: sex / age supplement (optional) ───────
    sex_age_f = DATA / "clinical_metadata" / "sex_age_supplement.csv"
    sex_age = None
    if sex_age_f.exists():
        sex_age = pd.read_csv(sex_age_f, index_col=0)  # index = patient_ID
        print(f"[COVARIATES] Loaded sex/age for {len(sex_age)} patients")
    else:
        print("[COVARIATES] sex_age_supplement.csv not found — skipping sex/age covariates")

    # ── Find shared genes & concatenate ────────────────────────────
    print("\n[MERGE] Finding shared genes...")
    gene_sets = [set(ad.var_names) for ad in datasets.values()]
    shared_genes = sorted(gene_sets[0].intersection(*gene_sets[1:]))
    print(f"  Shared genes: {len(shared_genes)}")

    subsets = []
    for name, ad in datasets.items():
        sub = ad[:, shared_genes].copy()
        # Ensure cell_type_raw column exists on every subset
        if "cell_type_raw" not in sub.obs.columns:
            sub.obs["cell_type_raw"] = UNLABELED
        subsets.append(sub)
        print(f"  {name}: {sub.shape[0]:,} cells")

    adata = sc.concat(subsets, join="outer")   # outer: keep all obs columns (response, patient_ID, tissue…)
    # Trim vars back to shared_genes (outer join may add NaN-filled genes if datasets diverged)
    adata = adata[:, shared_genes].copy()
    adata.obs_names_make_unique()
    print(f"\n  Combined: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

    # Attach sex/age to cells if supplement present
    if sex_age is not None and "patient_ID" in adata.obs.columns:
        if "sex" in sex_age.columns:
            adata.obs["sex"] = adata.obs["patient_ID"].map(
                sex_age["sex"].to_dict()).fillna("Unknown").astype(str)
        if "age" in sex_age.columns:
            adata.obs["age"] = pd.to_numeric(
                adata.obs["patient_ID"].map(sex_age["age"].to_dict()), errors="coerce")

    # ── Normalize ──────────────────────────────────────────────────
    print("\n[NORMALIZE] ...")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    # ── HVGs ───────────────────────────────────────────────────────
    print("[HVG] Selecting highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="dataset")
    print(f"  HVGs: {adata.var['highly_variable'].sum()}")

    # ── PCA (context; scVI will build its own latent space) ────────
    print("[PCA] ...")
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True)

    # ── Determine covariates to pass to scVI ───────────────────────
    cat_covariates = []
    cont_covariates = []
    if "disease_etiology" in adata.obs.columns:
        n_valid = (adata.obs["disease_etiology"].notna() &
                   (adata.obs["disease_etiology"] != "Unknown")).sum()
        if n_valid > 0:
            adata.obs["disease_etiology"] = (
                adata.obs["disease_etiology"].fillna("Unknown").astype(str))
            cat_covariates.append("disease_etiology")
            print(f"  Covariate (cat):  disease_etiology ({n_valid:,} annotated cells)")
    if "age" in adata.obs.columns:
        n_age = adata.obs["age"].notna().sum()
        if n_age / len(adata) > 0.5:   # only use if >50% non-NaN
            adata.obs["age"] = adata.obs["age"].fillna(adata.obs["age"].median())
            cont_covariates.append("age")
            print(f"  Covariate (cont): age ({n_age:,} annotated cells)")

    # ── scVI batch integration ─────────────────────────────────────
    scvi_dir = PROC / "scvi_model"
    if scvi_dir.exists() and not retrain:
        print(f"[scVI] Loading checkpoint from {scvi_dir}...")
        model = scvi.model.SCVI.load(str(scvi_dir), adata)
    else:
        print("[scVI] Setting up and training scVI model...")
        scvi.model.SCVI.setup_anndata(
            adata,
            layer="counts",
            batch_key="dataset",
            categorical_covariate_keys=cat_covariates or None,
            continuous_covariate_keys=cont_covariates or None,
        )
        model = scvi.model.SCVI(
            adata,
            n_latent=30,
            n_layers=2,
            gene_likelihood="nb",
        )
        max_epochs = max(60, min(200, int(2_000_000 / max(adata.n_obs, 1))))
        print(f"  max_epochs={max_epochs}  accelerator={accelerator}")
        model.train(max_epochs=max_epochs, accelerator=accelerator)
        model.save(str(scvi_dir), overwrite=True)
        print(f"  Saved scVI model → {scvi_dir}")

    adata.obsm["X_scVI"] = model.get_latent_representation()

    # ── scANVI semi-supervised label transfer ─────────────────────
    # cell_type_raw = known labels for GSE206325/GSE140228/GSE149614;
    # UNLABELED ("Unknown") for all other cells — scANVI infers those.
    print("\n[scANVI] Semi-supervised cell-type label transfer...")
    n_labeled = (adata.obs["cell_type_raw"] != UNLABELED).sum()
    label_frac = n_labeled / len(adata)
    print(f"  Labeled cells:   {n_labeled:,} ({label_frac:.1%})")
    print(f"  Unlabeled cells: {len(adata) - n_labeled:,}")

    scanvi_dir = PROC / "scanvi_model"
    if scanvi_dir.exists() and not retrain:
        print(f"  Loading checkpoint from {scanvi_dir}...")
        scanvi_model = scvi.model.SCANVI.load(str(scanvi_dir), adata)
    else:
        scanvi_model = scvi.model.SCANVI.from_scvi_model(
            model,
            unlabeled_category=UNLABELED,
            labels_key="cell_type_raw",
        )
        # Lower LR (default ~1e-3 can cause NaN divergence) + early stopping
        scanvi_model.train(
            max_epochs=20,
            accelerator=accelerator,
            plan_kwargs={"lr": 5e-4},
            early_stopping=True,
        )
        scanvi_model.save(str(scanvi_dir), overwrite=True)
        print(f"  Saved scANVI model → {scanvi_dir}")

    adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()
    adata.obs["cell_type"] = scanvi_model.predict()        # hard label, all cells
    soft = scanvi_model.predict(soft=True)                 # DataFrame: cells × labels
    adata.obs["cell_type_confidence"] = soft.max(axis=1).values

    print(f"  cell_type distribution (top 15):")
    print(adata.obs["cell_type"].value_counts().head(15).to_string())
    print(f"  Confidence: median={adata.obs['cell_type_confidence'].median():.3f}  "
          f"mean={adata.obs['cell_type_confidence'].mean():.3f}")

    # Save cell-type table for audit
    ct_counts = adata.obs[["dataset", "cell_type"]].value_counts().unstack(fill_value=0)
    ct_counts.to_csv(TAB / "cell_type_counts_by_dataset.csv")

    # ── Neighbors & UMAP on X_scANVI ─────────────────────────────
    print("\n[UMAP] Building graph on X_scANVI...")
    sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=15)
    sc.tl.umap(adata, min_dist=0.3)

    # ── Leiden clustering ─────────────────────────────────────────
    print("[LEIDEN] ...")
    sc.tl.leiden(adata, resolution=1.0, key_added="leiden")
    print(f"  Clusters: {adata.obs['leiden'].nunique()}")

    # ── QC summary table ─────────────────────────────────────────
    qc_summary = adata.obs.groupby("dataset").agg(
        n_cells=("n_genes_by_counts", "count"),
        median_genes=("n_genes_by_counts", "median"),
        median_counts=("total_counts", "median"),
        median_mito=("pct_counts_mt", "median"),
    ).round(1)
    qc_summary.to_csv(TAB / "qc_summary.csv")
    print("\n  QC Summary:")
    print(qc_summary.to_string())

    if "response" in adata.obs.columns:
        resp_counts = (
            adata.obs.groupby(["dataset", "response"], dropna=False)
            .size().unstack(fill_value=0)
        )
        resp_counts.to_csv(TAB / "response_counts_by_dataset.csv")
        print("\n  Response counts:")
        print(resp_counts.to_string())

    # ═══════════════════════════════════════════════════════════════
    # FIGURES
    # ═══════════════════════════════════════════════════════════════
    print("\n[FIGURES] Generating QC & integration plots...")

    # ── Fig 1: QC violins ─────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, metric, title in zip(
            axes,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            ["Genes / Cell", "UMI Counts / Cell", "% Mitochondrial"]):
        data_list, labels = [], []
        for ds in sorted(adata.obs["dataset"].unique()):
            data_list.append(adata.obs.loc[adata.obs["dataset"] == ds, metric].values)
            labels.append(ds)
        vp = ax.violinplot(data_list, showmedians=True)
        for body in vp["bodies"]:
            body.set_alpha(0.7)
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(
            [l.replace("GSE", "").replace("_", "\n") for l in labels], fontsize=8)
        ax.set_title(title, fontweight="bold")
        ax.set_ylabel(title)
    plt.suptitle(
        f"QC Metrics — {adata.shape[0]:,} cells, {len(datasets)} datasets",
        fontweight="bold")
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"01_qc_violins.{ext}", dpi=200, bbox_inches="tight")
    plt.close()

    # ── Fig 2: Integration UMAP (4-panel) ─────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(14, 13))

    ax = axes[0, 0]
    ds_list = sorted(adata.obs["dataset"].unique())
    cmap_ds = plt.cm.Set2(np.linspace(0, 1, max(len(ds_list), 3)))
    for i, ds in enumerate(ds_list):
        mask = adata.obs["dataset"] == ds
        ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                   s=0.5, alpha=0.2, c=[cmap_ds[i]], label=ds, rasterized=True)
    ax.legend(fontsize=7, markerscale=10)
    ax.set_title("A. Dataset (scVI batch correction)", fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    ax = axes[0, 1]
    if "tissue" in adata.obs.columns:
        tissues = sorted(adata.obs["tissue"].dropna().unique())
        cmap_ti = plt.cm.Paired(np.linspace(0, 1, max(len(tissues), 2)))
        for i, ti in enumerate(tissues):
            mask = adata.obs["tissue"] == ti
            ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                       s=0.5, alpha=0.2, c=[cmap_ti[i]], label=ti, rasterized=True)
        ax.legend(fontsize=8, markerscale=10)
    ax.set_title("B. Tissue", fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    ax = axes[1, 0]
    if "response" in adata.obs.columns:
        resp_colors = {"Responder": "#2ecc71", "Non-Responder": "#e74c3c"}
        for resp, color in resp_colors.items():
            mask = adata.obs["response"] == resp
            if mask.sum() > 0:
                ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                           s=0.5, alpha=0.2, c=color,
                           label=f"{resp} ({mask.sum():,})", rasterized=True)
        mask_na = ~adata.obs["response"].isin(["Responder", "Non-Responder"])
        if mask_na.sum() > 0:
            ax.scatter(adata.obsm["X_umap"][mask_na, 0], adata.obsm["X_umap"][mask_na, 1],
                       s=0.3, alpha=0.05, c="gray",
                       label=f"No label ({mask_na.sum():,})", rasterized=True)
        ax.legend(fontsize=8, markerscale=10)
    ax.set_title("C. IO Response", fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    ax = axes[1, 1]
    n_clusters = adata.obs["leiden"].nunique()
    cmap_cl = plt.cm.tab20(np.linspace(0, 1, min(n_clusters, 20)))
    for i, cl in enumerate(sorted(adata.obs["leiden"].unique(), key=int)):
        mask = adata.obs["leiden"] == cl
        ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                   s=0.5, alpha=0.3, c=[cmap_cl[i % 20]], rasterized=True)
        cx, cy = adata.obsm["X_umap"][mask].mean(axis=0)
        ax.text(cx, cy, cl, fontsize=6, ha="center", va="center", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="gray", alpha=0.7))
    ax.set_title(f"D. Leiden (n={n_clusters})", fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    plt.suptitle(
        f"scVI/scANVI Atlas — {adata.shape[0]:,} cells, {len(datasets)} datasets",
        fontweight="bold", fontsize=13)
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"01_integration_umap.{ext}", dpi=200, bbox_inches="tight")
    plt.close()

    # ── Fig 3: scANVI cell-type annotation (all cells) ─────────────
    cts = [ct for ct in adata.obs["cell_type"].value_counts().index
           if ct not in {UNLABELED, "nan", ""}]
    n_ct = len(cts)
    ct_colors = plt.cm.tab20(np.linspace(0, 1, min(n_ct, 20)))
    ct_color_map = {ct: ct_colors[i % 20] for i, ct in enumerate(cts)}

    fig, axes = plt.subplots(1, 2, figsize=(22, 9))
    ax = axes[0]
    ax.scatter(
        adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1],
        s=0.2, alpha=0.1, c="lightgray", rasterized=True, zorder=0)
    for ct in cts:
        mask = adata.obs["cell_type"] == ct
        ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                   s=0.4, alpha=0.3, c=[ct_color_map[ct]], rasterized=True,
                   label=f"{ct} ({mask.sum():,})")
    ax.legend(fontsize=7, markerscale=8, ncol=2,
              bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.set_title("Cell-Type Annotations (scANVI)", fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    ax = axes[1]
    top_cts = adata.obs["cell_type"].value_counts().head(20)
    bar_colors = [ct_color_map.get(ct, "gray") for ct in top_cts.index]
    ax.barh(range(len(top_cts)), top_cts.values, color=bar_colors)
    ax.set_yticks(range(len(top_cts)))
    ax.set_yticklabels(top_cts.index, fontsize=9)
    ax.invert_yaxis()
    for j, v in enumerate(top_cts.values):
        ax.text(v + top_cts.values.max() * 0.01, j, f"{v:,}", va="center", fontsize=8)
    ax.set_title("Cell-Type Abundance (Top 20)", fontweight="bold")
    ax.set_xlabel("Number of Cells")

    plt.suptitle(
        f"HCC IO Atlas Cell-Type Map — {adata.shape[0]:,} cells",
        fontweight="bold", fontsize=13)
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"01_celltype_umap.{ext}", dpi=200, bbox_inches="tight")
    plt.close()

    # ── Fig 4: scANVI confidence map ──────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 8))
    sc_conf = ax.scatter(
        adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1],
        c=adata.obs["cell_type_confidence"].values,
        s=0.3, alpha=0.4, cmap="viridis", vmin=0, vmax=1, rasterized=True)
    plt.colorbar(sc_conf, ax=ax, label="scANVI prediction confidence")
    ax.set_title("scANVI Label Transfer Confidence", fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"01_scanvi_confidence.{ext}", dpi=200, bbox_inches="tight")
    plt.close()

    # ── Save ──────────────────────────────────────────────────────
    out_path = PROC / "integrated_io_atlas.h5ad"
    print(f"\n[SAVE] Writing {out_path}...")
    # Sanitize obs: h5py requires pure str, not mixed str/float (NaN) columns
    for _col in adata.obs.select_dtypes(include="object").columns:
        adata.obs[_col] = adata.obs[_col].astype(str)
    adata.write(out_path)

    print(f"\n{'='*70}")
    print(f"[DONE] Integrated IO Atlas:")
    print(f"  Cells:    {adata.shape[0]:,}")
    print(f"  Genes:    {adata.shape[1]:,}")
    print(f"  Datasets: {adata.obs['dataset'].nunique()}")
    print(f"  Clusters: {adata.obs['leiden'].nunique()}")
    print(f"  Cell types (scANVI): {adata.obs['cell_type'].nunique()}")
    if "response" in adata.obs.columns:
        resp = adata.obs["response"].value_counts().to_dict()
        print(f"  Response: {resp}")
    print(f"  scVI model:   {scvi_dir}")
    print(f"  scANVI model: {scanvi_dir}")
    print(f"  Output:   {out_path}")
    print(f"  Figures:  {FIG}")
    print(f"={'='*70}")


if __name__ == "__main__":
    main()
