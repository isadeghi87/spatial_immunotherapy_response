#!/usr/bin/env python3
"""
Spatial Immunotherapy Response — Public Data Analysis
Uses public scRNA-seq datasets from HCC anti-PD1 clinical trials:
  - GSE125449 (Ma et al., 2019): HCC scRNA-seq pre/post anti-PD1
  - GSE234241 (Xue et al., 2023): HCC combination immunotherapy
  - GSE206325 (Liu et al., 2022): HCC anti-PD1 treated
  - GSE202069 (Zhu et al., 2023): HCC anti-PDL1 treated (bulk)

Analyses:
  1. Load & preprocess public scRNA-seq datasets
  2. Score IO signatures per cell
  3. Compare responder vs non-responder TME composition
  4. Identify spatial biomarkers of IO response
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
from scipy.stats import mannwhitneyu, kruskal
from collections import Counter

PROJECT = Path(__file__).resolve().parents[2]
DATA = PROJECT / "data"
OUT = PROJECT / "out"
FIG_DIR = OUT / "figures"
TABLE_DIR = OUT / "tables"
for d in [FIG_DIR, TABLE_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ── IO-relevant gene signatures ────────────────────────────────────────
IO_SIGNATURES = {
    "T_cell_exhaustion": ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX", "ENTPD1"],
    "T_cell_activation": ["IFNG", "GZMB", "PRF1", "TNF", "CD69", "IL2"],
    "T_cell_cytotoxicity": ["GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY", "NKG7"],
    "Treg_signature": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"],
    "M1_macrophage": ["NOS2", "IL1B", "IL6", "TNF", "CD80", "CD86"],
    "M2_macrophage": ["CD163", "MRC1", "MSR1", "CD200R1", "IL10"],
    "IFN_gamma": ["IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "IDO1"],
    "TGF_beta": ["TGFB1", "TGFB2", "TGFB3", "ACTA2", "COL1A1", "SERPINE1"],
    "Angiogenesis": ["VEGFA", "PECAM1", "KDR", "FLT1", "ANGPT2", "FGF2"],
    "Antigen_presentation": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB9"],
    "WNT_beta_catenin": ["CTNNB1", "APC", "AXIN1", "AXIN2", "LGR5", "MYC"],
    "Immune_checkpoint": ["CD274", "PDCD1LG2", "PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "VSIR", "SIGLEC15"],
}

CELL_TYPE_MARKERS = {
    "T_CD8": ["CD8A", "CD8B", "GZMB", "PRF1"],
    "T_CD4": ["CD4", "IL7R", "CCR7"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRB1"],
    "B_cell": ["CD19", "CD79A", "MS4A1", "CD79B"],
    "Macrophage": ["CD68", "CD163", "CSF1R", "MARCO"],
    "DC": ["CLEC9A", "XCR1", "BATF3", "IRF8"],
    "Hepatocyte": ["ALB", "APOB", "HP", "TF"],
    "Endothelial": ["PECAM1", "CDH5", "VWF", "ERG"],
    "Fibroblast": ["COL1A1", "COL3A1", "ACTA2", "DCN"],
}


def load_10x_mtx(data_dir, prefix=""):
    """Load 10x-style matrix (barcodes, genes/features, matrix)."""
    barcodes_files = sorted(data_dir.glob(f"{prefix}*barcodes*"))
    matrix_files = sorted(data_dir.glob(f"{prefix}*matrix*"))
    gene_files = sorted(data_dir.glob(f"{prefix}*genes*")) + sorted(data_dir.glob(f"{prefix}*features*"))

    if not (barcodes_files and matrix_files and gene_files):
        return None

    barcodes = pd.read_csv(barcodes_files[0], header=None, sep="\t")
    genes = pd.read_csv(gene_files[0], header=None, sep="\t")
    mat = sc.read_mtx(matrix_files[0]).T

    mat.obs_names = barcodes[0].values.astype(str)
    if genes.shape[1] >= 2:
        mat.var_names = genes[1].values.astype(str)
        mat.var["gene_ids"] = genes[0].values.astype(str)
    else:
        mat.var_names = genes[0].values.astype(str)

    # Deduplicate var names
    mat.var_names_make_unique()
    return mat


def annotate_cells(adata):
    """Simple marker-based cell type annotation."""
    gene_set = set(adata.var_names)
    scores = {}
    for ct, markers in CELL_TYPE_MARKERS.items():
        present = [m for m in markers if m in gene_set]
        if present:
            sc.tl.score_genes(adata, present, score_name=f"score_{ct}", use_raw=False)
            scores[ct] = f"score_{ct}"

    if not scores:
        adata.obs["cell_type"] = "Unknown"
        return

    # Assign based on max score
    score_df = adata.obs[[v for v in scores.values()]]
    adata.obs["cell_type"] = score_df.idxmax(axis=1).str.replace("score_", "")
    # Low-confidence → Unknown
    max_scores = score_df.max(axis=1)
    adata.obs.loc[max_scores < 0.1, "cell_type"] = "Unknown"


def score_io_signatures(adata):
    """Score all IO signatures."""
    gene_set = set(adata.var_names)
    scored = []
    for sig_name, sig_genes in IO_SIGNATURES.items():
        present = [g for g in sig_genes if g in gene_set]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=sig_name, use_raw=False)
            scored.append(sig_name)
            print(f"    {sig_name}: {len(present)}/{len(sig_genes)} genes present")
    return scored


def main():
    print("=" * 70)
    print("Spatial Immunotherapy Response — Public Data Analysis")
    print("=" * 70)

    all_datasets = []

    # ─── Dataset 1: GSE125449 (Ma et al., 2019) ──────────────────────────
    print("\n[1] Loading GSE125449 (Ma et al., 2019 — HCC anti-PD1 scRNA-seq)...")
    gse125449_dir = DATA / "pre_treatment" / "GSE125449_Ma2019_antiPD1"

    for set_name in ["Set1", "Set2"]:
        barcodes_f = gse125449_dir / f"GSE125449_{set_name}_barcodes.tsv.gz"
        genes_f = gse125449_dir / f"GSE125449_{set_name}_genes.tsv.gz"
        matrix_f = gse125449_dir / f"GSE125449_{set_name}_matrix.mtx.gz"
        samples_f = gse125449_dir / f"GSE125449_{set_name}_samples.txt.gz"

        if not all(f.exists() for f in [barcodes_f, genes_f, matrix_f]):
            print(f"  [SKIP] {set_name} files not found")
            continue

        barcodes = pd.read_csv(barcodes_f, header=None, sep="\t")
        genes = pd.read_csv(genes_f, header=None, sep="\t")
        mat = sc.read_mtx(matrix_f).T

        mat.obs_names = barcodes[0].values.astype(str)
        if genes.shape[1] >= 2:
            mat.var_names = genes[1].values.astype(str)
        else:
            mat.var_names = genes[0].values.astype(str)
        mat.var_names_make_unique()

        # Load sample info if available
        if samples_f.exists():
            samples = pd.read_csv(samples_f, sep="\t")
            print(f"  {set_name}: samples metadata shape {samples.shape}")
            print(f"  Columns: {list(samples.columns)}")

        mat.obs["dataset"] = f"GSE125449_{set_name}"
        mat.obs["study"] = "Ma2019"
        print(f"  {set_name}: {mat.shape[0]} cells × {mat.shape[1]} genes")
        all_datasets.append(mat)

    # ─── Dataset 2: GSE234241 (Xue/Zhang et al., 2023) ───────────────────
    print("\n[2] Loading GSE234241 (Xue et al., 2023 — HCC combination IO)...")
    gse234241_dir = DATA / "pre_treatment" / "GSE234241_Xue2023_combo"
    adata234 = load_10x_mtx(gse234241_dir, prefix="GSE234241_")
    if adata234 is not None:
        # Load metadata
        meta_f = gse234241_dir / "GSE234241_Metadata.txt.gz"
        if meta_f.exists():
            meta = pd.read_csv(meta_f, sep="\t", index_col=0)
            print(f"  Metadata: {meta.shape[0]} cells, columns: {list(meta.columns)[:10]}")
            # Merge
            common = adata234.obs_names.intersection(meta.index)
            if len(common) > 0:
                adata234 = adata234[common].copy()
                for col in meta.columns:
                    adata234.obs[col] = meta.loc[common, col].values

        adata234.obs["dataset"] = "GSE234241"
        adata234.obs["study"] = "Xue2023"
        print(f"  Loaded: {adata234.shape[0]} cells × {adata234.shape[1]} genes")
        all_datasets.append(adata234)
    else:
        print("  [SKIP] Could not load GSE234241")

    # ─── Dataset 3: GSE202069 (Zhu et al., 2023 — Bulk RNA-seq) ──────────
    print("\n[3] Loading GSE202069 (Zhu et al., 2023 — bulk RNA-seq anti-PDL1)...")
    bulk_f = DATA / "pre_treatment" / "GSE202069_Zhu2023_antiPDL1" / "GSE202069_gene_tpm_expression.txt.gz"
    bulk_data = None
    if bulk_f.exists():
        bulk_data = pd.read_csv(bulk_f, sep="\t", index_col=0)
        print(f"  Bulk data: {bulk_data.shape[0]} genes × {bulk_data.shape[1]} samples")
    else:
        print("  [SKIP] Bulk data not found")

    # ─── Merge scRNA-seq datasets ─────────────────────────────────────────
    print("\n[4] Merging scRNA-seq datasets...")
    if len(all_datasets) == 0:
        print("  [ERROR] No datasets loaded!")
        return

    if len(all_datasets) == 1:
        adata = all_datasets[0].copy()
    else:
        # Find shared genes
        shared_genes = set(all_datasets[0].var_names)
        for ad in all_datasets[1:]:
            shared_genes = shared_genes.intersection(set(ad.var_names))
        shared_genes = sorted(shared_genes)
        print(f"  Shared genes across datasets: {len(shared_genes)}")

        # Subset and concatenate
        subsets = [ad[:, shared_genes].copy() for ad in all_datasets]
        adata = sc.concat(subsets, join="inner")

    print(f"  Combined: {adata.shape[0]} cells × {adata.shape[1]} genes")

    # ─── Preprocess ───────────────────────────────────────────────────────
    print("\n[5] Preprocessing...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=10)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print(f"  After filtering: {adata.shape}")

    # HVGs, PCA, neighbors, UMAP, Leiden
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="dataset" if "dataset" in adata.obs else None)
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
    sc.pp.neighbors(adata, n_pcs=30)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.8)
    print(f"  Leiden clusters: {adata.obs['leiden'].nunique()}")

    # ─── Cell type annotation ─────────────────────────────────────────────
    print("\n[6] Cell type annotation...")
    annotate_cells(adata)
    ct_counts = adata.obs["cell_type"].value_counts()
    print("  Cell type distribution:")
    for ct, n in ct_counts.items():
        print(f"    {ct}: {n} ({100*n/len(adata):.1f}%)")

    # ─── IO signature scoring ─────────────────────────────────────────────
    print("\n[7] Scoring IO signatures...")
    scored_sigs = score_io_signatures(adata)
    print(f"  Scored {len(scored_sigs)} signatures")

    # ─── Response group analysis ──────────────────────────────────────────
    print("\n[8] Analyzing response groups...")
    # Try to find response labels in metadata
    response_col = None
    for col in adata.obs.columns:
        col_lower = col.lower()
        if any(kw in col_lower for kw in ["response", "responder", "treatment", "therapy"]):
            print(f"  Found potential response column: {col}")
            print(f"    Values: {adata.obs[col].value_counts().to_dict()}")
            response_col = col
            break

    # If no explicit response labels, create pseudo-response groups based on immune activity
    if response_col is None:
        print("  No explicit response labels found. Creating immune-based phenotypes...")
        if "T_cell_activation" in scored_sigs and "Immune_checkpoint" in scored_sigs:
            # Classify by immune activity
            act_score = adata.obs["T_cell_activation"]
            ckpt_score = adata.obs["Immune_checkpoint"]
            act_med = act_score.median()
            ckpt_med = ckpt_score.median()

            conditions = [
                (act_score > act_med) & (ckpt_score > ckpt_med),
                (act_score > act_med) & (ckpt_score <= ckpt_med),
                (act_score <= act_med) & (ckpt_score > ckpt_med),
                (act_score <= act_med) & (ckpt_score <= ckpt_med),
            ]
            labels = ["Hot_Checkpoint_High", "Hot_Checkpoint_Low", "Cold_Checkpoint_High", "Cold_Checkpoint_Low"]
            adata.obs["io_phenotype"] = "Unclassified"
            for cond, label in zip(conditions, labels):
                adata.obs.loc[cond, "io_phenotype"] = label
            response_col = "io_phenotype"
            print(f"  IO phenotype distribution:")
            for label, n in adata.obs[response_col].value_counts().items():
                print(f"    {label}: {n}")

    # ─── Differential signature analysis between groups ───────────────────
    print("\n[9] Differential analysis between IO phenotypes...")
    if response_col and len(scored_sigs) > 0:
        groups = adata.obs[response_col].unique()
        diff_results = []

        for sig in scored_sigs:
            for i, g1 in enumerate(groups):
                for g2 in groups[i+1:]:
                    vals1 = adata.obs.loc[adata.obs[response_col] == g1, sig].dropna()
                    vals2 = adata.obs.loc[adata.obs[response_col] == g2, sig].dropna()
                    if len(vals1) >= 10 and len(vals2) >= 10:
                        stat, pval = mannwhitneyu(vals1, vals2, alternative="two-sided")
                        effect = vals1.mean() - vals2.mean()
                        diff_results.append({
                            "signature": sig,
                            "group1": g1, "group2": g2,
                            "mean_g1": vals1.mean(), "mean_g2": vals2.mean(),
                            "effect_size": effect, "pvalue": pval,
                            "n_g1": len(vals1), "n_g2": len(vals2),
                        })

        diff_df = pd.DataFrame(diff_results)
        if len(diff_df) > 0:
            diff_df = diff_df.sort_values("pvalue")
            diff_df.to_csv(TABLE_DIR / "differential_signatures.csv", index=False)
            print(f"  Top differential signatures:")
            for _, row in diff_df.head(10).iterrows():
                print(f"    {row['signature']}: {row['group1']} vs {row['group2']} "
                      f"(effect={row['effect_size']:.3f}, p={row['pvalue']:.2e})")

    # ─── Cell type composition by IO phenotype ────────────────────────────
    print("\n[10] Cell type composition by IO phenotype...")
    if response_col:
        comp_df = pd.crosstab(adata.obs[response_col], adata.obs["cell_type"], normalize="index")
        comp_df.to_csv(TABLE_DIR / "celltype_composition_by_phenotype.csv")
        print(comp_df.round(3).to_string())

    # ─── Bulk RNA-seq analysis (if available) ─────────────────────────────
    bulk_scored = {}
    if bulk_data is not None:
        print("\n[11] Scoring IO signatures on bulk RNA-seq...")
        gene_set_bulk = set(bulk_data.index)
        for sig_name, sig_genes in IO_SIGNATURES.items():
            present = [g for g in sig_genes if g in gene_set_bulk]
            if len(present) >= 2:
                bulk_scored[sig_name] = bulk_data.loc[present].mean(axis=0)
                print(f"    {sig_name}: {len(present)} genes")
        if bulk_scored:
            bulk_sig_df = pd.DataFrame(bulk_scored)
            bulk_sig_df.to_csv(TABLE_DIR / "bulk_io_signatures.csv")
            print(f"  Bulk signatures saved for {len(bulk_sig_df)} samples")

    # ─── FIGURES ──────────────────────────────────────────────────────────
    print("\n[12] Generating figures...")

    fig, axes = plt.subplots(2, 3, figsize=(20, 13))

    # Panel A: UMAP by dataset
    ax = axes[0, 0]
    datasets = adata.obs["dataset"].unique()
    colors_ds = plt.cm.Set2(np.linspace(0, 1, len(datasets)))
    for i, ds in enumerate(datasets):
        mask = adata.obs["dataset"] == ds
        ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                   s=1, alpha=0.3, c=[colors_ds[i]], label=ds, rasterized=True)
    ax.legend(fontsize=7, markerscale=5, loc="upper right")
    ax.set_title("A. UMAP — Dataset", fontsize=11)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    # Panel B: UMAP by cell type
    ax = axes[0, 1]
    cts = adata.obs["cell_type"].unique()
    colors_ct = plt.cm.tab20(np.linspace(0, 1, len(cts)))
    for i, ct in enumerate(cts):
        mask = adata.obs["cell_type"] == ct
        ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                   s=1, alpha=0.3, c=[colors_ct[i]], label=ct, rasterized=True)
    ax.legend(fontsize=6, markerscale=5, ncol=2, loc="upper right")
    ax.set_title("B. UMAP — Cell Types", fontsize=11)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    # Panel C: UMAP by IO phenotype
    ax = axes[0, 2]
    if response_col:
        groups = adata.obs[response_col].unique()
        colors_ph = plt.cm.RdYlBu(np.linspace(0, 1, len(groups)))
        for i, g in enumerate(groups):
            mask = adata.obs[response_col] == g
            ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
                       s=1, alpha=0.3, c=[colors_ph[i]], label=g, rasterized=True)
        ax.legend(fontsize=7, markerscale=5, loc="upper right")
    ax.set_title("C. UMAP — IO Phenotype", fontsize=11)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    # Panel D: IO signature heatmap
    ax = axes[1, 0]
    if response_col and scored_sigs:
        sig_means = adata.obs.groupby(response_col)[scored_sigs].mean()
        sns.heatmap(sig_means.T, ax=ax, cmap="RdBu_r", center=0,
                    xticklabels=True, yticklabels=True, cbar_kws={"shrink": 0.8})
        ax.set_title("D. IO Signatures by Phenotype", fontsize=11)
    else:
        ax.set_title("D. IO Signatures (no groups available)")

    # Panel E: Cell type composition
    ax = axes[1, 1]
    if response_col:
        comp_df_plot = pd.crosstab(adata.obs[response_col], adata.obs["cell_type"], normalize="index")
        comp_df_plot.plot(kind="bar", stacked=True, ax=ax, colormap="tab20", legend=False)
        ax.legend(fontsize=6, bbox_to_anchor=(1.05, 1), loc="upper left")
        ax.set_title("E. Cell Type Composition", fontsize=11)
        ax.set_ylabel("Fraction")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    # Panel F: Top differential signatures
    ax = axes[1, 2]
    if 'diff_df' in dir() and len(diff_df) > 0:
        top10 = diff_df.head(10).copy()
        top10["neg_log_p"] = -np.log10(top10["pvalue"].clip(lower=1e-300))
        top10["label"] = top10["signature"] + "\n" + top10["group1"].str[:8] + " vs " + top10["group2"].str[:8]
        colors_bar = ["#e74c3c" if e > 0 else "#3498db" for e in top10["effect_size"]]
        ax.barh(range(len(top10)), top10["neg_log_p"].values, color=colors_bar)
        ax.set_yticks(range(len(top10)))
        ax.set_yticklabels(top10["label"].values, fontsize=7)
        ax.set_xlabel("-log10(p-value)")
        ax.set_title("F. Top Differential Signatures", fontsize=11)
    else:
        ax.text(0.5, 0.5, "No differential results", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("F. Top Differential Signatures", fontsize=11)

    plt.suptitle("Immunotherapy Response Biomarker Analysis\n(Public scRNA-seq Data)", fontsize=14, y=1.02)
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG_DIR / f"io_biomarker_public.{ext}", dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved io_biomarker_public.png/pdf")

    # ─── Bulk heatmap if available ────────────────────────────────────────
    if bulk_scored:
        fig, ax = plt.subplots(figsize=(12, 6))
        bulk_sig_df_plot = pd.DataFrame(bulk_scored)
        # Cluster samples
        from scipy.cluster.hierarchy import linkage, leaves_list
        Z = linkage(bulk_sig_df_plot.values, method="ward")
        order = leaves_list(Z)
        sns.heatmap(bulk_sig_df_plot.iloc[order].T, ax=ax, cmap="RdBu_r", center=0,
                    xticklabels=True, yticklabels=True)
        ax.set_title("IO Signatures — Bulk RNA-seq (GSE202069)")
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG_DIR / f"bulk_io_heatmap.{ext}", dpi=200, bbox_inches="tight")
        plt.close()
        print(f"  Saved bulk_io_heatmap.png/pdf")

    # ─── Save processed data ─────────────────────────────────────────────
    adata.write(OUT / "processed_io_public.h5ad")
    print(f"\n  Processed data saved: {adata.shape}")

    print("\n" + "=" * 70)
    print("[DONE] Immunotherapy Response Analysis — Public Data")
    print(f"  Total cells: {adata.shape[0]}")
    print(f"  Datasets: {adata.obs['dataset'].nunique()}")
    print(f"  Cell types: {adata.obs['cell_type'].nunique()}")
    print(f"  IO signatures scored: {len(scored_sigs)}")
    print(f"  Figures: {FIG_DIR}")
    print(f"  Tables: {TABLE_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
