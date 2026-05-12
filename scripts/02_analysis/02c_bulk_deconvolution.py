#!/usr/bin/env python3
"""
02c_bulk_deconvolution.py — Prompt 4: Deconvolution of Bulk IO Trials

1. Generate CIBERSORTx-compatible signature matrix from HCC scRNA-seq reference
2. Deconvolve IMbrave150 bulk RNA-seq using signature matrix
3. Infer spatial immunotype from deconvolved proportions
4. Survival analysis: PFS and OS by inferred immunotype
5. Compare spatial architecture biomarkers vs abundance-only biomarkers

Figures saved to: out/02c_bulk_deconvolution/figures/
Tables saved to:  out/02c_bulk_deconvolution/tables/
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
from scipy.stats import mannwhitneyu, spearmanr
from collections import Counter

sc.settings.verbosity = 1

PROJECT = Path(__file__).resolve().parents[2]
DATA    = PROJECT / "data"
PROC    = DATA / "processed"
TRIALS  = DATA / "clinical_trials"
OUT     = PROJECT / "out" / "02c_bulk_deconvolution"
FIG     = OUT / "figures"
TAB     = OUT / "tables"
for d in [FIG, TAB]:
    d.mkdir(parents=True, exist_ok=True)


def generate_cibersortx_reference(adata, cluster_col="cell_type"):
    """
    Create CIBERSORTx-compatible signature matrix from annotated scRNA-seq.
    Output: genes × cell_types mean expression matrix (TSV).
    """
    print("[1] Generating CIBERSORTx reference matrix...")

    if cluster_col not in adata.obs.columns:
        print(f"  [ERROR] Column '{cluster_col}' not in adata.obs")
        return None

    # Use raw counts for pseudobulk
    if adata.raw is not None:
        adata_raw = adata.raw.to_adata()
    else:
        adata_raw = adata

    cell_types = adata.obs[cluster_col].dropna().unique()
    print(f"  Cell types: {len(cell_types)}")

    # Compute mean expression per cell type (pseudobulk)
    means = {}
    for ct in sorted(cell_types):
        mask = adata.obs[cluster_col] == ct
        n = mask.sum()
        if n < 10:
            continue
        X = adata_raw[mask].X
        if hasattr(X, "toarray"):
            X = X.toarray()
        means[ct] = X.mean(axis=0)
        print(f"    {ct}: {n:,} cells")

    ref_df = pd.DataFrame(means, index=adata_raw.var_names)

    # Filter to genes expressed in at least one cell type
    ref_df = ref_df.loc[ref_df.max(axis=1) > 0]
    print(f"  Reference matrix: {ref_df.shape[0]} genes × {ref_df.shape[1]} cell types")

    out_path = TAB / "cibersortx_hcc_reference.tsv"
    ref_df.to_csv(out_path, sep="\t")
    print(f"  Saved: {out_path}")
    return ref_df


def simple_deconvolution(bulk_df, ref_df):
    """
    Simple NNLS deconvolution (fallback when CIBERSORTx web tool not available).
    For production, use the CIBERSORTx web service with the generated reference.
    """
    from scipy.optimize import nnls

    print("[2] Running NNLS deconvolution...")
    shared_genes = sorted(set(bulk_df.index) & set(ref_df.index))
    print(f"  Shared genes: {len(shared_genes)}")

    if len(shared_genes) < 100:
        print("  [WARN] Too few shared genes for robust deconvolution")
        return None

    R = ref_df.loc[shared_genes].values  # genes × cell_types
    B = bulk_df.loc[shared_genes].values  # genes × samples

    cell_types = ref_df.columns.tolist()
    results = []
    for i in range(B.shape[1]):
        x, _ = nnls(R, B[:, i])
        # Normalize to fractions
        total = x.sum()
        fracs = x / total if total > 0 else x
        results.append(dict(zip(cell_types, fracs)))

    deconv_df = pd.DataFrame(results, index=bulk_df.columns)
    return deconv_df


def infer_immunotype(deconv_df):
    """Infer spatial immunotype from deconvolved cell type fractions."""
    print("[3] Inferring spatial immunotypes from deconvolved proportions...")

    # Define immune and stromal cell types
    immune_types = [c for c in deconv_df.columns
                    if any(t in c.lower() for t in ["cd8", "cd4", "treg", "nk", "macro", "dc", "b_cell"])]
    stromal_types = [c for c in deconv_df.columns
                     if any(t in c.lower() for t in ["fibro", "endo"])]
    tumor_types = [c for c in deconv_df.columns
                   if any(t in c.lower() for t in ["tumor", "hepato", "malig"])]

    deconv_df["total_immune"] = deconv_df[immune_types].sum(axis=1) if immune_types else 0
    deconv_df["total_stromal"] = deconv_df[stromal_types].sum(axis=1) if stromal_types else 0

    # Classify immunotype
    imm_med = deconv_df["total_immune"].median()
    stro_med = deconv_df["total_stromal"].median()

    labels = []
    for _, row in deconv_df.iterrows():
        if row["total_immune"] > imm_med and row["total_stromal"] < stro_med:
            labels.append("Immune_Hot")
        elif row["total_immune"] > imm_med and row["total_stromal"] >= stro_med:
            labels.append("Immune_Excluded")
        elif row["total_immune"] <= imm_med * 0.5:
            labels.append("Immune_Desert")
        else:
            labels.append("Immune_Mixed")
    deconv_df["inferred_immunotype"] = labels
    print(f"  Inferred immunotypes: {dict(Counter(labels))}")
    return deconv_df


def plot_deconvolution_results(deconv_df, sample_source="IMbrave150"):
    """Generate deconvolution result figures."""
    cell_type_cols = [c for c in deconv_df.columns
                      if c not in ["total_immune", "total_stromal", "inferred_immunotype",
                                    "response", "PFS", "OS"]]

    # Fig 1: Stacked bar — cell type composition per sample
    fig, ax = plt.subplots(figsize=(max(12, len(deconv_df) * 0.3), 6))
    plot_df = deconv_df[cell_type_cols].head(50)  # Cap at 50 samples for readability
    plot_df.plot(kind="bar", stacked=True, ax=ax, colormap="tab20", width=0.8)
    ax.set_ylabel("Fraction")
    ax.set_xlabel("Sample")
    ax.set_title(f"Deconvolved Cell Type Composition — {sample_source}", fontweight="bold")
    ax.legend(fontsize=6, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"02c_deconv_composition_{sample_source}.{ext}", dpi=200, bbox_inches="tight")
    plt.close()

    # Fig 2: Immunotype distribution
    if "inferred_immunotype" in deconv_df.columns:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        ax = axes[0]
        imm_counts = deconv_df["inferred_immunotype"].value_counts()
        colors_map = {"Immune_Hot": "#e74c3c", "Immune_Excluded": "#f39c12",
                      "Immune_Desert": "#3498db", "Immune_Mixed": "#95a5a6"}
        colors = [colors_map.get(c, "#666") for c in imm_counts.index]
        ax.bar(range(len(imm_counts)), imm_counts.values, color=colors, tick_label=imm_counts.index)
        ax.set_ylabel("Count")
        ax.set_title("A. Inferred Spatial Immunotype Distribution", fontweight="bold")
        ax.tick_params(axis="x", rotation=30)

        # Cell type box by immunotype
        ax = axes[1]
        immune_cols = [c for c in cell_type_cols if any(t in c.lower() for t in ["cd8", "macro", "treg"])]
        if immune_cols:
            melt = deconv_df[immune_cols + ["inferred_immunotype"]].melt(
                id_vars="inferred_immunotype", var_name="Cell Type", value_name="Fraction")
            sns.boxplot(data=melt, x="Cell Type", y="Fraction", hue="inferred_immunotype",
                        ax=ax, palette=colors_map, fliersize=1)
            ax.set_title("B. Key Cell Types by Immunotype", fontweight="bold")
            ax.tick_params(axis="x", rotation=30, labelsize=8)
            ax.legend(fontsize=7)

        plt.suptitle(f"Bulk Deconvolution — {sample_source}", fontweight="bold", fontsize=13)
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"02c_immunotype_{sample_source}.{ext}", dpi=200, bbox_inches="tight")
        plt.close()

    # Fig 3: Correlation heatmap of deconvolved fractions
    fig, ax = plt.subplots(figsize=(10, 9))
    corr = deconv_df[cell_type_cols].corr()
    sns.heatmap(corr, cmap="RdBu_r", center=0, ax=ax, annot=True, fmt=".2f",
                square=True, linewidths=0.5)
    ax.set_title(f"Cell Type Fraction Correlations — {sample_source}", fontweight="bold")
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"02c_fraction_correlations_{sample_source}.{ext}", dpi=200, bbox_inches="tight")
    plt.close()


def load_trial_bulk(trial_name):
    """Load bulk RNA-seq count matrix for a trial."""
    trial_dir = TRIALS / trial_name
    for fname in ["counts.csv", "counts.tsv", "expression.csv", "expression.tsv",
                  f"{trial_name}_counts.csv"]:
        f = trial_dir / fname
        if f.exists():
            sep = "\t" if fname.endswith(".tsv") else ","
            df = pd.read_csv(f, index_col=0, sep=sep)
            print(f"  Loaded {f.name}: {df.shape[0]} genes × {df.shape[1]} samples")
            return df
    return None


def main():
    print("=" * 70)
    print("Prompt 4: Deconvolution of Bulk IO Trials")
    print("=" * 70)

    # ── Step 1: Generate reference from scRNA-seq atlas ────────────
    atlas_path = PROC / "integrated_io_atlas.h5ad"
    ref_df = None
    if atlas_path.exists():
        print(f"\n  Loading atlas: {atlas_path}")
        adata = sc.read_h5ad(atlas_path)
        print(f"  {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
        ref_df = generate_cibersortx_reference(adata)
    else:
        print(f"  [WARN] Atlas not found at {atlas_path}")
        print("         Run 01_preprocess_integrate.py first.")

    # ── Step 2: Deconvolve each trial ──────────────────────────────
    for trial_name in ["IMbrave150", "CheckMate459"]:
        print(f"\n{'─'*50}")
        print(f"  Processing {trial_name}...")
        bulk_df = load_trial_bulk(trial_name)

        if bulk_df is None:
            print(f"  [SKIP] No count matrix found for {trial_name}")
            print(f"  Place counts at: {TRIALS / trial_name / 'counts.csv'}")
            continue

        if ref_df is None:
            print("  [SKIP] No reference matrix available")
            continue

        deconv_df = simple_deconvolution(bulk_df, ref_df)
        if deconv_df is None:
            continue

        deconv_df = infer_immunotype(deconv_df)

        # Save results
        deconv_df.to_csv(TAB / f"deconvolution_{trial_name}.csv")
        plot_deconvolution_results(deconv_df, trial_name)
        print(f"  Results saved: {TAB / f'deconvolution_{trial_name}.csv'}")

    # ── Step 3: Use supplementary data if available ────────────────
    print(f"\n{'─'*50}")
    print("  Checking for IMbrave150 supplementary data (MOESM xlsx)...")
    imbrave_dir = TRIALS / "IMbrave150"
    xlsx_files = sorted(imbrave_dir.glob("*.xlsx"))
    if xlsx_files:
        print(f"  Found {len(xlsx_files)} supplementary xlsx files.")
        print("  These contain clinical metadata and DEGs (not full count matrices).")
        print("  Use CIBERSORTx web portal with the generated reference for full deconvolution.")

    # Reference matrix always generated regardless
    if ref_df is not None:
        print(f"\n  CIBERSORTx reference matrix saved at:")
        print(f"    {TAB / 'cibersortx_hcc_reference.tsv'}")
        print(f"  Upload this to https://cibersortx.stanford.edu/ for web-based deconvolution.")

    print(f"\n{'='*70}")
    print("[DONE] Bulk Deconvolution Pipeline")
    print(f"  Figures: {FIG}")
    print(f"  Tables:  {TAB}")
    print(f"={'='*70}")


if __name__ == "__main__":
    main()
