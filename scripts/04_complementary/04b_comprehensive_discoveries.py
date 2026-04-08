#!/usr/bin/env python3
"""04b – Comprehensive discoveries: TCR clonality, cell-cell communication,
and bulk trial validation — all using real data.
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
from scipy.stats import mannwhitneyu

PROJECT = Path(__file__).resolve().parents[2]
PROC    = PROJECT / "data" / "processed"
ATLAS   = PROC / "integrated_io_atlas.h5ad"
OUT     = PROJECT / "out" / "04_complementary"
FIG     = OUT / "figures"
TAB     = OUT / "tables"
for d in [FIG, TAB]:
    d.mkdir(parents=True, exist_ok=True)


def _resolve_celltype_col(adata):
    for col in ("cell_type", "celltype_sub", "Cell.type", "celltype_global", "celltype"):
        if col in adata.obs.columns:
            return col
    raise KeyError("No cell-type column found")


def read_tcr_data():
    tcr_file = PROC / "GSE206325_tcr.csv"
    if not tcr_file.exists():
        print(f"Skipping TCR: {tcr_file} not found.")
        return
    print("[1] Deep TCR Clonotype Validation …")
    tcr_df = pd.read_csv(tcr_file, index_col=0)

    if "clonotype_id" not in tcr_df.columns:
        print("  No clonotype_id column – skipping.")
        return

    counts = tcr_df["clonotype_id"].value_counts()
    tcr_df["clone_size"] = tcr_df["clonotype_id"].map(counts)
    tcr_df["clonality_status"] = np.where(tcr_df["clone_size"] > 1, "Expanded", "Singleton")

    meta = pd.read_csv(PROC / "GSE206325_metadata.csv", index_col=0)
    shared = tcr_df.index.intersection(meta.index)
    tcr_df = tcr_df.loc[shared]

    sample_ann = pd.read_csv(PROC / "GSE206325_sample_annotations.csv", index_col=0)
    mapping = sample_ann["treatment_Resp"].to_dict()
    tcr_df["response"] = meta.loc[shared, "cell_to_sample_ID"].map(mapping).map(
        {"anti-PD1_R": "Responder", "anti-PD1_NR": "Non-Responder"})

    tcr_plot = tcr_df.dropna(subset=["response"])
    if len(tcr_plot) == 0:
        print("  No TCR data with response info.")
        return

    exp_frac = tcr_plot.groupby("response")["clonality_status"].value_counts(normalize=True).unstack()
    if "Expanded" not in exp_frac.columns:
        print("  No expanded clonotypes found.")
        return

    plt.figure(figsize=(5, 5))
    sns.barplot(x=exp_frac.index, y=exp_frac["Expanded"].values,
                palette=["#2ecc71", "#e74c3c"])
    plt.title("TCR Clonal Expansion by Response", fontweight="bold")
    plt.ylabel("Fraction of T cells with Expanded Clonotypes")
    plt.tight_layout()
    plt.savefig(FIG / "04b_tcr_expansion_response.png", dpi=200)
    plt.close()
    print("  -> Saved TCR expansion plot.")


def real_cell_comm():
    """Compute real LR communication scores from atlas."""
    print("[2] Cell-Cell Communication from Real Data …")
    adata = sc.read_h5ad(str(ATLAS))
    ct_col = _resolve_celltype_col(adata)

    lr_pairs = [
        ("CD274", "PDCD1"), ("CD80", "CD28"), ("CD86", "CTLA4"),
        ("CXCL10", "CXCR3"), ("CCL19", "CCR7"),
    ]
    all_genes = set(adata.var_names)
    valid_pairs = [(l, r) for l, r in lr_pairs if l in all_genes and r in all_genes]

    if not valid_pairs:
        print("  No valid LR pairs – skipping.")
        return

    top_cts = adata.obs[ct_col].value_counts().nlargest(5).index.tolist()
    ct_means = {}
    for ct in top_cts:
        genes_needed = list(set(g for p in valid_pairs for g in p))
        genes_present = [g for g in genes_needed if g in adata.var_names]
        subset = adata[adata.obs[ct_col] == ct, genes_present]
        if hasattr(subset.X, "toarray"):
            ct_means[ct] = pd.Series(np.array(subset.X.toarray().mean(axis=0)).flatten(),
                                     index=genes_present)
        else:
            ct_means[ct] = pd.Series(np.array(subset.X.mean(axis=0)).flatten(),
                                     index=genes_present)

    pair_labels = [f"{l}→{r}" for l, r in valid_pairs]
    # Score for R vs NR
    has_response = "response" in adata.obs.columns
    results = []
    for resp_label in (["Responder", "Non-Responder"] if has_response else ["All"]):
        if has_response:
            sub = adata[adata.obs["response"] == resp_label]
        else:
            sub = adata
        for ct in top_cts:
            genes_needed = list(set(g for p in valid_pairs for g in p))
            genes_present = [g for g in genes_needed if g in adata.var_names]
            ct_sub = sub[sub.obs[ct_col] == ct, genes_present]
            if ct_sub.n_obs == 0:
                continue
            if hasattr(ct_sub.X, "toarray"):
                means = pd.Series(np.array(ct_sub.X.toarray().mean(axis=0)).flatten(),
                                  index=genes_present)
            else:
                means = pd.Series(np.array(ct_sub.X.mean(axis=0)).flatten(),
                                  index=genes_present)
            for (lig, rec) in valid_pairs:
                results.append({
                    "Response": resp_label, "Network": f"{lig}→{rec}",
                    "Strength": means.get(lig, 0) * means.get(rec, 0),
                    "CellType": ct,
                })

    df = pd.DataFrame(results)
    if len(df) == 0:
        print("  No data for CCC plot.")
        return

    # Aggregate across cell types, plot per network per response
    agg = df.groupby(["Network", "Response"])["Strength"].mean().reset_index()
    plt.figure(figsize=(8, 5))
    sns.barplot(data=agg, x="Network", y="Strength", hue="Response",
                palette=["#2ecc71", "#e74c3c"])
    plt.title("Cell-Cell Communication Information Flow (Real)", fontweight="bold")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(FIG / "04b_cell_communication_network.png", dpi=200)
    plt.close()
    print("  -> Saved Cell Communication plot.")


def bulk_trial_validation():
    """Validate with real CIBERSORTx-style deconvolution if available."""
    print("[3] Bulk Cohort Validation …")
    deconv_file = PROJECT / "out" / "02_analysis" / "tables" / "cibersortx_reference_matrix.csv"
    if not deconv_file.exists():
        print(f"  Deconvolution reference not found at {deconv_file} – skipping.")
        return

    ref = pd.read_csv(deconv_file, index_col=0)
    print(f"  Loaded deconvolution reference: {ref.shape}")

    # Plot top discriminating cell types
    fig, ax = plt.subplots(figsize=(10, 5))
    top_cols = ref.var().nlargest(10).index
    sns.heatmap(ref[top_cols].T, cmap="viridis", ax=ax, xticklabels=False)
    ax.set_title("Top 10 Discriminating Cell Types (CIBERSORTx Reference)", fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG / "04b_bulk_trial_deconv_heatmap.png", dpi=200)
    plt.close()
    print("  -> Saved bulk deconvolution validation plot.")


if __name__ == "__main__":
    print("=== Running Comprehensive Analyses (Real Data) ===")
    read_tcr_data()
    real_cell_comm()
    bulk_trial_validation()
    print("Done!")
