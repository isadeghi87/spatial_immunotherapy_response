#!/usr/bin/env python3
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

# Directories
PROJECT = Path(__file__).resolve().parents[2]
OUT_02A = PROJECT / "out" / "02a_response_tme"
OUT_02B = PROJECT / "out" / "02b_spatial_biomarkers"
OUT_02C = PROJECT / "out" / "02c_bulk_deconvolution"

def plot_deg_heatmap():
    deg_path = OUT_02A / "tables" / "deg_R_vs_NR_per_celltype.csv"
    if not deg_path.exists():
        return
    print("Generating EDA for 02a DEGs...")
    df = pd.read_csv(deg_path)
    # Count signficant DEGs up/down
    df_sig = df[(df["pvals_adj"] < 0.05) & (df["logfoldchanges"].abs() > 0.25)]
    df_sig["direction"] = np.where(df_sig["logfoldchanges"] > 0, "Up in Responder", "Down in Responder")
    counts = df_sig.groupby(["cell_type", "direction"]).size().unstack(fill_value=0)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    counts.plot(kind="bar", stacked=True, ax=ax, color=["#e74c3c", "#2ecc71"])
    ax.set_title("Number of Significant DEGs per Cell Type\n(Responder vs Non-Responder)", fontweight="bold")
    ax.set_ylabel("Number of DEGs")
    ax.set_xlabel("Cell Type")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(OUT_02A / "figures" / "eda_deg_counts_barplot.png", dpi=200, bbox_inches="tight")
    plt.savefig(OUT_02A / "figures" / "eda_deg_counts_barplot.pdf", dpi=200, bbox_inches="tight")
    plt.close()
    
    # Top 5 DEGs logFC heatmap
    top_genes = df_sig.groupby("cell_type").apply(lambda x: x.nlargest(5, "logfoldchanges")["names"]).explode().unique()
    if len(top_genes) > 0:
        pivot = df[df["names"].isin(top_genes)].pivot_table(index="names", columns="cell_type", values="logfoldchanges").fillna(0)
        plt.figure(figsize=(10, max(6, len(top_genes)*0.2)))
        sns.heatmap(pivot, cmap="RdBu_r", center=0, cbar_kws={'label': 'Log2FC'})
        plt.title("Top Responder DEGs across Cell Types", fontweight="bold")
        plt.tight_layout()
        plt.savefig(OUT_02A / "figures" / "eda_top_degs_heatmap.png", dpi=200, bbox_inches="tight")
        plt.close()
        
def plot_spatial_biomarkers_eda():
    h5ad_path = OUT_02B / "tables" / "GSE173905_pseudobulk_scored.h5ad"
    if not h5ad_path.exists():
        return
    print("Generating EDA for 02b Spatial Biomarkers...")
    adata = sc.read_h5ad(h5ad_path)
    
    # Heatmap of signature scores
    sig_cols = [c for c in adata.obs.columns if c in ["T_cell_exhaustion", "T_cell_activation", "T_cell_cytotoxicity", 
                                                      "Progenitor_CD8", "Terminal_exhaustion", "Treg_signature", 
                                                      "M1_macrophage", "M2_macrophage", "IFN_gamma", "TGF_beta", 
                                                      "mregDC", "Immune_checkpoint"]]
    if len(sig_cols) > 0:
        df = adata.obs[sig_cols + ["group"]]
        df_grouped = df.groupby("group").mean().T
        plt.figure(figsize=(6, 8))
        sns.heatmap(df_grouped, cmap="viridis", annot=True, fmt=".2f", linewidths=0.5)
        plt.title("Mean IO Signature Scores per Disease Group", fontweight="bold")
        plt.tight_layout()
        plt.savefig(OUT_02B / "figures" / "eda_signature_mean_heatmap.png", dpi=200, bbox_inches="tight")
        plt.savefig(OUT_02B / "figures" / "eda_signature_mean_heatmap.pdf", dpi=200, bbox_inches="tight")
        plt.close()

def plot_bulk_deconv_eda():
    ref_path = OUT_02C / "tables" / "cibersortx_hcc_reference.tsv"
    if not ref_path.exists():
        return
    print("Generating EDA for 02c Bulk Deconvolution Reference...")
    df = pd.read_csv(ref_path, sep="\t", index_col=0)
    
    # Plot top variable genes in the reference
    var = df.var(axis=1).sort_values(ascending=False).head(50)
    df_var = df.loc[var.index]
    
    # z-score for heatmap
    df_z = (df_var - df_var.mean(axis=1).values.reshape(-1,1)) / df_var.std(axis=1).values.reshape(-1,1)
    
    plt.figure(figsize=(10, 10))
    # Clip extreme values for clear visualization
    sns.heatmap(df_z.clip(-2, 2), cmap="magma", yticklabels=True, xticklabels=True)
    plt.title("Top 50 Most Variable Genes in Immune Reference (Z-score)", fontweight="bold")
    plt.tight_layout()
    plt.savefig(OUT_02C / "figures" / "eda_highly_variable_genes_heatmap.png", dpi=200, bbox_inches="tight")
    plt.savefig(OUT_02C / "figures" / "eda_highly_variable_genes_heatmap.pdf", dpi=200, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    plot_deg_heatmap()
    plot_spatial_biomarkers_eda()
    plot_bulk_deconv_eda()
    print("EDA plots generation complete!")
