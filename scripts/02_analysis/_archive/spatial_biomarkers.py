#!/usr/bin/env python3
"""
Spatial Immunotherapy Response — IO Biomarker Discovery

Analyzes the spatial and single-cell transcriptomic landscape of HCC patients
treated with anti-PD1/anti-PDL1 immunotherapy to identify:
  1. Spatial biomarkers of response/resistance
  2. TME remodeling signatures
  3. Predictive immune features

Uses HCC_omics data + downloaded IO datasets.
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

PROJECT = Path(__file__).resolve().parents[2]
DATA = PROJECT / "data"
OUT = PROJECT / "out" / "analysis"
FIG_DIR = PROJECT / "out" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

HCC_OMICS = Path.home() / "projects" / "HCC_omics"

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

# Cell type markers for annotation
MARKERS = {
    "Hepatocyte":   ["ALB", "APOA1", "TF", "HP"],
    "T cell":       ["CD3D", "CD3E", "CD2"],
    "CD8 T":        ["CD8A", "CD8B", "GZMB"],
    "Treg":         ["FOXP3", "IL2RA", "CTLA4"],
    "NK":           ["NCAM1", "GNLY", "KLRD1"],
    "Macrophage":   ["CD68", "CD163", "C1QA"],
    "B cell":       ["CD79A", "MS4A1"],
    "DC":           ["CD1C", "CLEC9A"],
    "Endothelial":  ["PECAM1", "VWF"],
    "Fibroblast":   ["COL1A1", "DCN", "LUM"],
    "Tumor":        ["GPC3", "AFP", "GOLM1"],
}


def load_hcc_omics_data():
    """Try to load processed data from HCC_omics project."""
    # Try Cell2Location outputs
    c2l_dirs = sorted(HCC_OMICS.glob("out/hcc_spata2_tma/*/cell2location/visium_counts.h5ad"))
    
    datasets = []
    for h5_path in c2l_dirs[:5]:  # Take up to 5 samples
        sample = h5_path.parent.parent.name
        try:
            adata = sc.read_h5ad(str(h5_path))
            adata.obs["sample"] = sample
            datasets.append(adata)
            print(f"  Loaded HCC_omics sample: {sample} ({adata.n_obs} spots)")
        except Exception as e:
            print(f"  [SKIP] {sample}: {e}")
    
    return datasets


def create_io_analysis_dataset():
    """Create analysis dataset from available sources."""
    
    # First try HCC_omics Visium data
    print("[1] Attempting to load HCC_omics Visium data...")
    datasets = load_hcc_omics_data()
    
    if datasets:
        adata = sc.concat(datasets, join="outer", fill_value=0)
        adata.var_names_make_unique()
        adata.obs_names_make_unique()
        print(f"  Combined: {adata.n_obs} spots × {adata.n_vars} genes")
        return adata, "visium"
    
    # Fallback: use HCC TME evolution data
    tme_path = Path.home() / "projects" / "hcc_tme_evolution" / "data" / "processed" / "integrated_tme_annotated.h5ad"
    if tme_path.exists():
        print("[1] Loading HCC TME Evolution integrated data...")
        adata = sc.read_h5ad(str(tme_path))
        return adata, "scrna"
    
    print("[ERROR] No suitable data found")
    return None, None


def score_io_signatures(adata):
    """Score all IO-relevant gene signatures."""
    print("\n[2] Scoring IO signatures...")
    available = set(adata.var_names)
    
    scored = {}
    for sig_name, genes in IO_SIGNATURES.items():
        present = [g for g in genes if g in available]
        if len(present) >= 2:
            sc.tl.score_genes(adata, gene_list=present, score_name=sig_name)
            scored[sig_name] = len(present)
            print(f"  {sig_name}: {len(present)}/{len(genes)} genes")
    
    return scored


def classify_io_phenotype(adata):
    """Classify spots/cells into IO-relevant phenotypes."""
    print("\n[3] Classifying IO phenotypes...")
    
    # Determine if immune-hot or immune-cold based on T cell signatures
    has_tcell_act = "T_cell_activation" in adata.obs.columns
    has_ifng = "IFN_gamma" in adata.obs.columns
    has_tgfb = "TGF_beta" in adata.obs.columns
    has_exhaust = "T_cell_exhaustion" in adata.obs.columns
    
    if has_tcell_act and has_tgfb:
        # Classify into IO-relevant categories
        t_act = adata.obs["T_cell_activation"].values
        tgfb = adata.obs["TGF_beta"].values
        
        t_act_thresh = np.percentile(t_act, 66)
        tgfb_thresh = np.percentile(tgfb, 66)
        
        conditions = []
        for i in range(len(adata)):
            if t_act[i] > t_act_thresh and tgfb[i] < tgfb_thresh:
                conditions.append("Immune_Hot_Favorable")
            elif t_act[i] > t_act_thresh and tgfb[i] > tgfb_thresh:
                conditions.append("Immune_Excluded")
            elif t_act[i] < np.percentile(t_act, 33):
                conditions.append("Immune_Desert")
            else:
                conditions.append("Immune_Mixed")
        
        adata.obs["io_phenotype"] = conditions
        print(f"  IO phenotypes: {dict(Counter(conditions))}")
    
    return adata


def differential_signature_analysis(adata):
    """Compare IO signatures across phenotypes."""
    print("\n[4] Differential signature analysis...")
    
    results = []
    sig_cols = [c for c in adata.obs.columns if c in IO_SIGNATURES]
    
    if "io_phenotype" not in adata.obs.columns:
        return pd.DataFrame()
    
    phenotypes = adata.obs["io_phenotype"].unique()
    
    for sig in sig_cols:
        for i, p1 in enumerate(phenotypes):
            for p2 in phenotypes[i+1:]:
                vals1 = adata.obs.loc[adata.obs["io_phenotype"] == p1, sig].values
                vals2 = adata.obs.loc[adata.obs["io_phenotype"] == p2, sig].values
                
                if len(vals1) < 10 or len(vals2) < 10:
                    continue
                
                stat, pval = mannwhitneyu(vals1, vals2, alternative="two-sided")
                fc = np.mean(vals1) - np.mean(vals2)
                
                results.append({
                    "signature": sig,
                    "group1": p1,
                    "group2": p2,
                    "mean_diff": round(fc, 4),
                    "pvalue": pval,
                    "n1": len(vals1),
                    "n2": len(vals2),
                })
    
    if results:
        diff_df = pd.DataFrame(results)
        diff_df["padj"] = np.minimum(diff_df["pvalue"] * len(diff_df), 1.0)
        diff_df = diff_df.sort_values("pvalue")
        return diff_df
    return pd.DataFrame()


def plot_io_landscape(adata, diff_df):
    """Generate IO biomarker landscape figures."""
    print("\n[5] Generating figures...")
    
    sig_cols = [c for c in adata.obs.columns if c in IO_SIGNATURES]
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # 1. Signature correlation heatmap
    ax = axes[0, 0]
    if sig_cols:
        sig_corr = adata.obs[sig_cols].corr()
        sns.heatmap(sig_corr, cmap="RdBu_r", center=0, ax=ax, 
                    cbar_kws={"label": "Pearson r"}, square=True,
                    xticklabels=True, yticklabels=True)
        ax.set_title("IO Signature Correlations")
        ax.tick_params(labelsize=7)
    
    # 2. IO phenotype distribution
    ax = axes[0, 1]
    if "io_phenotype" in adata.obs.columns:
        pheno_counts = adata.obs["io_phenotype"].value_counts()
        colors = {"Immune_Hot_Favorable": "#e74c3c", "Immune_Excluded": "#f39c12",
                  "Immune_Desert": "#3498db", "Immune_Mixed": "#95a5a6"}
        pheno_counts.plot(kind="bar", ax=ax, 
                         color=[colors.get(p, "#666") for p in pheno_counts.index])
        ax.set_title("IO Phenotype Distribution")
        ax.set_ylabel("Count")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    
    # 3. Key signatures by phenotype
    ax = axes[0, 2]
    key_sigs = ["T_cell_activation", "T_cell_exhaustion", "IFN_gamma", "TGF_beta"]
    key_sigs = [s for s in key_sigs if s in sig_cols]
    if key_sigs and "io_phenotype" in adata.obs.columns:
        melt_data = []
        for sig in key_sigs:
            for _, row in adata.obs[["io_phenotype", sig]].iterrows():
                melt_data.append({"Phenotype": row["io_phenotype"], "Score": row[sig], "Signature": sig})
        melt_df = pd.DataFrame(melt_data)
        # Subsample for plotting
        if len(melt_df) > 10000:
            melt_df = melt_df.sample(10000, random_state=42)
        sns.boxplot(data=melt_df, x="Signature", y="Score", hue="Phenotype", ax=ax, 
                   palette="Set2", fliersize=0)
        ax.set_title("Signatures by IO Phenotype")
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6, title_fontsize=7)
    
    # 4. Exhaustion vs Activation scatter
    ax = axes[1, 0]
    if "T_cell_activation" in sig_cols and "T_cell_exhaustion" in sig_cols:
        ax.scatter(adata.obs["T_cell_activation"], adata.obs["T_cell_exhaustion"],
                   alpha=0.1, s=1, rasterized=True)
        r, p = spearmanr(adata.obs["T_cell_activation"], adata.obs["T_cell_exhaustion"])
        ax.set_xlabel("T Cell Activation Score")
        ax.set_ylabel("T Cell Exhaustion Score")
        ax.set_title(f"Activation vs Exhaustion (ρ={r:.3f})")
    
    # 5. Checkpoint ligand scores
    ax = axes[1, 1]
    if "Immune_checkpoint" in sig_cols:
        if "io_phenotype" in adata.obs.columns:
            sns.violinplot(data=adata.obs, x="io_phenotype", y="Immune_checkpoint",
                          ax=ax, palette="Set2", inner="quartile")
            ax.set_title("Immune Checkpoint Score")
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    
    # 6. Top differential signatures
    ax = axes[1, 2]
    if len(diff_df) > 0:
        top = diff_df.head(10)
        ax.barh(range(len(top)), top["mean_diff"], color=["#e74c3c" if d > 0 else "#3498db" for d in top["mean_diff"]])
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels([f"{r['signature']}\n({r['group1']} vs {r['group2']})" for _, r in top.iterrows()], fontsize=7)
        ax.set_xlabel("Mean Difference")
        ax.set_title("Top Differential Signatures")
    
    plt.suptitle("Spatial Immunotherapy Response — IO Biomarker Landscape", fontsize=14, y=1.02)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "io_biomarker_landscape.png", dpi=200)
    fig.savefig(FIG_DIR / "io_biomarker_landscape.pdf")
    print(f"  Saved: {FIG_DIR / 'io_biomarker_landscape.png'}")
    plt.close()


def main():
    print("=" * 70)
    print("Spatial Immunotherapy Response — IO Biomarker Discovery")
    print("=" * 70)
    
    # Load data
    adata, data_type = create_io_analysis_dataset()
    if adata is None:
        return
    
    print(f"\nData type: {data_type}")
    print(f"Shape: {adata.shape}")
    
    # Ensure normalized
    if adata.X.max() > 100:
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    
    # Score IO signatures
    scored = score_io_signatures(adata)
    
    # Classify IO phenotypes
    adata = classify_io_phenotype(adata)
    
    # Differential analysis
    diff_df = differential_signature_analysis(adata)
    if len(diff_df) > 0:
        diff_df.to_csv(OUT / "differential_signatures.tsv", sep="\t", index=False)
        print(f"\nTop differential signatures:")
        print(diff_df.head(10).to_string(index=False))
    
    # Generate figures
    plot_io_landscape(adata, diff_df)
    
    # Save results
    sig_summary = adata.obs[[c for c in adata.obs.columns if c in IO_SIGNATURES or c == "io_phenotype"]]
    sig_summary.to_csv(OUT / "io_signature_scores.tsv", sep="\t")
    
    print(f"\n[DONE] IO Biomarker analysis complete.")
    print(f"  Results: {OUT}")
    print(f"  Figures: {FIG_DIR}")


if __name__ == "__main__":
    main()
