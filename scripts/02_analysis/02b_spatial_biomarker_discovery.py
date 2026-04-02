#!/usr/bin/env python3
"""
02b_spatial_biomarker_discovery.py — Prompt 3: Spatial Biomarker Discovery

1. Apply spatial immunotyping to all available spatial HCC data
   (GSE173905 Visium + any HCC_omics Visium data)
2. Extract spatial features: EIES, TLS, CD8-tumor distance, border regularity
3. Correlate spatial features with response (where data available)
4. Build LASSO model: spatial features → response prediction
5. Compare with bulk-derived biomarkers (TIDE, TMB, PD-L1)

Figures saved to: out/02b_spatial_biomarkers/figures/
Tables saved to:  out/02b_spatial_biomarkers/tables/
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
OUT     = PROJECT / "out" / "02b_spatial_biomarkers"
FIG     = OUT / "figures"
TAB     = OUT / "tables"
for d in [FIG, TAB]:
    d.mkdir(parents=True, exist_ok=True)

# ── IO-relevant gene signatures ────────────────────────────────────────
IO_SIGNATURES = {
    "T_cell_exhaustion":    ["PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","TOX","ENTPD1"],
    "T_cell_activation":    ["IFNG","GZMB","PRF1","TNF","CD69","IL2"],
    "T_cell_cytotoxicity":  ["GZMA","GZMB","GZMH","GZMK","PRF1","GNLY","NKG7"],
    "Treg_signature":       ["FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"],
    "M1_macrophage":        ["NOS2","IL1B","IL6","TNF","CD80","CD86"],
    "M2_macrophage":        ["CD163","MRC1","MSR1","CD200R1","IL10"],
    "IFN_gamma":            ["IFNG","STAT1","IRF1","CXCL9","CXCL10","CXCL11","IDO1"],
    "TGF_beta":             ["TGFB1","TGFB2","TGFB3","ACTA2","COL1A1","SERPINE1"],
    "Angiogenesis":         ["VEGFA","PECAM1","KDR","FLT1","ANGPT2","FGF2"],
    "Antigen_presentation": ["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","PSMB9"],
    "Immune_checkpoint":    ["CD274","PDCD1LG2","PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","VSIR"],
    "WNT_beta_catenin":     ["CTNNB1","APC","AXIN2","LGR5","MYC"],
    "mregDC":               ["CCR7","LAMP3","CCL19","CCL22","CD80","CD83","FSCN1","MARCKSL1"],
    "TLS_signature":        ["CXCL13","CXCR5","CR2","CD79A","MS4A1","CD19","IGHG1"],
}

CELL_TYPE_MARKERS = {
    "CD8_T":        ["CD8A","CD8B","GZMB","PRF1"],
    "CD4_T":        ["CD4","IL7R","CCR7"],
    "Treg":         ["FOXP3","IL2RA","CTLA4"],
    "NK":           ["NKG7","GNLY","KLRD1"],
    "B_cell":       ["CD79A","MS4A1","CD19"],
    "Macrophage":   ["CD68","CD163","C1QA","CSF1R"],
    "DC":           ["CLEC9A","XCR1","BATF3"],
    "Hepatocyte":   ["ALB","APOB","HP","TF"],
    "Endothelial":  ["PECAM1","CDH5","VWF"],
    "Fibroblast":   ["COL1A1","COL3A1","DCN","LUM"],
    "Tumor":        ["GPC3","AFP","GOLM1","EPCAM"],
}


def load_gse173905():
    """Load GSE173905 (Wu et al. 2021) HCC bulk RNA-seq as pseudo-bulk samples.

    Each file contains gene-level counts & FPKM for one sample.
    Columns: gene_id, <sample>_count, <sample>_fpkm, gene_name, gene_chr, ...
    Groups: NC (non-cancer), ER (early recurrence), RS (recurrence-sensitive).
    We load FPKM values and use gene_name as var index.
    """
    print("[GSE173905] Loading Wu 2021 HCC pseudo-bulk RNA-seq data...")
    d = DATA / "spatial_data" / "GSE173905_Wu2021_Visium"
    gz_files = sorted(d.glob("*.xls.gz"))
    if not gz_files:
        print("  [SKIP] No .xls.gz files found")
        return []

    datasets = []
    for gz in gz_files:
        sample_name = gz.stem.replace("fpkm.xls", "").rstrip(".")
        try:
            df = pd.read_csv(gz, sep="\t")
            # Identify count and fpkm columns
            count_col = [c for c in df.columns if c.endswith("_count")][0]
            fpkm_col = [c for c in df.columns if c.endswith("_fpkm")][0]
            # Use gene_name as identifier, drop duplicates
            df = df.dropna(subset=["gene_name"])
            df = df.drop_duplicates(subset="gene_name", keep="first")
            gene_names = df["gene_name"].values
            counts = df[count_col].values.astype(np.float64)

            # Create AnnData with 1 "spot" = 1 pseudo-bulk sample
            adata = sc.AnnData(
                X=np.array([counts]),
                obs=pd.DataFrame(index=[sample_name]),
                var=pd.DataFrame(index=gene_names),
            )
            adata.var_names_make_unique()
            adata.obs["sample"] = sample_name
            adata.obs["dataset"] = "GSE173905"
            # Assign group from sample name prefix
            if "NC" in sample_name:
                adata.obs["group"] = "Non_Cancer"
            elif "ER" in sample_name:
                adata.obs["group"] = "Early_Recurrence"
            elif "RS" in sample_name:
                adata.obs["group"] = "Recurrence_Sensitive"
            datasets.append(adata)
            print(f"  {sample_name}: 1 pseudo-bulk × {adata.n_vars} genes")
        except Exception as e:
            print(f"  [WARN] Failed to load {gz.name}: {e}")
    return datasets


def load_hcc_omics_visium():
    """Try to load HCC_omics Visium data."""
    print("[HCC_omics] Checking for Visium data...")
    hcc_dir = Path("/mnt/projects/HCC_omics")
    c2l_dirs = sorted(hcc_dir.glob("out/hcc_spata2_tma/*/cell2location/visium_counts.h5ad"))

    datasets = []
    for h5_path in c2l_dirs[:5]:
        sample = h5_path.parent.parent.name
        try:
            adata = sc.read_h5ad(str(h5_path))
            adata.obs["sample"] = sample
            adata.obs["dataset"] = "HCC_omics"
            datasets.append(adata)
            print(f"  {sample}: {adata.n_obs} spots × {adata.n_vars} genes")
        except Exception as e:
            print(f"  [SKIP] {sample}: {e}")
    return datasets


def score_spatial_signatures(adata):
    """Score IO signatures on spatial data."""
    available = set(adata.var_names)
    scored = []
    for name, genes in IO_SIGNATURES.items():
        present = [g for g in genes if g in available]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=name, use_raw=False)
            scored.append(name)
    print(f"  Scored {len(scored)}/{len(IO_SIGNATURES)} signatures")
    return scored


def classify_spatial_immunotype(adata):
    """Classify spots into spatial immunotypes."""
    has_act  = "T_cell_activation" in adata.obs.columns
    has_tgfb = "TGF_beta" in adata.obs.columns

    if has_act and has_tgfb:
        t_act = adata.obs["T_cell_activation"].values
        tgfb  = adata.obs["TGF_beta"].values
        act_hi = np.percentile(t_act, 66)
        act_lo = np.percentile(t_act, 33)
        tgfb_hi = np.percentile(tgfb, 66)

        labels = []
        for a, t in zip(t_act, tgfb):
            if a > act_hi and t < tgfb_hi:
                labels.append("Immune_Hot")
            elif a > act_hi and t >= tgfb_hi:
                labels.append("Immune_Excluded")
            elif a < act_lo:
                labels.append("Immune_Desert")
            else:
                labels.append("Immune_Mixed")
        adata.obs["spatial_immunotype"] = labels
        print(f"  Immunotypes: {dict(Counter(labels))}")
    return adata


def compute_spatial_features(adata, sample_name):
    """Extract spatial features per sample: EIES, TLS proxy, etc."""
    features = {"sample": sample_name}
    sig_cols = [c for c in adata.obs.columns if c in IO_SIGNATURES]

    for sig in sig_cols:
        features[f"mean_{sig}"] = adata.obs[sig].mean()
        features[f"std_{sig}"]  = adata.obs[sig].std()

    # TLS proxy score
    if "TLS_signature" in sig_cols:
        features["TLS_score"] = adata.obs["TLS_signature"].mean()
        features["TLS_high_frac"] = (adata.obs["TLS_signature"] > adata.obs["TLS_signature"].quantile(0.75)).mean()

    # EIES (Effective Immune-Excluded Score): ratio of border vs interior immune
    if "T_cell_activation" in sig_cols and "TGF_beta" in sig_cols:
        features["EIES"] = adata.obs["TGF_beta"].mean() / (adata.obs["T_cell_activation"].mean() + 1e-6)

    # Immunotype fractions
    if "spatial_immunotype" in adata.obs.columns:
        imm_counts = adata.obs["spatial_immunotype"].value_counts(normalize=True)
        for imm_type, frac in imm_counts.items():
            features[f"frac_{imm_type}"] = frac

    return features


def plot_spatial_overview(adata, scored, sample_name):
    """Generate spatial overview plots per sample."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # Signature heatmap
    ax = axes[0, 0]
    sig_cols = [s for s in scored if s in adata.obs.columns]
    if sig_cols:
        sig_corr = adata.obs[sig_cols].corr()
        sns.heatmap(sig_corr, cmap="RdBu_r", center=0, ax=ax, square=True,
                    xticklabels=True, yticklabels=True, cbar_kws={"shrink": 0.8})
        ax.set_title("A. Signature Correlations")
        ax.tick_params(labelsize=6)

    # Immunotype pie
    ax = axes[0, 1]
    if "spatial_immunotype" in adata.obs.columns:
        imm_counts = adata.obs["spatial_immunotype"].value_counts()
        colors_map = {"Immune_Hot": "#e74c3c", "Immune_Excluded": "#f39c12",
                      "Immune_Desert": "#3498db", "Immune_Mixed": "#95a5a6"}
        colors = [colors_map.get(c, "#666") for c in imm_counts.index]
        ax.pie(imm_counts.values, labels=imm_counts.index, colors=colors,
               autopct="%1.1f%%", startangle=140)
        ax.set_title("B. Spatial Immunotype Distribution")

    # key signature boxplots by immunotype
    ax = axes[0, 2]
    key_sigs = ["T_cell_activation", "T_cell_exhaustion", "IFN_gamma", "TGF_beta"]
    key_sigs = [s for s in key_sigs if s in sig_cols]
    if key_sigs and "spatial_immunotype" in adata.obs.columns:
        melt = adata.obs[key_sigs + ["spatial_immunotype"]].melt(
            id_vars="spatial_immunotype", var_name="Signature", value_name="Score")
        if len(melt) > 10000:
            melt = melt.sample(10000, random_state=42)
        sns.boxplot(data=melt, x="Signature", y="Score", hue="spatial_immunotype",
                    ax=ax, palette="Set2", fliersize=0)
        ax.set_title("C. Signatures by Immunotype")
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6)

    # Exhaustion vs Activation
    ax = axes[1, 0]
    if "T_cell_activation" in sig_cols and "T_cell_exhaustion" in sig_cols:
        ax.scatter(adata.obs["T_cell_activation"], adata.obs["T_cell_exhaustion"],
                   alpha=0.2, s=3, rasterized=True)
        r, p = spearmanr(adata.obs["T_cell_activation"], adata.obs["T_cell_exhaustion"])
        ax.set_xlabel("Activation"); ax.set_ylabel("Exhaustion")
        ax.set_title(f"D. Act vs Exh (ρ={r:.3f})")

    # Checkpoint violin by immunotype
    ax = axes[1, 1]
    if "Immune_checkpoint" in sig_cols and "spatial_immunotype" in adata.obs.columns:
        sns.violinplot(data=adata.obs, x="spatial_immunotype", y="Immune_checkpoint",
                       ax=ax, palette="Set2", inner="quartile")
        ax.set_title("E. Checkpoint by Immunotype")
        ax.tick_params(axis="x", rotation=30, labelsize=8)

    # TLS score histogram
    ax = axes[1, 2]
    if "TLS_signature" in sig_cols:
        ax.hist(adata.obs["TLS_signature"], bins=50, color="#2ecc71", alpha=0.7, edgecolor="white")
        ax.axvline(adata.obs["TLS_signature"].quantile(0.75), color="red", linestyle="--",
                   label="Q75")
        ax.set_xlabel("TLS Score"); ax.set_ylabel("Count")
        ax.set_title("F. TLS Signature Distribution")
        ax.legend()

    plt.suptitle(f"Spatial Biomarker Landscape — {sample_name}", fontweight="bold", fontsize=13)
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FIG / f"02b_spatial_overview_{sample_name}.{ext}", dpi=200, bbox_inches="tight")
    plt.close()


def plot_sample_comparison(features_df, scored):
    """Compare spatial features across samples."""
    if len(features_df) < 2:
        return

    # Immunotype fractions across samples
    frac_cols = [c for c in features_df.columns if c.startswith("frac_")]
    if frac_cols:
        fig, ax = plt.subplots(figsize=(10, 6))
        features_df.set_index("sample")[frac_cols].plot(kind="bar", stacked=True, ax=ax,
                                                         colormap="Set2")
        ax.set_ylabel("Fraction")
        ax.set_title("Spatial Immunotype Composition by Sample", fontweight="bold")
        ax.legend(loc="upper right", fontsize=8)
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"02b_immunotype_comparison.{ext}", dpi=200, bbox_inches="tight")
        plt.close()

    # Signature heatmap across samples
    mean_cols = [c for c in features_df.columns if c.startswith("mean_")]
    if mean_cols:
        fig, ax = plt.subplots(figsize=(12, max(4, len(features_df) * 0.5)))
        plot_df = features_df.set_index("sample")[mean_cols]
        plot_df.columns = [c.replace("mean_", "") for c in plot_df.columns]
        sns.heatmap(plot_df, cmap="RdBu_r", center=0, ax=ax, annot=True, fmt=".2f",
                    linewidths=0.5)
        ax.set_title("Mean IO Signatures Across Samples", fontweight="bold")
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FIG / f"02b_signature_heatmap_samples.{ext}", dpi=200, bbox_inches="tight")
        plt.close()


def main():
    print("=" * 70)
    print("Prompt 3: Spatial Biomarker Discovery")
    print("=" * 70)

    # ── Load spatial / pseudo-bulk datasets ────────────────────────
    gse_datasets = load_gse173905()
    hcc_datasets = load_hcc_omics_visium()

    has_spatial = len(hcc_datasets) > 0  # true Visium data from HCC_omics
    has_pseudobulk = len(gse_datasets) > 0  # GSE173905 pseudo-bulk

    if not has_spatial and not has_pseudobulk:
        print("[ERROR] No spatial or pseudo-bulk datasets available.")
        return

    # ── Part A: True spatial analysis (HCC_omics Visium) ───────────
    if has_spatial:
        print(f"\n  True spatial datasets: {len(hcc_datasets)}")
        all_features = []
        for adata in hcc_datasets:
            sample_name = adata.obs["sample"].iloc[0]
            print(f"\n--- Analyzing {sample_name} ---")

            if adata.X.max() > 100:
                adata.layers["counts"] = adata.X.copy()
                sc.pp.normalize_total(adata, target_sum=1e4)
                sc.pp.log1p(adata)

            n_top = min(2000, adata.n_vars)
            if adata.n_obs >= 10:
                sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="seurat_v3",
                                            span=1.0 if adata.n_obs < 5000 else 0.3)
            n_comps = min(30, adata.n_obs - 1, adata.n_vars - 1)
            if n_comps >= 2:
                sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=True)
                sc.pp.neighbors(adata, n_pcs=min(20, n_comps))
                sc.tl.umap(adata)

            scored = score_spatial_signatures(adata)
            adata = classify_spatial_immunotype(adata)
            features = compute_spatial_features(adata, sample_name)
            all_features.append(features)
            plot_spatial_overview(adata, scored, sample_name)

        if all_features:
            features_df = pd.DataFrame(all_features)
            features_df.to_csv(TAB / "spatial_features_per_sample.csv", index=False)
            print(f"\n  Spatial features table: {features_df.shape}")
            plot_sample_comparison(features_df, scored)

    # ── Part B: Pseudo-bulk analysis (GSE173905) ───────────────────
    if has_pseudobulk:
        print(f"\n  Pseudo-bulk datasets (GSE173905): {len(gse_datasets)}")

        # Concatenate into one AnnData
        combined = sc.concat(gse_datasets, join="inner")
        combined.obs_names_make_unique()
        print(f"  Combined: {combined.n_obs} samples × {combined.n_vars} genes")

        # Normalize
        if combined.X.max() > 100:
            combined.layers["counts"] = combined.X.copy()
            sc.pp.normalize_total(combined, target_sum=1e4)
            sc.pp.log1p(combined)

        # Score IO signatures
        scored = score_spatial_signatures(combined)

        # Build group comparison figure
        sig_cols = [s for s in scored if s in combined.obs.columns]
        if sig_cols and "group" in combined.obs.columns:
            groups = combined.obs["group"].unique()
            print(f"\n  Groups: {list(groups)}")

            # Fig 1: Heatmap of IO signatures across samples
            fig, axes = plt.subplots(1, 2, figsize=(16, 6))

            ax = axes[0]
            plot_df = combined.obs.set_index(combined.obs.index)[sig_cols]
            plot_df["group"] = combined.obs["group"].values
            plot_df = plot_df.sort_values("group")
            plot_mat = plot_df[sig_cols]
            sns.heatmap(plot_mat.T, cmap="RdBu_r", center=0, ax=ax, annot=True, fmt=".2f",
                        xticklabels=plot_df.index, yticklabels=[s.replace("_", " ") for s in sig_cols],
                        linewidths=0.5)
            ax.set_title("A. IO Signature Scores by Sample", fontweight="bold")
            ax.tick_params(labelsize=7)

            # Fig 2: Group means barplot
            ax = axes[1]
            group_means = combined.obs.groupby("group")[sig_cols].mean().T
            group_means.columns = [c.replace("_", " ") for c in group_means.columns]
            group_means.index = [s.replace("_", " ") for s in group_means.index]
            group_means.plot(kind="barh", ax=ax, colormap="Set2")
            ax.set_title("B. Mean IO Signatures by Group", fontweight="bold")
            ax.set_xlabel("Signature Score")
            ax.legend(fontsize=8)

            plt.suptitle("GSE173905: HCC Pseudo-Bulk IO Signature Analysis", fontweight="bold", fontsize=13)
            plt.tight_layout()
            for ext in ["png", "pdf"]:
                plt.savefig(FIG / f"02b_pseudobulk_io_signatures.{ext}", dpi=200, bbox_inches="tight")
            plt.close()

            # Group comparison stats (Mann-Whitney where possible)
            stat_results = []
            if len(groups) >= 2:
                # Compare Early_Recurrence vs Recurrence_Sensitive
                g1, g2 = "Early_Recurrence", "Recurrence_Sensitive"
                mask1 = combined.obs["group"] == g1
                mask2 = combined.obs["group"] == g2
                if mask1.sum() >= 2 and mask2.sum() >= 2:
                    for sig in sig_cols:
                        u, p = mannwhitneyu(
                            combined.obs.loc[mask1, sig].values,
                            combined.obs.loc[mask2, sig].values,
                            alternative="two-sided"
                        )
                        stat_results.append({
                            "signature": sig, "group1": g1, "group2": g2,
                            "mean_g1": combined.obs.loc[mask1, sig].mean(),
                            "mean_g2": combined.obs.loc[mask2, sig].mean(),
                            "U_stat": u, "p_value": p
                        })

            if stat_results:
                stats_df = pd.DataFrame(stat_results).sort_values("p_value")
                stats_df.to_csv(TAB / "pseudobulk_group_comparison.csv", index=False)
                print(f"\n  Group comparison stats saved ({len(stats_df)} tests)")

            # ------------- ADD UMAP / PCA ---------------
            print("\n[UMAP] Generating Pseudo-Bulk Sample Embeddings...")
            if len(combined) >= 4:
                sc.pp.pca(combined, n_comps=min(15, len(combined)-1))
                sc.pp.neighbors(combined, use_rep="X_pca", n_neighbors=min(15, len(combined)-1))
                sc.tl.umap(combined)
                
                fig, axes = plt.subplots(1, 2, figsize=(11, 5))
                sc.pl.pca(combined, color="group", size=200, ax=axes[0], show=False)
                axes[0].set_title("A. PCA by Group", fontweight="bold")
                
                sc.pl.umap(combined, color="group", size=200, ax=axes[1], show=False)
                axes[1].set_title("B. UMAP by Group", fontweight="bold")
                
                plt.suptitle("Sample Embeddings: GSE173905 Pseudo-Bulk", fontweight="bold", fontsize=12)
                plt.tight_layout()
                for ext in ["png", "pdf"]:
                    plt.savefig(FIG / f"02b_pseudobulk_embeddings.{ext}", dpi=200, bbox_inches="tight")
                plt.close()
            # --------------------------------------------

        # Save scored pseudo-bulk data
        combined.write(TAB / "GSE173905_pseudobulk_scored.h5ad")
        print(f"  Scored pseudo-bulk saved to {TAB / 'GSE173905_pseudobulk_scored.h5ad'}")

    print(f"\n{'='*70}")
    print("[DONE] Spatial Biomarker Discovery")
    print(f"  Figures: {FIG}")
    print(f"  Tables:  {TAB}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
