#!/usr/bin/env python3
"""
02b_visium_downstream.py — Downstream Spatial Analyses for HCC IO Visium

Extends 02b_spatial_visium_analysis.py with:
  1. Sample-level permutation tests (honest n=4 vs n=3)
  2. Leave-one-out sensitivity analysis
  3. Spatially variable genes (SVGs) per section
  4. Spatial co-occurrence of cell-type marker modules
  5. Spatial clustering (Leiden) per section with R vs NR comparison
  6. Immune infiltration gradient analysis (distance to tumour core)
  7. Spatial entropy / heterogeneity per section
  8. Ligand–receptor co-localization (top L-R pairs)

Outputs → out/02b_spatial_biomarkers/figures/ and tables/
"""

import warnings, sys
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from pathlib import Path
from scipy.stats import mannwhitneyu, permutation_test, spearmanr
from scipy.sparse import issparse
from scipy.spatial.distance import cdist
from collections import OrderedDict
from itertools import combinations
import gc

sc.settings.verbosity = 1

# ── Paths ──────────────────────────────────────────────────────────────
PROJECT = Path(__file__).resolve().parents[2]
DATA    = PROJECT / "data" / "spatial_data"
OUT     = PROJECT / "out" / "02b_spatial_biomarkers"
FIG     = OUT / "figures"
TAB     = OUT / "tables"
for d in [FIG, TAB]:
    d.mkdir(parents=True, exist_ok=True)

# ── IO Signatures (identical to main script) ───────────────────────────
IO_SIGS = OrderedDict({
    "T_cell_exhaustion":    ["PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","TOX","ENTPD1"],
    "T_cell_activation":    ["IFNG","GZMB","PRF1","TNF","CD69","IL2"],
    "T_cell_cytotoxicity":  ["GZMA","GZMB","GZMH","GZMK","PRF1","GNLY","NKG7"],
    "Progenitor_CD8":       ["TCF7","SLAMF6","LEF1","CCR7","IL7R","SELL"],
    "Terminal_exhaustion":  ["TOX","ENTPD1","HAVCR2","LAYN","TIGIT","CXCL13"],
    "Treg_signature":       ["FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"],
    "M2_macrophage":        ["CD163","MRC1","MSR1","CD200R1","IL10"],
    "IFN_gamma":            ["IFNG","STAT1","IRF1","CXCL9","CXCL10","CXCL11","IDO1"],
    "TGF_beta":             ["TGFB1","TGFB2","TGFB3","ACTA2","COL1A1","SERPINE1"],
    "Antigen_presentation": ["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","PSMB9"],
    "Immune_checkpoint":    ["CD274","PDCD1LG2","PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","VSIR"],
    "Angiogenesis":         ["VEGFA","PECAM1","KDR","FLT1","ANGPT2","FGF2"],
    "WNT_beta_catenin":     ["CTNNB1","APC","AXIN2","LGR5","MYC"],
    "mregDC":               ["CCR7","LAMP3","CCL19","CCL22","CD80","CD83","FSCN1","MARCKSL1"],
    "TLS_signature":        ["CXCL13","CXCR5","CR2","CD79A","MS4A1","CD19","IGHG1"],
})

CELL_MARKERS = {
    "CD8_T":       ["CD8A","CD8B","GZMB","PRF1"],
    "CD4_T":       ["CD4","IL7R","CCR7"],
    "Treg":        ["FOXP3","IL2RA","CTLA4"],
    "NK":          ["NKG7","GNLY","KLRD1"],
    "B_cell":      ["CD79A","MS4A1","CD19"],
    "Macrophage":  ["CD68","CD163","C1QA","CSF1R"],
    "DC":          ["CLEC9A","XCR1","BATF3"],
    "Hepatocyte":  ["ALB","APOB","HP","TF"],
    "Endothelial": ["PECAM1","CDH5","VWF"],
    "Fibroblast":  ["COL1A1","COL3A1","DCN","LUM"],
    "Tumor":       ["GPC3","AFP","EPCAM"],
}

# Key L-R pairs relevant to IO and TME
LR_PAIRS = [
    ("CD274", "PDCD1"),     # PD-L1 / PD-1
    ("PDCD1LG2", "PDCD1"),  # PD-L2 / PD-1
    ("CD80", "CTLA4"),      # B7-1 / CTLA-4
    ("CD86", "CTLA4"),      # B7-2 / CTLA-4
    ("LGALS9", "HAVCR2"),   # Galectin-9 / TIM-3
    ("CXCL13", "CXCR5"),   # TLS axis
    ("CXCL9", "CXCR3"),    # IFNg-driven T recruitment
    ("CXCL10", "CXCR3"),
    ("CCL19", "CCR7"),      # mregDC / T cell migration
    ("TGFB1", "TGFBR2"),   # TGF-β axis
    ("VEGFA", "KDR"),       # Angiogenesis
    ("SPP1", "CD44"),       # Macrophage–tumor
    ("TNF", "TNFRSF1A"),   # TNF signalling
    ("IFNG", "IFNGR1"),    # IFNg signalling
    ("IL10", "IL10RA"),     # Immunosuppression
]

RESPONSE_PAL = {"Responder": "#2ecc71", "Non_Responder": "#e74c3c"}


# ═══════════════════════════════════════════════════════════════════════
#  DATA LOADING (reuse from main script)
# ═══════════════════════════════════════════════════════════════════════

def load_visium_h5ad(path, sample_name):
    adata = sc.read_h5ad(str(path))
    adata.var_names_make_unique()
    adata.obs_names = [f"{sample_name}_{b}" for b in adata.obs_names]
    adata.obsm["spatial"] = adata.obs[["pxl_col", "pxl_row"]].values.astype(float)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
    adata.obs["sample"] = sample_name
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def score_all(adata):
    avail = set(adata.var_names)
    scored_io = []
    for name, genes in IO_SIGS.items():
        present = [g for g in genes if g in avail]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=name, use_raw=False)
            scored_io.append(name)
    scored_ct = []
    for name, genes in CELL_MARKERS.items():
        present = [g for g in genes if g in avail]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=f"ct_{name}", use_raw=False)
            scored_ct.append(f"ct_{name}")
    return scored_io, scored_ct


def load_gse238264():
    d = DATA / "GSE238264_HCC_IO_Visium"
    samples = {
        "HCC1R": "Responder", "HCC2R": "Responder",
        "HCC3R": "Responder", "HCC4R": "Responder",
        "HCC5NR": "Non_Responder", "HCC6NR": "Non_Responder",
        "HCC7NR": "Non_Responder",
    }
    datasets = {}
    for name, resp in samples.items():
        path = d / f"{name}.h5ad"
        if not path.exists():
            continue
        print(f"  Loading {name} ({resp})...")
        adata = load_visium_h5ad(path, name)
        adata.obs["response"] = resp
        score_all(adata)
        sq.gr.spatial_neighbors(adata, n_neighs=6, coord_type="generic")
        datasets[name] = adata
        print(f"    {adata.n_obs} spots × {adata.n_vars} genes")
    return datasets


# ═══════════════════════════════════════════════════════════════════════
#  1. SAMPLE-LEVEL PERMUTATION TESTS  (honest n=4 vs n=3)
# ═══════════════════════════════════════════════════════════════════════

def sample_level_permutation_tests(datasets):
    """Mann-Whitney + permutation test on sample means (n=4 vs n=3)."""
    print("\n" + "=" * 70)
    print("[1] Sample-Level Permutation Tests (n=4 R vs n=3 NR)")
    print("=" * 70)

    scored = [s for s in IO_SIGS.keys() if s in list(datasets.values())[0].obs.columns]

    # Compute sample-level means
    rows = []
    for sname, adata in datasets.items():
        row = {"sample": sname, "response": adata.obs["response"].iloc[0]}
        for sig in scored:
            row[sig] = adata.obs[sig].mean()
        rows.append(row)
    sample_df = pd.DataFrame(rows)

    r_mask = sample_df["response"] == "Responder"

    results = []
    for sig in scored:
        r_vals = sample_df.loc[r_mask, sig].values
        nr_vals = sample_df.loc[~r_mask, sig].values

        # Mann-Whitney (n=4 vs 3)
        try:
            u, p_mw = mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
        except Exception:
            u, p_mw = np.nan, np.nan

        # Exact permutation test on mean difference
        def stat_func(x, y, axis):
            return np.mean(x, axis=axis) - np.mean(y, axis=axis)

        try:
            perm_res = permutation_test(
                (r_vals, nr_vals), stat_func, n_resamples=9999,
                alternative="two-sided"
            )
            p_perm = perm_res.pvalue
        except Exception:
            p_perm = np.nan

        effect = np.mean(r_vals) - np.mean(nr_vals)
        results.append({
            "signature": sig,
            "mean_R_sample": np.mean(r_vals),
            "mean_NR_sample": np.mean(nr_vals),
            "effect_size": effect,
            "U_stat": u,
            "p_MannWhitney_n7": p_mw,
            "p_permutation_9999": p_perm,
        })

    res_df = pd.DataFrame(results).sort_values("p_permutation_9999")
    res_df.to_csv(TAB / "gse238264_sample_level_permutation.csv", index=False)

    # ── Figure: sample-level dot plot with permutation p ───────────────
    fig, axes = plt.subplots(3, 5, figsize=(22, 12))
    axes_flat = axes.flatten()

    for idx, sig in enumerate(scored[:15]):
        ax = axes_flat[idx]
        for _, row in sample_df.iterrows():
            color = RESPONSE_PAL[row["response"]]
            marker = "o" if row["response"] == "Responder" else "s"
            ax.scatter(row["response"].replace("_", "\n"), row[sig],
                       c=color, marker=marker, s=80, edgecolors="black",
                       linewidth=0.5, zorder=3)

        # Means
        for grp, col in RESPONSE_PAL.items():
            mean_v = sample_df.loc[sample_df["response"] == grp, sig].mean()
            ax.hlines(mean_v, -0.3, 0.3 if grp == "Responder" else 1.3,
                      transform=ax.get_yaxis_transform(), colors=col,
                      linewidth=2, linestyle="--", alpha=0.7)

        row_res = res_df[res_df["signature"] == sig].iloc[0]
        p_perm = row_res["p_permutation_9999"]
        stars = "***" if p_perm < 0.001 else ("**" if p_perm < 0.01 else
                ("*" if p_perm < 0.05 else ("†" if p_perm < 0.1 else "ns")))
        ax.set_title(f"{sig.replace('_', ' ')}\nperm p={p_perm:.3f} {stars}",
                     fontsize=8, fontweight="bold")
        ax.tick_params(labelsize=8)

    for idx in range(len(scored[:15]), 15):
        axes_flat[idx].axis("off")

    plt.suptitle("GSE238264: Sample-Level IO Signatures (n=4 R vs n=3 NR)\n"
                 "Permutation test (9999 resamples) — dots = individual samples",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG / "02b_ds_sample_level_permutation.png", dpi=200, bbox_inches="tight")
    plt.close()

    print(f"  Table: {res_df.shape[0]} signatures")
    print("  Significant (p<0.05):", res_df[res_df["p_permutation_9999"] < 0.05]["signature"].tolist())
    return sample_df, res_df


# ═══════════════════════════════════════════════════════════════════════
#  2. LEAVE-ONE-OUT SENSITIVITY ANALYSIS
# ═══════════════════════════════════════════════════════════════════════

def leave_one_out_sensitivity(datasets):
    """Check if results hold when excluding each sample."""
    print("\n" + "=" * 70)
    print("[2] Leave-One-Out Sensitivity Analysis")
    print("=" * 70)

    scored = [s for s in IO_SIGS.keys() if s in list(datasets.values())[0].obs.columns]
    all_names = list(datasets.keys())

    # Compute sample means
    sample_means = {}
    for sname, adata in datasets.items():
        sample_means[sname] = {sig: adata.obs[sig].mean() for sig in scored}
        sample_means[sname]["response"] = adata.obs["response"].iloc[0]

    loo_rows = []
    for excluded in all_names:
        remaining = {k: v for k, v in sample_means.items() if k != excluded}
        r_names = [k for k, v in remaining.items() if v["response"] == "Responder"]
        nr_names = [k for k, v in remaining.items() if v["response"] == "Non_Responder"]

        if len(r_names) < 2 or len(nr_names) < 2:
            continue

        for sig in scored:
            r_vals = [remaining[k][sig] for k in r_names]
            nr_vals = [remaining[k][sig] for k in nr_names]
            eff = np.mean(r_vals) - np.mean(nr_vals)
            try:
                u, p = mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
            except Exception:
                p = np.nan
            loo_rows.append({
                "excluded": excluded,
                "excluded_response": sample_means[excluded]["response"],
                "signature": sig,
                "effect_size": eff,
                "p_value": p,
                "n_R": len(r_names),
                "n_NR": len(nr_names),
            })

    loo_df = pd.DataFrame(loo_rows)
    loo_df.to_csv(TAB / "gse238264_leave_one_out.csv", index=False)

    # Figure: effect sizes with LOO
    fig, ax = plt.subplots(figsize=(14, 8))
    pivot = loo_df.pivot_table(index="signature", columns="excluded", values="effect_size")

    # Add "Full" column
    full_eff = {}
    r_names_full = [k for k, v in sample_means.items() if v["response"] == "Responder"]
    nr_names_full = [k for k, v in sample_means.items() if v["response"] == "Non_Responder"]
    for sig in scored:
        full_eff[sig] = np.mean([sample_means[k][sig] for k in r_names_full]) - \
                        np.mean([sample_means[k][sig] for k in nr_names_full])
    pivot["Full (n=7)"] = pd.Series(full_eff)

    # Sort by full effect
    pivot = pivot.loc[pd.Series(full_eff).sort_values(ascending=True).index]

    pivot_plot = pivot.copy()
    pivot_plot.index = [s.replace("_", " ") for s in pivot_plot.index]

    sns.heatmap(pivot_plot, cmap="RdBu_r", center=0, annot=True, fmt=".3f",
                linewidths=0.5, ax=ax,
                cbar_kws={"label": "Effect Size (R − NR)", "shrink": 0.7})
    ax.set_title("Leave-One-Out Sensitivity: Effect Sizes (sample-level means)\n"
                 "Each column excludes one sample — 'Full' uses all 7",
                 fontweight="bold", fontsize=11)
    ax.tick_params(labelsize=8)

    # Highlight cells where sign flips vs Full
    for i, sig in enumerate(pivot.index):
        full_sign = np.sign(pivot.loc[sig, "Full (n=7)"])
        for j, col in enumerate(pivot.columns):
            if col == "Full (n=7)":
                continue
            val_sign = np.sign(pivot.loc[sig, col])
            if val_sign != full_sign and full_sign != 0:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False,
                             edgecolor="black", linewidth=2.5))

    plt.tight_layout()
    plt.savefig(FIG / "02b_ds_leave_one_out.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Count sign flips
    n_flips = 0
    for sig in pivot.index:
        full_sign = np.sign(pivot.loc[sig, "Full (n=7)"])
        for col in [c for c in pivot.columns if c != "Full (n=7)"]:
            if np.sign(pivot.loc[sig, col]) != full_sign and full_sign != 0:
                n_flips += 1
    print(f"  Sign flips: {n_flips} (boxed in figure)")
    return loo_df


# ═══════════════════════════════════════════════════════════════════════
#  3. SPATIALLY VARIABLE GENES (SVGs) per section
# ═══════════════════════════════════════════════════════════════════════

def spatially_variable_genes(datasets):
    """Identify top SVGs per section using Moran's I on individual genes."""
    print("\n" + "=" * 70)
    print("[3] Spatially Variable Genes (Moran's I per section)")
    print("=" * 70)

    top_n = 20
    all_svg_rows = []

    for sname, adata in datasets.items():
        print(f"  {sname}: computing SVGs...")
        # Use HVGs for speed
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        hvg_names = adata.var_names[adata.var["highly_variable"]].tolist()

        try:
            sq.gr.spatial_autocorr(adata, mode="moran", genes=hvg_names,
                                   n_perms=99, n_jobs=1)
            mi = adata.uns["moranI"].copy()
            mi["sample"] = sname
            mi["response"] = adata.obs["response"].iloc[0]
            mi["gene"] = mi.index
            mi = mi.sort_values("I", ascending=False)
            all_svg_rows.append(mi.head(top_n))
        except Exception as e:
            print(f"    [WARN] SVG for {sname} failed: {e}")

    if not all_svg_rows:
        print("  [WARN] No SVGs computed")
        return None

    svg_df = pd.concat(all_svg_rows, ignore_index=True)
    svg_df.to_csv(TAB / "gse238264_top_SVGs_per_section.csv", index=False)
    print(f"  Total SVGs: {svg_df.shape[0]} rows")

    # ── Figure: top SVGs heatmap ───────────────────────────────────────
    # Pivot: gene × sample for top SVGs (union of top-10 per sample)
    top10_genes = set()
    for sname in datasets:
        sub = svg_df[svg_df["sample"] == sname].head(10)
        top10_genes.update(sub["gene"].values)

    pivot = svg_df[svg_df["gene"].isin(top10_genes)].pivot_table(
        index="gene", columns="sample", values="I", aggfunc="first"
    ).fillna(0)

    # Sort by mean Moran's I
    pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=False).index]

    if pivot.shape[0] > 40:
        pivot = pivot.head(40)

    fig, ax = plt.subplots(figsize=(10, max(6, 0.3 * pivot.shape[0])))
    sns.heatmap(pivot, cmap="YlOrRd", annot=True, fmt=".2f",
                linewidths=0.3, ax=ax,
                cbar_kws={"label": "Moran's I", "shrink": 0.7})
    ax.set_title("Top Spatially Variable Genes per Visium Section (Moran's I)\n"
                 "Union of top-10 SVGs per section", fontweight="bold", fontsize=10)
    ax.tick_params(labelsize=7)

    # Color sample labels by response
    for tick_l in ax.get_xticklabels():
        sn = tick_l.get_text()
        if sn in datasets:
            resp = datasets[sn].obs["response"].iloc[0]
            tick_l.set_color(RESPONSE_PAL.get(resp, "black"))
            tick_l.set_fontweight("bold")

    plt.tight_layout()
    plt.savefig(FIG / "02b_ds_spatially_variable_genes.png", dpi=200, bbox_inches="tight")
    plt.close()

    # ── R-specific vs NR-specific SVGs ─────────────────────────────────
    r_svgs = svg_df[svg_df["response"] == "Responder"]
    nr_svgs = svg_df[svg_df["response"] == "Non_Responder"]

    r_top = r_svgs.groupby("gene")["I"].mean().nlargest(30)
    nr_top = nr_svgs.groupby("gene")["I"].mean().nlargest(30)

    r_specific = set(r_top.index) - set(nr_top.index)
    nr_specific = set(nr_top.index) - set(r_top.index)
    shared = set(r_top.index) & set(nr_top.index)

    summary = {
        "R-specific top SVGs": sorted(r_specific),
        "NR-specific top SVGs": sorted(nr_specific),
        "Shared top SVGs": sorted(shared),
    }
    for k, v in summary.items():
        print(f"  {k}: {len(v)} genes")

    return svg_df


# ═══════════════════════════════════════════════════════════════════════
#  4. SPATIAL CO-OCCURRENCE OF CELL-TYPE MODULES
# ═══════════════════════════════════════════════════════════════════════

def spatial_cooccurrence(datasets):
    """Spatial co-occurrence of cell-type marker modules using squidpy."""
    print("\n" + "=" * 70)
    print("[4] Spatial Co-Occurrence of Cell-Type Modules")
    print("=" * 70)

    ct_cols = [f"ct_{k}" for k in CELL_MARKERS.keys()]

    # Find common cell-type columns across all datasets
    common_cts = None
    for sname, adata in datasets.items():
        sample_cts = set(c for c in ct_cols if c in adata.obs.columns)
        common_cts = sample_cts if common_cts is None else common_cts & sample_cts
    common_cts = sorted(common_cts)

    # Per-sample co-occurrence heatmap for R vs NR representative
    r_sample = [s for s, a in datasets.items() if a.obs["response"].iloc[0] == "Responder"][0]
    nr_sample = [s for s, a in datasets.items() if a.obs["response"].iloc[0] == "Non_Responder"][0]

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    for ax_idx, (sname, label) in enumerate([(r_sample, "Responder"), (nr_sample, "Non-Responder")]):
        adata = datasets[sname]
        ax = axes[ax_idx]

        # Compute pairwise Spearman correlations of cell-type scores across spatial neighbours
        ct_data = adata.obs[common_cts]
        corr = ct_data.corr(method="spearman")

        labels = [c.replace("ct_", "") for c in corr.columns]
        sns.heatmap(corr.values, xticklabels=labels, yticklabels=labels,
                    cmap="RdBu_r", center=0, annot=True, fmt=".2f",
                    linewidths=0.5, ax=ax, square=True,
                    cbar_kws={"label": "Spearman ρ", "shrink": 0.7})
        ax.set_title(f"{sname} ({label})\nCell-Type Score Correlations",
                     fontweight="bold", fontsize=10)
        ax.tick_params(labelsize=7)

    plt.suptitle("Spatial Co-Occurrence: Cell-Type Module Correlations (R vs NR)",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG / "02b_ds_celltype_cooccurrence.png", dpi=200, bbox_inches="tight")
    plt.close()

    # ── Aggregated R vs NR co-occurrence difference ────────────────────
    r_corrs, nr_corrs = [], []
    for sname, adata in datasets.items():
        corr = adata.obs[common_cts].corr(method="spearman")
        if adata.obs["response"].iloc[0] == "Responder":
            r_corrs.append(corr.values)
        else:
            nr_corrs.append(corr.values)

    if r_corrs and nr_corrs:
        mean_r = np.mean(r_corrs, axis=0)
        mean_nr = np.mean(nr_corrs, axis=0)
        diff = mean_r - mean_nr

        labels = [c.replace("ct_", "") for c in common_cts]

        fig, axes = plt.subplots(1, 3, figsize=(21, 6))
        for ax, data, title in zip(axes, [mean_r, mean_nr, diff],
                                    ["Mean R Correlation", "Mean NR Correlation",
                                     "Difference (R − NR)"]):
            cmap = "RdBu_r"
            sns.heatmap(data, xticklabels=labels, yticklabels=labels,
                        cmap=cmap, center=0, annot=True, fmt=".2f",
                        linewidths=0.5, ax=ax, square=True,
                        cbar_kws={"shrink": 0.7})
            ax.set_title(title, fontweight="bold", fontsize=10)
            ax.tick_params(labelsize=7)

        plt.suptitle("Cell-Type Co-Occurrence: Responder vs Non-Responder (mean across sections)",
                     fontsize=12, fontweight="bold")
        plt.tight_layout()
        plt.savefig(FIG / "02b_ds_cooccurrence_RvsNR_diff.png", dpi=200, bbox_inches="tight")
        plt.close()

    print("  Done.")


# ═══════════════════════════════════════════════════════════════════════
#  5. SPATIAL LEIDEN CLUSTERING PER SECTION
# ═══════════════════════════════════════════════════════════════════════

def spatial_clustering(datasets):
    """Leiden clustering on spatial graph, compare cluster composition R vs NR."""
    print("\n" + "=" * 70)
    print("[5] Spatial Leiden Clustering")
    print("=" * 70)

    cluster_compositions = []
    for sname, adata in datasets.items():
        resp = adata.obs["response"].iloc[0]
        # PCA on expression
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        sc.tl.pca(adata, n_comps=30)
        # Cluster on spatial graph
        sc.tl.leiden(adata, adjacency=adata.obsp["spatial_connectivities"],
                     resolution=0.5, key_added="spatial_cluster")
        n_clusters = adata.obs["spatial_cluster"].nunique()
        print(f"  {sname}: {n_clusters} spatial clusters")

        # Characterise clusters by IO signature means
        scored = [s for s in IO_SIGS.keys() if s in adata.obs.columns]
        cluster_means = adata.obs.groupby("spatial_cluster")[scored].mean()
        cluster_means["sample"] = sname
        cluster_means["response"] = resp
        cluster_means["cluster"] = cluster_means.index
        cluster_compositions.append(cluster_means)

    comp_df = pd.concat(cluster_compositions, ignore_index=True)
    comp_df.to_csv(TAB / "gse238264_spatial_cluster_profiles.csv", index=False)

    # Figure: spatial cluster maps for representative R and NR
    r_sample = [s for s, a in datasets.items() if a.obs["response"].iloc[0] == "Responder"][0]
    nr_sample = [s for s, a in datasets.items() if a.obs["response"].iloc[0] == "Non_Responder"][-1]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, sname, label in zip(axes, [r_sample, nr_sample], ["Responder", "Non-Responder"]):
        adata = datasets[sname]
        coords = adata.obsm["spatial"]
        clusters = adata.obs["spatial_cluster"].astype(int).values
        n_clust = len(np.unique(clusters))
        cmap = plt.cm.get_cmap("tab10", n_clust)

        for c in range(n_clust):
            mask = clusters == c
            ax.scatter(coords[mask, 0], coords[mask, 1],
                       c=[cmap(c)], s=8, label=f"C{c}", edgecolors="none",
                       rasterized=True)
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{sname} ({label})\n{n_clust} spatial clusters", fontweight="bold")
        ax.legend(fontsize=7, loc="lower right", framealpha=0.8)

    plt.suptitle("Spatial Leiden Clusters on Visium Sections",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG / "02b_ds_spatial_clusters.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Figure: cluster signature profiles heatmap
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    for ax, sname, label in zip(axes, [r_sample, nr_sample], ["Responder", "Non-Responder"]):
        sub = comp_df[comp_df["sample"] == sname]
        scored = [s for s in IO_SIGS.keys() if s in sub.columns]
        heat_data = sub.set_index("cluster")[scored].T
        heat_data.index = [s.replace("_", " ") for s in heat_data.index]
        sns.heatmap(heat_data, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
                    linewidths=0.3, ax=ax, cbar_kws={"label": "Mean Score", "shrink": 0.7})
        ax.set_title(f"{sname} ({label})\nCluster IO Signature Profiles", fontweight="bold", fontsize=10)
        ax.set_xlabel("Spatial Cluster")
        ax.tick_params(labelsize=7)

    plt.suptitle("IO Signature Profiles per Spatial Cluster",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG / "02b_ds_cluster_profiles.png", dpi=200, bbox_inches="tight")
    plt.close()

    return comp_df


# ═══════════════════════════════════════════════════════════════════════
#  6. IMMUNE INFILTRATION GRADIENT (distance to tumour core)
# ═══════════════════════════════════════════════════════════════════════

def immune_gradient_analysis(datasets):
    """Analyse immune signature gradient from tumour core to periphery."""
    print("\n" + "=" * 70)
    print("[6] Immune Infiltration Gradient Analysis")
    print("=" * 70)

    gradient_rows = []
    n_bins = 5

    for sname, adata in datasets.items():
        resp = adata.obs["response"].iloc[0]
        coords = adata.obsm["spatial"]

        # Identify tumour core: spots with highest tumour marker score
        if "ct_Tumor" in adata.obs:
            tumor_scores = adata.obs["ct_Tumor"].values
        elif "ct_Hepatocyte" in adata.obs:
            tumor_scores = adata.obs["ct_Hepatocyte"].values
        else:
            continue

        # Top 20% tumour-score spots = core
        threshold = np.percentile(tumor_scores, 80)
        core_mask = tumor_scores >= threshold
        core_coords = coords[core_mask]

        if core_coords.shape[0] < 5:
            continue

        # Distance to nearest core spot
        centroid = np.median(core_coords, axis=0)
        dists = np.sqrt(np.sum((coords - centroid) ** 2, axis=1))

        # Normalise to [0, 1]
        dists_norm = (dists - dists.min()) / (dists.max() - dists.min() + 1e-12)
        adata.obs["dist_to_core"] = dists_norm

        # Bin into distance zones
        bins = np.linspace(0, 1, n_bins + 1)
        adata.obs["dist_zone"] = pd.cut(dists_norm, bins=bins, labels=[f"Z{i+1}" for i in range(n_bins)])

        scored = [s for s in IO_SIGS.keys() if s in adata.obs.columns]
        for zone_label in [f"Z{i+1}" for i in range(n_bins)]:
            mask = adata.obs["dist_zone"] == zone_label
            if mask.sum() < 5:
                continue
            row = {"sample": sname, "response": resp, "zone": zone_label,
                   "n_spots": int(mask.sum()),
                   "mean_dist": dists_norm[mask].mean()}
            for sig in scored:
                row[f"mean_{sig}"] = adata.obs.loc[mask, sig].mean()
            gradient_rows.append(row)

    grad_df = pd.DataFrame(gradient_rows)
    grad_df.to_csv(TAB / "gse238264_immune_gradient.csv", index=False)

    if grad_df.empty:
        print("  [WARN] No gradient data")
        return

    # Figure: gradient line plots for key signatures
    key_sigs = ["T_cell_cytotoxicity", "TLS_signature", "TGF_beta",
                "Antigen_presentation", "IFN_gamma", "M2_macrophage"]
    key_sigs = [s for s in key_sigs if f"mean_{s}" in grad_df.columns]

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes_flat = axes.flatten()

    for idx, sig in enumerate(key_sigs[:6]):
        ax = axes_flat[idx]
        for resp, color in RESPONSE_PAL.items():
            sub = grad_df[grad_df["response"] == resp]
            means = sub.groupby("zone")[f"mean_{sig}"].agg(["mean", "std"]).reset_index()
            x = range(len(means))
            ax.errorbar(x, means["mean"], yerr=means["std"], color=color,
                        marker="o", capsize=4, linewidth=2, markersize=6,
                        label=resp.replace("_", " "), alpha=0.85)
        ax.set_xticks(range(n_bins))
        ax.set_xticklabels([f"Z{i+1}\n(core)" if i == 0 else
                           (f"Z{i+1}\n(edge)" if i == n_bins - 1 else f"Z{i+1}")
                           for i in range(n_bins)], fontsize=8)
        ax.set_ylabel("Mean Score", fontsize=9)
        ax.set_title(sig.replace("_", " "), fontweight="bold", fontsize=10)
        ax.legend(fontsize=7)

        # Spearman correlation with distance for R and NR
        for resp, color in RESPONSE_PAL.items():
            sub = grad_df[grad_df["response"] == resp]
            if len(sub) > 3:
                rho, p = spearmanr(sub["mean_dist"], sub[f"mean_{sig}"])
                ax.annotate(f"{resp[:1]}: ρ={rho:.2f}", xy=(0.02, 0.95 if resp == "Responder" else 0.87),
                            xycoords="axes fraction", fontsize=7, color=color)

    for idx in range(len(key_sigs[:6]), 6):
        axes_flat[idx].axis("off")

    plt.suptitle("Immune Infiltration Gradient: Core → Periphery\n"
                 "Mean ± SD across samples per distance zone (R vs NR)",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG / "02b_ds_immune_gradient.png", dpi=200, bbox_inches="tight")
    plt.close()

    # ── Spatial gradient map for one R and one NR ──────────────────────
    r_name = [s for s, a in datasets.items() if a.obs["response"].iloc[0] == "Responder"][0]
    nr_name = [s for s, a in datasets.items() if a.obs["response"].iloc[0] == "Non_Responder"][-1]

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    for row_idx, (sname, label) in enumerate([(r_name, "Responder"), (nr_name, "Non-Responder")]):
        adata = datasets[sname]
        coords = adata.obsm["spatial"]
        for col_idx, (var, cmap_name) in enumerate([
            ("dist_to_core", "viridis"),
            ("T_cell_cytotoxicity", "RdYlBu_r"),
            ("TLS_signature", "RdYlBu_r"),
        ]):
            ax = axes[row_idx, col_idx]
            if var in adata.obs:
                sp = ax.scatter(coords[:, 0], coords[:, 1],
                                c=adata.obs[var].values, cmap=cmap_name,
                                s=8, edgecolors="none", rasterized=True)
                plt.colorbar(sp, ax=ax, shrink=0.6, pad=0.02)
            ax.set_aspect("equal"); ax.invert_yaxis()
            ax.set_xticks([]); ax.set_yticks([])
            if row_idx == 0:
                ax.set_title(var.replace("_", " "), fontsize=10, fontweight="bold")
            if col_idx == 0:
                ax.set_ylabel(f"{sname} ({label})", fontsize=9, fontweight="bold")

    plt.suptitle("Spatial Gradient: Distance to Tumour Core & Immune Signatures",
                 fontsize=12, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(FIG / "02b_ds_gradient_spatial_maps.png", dpi=200, bbox_inches="tight")
    plt.close()

    print(f"  Gradient table: {grad_df.shape}")


# ═══════════════════════════════════════════════════════════════════════
#  7. SPATIAL ENTROPY / HETEROGENEITY
# ═══════════════════════════════════════════════════════════════════════

def spatial_entropy(datasets):
    """Shannon entropy of IO signatures per section as heterogeneity metric."""
    print("\n" + "=" * 70)
    print("[7] Spatial Entropy / Heterogeneity")
    print("=" * 70)

    scored = [s for s in IO_SIGS.keys() if s in list(datasets.values())[0].obs.columns]

    entropy_rows = []
    for sname, adata in datasets.items():
        resp = adata.obs["response"].iloc[0]
        row = {"sample": sname, "response": resp}

        for sig in scored:
            vals = adata.obs[sig].values
            # Discretise into 10 bins for entropy
            counts, _ = np.histogram(vals, bins=10)
            probs = counts / counts.sum()
            probs = probs[probs > 0]
            entropy = -np.sum(probs * np.log2(probs))
            row[f"entropy_{sig}"] = entropy

            # Also compute coefficient of variation
            row[f"cv_{sig}"] = np.std(vals) / (np.abs(np.mean(vals)) + 1e-12)

        entropy_rows.append(row)

    ent_df = pd.DataFrame(entropy_rows)
    ent_df.to_csv(TAB / "gse238264_spatial_entropy.csv", index=False)

    # Figure: entropy comparison R vs NR
    ent_cols = [c for c in ent_df.columns if c.startswith("entropy_")]

    fig, axes = plt.subplots(3, 5, figsize=(22, 12))
    axes_flat = axes.flatten()

    for idx, col in enumerate(ent_cols[:15]):
        ax = axes_flat[idx]
        sig_name = col.replace("entropy_", "")
        for _, row in ent_df.iterrows():
            color = RESPONSE_PAL[row["response"]]
            marker = "o" if row["response"] == "Responder" else "s"
            ax.scatter(row["response"].replace("_", "\n"), row[col],
                       c=color, marker=marker, s=80, edgecolors="black",
                       linewidth=0.5, zorder=3)

        # Means
        for grp, col_c in RESPONSE_PAL.items():
            mean_v = ent_df.loc[ent_df["response"] == grp, col].mean()
            ax.axhline(mean_v, color=col_c, linewidth=1.5, linestyle="--", alpha=0.5)

        ax.set_title(sig_name.replace("_", " "), fontsize=8, fontweight="bold")
        ax.set_ylabel("Shannon Entropy", fontsize=7)

    for idx in range(len(ent_cols[:15]), 15):
        axes_flat[idx].axis("off")

    plt.suptitle("Spatial Heterogeneity: Shannon Entropy of IO Signature Distributions\n"
                 "Higher entropy = more spatially heterogeneous signal",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG / "02b_ds_spatial_entropy.png", dpi=200, bbox_inches="tight")
    plt.close()

    print(f"  Entropy table: {ent_df.shape}")
    return ent_df


# ═══════════════════════════════════════════════════════════════════════
#  8. LIGAND–RECEPTOR CO-LOCALIZATION
# ═══════════════════════════════════════════════════════════════════════

def ligand_receptor_analysis(datasets):
    """Spatial co-expression of L-R pairs with neighbourhood-level correlation."""
    print("\n" + "=" * 70)
    print("[8] Ligand–Receptor Spatial Co-Localization")
    print("=" * 70)

    lr_rows = []
    for sname, adata in datasets.items():
        resp = adata.obs["response"].iloc[0]
        avail = set(adata.var_names)

        for lig, rec in LR_PAIRS:
            if lig not in avail or rec not in avail:
                continue

            # Get expression
            if issparse(adata.X):
                lig_expr = np.asarray(adata[:, lig].X.todense()).flatten()
                rec_expr = np.asarray(adata[:, rec].X.todense()).flatten()
            else:
                lig_expr = adata[:, lig].X.flatten()
                rec_expr = adata[:, rec].X.flatten()

            # Global Spearman
            rho_global, p_global = spearmanr(lig_expr, rec_expr)

            # Neighbour-averaged co-expression (ligand in spot, receptor in neighbours)
            adj = adata.obsp["spatial_connectivities"]
            # Average receptor expression in neighbours
            rec_neigh = np.asarray(adj.dot(rec_expr.reshape(-1, 1)) /
                                   (adj.sum(axis=1) + 1e-12)).flatten()
            rho_spatial, p_spatial = spearmanr(lig_expr, rec_neigh)

            # Fraction of spots co-expressing (both > 0)
            coexpr_frac = np.mean((lig_expr > 0) & (rec_expr > 0))

            lr_rows.append({
                "sample": sname, "response": resp,
                "ligand": lig, "receptor": rec,
                "rho_global": rho_global, "p_global": p_global,
                "rho_spatial": rho_spatial, "p_spatial": p_spatial,
                "coexpr_fraction": coexpr_frac,
                "mean_ligand": lig_expr.mean(), "mean_receptor": rec_expr.mean(),
            })

    lr_df = pd.DataFrame(lr_rows)
    lr_df.to_csv(TAB / "gse238264_ligand_receptor.csv", index=False)
    print(f"  L-R table: {lr_df.shape[0]} rows")

    if lr_df.empty:
        return

    # Figure: spatial L-R correlation heatmap (R vs NR)
    lr_df["pair"] = lr_df["ligand"] + " → " + lr_df["receptor"]

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    for ax, resp, label in zip(axes, ["Responder", "Non_Responder"],
                                ["Responder", "Non-Responder"]):
        sub = lr_df[lr_df["response"] == resp]
        pivot = sub.pivot_table(index="pair", columns="sample", values="rho_spatial")
        pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=False).index]

        sns.heatmap(pivot, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
                    linewidths=0.3, ax=ax, cbar_kws={"label": "Spatial ρ", "shrink": 0.7})
        ax.set_title(f"{label}\nL-R Spatial Co-Localization (neighbour ρ)",
                     fontweight="bold", fontsize=10)
        ax.tick_params(labelsize=7)

    plt.suptitle("Ligand–Receptor Spatial Co-Localization: R vs NR",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG / "02b_ds_ligand_receptor_heatmap.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Figure: R vs NR mean spatial rho comparison
    mean_lr = lr_df.groupby(["pair", "response"])["rho_spatial"].mean().unstack(fill_value=0)
    if "Responder" in mean_lr.columns and "Non_Responder" in mean_lr.columns:
        mean_lr["diff"] = mean_lr["Responder"] - mean_lr["Non_Responder"]
        mean_lr = mean_lr.sort_values("diff")

        fig, ax = plt.subplots(figsize=(10, 8))
        y = np.arange(len(mean_lr))
        w = 0.35
        ax.barh(y - w/2, mean_lr["Responder"], w, color=RESPONSE_PAL["Responder"],
                label="R", alpha=0.85, edgecolor="black", linewidth=0.3)
        ax.barh(y + w/2, mean_lr["Non_Responder"], w, color=RESPONSE_PAL["Non_Responder"],
                label="NR", alpha=0.85, edgecolor="black", linewidth=0.3)
        ax.set_yticks(y)
        ax.set_yticklabels(mean_lr.index, fontsize=8)
        ax.set_xlabel("Mean Spatial ρ (ligand–neighbour receptor)")
        ax.axvline(0, color="black", linewidth=0.5)
        ax.legend(fontsize=9)
        ax.set_title("L-R Spatial Co-Localization: Responder vs Non-Responder\n"
                     "Mean across sections", fontweight="bold")
        plt.tight_layout()
        plt.savefig(FIG / "02b_ds_LR_RvsNR_comparison.png", dpi=200, bbox_inches="tight")
        plt.close()


# ═══════════════════════════════════════════════════════════════════════
#  9. COMPREHENSIVE DOWNSTREAM SUMMARY FIGURE
# ═══════════════════════════════════════════════════════════════════════

def downstream_summary(perm_df, loo_df, ent_df, datasets):
    """Multi-panel summary of all downstream analyses."""
    print("\n" + "=" * 70)
    print("[9] Downstream Summary Figure")
    print("=" * 70)

    fig = plt.figure(figsize=(20, 14))
    gs = gridspec.GridSpec(2, 3, hspace=0.35, wspace=0.35)

    # Panel A: sample-level permutation volcano
    ax = fig.add_subplot(gs[0, 0])
    if perm_df is not None:
        perm_df["neg_log10_p"] = -np.log10(perm_df["p_permutation_9999"].clip(lower=1e-6))
        colors = ["#2ecc71" if e > 0 else "#e74c3c" for e in perm_df["effect_size"]]
        ax.scatter(perm_df["effect_size"], perm_df["neg_log10_p"],
                   c=colors, s=60, edgecolors="black", linewidth=0.5)
        ax.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.8)
        ax.axvline(0, color="gray", linewidth=0.5)
        for _, row in perm_df.iterrows():
            if row["p_permutation_9999"] < 0.1:
                ax.annotate(row["signature"].replace("_", " "),
                           (row["effect_size"], row["neg_log10_p"]),
                           fontsize=6, ha="center", va="bottom")
        ax.set_xlabel("Effect Size (R − NR)")
        ax.set_ylabel("−log₁₀(permutation p)")
        ax.set_title("A. Sample-Level Permutation Tests\n(n=4 R vs n=3 NR)", fontweight="bold", fontsize=10)

    # Panel B: entropy R vs NR
    ax = fig.add_subplot(gs[0, 1])
    if ent_df is not None:
        ent_cols = [c for c in ent_df.columns if c.startswith("entropy_")]
        # Mean entropy across all sigs
        ent_df["mean_entropy"] = ent_df[ent_cols].mean(axis=1)
        for _, row in ent_df.iterrows():
            color = RESPONSE_PAL[row["response"]]
            ax.scatter(row["response"].replace("_", " "), row["mean_entropy"],
                       c=color, s=100, edgecolors="black", linewidth=0.5, zorder=3)
        ax.set_ylabel("Mean Shannon Entropy")
        ax.set_title("B. Spatial Heterogeneity\n(mean entropy across IO sigs)", fontweight="bold", fontsize=10)

    # Panel C: LOO summary (number of sign flips per signature)
    ax = fig.add_subplot(gs[0, 2])
    if loo_df is not None and not loo_df.empty:
        # Compute full effect signs
        scored = loo_df["signature"].unique()
        sample_means_dict = {}
        for sname, adata in datasets.items():
            resp = adata.obs["response"].iloc[0]
            for sig in scored:
                if sig in adata.obs:
                    sample_means_dict.setdefault(sig, []).append(
                        (resp, adata.obs[sig].mean()))

        flip_counts = {}
        for sig in scored:
            pairs = sample_means_dict.get(sig, [])
            r_mean = np.mean([v for r, v in pairs if r == "Responder"])
            nr_mean = np.mean([v for r, v in pairs if r == "Non_Responder"])
            full_sign = np.sign(r_mean - nr_mean)

            sub = loo_df[loo_df["signature"] == sig]
            flips = sum(1 for _, row in sub.iterrows()
                        if np.sign(row["effect_size"]) != full_sign and full_sign != 0)
            flip_counts[sig] = flips

        flip_s = pd.Series(flip_counts).sort_values(ascending=False)
        colors = ["#e74c3c" if v > 0 else "#2ecc71" for v in flip_s.values]
        ax.barh(range(len(flip_s)), flip_s.values, color=colors, height=0.5,
                edgecolor="black", linewidth=0.3)
        ax.set_yticks(range(len(flip_s)))
        ax.set_yticklabels([s.replace("_", " ") for s in flip_s.index], fontsize=7)
        ax.set_xlabel("Sign Flips (out of 7 LOO)")
        ax.set_title("C. LOO Sensitivity\n(0 flips = robust)", fontweight="bold", fontsize=10)

    # Panel D: LR co-localization top pairs
    ax = fig.add_subplot(gs[1, 0])
    lr_file = TAB / "gse238264_ligand_receptor.csv"
    if lr_file.exists():
        lr_df = pd.read_csv(lr_file)
        lr_df["pair"] = lr_df["ligand"] + "→" + lr_df["receptor"]
        mean_lr = lr_df.groupby(["pair", "response"])["rho_spatial"].mean().unstack(fill_value=0)
        if "Responder" in mean_lr.columns and "Non_Responder" in mean_lr.columns:
            mean_lr["diff"] = mean_lr["Responder"] - mean_lr["Non_Responder"]
            mean_lr = mean_lr.sort_values("diff")
            top_diff = pd.concat([mean_lr.head(5), mean_lr.tail(5)])
            y = np.arange(len(top_diff))
            ax.barh(y, top_diff["diff"], color=["#e74c3c" if d < 0 else "#2ecc71" for d in top_diff["diff"]],
                    height=0.5, edgecolor="black", linewidth=0.3)
            ax.set_yticks(y)
            ax.set_yticklabels(top_diff.index, fontsize=7)
            ax.axvline(0, color="black", linewidth=0.5)
            ax.set_xlabel("Δ Spatial ρ (R − NR)")
            ax.set_title("D. L-R Co-Localization Differences\n(top R-enriched and NR-enriched)", fontweight="bold", fontsize=10)

    # Panel E: immune gradient key sigs
    ax = fig.add_subplot(gs[1, 1])
    grad_file = TAB / "gse238264_immune_gradient.csv"
    if grad_file.exists():
        grad_df = pd.read_csv(grad_file)
        sig = "T_cell_cytotoxicity"
        col = f"mean_{sig}"
        if col in grad_df.columns:
            for resp, color in RESPONSE_PAL.items():
                sub = grad_df[grad_df["response"] == resp]
                means = sub.groupby("zone")[col].mean()
                ax.plot(range(len(means)), means.values, color=color,
                        marker="o", linewidth=2, markersize=6,
                        label=resp.replace("_", " "))
            ax.set_xticks(range(5))
            ax.set_xticklabels(["Core", "Z2", "Z3", "Z4", "Edge"], fontsize=8)
            ax.set_ylabel("Mean Score")
            ax.legend(fontsize=8)
            ax.set_title("E. T-Cell Cytotoxicity Gradient\n(core → periphery)", fontweight="bold", fontsize=10)

    # Panel F: SVG overlap
    ax = fig.add_subplot(gs[1, 2])
    svg_file = TAB / "gse238264_top_SVGs_per_section.csv"
    if svg_file.exists():
        svg_df = pd.read_csv(svg_file)
        r_svgs = set(svg_df[svg_df["response"] == "Responder"]["gene"])
        nr_svgs = set(svg_df[svg_df["response"] == "Non_Responder"]["gene"])
        shared = r_svgs & nr_svgs
        r_only = r_svgs - nr_svgs
        nr_only = nr_svgs - r_svgs

        try:
            from matplotlib_venn import venn2
            venn2([r_svgs, nr_svgs], set_labels=("R SVGs", "NR SVGs"),
                  set_colors=(RESPONSE_PAL["Responder"], RESPONSE_PAL["Non_Responder"]),
                  alpha=0.6, ax=ax)
            ax.set_title("F. Spatially Variable Gene Overlap\n(top 20 per section)", fontweight="bold", fontsize=10)
        except Exception:
            ax.bar(["R-only", "Shared", "NR-only"],
                   [len(r_only), len(shared), len(nr_only)],
                   color=[RESPONSE_PAL["Responder"], "#95a5a6", RESPONSE_PAL["Non_Responder"]],
                   width=0.4)
            ax.set_title("F. SVG Overlap (R vs NR)", fontweight="bold", fontsize=10)

    plt.suptitle("Downstream Spatial Analyses — GSE238264 HCC IO Visium\n"
                 "7 sections (4 R, 3 NR) | 44,516 total spots",
                 fontsize=13, fontweight="bold")
    plt.savefig(FIG / "02b_ds_comprehensive_summary.png", dpi=200, bbox_inches="tight")
    plt.close()
    print("  Summary figure saved.")


# ═══════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("02b — Downstream Visium Spatial Analyses")
    print("=" * 70)

    # Load GSE238264
    print("\n[0] Loading GSE238264...")
    datasets = load_gse238264()

    # 1. Sample-level permutation tests
    sample_df, perm_df = sample_level_permutation_tests(datasets)

    # 2. Leave-one-out sensitivity
    loo_df = leave_one_out_sensitivity(datasets)

    # 3. Spatially variable genes
    svg_df = spatially_variable_genes(datasets)

    # 4. Spatial co-occurrence
    spatial_cooccurrence(datasets)

    # 5. Spatial clustering
    cluster_df = spatial_clustering(datasets)

    # 6. Immune gradient
    immune_gradient_analysis(datasets)

    # 7. Spatial entropy
    ent_df = spatial_entropy(datasets)

    # 8. Ligand–receptor analysis
    ligand_receptor_analysis(datasets)

    # 9. Summary figure
    downstream_summary(perm_df, loo_df, ent_df, datasets)

    print("\n" + "=" * 70)
    print("[DONE] All downstream analyses complete")
    print(f"  Figures: {FIG}")
    print(f"  Tables:  {TAB}")
    print("=" * 70)


if __name__ == "__main__":
    main()
