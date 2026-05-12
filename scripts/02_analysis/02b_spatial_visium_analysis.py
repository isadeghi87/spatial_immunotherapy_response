#!/usr/bin/env python3
"""
02b_spatial_visium_analysis.py — True Visium Spatial Analysis for HCC IO Response

Datasets:
  GSE238264  7 HCC Visium sections (4 Responder, 3 Non-Responder)
  GSE281759  17 HCC multi-region Visium sections (Normal/Tumor/Primary/Recurrence)

Analyses:
  Part 1 — Per-section QC, normalisation, IO-signature scoring
  Part 2 — GSE238264 R vs NR spatial biomarker comparison
  Part 3 — Spatial autocorrelation (Moran's I) for IO signatures
  Part 4 — Spatial immunotype classification & neighbourhood enrichment
  Part 5 — GSE281759 Tumor vs Normal, Primary vs Recurrence
  Part 6 — Cross-dataset summary

Outputs → out/02b_spatial_biomarkers/figures/ and tables/
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import sys; sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
from config.publication_style import apply_style, PAL, CMAP, label_clean
apply_style()
from pathlib import Path
from scipy.stats import mannwhitneyu, spearmanr, kruskal
from scipy.sparse import issparse
from collections import OrderedDict
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

# ── IO Signatures ──────────────────────────────────────────────────────
IO_SIGS = OrderedDict({
    "T_cell_exhaustion":    ["PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","TOX","ENTPD1"],
    "T_cell_activation":    ["IFNG","GZMB","PRF1","TNF","CD69","IL2"],
    "T_cell_cytotoxicity":  ["GZMA","GZMB","GZMH","GZMK","PRF1","GNLY","NKG7"],
    "Progenitor_CD8":       ["TCF7","SLAMF6","LEF1","CCR7","IL7R","SELL"],
    "Terminal_exhaustion":  ["TOX","ENTPD1","HAVCR2","LAYN","TIGIT","CXCL13"],
    "Treg_signature":       ["FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"],
    "M1_macrophage":        ["NOS2","IL1B","IL6","TNF","CD80","CD86"],
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

# ── Palette ────────────────────────────────────────────────────────────
RESPONSE_PAL = {"Responder": "#0077BB", "Non_Responder": "#CC3311"}
IMMUNOTYPE_PAL = {
    "Immune_Hot": "#CC3311", "Immune_Excluded": "#f39c12",
    "Immune_Desert": "#3498db", "Immune_Mixed": "#95a5a6",
}
TISSUE_PAL = {"Normal": "#3498db", "Tumor": "#CC3311",
              "Primary": "#CC3311", "Recurrence": "#9b59b6"}

# ═══════════════════════════════════════════════════════════════════════
#  HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════

def load_visium_h5ad(path: Path, sample_name: str) -> sc.AnnData:
    """Load a Visium h5ad, set up spatial coords in obsm, normalise."""
    adata = sc.read_h5ad(str(path))
    adata.var_names_make_unique()
    adata.obs_names = [f"{sample_name}_{b}" for b in adata.obs_names]

    # Build obsm['spatial'] from pixel coords
    adata.obsm["spatial"] = adata.obs[["pxl_col", "pxl_row"]].values.astype(float)

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()

    # Basic QC
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
    adata.obs["sample"] = sample_name

    # Normalise
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    return adata


def score_io_signatures(adata):
    """Score all IO signatures. Returns list of scored names."""
    avail = set(adata.var_names)
    scored = []
    for name, genes in IO_SIGS.items():
        present = [g for g in genes if g in avail]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=name, use_raw=False)
            scored.append(name)
    return scored


def score_cell_markers(adata):
    """Score cell-type marker modules."""
    avail = set(adata.var_names)
    scored = []
    for name, genes in CELL_MARKERS.items():
        present = [g for g in genes if g in avail]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=f"ct_{name}", use_raw=False)
            scored.append(f"ct_{name}")
    return scored


def classify_immunotype(adata):
    """Classify spots into spatial immunotypes based on activation vs TGF-β."""
    if "T_cell_activation" not in adata.obs or "TGF_beta" not in adata.obs:
        return
    act = adata.obs["T_cell_activation"].values
    tgfb = adata.obs["TGF_beta"].values
    p66_act, p33_act = np.percentile(act, 66), np.percentile(act, 33)
    p66_tgfb = np.percentile(tgfb, 66)

    labels = []
    for a, t in zip(act, tgfb):
        if a > p66_act and t < p66_tgfb:
            labels.append("Immune_Hot")
        elif a > p66_act and t >= p66_tgfb:
            labels.append("Immune_Excluded")
        elif a < p33_act:
            labels.append("Immune_Desert")
        else:
            labels.append("Immune_Mixed")
    adata.obs["spatial_immunotype"] = pd.Categorical(labels,
        categories=["Immune_Hot","Immune_Excluded","Immune_Desert","Immune_Mixed"])


def build_spatial_graph(adata, n_neighs=6):
    """Build spatial connectivity graph for squidpy."""
    sq.gr.spatial_neighbors(adata, n_neighs=n_neighs, coord_type="generic")


# ═══════════════════════════════════════════════════════════════════════
#  PART 1 — Load all Visium data
# ═══════════════════════════════════════════════════════════════════════

def load_gse238264():
    """Load GSE238264 (7 HCC IO Visium: 4R + 3NR)."""
    d = DATA / "GSE238264_HCC_IO_Visium"
    samples = {
        "HCC1R": "Responder", "HCC2R": "Responder",
        "HCC3R": "Responder", "HCC4R": "Responder",
        "HCC5NR": "Non_Responder", "HCC6NR": "Non_Responder",
        "HCC7NR": "Non_Responder",
    }
    datasets = []
    for name, resp in samples.items():
        path = d / f"{name}.h5ad"
        if not path.exists():
            print(f"  [SKIP] {name}: file not found")
            continue
        print(f"  Loading {name} ({resp})...")
        adata = load_visium_h5ad(path, name)
        adata.obs["response"] = resp
        adata.obs["dataset"] = "GSE238264"
        datasets.append(adata)
        print(f"    {adata.n_obs} spots × {adata.n_vars} genes")
    return datasets


def load_gse281759():
    """Load GSE281759 (17 HCC multi-region Visium)."""
    d = DATA / "GSE281759_HCC_MultiRegion"
    datasets = []
    for h5 in sorted(d.glob("*.h5ad")):
        gsm_name = h5.stem
        # Infer tissue type from sample name
        short = gsm_name.split("_", 1)[1] if "_" in gsm_name else gsm_name
        if short.endswith("N"):
            tissue = "Normal"
        elif short.endswith("P"):
            tissue = "Primary"
        elif short.endswith("R"):
            tissue = "Recurrence"
        else:
            tissue = "Tumor"

        print(f"  Loading {short} ({tissue})...")
        adata = load_visium_h5ad(h5, short)
        adata.obs["tissue_type"] = tissue
        adata.obs["patient"] = short[:4]  # e.g. A009, B001, C004, C008
        adata.obs["dataset"] = "GSE281759"
        datasets.append(adata)
        print(f"    {adata.n_obs} spots × {adata.n_vars} genes")
    return datasets


# ═══════════════════════════════════════════════════════════════════════
#  PART 2 — GSE238264: R vs NR Spatial Comparison
# ═══════════════════════════════════════════════════════════════════════

def analyse_gse238264(datasets):
    """Full spatial analysis for IO-response Visium cohort."""
    print("\n" + "=" * 70)
    print("PART 2: GSE238264 — Responder vs Non-Responder Spatial Analysis")
    print("=" * 70)

    # Score signatures per section, then concatenate
    all_scored = None
    per_sample_stats = []

    for adata in datasets:
        scored = score_io_signatures(adata)
        score_cell_markers(adata)
        classify_immunotype(adata)
        build_spatial_graph(adata)

        sample = adata.obs["sample"].iloc[0]
        resp = adata.obs["response"].iloc[0]

        # Per-sample summary stats
        row = {"sample": sample, "response": resp, "n_spots": adata.n_obs}
        for sig in scored:
            row[f"mean_{sig}"] = adata.obs[sig].mean()
            row[f"std_{sig}"] = adata.obs[sig].std()
        if "spatial_immunotype" in adata.obs:
            vc = adata.obs["spatial_immunotype"].value_counts(normalize=True)
            for it in IMMUNOTYPE_PAL:
                row[f"frac_{it}"] = vc.get(it, 0.0)
        per_sample_stats.append(row)

    stats_df = pd.DataFrame(per_sample_stats)
    stats_df.to_csv(TAB / "gse238264_per_sample_stats.csv", index=False)
    print(f"\n  Per-sample stats: {stats_df.shape}")

    # Concatenate all sections
    combined = sc.concat(datasets, join="inner")
    combined.obs_names_make_unique()
    build_spatial_graph(combined)
    scored = score_io_signatures(combined)
    classify_immunotype(combined)

    print(f"  Combined: {combined.n_obs} spots × {combined.n_vars} genes")
    print(f"  Scored: {len(scored)} signatures")

    # ── Figure 1: Spatial signature maps per sample ────────────────────
    print("\n[Fig 1] Spatial signature maps...")
    key_sigs = ["T_cell_cytotoxicity", "T_cell_exhaustion", "TLS_signature",
                "M1_macrophage", "TGF_beta", "IFN_gamma"]
    key_sigs = [s for s in key_sigs if s in scored]

    n_samples = len(datasets)
    n_sigs = len(key_sigs)
    fig, axes = plt.subplots(n_samples, n_sigs, figsize=(4 * n_sigs, 3.5 * n_samples))
    if n_samples == 1:
        axes = axes[np.newaxis, :]

    for i, adata in enumerate(datasets):
        sname = adata.obs["sample"].iloc[0]
        resp = adata.obs["response"].iloc[0]
        for j, sig in enumerate(key_sigs):
            ax = axes[i, j]
            if sig in adata.obs:
                sc_vals = adata.obs[sig].values
                coords = adata.obsm["spatial"]
                sc_plot = ax.scatter(coords[:, 0], coords[:, 1],
                                     c=sc_vals, cmap="RdYlBu_r", s=6,
                                     edgecolors="none", rasterized=True)
                plt.colorbar(sc_plot, ax=ax, shrink=0.6, pad=0.02)
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.set_xticks([]); ax.set_yticks([])
            if i == 0:
                ax.set_title(sig.replace("_", " "), fontsize=9, fontweight="bold")
            if j == 0:
                resp_tag = "R" if resp == "Responder" else "NR"
                ax.set_ylabel(f"{sname} ({resp_tag})", fontsize=9, fontweight="bold")

    plt.suptitle("GSE238264: Spatial IO Signature Maps (Visium)", fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(FIG / "02b_visium_spatial_signature_maps.png", dpi=200, bbox_inches="tight")
    plt.close()

    # ── Figure 2: R vs NR violin comparison ────────────────────────────
    print("[Fig 2] R vs NR signature violins...")
    fig, axes = plt.subplots(4, 4, figsize=(20, 16))
    axes_flat = axes.flatten()
    stat_rows = []
    for idx, sig in enumerate(scored[:16]):
        ax = axes_flat[idx]
        r_vals = combined.obs.loc[combined.obs["response"] == "Responder", sig].dropna()
        nr_vals = combined.obs.loc[combined.obs["response"] == "Non_Responder", sig].dropna()

        parts = ax.violinplot([r_vals.values, nr_vals.values], positions=[0, 1], showmedians=True)
        for pc, color in zip(parts["bodies"], [RESPONSE_PAL["Responder"], RESPONSE_PAL["Non_Responder"]]):
            pc.set_facecolor(color)
            pc.set_alpha(0.7)

        u_stat, p_val = mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
        effect = r_vals.mean() - nr_vals.mean()
        stat_rows.append({"signature": sig, "mean_R": r_vals.mean(), "mean_NR": nr_vals.mean(),
                          "effect_size": effect, "U_stat": u_stat, "p_value": p_val})

        stars = "***" if p_val < 0.001 else ("**" if p_val < 0.01 else ("*" if p_val < 0.05 else "ns"))
        ax.set_title(f"{sig.replace('_', ' ')}\n{stars} (p={p_val:.2e})", fontsize=8, fontweight="bold")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["R", "NR"], fontsize=9)

    for idx in range(len(scored[:16]), 16):
        axes_flat[idx].axis("off")

    plt.suptitle("GSE238264: IO Signature Scores — Responder vs Non-Responder (spot-level)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG / "02b_visium_RvsNR_violins.png", dpi=200, bbox_inches="tight")
    plt.close()

    stat_df = pd.DataFrame(stat_rows).sort_values("p_value")
    stat_df.to_csv(TAB / "gse238264_RvsNR_signature_stats.csv", index=False)
    print(f"  Stats table: {stat_df.shape[0]} signatures tested")

    # ── Figure 3: Spatial immunotype composition R vs NR ───────────────
    print("[Fig 3] Immunotype composition...")
    if "spatial_immunotype" in combined.obs:
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        # 3a: stacked bar per sample
        ax = axes[0]
        comp = combined.obs.groupby(["sample", "spatial_immunotype"]).size().unstack(fill_value=0)
        comp_norm = comp.div(comp.sum(axis=1), axis=0)
        comp_norm.plot(kind="barh", stacked=True, ax=ax,
                       color=[IMMUNOTYPE_PAL[c] for c in comp_norm.columns if c in IMMUNOTYPE_PAL],
                       width=0.5)
        ax.set_xlabel("Fraction")
        ax.set_title("A. Immunotype Composition per Section", fontweight="bold")
        ax.legend(fontsize=7, loc="lower right")

        # 3b: grouped comparison R vs NR
        ax = axes[1]
        fracs = combined.obs.groupby(["response", "spatial_immunotype"]).size().unstack(fill_value=0)
        fracs = fracs.div(fracs.sum(axis=1), axis=0)
        fracs.plot(kind="bar", ax=ax,
                   color=[IMMUNOTYPE_PAL.get(c, "#666") for c in fracs.columns],
                   width=0.5)
        ax.set_title("B. Immunotype Distribution: R vs NR", fontweight="bold")
        ax.set_ylabel("Fraction")
        ax.legend(fontsize=7)
        ax.tick_params(axis="x", rotation=0)

        # 3c: spatial immunotype maps for two representative samples
        ax = axes[2]
        imm_counts = combined.obs.groupby("response")["spatial_immunotype"].value_counts(normalize=True).unstack()
        if imm_counts is not None:
            r_hot = imm_counts.loc["Responder", "Immune_Hot"] if "Immune_Hot" in imm_counts.columns else 0
            nr_hot = imm_counts.loc["Non_Responder", "Immune_Hot"] if "Immune_Hot" in imm_counts.columns else 0
            r_desert = imm_counts.loc["Responder", "Immune_Desert"] if "Immune_Desert" in imm_counts.columns else 0
            nr_desert = imm_counts.loc["Non_Responder", "Immune_Desert"] if "Immune_Desert" in imm_counts.columns else 0
            labels = ["Hot (R)", "Hot (NR)", "Desert (R)", "Desert (NR)"]
            values = [r_hot, nr_hot, r_desert, nr_desert]
            colors = [IMMUNOTYPE_PAL["Immune_Hot"], IMMUNOTYPE_PAL["Immune_Hot"],
                      IMMUNOTYPE_PAL["Immune_Desert"], IMMUNOTYPE_PAL["Immune_Desert"]]
            alphas = [0.9, 0.4, 0.9, 0.4]
            for k, (lbl, val, col, alp) in enumerate(zip(labels, values, colors, alphas)):
                ax.bar(k, val, color=col, alpha=alp, width=0.4)
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels)
            ax.set_ylabel("Fraction of Spots")
            ax.set_title("C. Immune Hot & Desert: R vs NR", fontweight="bold")

        plt.suptitle("GSE238264: Spatial Immunotype Analysis", fontsize=13, fontweight="bold")
        plt.tight_layout()
        plt.savefig(FIG / "02b_visium_immunotype_composition.png", dpi=200, bbox_inches="tight")
        plt.close()

    # ── Figure 4: Spatial autocorrelation (Moran's I) ──────────────────
    print("[Fig 4] Moran's I spatial autocorrelation...")
    morans_rows = []
    for adata in datasets:
        sname = adata.obs["sample"].iloc[0]
        resp = adata.obs["response"].iloc[0]
        for sig in scored:
            if sig in adata.obs:
                try:
                    sq.gr.spatial_autocorr(adata, mode="moran", genes=[sig],
                                           attr="obs", n_perms=99)
                    if hasattr(adata, 'uns') and 'moranI' in adata.uns:
                        mi_df = adata.uns["moranI"]
                        if sig in mi_df.index:
                            row = mi_df.loc[sig]
                            morans_rows.append({"sample": sname, "response": resp,
                                                "signature": sig, "I": row.get("I", np.nan),
                                                "pval_norm": row.get("pval_norm", np.nan)})
                except Exception:
                    pass

    if morans_rows:
        morans_df = pd.DataFrame(morans_rows)
        morans_df.to_csv(TAB / "gse238264_morans_I.csv", index=False)

        fig, axes = plt.subplots(1, 2, figsize=(16, 6))

        # Heatmap of Moran's I
        ax = axes[0]
        pivot = morans_df.pivot_table(index="signature", columns="sample", values="I")
        sns.heatmap(pivot, cmap="YlOrRd", annot=True, fmt=".2f", ax=ax,
                    linewidths=0.5, cbar_kws={"label": "Moran's I"})
        ax.set_title("A. Moran's I per Signature × Section", fontweight="bold", fontsize=10)
        ax.tick_params(labelsize=7)

        # R vs NR comparison of mean Moran's I
        ax = axes[1]
        mean_mi = morans_df.groupby(["response", "signature"])["I"].mean().unstack(level=0)
        if "Responder" in mean_mi.columns and "Non_Responder" in mean_mi.columns:
            x = np.arange(len(mean_mi))
            w = 0.35
            ax.barh(x - w/2, mean_mi["Responder"], w, color=RESPONSE_PAL["Responder"],
                    label="Responder", alpha=0.85)
            ax.barh(x + w/2, mean_mi["Non_Responder"], w, color=RESPONSE_PAL["Non_Responder"],
                    label="Non-Responder", alpha=0.85)
            ax.set_yticks(x)
            ax.set_yticklabels(mean_mi.index, fontsize=7)
            ax.set_xlabel("Mean Moran's I")
            ax.set_title("B. Mean Spatial Autocorrelation: R vs NR", fontweight="bold", fontsize=10)
            ax.legend(fontsize=8)

        plt.suptitle("GSE238264: Spatial Autocorrelation of IO Signatures",
                     fontsize=13, fontweight="bold")
        plt.tight_layout()
        plt.savefig(FIG / "02b_visium_morans_I.png", dpi=200, bbox_inches="tight")
        plt.close()
        print(f"  Moran's I table: {morans_df.shape}")
    else:
        print("  [WARN] Could not compute Moran's I")

    # ── Figure 5: Neighbourhood enrichment ─────────────────────────────
    print("[Fig 5] Neighbourhood enrichment by immunotype...")
    if "spatial_immunotype" in combined.obs:
        try:
            build_spatial_graph(combined)
            sq.gr.nhood_enrichment(combined, cluster_key="spatial_immunotype")
            fig, ax = plt.subplots(figsize=(6, 5))
            sq.pl.nhood_enrichment(combined, cluster_key="spatial_immunotype", ax=ax)
            ax.set_title("Neighbourhood Enrichment: Spatial Immunotypes", fontweight="bold")
            plt.tight_layout()
            plt.savefig(FIG / "02b_visium_nhood_enrichment.png", dpi=200, bbox_inches="tight")
            plt.close()
        except Exception as e:
            print(f"  [WARN] Nhood enrichment failed: {e}")

    # ── Figure 6: Cell-type marker score spatial maps ──────────────────
    print("[Fig 6] Cell-type marker spatial maps (representative R & NR)...")
    ct_sigs = ["ct_CD8_T", "ct_Macrophage", "ct_Tumor", "ct_B_cell"]
    r_sample = datasets[0]  # first R
    nr_sample = datasets[-1]  # last NR

    fig, axes = plt.subplots(2, len(ct_sigs), figsize=(4 * len(ct_sigs), 7))
    for row_idx, (adata, label) in enumerate([(r_sample, "Responder"), (nr_sample, "Non-Responder")]):
        score_cell_markers(adata)
        for col_idx, ct in enumerate(ct_sigs):
            ax = axes[row_idx, col_idx]
            if ct in adata.obs:
                coords = adata.obsm["spatial"]
                sc_plot = ax.scatter(coords[:, 0], coords[:, 1],
                                     c=adata.obs[ct].values, cmap="magma", s=6,
                                     edgecolors="none", rasterized=True)
                plt.colorbar(sc_plot, ax=ax, shrink=0.6, pad=0.02)
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.set_xticks([]); ax.set_yticks([])
            if row_idx == 0:
                ax.set_title(ct.replace("ct_", ""), fontsize=10, fontweight="bold")
            if col_idx == 0:
                sname = adata.obs["sample"].iloc[0]
                ax.set_ylabel(f"{sname}\n({label})", fontsize=9, fontweight="bold")

    plt.suptitle("GSE238264: Cell-Type Marker Scores — Responder vs Non-Responder",
                 fontsize=12, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(FIG / "02b_visium_celltype_spatial_maps.png", dpi=200, bbox_inches="tight")
    plt.close()

    # ── Figure 7: Signature correlation heatmap (all spots) ────────────
    print("[Fig 7] Signature correlation heatmap...")
    sig_cols = [s for s in scored if s in combined.obs.columns]
    if sig_cols:
        corr = combined.obs[sig_cols].corr(method="spearman")
        fig, ax = plt.subplots(figsize=(10, 9))
        sns.heatmap(corr, cmap="RdBu_r", center=0, square=True, ax=ax,
                    annot=True, fmt=".2f", linewidths=0.3,
                    xticklabels=[s.replace("_", " ") for s in sig_cols],
                    yticklabels=[s.replace("_", " ") for s in sig_cols],
                    cbar_kws={"label": "Spearman ρ", "shrink": 0.8})
        ax.set_title("IO Signature Correlations Across All Visium Spots (GSE238264)",
                     fontweight="bold", fontsize=11)
        ax.tick_params(labelsize=7)
        plt.tight_layout()
        plt.savefig(FIG / "02b_visium_signature_correlation.png", dpi=200, bbox_inches="tight")
        plt.close()

    # ── Figure 8: Sample-level mean signature heatmap ──────────────────
    print("[Fig 8] Sample-level signature heatmap...")
    if sig_cols:
        sample_means = combined.obs.groupby("sample")[sig_cols].mean()
        sample_means["response"] = [combined.obs.loc[combined.obs["sample"] == s, "response"].iloc[0]
                                    for s in sample_means.index]
        sample_means = sample_means.sort_values("response")
        resp_labels = sample_means.pop("response")

        fig, ax = plt.subplots(figsize=(12, 5))
        sns.heatmap(sample_means.T, cmap="RdBu_r", center=0, ax=ax,
                    annot=True, fmt=".2f", linewidths=0.5,
                    yticklabels=[s.replace("_", " ") for s in sig_cols],
                    cbar_kws={"label": "Mean Score", "shrink": 0.8})
        # Color sample labels by response
        for tick_label in ax.get_xticklabels():
            sname = tick_label.get_text()
            if sname in resp_labels.index:
                tick_label.set_color(RESPONSE_PAL.get(resp_labels[sname], "black"))
                tick_label.set_fontweight("bold")
        ax.set_title("Mean IO Signature per Section (green=R, red=NR)", fontweight="bold")
        ax.tick_params(labelsize=8)
        plt.tight_layout()
        plt.savefig(FIG / "02b_visium_sample_signature_heatmap.png", dpi=200, bbox_inches="tight")
        plt.close()

    return combined, stat_df


# ═══════════════════════════════════════════════════════════════════════
#  PART 3 — GSE281759: Multi-Region Spatial Analysis
# ═══════════════════════════════════════════════════════════════════════

def analyse_gse281759(datasets):
    """Multi-region HCC Visium: Normal vs Tumor, Primary vs Recurrence."""
    print("\n" + "=" * 70)
    print("PART 3: GSE281759 — Multi-Region HCC Spatial Analysis")
    print("=" * 70)

    for adata in datasets:
        score_io_signatures(adata)
        score_cell_markers(adata)
        classify_immunotype(adata)

    combined = sc.concat(datasets, join="inner")
    combined.obs_names_make_unique()
    scored = score_io_signatures(combined)
    classify_immunotype(combined)

    print(f"  Combined: {combined.n_obs} spots × {combined.n_vars} genes")
    print(f"  Tissue types: {dict(combined.obs['tissue_type'].value_counts())}")

    sig_cols = [s for s in scored if s in combined.obs.columns]

    # ── Figure 9: Tumor vs Normal signature comparison ─────────────────
    print("[Fig 9] Tumor vs Normal signatures...")
    tumor = combined.obs[combined.obs["tissue_type"].isin(["Tumor", "Primary", "Recurrence"])]
    normal = combined.obs[combined.obs["tissue_type"] == "Normal"]

    if len(normal) > 0 and len(tumor) > 0:
        tn_stats = []
        for sig in sig_cols:
            u, p = mannwhitneyu(tumor[sig].dropna(), normal[sig].dropna(), alternative="two-sided")
            tn_stats.append({"signature": sig,
                             "mean_Tumor": tumor[sig].mean(), "mean_Normal": normal[sig].mean(),
                             "effect": tumor[sig].mean() - normal[sig].mean(),
                             "U": u, "p_value": p})
        tn_df = pd.DataFrame(tn_stats).sort_values("p_value")
        tn_df.to_csv(TAB / "gse281759_tumor_vs_normal_stats.csv", index=False)

        fig, axes = plt.subplots(4, 4, figsize=(20, 16))
        axes_flat = axes.flatten()
        for idx, sig in enumerate(sig_cols[:16]):
            ax = axes_flat[idx]
            data_plot = pd.DataFrame({
                "score": pd.concat([tumor[sig], normal[sig]]),
                "group": ["Tumor"] * len(tumor) + ["Normal"] * len(normal),
            })
            sns.violinplot(data=data_plot, x="group", y="score", ax=ax,
                           palette=TISSUE_PAL, inner="box", cut=0, linewidth=0.8)
            row = tn_df[tn_df["signature"] == sig].iloc[0]
            stars = "***" if row.p_value < 0.001 else ("**" if row.p_value < 0.01 else
                    ("*" if row.p_value < 0.05 else "ns"))
            ax.set_title(f"{sig.replace('_', ' ')}\n{stars}", fontsize=8, fontweight="bold")
            ax.set_xlabel("")

        for idx in range(len(sig_cols[:16]), 16):
            axes_flat[idx].axis("off")

        plt.suptitle("GSE281759: IO Signatures — Tumor vs Normal (spot-level)",
                     fontsize=13, fontweight="bold")
        plt.tight_layout()
        plt.savefig(FIG / "02b_visium_gse281759_tumor_vs_normal.png", dpi=200, bbox_inches="tight")
        plt.close()

    # ── Figure 10: Primary vs Recurrence ───────────────────────────────
    print("[Fig 10] Primary vs Recurrence...")
    primary = combined.obs[combined.obs["tissue_type"] == "Primary"]
    recurrence = combined.obs[combined.obs["tissue_type"] == "Recurrence"]

    if len(primary) > 0 and len(recurrence) > 0:
        pr_stats = []
        for sig in sig_cols:
            u, p = mannwhitneyu(primary[sig].dropna(), recurrence[sig].dropna(), alternative="two-sided")
            pr_stats.append({"signature": sig,
                             "mean_Primary": primary[sig].mean(),
                             "mean_Recurrence": recurrence[sig].mean(),
                             "effect": primary[sig].mean() - recurrence[sig].mean(),
                             "U": u, "p_value": p})
        pr_df = pd.DataFrame(pr_stats).sort_values("p_value")
        pr_df.to_csv(TAB / "gse281759_primary_vs_recurrence_stats.csv", index=False)

        fig, ax = plt.subplots(figsize=(12, 6))
        plot_data = []
        for sig in sig_cols:
            for _, row in primary[[sig]].dropna().iterrows():
                plot_data.append({"Signature": sig.replace("_", " "), "Score": row[sig], "Group": "Primary"})
            for _, row in recurrence[[sig]].dropna().iterrows():
                plot_data.append({"Signature": sig.replace("_", " "), "Score": row[sig], "Group": "Recurrence"})
        plot_df = pd.DataFrame(plot_data)
        # Subsample for plotting speed
        if len(plot_df) > 50000:
            plot_df = plot_df.sample(50000, random_state=42)

        sns.boxplot(data=plot_df, x="Signature", y="Score", hue="Group",
                    palette={"Primary": "#CC3311", "Recurrence": "#9b59b6"},
                    ax=ax, fliersize=0, linewidth=0.8)
        ax.set_title("GSE281759: IO Signatures — Primary vs Recurrence", fontweight="bold")
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(FIG / "02b_visium_gse281759_primary_vs_recurrence.png", dpi=200, bbox_inches="tight")
        plt.close()

    # ── Figure 11: Per-patient signature heatmap ───────────────────────
    print("[Fig 11] Per-patient signature heatmap...")
    patient_means = combined.obs.groupby(["patient", "tissue_type"])[sig_cols].mean()
    patient_means = patient_means.reset_index()
    patient_means["label"] = patient_means["patient"] + "_" + patient_means["tissue_type"]
    patient_means = patient_means.set_index("label")

    fig, ax = plt.subplots(figsize=(14, 6))
    sns.heatmap(patient_means[sig_cols].T, cmap="RdBu_r", center=0, ax=ax,
                annot=True, fmt=".2f", linewidths=0.3,
                yticklabels=[s.replace("_", " ") for s in sig_cols],
                cbar_kws={"label": "Mean Score", "shrink": 0.8})
    ax.set_title("GSE281759: Mean IO Signatures per Patient × Tissue", fontweight="bold")
    ax.tick_params(labelsize=7)
    plt.tight_layout()
    plt.savefig(FIG / "02b_visium_gse281759_patient_heatmap.png", dpi=200, bbox_inches="tight")
    plt.close()

    # ── Figure 12: Spatial maps for representative sections ────────────
    print("[Fig 12] Spatial maps for representative sections...")
    # Pick one Normal and one Tumor section
    repr_sections = []
    for adata in datasets:
        tt = adata.obs["tissue_type"].iloc[0]
        if tt == "Normal" and not any(r[1] == "Normal" for r in repr_sections):
            repr_sections.append((adata, "Normal"))
        elif tt == "Tumor" and not any(r[1] == "Tumor" for r in repr_sections):
            repr_sections.append((adata, "Tumor"))
        if len(repr_sections) == 2:
            break

    if repr_sections:
        show_sigs = ["T_cell_cytotoxicity", "TLS_signature", "TGF_beta", "ct_Tumor"]
        show_sigs = [s for s in show_sigs if s in combined.obs.columns or
                     s in repr_sections[0][0].obs.columns]

        fig, axes = plt.subplots(len(repr_sections), len(show_sigs),
                                 figsize=(4 * len(show_sigs), 3.5 * len(repr_sections)))
        if len(repr_sections) == 1:
            axes = axes[np.newaxis, :]

        for i, (adata, tt) in enumerate(repr_sections):
            for j, sig in enumerate(show_sigs):
                ax = axes[i, j]
                if sig in adata.obs:
                    coords = adata.obsm["spatial"]
                    sp = ax.scatter(coords[:, 0], coords[:, 1],
                                    c=adata.obs[sig].values, cmap="inferno", s=8,
                                    edgecolors="none", rasterized=True)
                    plt.colorbar(sp, ax=ax, shrink=0.6, pad=0.02)
                ax.set_aspect("equal"); ax.invert_yaxis()
                ax.set_xticks([]); ax.set_yticks([])
                if i == 0:
                    ax.set_title(sig.replace("_", " ").replace("ct ", ""), fontsize=9, fontweight="bold")
                if j == 0:
                    sname = adata.obs["sample"].iloc[0]
                    ax.set_ylabel(f"{sname} ({tt})", fontsize=9, fontweight="bold")

        plt.suptitle("GSE281759: Representative Spatial Maps — Normal vs Tumor",
                     fontsize=12, fontweight="bold", y=1.01)
        plt.tight_layout()
        plt.savefig(FIG / "02b_visium_gse281759_spatial_maps.png", dpi=200, bbox_inches="tight")
        plt.close()

    return combined


# ═══════════════════════════════════════════════════════════════════════
#  PART 4 — Cross-Dataset Summary
# ═══════════════════════════════════════════════════════════════════════

def cross_dataset_summary(gse238_combined, gse281_combined):
    """Final cross-dataset summary figure."""
    print("\n" + "=" * 70)
    print("PART 4: Cross-Dataset Summary")
    print("=" * 70)

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Panel A: total spots summary
    ax = axes[0, 0]
    datasets_info = {
        "GSE238264\n(IO Response)": gse238_combined.n_obs,
        "GSE281759\n(Multi-Region)": gse281_combined.n_obs,
    }
    bars = ax.bar(datasets_info.keys(), datasets_info.values(),
                  color=["#0077BB", "#3498db"], width=0.4, edgecolor="black", linewidth=0.5)
    for bar, val in zip(bars, datasets_info.values()):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 100,
                f"{val:,}", ha="center", fontweight="bold")
    ax.set_ylabel("Number of Visium Spots")
    ax.set_title("A. Dataset Size", fontweight="bold")

    # Panel B: R vs NR key signature effect sizes (GSE238264)
    ax = axes[0, 1]
    eff_file = TAB / "gse238264_RvsNR_signature_stats.csv"
    if eff_file.exists():
        eff_df = pd.read_csv(eff_file).sort_values("effect_size")
        colors = ["#0077BB" if x > 0 else "#CC3311" for x in eff_df["effect_size"]]
        ax.barh(eff_df["signature"].str.replace("_", " "), eff_df["effect_size"],
                color=colors, edgecolor="black", linewidth=0.3, height=0.5)
        ax.axvline(0, color="black", linewidth=0.5)
        ax.set_xlabel("Effect Size (R − NR)")
        ax.set_title("B. IO Signature Effect Sizes (R vs NR)", fontweight="bold")
        ax.tick_params(labelsize=7)

    # Panel C: Immunotype fractions summary
    ax = axes[1, 0]
    if "spatial_immunotype" in gse238_combined.obs:
        fracs_238 = gse238_combined.obs.groupby("response")["spatial_immunotype"].value_counts(normalize=True)
        fracs_238 = fracs_238.unstack(fill_value=0)
        fracs_238.plot(kind="bar", ax=ax,
                       color=[IMMUNOTYPE_PAL.get(c, "#666") for c in fracs_238.columns],
                       width=0.5, edgecolor="black", linewidth=0.3)
        ax.set_title("C. Spatial Immunotype: R vs NR (GSE238264)", fontweight="bold")
        ax.set_ylabel("Fraction of Spots"); ax.tick_params(axis="x", rotation=0)
        ax.legend(fontsize=7)

    # Panel D: Tumor vs Normal key signatures (GSE281759)
    ax = axes[1, 1]
    tn_file = TAB / "gse281759_tumor_vs_normal_stats.csv"
    if tn_file.exists():
        tn_df = pd.read_csv(tn_file).sort_values("effect")
        colors = ["#CC3311" if x > 0 else "#3498db" for x in tn_df["effect"]]
        ax.barh(tn_df["signature"].str.replace("_", " "), tn_df["effect"],
                color=colors, edgecolor="black", linewidth=0.3, height=0.5)
        ax.axvline(0, color="black", linewidth=0.5)
        ax.set_xlabel("Effect Size (Tumor − Normal)")
        ax.set_title("D. IO Signature Effects: Tumor vs Normal (GSE281759)", fontweight="bold")
        ax.tick_params(labelsize=7)

    plt.suptitle("Visium Spatial Transcriptomics — Cross-Dataset Summary\n"
                 f"GSE238264: {gse238_combined.n_obs:,} spots | GSE281759: {gse281_combined.n_obs:,} spots",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG / "02b_visium_cross_dataset_summary.png", dpi=200, bbox_inches="tight")
    plt.close()


# ═══════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("02b — True Visium Spatial Analysis for HCC Immunotherapy Response")
    print("=" * 70)

    # ── Load datasets ──────────────────────────────────────────────────
    print("\n[1] Loading GSE238264 (HCC IO Response Visium)...")
    gse238 = load_gse238264()

    print(f"\n[2] Loading GSE281759 (HCC Multi-Region Visium)...")
    gse281 = load_gse281759()

    # ── Run analyses ───────────────────────────────────────────────────
    gse238_combined, gse238_stats = analyse_gse238264(gse238)
    gse281_combined = analyse_gse281759(gse281)

    # ── Cross-dataset summary ──────────────────────────────────────────
    cross_dataset_summary(gse238_combined, gse281_combined)

    # ── Print summary ─────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("[DONE] Visium Spatial Analysis Complete")
    print(f"  GSE238264: {gse238_combined.n_obs:,} spots (4R + 3NR)")
    print(f"  GSE281759: {gse281_combined.n_obs:,} spots (17 multi-region)")
    print(f"  Figures: {FIG}")
    print(f"  Tables:  {TAB}")
    print("=" * 70)


if __name__ == "__main__":
    main()
