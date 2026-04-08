#!/usr/bin/env python3
"""02h – CNV ↔ TME ↔ Response integration dashboard.

Uses real data from the integrated atlas:
  A. Spearman correlation heatmap between cnv_leiden clones and immune cell fractions
  B. Volcano plot of DEGs in high-CNV (malignant) vs normal cells
  C. Immune fraction bar plot by responder status
  D. Kaplan-Meier-style survival plot (if clinical data available, else
     CNV-score boxplot split by response)
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
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu

PROJECT = Path(__file__).resolve().parents[2]
ATLAS   = PROJECT / "data" / "processed" / "integrated_io_atlas.h5ad"
OUT     = PROJECT / "out" / "02_analysis" / "figures"
OUT.mkdir(parents=True, exist_ok=True)


def _resolve_celltype_col(adata):
    for col in ("cell_type", "celltype_sub", "Cell.type", "celltype_global", "celltype", "leiden"):
        if col in adata.obs.columns:
            return col
    raise KeyError("No cell-type column found")


def run_integration(adata_path=None):
    if adata_path is None:
        adata_path = str(ATLAS)

    print("Loading atlas …")
    adata = sc.read_h5ad(adata_path)
    print(f"  Shape: {adata.shape}")

    ct_col = _resolve_celltype_col(adata)
    has_cnv = "cnv_score" in adata.obs.columns and "cnv_leiden" in adata.obs.columns
    has_malig = "malignant_status" in adata.obs.columns
    has_response = "response" in adata.obs.columns

    # If InferCNV hasn't been run yet, compute cnv_score from 02g first
    if not has_cnv:
        print("  cnv_score / cnv_leiden not found – running InferCNV first …")
        try:
            import infercnvpy as cnv_mod
            from scripts.analysis_02g import _pick_reference_cats  # noqa
        except ImportError:
            pass
        # Try to load from a saved version
        cnv_path = PROJECT / "data" / "processed" / "integrated_io_atlas_cnv.h5ad"
        if cnv_path.exists():
            adata = sc.read_h5ad(str(cnv_path))
            has_cnv = True
            has_malig = "malignant_status" in adata.obs.columns

    if not has_cnv:
        raise RuntimeError("Run 02g_infercnv.py first to generate cnv_score and cnv_leiden columns")

    # ── aesthetics ───────────────────────────────────────────────────
    plt.style.use("seaborn-v0_8-white")
    font_title = {"fontweight": "bold", "fontsize": 16}
    font_label = {"fontweight": "bold", "fontsize": 12}
    hc_cmap = LinearSegmentedColormap.from_list("hc_diverge", ["#00E5FF", "#111111", "#FF004D"])

    fig = plt.figure(figsize=(20, 16))
    gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.30)

    # ── A. Clone ↔ Immune Correlation ────────────────────────────────
    axA = fig.add_subplot(gs[0, 0])
    print("  Panel A: Clone × immune-cell correlation …")

    # Build per-patient fraction table
    if "patient_ID" in adata.obs.columns:
        pid_col = "patient_ID"
    elif "Donor" in adata.obs.columns:
        pid_col = "Donor"
    else:
        pid_col = "sample_id"

    obs = adata.obs[[pid_col, ct_col, "cnv_leiden"]].copy()
    obs["cnv_leiden"] = obs["cnv_leiden"].astype(str)

    # Cell-type fractions per patient
    ct_frac = obs.groupby([pid_col, ct_col]).size().unstack(fill_value=0)
    ct_frac = ct_frac.div(ct_frac.sum(axis=1), axis=0)

    # Clone fractions per patient
    clone_frac = obs.groupby([pid_col, "cnv_leiden"]).size().unstack(fill_value=0)
    clone_frac = clone_frac.div(clone_frac.sum(axis=1), axis=0)

    # Keep top clones (max 6) and top cell types (max 8) for readability
    top_clones = clone_frac.mean().nlargest(6).index.tolist()
    top_cts = ct_frac.mean().nlargest(8).index.tolist()

    shared_patients = ct_frac.index.intersection(clone_frac.index)
    corr_matrix = pd.DataFrame(index=top_clones, columns=top_cts, dtype=float)
    for cl in top_clones:
        for ct in top_cts:
            rho, _ = spearmanr(clone_frac.loc[shared_patients, cl],
                               ct_frac.loc[shared_patients, ct])
            corr_matrix.loc[cl, ct] = rho

    corr_matrix = corr_matrix.astype(float)
    sns.heatmap(corr_matrix, annot=True, fmt=".2f", cmap=hc_cmap, center=0,
                ax=axA, cbar_kws={"label": "Spearman ρ", "shrink": 0.7},
                linewidths=1, linecolor="white")
    axA.set_title("A. TME Immune Enrichment by CNV Clone", **font_title, pad=15)
    axA.set_xticklabels(axA.get_xticklabels(), rotation=45, ha="right", fontweight="bold", fontsize=9)
    axA.set_yticklabels([f"Clone {c}" for c in top_clones], fontweight="bold", fontsize=10)

    # ── B. Volcano: Malignant vs Normal DEGs ─────────────────────────
    axB = fig.add_subplot(gs[0, 1])
    axB.set_facecolor("#1A1A1A")
    print("  Panel B: DE analysis (malignant vs normal) …")

    if has_malig:
        adata_de = adata[adata.obs["malignant_status"].isin(["Malignant", "Normal"])].copy()
        sc.tl.rank_genes_groups(adata_de, groupby="malignant_status", reference="Normal",
                                method="wilcoxon", n_genes=adata_de.n_vars)
        result = sc.get.rank_genes_groups_df(adata_de, group="Malignant")
        result["neg_log10_pval"] = -np.log10(result["pvals_adj"].clip(lower=1e-300))

        sig_up = (result["logfoldchanges"] > 1) & (result["neg_log10_pval"] > 2)
        sig_dn = (result["logfoldchanges"] < -1) & (result["neg_log10_pval"] > 2)
        ns = ~(sig_up | sig_dn)

        axB.scatter(result.loc[ns, "logfoldchanges"], result.loc[ns, "neg_log10_pval"],
                    color="#525252", alpha=0.5, s=12, label="NS")
        axB.scatter(result.loc[sig_dn, "logfoldchanges"], result.loc[sig_dn, "neg_log10_pval"],
                    color="#00E5FF", alpha=0.9, s=30, label="Down")
        axB.scatter(result.loc[sig_up, "logfoldchanges"], result.loc[sig_up, "neg_log10_pval"],
                    color="#FF004D", alpha=0.9, s=30, label="Up")

        # Label top genes
        for direction, subset in [("up", sig_up), ("down", sig_dn)]:
            top = result.loc[subset].nlargest(5, "neg_log10_pval")
            for _, row in top.iterrows():
                col = "#FF004D" if direction == "up" else "#00E5FF"
                axB.annotate(row["names"], (row["logfoldchanges"], row["neg_log10_pval"]),
                             color="white", fontweight="bold", fontsize=8,
                             xytext=(5, 5), textcoords="offset points",
                             bbox=dict(boxstyle="round,pad=0.2", fc="#333", ec=col, lw=1))

        axB.axvline(1, color="white", ls="--", lw=0.8, alpha=0.5)
        axB.axvline(-1, color="white", ls="--", lw=0.8, alpha=0.5)
        axB.axhline(2, color="white", ls="--", lw=0.8, alpha=0.5)
    else:
        axB.text(0.5, 0.5, "Run 02g first for\nmalignant_status", ha="center", va="center",
                 transform=axB.transAxes, color="white", fontsize=14)

    axB.set_title("B. Malignant vs Normal DEGs", **font_title, pad=15)
    axB.set_xlabel("Log2 Fold Change", **font_label)
    axB.set_ylabel("-Log10(p-adj)", **font_label)
    axB.legend(loc="upper right", facecolor="#1A1A1A", edgecolor="white", labelcolor="white")

    # ── C. Immune Fraction by Response ───────────────────────────────
    axC = fig.add_subplot(gs[1, 0])
    print("  Panel C: Immune fractions by response …")

    if has_response:
        resp_ct = adata.obs.groupby(["response", ct_col]).size().unstack(fill_value=0)
        resp_ct = resp_ct.div(resp_ct.sum(axis=1), axis=0)
        top_cts_bar = resp_ct.mean().nlargest(6).index.tolist()
        resp_ct_top = resp_ct[top_cts_bar]

        palette = ["#00E5FF", "#FF004D", "#FFD700", "#39FF14", "#FF6EC7", "#A855F7"]
        x = np.arange(len(resp_ct_top.index))
        width = 0.8 / len(top_cts_bar)
        for j, ct_name in enumerate(top_cts_bar):
            axC.bar(x + j * width, resp_ct_top[ct_name], width,
                    label=ct_name, color=palette[j % len(palette)], edgecolor="black", lw=0.8)
        axC.set_xticks(x + width * len(top_cts_bar) / 2)
        axC.set_xticklabels(resp_ct_top.index, **font_label)
    axC.set_ylabel("Cell Type Fraction", **font_label)
    axC.set_title("C. Immune Ecosystem by Response", **font_title, pad=15)
    axC.legend(title="Cell Type", fontsize=8, title_fontproperties={"weight": "bold"},
               bbox_to_anchor=(1.0, 1), loc="upper left")
    axC.grid(axis="y", ls="--", alpha=0.5)

    # ── D. CNV Score by Response ─────────────────────────────────────
    axD = fig.add_subplot(gs[1, 1])
    axD.set_facecolor("#F9F9F9")
    print("  Panel D: CNV score distribution by response …")

    if has_response and has_cnv:
        resp_groups = adata.obs["response"].dropna().unique()
        data_d = []
        labels_d = []
        for rg in sorted(resp_groups):
            vals = adata.obs.loc[adata.obs["response"] == rg, "cnv_score"].dropna().values
            data_d.append(vals)
            labels_d.append(rg)

        colors_d = ["#2ECC71", "#E74C3C"][:len(data_d)]
        parts = axD.violinplot(data_d, showmeans=False, showmedians=False, showextrema=False)
        for pc, c in zip(parts["bodies"], colors_d):
            pc.set_facecolor(c); pc.set_edgecolor("black"); pc.set_alpha(0.85)
        for i, d in enumerate(data_d):
            q1, med, q3 = np.percentile(d, [25, 50, 75])
            axD.scatter(i + 1, med, marker="o", color="white", s=30, zorder=3)
            axD.vlines(i + 1, q1, q3, color="#333", lw=5)
            axD.vlines(i + 1, d.min(), d.max(), color="#333", lw=1)
        axD.set_xticks(range(1, len(labels_d) + 1))
        axD.set_xticklabels(labels_d, **font_label)

        # Statistical test
        if len(data_d) == 2:
            stat, pval = mannwhitneyu(data_d[0], data_d[1], alternative="two-sided")
            axD.text(0.5, 0.95, f"Mann-Whitney p = {pval:.2e}", transform=axD.transAxes,
                     ha="center", fontsize=12, fontweight="bold",
                     bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.4"))

    axD.set_ylabel("CNV Score", **font_label)
    axD.set_title("D. CNV Burden by IO Response", **font_title, pad=15)
    axD.grid(True, ls=":", alpha=0.6)

    plt.tight_layout()
    out_path = OUT / "02h_cnv_tme_high_contrast_dashboard.png"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {out_path}")


if __name__ == "__main__":
    run_integration()
