#!/usr/bin/env python3
"""
06_compile_extended_data.py — Extended Data figures that capture
core-story analyses not shown in the six main figures.

  Extended Data Fig. 1 — Genomic instability (CNV / InferCNV / CNV–TME)
  Extended Data Fig. 2 — Regulatory networks & cell-cell communication
  Extended Data Fig. 3 — High-contrast spatial niche maps
  Extended Data Fig. 4 — Clinical-trial bulk validation & deconvolution

All figures saved at 300 DPI as PNG + PDF in out/06_manuscript_figures/.
"""
import warnings
warnings.filterwarnings("ignore")

import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.image import imread
from matplotlib.gridspec import GridSpec
from config.publication_style import apply_style

apply_style()

OUT = PROJECT / "out" / "06_manuscript_figures"
OUT.mkdir(parents=True, exist_ok=True)

P = {
    # CNV
    "cnv_dashboard":   PROJECT / "out" / "02_analysis" / "figures" / "02g_infercnv_attractive_dashboard.png",
    "cnv_heatmap":     PROJECT / "out" / "02_analysis" / "figures" / "02g_infercnv_chromosome_heatmap.png",
    "cnv_umap":        PROJECT / "out" / "02_analysis" / "figures" / "02g_cnv_score_umap.png",
    "cnv_tme":         PROJECT / "out" / "02_analysis" / "figures" / "02h_cnv_tme_high_contrast_dashboard.png",
    # Regulon + cell-cell
    "regulon":         PROJECT / "out" / "02a_response_tme" / "figures" / "02f_tf_regulon_activity.png",
    "ccc_network":     PROJECT / "out" / "04_complementary" / "figures" / "04b_cell_communication_network.png",
    "lr_dotplot":      PROJECT / "out" / "04_complementary" / "figures" / "04d_ccc_lr_dotplot.png",
    "lr_heatmap":      PROJECT / "out" / "04_complementary" / "figures" / "04_ligand_receptor.png",
    # High-contrast spatial
    "hcs_dashboard":   PROJECT / "out" / "04_complementary" / "figures" / "04h_spatial_high_contrast_maps.png",
    "hcs_tissue":      PROJECT / "out" / "04_complementary" / "figures" / "04h_tissue_stratified_spatial.png",
    "hcs_hypoxia":     PROJECT / "out" / "04_complementary" / "figures" / "spatial_1_hypoxia_gradient.png",
    "hcs_cd8exh":      PROJECT / "out" / "04_complementary" / "figures" / "spatial_2_cd8_exhaustion_hubs.png",
    "hcs_macpol":      PROJECT / "out" / "04_complementary" / "figures" / "spatial_3_macrophage_polarization.png",
    "hcs_iopred":      PROJECT / "out" / "04_complementary" / "figures" / "spatial_4_io_response_prediction.png",
    # Bulk trial validation
    "trial_surv":      PROJECT / "out" / "04_complementary" / "figures" / "04b_bulk_trial_survival_validation.png",
    "bulk_io":         PROJECT / "out" / "02c_bulk_deconvolution" / "figures" / "bulk_io_heatmap.png",
    "ref_heatmap":     PROJECT / "out" / "02c_bulk_deconvolution" / "figures" / "02c_reference_heatmap.png",
    "ref_corr":        PROJECT / "out" / "02c_bulk_deconvolution" / "figures" / "02c_reference_correlation.png",
    "ref_counts":      PROJECT / "out" / "02c_bulk_deconvolution" / "figures" / "02c_reference_cellcounts.png",
    "io_pub":          PROJECT / "out" / "03_validation" / "figures" / "io_biomarker_public.png",
    # Cappuyns external validation
    "cap_val_roc":     PROJECT / "out" / "03_validation" / "figures" / "03c_validation_roc.png",
    "cap_val_water":   PROJECT / "out" / "03_validation" / "figures" / "03c_validation_waterfall.png",
    "cap_val_layers":  PROJECT / "out" / "03_validation" / "figures" / "03c_validation_perlayer.png",
}


def _place(ax, key, label=None, fontsize_label=14):
    path = P.get(key)
    if path and path.exists():
        img = imread(str(path)); ax.imshow(img)
    else:
        ax.text(0.5, 0.5, f"[{key}]\nnot found", ha="center", va="center",
                fontsize=8, color="grey", transform=ax.transAxes)
    ax.set_axis_off()
    if label:
        ax.text(-0.02, 1.02, label, transform=ax.transAxes,
                fontsize=fontsize_label, fontweight="bold",
                va="bottom", ha="right")


def _save(fig, stem):
    for fmt in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{fmt}", dpi=300,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {stem}  saved.")


# ── Extended Data Fig. 1 — Genomic instability ────────────────────────
def ed1():
    fig = plt.figure(figsize=(20, 14), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.1])
    _place(fig.add_subplot(gs[0, 0]), "cnv_umap",      "A")
    _place(fig.add_subplot(gs[0, 1]), "cnv_dashboard", "B")
    _place(fig.add_subplot(gs[1, 0]), "cnv_heatmap",   "C")
    _place(fig.add_subplot(gs[1, 1]), "cnv_tme",       "D")
    fig.suptitle(
        "Extended Data Fig. 1 — Copy-number variation, malignant clones, and CNV–TME cross-modal integration",
        fontsize=14, fontweight="bold", y=1.01)
    _save(fig, "ExtendedData_Fig1_CNV_genomic_instability")


# ── Extended Data Fig. 2 — Regulatory & cell-cell communication ──────
def ed2():
    fig = plt.figure(figsize=(20, 14), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig)
    _place(fig.add_subplot(gs[0, 0]), "regulon",     "A")
    _place(fig.add_subplot(gs[0, 1]), "ccc_network", "B")
    _place(fig.add_subplot(gs[1, 0]), "lr_dotplot",  "C")
    _place(fig.add_subplot(gs[1, 1]), "lr_heatmap",  "D")
    fig.suptitle(
        "Extended Data Fig. 2 — TF regulons (SCENIC) and cell-cell communication (CellPhoneDB)",
        fontsize=14, fontweight="bold", y=1.01)
    _save(fig, "ExtendedData_Fig2_regulon_communication")


# ── Extended Data Fig. 3 — High-contrast spatial niches ──────────────
def ed3():
    fig = plt.figure(figsize=(22, 14), constrained_layout=True)
    gs = GridSpec(2, 3, figure=fig)
    _place(fig.add_subplot(gs[0, 0]), "hcs_hypoxia", "A")
    _place(fig.add_subplot(gs[0, 1]), "hcs_cd8exh",  "B")
    _place(fig.add_subplot(gs[0, 2]), "hcs_macpol",  "C")
    _place(fig.add_subplot(gs[1, 0]), "hcs_iopred",  "D")
    _place(fig.add_subplot(gs[1, 1]), "hcs_dashboard","E")
    _place(fig.add_subplot(gs[1, 2]), "hcs_tissue",  "F")
    fig.suptitle(
        "Extended Data Fig. 3 — High-contrast spatial niche maps "
        "(hypoxia, CD8 exhaustion, M1/M2 polarisation, IO-response gradient)",
        fontsize=14, fontweight="bold", y=1.01)
    _save(fig, "ExtendedData_Fig3_high_contrast_spatial")


# ── Extended Data Fig. 4 — Clinical-trial bulk validation ────────────
def ed4():
    fig = plt.figure(figsize=(20, 14), constrained_layout=True)
    gs = GridSpec(2, 3, figure=fig, height_ratios=[0.9, 1.1])
    _place(fig.add_subplot(gs[0, 0]), "ref_counts",  "A")
    _place(fig.add_subplot(gs[0, 1]), "ref_heatmap", "B")
    _place(fig.add_subplot(gs[0, 2]), "ref_corr",    "C")
    _place(fig.add_subplot(gs[1, 0]), "bulk_io",     "D")
    _place(fig.add_subplot(gs[1, 1]), "io_pub",      "E")
    _place(fig.add_subplot(gs[1, 2]), "trial_surv",  "F")
    fig.suptitle(
        "Extended Data Fig. 4 — CIBERSORTx reference & bulk clinical-trial validation "
        "(IMbrave150, CheckMate-459 proxies)",
        fontsize=14, fontweight="bold", y=1.01)
    _save(fig, "ExtendedData_Fig4_bulk_trial_validation")


# ── Extended Data Fig. 5 — Held-out Cappuyns validation ─────────────
def ed5():
    fig = plt.figure(figsize=(20, 12), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.0])
    _place(fig.add_subplot(gs[0, 0]), "cap_val_roc",    "A")
    _place(fig.add_subplot(gs[0, 1]), "cap_val_layers", "B")
    _place(fig.add_subplot(gs[1, :]), "cap_val_water",  "C")
    fig.suptitle(
        "Extended Data Fig. 5 — Fully held-out Cappuyns external validation "
        "of discovery-derived predictors",
        fontsize=14, fontweight="bold", y=1.01)
    _save(fig, "ExtendedData_Fig5_cappuyns_validation")


if __name__ == "__main__":
    print("=" * 70)
    print("Compiling Extended Data Figures")
    print("=" * 70)
    ed1()
    ed2()
    ed3()
    ed4()
    ed5()
    print("=" * 70)
    print(f"Done. Extended Data figures in: {OUT}")
    print("=" * 70)
