#!/usr/bin/env python3
"""
06_compile_main_figures.py — Assemble named manuscript figures from
individual analysis outputs, preserving native aspect ratios of panels.

Redesigned around the restructured narrative:

    Tumour-intrinsic (L1 CNV)  ×  Regulatory state (L2)
       × Immune engagement (L3-L4)  × Spatial geometry
          × NMF integrator (L5)  × Multi-layer predictor.

  Figure 1 — Integrated HCC IO atlas landscape
  Figure 2 — Tumour-intrinsic axis          (CNV, clonality)
  Figure 3 — Immune engagement axis         (composition + IO sigs)
  Figure 4 — Spatial axis                   (Visium + niches)
  Figure 5 — Cross-modality NMF integrator  (sc × bulk × spatial)
  Figure 6 — Multi-layer predictor + clinical translation

Each composite figure is saved at 300 DPI as both PNG and PDF.
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

# ── Paths to individual figure panels ──────────────────────────────────
P = {
    # 01 — preprocessing / atlas
    "umap_celltype":       PROJECT / "out" / "01_preprocessing" / "figures" / "01_celltype_annotation_umap.png",
    "umap_integration":    PROJECT / "out" / "01_preprocessing" / "figures" / "01_integration_umap.png",
    "qc_violins":          PROJECT / "out" / "01_preprocessing" / "figures" / "01_qc_violins.png",
    "marker_features":     PROJECT / "out" / "01_preprocessing" / "figures" / "01_marker_feature_plots.png",
    "immune_umaps":        PROJECT / "out" / "01_preprocessing" / "figures" / "01_immune_umaps.png",
    # 02a — TME
    "celltype_response":   PROJECT / "out" / "02a_response_tme" / "figures" / "02a_celltype_response_umap.png",
    "abundance_bar":       PROJECT / "out" / "02a_response_tme" / "figures" / "02a_abundance_barplot.png",
    "signature_violins":   PROJECT / "out" / "02a_response_tme" / "figures" / "02a_signature_violins.png",
    "volcano":             PROJECT / "out" / "02a_response_tme" / "figures" / "02a_volcano_deg.png",
    "marker_dotplot":      PROJECT / "out" / "02a_response_tme" / "figures" / "02a_marker_dotplot.png",
    "tissue_heatmap":      PROJECT / "out" / "02a_response_tme" / "figures" / "02a_tissue_signature_heatmap.png",
    "signature_radar":     PROJECT / "out" / "02a_response_tme" / "figures" / "02a_signature_radar.png",
    "composition_bar":     PROJECT / "out" / "02a_response_tme" / "figures" / "02a_composition_stacked_bar.png",
    # 02b — Spatial
    "spatial_ct_maps":     PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_ct_spatial_celltype_maps.png",
    "spatial_sig_R":       PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_visium_spatial_signature_maps_responder.png",
    "spatial_sig_NR":      PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_visium_spatial_signature_maps_non_responder.png",
    "spatial_nhood":       PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_visium_nhood_enrichment.png",
    "spatial_morans":      PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_visium_morans_I.png",
    "spatial_rvnr_violin": PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_visium_RvsNR_violins.png",
    "spatial_summary":     PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_ds_comprehensive_summary.png",
    # 02i — Cappuyns
    "cap_heatmap":         PROJECT / "out" / "02i_cappuyns_signatures" / "figures" / "cappuyns_signature_heatmap.png",
    "cap_key_sigs":        PROJECT / "out" / "02i_cappuyns_signatures" / "figures" / "cappuyns_key_signatures_response.png",
    "cap_subsets":         PROJECT / "out" / "02i_cappuyns_signatures" / "figures" / "cappuyns_molecular_subsets.png",
    "cap_immune_status":   PROJECT / "out" / "02i_cappuyns_signatures" / "figures" / "cappuyns_immune_status_comparison.png",
    "concordance":         PROJECT / "out" / "02j_concordance" / "figures" / "cappuyns_internal_concordance.png",
    "cap_spatial":         PROJECT / "out" / "02k_cappuyns_spatial" / "figures" / "cappuyns_spatial_heatmaps.png",
    "cap_spatial_comp":    PROJECT / "out" / "02k_cappuyns_spatial" / "figures" / "cappuyns_spatial_section_comparison.png",
    # 03 — Validation
    "roc_lasso":           PROJECT / "out" / "03_validation" / "figures" / "03_roc_lasso.png",
    "rf_importance":       PROJECT / "out" / "03_validation" / "figures" / "03_rf_importance.png",
    "bootstrap_ci":        PROJECT / "out" / "03_validation" / "figures" / "03_bootstrap_ci.png",
    "bulk_clustermap":     PROJECT / "out" / "03_validation" / "figures" / "03_bulk_clustermap.png",
    "pathway_enrich":      PROJECT / "out" / "03_validation" / "figures" / "03_pathway_enrichment.png",
    "survival":            PROJECT / "out" / "03_validation" / "figures" / "03b_clinical_survival_validation.png",
    # 04 — Complementary
    "trajectory":          PROJECT / "out" / "04_complementary" / "figures" / "04c_tcell_pseudotime_trajectories.png",
    "metabolic":           PROJECT / "out" / "04_complementary" / "figures" / "04f_metabolic_reprogramming.png",
    "tls":                 PROJECT / "out" / "04_complementary" / "figures" / "04g_tls_spatial_hubs.png",
    "cd8_states":          PROJECT / "out" / "04_complementary" / "figures" / "04_cd8_states.png",
    "lr_dotplot":          PROJECT / "out" / "04_complementary" / "figures" / "04d_ccc_lr_dotplot.png",
    "tcr_alluvial":        PROJECT / "out" / "04_complementary" / "figures" / "04e_tcr_clonotype_alluvial.png",
    "immune_ratios":       PROJECT / "out" / "04_complementary" / "figures" / "04_immune_ratios.png",
    "patient_subtypes":    PROJECT / "out" / "04_complementary" / "figures" / "04_patient_subtypes.png",
    # 02 — CNV / tumour-intrinsic axis
    "cnv_umap":            PROJECT / "out" / "02_analysis" / "figures" / "02g_cnv_score_umap.png",
    "cnv_dashboard":       PROJECT / "out" / "02_analysis" / "figures" / "02g_infercnv_attractive_dashboard.png",
    "cnv_chrom":           PROJECT / "out" / "02_analysis" / "figures" / "02g_infercnv_chromosome_heatmap.png",
    "cnv_tme":             PROJECT / "out" / "02_analysis" / "figures" / "02h_cnv_tme_high_contrast_dashboard.png",
    # 05 — multi-layer polished panels
    "cnv_RvsNR":           PROJECT / "out" / "05_multilayer" / "figures" / "05e_cnv_RvsNR.png",
    "nmf_forest":          PROJECT / "out" / "05_multilayer" / "figures" / "05e_nmf_program_vs_response.png",
    "nmf_topgenes":        PROJECT / "out" / "05_multilayer" / "figures" / "05e_nmf_top_genes_dot.png",
    "nmf_bulk_vs_sc":      PROJECT / "out" / "05_multilayer" / "figures" / "05e_nmf_bulk_vs_sc_concordance.png",
    "nmf_spatial":         PROJECT / "out" / "05_multilayer" / "figures" / "05e_nmf_spatial_heatmap.png",
    "layer_auc":           PROJECT / "out" / "05_multilayer" / "figures" / "05e_auc_bars.png",
    "layer_roc":           PROJECT / "out" / "05_multilayer" / "figures" / "05e_roc_by_layer.png",
    "feature_effects":     PROJECT / "out" / "05_multilayer" / "figures" / "05e_feature_effects.png",
    "variance_partition":  PROJECT / "out" / "05_multilayer" / "figures" / "05e_variance_partition.png",
}


# ── Helper ─────────────────────────────────────────────────────────────
def _place(ax, key, label=None, label_xy=(-0.02, 1.02), fontsize_label=14):
    """Place a panel image on *ax* preserving native aspect ratio."""
    path = P.get(key)
    if path and path.exists():
        img = imread(str(path))
        ax.imshow(img)  # default aspect='equal' preserves native ratio
    else:
        ax.text(0.5, 0.5, f"[{key}]\nnot found", ha="center", va="center",
                fontsize=8, color="grey", transform=ax.transAxes)
    ax.set_axis_off()
    if label:
        ax.text(label_xy[0], label_xy[1], label, transform=ax.transAxes,
                fontsize=fontsize_label, fontweight="bold",
                va="bottom", ha="right")


# Uniform canvas so every main figure shares width/height/DPI and
# can be stacked or placed side-by-side without rescaling artefacts.
MAIN_FIGSIZE = (19.0, 12.0)   # inches  → 5700 × 3600 px @ 300 dpi
MAIN_DPI     = 350            # print-quality


def _save(fig, stem):
    for fmt in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{fmt}", dpi=MAIN_DPI,
                    facecolor="white")
    plt.close(fig)
    print(f"  {stem}  saved.")


# Canonical canvas all main figures are padded to (identical pixel
# dimensions make side-by-side / stacked composition trivial downstream).
CANVAS_PX = (int(MAIN_FIGSIZE[0] * MAIN_DPI),   # width
             int(MAIN_FIGSIZE[1] * MAIN_DPI))   # height  → 6650 × 4200


def _normalize_canvas():
    """Pad every Figure*.png/pdf saved above to CANVAS_PX on a white
    background so widths and heights are identical across figures.
    PDFs are re-exported by converting the padded PNG back to PDF."""
    from PIL import Image
    for png in sorted(OUT.glob("Figure?_*.png")):
        im = Image.open(png).convert("RGB")
        # scale down if larger than canvas while preserving aspect
        w, h = im.size
        tw, th = CANVAS_PX
        scale = min(tw / w, th / h, 1.0)
        if scale < 1.0:
            im = im.resize((int(w * scale), int(h * scale)),
                           Image.LANCZOS)
        w, h = im.size
        canvas = Image.new("RGB", CANVAS_PX, "white")
        canvas.paste(im, ((tw - w) // 2, (th - h) // 2))
        canvas.save(png, dpi=(MAIN_DPI, MAIN_DPI), optimize=True)
        canvas.save(png.with_suffix(".pdf"), "PDF",
                    resolution=MAIN_DPI)
    print(f"  normalised all main figures to {CANVAS_PX[0]}\u00d7{CANVAS_PX[1]} px")


# ====================================================================
# Figure 1 — Integrated HCC IO atlas landscape
# ====================================================================
def figure1():
    fig = plt.figure(figsize=MAIN_FIGSIZE, constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 0.9])
    _place(fig.add_subplot(gs[0, 0]), "umap_celltype",    "A")
    _place(fig.add_subplot(gs[0, 1]), "umap_integration", "B")
    _place(fig.add_subplot(gs[1, 0]), "marker_dotplot",   "C")
    _place(fig.add_subplot(gs[1, 1]), "composition_bar",  "D")
    fig.suptitle(
        "Figure 1  |  Integrated HCC IO atlas  —  ~970 000 cells, "
        "5 cohorts, 24 Responder / Non-Responder patients",
        fontsize=14, fontweight="bold", y=1.01)
    _save(fig, "Figure1_atlas_landscape")


# ====================================================================
# Figure 2 — Tumour-intrinsic axis  (L1: CNV & clonal heterogeneity)
# ====================================================================
def figure2():
    fig = plt.figure(figsize=MAIN_FIGSIZE, constrained_layout=True)
    gs = GridSpec(2, 3, figure=fig, height_ratios=[0.95, 1.0])
    _place(fig.add_subplot(gs[0, 0]), "cnv_umap",      "A")
    _place(fig.add_subplot(gs[0, 1]), "cnv_chrom",     "B")
    _place(fig.add_subplot(gs[0, 2]), "cnv_dashboard", "C")
    _place(fig.add_subplot(gs[1, :2]), "cnv_RvsNR",    "D")
    _place(fig.add_subplot(gs[1, 2]),  "cnv_tme",      "E")
    fig.suptitle(
        "Figure 2  |  Tumour-intrinsic axis  —  "
        "CNV burden is similar between R and NR, but clonal "
        "heterogeneity (SD) is elevated in Non-Responders (p=0.028)",
        fontsize=13, fontweight="bold", y=1.01)
    _save(fig, "Figure2_tumour_intrinsic")


# ====================================================================
# Figure 3 — Immune engagement axis  (L3 composition + L4 IO sigs)
# ====================================================================
def figure3():
    fig = plt.figure(figsize=MAIN_FIGSIZE, constrained_layout=True)
    gs = GridSpec(2, 3, figure=fig)
    _place(fig.add_subplot(gs[0, 0]), "celltype_response", "A")
    _place(fig.add_subplot(gs[0, 1]), "abundance_bar",     "B")
    _place(fig.add_subplot(gs[0, 2]), "signature_radar",   "C")
    _place(fig.add_subplot(gs[1, 0]), "signature_violins", "D")
    _place(fig.add_subplot(gs[1, 1]), "tissue_heatmap",    "E")
    _place(fig.add_subplot(gs[1, 2]), "volcano",           "F")
    fig.suptitle(
        "Figure 3  |  Immune engagement axis  —  Responders show "
        "Tfh-like CD4, dysfunctional-effector CD8, B/Plasma enrichment, "
        "and an exhaustion-cytotoxicity paradox",
        fontsize=13, fontweight="bold", y=1.01)
    _save(fig, "Figure3_immune_engagement")


# ====================================================================
# Figure 4 — Spatial axis
# ====================================================================
def figure4():
    fig = plt.figure(figsize=MAIN_FIGSIZE, constrained_layout=True)
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1.15, 0.85])
    _place(fig.add_subplot(gs[0, 0]), "spatial_sig_R",       "A")
    _place(fig.add_subplot(gs[0, 1]), "spatial_sig_NR",      "B")
    _place(fig.add_subplot(gs[0, 2]), "spatial_morans",      "C")
    _place(fig.add_subplot(gs[1, 0]), "spatial_rvnr_violin", "D")
    _place(fig.add_subplot(gs[1, 1]), "spatial_nhood",       "E")
    _place(fig.add_subplot(gs[1, 2]), "tls",                 "F")
    fig.suptitle(
        "Figure 4  |  Spatial axis  —  Visium (GSE238264; 4 R, 3 NR): "
        "Responder tumours show spatially organised TLS-like niches "
        "and higher Moran's I for immune signatures",
        fontsize=13, fontweight="bold", y=1.01)
    _save(fig, "Figure4_spatial_axis")


# ====================================================================
# Figure 5 — Cross-modality NMF integrator (L5)
# ====================================================================
def figure5():
    fig = plt.figure(figsize=MAIN_FIGSIZE, constrained_layout=True)
    gs = GridSpec(2, 3, figure=fig,
                  width_ratios=[0.85, 1.15, 1.0],
                  height_ratios=[1.0, 0.9])
    _place(fig.add_subplot(gs[:, 0]),  "nmf_topgenes",    "A")
    _place(fig.add_subplot(gs[0, 1]),  "nmf_forest",      "B")
    _place(fig.add_subplot(gs[0, 2]),  "nmf_bulk_vs_sc",  "C")
    _place(fig.add_subplot(gs[1, 1:]), "nmf_spatial",     "D")
    fig.suptitle(
        "Figure 5  |  Cross-modality NMF integrator  —  "
        "NMF9 (cytotoxic T: CCL5/NKG7/IL32) sig. higher in Responders; "
        "NMF3 (hepatocyte: ALB/APOA2) trends higher in Non-Responders; "
        "patterns recur across bulk and Visium",
        fontsize=13, fontweight="bold", y=1.01)
    _save(fig, "Figure5_nmf_integrator")


# ====================================================================
# Figure 6 — Multi-layer predictor & clinical translation
# ====================================================================
def figure6():
    fig = plt.figure(figsize=MAIN_FIGSIZE, constrained_layout=True)
    gs = GridSpec(2, 3, figure=fig,
                  height_ratios=[1.0, 1.05],
                  width_ratios=[1.0, 1.0, 1.1])
    _place(fig.add_subplot(gs[0, 0]), "layer_auc",          "A")
    _place(fig.add_subplot(gs[0, 1]), "layer_roc",          "B")
    _place(fig.add_subplot(gs[0, 2]), "variance_partition", "C")
    _place(fig.add_subplot(gs[1, 0]), "feature_effects",    "D")
    _place(fig.add_subplot(gs[1, 1]), "bulk_clustermap",    "E")
    _place(fig.add_subplot(gs[1, 2]), "survival",           "F")
    fig.suptitle(
        "Figure 6  |  Multi-layer predictor & clinical translation  —  "
        "balanced top-3 features per layer; composition dominates "
        "(AUC 0.94), CNV and NMF add orthogonal signal; "
        "validated on IMbrave150",
        fontsize=13, fontweight="bold", y=1.01)
    _save(fig, "Figure6_multilayer_predictor")


# ====================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("Compiling Manuscript Main Figures (publication-ready)")
    print("=" * 70)
    figure1()
    figure2()
    figure3()
    figure4()
    figure5()
    figure6()
    _normalize_canvas()
    print("=" * 70)
    print(f"Done. All main figures in: {OUT}")
    print("=" * 70)
