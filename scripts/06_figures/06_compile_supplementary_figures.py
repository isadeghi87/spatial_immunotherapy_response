#!/usr/bin/env python3
"""
06_compile_supplementary_figures.py — Supplementary Fig. S1-S8 composites
from individual QC/diagnostic/auxiliary outputs.

  Fig. S1 — Atlas QC & integration diagnostics
  Fig. S2 — Extended per-cohort & immune-subtype UMAPs
  Fig. S3 — Marker-based cell-type annotation support
  Fig. S4 — Spatial transcriptomics — downstream diagnostics
  Fig. S5 — Cross-cohort Cappuyns spatial consistency
  Fig. S6 — Signature correlation & bootstrap stability
  Fig. S7 — Cross-dataset reproducibility
  Fig. S8 — Patient-level immune subtyping

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
    # 01 — atlas QC
    "qc_violins":      PROJECT / "out" / "01_preprocessing" / "figures" / "01_qc_violins.png",
    "leiden":          PROJECT / "out" / "01_preprocessing" / "figures" / "01_leiden_clusters.png",
    "integration":     PROJECT / "out" / "01_preprocessing" / "figures" / "01_integration_umap.png",
    "celltype":        PROJECT / "out" / "01_preprocessing" / "figures" / "01_celltype_annotation_umap.png",
    "markerfeat":      PROJECT / "out" / "01_preprocessing" / "figures" / "01_marker_feature_plots.png",
    "immune_umaps":    PROJECT / "out" / "01_preprocessing" / "figures" / "01_immune_umaps.png",
    # 02a — marker evidence
    "dotplot":         PROJECT / "out" / "02a_response_tme" / "figures" / "02a_marker_dotplot.png",
    "hexbin":          PROJECT / "out" / "02a_response_tme" / "figures" / "02a_marker_feature_hexbin.png",
    "ct_enhanced":     PROJECT / "out" / "02a_response_tme" / "figures" / "02a_celltype_annotation_enhanced.png",
    "composition_bar": PROJECT / "out" / "02a_response_tme" / "figures" / "02a_composition_stacked_bar.png",
    "actexh_scatter":  PROJECT / "out" / "02a_response_tme" / "figures" / "02a_activation_exhaustion_scatter.png",
    "eda_degs":        PROJECT / "out" / "02a_response_tme" / "figures" / "eda_deg_counts_barplot.png",
    "top_degs_hmap":   PROJECT / "out" / "02a_response_tme" / "figures" / "eda_top_degs_heatmap.png",
    "eda_sig":         PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "eda_signature_mean_heatmap.png",
    # 02b — spatial downstream
    "pseudobulk_emb":  PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_pseudobulk_embeddings.png",
    "pseudobulk_io":   PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_pseudobulk_io_signatures.png",
    "ds_summary":      PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_ds_comprehensive_summary.png",
    "spatial_entropy": PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_ds_spatial_entropy.png",
    "spatial_grad":    PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_ds_gradient_spatial_maps.png",
    "spatial_topo":    PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02e_spatial_topology_rgb.png",
    "coocc":           PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_ds_celltype_cooccurrence.png",
    "coocc_diff":      PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_ds_cooccurrence_RvsNR_diff.png",
    "ct_ct_umap_resp": PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_ct_umap_celltype_response.png",
    "ct_ct_umap_imm":  PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_ct_umap_immune_subtypes.png",
    "ct_ct_umap_sam":  PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_ct_umap_sample.png",
    "ct_spmaps":       PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_ct_spatial_celltype_maps.png",
    "io_landscape":    PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "io_biomarker_landscape.png",
    "svg":             PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_ds_spatially_variable_genes.png",
    "loo":             PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_ds_leave_one_out.png",
    "sample_perm":     PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_ds_sample_level_permutation.png",
    # 02b — second Visium cohort
    "v_primrec":       PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_visium_gse281759_primary_vs_recurrence.png",
    "v_tumnorm":       PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_visium_gse281759_tumor_vs_normal.png",
    "v_phmap":         PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_visium_gse281759_patient_heatmap.png",
    "v_spmaps":        PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_visium_gse281759_spatial_maps.png",
    "v_imm":           PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_visium_immunotype_composition.png",
    "v_crossds":       PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_visium_cross_dataset_summary.png",
    "v_samp_heatmap":  PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_visium_sample_signature_heatmap.png",
    "v_sigcorr":       PROJECT / "out" / "02b_spatial_biomarkers" / "figures" / "02b_visium_signature_correlation.png",
    # Cappuyns auxiliaries
    "cap_imm_status":  PROJECT / "out" / "02i_cappuyns_signatures" / "figures" / "cappuyns_immune_status_comparison.png",
    # Validation support
    "sig_corr":        PROJECT / "out" / "03_validation" / "figures" / "03_signature_correlation.png",
    "bootstrap":       PROJECT / "out" / "03_validation" / "figures" / "03_bootstrap_ci.png",
    "bulk_cm":         PROJECT / "out" / "03_validation" / "figures" / "03_bulk_clustermap.png",
    "pathways":        PROJECT / "out" / "03_validation" / "figures" / "03_pathway_enrichment.png",
    # Complementary
    "hallmark":        PROJECT / "out" / "04_complementary" / "figures" / "04_hallmark_pathways.png",
    "crossds":         PROJECT / "out" / "04_complementary" / "figures" / "04_cross_dataset_signatures.png",
    "patient_sub":     PROJECT / "out" / "04_complementary" / "figures" / "04_patient_subtypes.png",
    "cd8_states":      PROJECT / "out" / "04_complementary" / "figures" / "04_cd8_states.png",
    "sum_multi":       PROJECT / "out" / "04_complementary" / "figures" / "04_summary_multipanel.png",
    "lr_heatmap":      PROJECT / "out" / "04_complementary" / "figures" / "04_ligand_receptor.png",
    "lr_dot_classic":  PROJECT / "out" / "04_complementary" / "figures" / "04_lr_dotplot.png",
    "imm_ratios":      PROJECT / "out" / "04_complementary" / "figures" / "04_immune_ratios.png",
}


def _place(ax, key, label=None):
    path = P.get(key)
    if path and path.exists():
        img = imread(str(path)); ax.imshow(img)
    else:
        ax.text(0.5, 0.5, f"[{key}]\nnot found", ha="center", va="center",
                fontsize=8, color="grey", transform=ax.transAxes)
    ax.set_axis_off()
    if label:
        ax.text(-0.02, 1.02, label, transform=ax.transAxes,
                fontsize=13, fontweight="bold", va="bottom", ha="right")


def _save(fig, stem):
    for fmt in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{fmt}", dpi=300,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {stem}  saved.")


# ── Fig. S1 — Atlas QC & integration diagnostics ───────────────────────
def s1():
    fig = plt.figure(figsize=(18, 14), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig)
    _place(fig.add_subplot(gs[0, 0]), "qc_violins",  "A")
    _place(fig.add_subplot(gs[0, 1]), "leiden",      "B")
    _place(fig.add_subplot(gs[1, 0]), "integration", "C")
    _place(fig.add_subplot(gs[1, 1]), "celltype",    "D")
    fig.suptitle("Supplementary Figure S1 — Atlas QC & integration diagnostics",
                 fontsize=13, fontweight="bold", y=1.01)
    _save(fig, "SupplementaryFig_S1_atlas_QC")


# ── Fig. S2 — Immune subsets UMAPs ─────────────────────────────────────
def s2():
    fig = plt.figure(figsize=(16, 12), constrained_layout=True)
    gs = GridSpec(1, 1, figure=fig)
    _place(fig.add_subplot(gs[0, 0]), "immune_umaps", None)
    fig.suptitle("Supplementary Figure S2 — Per-cohort and immune-subtype UMAPs",
                 fontsize=13, fontweight="bold", y=1.01)
    _save(fig, "SupplementaryFig_S2_immune_UMAPs")


# ── Fig. S3 — Marker-based cell-type annotation support ────────────────
def s3():
    fig = plt.figure(figsize=(20, 14), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig)
    _place(fig.add_subplot(gs[0, 0]), "dotplot",      "A")
    _place(fig.add_subplot(gs[0, 1]), "ct_enhanced",  "B")
    _place(fig.add_subplot(gs[1, 0]), "markerfeat",   "C")
    _place(fig.add_subplot(gs[1, 1]), "hexbin",       "D")
    fig.suptitle("Supplementary Figure S3 — Marker-based cell-type annotation support",
                 fontsize=13, fontweight="bold", y=1.01)
    _save(fig, "SupplementaryFig_S3_marker_annotation")


# ── Fig. S4 — Spatial transcriptomics diagnostics ─────────────────────
def s4():
    fig = plt.figure(figsize=(22, 14), constrained_layout=True)
    gs = GridSpec(2, 3, figure=fig)
    _place(fig.add_subplot(gs[0, 0]), "pseudobulk_emb",  "A")
    _place(fig.add_subplot(gs[0, 1]), "pseudobulk_io",   "B")
    _place(fig.add_subplot(gs[0, 2]), "spatial_topo",    "C")
    _place(fig.add_subplot(gs[1, 0]), "spatial_entropy", "D")
    _place(fig.add_subplot(gs[1, 1]), "coocc",           "E")
    _place(fig.add_subplot(gs[1, 2]), "coocc_diff",      "F")
    fig.suptitle("Supplementary Figure S4 — Spatial transcriptomics — downstream diagnostics",
                 fontsize=13, fontweight="bold", y=1.01)
    _save(fig, "SupplementaryFig_S4_spatial_diagnostics")


# ── Fig. S5 — Cross-cohort Visium consistency ─────────────────────────
def s5():
    fig = plt.figure(figsize=(20, 14), constrained_layout=True)
    gs = GridSpec(2, 3, figure=fig)
    _place(fig.add_subplot(gs[0, 0]), "v_primrec",  "A")
    _place(fig.add_subplot(gs[0, 1]), "v_tumnorm",  "B")
    _place(fig.add_subplot(gs[0, 2]), "v_phmap",    "C")
    _place(fig.add_subplot(gs[1, 0]), "v_spmaps",   "D")
    _place(fig.add_subplot(gs[1, 1]), "v_imm",      "E")
    _place(fig.add_subplot(gs[1, 2]), "v_crossds",  "F")
    fig.suptitle("Supplementary Figure S5 — Cross-cohort Visium consistency (GSE238264 & GSE281759)",
                 fontsize=13, fontweight="bold", y=1.01)
    _save(fig, "SupplementaryFig_S5_crosscohort_visium")


# ── Fig. S6 — Signature correlation & bootstrap stability ──────────────
def s6():
    fig = plt.figure(figsize=(20, 14), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig)
    _place(fig.add_subplot(gs[0, 0]), "sig_corr",   "A")
    _place(fig.add_subplot(gs[0, 1]), "bootstrap",  "B")
    _place(fig.add_subplot(gs[1, 0]), "bulk_cm",    "C")
    _place(fig.add_subplot(gs[1, 1]), "pathways",   "D")
    fig.suptitle("Supplementary Figure S6 — Signature correlation, bootstrap stability, and pathway enrichment",
                 fontsize=13, fontweight="bold", y=1.01)
    _save(fig, "SupplementaryFig_S6_signature_stability")


# ── Fig. S7 — Cross-dataset reproducibility ───────────────────────────
def s7():
    fig = plt.figure(figsize=(20, 14), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig)
    _place(fig.add_subplot(gs[0, 0]), "hallmark",   "A")
    _place(fig.add_subplot(gs[0, 1]), "crossds",    "B")
    _place(fig.add_subplot(gs[1, 0]), "sum_multi",  "C")
    _place(fig.add_subplot(gs[1, 1]), "lr_heatmap", "D")
    fig.suptitle("Supplementary Figure S7 — Cross-dataset reproducibility & pathway summary",
                 fontsize=13, fontweight="bold", y=1.01)
    _save(fig, "SupplementaryFig_S7_reproducibility")


# ── Fig. S8 — Patient-level immune subtyping ──────────────────────────
def s8():
    fig = plt.figure(figsize=(20, 14), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig)
    _place(fig.add_subplot(gs[0, 0]), "patient_sub",    "A")
    _place(fig.add_subplot(gs[0, 1]), "cd8_states",     "B")
    _place(fig.add_subplot(gs[1, 0]), "imm_ratios",     "C")
    _place(fig.add_subplot(gs[1, 1]), "cap_imm_status", "D")
    fig.suptitle("Supplementary Figure S8 — Patient-level immune subtyping & ratio biomarkers",
                 fontsize=13, fontweight="bold", y=1.01)
    _save(fig, "SupplementaryFig_S8_patient_subtyping")


if __name__ == "__main__":
    print("=" * 70)
    print("Compiling Supplementary Figures S1–S8")
    print("=" * 70)
    s1(); s2(); s3(); s4(); s5(); s6(); s7(); s8()
    print("=" * 70)
    print(f"Done. Supplementary figures in: {OUT}")
    print("=" * 70)
