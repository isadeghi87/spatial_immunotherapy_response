#!/usr/bin/env python3
"""Generate biological figures_interpretation.md (+ .docx via pandoc) for every
subfolder under `out/` that contains PNG figures.

Descriptions are grounded in the manuscript findings (docs/manuscript.md, May 2026
revision) and the corrected results recorded in /memories/repo/spatial_io_project.md.

Usage:
    python scripts/06_figures/generate_figure_interpretations.py
"""
from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path
from textwrap import dedent

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "out"

# ---------------------------------------------------------------------------
# Per-folder thematic introduction
# ---------------------------------------------------------------------------
FOLDER_INTROS: dict[str, str] = {
    "00_pipeline_audit": (
        "Pipeline-audit logs and YAML summarising the reframe-stage rerun of the "
        "HCC immunotherapy atlas. These artefacts document which scripts ran, "
        "which outputs were regenerated, and which QC checks passed (no figures)."
    ),
    "01_preprocessing": (
        "Stage-1 figures cover atlas QC, per-dataset integration, scVI/Harmony "
        "benchmarking, and reference cell-type annotation across the 1,058,964-cell "
        "integrated HCC scRNA-seq atlas (5 HCC cohorts + non-tumour liver reference)."
    ),
    "02_analysis": (
        "Composite analysis outputs (tables + figures) for the harmonised atlas "
        "before stage-specific stratification."
    ),
    "02a_response_tme": (
        "Tumour-microenvironment composition and signature contrasts between "
        "Responders (R, n=7) and Non-Responders (NR, n=17) in the GSE206325 "
        "anti-PD-1 discovery cohort. Patient-level Mann–Whitney U with BH-FDR; "
        "effect sizes reported as Cliff's δ (positive = higher in R)."
    ),
    "02b_spatial_biomarkers": (
        "Visium spatial transcriptomics analyses across GSE173905, GSE238264 "
        "(4R/3NR) and GSE281759, combined with a 24-section scVI-integrated "
        "spatial atlas (43,277 spots × 30 latent dimensions) and per-sample "
        "InferCNV-based tumour delineation."
    ),
    "02c_bulk_deconvolution": (
        "Pseudobulk and CIBERSORTx-style deconvolution of the integrated atlas "
        "(79 patients × 15,697 genes), with cross-trial projection onto IMbrave150."
    ),
    "02i_cappuyns_signatures": (
        "Signature-scoring outputs on the Cappuyns et al. 2023 atezolizumab + "
        "bevacizumab external validation cohort (32 response-evaluable patients, "
        "19 R / 13 NR; pre-treatment week-0 biopsies)."
    ),
    "02j_concordance": (
        "Internal concordance checks between Cappuyns-derived signatures and the "
        "discovery cohort's molecular sub-typing."
    ),
    "02k_cappuyns_spatial": (
        "Section-level signature heatmaps from the Cappuyns spatial sub-cohort, "
        "used to test cross-platform program coherence."
    ),
    "03_response_modeling": (
        "Logistic-regression / LASSO hyper-parameter sensitivity for the multi-layer "
        "predictor."
    ),
    "03_validation": (
        "Cross-cohort validation and biomarker robustness: LASSO/RF ROC, bootstrap "
        "CIs, pathway enrichment, bulk RNA-seq concordance, public-data IO "
        "biomarker landscape, and spatial spot QC."
    ),
    "04_complementary": (
        "Complementary biological axes: CD8 T-cell sub-states and pseudotime, "
        "TLS spatial hubs, ligand–receptor cross-talk, TCR clonotype dynamics, "
        "metabolic reprogramming, hallmark pathway shifts, and IMbrave150 real "
        "per-patient survival validation."
    ),
    "05_multilayer": (
        "Five-layer predictive framework: L1 tumour-intrinsic CNV, L2 regulons, "
        "L3 cell-type composition, L4 IO signatures, L5 NMF programs. Includes "
        "balanced top-3-features-per-layer LOOCV (24-patient discovery), 1,000-"
        "iteration permutation null, 2,000-iteration bootstrap CIs, SHAP "
        "attribution, NMF program forest plots, CNV chromosome/pathway analyses, "
        "Geneformer benchmark comparisons, and held-out Cappuyns transfer."
    ),
    "06_manuscript_figures": (
        "Final, manuscript-grade composite figures (main Figures 1–6, Extended "
        "Data Fig. 1–5, Supplementary Fig. S1–S9) assembled with strict "
        "publication-style typography, 6650×4200 master canvas, and shared "
        "colour palette (Responder #0077BB / Non_Responder #CC3311)."
    ),
    "06_reframe": (
        "Reframe-stage notes and tables that justify the May 2026 corrected "
        "narrative (cnv_clone_diversity null in CNV atlas; NMF10 IgG/antigen-"
        "presentation programme as the BH-significant responder programme; "
        "discovery-leading composition layer fails Cappuyns transfer)."
    ),
}

# ---------------------------------------------------------------------------
# Figure-specific biological interpretation registry
# Keys are matched as substrings against the PNG basename (without extension).
# The first match wins, so most specific keys should come first.
# ---------------------------------------------------------------------------
FIGURE_DESCRIPTIONS: list[tuple[str, str]] = [
    # ── Main manuscript composites ────────────────────────────────────────
    ("Figure1_atlas_landscape", "Integrated atlas landscape. (a) Harmony-corrected UMAP of 1,058,964 cells coloured by manual cell-type annotation, with 44 sub-lineages spanning malignant hepatocytes, six CD8 states, CD4 helpers/Tregs, NK, mregDC, conventional DC, M1/M2 macrophages, B/plasma cells, endothelium and fibroblasts. (b) Same embedding coloured by dataset of origin — overlap confirms biology, not batch, dominates the manifold (dataset silhouette 0.41→0.04 post-Harmony; cell-type silhouette preserved 0.27→0.26). (c) Marker dot plot for lineage-canonical genes (CD3D/CD8A/CD4/FOXP3/CD79A/KLRD1/LYZ/COL1A1/ALB). (d) Per-patient stacked composition for the 24 GSE206325 anti-PD-1-treated patients grouped by response label."),
    ("Figure2_tumour_intrinsic", "Tumour-intrinsic axis. (a) Per-cell InferCNV score on UMAP — hepatocyte clusters carry the expected elevated burden. (b) Per-patient chromosome heatmap with recurrent 1q/8q gains and 4q/8p/17p losses. (c) InferCNV dashboard summarising per-patient burden, malignant fraction and clone diversity. (d) Patient-level R-vs-NR violins for `cnv_score_mean` (p=0.38, n.s.), `cnv_score_sd` (p=0.065, borderline), malignant fraction (p=0.90, n.s.) and `cnv_clone_diversity` (Cliff's δ=−1.03, p=0.026) — the dominant tumour-intrinsic non-response axis is *clonal fragmentation*, not mean burden. (e) Higher-clonality patients show depleted lymphoid and enriched myeloid signatures, linking genomic instability to immune-excluded TME structure."),
    ("Figure3_immune_engagement", "Immune engagement axis. (a) Immune-only UMAP coloured by response. (b) Per-cell-type abundance bars (mean ± SE of per-patient fraction) with BH-adjusted Mann–Whitney p-values; Dysfunctional-effector CD8 (δ=+0.78, p=0.010), Cytotoxic and Effector CD8, and plasma cells expand in R, while Naïve/Memory CD8 contract. (c) Radar plot of ten IO signatures per response group. (d) Violins for cytotoxicity, exhaustion, IFN-γ, Tfh-like and TLS signatures. (e) Tissue-level z-scored heatmap of twenty hallmark + IO signatures. (f) Pseudo-bulk DE volcano (DESeq2 on patient pseudobulks)."),
    ("Figure4_spatial_axis", "Spatial projection of responder programmes onto Visium. (a) Responder GSE238264 sections show spatially clustered immune/TLS hubs at tumour–stroma interfaces; NMF10 (mean section activity 0.274 vs 0.058 in NR) and cytotoxicity co-localise with TLS-defining IGHG1/CD79A spots. (b) Non-Responder sections are dominated by NMF3 hepatocyte-differentiation activity (0.944 vs 0.740 in R) and diffuse, low-autocorrelation cytotoxicity. (c) Per-section Moran's I — responder signatures show higher spatial structure (0.12–0.41 vs 0.03–0.09; p=0.042). (d) Per-spot R-vs-NR signature distributions with permutation p-values. (e) Neighbourhood-enrichment heatmap between deconvolved cell types. (f) Top TLS-score tiles: median TLS-hub count 18 (R) vs 6 (NR); p=0.057 at n=7, hypothesis-generating."),
    ("Figure5_nmf_integrator", "NMF cross-modality integrator. (a) Top-8 gene loadings per NMF program (k=10). (b) Forest plot of Cliff's δ for per-patient program activity (R − NR), annotated with top genes: NMF10 (IGKC/IGHG1/SSR4/CD74/MZB1; δ=+0.85, BH-p=0.0054) is the BH-significant responder programme; NMF8 (CCL4/NKG7/CD8A/CTSW/PRF1; nominal R) and NMF3 (APOA2/ALB/APOC3/APOC1; nominal NR) are secondary axes. (c) Cross-modality concordance — sc-derived NMF programmes projected onto bulk RNA-seq. (d) Per-section NMF activity heatmap (HCC1R–HCC7NR); NMF10 enriched in R, NMF3 in NR sections."),
    ("Figure6_multilayer_predictor", "Balanced five-layer comparison and external transfer. (a) Per-layer LOOCV ROC-AUC / PR-AUC / balanced-accuracy bars with 2,000-iteration bootstrap 95% CIs and the 1,000-iteration label-shuffle permutation 95% null. L3=0.96 (perm p=0.014); L3-noTech=0.95 (perm p=0.022); L5=0.95; L4=0.90; L1=0.78; L2=0.55; ensemble=0.96 (perm p=0.015). (b) ROC curves per layer. (c) Variance-partition: removing L5 produces the largest ensemble AUC drop. (d) Coefficient plot. (e) IMbrave150 per-patient KM by intratumoural CD8 IHC, within-arm: atezo+bev PFS log-rank p=0.009; sorafenib p=0.90. Cox arm × CD8 interaction HR=0.57 PFS / 0.36 OS — predictive, not purely prognostic. (f) Bulk-RNA-seq clustermap."),
    ("Fig6A_balanced_layer_metrics", "Updated Figure 6 panel A — balanced LOOCV metrics with ROC-AUC, PR-AUC and balanced-accuracy bars per layer plus the 15-feature ensemble. Permutation-null 95-th percentiles are marked in red for L3 and the ensemble (the only entries exceeding their null at p<0.05). PR-AUC prevalence baseline = 0.29."),
    # ── Extended Data ────────────────────────────────────────────────────
    ("ExtendedData_Fig1_CNV_genomic_instability", "Per-patient InferCNV summaries and clonal-heterogeneity distributions — extends the patient-level CNV story of Figure 2 with per-clone violins and arm-level Responder − Non-Responder residuals."),
    ("ExtendedData_Fig2_regulon_communication", "L2 SCENIC-style AUCell activity for 12 master TFs (IRF1/STAT1/TCF7/TOX/EOMES/FOXP3/MYC/SPI1/IRF8/NR4A1/BATF/PRDM1) and LIANA ligand–receptor pair activity — confirms L2 is informative biologically (CXCL9/CXCR3, CCL5/CCR5, HLA-A/CD8A leading) even though it sits near chance as a stand-alone classifier."),
    ("ExtendedData_Fig3_high_contrast_spatial", "High-contrast spatial maps for hypoxia, CD8 exhaustion, M1−M2 polarity and IO response on both the UMAP and Visium coordinates, emphasising spatial gradients that survive the publication-style colour normalisation."),
    ("ExtendedData_Fig4_bulk_trial_validation", "IMbrave150 trial-level KM stratified by intratumoural CD8 IHC (median split), within each arm. Within-arm atezo+bev PFS p=0.009; sorafenib p=0.90 — predictive, not purely prognostic, and based on real per-patient MOESM4 data (n=177)."),
    ("ExtendedData_Fig4_imbrave150_real", "IMbrave150 real per-patient KM curves with Cox arm × CD8-IHC interaction (HR=0.57 PFS, HR=0.36 OS) — the figure used in main-text §3.10 to anchor the dominant L3 composition signal in registration-trial data."),
    ("ExtendedData_Fig5_cappuyns_validation", "Held-out external validation on Cappuyns 2023 (32 atezo+bev patients, 19 R / 13 NR). Discovery-vs-validation ROC overlay shows the L3 composition model collapses to chance (validation ROC-AUC=0.50, CI 0.50–0.50), whereas L5 NMF (0.69, 0.49–0.86) and L2 regulons (0.68) retain modest transfer. Per-patient waterfall illustrates the failure of composition to discriminate at the validation cohort."),
    # ── Supplementary ────────────────────────────────────────────────────
    ("SupplementaryFig_S1_atlas_QC", "Atlas QC and integration diagnostics: per-cell n_genes / n_UMI / pct_mt violins; dataset and cell-type silhouette pre/post Harmony; cluster-stability bars."),
    ("SupplementaryFig_S2_immune_UMAPs", "Per-lineage immune UMAPs (CD8, CD4, NK, Treg, B/plasma, macrophage/monocyte, DC) coloured by donor, dataset and response — supports the visual separation of dysfunctional-effector CD8 and Tfh-like CD4 niches."),
    ("SupplementaryFig_S3_marker_annotation", "Curated lineage-marker enrichment heatmap and dot plot — the basis for the 44-sub-lineage `cell_type_original` column."),
    ("SupplementaryFig_S4_spatial_diagnostics", "Per-section Visium QC: spot counts, gene-detection, mitochondrial fraction, and within-section Moran's I null distributions."),
    ("SupplementaryFig_S5_crosscohort_visium", "Cross-cohort spatial reproducibility across GSE173905 / GSE238264 / GSE281759 — same IO signature ranks recovered in independent platforms."),
    ("SupplementaryFig_S6_signature_stability", "Signature-stability and NMF-rank selection diagnostics: reconstruction error vs k ∈ {5,7,10,12,15,20}; cophenetic / silhouette curves identifying k=10."),
    ("SupplementaryFig_S7_reproducibility", "Leave-one-cohort-out reproducibility of the L3 / L5 layer metrics — none of GSE125449, GSE140228, GSE149614 dominates the discovery signal."),
    ("SupplementaryFig_S8_patient_subtyping", "Hierarchical clustering of the 24 GSE206325 patients by 16 IO signatures into S1 immune-active (6R/3NR), S2 intermediate (0R/2NR), S3 immune-cold (1R/12NR); χ² p=0.0092."),
    ("SupplementaryFig_S9_patient_demographics", "Patient and cohort demographics for GSE206325 + Cappuyns 2023 — age, sex, etiology, BCLC stage, viral status, prior treatment, and per-patient cell counts."),
    # ── 01_preprocessing ────────────────────────────────────────────────
    ("01_qc_violins", "Per-cell QC violins (n_genes_by_counts, total_counts, pct_counts_mt) before and after Scanpy QC thresholds (200 ≤ n_genes ≤ 8000, pct_mt < 15%)."),
    ("01_integration_umap", "Pre- vs post-Harmony UMAP coloured by dataset and donor — visual confirmation that Harmony removes batch while preserving cell-type structure."),
    ("01_celltype_umap", "UMAP coloured by the 44-class `cell_type_original` annotation used throughout the manuscript."),
    ("01_immune_umaps", "Immune-restricted UMAPs split by lineage."),
    ("01_leiden_clusters", "Leiden clustering at resolution 0.8 on the harmonised 15-NN graph."),
    ("01_marker_feature_plots", "Marker-gene feature plots overlaying canonical genes on the harmonised UMAP."),
    ("01_scanvi_confidence", "scANVI label-transfer confidence per cell — supports the manual annotation against an automated baseline."),
    ("01c_batch_metrics_summary", "Bar summary of ASW-batch / ASW-bio / iLISI / cLISI / kBET for PCA vs Harmony (matches §2.16 of M&M)."),
    ("01c_umap_integration_comparison", "Side-by-side UMAPs from PCA, Harmony, and scVI embeddings."),
    ("01c_lisi_histograms", "Local Inverse-Simpson histograms: iLISI (dataset) and cLISI (cell type) before and after integration."),
    ("01d_scvi_umap", "scVI latent UMAP used as an orthogonal integration cross-check."),
    ("01d_scvi_metric_bars", "scVI metric summary bars vs Harmony baseline."),
    ("01y_umap_dataset", "UMAP coloured by dataset across all cohorts in the integrated atlas."),
    ("01y_umap_celltype", "UMAP coloured by cell-type after Harmony."),
    ("01y_umap_response", "UMAP coloured by RECIST response across the 24-patient discovery subset."),
    ("01y_visium_celltype", "Visium spot-level cell-type calls from Cell2location deconvolution."),
    ("01y_visium_immunotype", "Per-spot immunotype assignment derived from joint IO-signature thresholds."),
    ("01y_proportion_by_dataset", "Per-dataset cell-type proportion bars — sanity check for relative lineage capture across cohorts."),
    ("01y_proportion_cappuyns", "Cell-type proportion bars in the Cappuyns 2023 external validation cohort."),
    ("01y_proportion_gse206325", "Cell-type proportion bars in the GSE206325 anti-PD-1 discovery cohort."),
    # ── 02a_response_tme ────────────────────────────────────────────────
    ("02a_abundance_barplot", "Per-cell-type abundance bars for R vs NR with BH-adjusted Mann–Whitney p-values."),
    ("02_abundance_barplot", "Legacy duplicate of the per-cell-type abundance bars (kept for traceability)."),
    ("02a_celltype_annotation_enhanced", "Enhanced annotation panel combining marker dot plot, UMAP and per-patient prevalence."),
    ("02a_celltype_response_umap", "Cell-type UMAP coloured by R/NR; visually highlights expansion of Dysfunctional-effector CD8 in R."),
    ("02_celltype_response_umap", "Legacy duplicate of the response-coloured cell-type UMAP."),
    ("02a_composition_stacked_bar", "Per-patient stacked composition bars ordered by response group."),
    ("02_composition_stacked_bar", "Legacy duplicate of the stacked composition bar plot."),
    ("02a_marker_dotplot", "Lineage-marker dot plot used for the annotation pass on GSE206325."),
    ("02a_marker_feature_hexbin", "Hexbin marker feature plots — handles high cell density in the immune compartment."),
    ("02a_signature_radar", "Radar plot of 10 IO signatures for R vs NR."),
    ("02_signature_radar", "Legacy duplicate of the IO signature radar plot."),
    ("02a_signature_violins", "Per-patient violins of cytotoxicity, exhaustion, IFN-γ, Tfh-like and TLS signatures."),
    ("02_signature_violins", "Legacy duplicate of the per-patient signature violins."),
    ("02a_tissue_signature_heatmap", "Tissue-level z-scored heatmap of 20 hallmark + IO signatures across R/NR samples."),
    ("02_tissue_signature_heatmap", "Legacy duplicate of the tissue-level signature heatmap."),
    ("02a_volcano_deg", "Pseudo-bulk DE volcano (DESeq2 on patient pseudobulks; |log2FC|≥1, padj<0.05)."),
    ("02_volcano_deg", "Legacy duplicate of the pseudo-bulk DE volcano."),
    ("02_activation_exhaustion_scatter", "Per-patient scatter of activation vs exhaustion signature means — Responders cluster in a 'high-activation-with-productive-dysfunction' region."),
    ("02f_tf_regulon_activity", "L2 SCENIC-style TF regulon activity heatmap per patient; supports Extended Data Fig. 2."),
    ("eda_deg_counts_barplot", "Exploratory DEG-count bars per cell type from the pseudo-bulk DE pass."),
    ("eda_top_degs_heatmap", "Top-DEG heatmap from pseudo-bulk DE — exploratory companion to the volcano plot."),
    # ── 02b_spatial_biomarkers ──────────────────────────────────────────
    ("02b_ct_marker_dotplot", "Spatial-atlas cell-type marker dot plot from the 24-section scVI-integrated Visium analysis."),
    ("02b_ct_marker_heatmap", "Spatial-atlas marker heatmap supporting the spot-level annotation."),
    ("02b_ct_sample_annotation_support", "Per-sample support bars for each cell type — confirms cross-section coverage of every annotated class."),
    ("02b_ct_spatial_celltype_maps", "Per-section spot maps coloured by dominant cell-type call."),
    ("02b_ct_umap_celltype_response", "scVI-spatial UMAP of all 43,277 spots coloured by cell type and response."),
    ("02b_ct_umap_immune_subtypes", "Spatial-UMAP focus on immune sub-types only."),
    ("02b_ct_umap_leiden", "Leiden clusters on the scVI-spatial latent embedding (within-sample 20-NN purity 0.576; random ≈ 0.042)."),
    ("02b_ct_umap_sample", "Spatial-UMAP coloured by section ID — confirms cross-section integration."),
    ("02b_ds_celltype_cooccurrence", "Pairwise neighbour cell-type co-occurrence in Visium spots."),
    ("02b_ds_cluster_profiles", "Per-cluster cell-type composition profiles from the spatial atlas."),
    ("02b_ds_comprehensive_summary", "Composite spatial-discovery summary: cell-type composition + IO signatures + spatial clusters per section."),
    ("02b_ds_cooccurrence_RvsNR_diff", "Difference of pairwise neighbour co-occurrence (R − NR) across 24 sections — emphasises gain of Dysfunctional-effector CD8 ↔ Tfh-like CD4 contacts in R."),
    ("02b_ds_gradient_spatial_maps", "Per-section spatial gradient maps of cytotoxicity vs hepatocyte differentiation."),
    ("02b_ds_immune_gradient", "Tumour-edge → tumour-core immune gradient per section."),
    ("02b_ds_leave_one_out", "Leave-one-section-out reproducibility of the spatial layer."),
    ("02b_ds_ligand_receptor_heatmap", "LIANA ligand–receptor scoring heatmap for the Visium cohort."),
    ("02b_ds_LR_RvsNR_comparison", "R vs NR comparison of top ligand–receptor pairs (CXCL9/CXCR3, CCL5/CCR5, HLA-A/CD8A)."),
    ("02b_ds_sample_level_permutation", "Per-section permutation null for IO signature R/NR contrasts."),
    ("02b_ds_spatial_clusters", "Per-section spatial Leiden cluster maps."),
    ("02b_ds_spatial_entropy", "Per-section spatial Shannon entropy of cell-type composition."),
    ("02b_ds_spatially_variable_genes", "Top spatially variable genes (Squidpy Moran's I, 999-permutation null) per section."),
    ("02b_pseudobulk_embeddings", "Spatial pseudobulk UMAP/PCA of sections, coloured by response."),
    ("02b_pseudobulk_io_signatures", "Spatial pseudobulk IO signature heatmap across sections."),
    ("02b_visium_celltype_spatial_maps", "Per-section Cell2location cell-type spatial maps."),
    ("02b_visium_cross_dataset_summary", "Cross-dataset spatial summary across GSE173905 / GSE238264 / GSE281759."),
    ("02b_visium_gse281759_patient_heatmap", "GSE281759 multi-region per-patient signature heatmap (primary vs recurrence)."),
    ("02b_visium_gse281759_primary_vs_recurrence", "Primary vs recurrence R/NR contrast on GSE281759 sections."),
    ("02b_visium_gse281759_spatial_maps", "Per-section signature spatial maps on GSE281759."),
    ("02b_visium_gse281759_tumor_vs_normal", "Tumour vs adjacent-normal contrast on GSE281759 sections."),
    ("02b_visium_immunotype_composition", "Per-section immunotype composition bars."),
    ("02b_visium_morans_I", "Section-level Moran's I per IO signature (responder-associated signatures show 0.12–0.41 vs 0.03–0.09 in NR)."),
    ("02b_visium_neighbor_enrichment_heatmap", "Neighbourhood-enrichment heatmap between cell-type categories."),
    ("02b_visium_nhood_enrichment", "Squidpy nhood-enrichment z-scores."),
    ("02b_visium_RvsNR_violins", "Per-spot R vs NR signature violins aggregated across sections."),
    ("02b_visium_sample_signature_heatmap", "Per-sample mean signature heatmap (HCC1R–HCC7NR rows)."),
    ("02b_visium_signature_correlation", "Pairwise correlation of 16 IO signatures across spots."),
    ("02b_visium_spatial_immunotype_maps", "Per-section immunotype spatial maps."),
    ("02b_visium_spatial_signature_maps", "Per-section IO-signature spatial overlays."),
    ("02b_visium_spatial_signature_maps_responder", "Responder-only signature spatial maps — TLS / cytotoxicity hubs at tumour–stroma interfaces."),
    ("02b_visium_spatial_signature_maps_non_responder", "Non-responder-only signature spatial maps — NMF3 hepatocyte dominance."),
    ("eda_signature_mean_heatmap", "Exploratory mean-signature heatmap from the spatial atlas pass."),
    ("spatial_scvi_umap_composite", "Six-panel scVI spatial UMAP composite (sample, dataset, response, cell type, CNV-derived tumour status, immunotype)."),
    ("io_biomarker_landscape", "Cross-cohort IO biomarker landscape — projects atlas-derived signatures onto public IO trial transcriptomes."),
    # ── 02b qc panels ────────────────────────────────────────────────────
    ("_qc_panel", "Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail."),
    ("02b_qc_overview", "Cross-section Visium QC overview."),
    # ── 02b CNV ─────────────────────────────────────────────────────────
    ("spatial_cnv_universal_heatmap", "Universal infercnv heatmap across 23/24 sections (43,196 spots × 2,115 genome windows) — HCC5NR excluded (0 reference cells)."),
    ("spatial_cnv_", "Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline."),
    ("infercnv.preliminary", "Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only."),
    ("infercnv.png", "Final infercnv heatmap (post-HMM, denoised) for this section."),
    # ── 02c bulk ─────────────────────────────────────────────────────────
    ("02c_reference_cellcounts", "CIBERSORTx-style reference: cell counts per signature cell type used to build the 43-class reference matrix."),
    ("02c_reference_correlation", "Reference-matrix cell-type correlation heatmap — confirms low cross-cell-type collinearity."),
    ("02c_reference_heatmap", "Reference signature matrix heatmap (15,853 genes × 60 sub-types)."),
    ("bulk_deconv_advanced_composite", "Nine-panel bulk deconvolution composite: UMAP, stacked composition, box plots, volcano, heatmap, radar, co-abundance, deconv-vs-truth, KM (§P3 of the project memory)."),
    ("bulk_io_heatmap", "Bulk IO-signature heatmap projected from the atlas onto pseudobulk samples (79 patients × 15,697 genes)."),
    ("eda_highly_variable_genes_heatmap", "Exploratory HVG heatmap supporting the pseudobulk normalisation."),
    # ── 02i / 02j / 02k Cappuyns ────────────────────────────────────────
    ("cappuyns_immune_status_comparison", "Cappuyns 2023 immune-status comparison: per-patient IO signature means across R/NR/excluded subsets."),
    ("cappuyns_key_signatures_response", "Per-signature R-vs-NR boxplots in the Cappuyns cohort."),
    ("cappuyns_molecular_subsets", "Cappuyns molecular sub-set assignment via the discovery signature panel."),
    ("cappuyns_signature_heatmap", "Per-patient signature heatmap on the Cappuyns external cohort."),
    ("cappuyns_internal_concordance", "Internal concordance check: Cappuyns scoring against the discovery subtyping (Cohen's κ ≈ 0.55)."),
    ("cappuyns_spatial_heatmaps", "Cappuyns spatial sub-cohort signature heatmaps."),
    ("cappuyns_spatial_section_comparison", "Per-section comparison of immune programmes in the Cappuyns spatial subset."),
    # ── 03 modeling / validation ────────────────────────────────────────
    ("05m_lasso_C_sensitivity", "LASSO penalty (C) sensitivity scan for the response classifier — selects C=0.5 for the 15-feature ensemble."),
    ("03_signature_correlation", "Cross-cohort signature correlation — IFN-γ and cytotoxicity correlate strongly across datasets, supporting cross-cohort transfer."),
    ("03_roc_lasso", "ROC for the LASSO classifier on the discovery cohort (AUC≈0.97 internal)."),
    ("03_rf_importance", "Random-forest feature importance — top contributors are CD8 Dysfunctional-effector, Tfh-like, IFN-γ score and NMF10."),
    ("03_bootstrap_ci", "2,000-iteration bootstrap CI distributions for ROC-AUC by classifier."),
    ("03_bulk_clustermap", "Bulk RNA-seq sample clustermap scored with atlas-derived signatures."),
    ("03_pathway_enrichment", "Hallmark pathway enrichment of R-vs-NR DEGs (IFN-α, IFN-γ, TNF-NFκB elevated in R)."),
    ("03b_clinical_survival_validation", "Clinical-survival validation panel: KM stratification by composition-derived risk score."),
    ("03c_validation_perlayer", "Per-layer external-validation AUC bars on Cappuyns — L5 NMF (0.69) and L2 regulon (0.68) retain transfer; L3 collapses to 0.50."),
    ("03c_validation_roc", "Per-layer external-validation ROC curves on Cappuyns."),
    ("03c_validation_waterfall", "Per-patient L3 prediction-score waterfall on Cappuyns illustrating the failure of cross-regimen transfer."),
    ("03d_spatial_spot_qc", "Spatial-spot QC summary on the validation Visium subset."),
    ("io_biomarker_public", "Public-data IO biomarker landscape projection."),
    # ── 04 complementary ────────────────────────────────────────────────
    ("04_cd8_states", "CD8 T-cell sub-state panel (naïve-like, memory, effector, cytotoxic, dysfunctional-effector, dysfunctional-terminal/progenitor) with per-patient R/NR enrichment."),
    ("04_cross_dataset_signatures", "Cross-dataset IO-signature reproducibility across the five HCC cohorts."),
    ("04_hallmark_pathways", "Hallmark pathway R-vs-NR shift summary — IFN-α/IFN-γ/TNF-NFκB ↑R; WNT/β-catenin / oxphos / coagulation / TGF-β ↑NR."),
    ("04_immune_ratios", "Composite immune ratios: T/Tumor AUC=0.87, CD8/Treg AUC=0.82, DC/Mac=0.49, Effector/Exhausted=0.44."),
    ("04_ligand_receptor", "Ligand–receptor heatmap from LIANA across the discovery cohort."),
    ("04_lr_dotplot", "Top R-vs-NR ligand–receptor pairs as a dot plot (size = expression, colour = R−NR)."),
    ("04_patient_subtypes", "Internal patient subtyping (S1 immune-active 6R/3NR, S2 intermediate 0R/2NR, S3 immune-cold 1R/12NR; χ² p=0.0092)."),
    ("04_summary_multipanel", "Multi-panel complementary summary tying CD8 states, signatures, LR pairs and patient subtypes into one figure."),
    ("04b_bulk_trial_survival_validation", "Bulk trial-survival validation panel — same atezo+bev arm-specific PFS stratification as Figure 6e."),
    ("04b_cell_communication_network", "CellPhoneDB-style cell–cell communication network — mregDC orchestrates dysfunctional-effector CD8 interactions in R."),
    ("04c_tcell_pseudotime_trajectories", "CD8 diffusion-pseudotime trajectories per patient; root = naïve-like cluster, rescaled to [0,1]; patient-level R-vs-NR median 0.092 vs 0.091 (δ=+0.38, p=0.166, under-powered)."),
    ("04d_ccc_lr_dotplot", "Cell–cell communication LR dot plot (per cell-type pair)."),
    ("04e_tcr_clonotype_alluvial", "TCR clonotype alluvial plot tracking clonal expansion across CD8 states."),
    ("04f_metabolic_reprogramming", "Metabolic reprogramming axes (oxphos vs glycolysis vs FAO) per patient."),
    ("04g_tls_spatial_hubs", "TLS spatial hub detection (12-gene panel CXCL13/CCL19/CCL21/CXCR5/CR2/CD79A/MS4A1/LAMP3/CD83/CCR7/IGHG1/IGHG2); median hub count 18 (R) vs 6 (NR), p=0.057 at n=7."),
    ("04h_tissue_stratified_spatial", "Tissue-stratified spatial signature maps (tumour core / tumour edge / stromal)."),
    ("04h_spatial_high_contrast_maps", "High-contrast spatial maps for the manuscript companion (Extended Data Fig. 3)."),
    ("04i_trajectory_response", "Pseudotime → response link: per-patient mean pseudotime vs R/NR (under-powered direction-consistent shift)."),
    ("spatial_1_hypoxia_gradient", "Spatial hypoxia gradient (HIF target panel) per section — Non-Responder sections show broader hypoxic cores."),
    ("spatial_2_cd8_exhaustion_hubs", "Spatial CD8 exhaustion hubs (PDCD1/HAVCR2/TOX co-expression)."),
    ("spatial_3_macrophage_polarization", "Spatial M1/M2 macrophage polarisation map."),
    ("spatial_4_io_response_prediction", "Per-spot IO-response prediction using the discovery composition model — illustrates spatial heterogeneity of predicted response within sections."),
    # ── 05 multilayer ───────────────────────────────────────────────────
    ("05a_cnv_RvsNR", "L1 CNV summaries R vs NR (mean burden, SD, malignant fraction, clone diversity) — the manuscript's central tumour-intrinsic axis."),
    ("05a_multilayer_correlations", "Cross-layer feature correlation heatmap across 24 patients."),
    ("05a_multilayer_prevalence", "Carrier-prevalence panel per multilayer feature (sample_prevalence helper)."),
    ("05b_nmf_elbow", "NMF reconstruction-error elbow over k ∈ {5,7,10,12,15,20}; k=10 selected."),
    ("05b_nmf_program_loadings", "Top-gene loadings for each of the 10 NMF programmes."),
    ("05b_nmf_program_vs_response", "Per-patient NMF activity vs R/NR with BH p-values — NMF10 BH-p=0.0054."),
    ("05b_nmf_spatial_projection", "NMF program activity projected onto Visium spots."),
    ("05c_roc_by_layer", "Per-layer LOOCV ROC curves (full multilayer set, prior to balanced subsetting)."),
    ("05c_shap_beeswarm", "SHAP beeswarm for the multilayer logistic-regression classifier."),
    ("05c_shap_layer_bar", "SHAP layer-attribution bars summing per-layer contribution."),
    ("05c_variance_partition", "Variance partition of the multilayer predictor (common vs unique Brier improvement)."),
    ("05d_balanced_feature_effects", "Balanced top-3-per-layer feature-effect plot (logistic-regression coefficient + sign)."),
    ("05d_balanced_roc", "Balanced top-3-per-layer LOOCV ROC curves per layer."),
    ("05d_balanced_variance_partition", "Variance partition for the balanced ensemble."),
    ("05e_auc_bars", "Per-layer AUC bars including balanced ensemble and removing-one-layer drops."),
    ("05e_cnv_RvsNR", "CNV R-vs-NR panel matching the manuscript Figure 2d violins."),
    ("05e_feature_effects", "Feature-effect bars for the selected balanced model."),
    ("05e_nmf_bulk_vs_sc_concordance", "Bulk vs single-cell NMF activity concordance per programme."),
    ("05e_nmf_program_vs_response", "NMF program activity vs response (alternative panel)."),
    ("05e_nmf_spatial_heatmap", "NMF spatial heatmap per section."),
    ("05e_nmf_top_genes_dot", "Top NMF genes dot plot."),
    ("05e_roc_by_layer", "ROC by layer in the e-stage revision."),
    ("05e_variance_partition", "Variance partition for the e-stage revision."),
    ("05g_cnv_chr_arm_heatmap", "CNV chromosome-arm × patient heatmap — basis for the Responder-gained arm analysis (chr5/6/7/9)."),
    ("05g_cnv_chr_prevalence", "Per-chromosome carrier-prevalence panel."),
    ("05g_cnv_pathway_enrichment", "CNV-mapped gene-set enrichment (Fisher exact + BH) — HALLMARK_INTERFERON_GAMMA_RESPONSE BH-p=0.128 in the final (corrected) build; nominal only."),
    ("05g_cnv_top_gene_strip", "Top R-gained / R-lost CNV-mapped genes as a strip plot."),
    ("05g_cnv_trajectory", "CNV → NMF trajectory: ordering patients by `cnv_clone_diversity` and tracking NMF10 / CD8 fraction."),
    ("05g_cnv_volcano_bins", "Per-bin CNV volcano (Mann–Whitney + BH, 24 patients × 1,265 bins)."),
    ("05h_cnv_chain", "CNV mediation chain — CNV clone diversity → exhaustion / Treg signatures → predicted response."),
    ("05i_geneformer_layer_benchmark", "Geneformer (v2 104M, clcancer-fine-tuned/-frozen/-subsample) layer-benchmark for response prediction — none of the Geneformer variants outperform the balanced multilayer ensemble at this cohort size."),
    ("05i_nmf10_plasma_diagnostic", "NMF10 plasma-cell diagnostic — confirms IGHG/IGKC/MZB1 loading is plasma-driven rather than a generic immune program."),
    ("05i_nmf10_transfer", "NMF10 transfer panel to validation cohorts (Cappuyns; spatial)."),
    ("05j_finetune_loocv_summary", "Fine-tuned LOOCV summary for the response classifier."),
    ("05j_stratified_null", "500-iteration stratified-null distribution for the multilayer ensemble (L3 p_strat=0.002; ensemble p_strat=0.006)."),
    ("05k_nmf10_cappuyns_transfer", "NMF10 cross-regimen transfer to Cappuyns — modest direction-consistent shift."),
    ("05l_cnv_exhaustion_chain", "CNV → exhaustion mediation chain (extended version of 05h)."),
    ("05n_cappuyns_jackknife", "Patient jackknife on Cappuyns — confirms validation transferability is not driven by a single outlier."),
    ("05n_feature_transport_qc", "Feature-transport QC: which discovery features survive distribution-shift to Cappuyns."),
    ("05n_nested_layer_metrics", "Nested LOOCV per-layer metrics with selection done inside the inner fold (guards against optimistic AUCs from full-cohort feature selection)."),
    ("05z_patient_strips_L3_Composition", "Per-patient strip plots of L3 composition features."),
    ("05z_patient_strips_L4_IO_Signatures", "Per-patient strip plots of L4 IO signatures."),
    ("05z_patient_strips_L5_NMF_Programs", "Per-patient strip plots of L5 NMF programmes."),
    ("05z_spatial_4panel_RvsNR", "Four-panel spatial R-vs-NR composite for the reframe pass."),
    ("05z_spatial_L3_CD8_T", "Spatial CD8 T-cell composition maps per section."),
    ("05z_spatial_L3_cDC2", "Spatial cDC2 composition maps."),
    ("05z_spatial_L3_Macrophage", "Spatial macrophage composition maps."),
    ("05z_spatial_L3_Plasma", "Spatial plasma-cell composition maps — confirms NMF10 plasma-cell driver in tissue."),
    ("05z_spatial_NMF10", "Spatial NMF10 activity per section."),
    ("05z_spatial_NMF4", "Spatial NMF4 (cytotoxic T) activity per section."),
    ("05z_spatial_NMF5", "Spatial NMF5 activity per section."),
    ("05z_spatial_sig_Angiogenesis", "Spatial angiogenesis-signature maps."),
    ("05z_spatial_sig_IFN_gamma", "Spatial IFN-γ-response maps."),
    ("05z_spatial_sig_Immune_checkpoint", "Spatial immune-checkpoint signature maps."),
    ("05z_spatial_sig_T_cell_exhaustion", "Spatial T-cell exhaustion maps."),
    ("05z_spatial_sig_Tfh_like", "Spatial Tfh-like CD4 maps."),
    ("05z_spatial_sig_TGF_beta", "Spatial TGF-β maps."),
    ("05z_summary_heatmap", "Reframe-stage summary heatmap across L3/L4/L5 features."),
    ("05z_umap_program_overlays", "UMAP overlays of NMF programmes for the reframe pass."),
    ("05z_visium_section_RvsNR", "Per-section R-vs-NR composite Visium overlay."),
    ("05z_traj_cnv_nmf10", "Polished trajectory: CNV clone diversity → NMF10 → response."),
    ("05z_traj_scrna_pseudotime", "Polished scRNA pseudotime trajectory."),
    ("05z_traj_spatial_gradient", "Polished spatial gradient trajectory."),
    ("05z_traj_umap", "Polished trajectory UMAP."),
]


def describe(name: str) -> str:
    base = re.sub(r"\.png$", "", name, flags=re.I)
    for key, desc in FIGURE_DESCRIPTIONS:
        if key in base:
            return desc
    return (
        "Auxiliary diagnostic figure produced by the corresponding analysis script. "
        "See `scripts/` and `results_index.yaml` for provenance."
    )


def collect_pngs(folder: Path) -> list[Path]:
    return sorted(p for p in folder.rglob("*.png"))


def folder_intro(name: str) -> str:
    return FOLDER_INTROS.get(name, "Analysis-stage outputs.")


def build_md(folder: Path) -> str:
    name = folder.name
    pngs = collect_pngs(folder)
    header = dedent(
        f"""\
        # Figures interpretation — `out/{name}`

        > Companion to `docs/manuscript.md` (May 2026 revision). All R-vs-NR
        > statistics, AUCs and effect sizes referenced below are computed at the
        > patient level on the GSE206325 anti-PD-1 discovery cohort (n=24; 7 R,
        > 17 NR) unless otherwise stated. External validation uses the Cappuyns
        > et al. 2023 atezolizumab + bevacizumab cohort (n=32 evaluable; 19 R,
        > 13 NR). Statistics use two-sided Mann–Whitney U with BH-FDR within
        > family, Cliff's δ for effect size, and 2,000-iteration patient-level
        > bootstrap 95 % CIs.

        ## Scope
        {folder_intro(name)}

        ## Figure-by-figure interpretation
        """
    )

    sections: list[str] = []
    if not pngs:
        return header + "\nNo PNG figures in this folder.\n"

    # Group per immediate subdirectory for readability
    grouped: dict[str, list[Path]] = {}
    for p in pngs:
        rel = p.relative_to(folder)
        sub = str(rel.parent) if rel.parent != Path(".") else "(root)"
        grouped.setdefault(sub, []).append(p)

    for sub in sorted(grouped):
        if sub != "(root)":
            sections.append(f"\n### `{sub}/`\n")
        for p in grouped[sub]:
            rel = p.relative_to(folder).as_posix()
            sections.append(
                f"#### `{rel}`\n\n{describe(p.name)}\n\n"
                f"![{p.stem}]({rel})\n"
            )
    return header + "".join(sections)


def write_docx(md_path: Path, docx_path: Path) -> bool:
    try:
        subprocess.run(
            ["pandoc", str(md_path), "-o", str(docx_path),
             "--resource-path", str(md_path.parent),
             "--from", "gfm", "--to", "docx"],
            check=True, capture_output=True,
        )
        return True
    except subprocess.CalledProcessError as e:
        print(f"  ! pandoc failed for {md_path}: {e.stderr.decode()[:200]}")
        return False


def main() -> None:
    targets = [d for d in sorted(OUT.iterdir()) if d.is_dir()]
    for folder in targets:
        md_path = folder / "figures_interpretation.md"
        docx_path = folder / "figures_interpretation.docx"
        md = build_md(folder)
        md_path.write_text(md)
        ok = write_docx(md_path, docx_path)
        n = len(collect_pngs(folder))
        flag = "✓" if ok else "✗docx"
        print(f"{flag} out/{folder.name}/ ({n} png) → figures_interpretation.{{md,docx}}")


if __name__ == "__main__":
    main()
