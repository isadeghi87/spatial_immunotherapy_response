# Figures interpretation — `out/05_multilayer`

> Companion to `docs/manuscript.md` (May 2026 revision). All R-vs-NR
> statistics, AUCs and effect sizes referenced below are computed at the
> patient level on the GSE206325 anti-PD-1 discovery cohort (n=24; 7 R,
> 17 NR) unless otherwise stated. External validation uses the Cappuyns
> et al. 2023 atezolizumab + bevacizumab cohort (n=32 evaluable; 19 R,
> 13 NR). Statistics use two-sided Mann–Whitney U with BH-FDR within
> family, Cliff's δ for effect size, and 2,000-iteration patient-level
> bootstrap 95 % CIs.

## Scope
Five-layer predictive framework: L1 tumour-intrinsic CNV, L2 regulons, L3 cell-type composition, L4 IO signatures, L5 NMF programs. Includes balanced top-3-features-per-layer LOOCV (24-patient discovery), 1,000-iteration permutation null, 2,000-iteration bootstrap CIs, SHAP attribution, NMF program forest plots, CNV chromosome/pathway analyses, Geneformer benchmark comparisons, and held-out Cappuyns transfer.

## Figure-by-figure interpretation

### `figures/`
#### `figures/05a_cnv_RvsNR.png`

L1 CNV summaries R vs NR (mean burden, SD, malignant fraction, clone diversity) — the manuscript's central tumour-intrinsic axis.

![05a_cnv_RvsNR](figures/05a_cnv_RvsNR.png)
#### `figures/05a_multilayer_correlations.png`

Cross-layer feature correlation heatmap across 24 patients.

![05a_multilayer_correlations](figures/05a_multilayer_correlations.png)
#### `figures/05a_multilayer_prevalence.png`

Carrier-prevalence panel per multilayer feature (sample_prevalence helper).

![05a_multilayer_prevalence](figures/05a_multilayer_prevalence.png)
#### `figures/05b_nmf_elbow.png`

NMF reconstruction-error elbow over k ∈ {5,7,10,12,15,20}; k=10 selected.

![05b_nmf_elbow](figures/05b_nmf_elbow.png)
#### `figures/05b_nmf_program_loadings.png`

Top-gene loadings for each of the 10 NMF programmes.

![05b_nmf_program_loadings](figures/05b_nmf_program_loadings.png)
#### `figures/05b_nmf_program_vs_response.png`

Per-patient NMF activity vs R/NR with BH p-values — NMF10 BH-p=0.0054.

![05b_nmf_program_vs_response](figures/05b_nmf_program_vs_response.png)
#### `figures/05b_nmf_spatial_projection.png`

NMF program activity projected onto Visium spots.

![05b_nmf_spatial_projection](figures/05b_nmf_spatial_projection.png)
#### `figures/05c_roc_by_layer.png`

Per-layer LOOCV ROC curves (full multilayer set, prior to balanced subsetting).

![05c_roc_by_layer](figures/05c_roc_by_layer.png)
#### `figures/05c_shap_beeswarm.png`

SHAP beeswarm for the multilayer logistic-regression classifier.

![05c_shap_beeswarm](figures/05c_shap_beeswarm.png)
#### `figures/05c_shap_layer_bar.png`

SHAP layer-attribution bars summing per-layer contribution.

![05c_shap_layer_bar](figures/05c_shap_layer_bar.png)
#### `figures/05c_variance_partition.png`

Variance partition of the multilayer predictor (common vs unique Brier improvement).

![05c_variance_partition](figures/05c_variance_partition.png)
#### `figures/05d_balanced_feature_effects.png`

Balanced top-3-per-layer feature-effect plot (logistic-regression coefficient + sign).

![05d_balanced_feature_effects](figures/05d_balanced_feature_effects.png)
#### `figures/05d_balanced_roc.png`

Balanced top-3-per-layer LOOCV ROC curves per layer.

![05d_balanced_roc](figures/05d_balanced_roc.png)
#### `figures/05d_balanced_variance_partition.png`

Variance partition for the balanced ensemble.

![05d_balanced_variance_partition](figures/05d_balanced_variance_partition.png)
#### `figures/05e_auc_bars.png`

Per-layer AUC bars including balanced ensemble and removing-one-layer drops.

![05e_auc_bars](figures/05e_auc_bars.png)
#### `figures/05e_cnv_RvsNR.png`

CNV R-vs-NR panel matching the manuscript Figure 2d violins.

![05e_cnv_RvsNR](figures/05e_cnv_RvsNR.png)
#### `figures/05e_feature_effects.png`

Feature-effect bars for the selected balanced model.

![05e_feature_effects](figures/05e_feature_effects.png)
#### `figures/05e_nmf_bulk_vs_sc_concordance.png`

Bulk vs single-cell NMF activity concordance per programme.

![05e_nmf_bulk_vs_sc_concordance](figures/05e_nmf_bulk_vs_sc_concordance.png)
#### `figures/05e_nmf_program_vs_response.png`

NMF program activity vs response (alternative panel).

![05e_nmf_program_vs_response](figures/05e_nmf_program_vs_response.png)
#### `figures/05e_nmf_spatial_heatmap.png`

NMF spatial heatmap per section.

![05e_nmf_spatial_heatmap](figures/05e_nmf_spatial_heatmap.png)
#### `figures/05e_nmf_top_genes_dot.png`

Top NMF genes dot plot.

![05e_nmf_top_genes_dot](figures/05e_nmf_top_genes_dot.png)
#### `figures/05e_roc_by_layer.png`

ROC by layer in the e-stage revision.

![05e_roc_by_layer](figures/05e_roc_by_layer.png)
#### `figures/05e_variance_partition.png`

Variance partition for the e-stage revision.

![05e_variance_partition](figures/05e_variance_partition.png)
#### `figures/05g_cnv_chr_arm_heatmap.png`

CNV chromosome-arm × patient heatmap — basis for the Responder-gained arm analysis (chr5/6/7/9).

![05g_cnv_chr_arm_heatmap](figures/05g_cnv_chr_arm_heatmap.png)
#### `figures/05g_cnv_chr_prevalence.png`

Per-chromosome carrier-prevalence panel.

![05g_cnv_chr_prevalence](figures/05g_cnv_chr_prevalence.png)
#### `figures/05g_cnv_pathway_enrichment.png`

CNV-mapped gene-set enrichment (Fisher exact + BH) — HALLMARK_INTERFERON_GAMMA_RESPONSE BH-p=0.128 in the final (corrected) build; nominal only.

![05g_cnv_pathway_enrichment](figures/05g_cnv_pathway_enrichment.png)
#### `figures/05g_cnv_top_gene_strip.png`

Top R-gained / R-lost CNV-mapped genes as a strip plot.

![05g_cnv_top_gene_strip](figures/05g_cnv_top_gene_strip.png)
#### `figures/05g_cnv_trajectory.png`

CNV → NMF trajectory: ordering patients by `cnv_clone_diversity` and tracking NMF10 / CD8 fraction.

![05g_cnv_trajectory](figures/05g_cnv_trajectory.png)
#### `figures/05g_cnv_volcano_bins.png`

Per-bin CNV volcano (Mann–Whitney + BH, 24 patients × 1,265 bins).

![05g_cnv_volcano_bins](figures/05g_cnv_volcano_bins.png)
#### `figures/05h_cnv_chain.png`

CNV mediation chain — CNV clone diversity → exhaustion / Treg signatures → predicted response.

![05h_cnv_chain](figures/05h_cnv_chain.png)
#### `figures/05i_geneformer_layer_benchmark_geneformer_v2_104m_clcancer_finetuned.png`

Geneformer (v2 104M, clcancer-fine-tuned/-frozen/-subsample) layer-benchmark for response prediction — none of the Geneformer variants outperform the balanced multilayer ensemble at this cohort size.

![05i_geneformer_layer_benchmark_geneformer_v2_104m_clcancer_finetuned](figures/05i_geneformer_layer_benchmark_geneformer_v2_104m_clcancer_finetuned.png)
#### `figures/05i_geneformer_layer_benchmark_geneformer_v2_104m_clcancer_frozen_200cpp.png`

Geneformer (v2 104M, clcancer-fine-tuned/-frozen/-subsample) layer-benchmark for response prediction — none of the Geneformer variants outperform the balanced multilayer ensemble at this cohort size.

![05i_geneformer_layer_benchmark_geneformer_v2_104m_clcancer_frozen_200cpp](figures/05i_geneformer_layer_benchmark_geneformer_v2_104m_clcancer_frozen_200cpp.png)
#### `figures/05i_geneformer_layer_benchmark_geneformer_v2_104m_clcancer_subsample_200cpp.png`

Geneformer (v2 104M, clcancer-fine-tuned/-frozen/-subsample) layer-benchmark for response prediction — none of the Geneformer variants outperform the balanced multilayer ensemble at this cohort size.

![05i_geneformer_layer_benchmark_geneformer_v2_104m_clcancer_subsample_200cpp](figures/05i_geneformer_layer_benchmark_geneformer_v2_104m_clcancer_subsample_200cpp.png)
#### `figures/05i_geneformer_layer_benchmark_geneformer_v2_104m_clcancer_subsample_2cpp.png`

Geneformer (v2 104M, clcancer-fine-tuned/-frozen/-subsample) layer-benchmark for response prediction — none of the Geneformer variants outperform the balanced multilayer ensemble at this cohort size.

![05i_geneformer_layer_benchmark_geneformer_v2_104m_clcancer_subsample_2cpp](figures/05i_geneformer_layer_benchmark_geneformer_v2_104m_clcancer_subsample_2cpp.png)
#### `figures/05i_nmf10_plasma_diagnostic.png`

NMF10 plasma-cell diagnostic — confirms IGHG/IGKC/MZB1 loading is plasma-driven rather than a generic immune program.

![05i_nmf10_plasma_diagnostic](figures/05i_nmf10_plasma_diagnostic.png)
#### `figures/05i_nmf10_transfer.png`

NMF10 transfer panel to validation cohorts (Cappuyns; spatial).

![05i_nmf10_transfer](figures/05i_nmf10_transfer.png)
#### `figures/05j_finetune_loocv_summary.png`

Fine-tuned LOOCV summary for the response classifier.

![05j_finetune_loocv_summary](figures/05j_finetune_loocv_summary.png)
#### `figures/05j_stratified_null.png`

500-iteration stratified-null distribution for the multilayer ensemble (L3 p_strat=0.002; ensemble p_strat=0.006).

![05j_stratified_null](figures/05j_stratified_null.png)
#### `figures/05k_nmf10_cappuyns_transfer.png`

NMF10 cross-regimen transfer to Cappuyns — modest direction-consistent shift.

![05k_nmf10_cappuyns_transfer](figures/05k_nmf10_cappuyns_transfer.png)
#### `figures/05l_cnv_exhaustion_chain.png`

CNV → exhaustion mediation chain (extended version of 05h).

![05l_cnv_exhaustion_chain](figures/05l_cnv_exhaustion_chain.png)
#### `figures/05n_cappuyns_jackknife.png`

Patient jackknife on Cappuyns — confirms validation transferability is not driven by a single outlier.

![05n_cappuyns_jackknife](figures/05n_cappuyns_jackknife.png)
#### `figures/05n_feature_transport_qc.png`

Feature-transport QC: which discovery features survive distribution-shift to Cappuyns.

![05n_feature_transport_qc](figures/05n_feature_transport_qc.png)
#### `figures/05n_nested_layer_metrics.png`

Nested LOOCV per-layer metrics with selection done inside the inner fold (guards against optimistic AUCs from full-cohort feature selection).

![05n_nested_layer_metrics](figures/05n_nested_layer_metrics.png)
#### `figures/Fig6A_balanced_layer_metrics.png`

Updated Figure 6 panel A — balanced LOOCV metrics with ROC-AUC, PR-AUC and balanced-accuracy bars per layer plus the 15-feature ensemble. Permutation-null 95-th percentiles are marked in red for L3 and the ensemble (the only entries exceeding their null at p<0.05). PR-AUC prevalence baseline = 0.29.

![Fig6A_balanced_layer_metrics](figures/Fig6A_balanced_layer_metrics.png)

### `figures/05z_expression_spatial/`
#### `figures/05z_expression_spatial/05z_patient_strips_L3_Composition.png`

Per-patient strip plots of L3 composition features.

![05z_patient_strips_L3_Composition](figures/05z_expression_spatial/05z_patient_strips_L3_Composition.png)
#### `figures/05z_expression_spatial/05z_patient_strips_L4_IO_Signatures.png`

Per-patient strip plots of L4 IO signatures.

![05z_patient_strips_L4_IO_Signatures](figures/05z_expression_spatial/05z_patient_strips_L4_IO_Signatures.png)
#### `figures/05z_expression_spatial/05z_patient_strips_L5_NMF_Programs.png`

Per-patient strip plots of L5 NMF programmes.

![05z_patient_strips_L5_NMF_Programs](figures/05z_expression_spatial/05z_patient_strips_L5_NMF_Programs.png)
#### `figures/05z_expression_spatial/05z_spatial_4panel_RvsNR.png`

Four-panel spatial R-vs-NR composite for the reframe pass.

![05z_spatial_4panel_RvsNR](figures/05z_expression_spatial/05z_spatial_4panel_RvsNR.png)
#### `figures/05z_expression_spatial/05z_spatial_L3_CD8_T.png`

Spatial CD8 T-cell composition maps per section.

![05z_spatial_L3_CD8_T](figures/05z_expression_spatial/05z_spatial_L3_CD8_T.png)
#### `figures/05z_expression_spatial/05z_spatial_L3_Macrophage.png`

Spatial macrophage composition maps.

![05z_spatial_L3_Macrophage](figures/05z_expression_spatial/05z_spatial_L3_Macrophage.png)
#### `figures/05z_expression_spatial/05z_spatial_L3_Plasma.png`

Spatial plasma-cell composition maps — confirms NMF10 plasma-cell driver in tissue.

![05z_spatial_L3_Plasma](figures/05z_expression_spatial/05z_spatial_L3_Plasma.png)
#### `figures/05z_expression_spatial/05z_spatial_L3_cDC2.png`

Spatial cDC2 composition maps.

![05z_spatial_L3_cDC2](figures/05z_expression_spatial/05z_spatial_L3_cDC2.png)
#### `figures/05z_expression_spatial/05z_spatial_NMF10.png`

Spatial NMF10 activity per section.

![05z_spatial_NMF10](figures/05z_expression_spatial/05z_spatial_NMF10.png)
#### `figures/05z_expression_spatial/05z_spatial_NMF4.png`

Spatial NMF4 (cytotoxic T) activity per section.

![05z_spatial_NMF4](figures/05z_expression_spatial/05z_spatial_NMF4.png)
#### `figures/05z_expression_spatial/05z_spatial_NMF5.png`

Spatial NMF5 activity per section.

![05z_spatial_NMF5](figures/05z_expression_spatial/05z_spatial_NMF5.png)
#### `figures/05z_expression_spatial/05z_spatial_sig_Angiogenesis.png`

Spatial angiogenesis-signature maps.

![05z_spatial_sig_Angiogenesis](figures/05z_expression_spatial/05z_spatial_sig_Angiogenesis.png)
#### `figures/05z_expression_spatial/05z_spatial_sig_IFN_gamma.png`

Spatial IFN-γ-response maps.

![05z_spatial_sig_IFN_gamma](figures/05z_expression_spatial/05z_spatial_sig_IFN_gamma.png)
#### `figures/05z_expression_spatial/05z_spatial_sig_Immune_checkpoint.png`

Spatial immune-checkpoint signature maps.

![05z_spatial_sig_Immune_checkpoint](figures/05z_expression_spatial/05z_spatial_sig_Immune_checkpoint.png)
#### `figures/05z_expression_spatial/05z_spatial_sig_TGF_beta.png`

Spatial TGF-β maps.

![05z_spatial_sig_TGF_beta](figures/05z_expression_spatial/05z_spatial_sig_TGF_beta.png)
#### `figures/05z_expression_spatial/05z_spatial_sig_T_cell_exhaustion.png`

Spatial T-cell exhaustion maps.

![05z_spatial_sig_T_cell_exhaustion](figures/05z_expression_spatial/05z_spatial_sig_T_cell_exhaustion.png)
#### `figures/05z_expression_spatial/05z_spatial_sig_Tfh_like.png`

Spatial Tfh-like CD4 maps.

![05z_spatial_sig_Tfh_like](figures/05z_expression_spatial/05z_spatial_sig_Tfh_like.png)
#### `figures/05z_expression_spatial/05z_summary_heatmap.png`

Reframe-stage summary heatmap across L3/L4/L5 features.

![05z_summary_heatmap](figures/05z_expression_spatial/05z_summary_heatmap.png)
#### `figures/05z_expression_spatial/05z_umap_program_overlays.png`

UMAP overlays of NMF programmes for the reframe pass.

![05z_umap_program_overlays](figures/05z_expression_spatial/05z_umap_program_overlays.png)
#### `figures/05z_expression_spatial/05z_visium_section_RvsNR.png`

Per-section R-vs-NR composite Visium overlay.

![05z_visium_section_RvsNR](figures/05z_expression_spatial/05z_visium_section_RvsNR.png)

### `figures/05z_trajectory_polished/`
#### `figures/05z_trajectory_polished/05z_traj_cnv_nmf10.png`

Polished trajectory: CNV clone diversity → NMF10 → response.

![05z_traj_cnv_nmf10](figures/05z_trajectory_polished/05z_traj_cnv_nmf10.png)
#### `figures/05z_trajectory_polished/05z_traj_scrna_pseudotime.png`

Polished scRNA pseudotime trajectory.

![05z_traj_scrna_pseudotime](figures/05z_trajectory_polished/05z_traj_scrna_pseudotime.png)
#### `figures/05z_trajectory_polished/05z_traj_spatial_gradient.png`

Polished spatial gradient trajectory.

![05z_traj_spatial_gradient](figures/05z_trajectory_polished/05z_traj_spatial_gradient.png)
#### `figures/05z_trajectory_polished/05z_traj_umap.png`

Polished trajectory UMAP.

![05z_traj_umap](figures/05z_trajectory_polished/05z_traj_umap.png)
