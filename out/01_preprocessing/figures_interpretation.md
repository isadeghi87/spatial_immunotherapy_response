# Figures interpretation — `out/01_preprocessing`

> Companion to `docs/manuscript.md` (May 2026 revision). All R-vs-NR
> statistics, AUCs and effect sizes referenced below are computed at the
> patient level on the GSE206325 anti-PD-1 discovery cohort (n=24; 7 R,
> 17 NR) unless otherwise stated. External validation uses the Cappuyns
> et al. 2023 atezolizumab + bevacizumab cohort (n=32 evaluable; 19 R,
> 13 NR). Statistics use two-sided Mann–Whitney U with BH-FDR within
> family, Cliff's δ for effect size, and 2,000-iteration patient-level
> bootstrap 95 % CIs.

## Scope
Stage-1 figures cover atlas QC, per-dataset integration, scVI/Harmony benchmarking, and reference cell-type annotation across the 1,058,964-cell integrated HCC scRNA-seq atlas (5 HCC cohorts + non-tumour liver reference).

## Figure-by-figure interpretation

### `figures/`
#### `figures/01_celltype_umap.png`

UMAP coloured by the 44-class `cell_type_original` annotation used throughout the manuscript.

![01_celltype_umap](figures/01_celltype_umap.png)
#### `figures/01_immune_umaps.png`

Immune-restricted UMAPs split by lineage.

![01_immune_umaps](figures/01_immune_umaps.png)
#### `figures/01_integration_umap.png`

Pre- vs post-Harmony UMAP coloured by dataset and donor — visual confirmation that Harmony removes batch while preserving cell-type structure.

![01_integration_umap](figures/01_integration_umap.png)
#### `figures/01_leiden_clusters.png`

Leiden clustering at resolution 0.8 on the harmonised 15-NN graph.

![01_leiden_clusters](figures/01_leiden_clusters.png)
#### `figures/01_marker_feature_plots.png`

Marker-gene feature plots overlaying canonical genes on the harmonised UMAP.

![01_marker_feature_plots](figures/01_marker_feature_plots.png)
#### `figures/01_qc_violins.png`

Per-cell QC violins (n_genes_by_counts, total_counts, pct_counts_mt) before and after Scanpy QC thresholds (200 ≤ n_genes ≤ 8000, pct_mt < 15%).

![01_qc_violins](figures/01_qc_violins.png)
#### `figures/01_scanvi_confidence.png`

scANVI label-transfer confidence per cell — supports the manual annotation against an automated baseline.

![01_scanvi_confidence](figures/01_scanvi_confidence.png)
#### `figures/01c_batch_metrics_summary.png`

Bar summary of ASW-batch / ASW-bio / iLISI / cLISI / kBET for PCA vs Harmony (matches §2.16 of M&M).

![01c_batch_metrics_summary](figures/01c_batch_metrics_summary.png)
#### `figures/01c_lisi_histograms.png`

Local Inverse-Simpson histograms: iLISI (dataset) and cLISI (cell type) before and after integration.

![01c_lisi_histograms](figures/01c_lisi_histograms.png)
#### `figures/01c_umap_integration_comparison.png`

Side-by-side UMAPs from PCA, Harmony, and scVI embeddings.

![01c_umap_integration_comparison](figures/01c_umap_integration_comparison.png)
#### `figures/01d_scvi_metric_bars.png`

scVI metric summary bars vs Harmony baseline.

![01d_scvi_metric_bars](figures/01d_scvi_metric_bars.png)
#### `figures/01d_scvi_umap.png`

scVI latent UMAP used as an orthogonal integration cross-check.

![01d_scvi_umap](figures/01d_scvi_umap.png)
#### `figures/01y_proportion_by_dataset.png`

Per-dataset cell-type proportion bars — sanity check for relative lineage capture across cohorts.

![01y_proportion_by_dataset](figures/01y_proportion_by_dataset.png)
#### `figures/01y_proportion_cappuyns.png`

Cell-type proportion bars in the Cappuyns 2023 external validation cohort.

![01y_proportion_cappuyns](figures/01y_proportion_cappuyns.png)
#### `figures/01y_proportion_gse206325.png`

Cell-type proportion bars in the GSE206325 anti-PD-1 discovery cohort.

![01y_proportion_gse206325](figures/01y_proportion_gse206325.png)
#### `figures/01y_umap_celltype.png`

UMAP coloured by cell-type after Harmony.

![01y_umap_celltype](figures/01y_umap_celltype.png)
#### `figures/01y_umap_dataset.png`

UMAP coloured by dataset across all cohorts in the integrated atlas.

![01y_umap_dataset](figures/01y_umap_dataset.png)
#### `figures/01y_umap_response.png`

UMAP coloured by RECIST response across the 24-patient discovery subset.

![01y_umap_response](figures/01y_umap_response.png)
#### `figures/01y_visium_celltype.png`

Visium spot-level cell-type calls from Cell2location deconvolution.

![01y_visium_celltype](figures/01y_visium_celltype.png)
#### `figures/01y_visium_immunotype.png`

Per-spot immunotype assignment derived from joint IO-signature thresholds.

![01y_visium_immunotype](figures/01y_visium_immunotype.png)
