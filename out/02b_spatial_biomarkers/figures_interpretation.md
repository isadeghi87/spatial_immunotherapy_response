# Figures interpretation — `out/02b_spatial_biomarkers`

> Companion to `docs/manuscript.md` (May 2026 revision). All R-vs-NR
> statistics, AUCs and effect sizes referenced below are computed at the
> patient level on the GSE206325 anti-PD-1 discovery cohort (n=24; 7 R,
> 17 NR) unless otherwise stated. External validation uses the Cappuyns
> et al. 2023 atezolizumab + bevacizumab cohort (n=32 evaluable; 19 R,
> 13 NR). Statistics use two-sided Mann–Whitney U with BH-FDR within
> family, Cliff's δ for effect size, and 2,000-iteration patient-level
> bootstrap 95 % CIs.

## Scope
Visium spatial transcriptomics analyses across GSE173905, GSE238264 (4R/3NR) and GSE281759, combined with a 24-section scVI-integrated spatial atlas (43,277 spots × 30 latent dimensions) and per-sample InferCNV-based tumour delineation.

## Figure-by-figure interpretation

### `cnv/figures/`
#### `cnv/figures/02g_spatial_cnv_arm_heatmap.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![02g_spatial_cnv_arm_heatmap](cnv/figures/02g_spatial_cnv_arm_heatmap.png)
#### `cnv/figures/02g_spatial_cnv_genome_heatmap_celltype.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![02g_spatial_cnv_genome_heatmap_celltype](cnv/figures/02g_spatial_cnv_genome_heatmap_celltype.png)
#### `cnv/figures/spatial_arm_burden_RvsNR.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![spatial_arm_burden_RvsNR](cnv/figures/spatial_arm_burden_RvsNR.png)
#### `cnv/figures/spatial_cnv_GSM8627601_A009N.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627601_A009N](cnv/figures/spatial_cnv_GSM8627601_A009N.png)
#### `cnv/figures/spatial_cnv_GSM8627602_A009T1.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627602_A009T1](cnv/figures/spatial_cnv_GSM8627602_A009T1.png)
#### `cnv/figures/spatial_cnv_GSM8627603_A009T2.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627603_A009T2](cnv/figures/spatial_cnv_GSM8627603_A009T2.png)
#### `cnv/figures/spatial_cnv_GSM8627604_A009T3.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627604_A009T3](cnv/figures/spatial_cnv_GSM8627604_A009T3.png)
#### `cnv/figures/spatial_cnv_GSM8627605_A009T4.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627605_A009T4](cnv/figures/spatial_cnv_GSM8627605_A009T4.png)
#### `cnv/figures/spatial_cnv_GSM8627606_A009T5.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627606_A009T5](cnv/figures/spatial_cnv_GSM8627606_A009T5.png)
#### `cnv/figures/spatial_cnv_GSM8627607_B0013T1.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627607_B0013T1](cnv/figures/spatial_cnv_GSM8627607_B0013T1.png)
#### `cnv/figures/spatial_cnv_GSM8627608_B0013T2.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627608_B0013T2](cnv/figures/spatial_cnv_GSM8627608_B0013T2.png)
#### `cnv/figures/spatial_cnv_GSM8627609_B0013T3.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627609_B0013T3](cnv/figures/spatial_cnv_GSM8627609_B0013T3.png)
#### `cnv/figures/spatial_cnv_GSM8627610_B0013T4.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627610_B0013T4](cnv/figures/spatial_cnv_GSM8627610_B0013T4.png)
#### `cnv/figures/spatial_cnv_GSM8627611_B013T5.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627611_B013T5](cnv/figures/spatial_cnv_GSM8627611_B013T5.png)
#### `cnv/figures/spatial_cnv_GSM8627612_C004T1P.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627612_C004T1P](cnv/figures/spatial_cnv_GSM8627612_C004T1P.png)
#### `cnv/figures/spatial_cnv_GSM8627613_C004T1R.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627613_C004T1R](cnv/figures/spatial_cnv_GSM8627613_C004T1R.png)
#### `cnv/figures/spatial_cnv_GSM8627614_C008T1.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627614_C008T1](cnv/figures/spatial_cnv_GSM8627614_C008T1.png)
#### `cnv/figures/spatial_cnv_GSM8627615_C008T2.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627615_C008T2](cnv/figures/spatial_cnv_GSM8627615_C008T2.png)
#### `cnv/figures/spatial_cnv_GSM8627616_C008T3.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627616_C008T3](cnv/figures/spatial_cnv_GSM8627616_C008T3.png)
#### `cnv/figures/spatial_cnv_GSM8627617_C008T4.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_GSM8627617_C008T4](cnv/figures/spatial_cnv_GSM8627617_C008T4.png)
#### `cnv/figures/spatial_cnv_HCC1R.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_HCC1R](cnv/figures/spatial_cnv_HCC1R.png)
#### `cnv/figures/spatial_cnv_HCC2R.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_HCC2R](cnv/figures/spatial_cnv_HCC2R.png)
#### `cnv/figures/spatial_cnv_HCC3R.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_HCC3R](cnv/figures/spatial_cnv_HCC3R.png)
#### `cnv/figures/spatial_cnv_HCC4R.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_HCC4R](cnv/figures/spatial_cnv_HCC4R.png)
#### `cnv/figures/spatial_cnv_HCC6NR.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_HCC6NR](cnv/figures/spatial_cnv_HCC6NR.png)
#### `cnv/figures/spatial_cnv_HCC7NR.png`

Per-section infercnv spatial map (X_cnv per spot, projected back to Visium coordinates) — supports the per-patient classification of malignant vs reference spots used in §P1.2 of the spatial pipeline.

![spatial_cnv_HCC7NR](cnv/figures/spatial_cnv_HCC7NR.png)
#### `cnv/figures/spatial_cnv_universal_heatmap.png`

Universal infercnv heatmap across 23/24 sections (43,196 spots × 2,115 genome windows) — HCC5NR excluded (0 reference cells).

![spatial_cnv_universal_heatmap](cnv/figures/spatial_cnv_universal_heatmap.png)

### `cnv/rds/GSM8627601_A009N/`
#### `cnv/rds/GSM8627601_A009N/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627601_A009N/infercnv.png)
#### `cnv/rds/GSM8627601_A009N/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627601_A009N/infercnv.preliminary.png)

### `cnv/rds/GSM8627602_A009T1/`
#### `cnv/rds/GSM8627602_A009T1/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627602_A009T1/infercnv.png)
#### `cnv/rds/GSM8627602_A009T1/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627602_A009T1/infercnv.preliminary.png)

### `cnv/rds/GSM8627603_A009T2/`
#### `cnv/rds/GSM8627603_A009T2/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627603_A009T2/infercnv.png)
#### `cnv/rds/GSM8627603_A009T2/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627603_A009T2/infercnv.preliminary.png)

### `cnv/rds/GSM8627604_A009T3/`
#### `cnv/rds/GSM8627604_A009T3/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627604_A009T3/infercnv.png)
#### `cnv/rds/GSM8627604_A009T3/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627604_A009T3/infercnv.preliminary.png)

### `cnv/rds/GSM8627605_A009T4/`
#### `cnv/rds/GSM8627605_A009T4/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627605_A009T4/infercnv.png)
#### `cnv/rds/GSM8627605_A009T4/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627605_A009T4/infercnv.preliminary.png)

### `cnv/rds/GSM8627606_A009T5/`
#### `cnv/rds/GSM8627606_A009T5/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627606_A009T5/infercnv.png)
#### `cnv/rds/GSM8627606_A009T5/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627606_A009T5/infercnv.preliminary.png)

### `cnv/rds/GSM8627607_B0013T1/`
#### `cnv/rds/GSM8627607_B0013T1/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627607_B0013T1/infercnv.png)
#### `cnv/rds/GSM8627607_B0013T1/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627607_B0013T1/infercnv.preliminary.png)

### `cnv/rds/GSM8627608_B0013T2/`
#### `cnv/rds/GSM8627608_B0013T2/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627608_B0013T2/infercnv.png)
#### `cnv/rds/GSM8627608_B0013T2/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627608_B0013T2/infercnv.preliminary.png)

### `cnv/rds/GSM8627609_B0013T3/`
#### `cnv/rds/GSM8627609_B0013T3/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627609_B0013T3/infercnv.png)
#### `cnv/rds/GSM8627609_B0013T3/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627609_B0013T3/infercnv.preliminary.png)

### `cnv/rds/GSM8627610_B0013T4/`
#### `cnv/rds/GSM8627610_B0013T4/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627610_B0013T4/infercnv.png)
#### `cnv/rds/GSM8627610_B0013T4/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627610_B0013T4/infercnv.preliminary.png)

### `cnv/rds/GSM8627611_B013T5/`
#### `cnv/rds/GSM8627611_B013T5/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627611_B013T5/infercnv.png)
#### `cnv/rds/GSM8627611_B013T5/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627611_B013T5/infercnv.preliminary.png)

### `cnv/rds/GSM8627612_C004T1P/`
#### `cnv/rds/GSM8627612_C004T1P/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627612_C004T1P/infercnv.png)
#### `cnv/rds/GSM8627612_C004T1P/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627612_C004T1P/infercnv.preliminary.png)

### `cnv/rds/GSM8627613_C004T1R/`
#### `cnv/rds/GSM8627613_C004T1R/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627613_C004T1R/infercnv.png)
#### `cnv/rds/GSM8627613_C004T1R/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627613_C004T1R/infercnv.preliminary.png)

### `cnv/rds/GSM8627614_C008T1/`
#### `cnv/rds/GSM8627614_C008T1/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627614_C008T1/infercnv.png)
#### `cnv/rds/GSM8627614_C008T1/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627614_C008T1/infercnv.preliminary.png)

### `cnv/rds/GSM8627615_C008T2/`
#### `cnv/rds/GSM8627615_C008T2/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627615_C008T2/infercnv.png)
#### `cnv/rds/GSM8627615_C008T2/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627615_C008T2/infercnv.preliminary.png)

### `cnv/rds/GSM8627616_C008T3/`
#### `cnv/rds/GSM8627616_C008T3/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627616_C008T3/infercnv.png)
#### `cnv/rds/GSM8627616_C008T3/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627616_C008T3/infercnv.preliminary.png)

### `cnv/rds/GSM8627617_C008T4/`
#### `cnv/rds/GSM8627617_C008T4/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/GSM8627617_C008T4/infercnv.png)
#### `cnv/rds/GSM8627617_C008T4/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/GSM8627617_C008T4/infercnv.preliminary.png)

### `cnv/rds/HCC1R/`
#### `cnv/rds/HCC1R/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/HCC1R/infercnv.png)
#### `cnv/rds/HCC1R/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/HCC1R/infercnv.preliminary.png)

### `cnv/rds/HCC2R/`
#### `cnv/rds/HCC2R/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/HCC2R/infercnv.png)
#### `cnv/rds/HCC2R/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/HCC2R/infercnv.preliminary.png)

### `cnv/rds/HCC3R/`
#### `cnv/rds/HCC3R/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/HCC3R/infercnv.png)
#### `cnv/rds/HCC3R/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/HCC3R/infercnv.preliminary.png)

### `cnv/rds/HCC4R/`
#### `cnv/rds/HCC4R/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/HCC4R/infercnv.png)
#### `cnv/rds/HCC4R/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/HCC4R/infercnv.preliminary.png)

### `cnv/rds/HCC5NR/`
#### `cnv/rds/HCC5NR/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/HCC5NR/infercnv.png)
#### `cnv/rds/HCC5NR/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/HCC5NR/infercnv.preliminary.png)

### `cnv/rds/HCC6NR/`
#### `cnv/rds/HCC6NR/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/HCC6NR/infercnv.png)
#### `cnv/rds/HCC6NR/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/HCC6NR/infercnv.preliminary.png)

### `cnv/rds/HCC7NR/`
#### `cnv/rds/HCC7NR/infercnv.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![infercnv](cnv/rds/HCC7NR/infercnv.png)
#### `cnv/rds/HCC7NR/infercnv.preliminary.png`

Preliminary infercnv heatmap (pre-HMM denoising) for this section — diagnostic only.

![infercnv.preliminary](cnv/rds/HCC7NR/infercnv.preliminary.png)

### `figures/`
#### `figures/02b_ct_marker_dotplot.png`

Spatial-atlas cell-type marker dot plot from the 24-section scVI-integrated Visium analysis.

![02b_ct_marker_dotplot](figures/02b_ct_marker_dotplot.png)
#### `figures/02b_ct_marker_heatmap.png`

Spatial-atlas marker heatmap supporting the spot-level annotation.

![02b_ct_marker_heatmap](figures/02b_ct_marker_heatmap.png)
#### `figures/02b_ct_sample_annotation_support.png`

Per-sample support bars for each cell type — confirms cross-section coverage of every annotated class.

![02b_ct_sample_annotation_support](figures/02b_ct_sample_annotation_support.png)
#### `figures/02b_ct_spatial_celltype_maps.png`

Per-section spot maps coloured by dominant cell-type call.

![02b_ct_spatial_celltype_maps](figures/02b_ct_spatial_celltype_maps.png)
#### `figures/02b_ct_umap_celltype_response.png`

scVI-spatial UMAP of all 43,277 spots coloured by cell type and response.

![02b_ct_umap_celltype_response](figures/02b_ct_umap_celltype_response.png)
#### `figures/02b_ct_umap_immune_subtypes.png`

Spatial-UMAP focus on immune sub-types only.

![02b_ct_umap_immune_subtypes](figures/02b_ct_umap_immune_subtypes.png)
#### `figures/02b_ct_umap_leiden.png`

Leiden clusters on the scVI-spatial latent embedding (within-sample 20-NN purity 0.576; random ≈ 0.042).

![02b_ct_umap_leiden](figures/02b_ct_umap_leiden.png)
#### `figures/02b_ct_umap_sample.png`

Spatial-UMAP coloured by section ID — confirms cross-section integration.

![02b_ct_umap_sample](figures/02b_ct_umap_sample.png)
#### `figures/02b_ds_LR_RvsNR_comparison.png`

R vs NR comparison of top ligand–receptor pairs (CXCL9/CXCR3, CCL5/CCR5, HLA-A/CD8A).

![02b_ds_LR_RvsNR_comparison](figures/02b_ds_LR_RvsNR_comparison.png)
#### `figures/02b_ds_celltype_cooccurrence.png`

Pairwise neighbour cell-type co-occurrence in Visium spots.

![02b_ds_celltype_cooccurrence](figures/02b_ds_celltype_cooccurrence.png)
#### `figures/02b_ds_cluster_profiles.png`

Per-cluster cell-type composition profiles from the spatial atlas.

![02b_ds_cluster_profiles](figures/02b_ds_cluster_profiles.png)
#### `figures/02b_ds_comprehensive_summary.png`

Composite spatial-discovery summary: cell-type composition + IO signatures + spatial clusters per section.

![02b_ds_comprehensive_summary](figures/02b_ds_comprehensive_summary.png)
#### `figures/02b_ds_cooccurrence_RvsNR_diff.png`

Difference of pairwise neighbour co-occurrence (R − NR) across 24 sections — emphasises gain of Dysfunctional-effector CD8 ↔ Tfh-like CD4 contacts in R.

![02b_ds_cooccurrence_RvsNR_diff](figures/02b_ds_cooccurrence_RvsNR_diff.png)
#### `figures/02b_ds_gradient_spatial_maps.png`

Per-section spatial gradient maps of cytotoxicity vs hepatocyte differentiation.

![02b_ds_gradient_spatial_maps](figures/02b_ds_gradient_spatial_maps.png)
#### `figures/02b_ds_immune_gradient.png`

Tumour-edge → tumour-core immune gradient per section.

![02b_ds_immune_gradient](figures/02b_ds_immune_gradient.png)
#### `figures/02b_ds_leave_one_out.png`

Leave-one-section-out reproducibility of the spatial layer.

![02b_ds_leave_one_out](figures/02b_ds_leave_one_out.png)
#### `figures/02b_ds_ligand_receptor_heatmap.png`

LIANA ligand–receptor scoring heatmap for the Visium cohort.

![02b_ds_ligand_receptor_heatmap](figures/02b_ds_ligand_receptor_heatmap.png)
#### `figures/02b_ds_sample_level_permutation.png`

Per-section permutation null for IO signature R/NR contrasts.

![02b_ds_sample_level_permutation](figures/02b_ds_sample_level_permutation.png)
#### `figures/02b_ds_spatial_clusters.png`

Per-section spatial Leiden cluster maps.

![02b_ds_spatial_clusters](figures/02b_ds_spatial_clusters.png)
#### `figures/02b_ds_spatial_entropy.png`

Per-section spatial Shannon entropy of cell-type composition.

![02b_ds_spatial_entropy](figures/02b_ds_spatial_entropy.png)
#### `figures/02b_ds_spatially_variable_genes.png`

Top spatially variable genes (Squidpy Moran's I, 999-permutation null) per section.

![02b_ds_spatially_variable_genes](figures/02b_ds_spatially_variable_genes.png)
#### `figures/02b_pseudobulk_embeddings.png`

Spatial pseudobulk UMAP/PCA of sections, coloured by response.

![02b_pseudobulk_embeddings](figures/02b_pseudobulk_embeddings.png)
#### `figures/02b_pseudobulk_io_signatures.png`

Spatial pseudobulk IO signature heatmap across sections.

![02b_pseudobulk_io_signatures](figures/02b_pseudobulk_io_signatures.png)
#### `figures/02b_visium_RvsNR_violins.png`

Per-spot R vs NR signature violins aggregated across sections.

![02b_visium_RvsNR_violins](figures/02b_visium_RvsNR_violins.png)
#### `figures/02b_visium_celltype_spatial_maps.png`

Per-section Cell2location cell-type spatial maps.

![02b_visium_celltype_spatial_maps](figures/02b_visium_celltype_spatial_maps.png)
#### `figures/02b_visium_cross_dataset_summary.png`

Cross-dataset spatial summary across GSE173905 / GSE238264 / GSE281759.

![02b_visium_cross_dataset_summary](figures/02b_visium_cross_dataset_summary.png)
#### `figures/02b_visium_gse281759_patient_heatmap.png`

GSE281759 multi-region per-patient signature heatmap (primary vs recurrence).

![02b_visium_gse281759_patient_heatmap](figures/02b_visium_gse281759_patient_heatmap.png)
#### `figures/02b_visium_gse281759_primary_vs_recurrence.png`

Primary vs recurrence R/NR contrast on GSE281759 sections.

![02b_visium_gse281759_primary_vs_recurrence](figures/02b_visium_gse281759_primary_vs_recurrence.png)
#### `figures/02b_visium_gse281759_spatial_maps.png`

Per-section signature spatial maps on GSE281759.

![02b_visium_gse281759_spatial_maps](figures/02b_visium_gse281759_spatial_maps.png)
#### `figures/02b_visium_gse281759_tumor_vs_normal.png`

Tumour vs adjacent-normal contrast on GSE281759 sections.

![02b_visium_gse281759_tumor_vs_normal](figures/02b_visium_gse281759_tumor_vs_normal.png)
#### `figures/02b_visium_immunotype_composition.png`

Per-section immunotype composition bars.

![02b_visium_immunotype_composition](figures/02b_visium_immunotype_composition.png)
#### `figures/02b_visium_morans_I.png`

Section-level Moran's I per IO signature (responder-associated signatures show 0.12–0.41 vs 0.03–0.09 in NR).

![02b_visium_morans_I](figures/02b_visium_morans_I.png)
#### `figures/02b_visium_neighbor_enrichment_heatmap.png`

Neighbourhood-enrichment heatmap between cell-type categories.

![02b_visium_neighbor_enrichment_heatmap](figures/02b_visium_neighbor_enrichment_heatmap.png)
#### `figures/02b_visium_nhood_enrichment.png`

Squidpy nhood-enrichment z-scores.

![02b_visium_nhood_enrichment](figures/02b_visium_nhood_enrichment.png)
#### `figures/02b_visium_sample_signature_heatmap.png`

Per-sample mean signature heatmap (HCC1R–HCC7NR rows).

![02b_visium_sample_signature_heatmap](figures/02b_visium_sample_signature_heatmap.png)
#### `figures/02b_visium_signature_correlation.png`

Pairwise correlation of 16 IO signatures across spots.

![02b_visium_signature_correlation](figures/02b_visium_signature_correlation.png)
#### `figures/02b_visium_spatial_immunotype_maps.png`

Per-section immunotype spatial maps.

![02b_visium_spatial_immunotype_maps](figures/02b_visium_spatial_immunotype_maps.png)
#### `figures/02b_visium_spatial_signature_maps.png`

Per-section IO-signature spatial overlays.

![02b_visium_spatial_signature_maps](figures/02b_visium_spatial_signature_maps.png)
#### `figures/02b_visium_spatial_signature_maps_non_responder.png`

Per-section IO-signature spatial overlays.

![02b_visium_spatial_signature_maps_non_responder](figures/02b_visium_spatial_signature_maps_non_responder.png)
#### `figures/02b_visium_spatial_signature_maps_responder.png`

Per-section IO-signature spatial overlays.

![02b_visium_spatial_signature_maps_responder](figures/02b_visium_spatial_signature_maps_responder.png)
#### `figures/eda_signature_mean_heatmap.png`

Exploratory mean-signature heatmap from the spatial atlas pass.

![eda_signature_mean_heatmap](figures/eda_signature_mean_heatmap.png)
#### `figures/io_biomarker_landscape.png`

Cross-cohort IO biomarker landscape — projects atlas-derived signatures onto public IO trial transcriptomes.

![io_biomarker_landscape](figures/io_biomarker_landscape.png)
#### `figures/spatial_scvi_umap_composite.png`

Six-panel scVI spatial UMAP composite (sample, dataset, response, cell type, CNV-derived tumour status, immunotype).

![spatial_scvi_umap_composite](figures/spatial_scvi_umap_composite.png)

### `qc/`
#### `qc/GSM8627601_A009N_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627601_A009N_qc_panel](qc/GSM8627601_A009N_qc_panel.png)
#### `qc/GSM8627602_A009T1_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627602_A009T1_qc_panel](qc/GSM8627602_A009T1_qc_panel.png)
#### `qc/GSM8627603_A009T2_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627603_A009T2_qc_panel](qc/GSM8627603_A009T2_qc_panel.png)
#### `qc/GSM8627604_A009T3_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627604_A009T3_qc_panel](qc/GSM8627604_A009T3_qc_panel.png)
#### `qc/GSM8627605_A009T4_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627605_A009T4_qc_panel](qc/GSM8627605_A009T4_qc_panel.png)
#### `qc/GSM8627606_A009T5_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627606_A009T5_qc_panel](qc/GSM8627606_A009T5_qc_panel.png)
#### `qc/GSM8627607_B0013T1_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627607_B0013T1_qc_panel](qc/GSM8627607_B0013T1_qc_panel.png)
#### `qc/GSM8627608_B0013T2_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627608_B0013T2_qc_panel](qc/GSM8627608_B0013T2_qc_panel.png)
#### `qc/GSM8627609_B0013T3_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627609_B0013T3_qc_panel](qc/GSM8627609_B0013T3_qc_panel.png)
#### `qc/GSM8627610_B0013T4_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627610_B0013T4_qc_panel](qc/GSM8627610_B0013T4_qc_panel.png)
#### `qc/GSM8627611_B013T5_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627611_B013T5_qc_panel](qc/GSM8627611_B013T5_qc_panel.png)
#### `qc/GSM8627612_C004T1P_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627612_C004T1P_qc_panel](qc/GSM8627612_C004T1P_qc_panel.png)
#### `qc/GSM8627613_C004T1R_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627613_C004T1R_qc_panel](qc/GSM8627613_C004T1R_qc_panel.png)
#### `qc/GSM8627614_C008T1_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627614_C008T1_qc_panel](qc/GSM8627614_C008T1_qc_panel.png)
#### `qc/GSM8627615_C008T2_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627615_C008T2_qc_panel](qc/GSM8627615_C008T2_qc_panel.png)
#### `qc/GSM8627616_C008T3_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627616_C008T3_qc_panel](qc/GSM8627616_C008T3_qc_panel.png)
#### `qc/GSM8627617_C008T4_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627617_C008T4_qc_panel](qc/GSM8627617_C008T4_qc_panel.png)
#### `qc/HCC1R_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![HCC1R_qc_panel](qc/HCC1R_qc_panel.png)
#### `qc/HCC2R_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![HCC2R_qc_panel](qc/HCC2R_qc_panel.png)
#### `qc/HCC3R_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![HCC3R_qc_panel](qc/HCC3R_qc_panel.png)
#### `qc/HCC4R_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![HCC4R_qc_panel](qc/HCC4R_qc_panel.png)
#### `qc/HCC5NR_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![HCC5NR_qc_panel](qc/HCC5NR_qc_panel.png)
#### `qc/HCC6NR_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![HCC6NR_qc_panel](qc/HCC6NR_qc_panel.png)
#### `qc/HCC7NR_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![HCC7NR_qc_panel](qc/HCC7NR_qc_panel.png)

### `qc/plots/`
#### `qc/plots/02b_qc_overview.png`

Cross-section Visium QC overview.

![02b_qc_overview](qc/plots/02b_qc_overview.png)
#### `qc/plots/GSM8627601_A009N_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627601_A009N_qc_panel](qc/plots/GSM8627601_A009N_qc_panel.png)
#### `qc/plots/GSM8627602_A009T1_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627602_A009T1_qc_panel](qc/plots/GSM8627602_A009T1_qc_panel.png)
#### `qc/plots/GSM8627603_A009T2_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627603_A009T2_qc_panel](qc/plots/GSM8627603_A009T2_qc_panel.png)
#### `qc/plots/GSM8627604_A009T3_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627604_A009T3_qc_panel](qc/plots/GSM8627604_A009T3_qc_panel.png)
#### `qc/plots/GSM8627605_A009T4_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627605_A009T4_qc_panel](qc/plots/GSM8627605_A009T4_qc_panel.png)
#### `qc/plots/GSM8627606_A009T5_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627606_A009T5_qc_panel](qc/plots/GSM8627606_A009T5_qc_panel.png)
#### `qc/plots/GSM8627607_B0013T1_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627607_B0013T1_qc_panel](qc/plots/GSM8627607_B0013T1_qc_panel.png)
#### `qc/plots/GSM8627608_B0013T2_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627608_B0013T2_qc_panel](qc/plots/GSM8627608_B0013T2_qc_panel.png)
#### `qc/plots/GSM8627609_B0013T3_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627609_B0013T3_qc_panel](qc/plots/GSM8627609_B0013T3_qc_panel.png)
#### `qc/plots/GSM8627610_B0013T4_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627610_B0013T4_qc_panel](qc/plots/GSM8627610_B0013T4_qc_panel.png)
#### `qc/plots/GSM8627611_B013T5_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627611_B013T5_qc_panel](qc/plots/GSM8627611_B013T5_qc_panel.png)
#### `qc/plots/GSM8627612_C004T1P_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627612_C004T1P_qc_panel](qc/plots/GSM8627612_C004T1P_qc_panel.png)
#### `qc/plots/GSM8627613_C004T1R_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627613_C004T1R_qc_panel](qc/plots/GSM8627613_C004T1R_qc_panel.png)
#### `qc/plots/GSM8627614_C008T1_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627614_C008T1_qc_panel](qc/plots/GSM8627614_C008T1_qc_panel.png)
#### `qc/plots/GSM8627615_C008T2_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627615_C008T2_qc_panel](qc/plots/GSM8627615_C008T2_qc_panel.png)
#### `qc/plots/GSM8627616_C008T3_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627616_C008T3_qc_panel](qc/plots/GSM8627616_C008T3_qc_panel.png)
#### `qc/plots/GSM8627617_C008T4_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![GSM8627617_C008T4_qc_panel](qc/plots/GSM8627617_C008T4_qc_panel.png)
#### `qc/plots/HCC1R_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![HCC1R_qc_panel](qc/plots/HCC1R_qc_panel.png)
#### `qc/plots/HCC2R_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![HCC2R_qc_panel](qc/plots/HCC2R_qc_panel.png)
#### `qc/plots/HCC3R_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![HCC3R_qc_panel](qc/plots/HCC3R_qc_panel.png)
#### `qc/plots/HCC4R_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![HCC4R_qc_panel](qc/plots/HCC4R_qc_panel.png)
#### `qc/plots/HCC5NR_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![HCC5NR_qc_panel](qc/plots/HCC5NR_qc_panel.png)
#### `qc/plots/HCC6NR_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![HCC6NR_qc_panel](qc/plots/HCC6NR_qc_panel.png)
#### `qc/plots/HCC7NR_qc_panel.png`

Per-section Visium QC panel: spot counts, gene detection, mitochondrial fraction, and H&E thumbnail.

![HCC7NR_qc_panel](qc/plots/HCC7NR_qc_panel.png)
