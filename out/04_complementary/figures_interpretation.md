# Figures interpretation — `out/04_complementary`

> Companion to `docs/manuscript.md` (May 2026 revision). All R-vs-NR
> statistics, AUCs and effect sizes referenced below are computed at the
> patient level on the GSE206325 anti-PD-1 discovery cohort (n=24; 7 R,
> 17 NR) unless otherwise stated. External validation uses the Cappuyns
> et al. 2023 atezolizumab + bevacizumab cohort (n=32 evaluable; 19 R,
> 13 NR). Statistics use two-sided Mann–Whitney U with BH-FDR within
> family, Cliff's δ for effect size, and 2,000-iteration patient-level
> bootstrap 95 % CIs.

## Scope
Complementary biological axes: CD8 T-cell sub-states and pseudotime, TLS spatial hubs, ligand–receptor cross-talk, TCR clonotype dynamics, metabolic reprogramming, hallmark pathway shifts, and IMbrave150 real per-patient survival validation.

## Figure-by-figure interpretation

### `figures/`
#### `figures/04_cd8_states.png`

CD8 T-cell sub-state panel (naïve-like, memory, effector, cytotoxic, dysfunctional-effector, dysfunctional-terminal/progenitor) with per-patient R/NR enrichment.

![04_cd8_states](figures/04_cd8_states.png)
#### `figures/04_cross_dataset_signatures.png`

Cross-dataset IO-signature reproducibility across the five HCC cohorts.

![04_cross_dataset_signatures](figures/04_cross_dataset_signatures.png)
#### `figures/04_hallmark_pathways.png`

Hallmark pathway R-vs-NR shift summary — IFN-α/IFN-γ/TNF-NFκB ↑R; WNT/β-catenin / oxphos / coagulation / TGF-β ↑NR.

![04_hallmark_pathways](figures/04_hallmark_pathways.png)
#### `figures/04_immune_ratios.png`

Composite immune ratios: T/Tumor AUC=0.87, CD8/Treg AUC=0.82, DC/Mac=0.49, Effector/Exhausted=0.44.

![04_immune_ratios](figures/04_immune_ratios.png)
#### `figures/04_ligand_receptor.png`

Ligand–receptor heatmap from LIANA across the discovery cohort.

![04_ligand_receptor](figures/04_ligand_receptor.png)
#### `figures/04_lr_dotplot.png`

Top R-vs-NR ligand–receptor pairs as a dot plot (size = expression, colour = R−NR).

![04_lr_dotplot](figures/04_lr_dotplot.png)
#### `figures/04_patient_subtypes.png`

Internal patient subtyping (S1 immune-active 6R/3NR, S2 intermediate 0R/2NR, S3 immune-cold 1R/12NR; χ² p=0.0092).

![04_patient_subtypes](figures/04_patient_subtypes.png)
#### `figures/04_summary_multipanel.png`

Multi-panel complementary summary tying CD8 states, signatures, LR pairs and patient subtypes into one figure.

![04_summary_multipanel](figures/04_summary_multipanel.png)
#### `figures/04b_bulk_trial_survival_validation.png`

Bulk trial-survival validation panel — same atezo+bev arm-specific PFS stratification as Figure 6e.

![04b_bulk_trial_survival_validation](figures/04b_bulk_trial_survival_validation.png)
#### `figures/04b_cell_communication_network.png`

CellPhoneDB-style cell–cell communication network — mregDC orchestrates dysfunctional-effector CD8 interactions in R.

![04b_cell_communication_network](figures/04b_cell_communication_network.png)
#### `figures/04c_tcell_pseudotime_trajectories.png`

CD8 diffusion-pseudotime trajectories per patient; root = naïve-like cluster, rescaled to [0,1]; patient-level R-vs-NR median 0.092 vs 0.091 (δ=+0.38, p=0.166, under-powered).

![04c_tcell_pseudotime_trajectories](figures/04c_tcell_pseudotime_trajectories.png)
#### `figures/04d_ccc_lr_dotplot.png`

Cell–cell communication LR dot plot (per cell-type pair).

![04d_ccc_lr_dotplot](figures/04d_ccc_lr_dotplot.png)
#### `figures/04e_tcr_clonotype_alluvial.png`

TCR clonotype alluvial plot tracking clonal expansion across CD8 states.

![04e_tcr_clonotype_alluvial](figures/04e_tcr_clonotype_alluvial.png)
#### `figures/04f_metabolic_reprogramming.png`

Metabolic reprogramming axes (oxphos vs glycolysis vs FAO) per patient.

![04f_metabolic_reprogramming](figures/04f_metabolic_reprogramming.png)
#### `figures/04g_tls_spatial_hubs.png`

TLS spatial hub detection (12-gene panel CXCL13/CCL19/CCL21/CXCR5/CR2/CD79A/MS4A1/LAMP3/CD83/CCR7/IGHG1/IGHG2); median hub count 18 (R) vs 6 (NR), p=0.057 at n=7.

![04g_tls_spatial_hubs](figures/04g_tls_spatial_hubs.png)
#### `figures/04h_spatial_high_contrast_maps.png`

High-contrast spatial maps for the manuscript companion (Extended Data Fig. 3).

![04h_spatial_high_contrast_maps](figures/04h_spatial_high_contrast_maps.png)
#### `figures/04h_tissue_stratified_spatial.png`

Tissue-stratified spatial signature maps (tumour core / tumour edge / stromal).

![04h_tissue_stratified_spatial](figures/04h_tissue_stratified_spatial.png)
#### `figures/04i_trajectory_response.png`

Pseudotime → response link: per-patient mean pseudotime vs R/NR (under-powered direction-consistent shift).

![04i_trajectory_response](figures/04i_trajectory_response.png)
#### `figures/ExtendedData_Fig4_imbrave150_real.png`

IMbrave150 real per-patient KM curves with Cox arm × CD8-IHC interaction (HR=0.57 PFS, HR=0.36 OS) — the figure used in main-text §3.10 to anchor the dominant L3 composition signal in registration-trial data.

![ExtendedData_Fig4_imbrave150_real](figures/ExtendedData_Fig4_imbrave150_real.png)
#### `figures/spatial_1_hypoxia_gradient.png`

Spatial hypoxia gradient (HIF target panel) per section — Non-Responder sections show broader hypoxic cores.

![spatial_1_hypoxia_gradient](figures/spatial_1_hypoxia_gradient.png)
#### `figures/spatial_2_cd8_exhaustion_hubs.png`

Spatial CD8 exhaustion hubs (PDCD1/HAVCR2/TOX co-expression).

![spatial_2_cd8_exhaustion_hubs](figures/spatial_2_cd8_exhaustion_hubs.png)
#### `figures/spatial_3_macrophage_polarization.png`

Spatial M1/M2 macrophage polarisation map.

![spatial_3_macrophage_polarization](figures/spatial_3_macrophage_polarization.png)
#### `figures/spatial_4_io_response_prediction.png`

Per-spot IO-response prediction using the discovery composition model — illustrates spatial heterogeneity of predicted response within sections.

![spatial_4_io_response_prediction](figures/spatial_4_io_response_prediction.png)
