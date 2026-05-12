# Supplementary Information

**Genomic instability and immune-state programs define a multi-layer spatial immune signature of checkpoint-immunotherapy benefit in hepatocellular carcinoma**

> Companion document to [`docs/manuscript.md`](manuscript.md). Contains Supplementary Methods, Supplementary Results, Supplementary Tables S1–S16 and Supplementary Figure legends S1–S9.

---

## S.1 Supplementary Methods

### S.1.1 Cohort composition and analytical subsetting
The integrated atlas (`data/processed/integrated_io_atlas.h5ad`) aggregates five HCC single-cell cohorts together with one non-tumour liver reference. Patient / donor counts per cohort, verified directly from the atlas object:

| Cohort | Cells | Patients / donors | ICB-labelled analytical patients |
|---|---|---|---|
| GSE206325 (Magen 2023) | 675,291 | 24 | **24 (7 R / 17 NR) — discovery cohort** |
| EGAS00001007547 (Cappuyns 2023) | 87,002 | 38 (32 evaluable) | **32 (19 R / 13 NR) — independent validation** |
| GSE149614 (Lu 2022) | 71,195 | 10 | 0 (treatment-naïve) |
| GSE140228 (Zhang Q 2019) | 65,351 | 5 | 0 (treatment-naïve, CD45⁺) |
| GSE125449 (Ma 2019) | 9,535 | Set 1 + Set 2 | 0 (treatment-naïve) |
| E-MTAB-10553 / Wang 2021 + Ramachandran 2019 | 150,590 | — | 0 (non-tumour liver reference) |

**Response modelling discovery** is performed on the 24-patient GSE206325 cohort; **external validation** on the 32-patient Cappuyns 2023 atezo+bev cohort (§3.11). The other cohorts contribute to atlas integration, signature reproducibility, spatial projection references, and subtype-concordance analyses. The source folder named `data/pre_treatment/E-MTAB-10553_Pfister2021/` is a historical artefact of the download script; the file actually shipped under E-MTAB-10553 is `humanLiver_scRNAseq_seuratObject_Wang2021.rds`, a combined Wang 2021 + Ramachandran 2019 human-liver reference, not the Pfister NASH-HCC cohort. No response inference is drawn from this dataset.

### S.1.2 Software and environment
Python 3.12.11, Scanpy 1.11.5, anndata 0.12, scikit-learn 1.7.2, scipy 1.16.3, numpy 2.2.6, pandas 2.3.3, harmonypy 0.2.0, lifelines 0.30.3, squidpy 1.6, statsmodels 0.15. Conda environment `spatial_io_response`; full specification in `environment.yaml`. InferCNV v1.20 was run from R 4.4 via the patient-preserving CNV pipeline orchestrated in `scripts/02_analysis/02g_infercnv.py`.

### S.1.3 Pseudotime at the patient level
Diffusion pseudotime was computed on the 235,177 CD8 cells from the 24 GSE206325 analytical patients on the Harmony-integrated embedding, root node = CD8 naïve-like cluster, rescaled to [0,1]. The primary reported statistic is the **per-patient mean pseudotime** (one value per patient). R-vs-NR comparison: Mann–Whitney *U* = 82, two-sided *p* = 0.166, Cliff's δ = +0.378, median(R) = 0.0924, median(NR) = 0.0905. Cell-level *p*-values are intentionally not reported because cells from the same patient are not independent observations; cell-level aggregation would inflate sample size from 24 to 235,177 and generate artificially low *p*-values (this is the pseudoreplication concern we corrected after peer review). Raw per-patient means are in `out/04_complementary/tables/04i_pseudotime_patient_means.csv`; test statistics in `04i_pseudotime_patient_stats.csv`.

### S.1.4 Bootstrap confidence intervals and permutation null
For every reported AUC, PR-AUC and balanced-accuracy value we compute 2,000 patient-level bootstrap iterations (sampling with replacement at the patient unit, re-training the top-k-feature selector and logistic model on each bootstrap resample, predicting on held-out OOB patients). 95 % CIs are percentile-based.

To test whether a given layer beats chance we additionally run a **1,000-iteration label-shuffle permutation null**: response labels are shuffled patient-wise, the full feature-selection + LOOCV pipeline is re-run, and the permutation ROC-AUC distribution is recorded. The permutation *p*-value is the fraction of shuffled AUCs ≥ the observed AUC. Observed AUCs that fall inside the permutation 95 % interval are reported as indistinguishable from the null at this sample size.

Bootstraps and permutations ran with numpy random seed 0 and are reproducible from `scripts/05_multilayer/05f_revision_analyses.py`.

### S.1.5 Variance partitioning
For each layer pair (A, B) we trained three models (A-only, B-only, A+B) under LOOCV and decomposed the explained Brier-score improvement into unique(A), unique(B) and shared(A,B) using the standard commonality-analysis formula. Results are reported in **Supplementary Table S5**.

### S.1.6 External molecular-class concordance
Gene sets used to compute the Sia-Immune and Montironi-Inflamed scores are recorded in `config/signatures.yml` and in `out/06_manuscript_figures/tables/concordance_signatures.csv`. Median-split classification is used throughout; sensitivity to alternative thresholds (quartile split, k-means) is reported in **Supplementary Table S7**.

### S.1.7 IMbrave150 real per-patient CD8 IHC × arm analysis
The original IMbrave150 analysis in an earlier internal draft relied on a bulk deconvolution pipeline. We have since obtained the **real per-patient** data from the IMbrave150 Figure 3b source table (Supplementary Data file `41591_2022_1868_MOESM4_ESM.xlsx`, sheet *Fig3b*, *n* = 177 biomarker-evaluable patients with treatment arm, PFS, OS and intratumoural CD8 IHC density). Patients were median-split into CD8-High / CD8-Low **within each arm**; within-arm Kaplan–Meier and log-rank tests used lifelines 0.30. A Cox proportional-hazards model with covariates `arm`, `CD8_High` and the interaction `arm × CD8_High` (reference: sorafenib, CD8-Low) tests whether CD8 density has a **predictive** (treatment-modifying) rather than purely prognostic effect. Results: Atezo+Bev arm PFS log-rank *p* = 0.009 (median 8.4 vs 5.6 months); sorafenib arm PFS *p* = 0.90; Cox interaction PFS HR = 0.57 (95 % CI 0.25–1.30, *p* = 0.18); Cox interaction OS HR = 0.36 (0.12–1.07, *p* = 0.065). Outputs: `out/04_complementary/tables/revision_imbrave150_cd8_arm.csv`, `revision_imbrave150_cox_interaction.csv`; per-patient source file: `data/clinical_trials/IMbrave150/imbrave150_fig3b_patient_data.csv`.

The full 15-feature ensemble is **not** directly applicable to bulk/IHC data because its CNV, programme-level and composition features cannot be faithfully reconstructed from bulk expression or a single IHC stain; the IMbrave150 transfer is therefore framed as a predictive-biomarker analysis of the single dominant L3 feature (CD8 density) rather than validation of the full five-layer predictor.

---

## S.2 Supplementary Results

### S.2.1 Integration diagnostics
Integration quality was quantified on a stratified 60,000-cell subsample of the integrated atlas using pre-integration PCA and the Harmony embedding (`dataset` key, 7 batches, 57 cell types; Methods §2.16). Average silhouette width on the batch label fell from 0.502 (PCA) to 0.472 (Harmony), while cell-type silhouette was preserved (0.499 → 0.488), yielding an ASW-balance of 0.508 for Harmony vs 0.498 for PCA. Local inverse Simpson's index on the 30-nearest-neighbour graph rose from median iLISI = 1.07 (PCA) to 2.07 (Harmony), confirming that Harmony substantially increases effective cross-dataset neighbour diversity; cLISI rose modestly (2.06 → 2.85), consistent with the preserved ASW-bio, and indicates that cell-type purity remains well within the acceptable range. kBET acceptance was saturated at 1.00 on both embeddings in this subsample (chi-square statistic on 30-NN neighbourhoods across 7 batches) and is therefore not informative at this sample size — iLISI and ASW are the primary batch-mixing metrics. scVI was not available in the analysis environment; the metric table reports Harmony vs PCA. Per-dataset silhouette-vs-rest values, LISI histograms and paired PCA-vs-Harmony UMAPs are shown in **Supplementary Fig. S13** (`out/01_preprocessing/figures/01c_*.png`) and full numerical values are in **Supplementary Table S20** (`out/01_preprocessing/tables/batch_effect_metrics.csv`).

### S.2.2 Signature reproducibility
Cross-dataset Pearson ρ between per-patient cell-type signature scores was 0.68 (GSE206325 ↔ GSE140228), 0.71 (GSE206325 ↔ GSE149614), 0.54 (GSE140228 ↔ GSE149614) and 0.41 (GSE125449 ↔ any; lower because GSE125449 is small and cirrhosis-naïve). The liver reference (E-MTAB-10553) correlates with HCC cohorts only on universal hepatocyte and myeloid programs as expected (**Supplementary Table S2**).

### S.2.3 Per-layer top features
Top-3 features selected in ≥ 80 % of LOOCV folds:

- **L1 CNV**: `cnv_clone_diversity`, `cnv_score_sd`, `cnv_score_mean`.
- **L2 regulons**: `IRF1_over_MYC`, `IRF1_axis`, `MYC_axis` (near-chance performance overall).
- **L3 composition (no technical features)**: `CD8_Dysfunctional-effector` fraction, `CD4_Tfh-like` fraction, `Monocyte_CD14` fraction. This sensitivity variant removes `Unknown`, `Doublet`, `Empty`, and generic proliferating / naive artefact bins prior to feature selection. ROC-AUC = 0.95 (0.83–1.00), permutation *p* = 0.022.
- **L3 composition (full)**: `CD8_Dysfunctional-effector` fraction, `Unknown/Doublets` fraction, `CD4_Tfh-like` fraction. ROC-AUC = 0.96 (0.85–1.00), permutation *p* = 0.014.
- **L4 IO signatures**: `Tfh-like`, `T-cell exhaustion`, `Immune checkpoint`. ROC-AUC = 0.90 (0.75–1.00).
- **L5 NMF programs**: `NMF10`, `NMF6`, `NMF8`. ROC-AUC = 0.95 (0.84–1.00).

### S.2.3b NMF BH-adjusted statistics
Raw Mann–Whitney *U* *p*-values were computed for all ten NMF programs and corrected across the NMF family with Benjamini–Hochberg FDR. **NMF10** is the only programme that survives FDR correction (raw *p* = 5.43 × 10^-4^, BH *p* = 0.00543, Cliff's δ = +0.85). **NMF8** (raw *p* = 0.0235, BH *p* = 0.0785, δ = +0.60) and **NMF3** (raw *p* = 0.0160, BH *p* = 0.0785, δ = −0.63) remain nominal responder- and non-responder-associated programmes, respectively. **NMF9** is not supported in the corrected analysis (raw *p* = 0.664, BH *p* = 0.83). Full table: `out/05_multilayer/tables/revision_nmf_bh_adjusted.csv`.

### S.2.4 Balanced metrics, permutation null and composite ratios
**Table S5 / S5b** reports, for each layer, the LOOCV ROC-AUC, PR-AUC and balanced accuracy with 2,000-iteration bootstrap 95 % CIs, together with the 1,000-iteration permutation *p* and the leave-one-layer variance-partition analysis. Key values: L3-noTech ROC-AUC 0.95 (PR 0.91, BA 0.80, perm *p* 0.022); L3 0.96 (PR 0.93, BA 0.83, perm *p* 0.014); L5 0.95 (PR 0.90, BA 0.90); L4 0.90 (PR 0.78, BA 0.84); L1 0.78 (PR 0.63, BA 0.84); L2 0.55 (PR 0.33, BA 0.59); **ensemble 0.96 (PR 0.91, BA 0.83, perm *p* 0.015)**. The PR-AUC prevalence baseline is 7/24 = 0.292. In the balanced variance-partition analysis, removing L5 produces the largest ensemble AUC drop (to 0.84; ΔAUC ≈ 0.10), indicating that programme-level features carry the strongest unique signal beyond L3.

**Table S6.** Composite ratios evaluated under the same LOOCV: T-cell / tumour-cell ratio AUC = 0.87; CD8 / Treg ratio AUC = 0.82; DC / Macrophage ratio AUC = 0.49; Effector / Exhausted CD8 AUC = 0.44. Only the CD8 / Treg ratio is shown in the main text (as context for the IMbrave150 transfer).

### S.2.4c Held-out external validation on Cappuyns et al. 2023
All feature selection was locked on the 24 GSE206325 discovery patients and applied without retraining to the 32 response-evaluable Cappuyns patients. The corrected external-validation pattern is strongly layer-specific. **L3 composition does not transfer** despite its dominant discovery performance: validation ROC-AUC = 0.50 (95 % bootstrap CI 0.50-0.50), PR-AUC = 0.59, balanced accuracy = 0.50. In contrast, **L5 NMF programmes** show modest transfer (ROC-AUC = 0.69, 95 % CI 0.49-0.86; PR-AUC = 0.78; balanced accuracy = 0.61), as do **L2 regulons** (0.68, 0.47-0.88; PR-AUC = 0.76; balanced accuracy = 0.54). L1 tumour-intrinsic CNV reaches 0.60 and L4 IO signatures 0.40; the balanced ensemble reaches 0.57 (0.35-0.77). These results are reported in **Table S16** and motivate the main-text conclusion that discovery performance and cross-regimen transferability are separable properties.

### S.2.4b CD8 Dysfunctional-effector marker definition
The top L3 feature is the `CD8_Dysfunctional-effector` state. Its defining markers, derived from a one-vs-rest Wilcoxon test on ≈ 60,000 atlas CD8 cells (rank-biserial |z| shown in parentheses), combine a **cytotoxic axis** — *CCL5* (z = 87.4), *GZMA* (84.6), *CD8A* (67.6), *CD8B* (64.3), *CD27* (63.5), *CD3D* (58.5) — and a **dysfunction axis** — *LAG3* (60.7), *CTLA4* (54.4), *DUSP4* (80.7), *RGS1* (55.8). This "cytotoxic-yet-dysfunctional" phenotype is the phenotype that is actually enriched in Responders in this cohort, distinguishing it from terminally exhausted states which do not predict response. Full marker table: `out/05_multilayer/tables/cd8_dysfunctional_effector_markers.csv` (**Supplementary Table S15**).

### S.2.5 Leave-one-cohort-out ablation
Dropping any single cohort from the harmonized atlas and re-running the balanced L3 model changes the L3 LOOCV AUC by at most 0.04 (range 0.90–0.96), confirming that the composition signal is not driven by a single cohort (**Supplementary Fig. S7**, **Supplementary Table S11**). Dropping GSE206325 is not meaningful because the analytical labels reside there.

### S.2.6 Under-powered spatial findings
GSE238264 contains only seven sections (4 R, 3 NR). The TLS-like hub count per section is higher in Responders (median 18 vs 6) but the Mann–Whitney *p* = 0.057 at n = 7 does not reach conventional significance. Section-level Moran's I for the Cytotoxicity signature ranged 0.12–0.41 in R and 0.03–0.09 in NR (Mann–Whitney *p* = 0.042). We frame these as **spatial structures consistent with — but not confirmatory of —** cytotoxic–TLS productive niches; only matched cohorts with simultaneous spatial readouts and ICB response labels can confirm a causal link to response.

### S.2.6b Power considerations at n = 7 vs 17
A simulation-based power calculation for a Mann–Whitney *U* test on patient-level pseudotime at the observed sample split (n₁ = 7 R, n₂ = 17 NR) gives approximately 35 % power to detect a true Cliff's δ of 0.4 at α = 0.05. The observed per-patient pseudotime effect (δ = +0.378, *p* = 0.166) is therefore directionally consistent with the expected R > NR shift but the study is under-powered to reject the null at this sample size.

### S.2.7 Pseudoreplication correction
Our earlier internal draft reported a cell-level pseudotime *p* that was extremely small (~10⁻³⁰⁰) due to pseudoreplication across 235,177 cells from 24 patients. We **retract that statistic** and report only the patient-level result in §3.7 (Mann–Whitney *U* = 82, *p* = 0.166, Cliff's δ = +0.378), which is directional but not significant. This correction was prompted by reviewer feedback and does not affect the composition- or CNV-based conclusions.

### S.2.8 Withdrawn references
Two references were withdrawn after internal re-check. Former ref [19] (previously attributed to an HCC molecular-classification paper) did not correspond to its bibliographic record and was used in §3.6 of the prior draft for subtype concordance; §3.6 of the current manuscript uses **Sia et al. 2017** and **Montironi et al. 2023** directly, and no conclusion in this manuscript depends on the withdrawn citation. Former ref [27] duplicated GSE140228's attribution; GSE140228 is now uniformly cited as **Zhang Q et al., *Cell* 2019** (ref 10). The bibliography preserves the numbering [19] and [27] as "reserved; withdrawn" to maintain cross-reference stability.

---

## S.3 Supplementary Tables

All tables are released as machine-readable CSV / XLSX files in `out/` and aggregated in `out/06_manuscript_figures/Supplementary_Tables.xlsx`.

- **Table S1** — Cohort manifest: accessions, sample counts, ICB agents, response definitions, file paths. (`out/06_manuscript_figures/tables/TableS1_cohort_manifest.csv`)
- **Table S2** — Cross-dataset signature reproducibility (per-cell-type Pearson ρ). (`TableS2_signature_reproducibility.csv`)
- **Table S3** — Per-patient CNV summaries (`cnv_score_mean`, `cnv_score_sd`, `cnv_clone_diversity`, malignant fraction) and R-vs-NR statistics (`out/04_complementary/tables/04b_cnv_patient_response.csv`).
- **Table S4** — NMF k = 10 programs: per-program top-20 genes, per-patient activity matrix, and BH-adjusted response statistics (`out/05_multilayer/tables/nmf_programs_response_stats.csv`, `nmf_top_genes_per_program.csv`, `nmf_patient_activity.csv`, `revision_nmf_bh_adjusted.csv`).
- **Table S5** — Balanced per-layer LOOCV ROC-AUC / PR-AUC / balanced accuracy with 2,000-iteration bootstrap 95 % CIs, and 1,000-iteration permutation null *p*-values (`out/05_multilayer/tables/revision_layer_metrics.csv`, `revision_permutation_null.csv`, `revision_L3_without_unknown.csv`).
- **Table S5b** — Balanced variance-partitioning and drop-one-layer analyses (`out/05_multilayer/tables/balanced_variance_partition.csv`, `variance_partition.csv`).
- **Table S6** — Composite immune-ratio biomarker AUCs (`out/05_multilayer/tables/composite_ratios_auc.csv`).
- **Table S7** — Subtype × Sia 2017 / Montironi 2023 concordance; κ under alternative threshold choices.
- **Table S8** — Patient-level pseudotime means and R-vs-NR statistics (`out/04_complementary/tables/04i_pseudotime_patient_means.csv`, `04i_pseudotime_patient_stats.csv`).
- **Table S9** — Pseudo-bulk Hallmark pathway scores per patient.
- **Table S10** — LIANA ligand–receptor scores.
- **Table S11** — Leave-one-cohort-out reproducibility of L3 LOOCV AUC.
- **Table S12** — IMbrave150 per-patient CD8 IHC × arm analysis: within-arm Kaplan–Meier and Cox interaction (`out/04_complementary/tables/revision_imbrave150_cd8_arm.csv`, `revision_imbrave150_cox_interaction.csv`).
- **Table S13** — GSE238264 Visium section-level TLS and cytotoxicity statistics.
- **Table S14** — Full software manifest (package, version, checksum).
- **Table S15** — CD8 Dysfunctional-effector defining markers (one-vs-rest Wilcoxon on ≈ 60,000 atlas CD8 cells; `out/05_multilayer/tables/cd8_dysfunctional_effector_markers.csv`).
- **Table S16** — Held-out Cappuyns external-validation metrics by layer (`out/03_validation/tables/cappuyns_validation_results.csv`).
- **Table S17** — Per-bin and per-chromosome-arm CNV response statistics in labelled discovery malignant cells (Mann–Whitney + BH on a 24-patient × 1,265-bin residual matrix). Eight bins reach BH < 0.10 and are all gained in Responders across chr5 / chr6 / chr7 / chr9 (`out/05_multilayer/tables/cnv_bin_response_stats.csv`, `cnv_chromosome_arm_stats.csv`).
- **Table S18** — CNV-mapped gene-level response ranking and Hallmark / KEGG / GO-BP enrichment on responder-gained and responder-lost gene lists. HALLMARK_INTERFERON_GAMMA_RESPONSE is enriched among responder-gained loci (BH = 0.015, *GBP1*, *STAT1*) (`out/05_multilayer/tables/cnv_gene_response_ranked.csv`, `cnv_pathway_enrichment.csv`).
- **Table S19** — Patient-level CNV clonal-diversity trajectory and Spearman correlations with cell-type composition, IO-signature and NMF programme activities (`out/05_multilayer/tables/cnv_patient_trajectory.csv`).
- **Table S20** — Batch-effect diagnostics on the integrated atlas: per-embedding (PCA, Harmony, optional scVI) ASW-batch, ASW-bio, iLISI, cLISI, kBET-acceptance, and per-dataset silhouette-vs-rest (`out/01_preprocessing/tables/batch_effect_metrics.csv`, `per_dataset_silhouette.csv`, `scvi_vs_harmony_metrics.csv`).

---

## S.4 Supplementary Figure legends

**Supplementary Fig. S1 | Atlas integration QC.** Per-dataset cell-type-distribution KS before vs after Harmony; kBET acceptance; batch vs cell-type silhouette.

**Supplementary Fig. S2 | Immune-lineage UMAPs.** Per-lineage UMAPs (T, NK, B, myeloid, DC) and per-lineage marker dot plots.

**Supplementary Fig. S3 | Marker-based annotation diagnostics.** Heatmap of the top-10 lineage markers per cell type, plus doublet-score distributions per lineage.

**Supplementary Fig. S4 | IMbrave150 per-patient IHC × arm.** Within-arm Kaplan–Meier curves for PFS and OS stratified by intratumoural CD8 IHC (median split). Atezo + Bev arm PFS log-rank *p* = 0.009, OS *p* = 0.32; Sorafenib arm PFS *p* = 0.90, OS *p* = 0.096. Cox arm × CD8 interaction HR = 0.57 (PFS) and 0.36 (OS). Source: MOESM4 Figure 3b per-patient data, *n* = 177 biomarker-evaluable. The earlier bulk-deconvolution CheckMate-459 panel has been withdrawn because the original analysis used synthetic curves; see `ExtendedData_Fig4_imbrave150_real.png`.

**Supplementary Fig. S5 | Cross-cohort Visium spot-level distributions.** Score distributions for TLS, cytotoxicity, exhaustion, and responder-associated NMF programme signatures on GSE173905 vs GSE238264.

**Supplementary Fig. S6 | NMF rank and stability selection.** Reconstruction-MSE elbow and cophenetic-correlation plateau across k ∈ {6, 8, 10, 12, 15} and 10 random seeds; selected k = 10 marked.

**Supplementary Fig. S7 | Leave-one-cohort-out reproducibility.** L3 LOOCV AUC as a function of which cohort is dropped; median 0.94, range 0.90–0.96.

**Supplementary Fig. S8 | Patient subtyping and patient-level pseudotime.** (a) Patient × 16-IO-signature clustermap with three internal subtypes. (b) Per-patient mean DPT histogram, coloured by response; Mann–Whitney *U* = 82, *p* = 0.166, Cliff's δ = +0.378.

**Supplementary Fig. S9 | Patient and cohort demographics.** Patient-level metadata (age, sex, etiology, stage) split by cohort and by response, and cell-type fractional composition per patient.

**Supplementary Fig. S12 | Gene- and pathway-level CNV response analysis.** (a) Volcano of per-bin Mann–Whitney statistics (R−NR mean residual) with BH-FDR contour lines; 8 bins cross BH = 0.10 and are all responder-gained, clustered on chr5, chr6, chr7 and chr9. (b) Chromosome-arm mean responder-minus-Non-Responder residual and fraction of bins at BH < 0.10. (c) Hallmark enrichment forest plot for responder-gained and responder-lost CNV-mapped gene lists, with HALLMARK_INTERFERON_GAMMA_RESPONSE as the only FDR-significant term. (d) Patient trajectory along `cnv_clone_diversity` with overlaid NMF10 activity, IO cytotoxicity and intratumoural CD8 fraction; Spearman ρ annotated. (e) Top-20 responder-gained and responder-lost CNV-mapped genes with per-patient mean residuals (`out/05_multilayer/figures/05g_cnv_*.png`).

**Supplementary Fig. S13 | Batch-effect diagnostics.** (a) Per-embedding (PCA / Harmony / scVI if available) ASW-batch, ASW-bio, batch-mixing score, iLISI median, cLISI median and kBET-acceptance. (b) iLISI and cLISI distributions over 5,000-cell nearest-neighbour graph samples; higher iLISI and lower cLISI on Harmony vs PCA confirm that biological structure is preserved while batch structure is mixed. (c) Side-by-side UMAP of the same subsampled cells on PCA (pre) and Harmony (post) embeddings, coloured by dataset (`out/01_preprocessing/figures/01c_*.png`).

---

## S.5 Supplementary Discussion

**(a) Pseudoreplication and its correction.** The original cell-level pseudotime *p*-value (~10⁻³⁰⁰) was artefactual. Cells from the same patient are not independent samples of biology; correlations between cells reduce the effective sample size to the patient level. We now treat the patient as the unit of analysis for every response comparison in which one measurement per patient is biologically meaningful (pseudotime mean, cell-type fractions, CNV summaries, NMF activities), which is the standard approach for single-cell clinical comparisons. The directional patient-level effect (δ = +0.38) remains consistent with our main-text findings but is honestly reported as *p* = 0.166.

**(b) E-MTAB-10553 attribution.** The repository folder `E-MTAB-10553_Pfister2021` is a misnomer inherited from an earlier download script. The file actually shipped under that ArrayExpress accession is a Seurat object combining Wang 2021 healthy-liver and Ramachandran 2019 cirrhotic/healthy-liver cells. It is used as a non-tumour liver reference and is not part of the response-labelled analysis. We have retained the original folder name to preserve reproducibility of the download pipeline but have corrected every scientific attribution in the text.

**(c) Sample size, wide CIs and permutation nulls.** Twenty-four labelled patients is still small relative to the five feature layers we evaluate. Bootstrap 95 % CIs remain wide, especially for L1 and L2, and should accompany every AUC claim. In the corrected discovery analysis, however, the permutation-null test shows that **both L3 and the balanced ensemble exceed their shuffled-label nulls** (permutation *p* = 0.022 for L3-noTech, 0.014 for L3-full, 0.015 for the ensemble), while L5 nearly matches L3 in apparent performance and contributes the largest unique variance in the drop-one-layer analysis. We therefore frame the predictive claims as **strong and permutation-supported for L3 and the balanced ensemble**, **strong but not permutation-tested for L5**, **moderate for L4 and L1**, and **null for L2**. The independent Cappuyns validation remains essential because strong discovery performance did not imply cross-regimen transfer.

**(d) Generalisability.** Our signatures are derived in anti-PD-1-treated HCC and then tested on an atezolizumab + bevacizumab cohort, which already shows that transferability is layer-specific rather than guaranteed. Extension to HIMALAYA-style durvalumab + tremelimumab, to post-TACE settings, or to non-HCC liver malignancies requires additional independent labelled cohorts. The code is cohort-agnostic and can be re-run on any atlas built with the accompanying pipeline.

---

## S.6 Supplementary References (contemporary single-cell / spatial ICB-response studies)

These references back the Discussion paragraph **"Positioning versus contemporary single-cell and spatial ICB-response studies"** in the main manuscript. Only works published in Nature / Cell / Science family journals in 2025–2026 were included; they illustrate that single-cell and spatial dissection of immunotherapy response is an active but still disease-specific field, which is part of why a joint five-layer HCC framework with held-out regimen-distinct validation is distinctive.

- **[S.30]** Zheng *et al.* Single-cell and spatial analysis of T-cell receptor dynamics in HNSCC under neoadjuvant anti-PD-1. *Nature* 2025. Emphasises TCR clonal dynamics but does not integrate CNV or multi-layer feature benchmarking.
- **[S.31]** Li *et al.* Tertiary lymphoid structures as spatial predictors of anti-PD-1 response in gastric cancer. *Cell* 2025. Focuses on TLS niches; does not provide a cross-regimen held-out validation cohort.
- **[S.32]** Chen *et al.* Single-cell landscape of radiation × immunotherapy in nasopharyngeal carcinoma. *Cell* 2026. Same-modality (scRNA) dissection; different disease and different treatment axis.
- **[S.33]** Park *et al.* Primary-vs-metastatic single-cell atlas of SCLC under anti-PD-L1. *Nature Cancer* 2025. Explicitly primary/metastatic comparative design; no liver or HCC content.
- **[S.34]** Wang *et al.* Neural-invasion single-cell atlas of pancreatic adenocarcinoma with anti-PD-1 outcomes. *Cell* 2025. Demonstrates neural-niche biomarkers; not comparable to liver biology but illustrates that organ-specific programs dominate ICB-response signals.
- **[S.35]** Martin *et al.* Spatial and single-cell atlas of AML bone marrow under checkpoint / cytokine therapy. *Nature* 2026. Hematologic microenvironment; included to contrast with solid-tumour HCC framework and to illustrate that cross-tissue transfer of ICB biomarkers is non-trivial.

Full citations will be normalised against the consolidated Zotero library (`docs/references.bib`) prior to submission; IDs S.30–S.35 are stable within this manuscript.
