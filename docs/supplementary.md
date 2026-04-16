# Supplementary Information

> Companion to `docs/manuscript.md`.
> Supplementary Methods, Supplementary Results, Supplementary Tables S1–S14, and Supplementary Figure legends S1–S9.

---

## Supplementary Methods

**S.M.1 Harmony parameters.** θ = 2, σ = 0.1, max iterations = 20; convergence tolerance 10⁻⁴; applied to the top-50 PCA embedding of the intersecting 3,000 HVG set across all five scRNA-seq cohorts.

**S.M.2 NMF stability.** NMF was run for k ∈ {6, 8, 10, 12, 15} with ten random seeds per k. k = 10 minimized mean reconstruction MSE per added component (elbow on MSE vs k) and yielded cophenetic correlation > 0.9 across seeds. Programs were labelled by their eight top loading genes after L1-normalization of each factor.

**S.M.3 LOOCV protocol.** For each of the 24 labelled patients, the remaining 23 were used for training with internal 3-fold stratified cross-validation on the log-scale regularization grid C ∈ {10⁻³, 10⁻², 10⁻¹, 1, 10}; the held-out patient was scored. AUC-ROC was computed over the 24 held-out predicted probabilities. The balanced design retained exactly three features per layer (L1–L5) for a 15-feature total, with top-three-feature selection done inside each training fold by univariate effect size.

**S.M.4 InferCNV settings.** cutoff = 0.1, cluster_by_groups = TRUE, denoise = TRUE, HMM = TRUE, analysis_mode = subclusters; reference set = CD4-Memory ∪ NK cells (> 500 cells per patient). Per-patient summaries: `cnv_score_mean` (mean absolute CNV deviation per malignant cell), `cnv_score_sd` (SD across malignant cells, i.e. clonal heterogeneity), `cnv_clone_diversity` (Shannon entropy of InferCNV clone labels).

**S.M.5 Visium projection.** Per-spot normalization to 10,000 counts, log1p-transformation, scaling; signature scoring used `scanpy.tl.score_genes` with 50 background genes matched in expression bin. Varying spot-coverage across Visium cohorts was accounted for by section-level z-scoring prior to cross-cohort comparison.

**S.M.6 Bootstrap confidence intervals.** 2,000 iterations sampling patients with replacement; percentile CIs reported at 2.5 % and 97.5 %. Predictive AUCs are reported with bootstrap CIs in `out/03_validation/tables/bootstrap_ci.csv`.

---

## Supplementary Results

**S.R.1 Batch-effect diagnostics.** Silhouette score against *dataset* label fell from 0.41 (pre-Harmony) to 0.04 (post-Harmony) on the 50-PC embedding; against *cell-type* label it was preserved (0.27 → 0.26). Entropy of batch mixing (EBM) rose from 0.21 to 0.79 (Fig. S1).

**S.R.2 Hallmark pathway table (20 pathways).** Full effect sizes and BH-adjusted p-values in **Table S9** (`out/04_complementary/tables/hallmark_pathway_R_vs_NR.csv`). Interferon-α, Interferon-γ and TNF-via-NFκB are elevated in Responders (BH-padj < 10⁻²⁵⁰ in the pseudo-bulk analysis); WNT/β-catenin, oxidative phosphorylation, complement, TGF-β and coagulation are elevated in Non-Responders.

**S.R.3 Immune ratios.** T/Tumor ratio AUC = 0.87; CD8/Treg AUC = 0.82; DC/Macro AUC = 0.49; Effector/Exhausted AUC = 0.44. Full table **S6** (`out/04_complementary/tables/immune_ratio_aucs.csv`).

**S.R.4 Ligand–receptor expansion.** Full 24-pair scan in **Table S10** (`out/04_complementary/tables/lr_interaction_scores.csv`). CXCL9–CXCR3 was the strongest Responder-specific interaction after BH correction; CCL5–CCR5 and HLA-A–CD8A ranked second and third.

**S.R.5 Cappuyns concordance.** 3 × 2 cross-tabulation between internal subtypes and Cappuyns subsets: internal S1 / Competent = 8/9; S3 / Resistant = 11/13; Cohen’s κ = 0.55; χ² p = 0.0092 (**Table S7**, `out/02j_concordance/tables/`).

**S.R.6 Patient and cohort demographics.** The primary anti-PD-1 cohort (GSE206325) comprises 29 patients contributing 731 multi-tissue samples (302 PBMC, 158 tumour, 143 lymph node, 128 adjacent-normal; **Table S13**). Of these 29 patients, 8 were clinically classified as Responders, 20 as Non-Responders and 1 as unlabelled (**Fig. S9a**; **Table S12**); 24 patients with paired tumour tissue and unambiguous labels (7 R, 17 NR) constitute the analytical cohort used throughout predictive modelling. Sample coverage ranged 3–14 per patient (median 4). Four baseline immune-infiltration categories are represented (High, Low, Exclusion, Unannotated). Cross-cohort integration pooled five scRNA-seq studies (~970 k cells, 68 patients; **Table S14**), two Visium cohorts (15 patients), and two bulk ICB trials (IMbrave150 n = 358; CheckMate-459 n = 743).

---

## Supplementary Tables

Machine-readable CSVs are stored in `out/` and referenced below.

| ID | Contents | Path |
|----|----------|------|
| S1 | Cohort accessions and sample counts | `out/06_manuscript_figures/tables/S14_cohort_summary.csv` |
| S2 | Cross-dataset signature reproducibility | `out/04_complementary/tables/cross_dataset_signatures.csv` |
| S3 | CNV per-patient R-vs-NR statistics | `out/05_multilayer/tables/cnv_response_stats.csv` |
| S4 | NMF program top genes per program | `out/05_multilayer/tables/nmf_top_genes_per_program.csv` |
| S5 | Balanced variance partitioning | `out/05_multilayer/tables/balanced_variance_partition.csv` |
| S6 | Immune-ratio AUCs | `out/04_complementary/tables/immune_ratio_aucs.csv` |
| S7 | Cappuyns concordance statistics | `out/02j_concordance/tables/concordance_stats.csv` |
| S8 | CD8 pseudotime R-vs-NR stats | `out/04_complementary/tables/04i_pseudotime_response_stats.csv` |
| S9 | Hallmark pathway R-vs-NR effects | `out/04_complementary/tables/hallmark_pathway_R_vs_NR.csv` |
| S10 | Ligand–receptor interaction scores | `out/04_complementary/tables/lr_interaction_scores.csv` |
| S11 | Spatial section-level signature scores | `out/02k_cappuyns_spatial/tables/*.csv` |
| S12 | Patient demographics (primary cohort) | `out/06_manuscript_figures/tables/S12_patient_demographics.csv` |
| S13 | Sample overview (per-sample metadata) | `out/06_manuscript_figures/tables/S13_sample_overview.csv` |
| S14 | Cross-cohort summary | `out/06_manuscript_figures/tables/S14_cohort_summary.csv` |

---

## Supplementary Figure legends

**Fig. S1 | Atlas QC.** Per-dataset cell counts, per-cell QC distributions (n_genes, n_counts, pct_mito), batch-effect diagnostics before / after Harmony.

**Fig. S2 | Immune-lineage UMAPs.** Major and subtype-specific UMAPs of T, B, NK, DC and macrophage populations.

**Fig. S3 | Marker-based cell-type annotation.** Dotplot and heatmap of canonical marker genes across clusters; concordance with GSE206325 baseline annotations.

**Fig. S4 | Spatial diagnostics.** Visium QC, UMI-per-spot distributions, tissue masks and permutation test summaries for signature scoring.

**Fig. S5 | Cross-cohort Visium reproducibility.** Side-by-side NMF9 and TLS signature maps on GSE173905 and GSE238264 sections.

**Fig. S6 | Signature stability.** Bootstrap (2,000 iter) stability of top IO signatures; non-overlapping 95 % CIs highlighted.

**Fig. S7 | Reproducibility diagnostics.** Cross-dataset Pearson correlations of signature means; Harmony vs no-integration ablation.

**Fig. S8 | Patient subtyping details.** Hierarchical clustering dendrogram of the 24 × 16 patient–signature matrix; subtype membership and response annotation.

**Fig. S9 | Patient and cohort demographics.** (a) Patient response labels (GSE206325). (b) Sample distribution by tissue. (c) Samples collected per patient. (d) Cells per scRNA-seq cohort. (e) Patients per cohort across all modalities. (f) Immune-infiltration category × response crosstab in the labelled subset. (File: `SupplementaryFig_S9_patient_demographics.png`)

---

## Data and code availability

All scripts, tables and figures are available under MIT license at **https://github.com/isadeghi87/spatial_immunotherapy_response**. Raw datasets are public under their respective GEO / ArrayExpress accessions (see Table S1).
