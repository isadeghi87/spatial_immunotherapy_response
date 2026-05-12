# Materials and Methods

> Companion document to [`docs/manuscript.md`](manuscript.md). Mirrors §2 of the manuscript and expands on parameter choices.

## 1. Study cohorts and data acquisition
Public single-cell RNA-seq and paired Visium datasets were downloaded from GEO, ArrayExpress and the Lambrechts VIB portal. The integrated atlas aggregates five HCC scRNA-seq cohorts and one non-tumour liver reference:

- **GSE206325** — Magen et al., *Nat Med* 2023; anti-PD-1 / combination ICB HCC single-cell cohort. The primary analytical subset consists of the **24 patients** with paired tumour tissue and unambiguous RECIST labels (7 Responders, 17 Non-Responders), identified by the atlas column `patient_ID` × `response`. This cohort serves as the **discovery set** for all response modelling.
- **EGAS00001007547 (Cappuyns 2023)** — Cappuyns et al., *Nat Commun* 2023; atezolizumab + bevacizumab HCC single-cell cohort (≈98,000 cells, 38 patients). Of the 38 patients, 32 carry evaluable response labels (19 Responders, 13 Non-Responders; 4 excluded: DeathBeforeImaging; 2: NA). Week-0 (pre-treatment) tumour biopsies only. This cohort serves as the **independent external validation set**.
- **GSE125449** — Ma et al., *Cancer Cell* 2019; treatment-naïve HCC landscape (≈9,500 cells across Set 1 and Set 2).
- **GSE140228** — Zhang Q et al., *Cell* 2019; CD45⁺ immune-enriched HCC (≈65,000 cells, five donors). All GSE140228 attribution in this manuscript is to Zhang Q et al., *Cell* 2019.
- **GSE149614** — Lu et al., *Nat Commun* 2022; multi-region HCC ecosystem (≈71,000 cells, ten patients).
- **E-MTAB-10553** — Wang et al. 2021 and Ramachandran et al. 2019 human liver single-cell reference. The file shipped under this accession is a Seurat object (`humanLiver_scRNAseq_seuratObject_Wang2021.rds`) that combines Wang 2021 healthy-liver cells with Ramachandran 2019 healthy and cirrhotic liver cells. We use it **exclusively as a non-tumour liver reference** and never for response inference. The name of the source directory in our repository (`E-MTAB-10553_Pfister2021`) is historical and does not correspond to the Pfister 2021 NASH-HCC cohort.

Spatial cohorts (projection only): **GSE173905** (Wu et al., *Sci Adv* 2021), **GSE238264** (Liu et al., *Cell Rep Med* 2023; 4 R + 3 NR sections). Trial-adjacent external transfer used **IMbrave150** per-patient intratumoural CD8 IHC data (Finn et al., *NEJM* 2020, supplementary matrix); **CheckMate-459** was not analysed because no equivalent per-patient table was publicly available.

## 2. Single-cell preprocessing and integration
Processed with Scanpy v1.11.5 (Wolf et al. 2018). QC: `200 ≤ n_genes ≤ 8,000`, `pct_mt < 15 %`. Normalized to 10,000 counts per cell, log1p, 3,000 highly variable genes (Seurat v3 flavour), PCA (50 components). Batch correction with Harmony (Korsunsky et al. 2019; θ = 2, σ = 0.1, 20 iterations). k-NN (k = 15) on the harmonized embedding, UMAP, Leiden clustering (resolution = 0.8). Final atlas stored as `data/processed/integrated_io_atlas.h5ad` with **1,058,964 cells × 15,697 genes**.

## 3. Cell-type annotation and IO signatures
Cell types were annotated by marker-gene enrichment cross-referenced against GSE206325 baseline labels. The canonical column `cell_type_original` resolves 44 sub-lineages including six CD8 T-cell states (naïve-like, memory, effector, cytotoxic, dysfunctional-effector, dysfunctional-terminal / progenitor), CD4 helpers (memory, Tfh-like, Th1-GZMK), Tregs, NK, mregDC, conventional DC, macrophage M1 / M2, B cells, plasma cells, endothelium, fibroblasts and malignant hepatocytes. Sixteen curated IO signatures (cytotoxicity, exhaustion, terminal exhaustion, progenitor CD8, activation, Treg, M1/M2 macrophage, IFN-γ / IFN-α response, antigen presentation, checkpoints, mregDC, Tfh, TLS, WNT/β-catenin, TGF-β, angiogenesis) were scored per cell with `scanpy.tl.score_genes` (50 mean-expression-bin matched background genes).

## 4. Differential abundance and signature testing
Per-patient fractions and mean signature scores were computed on tumour-compartment cells. R-vs-NR comparisons used two-sided Mann–Whitney U tests (scipy 1.16) with Benjamini–Hochberg FDR correction within each analytical family. Effect sizes are reported as **Cliff's δ** throughout (rank-based, scale-free, sign convention: positive δ = higher in Responders).

## 5. Copy-number variation (InferCNV)
InferCNV was run on a **50,000-cell patient-preserving subset** of the integrated atlas (`data/processed/integrated_io_atlas_cnv.h5ad`) because the full million-cell atlas is computationally intractable for expression-based CNV calling. All non-tumour immune / stromal lineages were used as the reference set and hepatocyte / tumour-like compartments as the test set. Cutoff 0.1, HMM enabled, subclusters mode, denoise = TRUE. Three patient-level CNV summaries were retained: **`cnv_score_mean`** (mean absolute inferred deviation per malignant cell; chromosomal-burden proxy), **`cnv_score_sd`** (standard deviation across malignant cells) and **`cnv_clone_diversity`** (Shannon entropy of InferCNV clone labels).

## 6. NMF expression programs
`sklearn.decomposition.NMF` (scikit-learn 1.7) was fit on a **20,000-cell discovery-only subset** drawn from the labelled GSE206325 cells after HVG selection. Parameters: `init='nndsvda'`, `max_iter=400`, `tol=1e-3`. `k ∈ {5, 7, 10, 12, 15, 20}`; **k = 10** was chosen by the reconstruction-error elbow. The fitted discovery dictionary was then projected onto all labelled discovery and validation cells. Programs were labelled by top-8 loading genes, and per-patient programme activity was defined as the mean projected score across all labelled cells from that patient.

## 7. Regulons (L2)
AUCell-style scoring (decoupler, Aibar et al. 2017) on 12 master regulators: IRF1, STAT1, TCF7, TOX, EOMES, FOXP3, MYC, SPI1, IRF8, NR4A1, BATF, PRDM1. Per-patient regulon activity vectors are the L2 feature layer. L2 is included for transparency; in the 24-patient cohort it yields near-chance discrimination.

## 8. Balanced five-layer predictive framework
Exactly **three features per layer**, selected on the training fold by univariate |t| statistic against response (L1 CNV, L2 regulons, L3 composition, L4 IO signatures, L5 NMF). Classifier: `sklearn.linear_model.LogisticRegression` with L2 penalty (`C = 1.0` per layer, `C = 0.5` for the 15-feature ensemble), `class_weight='balanced'`, `max_iter=3000`. Evaluation: **leave-one-patient-out cross-validation** on the 24-patient cohort.

For each layer (and the ensemble) we report three metrics: **ROC-AUC**, **precision–recall AUC** (prevalence baseline 7 / 24 = 0.292) and **balanced accuracy**. All three are reported with **2,000-iteration patient-level bootstrap 95 % confidence intervals**. To test whether a layer beats chance we additionally run a **1,000-iteration label-shuffle permutation null**: response labels are shuffled patient-wise, the full feature-selection + LOOCV pipeline is re-run, and the permutation *p* is the fraction of shuffled AUCs ≥ the observed AUC.

For L3 (cell-type composition) we additionally report an **L3-noTech** sensitivity variant in which `Unknown`, `Doublet`, `Empty`, and generic proliferating / naive artefact bins are removed prior to feature selection. L3-noTech ROC-AUC = 0.95 (0.83–1.00), permutation *p* = 0.022; L3 (full) ROC-AUC = 0.96 (0.85–1.00), permutation *p* = 0.014; L5 ROC-AUC = 0.95 (0.84–1.00); ensemble ROC-AUC = 0.96 (0.87–1.00), **permutation *p* = 0.015**. Leave-one-layer analysis further showed that removing L5 produced the largest ensemble AUC drop (to 0.84). All results are reproducible from `scripts/05_multilayer/05f_revision_analyses.py`.

## 9. Spatial projection and TLS scoring
Each Visium section was processed identically to the scRNA-seq atlas (normalization, log1p, HVG selection). Sixteen IO signatures and ten NMF programs were projected per spot with `scanpy.tl.score_genes`. TLS score used a 12-gene panel from Cabrita 2020 / Helmink 2020 (*CXCL13, CCL19, CCL21, CXCR5, CR2, CD79A, MS4A1, LAMP3, CD83, CCR7, IGHG1, IGHG2*). TLS-like hubs = spots with joint TLS and cytotoxicity z-scores > 1. Section-level comparisons used Mann–Whitney U; spatial autocorrelation used Moran's I with 999 permutations. GSE238264 contains 7 sections (4 R, 3 NR), so hub-count comparisons are explicitly under-powered and reported with their exact *p*-values as hypothesis-generating.

## 10. IMbrave150 real per-patient IHC × arm external transfer
External validation of the dominant L3 feature (intratumoural CD8 density) used the **real per-patient IMbrave150 Figure 3b source table** (`data/clinical_trials/IMbrave150/41591_2022_1868_MOESM4_ESM.xlsx`, sheet *Fig3b*, *n* = 177 biomarker-evaluable patients with treatment arm, PFS, OS and intratumoural CD8 IHC %; note that the MOESM convention is `censor = 1, event = 0`, so `event_observed = 1 − censor`). Patients were median-split into CD8-High / CD8-Low **within each arm**; within-arm Kaplan–Meier and log-rank tests used lifelines 0.30. A Cox proportional-hazards model with covariates `arm`, `CD8_High` and the interaction `arm × CD8_High` (reference: sorafenib, CD8-Low) tests whether CD8 density has a **predictive** (treatment-modifying) rather than purely prognostic effect.

Results: **Atezo + Bev arm** PFS log-rank *p* = 0.009 (median 8.4 vs 5.6 months), OS *p* = 0.32; **sorafenib arm** PFS *p* = 0.90 (median 4.2 vs 5.2 months), OS *p* = 0.096; **Cox arm × CD8 interaction** HR = 0.57 (95 % CI 0.25–1.30, *p* = 0.18) for PFS and HR = 0.36 (0.12–1.07, *p* = 0.065) for OS. The full 15-feature ensemble is **not** directly applicable to bulk/IHC data because its CNV, programme-level and composition features cannot be faithfully reconstructed from bulk expression or a single IHC stain; we therefore frame this analysis as an external transfer of the dominant composition signal rather than validation of the full five-layer predictor. Outputs: `out/04_complementary/tables/revision_imbrave150_cd8_arm.csv`, `revision_imbrave150_cox_interaction.csv`; reproducible from `scripts/05_multilayer/05f_revision_analyses.py`.

## 11. External molecular-class concordance
For concordance with published HCC molecular-immune classifications we used **Sia et al. 2017** Immune / Immune-excluded / Immune-cold signatures and **Montironi et al. 2023** Inflamed / Non-inflamed panels (genes listed in `config/signatures.yml`). Patients were classified by median split of their mean Sia-Immune score and evaluated against our internal subtypes with χ² and Cohen's κ. The reviewer's earlier concern about an incorrectly attributed molecular-class reference has been resolved by withdrawing that citation (formerly ref [19]); no claim in this manuscript depends on it.

## 12. Pseudotime
Diffusion pseudotime (Scanpy) on the Harmony-integrated CD8 compartment of the 24 GSE206325 analytical patients (235,177 cells). Root = CD8 naïve-like cluster. DPT values rescaled to [0,1]. Comparisons of R vs NR use **per-patient mean pseudotime** (one observation per patient) and Mann–Whitney U with Cliff's δ. Cell-level *p*-values are not reported because they pseudoreplicate across cells from the same donor.

## 13. Statistical and reproducibility notes
All continuous comparisons used two-sided non-parametric tests (Mann–Whitney U, Kolmogorov–Smirnov) with BH correction within each analysis family. AUC 95 % CIs: 2,000-iteration patient-level bootstrap. No imputation was performed on scRNA-seq counts.

**Software.** Python 3.12.11, Scanpy 1.11.5, scikit-learn 1.7, scipy 1.16, pandas 2.3, numpy 2.2, harmonypy 0.2, lifelines 0.30, squidpy 1.6.

**Pipeline.** SLURM + conda environment in `environment.yaml`; scripts under `scripts/`, orchestrated by `slurm/run_full_pipeline.sh`. Repository: **https://github.com/isadeghi87/spatial_immunotherapy_response** (MIT licence).

## 14. Per-datatype processing summary
Droplet scRNA-seq (GSE206325, GSE125449, GSE140228, GSE149614, Cappuyns 2023): Scanpy QC (`n_genes ≥ 300`, `pct_mt ≤ 20 %`, Scrublet), log1p-CPM, 3,000 HVGs per batch, Harmony on top 50 PCs keyed by `dataset`, marker-based cell-type annotation. Non-tumour liver reference (E-MTAB-10553): published Seurat combined Wang 2021 + Ramachandran 2019 object converted via SeuratDisk. Visium spatial (GSE173905, GSE238264): Squidpy spot-level QC (`≥ 500 counts`, `≥ 200 genes`), CPM + log1p, per-section signature scoring with atlas-derived gene lists; no cross-section integration. Bulk deconvolution replaced by IMbrave150 real per-patient CD8 IHC × arm analysis (MOESM4, Finn 2022 Nat Med). CheckMate-459 bulk deconvolution withdrawn (synthetic curves). Trial metadata analysed with lifelines 0.30.

## 15. Atlas integration diagnostics
Pre-integration PCA and Harmony embeddings were benchmarked on a 60,000-cell subsample using ASW-batch, ASW-bio, iLISI, cLISI, and kBET on the 30-NN graph. scVI (scvi-tools) was attempted as an orthogonal integration but was skipped when not available in the environment. Measured values for the current atlas (Harmony vs PCA): ASW-batch 0.472 vs 0.502, ASW-bio 0.488 vs 0.499, iLISI median 2.07 vs 1.07, cLISI median 2.85 vs 2.06, kBET-acceptance 1.00 on both. Reported in Supplementary Table S20 and Supplementary Fig. S13.

## 16. Comprehensive statistical analysis catalog
Patient-level two-sided Mann–Whitney U + BH for all R-vs-NR group comparisons (NMF programmes, IO signatures, cell-type fractions, CNV summaries). Per-bin CNV R-vs-NR: Mann–Whitney + BH on patient × 1,265-bin residual matrix. CNV-mapped gene-set enrichment: Fisher exact + BH vs MSigDB Hallmark. Classifier performance: ROC-AUC / PR-AUC / balanced accuracy under LOOCV with 2,000-iteration patient bootstrap 95 % CIs and 1,000-iteration label-shuffle permutation null. Variance partitioning: common/unique decomposition on Brier-score improvement. Survival: within-arm Kaplan–Meier, log-rank, and Cox arm × CD8 interaction (lifelines 0.30). Cross-cohort concordance: Cohen's κ. Integration QC: ASW, iLISI, cLISI, kBET. Pseudoreplication avoided by patient-level unit of analysis throughout.

## 17. Software and environment manifest
Python 3.12.7, Scanpy 1.11.5, anndata 0.10.9, scikit-learn 1.7.2, SciPy 1.16.2, NumPy 2.1.3, pandas 2.3.1, matplotlib 3.10.1, seaborn 0.13.2, harmonypy 0.2.5, InferCNV 1.22 (R), Squidpy 1.6.4, lifelines 0.30.0, statsmodels 0.14.4, scvi-tools 1.2.1 (optional), gseapy 1.1.5 (optional; curated Hallmark fallback), WeasyPrint 68.1, Pandoc 3.8.2.1. Released as Supplementary Table S14.

## References
Full reference list in [`docs/manuscript.md`](manuscript.md) §6.
