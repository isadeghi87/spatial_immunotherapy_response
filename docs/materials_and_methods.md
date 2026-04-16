# Materials and Methods

> Companion document to `docs/manuscript.md`.
> This file mirrors Section 2 of the manuscript and expands on parameter choices.

## 1. Study cohorts and data acquisition
Public single-cell RNA-seq and paired Visium datasets were downloaded from GEO and ArrayExpress. The discovery cohort aggregates five scRNA-seq studies:

- **GSE206325** — Magen et al., *Nat Med* 2023; anti-PD-1, 30 patients, primary cohort with R / NR labels.
- **GSE125449** — Ma et al., *Cancer Cell* 2019; treatment-naïve HCC landscape, ~10 k cells.
- **GSE140228** — Zhang et al., *Cell* 2019; CD45+ immune-enriched HCC, ~67 k cells.
- **GSE149614** — Lu et al., *Nat Commun* 2022; multi-region HCC ecosystem, ~70 k cells.
- **E-MTAB-10553** — Pfister et al., *Nature* 2021; NASH-driven HCC with anti-PD-1, ~150 k cells.

Spatial validation used **GSE173905** (Wu et al., *Sci Adv* 2021) and **GSE238264** (Liu et al., *Cell Rep Med* 2023). Bulk validation used IMbrave150 (Finn et al., *NEJM* 2020) and CheckMate-459 (Yau et al., *Lancet Oncol* 2022) supplementary matrices.

## 2. Single-cell preprocessing and integration
Processed with Scanpy v1.11 (Wolf et al. 2018). QC: 200 < n_genes < 8,000; pct_mt < 15 %. Normalized to 10,000 counts per cell, log1p, 3,000 HVGs (Seurat v3), PCA (50 PCs). Batch correction with Harmony (Korsunsky et al. 2019; θ = 2, σ = 0.1, 20 iterations). KNN (k = 15), UMAP, Leiden (resolution = 0.8). Final atlas: 970,632 cells × 15,853 genes.

## 3. Cell-type annotation and IO signatures
Cell types annotated by marker-gene enrichment against GSE206325 baseline. Sixteen curated IO signatures scored per cell with `scanpy.tl.score_genes` (50 background genes, expression-bin matched).

## 4. Differential abundance and signature testing
Two-sided Mann–Whitney U with BH correction (scipy v1.12). Effect sizes = Cohen’s d or (mean_R − mean_NR)/pooled SD.

## 5. Copy-number variation (InferCNV)
InferCNV (Patel et al. 2014). Reference cells = CD4-Memory ∪ NK. Cutoff 0.1, HMM enabled, subclusters mode, denoise TRUE. Per-patient summaries: `cnv_score_mean`, `cnv_score_sd`, `cnv_clone_diversity`.

## 6. NMF expression programs
`sklearn.decomposition.NMF(k=10, init='nndsvd', random_state=0)` on log-normalized HVG matrix; k chosen by elbow and cophenetic > 0.9 across seeds. Program labels = top-8 loading genes.

## 7. Regulons (L2)
SCENIC-style AUCell (Aibar et al. 2017) on 12 master regulators (IRF1, STAT1, TCF7, TOX, EOMES, FOXP3, MYC, SPI1, IRF8, NR4A1, BATF, PRDM1).

## 8. Five-layer predictor
Exactly three features per layer (top univariate effect size on each training fold). LASSO under LOOCV (n = 24) with inner 3-fold stratified CV over C ∈ [1e-3, 10]. 2,000-iter bootstrap CIs.

## 9. Spatial projection
Visium sections processed identically to scRNA. Signatures scored per spot. TLS panel = Cabrita 2020 / Helmink 2020 (12 genes).

## 10. Bulk deconvolution (trial validation)
60-cell-type × 15,853-gene reference built CIBERSORTx-style (Newman et al. 2019). NNLS deconvolution; KM / log-rank (lifelines v0.28).

## 11. Cappuyns concordance
21 signatures from Cappuyns et al. 2023 Supplementary Table S2, scored per cell, averaged per patient, median-split into Immune-competent / Resistant. χ² and Cohen’s κ vs internal subtypes.

## 12. Reproducibility
Python 3.11, SLURM, conda env in `environment.yaml`. Scripts in `scripts/`, orchestrated by `slurm/run_full_pipeline.sh`. Publication colours in `config/publication_style.py`.

## References
Full reference list in `docs/manuscript.md` section 6.
