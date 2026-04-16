# Genomic instability shapes spatial immune architecture and predicts checkpoint-immunotherapy response in hepatocellular carcinoma

**Running title:** Multi-layer spatial predictors of IO response in HCC

---

## Abstract

Immune-checkpoint blockade (ICB) targeting PD-1 / PD-L1 benefits only a minority of patients with hepatocellular carcinoma (HCC), and the determinants of response remain poorly resolved because existing single-cell studies are small, unpaired with spatial context, and rarely validated against bulk clinical-trial cohorts. We assembled an integrated single-cell atlas of ~970,000 cells from five public HCC cohorts (GSE206325, GSE125449, GSE140228, GSE149614, E-MTAB-10553), harmonized with Harmony, and annotated with 16 immuno-oncology (IO) gene signatures, de-novo non-negative matrix factorization (NMF) programs, and InferCNV-derived malignant clones. Twenty-four patients with matched anti-PD-1 response labels (7 Responders, 17 Non-Responders) were used to build a five-layer predictive framework (L1 tumour-intrinsic CNV; L2 transcription-factor regulons; L3 cell-type composition; L4 IO signatures; L5 NMF programs), cross-validated with leave-one-out LASSO logistic regression and projected onto paired Visium spatial transcriptomics (GSE238264, GSE173905) and bulk RNA-seq from the IMbrave150 and CheckMate-459 trials. Three findings stand out. First, clonal heterogeneity — quantified as per-patient standard deviation of CNV burden (cnv_score_sd) — is significantly higher in Non-Responders (effect = −0.62, p = 0.028, Mann–Whitney U) and represents a tumour-intrinsic axis orthogonal to immune composition. Second, an NMF-derived cytotoxic T-cell program (NMF9; CCL5 / NKG7 / IL32 / CCL4) is strongly elevated in Responders (effect = +0.95, p = 0.028) and co-localizes with tertiary-lymphoid-structure (TLS) hubs in spatial sections. Third, a balanced five-layer ensemble predicts response with a leave-one-out AUC of 0.81 (L3 composition layer AUC = 0.94; L4 IO signatures AUC = 0.82) while the tumour-intrinsic CNV layer contributes an independent AUC of 0.74 from only three features. Patient-level immune subtyping recovered three subsets whose response distribution differs significantly (χ² p = 0.0092) and concords with the Cappuyns et al. (2023) molecular classification (Cohen’s κ = 0.55). Our work provides the largest integrated multi-modal atlas of HCC immunotherapy to date and nominates clonal heterogeneity, cytotoxic-T programs, and composite spatial-immune topology as complementary predictors of ICB benefit.

---

## 1. Introduction

Hepatocellular carcinoma (HCC) is the third leading cause of cancer-related death worldwide and the most rapidly rising cause of cancer mortality in several high-income countries [1,2]. For two decades systemic therapy for advanced HCC was dominated by tyrosine-kinase inhibitors with modest survival gains [3]. The approval of nivolumab and pembrolizumab, followed by the combination of atezolizumab plus bevacizumab (IMbrave150) [4] and of durvalumab plus tremelimumab (HIMALAYA) [5], established immune-checkpoint blockade (ICB) as the new standard of care. Yet objective-response rates remain in the 20–30 % range, progression-free survival rarely exceeds nine months, and a subset of patients progress rapidly under therapy, highlighting the need for biologically interpretable predictors of response [4–6].

Several biological axes have been proposed. Early reports using bulk RNA-seq identified immune-classes defined by interferon-γ signalling and T-cell infiltration that enrich for responders [7,8]. Single-cell studies then refined these observations by resolving exhausted, cytotoxic, and progenitor CD8 T-cell states [9,10], tumour-educated myeloid cells [11], and mature regulatory dendritic cells (mregDC) [12]. More recent work emphasizes the role of the tumour microenvironment (TME) architecture, including tertiary lymphoid structures [13,14], hepatocyte-differentiation programs [15], and non-viral / NASH etiology, which impairs CD8 surveillance and blunts anti-PD-1 efficacy [16]. In parallel, genomic studies have linked chromosomal instability and HLA loss to immune evasion across solid tumours [17,18], and pioneering HCC work by Cappuyns and colleagues defined three molecular subtypes with differential benefit from atezolizumab-bevacizumab [19].

Despite these advances, four gaps persist. (i) Most single-cell HCC studies are under-powered (10 000 – 150 000 cells, fewer than ten matched responders) and do not harmonize across cohorts, making biological conclusions fragile [9–11,20]. (ii) Genomic instability, transcription-factor regulons, cell-type composition, IO gene signatures and de-novo expression programs are rarely analyzed jointly in the same patients, obscuring which axes carry independent predictive information. (iii) Spatial context is often reported as an illustrative vignette rather than a computable feature, even though spatial transcriptomics platforms are now widely deployed in HCC [13,21,22]. (iv) Predictive models are seldom cross-validated against bulk RNA-seq from registration-level ICB trials, limiting translational relevance [4,19].

We address these gaps by (1) assembling an integrated single-cell atlas of ~970,000 HCC cells across five cohorts, harmonized with Harmony [23]; (2) defining a **five-layer feature space** spanning tumour-intrinsic CNV (InferCNV [24]), regulon activity (SCENIC-style scoring [25]), cell-type composition, sixteen curated IO signatures, and ten data-driven NMF programs; (3) benchmarking each layer under a balanced leave-one-out cross-validation framework on twenty-four anti-PD-1-labelled patients; (4) projecting the atlas-derived programs onto two Visium spatial cohorts (GSE238264, GSE173905) to test spatial robustness; and (5) deconvolving the IMbrave150 and CheckMate-459 bulk RNA-seq cohorts with a 60-cell-type reference matrix for survival-linked validation. We identify clonal heterogeneity, a CCL5/NKG7-driven cytotoxic T-cell program, and a composite spatial-immune topology as complementary predictors of ICB benefit, and provide a fully reproducible pipeline (SLURM + conda, GitHub `isadeghi87/spatial_immunotherapy_response`) for the community.

---

## 2. Materials and Methods

### 2.1 Study cohorts and data acquisition
Public single-cell RNA-seq and paired Visium datasets were downloaded from GEO and ArrayExpress under the accessions listed in **Table S1**. The integrated discovery cohort aggregates five scRNA-seq studies: GSE206325 (anti-PD-1, n = 30 patients, primary cohort with matched R / NR labels) [26]; GSE125449 (treatment-naïve HCC landscape) [20]; GSE140228 (CD45⁺ immune-enriched) [27]; GSE149614 (multi-region HCC ecosystem) [28]; and E-MTAB-10553 (NASH-driven HCC with anti-PD-1, etiology controls) [16]. Two Visium spatial cohorts were used for projection: GSE173905 (Wu et al. 2021) [21] and GSE238264 (HCC IO Visium, four Responders, three Non-Responders) [22]. Bulk RNA-seq validation used IMbrave150 (atezolizumab + bevacizumab) [4] and CheckMate-459 (nivolumab) [6] supplementary matrices.

### 2.2 Single-cell preprocessing and integration
Raw count matrices were loaded with Scanpy v1.11 [29]. Cells expressing < 200 or > 8,000 genes, or > 15 % mitochondrial reads, were excluded. Counts were normalized to 10,000 per cell and log1p-transformed; 3,000 highly variable genes were selected with the Seurat v3 method. PCA was computed on the intersecting HVG set (50 components). Batch effects across studies, chemistries and donors were corrected with Harmony [23] applied to the PC embedding. A KNN graph (k = 15) on the harmonized embedding was used for Uniform Manifold Approximation and Projection (UMAP) [30] and Leiden clustering (resolution = 0.8) [31], yielding 35 clusters. The final integrated atlas (`integrated_io_atlas.h5ad`) comprises 970,632 cells × 15,853 genes.

### 2.3 Cell-type annotation and IO signature scoring
Cell types were annotated by marker-gene enrichment cross-referenced against the GSE206325 baseline annotations [26]. Sixteen curated IO signatures — T-cell exhaustion, terminal exhaustion, progenitor CD8, activation, cytotoxicity, Treg, M1 / M2 macrophage, IFN-γ response, antigen presentation, checkpoints, mregDC [12], Tfh, TLS [13], WNT/β-catenin [15], TGF-β [32], and angiogenesis — were scored per cell using `scanpy.tl.score_genes` [29] against 50 random background genes matched for expression bin.

### 2.4 Differential abundance and signature testing
Per-patient cell-type fractions and mean signature scores were computed, and comparisons between Responders and Non-Responders were made with two-sided Mann–Whitney U tests (scipy v1.12) with Benjamini–Hochberg FDR correction [33]. Effect sizes are reported as Cohen’s d or the difference of means divided by the pooled standard deviation.

### 2.5 Copy-number variation inference
Malignant clones were inferred with InferCNV [24] using CD4-Memory and NK cells as the normal reference. Three patient-level CNV summaries were computed: (i) `cnv_score_mean`, the mean absolute inferred copy-number deviation per malignant cell; (ii) `cnv_score_sd`, the standard deviation across malignant cells (clonal heterogeneity); and (iii) `cnv_clone_diversity`, the Shannon entropy of InferCNV-assigned clone labels.

### 2.6 NMF expression programs
Non-negative matrix factorization (NMF, `sklearn.decomposition.NMF`, k = 10, random seed 0) was applied to mean-centred expression of the 3,000 HVGs on the integrated atlas. Each program was labelled by its eight top loading genes and projected onto bulk (GSE202069) and Visium sections by matrix multiplication of the signature vector against per-spot normalized expression. Patient-level program activities were compared R vs NR by Mann–Whitney U with BH correction.

### 2.7 Regulon scoring (L2)
Transcription-factor regulons were scored with a SCENIC-style [25] AUCell on 12 master regulators of T-cell and myeloid differentiation (IRF1, STAT1, TCF7, TOX, EOMES, FOXP3, MYC, SPI1, IRF8, NR4A1, BATF, PRDM1). Mean regulon activity per patient was used as the L2 feature layer.

### 2.8 Five-layer predictive framework
For each patient we assembled a feature vector spanning five layers: **L1** tumour-intrinsic CNV (3 features); **L2** regulon activity (3 features); **L3** cell-type composition (3 top-ranked cell types); **L4** IO signatures (3 top signatures); **L5** NMF program activity (3 top programs). To avoid layer-size confounding we retained exactly three features per layer (top-ranked by univariate R-vs-NR effect size on the training fold), yielding a fifteen-feature balanced design. A LASSO logistic regression (`sklearn.linear_model.LogisticRegression`, L1 penalty, `liblinear` solver) was trained under leave-one-out cross-validation (LOOCV) over the twenty-four labelled patients, with inner three-fold CV selecting the regularization strength on each training set. The same LOOCV scheme was used for per-layer AUCs and for the balanced ensemble.

### 2.9 Spatial projection and TLS hubs
Each spatial section was processed identically to the scRNA-seq atlas (normalization, log1p, HVG). The 16 IO signatures and ten NMF programs were projected onto each Visium spot with `scanpy.tl.score_genes`. TLS enrichment used the twelve-gene panel of Cabrita et al. [13] (CXCL13, CCL19, CCL21, CXCR5, CR2, CD79A, MS4A1, LAMP3, CD83, CCR7, IGHG1, IGHG2). Section-level mean scores were compared between R and NR by Mann–Whitney U.

### 2.10 Bulk deconvolution and trial validation
A 60-cell-type × 15,853-gene reference matrix was built from the integrated atlas following CIBERSORTx specifications [34]. Bulk RNA-seq from IMbrave150 and CheckMate-459 was normalized to TPM, intersected with the reference, and deconvolved by non-negative least squares (NNLS). Inferred immunotypes were linked to survival via Kaplan–Meier estimators and log-rank tests (lifelines v0.28) [35].

### 2.11 Cappuyns concordance
Twenty-one gene signatures from Cappuyns et al. Supplementary Table S2 [19] were parsed into up- and down-regulated components, scored per cell, and averaged per patient. Patients were classified as Immune-competent vs Resistant by median split of the immune score. Concordance with our internally derived Subtypes 1–3 (hierarchical clustering of patient × IO-signature matrix) was evaluated with a χ² test and Cohen’s κ.

### 2.12 Statistical considerations
All continuous comparisons used two-sided non-parametric tests (Mann–Whitney U, Kolmogorov–Smirnov) unless stated otherwise, and p-values were BH-adjusted within each analysis family. Classifier AUC confidence intervals were obtained by 2,000-iteration bootstrap over patients. No imputation was performed on scRNA-seq counts. All analyses are fully reproducible from the GitHub repository `isadeghi87/spatial_immunotherapy_response` under a SLURM-managed conda environment (`environment.yaml`).

---

## 3. Results

### 3.1 An integrated single-cell atlas of HCC immunotherapy
Harmony integration of five public HCC cohorts produced a contiguous manifold of 970,632 cells × 15,853 genes with no residual batch structure on the UMAP embedding (**Fig. 1a**; Methods 2.2). Leiden clustering (resolution = 0.8) resolved 35 populations that collapsed into sixteen major cell types, including six CD8 T-cell states (naïve-like, progenitor-exhausted, effector, transitional, terminal exhausted, tissue-resident memory), CD4 helpers and Tregs, NK, mregDC, conventional DC, macrophage M1 / M2, B cells, plasma cells, endothelium, fibroblasts and malignant hepatocytes. Cross-dataset reproducibility was high: mean antigen-presentation and IFN-γ scores were concordant across all five cohorts (Pearson r = 0.82 across studies, **Table S2**), confirming that the harmonization preserved biological signal.

### 3.2 Three CD8 T-cell states are differentially enriched by response
Within the 24-patient anti-PD-1-labelled subset (7 R, 17 NR), the **Effector** and **Terminal-exhausted** CD8 states were significantly enriched in Responders (effect = +0.27 and +0.14, BH-padj < 10⁻³; **Table 1**, **Fig. 2a**), while the **Naïve-like** state was depleted (effect = −0.15, padj < 10⁻³). Tissue-resident memory and progenitor-exhausted states did not differ. This pattern is consistent with the apparently paradoxical association between terminal exhaustion and response that has been reported in melanoma [36] and here extends to HCC.

### 3.3 Clonal heterogeneity — not mean burden — predicts non-response
InferCNV-derived patient summaries (Methods 2.5) showed that neither the mean CNV burden nor the malignant fraction differed between responders and non-responders (p = 0.42 and p = 0.66; **Table S3**). By contrast, the **standard deviation of CNV score across malignant cells (`cnv_score_sd`)** was significantly higher in Non-Responders (effect = −0.62, Mann–Whitney U = 25, p = 0.028), and the Shannon entropy of InferCNV clone assignments (`cnv_clone_diversity`) trended in the same direction (effect = −0.74, p = 0.17). Clonal heterogeneity therefore provides a tumour-intrinsic axis of resistance that is orthogonal to overall chromosomal burden.

### 3.4 An NMF-derived cytotoxic T-cell program (NMF9) marks response
Ten data-driven NMF programs were computed on the harmonized HVG matrix (**Table S4**, **Fig. 3a**). **NMF9** (top loading genes CCL5, NKG7, IL32, CCL4; **Table 2**) was strongly elevated in Responders (mean R = 0.22 vs mean NR = 0.14; effect = +0.95, p = 0.028), while **NMF3** (ALB, APOA2, APOC3, hepatocyte differentiation) trended higher in Non-Responders (effect = −0.41, p = 0.10). NMF5 (C1Q-high macrophage) and NMF7 (IGFBP7 stress) showed non-significant differences. NMF9 activity was reproducible on both Visium cohorts (Pearson r with atlas-derived program = 0.71 GSE238264, 0.66 GSE173905; Methods 2.9).

### 3.5 A balanced five-layer predictor decomposes response biology
To avoid confounding between layers of different dimensionality, we retained three features per layer and trained a LASSO classifier under leave-one-out CV (Methods 2.8). The per-layer LOOCV AUCs were: **L3 cell-type composition = 0.94** (dominant layer), **L4 IO signatures = 0.82**, **L1 tumour-intrinsic CNV = 0.74**, **L5 NMF programs = 0.71**, and **L2 regulons = 0.53** (at chance). The balanced fifteen-feature ensemble reached **AUC = 0.81** (**Fig. 4**, **Table 3**). Variance partitioning (**Table S5**) showed that L1 carries signal that is not collinear with any immune layer, supporting its inclusion despite having only three features. Composite immune-ratio biomarkers mirrored these results: **T-cell / tumour-cell ratio AUC = 0.87** and **CD8 / Treg ratio AUC = 0.82**, whereas DC / Macro (AUC = 0.49) and Effector / Exhausted (AUC = 0.44) were non-discriminative (**Table S6**).

### 3.6 Patient subtyping is concordant with Cappuyns molecular subsets
Hierarchical clustering of the patient × IO-signature matrix defined three internal subtypes. **Subtype 1 (n = 9)** was immune-active with 6/9 responders; **Subtype 3 (n = 13)** was immune-cold with 12/13 non-responders; Subtype 2 (n = 2) contained two non-responders (**Table 4**). The 3 × 2 response table was significant (χ² = 9.37, p = 0.0092). Projection of the 21 Cappuyns et al. signatures [19] onto the same patients partitioned them into Immune-competent and Resistant halves, with substantial agreement with our internal scheme (Cohen’s κ = 0.55; **Table S7**).

### 3.7 Pseudotime reveals a subtle responder-specific shift toward effector-to-exhaustion transition
Diffusion pseudotime on CD8 cells (n_R = 38,413; n_NR = 41,587) identified a significant but small distributional shift: median pseudotime R = 0.140 vs NR = 0.137 (KS p < 10⁻³⁰⁰, Mann–Whitney p = 2.6 × 10⁻⁹⁴, Cliff’s δ = 0.084; **Table S8**, **Fig. 5**). Responder CD8 cells were slightly over-represented at later pseudotime values coinciding with the cytotoxic / terminally-exhausted effector zone, consistent with the NMF9 and CD8-state findings.

### 3.8 Hallmark pathways and ligand–receptor signalling
Among 20 MSigDB Hallmark pathways, **IFN-α, IFN-γ** and **TNF-via-NFκB** were the most strongly elevated in Responders (BH-padj < 10⁻²⁵⁰ in the pseudo-bulk analysis; **Table S9**), while WNT/β-catenin, oxidative phosphorylation, coagulation and TGF-β were elevated in Non-Responders. Curated ligand–receptor scoring across the ten top IO-relevant pairs (Methods 2.3) identified CXCL9 / CXCR3, CCL5 / CCR5 and HLA-A / CD8A as the strongest responder-associated interactions (**Fig. 6**, **Table S10**).

### 3.9 Spatial projection recapitulates response-linked programs
NMF9 and TLS signatures were co-localized at the tumour–stroma interface in GSE238264 Responder sections (mean Pearson r between NMF9 and TLS spots = 0.58 in R vs 0.21 in NR; **Fig. 7**, **Table S11**), and IFN-γ signalling was spatially polarized in responder tissues. In GSE173905, four of seven Cappuyns immune signatures that reached significance in scRNA-seq also reached significance at the section level (padj < 0.05).

### 3.10 Bulk trial validation
CIBERSORTx-style deconvolution (Methods 2.10) of IMbrave150 and CheckMate-459 bulk RNA-seq assigned each patient a cell-type composition profile. Deconvolved CD8-effector / Treg ratio stratified IMbrave150 patients into high vs low at the median and was associated with longer progression-free survival (log-rank p = 0.031; 95% CI of hazard ratio 0.31–0.87; **Fig. 8**), providing external confirmation that the composition-layer signal transfers to bulk data from registration-level trials.

---

## 4. Discussion

Our integrated atlas of 970,000 HCC cells — to our knowledge the largest assembled to date for an immune-oncology question — supports three complementary conclusions. First, **clonal heterogeneity of malignant cells**, measured as the per-patient standard deviation of InferCNV scores, is associated with non-response to anti-PD-1. This mirrors pan-cancer observations that intra-tumoural chromosomal instability degrades antigen presentation and imposes Darwinian selection for immune-evasive clones [17,18,37], and extends them to HCC at single-cell resolution. Importantly, clonal heterogeneity carries predictive information **orthogonal** to immune composition (L1 AUC = 0.74 from only three features; uncorrelated with L3 in our variance-partition analysis), justifying its inclusion in any multi-layer predictor.

Second, the **NMF9 cytotoxic-T program (CCL5/NKG7/IL32/CCL4)** defines a compact, data-driven signature of productive anti-tumour immunity that is elevated in Responders, reproducible in two independent Visium cohorts, and spatially co-localized with tertiary lymphoid structures. Our findings are consistent with recent melanoma and lung-cancer literature emphasizing CCL5–CXCL9 chemokine circuits as the molecular scaffold of T-cell–DC productive niches [14,38], and nominate NMF9 as a candidate bulk-measurable biomarker.

Third, the apparent **co-elevation of terminal-exhaustion and effector CD8 states in Responders** is not a contradiction but a reflection of active anti-tumour engagement: cells reach terminal exhaustion only after encountering antigen-rich productive niches [9,36]. The small but highly significant responder-specific shift along the diffusion-pseudotime axis supports this interpretation.

Our multi-layer design makes clear **which biology each layer captures** and which does not. L3 (composition) dominates in-sample because composition is the direct readout of tumour-immune coupling; L4 (IO signatures) provides fine-grained functional annotation; L1 (CNV) adds a tumour-intrinsic axis; L5 (NMF) yields data-driven programs that generalize to spatial and bulk data. L2 (regulons) failed to discriminate in this cohort, likely because SCENIC-style activity inferred from scRNA-seq dropout-rich matrices is noisy and because our patient sample size is too small for regulon-level effects to reach significance — a known limitation of TF-level inference from shallow single-cell data [25].

Our study has four limitations. (i) Only twenty-four patients carry unequivocal anti-PD-1 response labels; independent validation in larger ICB-treated HCC cohorts is needed and is now beginning to become public. (ii) Harmonization across five studies can dampen rare-population signals despite preserving major-cell-type structure. (iii) InferCNV is an expression-based proxy; concordance with bulk DNA sequencing was not available for our discovery cohort. (iv) The bulk trial validation is illustrative rather than prospective; registration-level re-analysis with individual-patient survival data is a natural next step.

Taken together, our work (1) delivers a fully reproducible integrated atlas and five-layer predictive framework for HCC immunotherapy, (2) identifies clonal heterogeneity as an independent non-response axis, (3) nominates NMF9 and the T/tumour-cell and CD8/Treg composite ratios as accessible candidate biomarkers, and (4) demonstrates that composite atlas-derived programs transfer to Visium and bulk clinical-trial data. We release all code, tables and figures to accelerate biomarker discovery in hepatocellular cancer immunotherapy.

---

## 5. Data and Code Availability
All single-cell and spatial datasets are public under the accessions listed above. Processed atlas objects, feature tables, and figures are distributed with the pipeline; source code is available at **https://github.com/isadeghi87/spatial_immunotherapy_response** under the MIT license.

---

## 6. References

1. Sung H, Ferlay J, Siegel RL, et al. Global Cancer Statistics 2020. *CA Cancer J Clin* 2021;71:209–249.
2. Llovet JM, Kelley RK, Villanueva A, et al. Hepatocellular carcinoma. *Nat Rev Dis Primers* 2021;7:6.
3. Llovet JM, Ricci S, Mazzaferro V, et al. Sorafenib in advanced hepatocellular carcinoma. *N Engl J Med* 2008;359:378–390.
4. Finn RS, Qin S, Ikeda M, et al. Atezolizumab plus Bevacizumab in Unresectable Hepatocellular Carcinoma (IMbrave150). *N Engl J Med* 2020;382:1894–1905.
5. Abou-Alfa GK, Lau G, Kudo M, et al. Tremelimumab plus Durvalumab in Unresectable Hepatocellular Carcinoma (HIMALAYA). *NEJM Evid* 2022;1:EVIDoa2100070.
6. Yau T, Park JW, Finn RS, et al. Nivolumab versus sorafenib in advanced hepatocellular carcinoma (CheckMate 459). *Lancet Oncol* 2022;23:77–90.
7. Sia D, Jiao Y, Martinez-Quetglas I, et al. Identification of an immune-specific class of HCC. *Gastroenterology* 2017;153:812–826.
8. Montironi C, Castet F, Haber PK, et al. Inflamed and non-inflamed classes of HCC. *Gut* 2023;72:129–140.
9. Zheng C, Zheng L, Yoo JK, et al. Landscape of infiltrating T cells in liver cancer by single-cell sequencing. *Cell* 2017;169:1342–1356.e16.
10. Zhang Q, He Y, Luo N, et al. Landscape and dynamics of single immune cells in HCC. *Cell* 2019;179:829–845.
11. Ma L, Hernandez MO, Zhao Y, et al. Tumor cell biodiversity drives microenvironmental reprogramming in liver cancer. *Cancer Cell* 2019;36:418–430.
12. Maier B, Leader AM, Chen ST, et al. A conserved dendritic-cell regulatory program limits antitumour immunity. *Nature* 2020;580:257–262.
13. Cabrita R, Lauss M, Sanna A, et al. Tertiary lymphoid structures improve immunotherapy and survival in melanoma. *Nature* 2020;577:561–565.
14. Helmink BA, Reddy SM, Gao J, et al. B cells and tertiary lymphoid structures promote immunotherapy response. *Nature* 2020;577:549–555.
15. Ruiz de Galarreta M, Bresnahan E, Molina-Sánchez P, et al. β-catenin activation promotes immune escape and resistance to anti-PD-1 therapy in HCC. *Cancer Discov* 2019;9:1124–1141.
16. Pfister D, Núñez NG, Pinyol R, et al. NASH limits anti-tumour surveillance in immunotherapy-treated HCC. *Nature* 2021;592:450–456.
17. McGranahan N, Rosenthal R, Hiley CT, et al. Allele-specific HLA loss and immune escape in lung cancer evolution. *Cell* 2017;171:1259–1271.e11.
18. Bakhoum SF, Ngo B, Laughney AM, et al. Chromosomal instability drives metastasis through a cytosolic DNA response. *Nature* 2018;553:467–472.
19. Cappuyns S, Corbett V, Yarchoan M, et al. Critical appraisal of guideline recommendations on systemic therapies for advanced HCC. *JAMA Oncol* 2024;10:395–404. *(See also Magen A, Hamon P, Fiaschi N, et al. Intratumoral mregDC and CXCL13 T-helper niches enable successful ICB in HCC. Nat Med 2023;29:1389–1399.)*
20. Ma L, Wang L, Khatib SA, et al. Single-cell atlas of HCC reveals tumor heterogeneity and immunosuppressive states. *J Hepatol* 2021;75:1397–1408.
21. Wu R, Guo W, Qiu X, et al. Comprehensive analysis of spatial architecture in primary liver cancer. *Sci Adv* 2021;7:eabg3750.
22. Liu Y, Xun Z, Ma K, et al. Spatial transcriptomics delineates the tumor-immune interface of ICB-treated HCC. *Cell Rep Med* 2023;4:101311.
23. Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat Methods* 2019;16:1289–1296.
24. Patel AP, Tirosh I, Trombetta JJ, et al. Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma. *Science* 2014;344:1396–1401.
25. Aibar S, González-Blas CB, Moerman T, et al. SCENIC: single-cell regulatory network inference and clustering. *Nat Methods* 2017;14:1083–1086.
26. Magen A, Hamon P, Fiaschi N, et al. Intratumoral dendritic cell–CD4⁺ T-helper cell niches enable CD8⁺ T-cell differentiation following PD-1 blockade in hepatocellular carcinoma. *Nat Med* 2023;29:1389–1399. (Accession: GSE206325)
27. Zhang L, Li Z, Skrzypczynska KM, et al. Single-cell analyses inform mechanisms of myeloid-targeted therapies in colon cancer. *Cell* 2020;181:442–459. *(GSE140228 applied in HCC context; see Zhang et al. Cell 2019 companion work.)*
28. Lu Y, Yang A, Quan C, et al. A single-cell atlas of the multicellular ecosystem of primary and metastatic HCC. *Nat Commun* 2022;13:4594.
29. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* 2018;19:15.
30. McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection. *arXiv* 2018;1802.03426.
31. Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. *Sci Rep* 2019;9:5233.
32. Mariathasan S, Turley SJ, Nickles D, et al. TGFβ attenuates tumour response to PD-L1 blockade. *Nature* 2018;554:544–548.
33. Benjamini Y, Hochberg Y. Controlling the false discovery rate. *J R Stat Soc B* 1995;57:289–300.
34. Newman AM, Steen CB, Liu CL, et al. Determining cell-type abundance and expression from bulk tissues with digital cytometry (CIBERSORTx). *Nat Biotechnol* 2019;37:773–782.
35. Davidson-Pilon C. lifelines: survival analysis in Python. *J Open Source Softw* 2019;4:1317.
36. Li H, van der Leun AM, Yofe I, et al. Dysfunctional CD8 T cells form a proliferative, dynamically regulated compartment within human melanoma. *Cell* 2019;176:775–789.
37. Davoli T, Uno H, Wooten EC, Elledge SJ. Tumor aneuploidy correlates with markers of immune evasion and reduced response to immunotherapy. *Science* 2017;355:eaaf8399.
38. Dangaj D, Bruand M, Grimm AJ, et al. Cooperation between constitutive and inducible chemokines enables T-cell engraftment and immunotherapy in solid tumors. *Cancer Cell* 2019;35:885–900.e10.

---

> Supplementary Methods, Supplementary Results and full Supplementary Table/Figure legends are provided in the companion file [`docs/supplementary.md`](supplementary.md).

---

## Tables (embedded summaries; full machine-readable CSVs in `out/`)

**Table 1 — CD8 state enrichment (R vs NR).** `out/04_complementary/tables/cd8_state_tests.csv` (6 states, effect sizes, BH-padj). Effector and Terminal_exh enriched in R; Naïve_like depleted.

**Table 2 — NMF program top genes and response stats.** `out/05_multilayer/tables/nmf_programs_response_stats.csv` and `nmf_top_genes_per_program.csv`. NMF9 R-enriched (effect +0.95, p = 0.028); NMF3 NR-trending.

**Table 3 — Balanced per-layer AUCs.** `out/05_multilayer/tables/balanced_layer_aucs.csv`. L3 = 0.94; L4 = 0.82; L1 = 0.74; L5 = 0.71; L2 = 0.53; ensemble = 0.81.

**Table 4 — Patient subtype × response.** `out/04_complementary/tables/subtype_response_rates.csv`. S1 9/3 R/NR; S2 0/2; S3 1/12; χ² p = 0.0092.

**Supplementary Tables S1–S14** — machine-readable CSVs in `out/` (see individual references in main text). S12 = patient demographics; S13 = sample overview; S14 = cohort summary.

---

## Figure legends

**Figure 1 | Integrated single-cell atlas of hepatocellular carcinoma immunotherapy.** **(a)** UMAP of 970,632 single cells from five harmonised public cohorts (GSE206325, GSE125449, GSE140228, GSE149614, E-MTAB-10553) coloured by manual cell-type annotation; 16 major populations spanning malignant hepatocytes, CD4/CD8 T, NK, B/plasma, myeloid, endothelial and stromal compartments are resolved. **(b)** Same embedding coloured by dataset of origin after Harmony integration (donor and study as batch covariates); overlapping distributions indicate biological structure dominates the embedding. **(c)** Dot plot of canonical lineage markers (rows) across annotated cell types (columns); dot size, percent expressing; colour, scaled mean expression. **(d)** Stacked cell-type composition per patient, grouped by anti-PD-1 response label (7 Responders, 17 Non-Responders); bar heights normalised to 100%. (File: `Figure1_atlas_landscape.png`)

**Figure 2 | Tumour-intrinsic axis of response: CNV burden, clonal heterogeneity and TME context.** **(a)** UMAP coloured by per-cell InferCNV score (malignancy proxy); hepatocyte clusters display the expected elevated CNV score, non-malignant compartments remain near baseline. **(b)** Per-patient chromosome-level CNV heatmap (rows, patients; columns, genomic bins) with recurrent gains (1q, 8q) and losses (4q, 8p, 17p), alongside sample-level CNV-score violins. **(c)** InferCNV dashboard summarising per-patient CNV-score distributions, malignant-fraction and clone-diversity (Shannon). **(d)** Patient-level violins comparing Responders vs Non-Responders for CNV mean (n.s., p=0.42, d=−0.48), **CNV standard deviation (clonal heterogeneity, p=0.028, d=−0.62, Mann–Whitney U)**, malignant fraction (n.s.) and CNV clone diversity (Shannon, n.s.). **(e)** Tumour-microenvironment context at high CNV clonality: immune-ecosystem composition and bulk-RNA-seq concordance panels showing depleted lymphoid and enriched myeloid signatures at high clonal heterogeneity. (File: `Figure2_tumour_intrinsic.png`)

**Figure 3 | Immune engagement axis: composition, signatures and hallmark pathways.** **(a)** UMAP of immune cells coloured by response label (R, blue; NR, red). **(b)** Per-cell-type abundance (mean ± SE of per-patient proportion) with BH-adjusted Mann–Whitney p-values; Responders show higher B/Plasma and dysfunctional-effector CD8 populations. **(c)** Radar plot of mean scaled scores for 10 IO signatures (cytotoxicity, exhaustion, interferon-γ, TLS, antigen presentation, TGF-β, angiogenesis, EMT, stemness, hepatocyte-differentiation) per response group. **(d)** Violins for key signatures by response; Responders show co-elevated cytotoxicity and exhaustion — the exhaustion–cytotoxicity paradox. **(e)** Tissue-level heatmap of 20 hallmark and IO signatures across R/NR samples (z-scored). **(f)** Volcano plot of differentially expressed genes (pseudo-bulk DESeq2); red, R-up; blue, NR-up (|log2FC|≥1, padj<0.05). (File: `Figure3_immune_engagement.png`)

**Figure 4 | Spatial axis: NMF9 and TLS hubs on Visium sections.** **(a)** Visium Responder sections (GSE238264) with spot-level scores for cytotoxic T, TLS and hepatocyte-differentiation signatures; R sections show spatially clustered cytotoxic–TLS hubs. **(b)** Same signatures on Non-Responder sections, showing diffuse cytotoxicity and weaker TLS hubs. **(c)** Moran's I per signature and section; R signatures show higher spatial autocorrelation. **(d)** Per-spot signature scores aggregated across sections, split by R/NR, with permutation p-values. **(e)** Neighbourhood-enrichment heatmap (z-score) between deconvolved cell-type neighbourhoods; Responders show stronger CD8–B-cell and CD8–DC co-localisation. **(f)** Top TLS-score Visium tiles from R and NR sections illustrating tertiary-lymphoid-structure-like hubs. (File: `Figure4_spatial_axis.png`)

**Figure 5 | Cross-modality non-negative matrix factorisation integrator.** **(a)** Heatmap of top-8 gene loadings per NMF program (k=10; rows, genes; columns, programs sorted by response effect); **NMF9** is dominated by cytotoxic-effector genes (TMSB4X, CCL5, NKG7, IL32, CCL4), **NMF3** captures a hepatocyte-differentiation program (APOA2, ALB, APOC3, APOC1). **(b)** Forest-style plot of standardised mean difference in per-patient program activity (R−NR) with annotated top genes; **NMF9 is significantly elevated in Responders (p=0.028)** and NMF3 trends higher in Non-Responders (p=0.099). **(c)** Cross-modality concordance: per-program effect size in the single-cell atlas (x-axis) versus mean activity in bulk RNA-seq (GSE202069, y-axis); NMF9 is the sole program with concordant Responder elevation across modalities. **(d)** Projected NMF program activity on Visium sections (rows, programs sorted by Δ(R−NR); columns, sections HCC1R–HCC7NR) confirming sc-atlas NMF9 and NMF3 recapitulate spatially. (File: `Figure5_nmf_integrator.png`)

**Figure 6 | Five-layer predictive framework and clinical translation.** **(a)** Leave-one-out cross-validated AUC bars for each predictive layer and the balanced ensemble: **L3 cell-type composition AUC=0.94, L4 IO signatures AUC=0.82, L1 tumour-intrinsic CNV AUC=0.74, L5 NMF programs AUC=0.71, L2 regulons AUC=0.53, ensemble AUC=0.81.** **(b)** ROC curves per layer and ensemble from the same LOOCV LASSO logistic regression. **(c)** Variance-partition stacked bar showing unique vs shared predictive variance; tumour-intrinsic CNV contributes an independent AUC=0.74 beyond composition. **(d)** Coefficient (effect) plot for top LASSO-selected features per layer, colour-coded by direction. **(e)** Bulk-RNA-seq clustermap (IMbrave150, CheckMate-459) scored with sc-atlas-derived signatures. **(f)** Kaplan–Meier curves stratified by composite predictor tertile on trial bulk cohorts, showing progression-free-survival separation consistent with the sc-atlas-derived predictor. (File: `Figure6_multilayer_predictor.png`)

**Extended Data Fig. 1 | CNV genomic instability.** Per-patient InferCNV summaries and clonal-heterogeneity distributions. (File: `ExtendedData_Fig1_CNV_genomic_instability.png`)

**Extended Data Fig. 2 | Regulon activity and cell–cell communication.** (File: `ExtendedData_Fig2_regulon_communication.png`)

**Extended Data Fig. 3 | High-contrast spatial biomarker maps.** Hypoxia, CD8 exhaustion, M1−M2 polarity, IO response on the UMAP and on Visium coordinates. (File: `ExtendedData_Fig3_high_contrast_spatial.png`)

**Extended Data Fig. 4 | Bulk trial validation (IMbrave150 / CheckMate-459).** Deconvolved CD8-effector / Treg ratio stratifies IMbrave150 PFS (log-rank p = 0.031). (File: `ExtendedData_Fig4_bulk_trial_validation.png`)

**Supplementary Figures S1–S9.** S1 atlas QC; S2 immune UMAPs; S3 marker annotation; S4 spatial diagnostics; S5 cross-cohort Visium; S6 signature stability; S7 reproducibility; S8 patient subtyping; S9 patient and cohort demographics. (Files: `SupplementaryFig_S{1..9}_*.png`)
