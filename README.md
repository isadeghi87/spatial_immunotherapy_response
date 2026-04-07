# Spatial Biomarkers of Immunotherapy Response in HCC

![Python](https://img.shields.io/badge/Python-3.11-blue)
![R](https://img.shields.io/badge/R-4.3%2B-276DC3)
![scanpy](https://img.shields.io/badge/scanpy-1.11-orange)
![License: MIT](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/Status-Active-brightgreen)

> Defining spatial immune topology biomarkers that predict response vs. resistance
> to anti-PD-1 / anti-PD-L1 checkpoint immunotherapy in hepatocellular carcinoma (HCC),
> using the largest integrated single-cell atlas of HCC immunotherapy to date (~970 k cells).

---

## Highlights

| What | Detail |
|------|--------|
| **Integrated atlas** | ~970,000 cells × 15,853 genes across 5 scRNA-seq studies, Harmony batch-corrected |
| **Response labels** | 30 patients with matched Responder / Non-Responder anti-PD-1 labels (GSE206325) |
| **16 IO signatures** | Exhaustion, activation, cytotoxicity, Treg, M1/M2, IFN-γ, TGF-β, checkpoints, mregDC, Tfh, angiogenesis, etc. |
| **Response prediction** | LASSO AUC = 0.966, Random Forest AUC = 0.941 (patient-level cross-validated) |
| **20 hallmark pathways** | GSEA-style pathway scoring comparing R vs NR |
| **24 ligand-receptor pairs** | Curated IO-relevant LR interaction analysis |
| **CD8 T cell states** | 6-state characterisation (naïve → progenitor → effector → terminal exhaustion) |
| **Bulk deconvolution** | CIBERSORTx-compatible reference matrix (15,853 genes × 60 cell types) |
| **Reproducibility** | Fully automated SLURM pipeline with dependency chains |

---



## Recent Updates (April 2026)

*   **Standardized Cell Type Annotations:** Cleaned up numerical prefixes and variation names across all immune subsets (converging variants into single categories like `NK` and `CD8-Memory`).
*   **Major/Subtype UMAP Representations:** Added major and subtype-specific immune UMAPs isolating non-immune tissues to high-contrast grey (`#E0E0E0`) within the `regen_umap.py` script.
*   **InferCNV Integration (`02g_infercnv.py`):** Added chromosomal copy number variation tracking using CD4-Memory and NK cells as references to identify and pinpoint malignant clones. Enhanced visualizations directly map horizontal genomic ranges (Chr 1-22) against malignant clone fractions.
*   **Multi-Modal CNV/TME Integration (`02h_cnv_tme_integration.py`):** Generated an integration dashboard cross-correlating isolated tumor clones against spatial TME shifts, highlighting expression pathways and translating burden events (like HLA-loss) to explicit progression probabilities using Kaplan-Meier hazard modeling.
*   **Trajectory Stream Dynamics (`04c_tcell_trajectories.py`):** Migrated T cell continuous trajectories to visualize topological progression using high-density vector arrows across the UMAP manifold and spline-interpolated (scipy.interpolate) pseudotime feature tracks.

## Repository Structure

```
spatial_immunotherapy_response/
│
├── scripts/
│   ├── 00_data_acquisition/          # Download & format conversion
│   │   ├── download_io_datasets.py        # GEO FTP downloader (4 datasets)
│   │   ├── download_direct.sh             # Shell-based GEO downloader
│   │   ├── download_additional_data.sh    # E-MTAB-10553 + GSE173905 + trial placeholders
│   │   ├── download_trial_data.sh         # IMbrave150 / CheckMate-459 XLSX
│   │   ├── convert_rda_to_h5ad.R          # GSE206325 SCE .Rda → MTX + CSV
│   │   └── 00b_convert_emtab10553.R       # E-MTAB-10553 Seurat .rds → MTX + CSV
│   │
│   ├── 01_preprocessing/
│   │   └── 01_preprocess_integrate.py     # QC → normalize → HVG → PCA → Harmony → UMAP → Leiden
│   │
│   ├── 02_analysis/
│   │   ├── 02a_response_tme_characterization.py  # Cell types, 16 IO signatures, DEGs, 8 figures
│   │   ├── 02b_spatial_biomarker_discovery.py     # Pseudo-bulk / Visium spatial signatures
│   │   └── 02c_bulk_deconvolution.py              # CIBERSORTx reference + NNLS deconvolution
│   │
│   ├── 03_validation/
│   │   └── 03_validation_downstream.py    # LASSO, RF, NMF, TCR, bulk validation, bootstrap
│   │
│   └── 04_complementary/
│       └── 04_complementary_analyses.py   # GSEA, LR pairs, CD8 states, diversity, ratio biomarkers
│
├── slurm/                            # SLURM job scripts (SBATCH)
│   ├── 00_setup_env.sh                    # Conda env creation
│   ├── 01_convert_data.sh                 # .Rda/.rds → MTX conversion
│   ├── 02_preprocess_integrate.sh         # Integration (16 CPUs, 200 GB, 12 h)
│   ├── 03a_response_tme.sh               # TME characterisation
│   ├── 03b_spatial_biomarkers.sh          # Spatial analysis
│   ├── 03c_bulk_deconv.sh                # Deconvolution
│   ├── 04_validation.sh                  # Validation & downstream
│   ├── 05_complementary.sh               # Complementary analyses
│   └── run_full_pipeline.sh              # Master dependency-chain orchestrator
│
├── data/
│   ├── pre_treatment/                # Raw scRNA-seq per study
│   │   ├── GSE206325_Liu2022_antiPD1/     # Primary IO dataset (682 k cells)
│   │   ├── GSE125449_Ma2019_antiPD1/      # HCC landscape (~10 k cells)
│   │   ├── GSE140228_Zhang2019_immune/    # CD45+ immune (~67 k cells)
│   │   ├── GSE149614_Lu2022_landscape/    # HCC ecosystem (~70 k cells)
│   │   └── GSE202069_Zhu2023_antiPDL1/    # Bulk RNA-seq validation
│   ├── spatial_data/                 # Visium / pseudo-bulk
│   ├── clinical_trials/              # IMbrave150, CheckMate-459 supplements
│   ├── processed/                    # Converted MTX, metadata CSVs
│   │   └── integrated_io_atlas.h5ad       # Final atlas (~46 GB)
│   └── clinical_metadata/
│
├── out/                              # All pipeline outputs (figures + tables)
│   ├── 01_preprocessing/
│   ├── 02a_response_tme/
│   ├── 02b_spatial_biomarkers/
│   ├── 02c_bulk_deconvolution/
│   ├── 03_validation/
│   └── 04_complementary/
│
├── prompts/                          # Structured analysis prompts
│   └── analysis_prompts.md
├── config/                           # (reserved for future YAML configs)
├── environment.yaml                  # Conda env specification
├── Makefile                          # make all / download / preprocess / ...
├── LICENSE                           # MIT
└── README.md
```

---

## Pipeline Overview

```
 ┌──────────────────┐
 │ 00  Data Download │ GEO / ArrayExpress / clinical trials
 │     & Conversion  │ .Rda/.rds → MTX + CSV
 └────────┬─────────┘
          ▼
 ┌──────────────────┐
 │ 01  Preprocessing │ QC → Normalize → 3 000 HVGs → PCA (50)
 │     & Integration │ Harmony batch correction → UMAP → Leiden (35 clusters)
 └────────┬─────────┘  Output: integrated_io_atlas.h5ad (970 k cells × 15 853 genes)
          ▼
 ┌──────────────────────────────────────────────────────┐
 │                   02  Analysis                       │
 │  02a  TME Characterisation                           │
 │       • Marker-based cell type annotation            │
 │       • 16 IO signature scores                       │
 │       • Cell abundance R vs NR (Mann-Whitney)        │
 │       • Per-cell-type DEGs (Wilcoxon rank-sum)       │
 │       • Tissue-specific (Tumor / Normal) signatures  │
 │                                                      │
 │  02b  Spatial / Pseudo-bulk Biomarkers               │
 │       • GSE173905 pseudo-bulk IO signature scoring   │
 │       • Spatial immunotype classification             │
 │                                                      │
 │  02c  Bulk Deconvolution                             │
 │       • CIBERSORTx reference (15 853 genes × 60 CT) │
 │       • NNLS deconvolution of trial bulk RNA-seq     │
 │       • Inferred immunotype from proportions         │
 └────────┬─────────────────────────────────────────────┘
          ▼
 ┌──────────────────┐
 │ 03  Validation   │ • LASSO logistic (CV AUC 0.966)
 │     & Downstream │ • Random Forest  (CV AUC 0.941)
 │                  │ • NMF programs, TCR clonality
 │                  │ • Bulk validation (GSE202069)
 │                  │ • Bootstrap CI, pathway enrichment
 └────────┬─────────┘
          ▼
 ┌──────────────────┐
 │ 04  Complementary│ • 20 Hallmark pathways (GSEA-style)
 │     Analyses     │ • 24 Ligand-receptor pairs
 │                  │ • CD8 T cell states (6 states)
 │                  │ • Clonotype diversity (Shannon / Simpson)
 │                  │ • Immune ratio biomarkers (CD8/Treg, M1/M2)
 │                  │ • Patient immune subtyping (hierarchical)
 │                  │ • Cross-dataset reproducibility
 └──────────────────┘
```

---

## Datasets

### Single-Cell RNA-seq (integrated into atlas)

| Study | Accession | Context | Patients | Cells | Key Feature |
|-------|-----------|---------|----------|-------|-------------|
| Magen et al. 2023 _Nat Med_ | GSE206325 | Anti-PD-1 HCC | 30 | ~682 k | **R vs NR labels**, scRNA + TCR |
| Lu et al. 2022 _Nat Commun_ | GSE149614 | HCC ecosystem | — | ~70 k | Primary & metastatic |
| Ma et al. 2019 | GSE125449 | HCC landscape | 8 | ~10 k | Treatment-naïve |
| Zhang et al. 2019 _Cell_ | GSE140228 | CD45+ immune | — | ~67 k | Immune-enriched |
| Pfister et al. 2021 _Nature_ | E-MTAB-10553 | NASH-HCC + anti-PD-1 | 12 | ~152 k | NASH etiology |

### Bulk RNA-seq (validation & deconvolution)

| Accession | Study | Use |
|-----------|-------|-----|
| GSE202069 | Zhu et al. 2023 | Anti-PD-L1 bulk validation (IO signature scoring) |
| GSE173905 | Wu et al. 2021 | Pseudo-bulk spatial analysis (NC / ER / RS samples) |

### Clinical Trial Supplements

| Trial | Treatment | Note |
|-------|-----------|------|
| IMbrave150 | Atezolizumab + Bevacizumab | Supplementary XLSX (DEGs, metadata) |
| CheckMate-459 | Nivolumab | Supplementary XLSX |

---

## Key Results

### Response Prediction (Patient-Level)

| Model | Cross-Validated AUC | Top Features |
|-------|---------------------|--------------|
| LASSO Logistic Regression | **0.966** | T cell cytotoxicity ↑R, T cell exhaustion ↑R, Terminal exhaustion ↑R, Progenitor CD8 ↓R |
| Random Forest | **0.941** | T cell cytotoxicity, T cell exhaustion, cell type fractions |

### IO Signature Differences (R vs NR)

All 16 IO signatures tested; top differentiators (by effect size):
- **T cell cytotoxicity** — strongly elevated in Responders
- **T cell exhaustion / Terminal exhaustion** — elevated in Responders (paradoxical but consistent with active anti-tumour immunity)
- **T cell activation** — elevated in Responders
- **Progenitor CD8** — reduced in Responders
- **M1 / M2 macrophage** — reduced in Responders

### Bootstrap Stability

All top signatures have non-overlapping 95% confidence intervals across 50 bootstrap iterations; effect directions are robust.

---

## Outputs

Each pipeline step writes to `out/<step>/figures/` (PNG + PDF) and `out/<step>/tables/` (CSV).

### Figures (30+)

| Step | Figures |
|------|---------|
| `01_preprocessing` | QC violins, 4-panel integration UMAP (dataset / study / tissue / response), Leiden clusters |
| `02a_response_tme` | Cell type + response UMAP, signature violins (8 top), abundance barplot, tissue×signature heatmap, volcano DEGs, composition stacked bar, activation vs exhaustion scatter, signature radar |
| `02b_spatial` | Pseudo-bulk IO signature heatmap + barplot |
| `02c_deconvolution` | Reference heatmap, reference correlation, cell count barplot |
| `03_validation` | ROC + LASSO coefficients, RF importance, NMF elbow + gene heatmap + response boxplots, bulk clustermap, signature correlation, bootstrap CI, pathway enrichment |
| `04_complementary` | Hallmark pathways, LR interaction barplot, LR dotplot, CD8 state composition + scores, TCR diversity, immune ratios, patient subtype clustermap, cross-dataset heatmap, multi-panel summary |

### Tables (20+)

QC summary, cell type abundance, DEGs per cell type, signature comparisons, tissue-specific effects, patient features, LASSO/RF coefficients, NMF loadings, TCR clonality, bulk signatures, bootstrap effects, pathway enrichment, hallmark pathways, LR scores, CD8 state fractions, immune ratios, patient subtypes, cross-dataset profiles.

---

## Installation & Usage

### Prerequisites

- **Python 3.11** with: `scanpy ≥ 1.9`, `anndata ≥ 0.10`, `harmonypy`, `scikit-learn`, `scikit-misc`, `leidenalg`, `lifelines`, `decoupler`, `infercnvpy`
- **R 4.3+** with: `Seurat` (v5), `SingleCellExperiment`, `Matrix`
- **SLURM** cluster (or run scripts directly with sufficient RAM ≥ 200 GB for the integration step)

### Setup

```bash
# Clone
git clone https://github.com/isadeghi87/spatial_immunotherapy_response.git
cd spatial_immunotherapy_response

# Create conda environment
conda env create -f environment.yaml
conda activate spatial_io_response

# Install additional pip packages
pip install scikit-misc
```

### Run with SLURM (recommended)

```bash
# Full pipeline with dependency chain:
bash slurm/run_full_pipeline.sh

# Or submit individual steps:
sbatch slurm/01_convert_data.sh
sbatch slurm/02_preprocess_integrate.sh              # depends on 01
sbatch --dependency=afterok:<preprocess_jobid> slurm/03a_response_tme.sh
sbatch --dependency=afterok:<tme_jobid> slurm/03c_bulk_deconv.sh
sbatch --dependency=afterok:<tme_jobid> slurm/04_validation.sh
sbatch --dependency=afterok:<tme_jobid> slurm/05_complementary.sh
sbatch slurm/03b_spatial_biomarkers.sh                # runs independently
```

### Run with Make

```bash
make all          # Full pipeline: download → convert → preprocess → analyze → validate
make preprocess   # Just preprocessing & integration
make clean        # Remove processed data and outputs
```

### Run Scripts Directly

```bash
python scripts/01_preprocessing/01_preprocess_integrate.py
python scripts/02_analysis/02a_response_tme_characterization.py
python scripts/02_analysis/02b_spatial_biomarker_discovery.py
python scripts/02_analysis/02c_bulk_deconvolution.py
python scripts/03_validation/03_validation_downstream.py
python scripts/04_complementary/04_complementary_analyses.py
```

---

## IO Signatures Scored

| Signature | Genes |
|-----------|-------|
| T cell exhaustion | PDCD1, CTLA4, HAVCR2, LAG3, TIGIT, TOX, ENTPD1 |
| T cell activation | IFNG, GZMB, PRF1, TNF, CD69, IL2 |
| T cell cytotoxicity | GZMA, GZMB, GZMH, GZMK, PRF1, GNLY, NKG7 |
| Progenitor CD8 | TCF7, SLAMF6, LEF1, CCR7, IL7R, SELL |
| Terminal exhaustion | TOX, ENTPD1, HAVCR2, LAYN, TIGIT, CXCL13 |
| Treg signature | FOXP3, IL2RA, CTLA4, IKZF2, TNFRSF18 |
| M1 macrophage | NOS2, IL1B, IL6, TNF, CD80, CD86 |
| M2 macrophage | CD163, MRC1, MSR1, CD200R1, IL10 |
| IFN-γ response | IFNG, STAT1, IRF1, CXCL9, CXCL10, CXCL11, IDO1 |
| TGF-β signalling | TGFB1, TGFB2, TGFB3, ACTA2, COL1A1, SERPINE1 |
| Antigen presentation | HLA-A, HLA-B, HLA-C, B2M, TAP1, TAP2, PSMB9 |
| Immune checkpoint | CD274, PDCD1LG2, PDCD1, CTLA4, HAVCR2, LAG3, TIGIT, VSIR |
| Angiogenesis | VEGFA, PECAM1, KDR, FLT1, ANGPT2, FGF2 |
| WNT / β-catenin | CTNNB1, APC, AXIN2, LGR5, MYC |
| mregDC | CCR7, LAMP3, CCL19, CCL22, CD80, CD83, FSCN1, MARCKSL1 |
| Tfh-like | CXCL13, CH25H, IL21, PDCD1, ICOS, BCL6, TOX2 |

---

## Citation

```bibtex
@article{sadeghi2026spatial,
  title   = {Spatial Biomarkers of Immunotherapy Response in Hepatocellular Carcinoma},
  author  = {Sadeghi, Iman},
  year    = {2026},
  note    = {Manuscript in preparation}
}
```

## License

MIT License — see [LICENSE](LICENSE).

## Contact

**Iman Sadeghi** — [GitHub](https://github.com/isadeghi87)


## Summary of Findings & Future Directions

### Current Findings

1. **Integrated TME Atlas Response Maps**: We successfully harmonized ~970k cells across major HCC immunotherapy cohorts. Our cell type abundance tests and DEG analysis revealed that Responders exhibit distinct T cell state shifts (e.g., increased cytotoxicity and late exhaustion markers paired paradoxically with high activation), while Non-Responders lean towards progenitor CD8 and M1/M2 macrophage dominancies.
2. **Spatial & Pseudo-bulk Vulnerabilities**: Analysis of the spatial/pseudo-bulk data (GSE173905) mapped key IO signature scores (T cell exhaustion, Treg signatures) to distinct clinical trajectories like *Early Recurrence* vs *Recurrence Sensitive*. 
3. **CIBERSORTx Deconvolution Setup**: A robust high-variance reference signature matrix (15,853 genes $\times$ 60 distinct immune/stromal subsets) was successfully generated, allowing us to deconvolute larger trial bulk RNA-seq cohorts.

### Advanced Implemented Elements
Based on recent analysis expansions, the pipeline now fully incorporates these critical dimensions:
- **True Spatial Transcriptomics (Visium)**: Full integration of structural architectures, including tertiary lymphoid structures (TLS) scoring via `04g_tls_spatial_hubs.py`, mapping spatial topology gradients (`02e_spatial_topology.py`) utilizing spatial deconvolution proxies (Tangram/Cell2location workflows) at the tumor-stroma boundaries.
- **Bulk Cohort Predictive Validation**: Implementation of CIBERSORTx signature deconvolution applied to CheckMate-459 and IMbrave150 data (`03b_clinical_survival.py`), validating R vs NR biomarkers via rigorous Kaplan-Meier estimation and time-dependent ROC metrics.
- **TCR/Clonality Dynamics**: Seamless mapping of paired scTCR-seq V(D)J repertoires tracking clonal expansion frequencies and lineage commitment in distinct response groups (`04e_tcr_dynamics.py`).
- **Cell-Cell Communication**: A sophisticated pipeline mapping inter-cellular cross-talk bridging interactions utilizing a CellPhoneDB topology inference methodology, explicitly capturing mregDC orchestration of CD8+ T cell exhaustion along spatially mapped boundaries (`04d_cell_communication.py`).

## Materials and Methods / Documentation

The fully referenced **Materials and Methods** for this project, suitable for manuscript inclusion (*Nature Cancer* submission draft), are available in the `docs/` directory:
- [Materials and Methods (Markdown)](docs/materials_and_methods.md)
- [Materials and Methods (Word / DOCX)](docs/materials_and_methods.docx)

These documents detail the study cohorts, scRNA-seq preprocessing, TME characterization, bulk RNA-seq deconvolution, and advanced T-cell trajectory, spatial, and regulatory network analyses used throughout these scripts.
