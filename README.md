# Spatial Biomarkers of Immunotherapy Response in HCC

![Python](https://img.shields.io/badge/Python-3.11-blue) ![R](https://img.shields.io/badge/R-4.3%2B-276DC3) ![License: MIT](https://img.shields.io/badge/License-MIT-green) ![Status](https://img.shields.io/badge/Status-Active-brightgreen)

> Defining spatial immune topology biomarkers that predict response vs. resistance to anti-PD1/anti-PD-L1 immunotherapy in hepatocellular carcinoma.

## Highlights

- Largest assembly of HCC immunotherapy scRNA-seq datasets (~1M cells integrated).
- TME characterization evaluating cell abundance, DEGs, and IO signature changes in Responder vs. Non-Responder.
- LASSO and Random Forest models predicting response.
- Hallmark pathway enrichment (GSEA-style), Ligand-Receptor communication, and CD8+ T cell states analysis.
- Deconvolution of bulk RNA-seq trials and patient subtyping.
- Fully automated SLURM-based pipeline.

## Pipeline Architecture

The analysis is broken down into modular steps running on SLURM:

1. **Data Acquisition & Conversion (`00_`):** Downloads and converts scRNA-seq datasets to `.h5ad` (e.g. Seurat v5 to AnnData).
2. **Preprocessing & Integration (`01_`):** QC, total-count normalization, log1p, highly variable genes, PCA, and `Harmony` batch integration across 5 scRNA-seq datasets. Produces the `integrated_io_atlas.h5ad` (~46GB).
3. **Response & TME Characterization (`02a_`):** Cell type annotation, IO signature scoring, cell abundance R vs NR (Mann-Whitney), and DEG discovery (Wilcoxon).
4. **Spatial Biomarkers (`02b_`):** Heatmaps and barplots on spatial/pseudo-bulk signatures.
5. **Bulk Deconvolution (`02c_`):** Generates a CIBERSORTx signature matrix from the scRNA-seq atlas to deconvolute trials (e.g. IMbrave150, CheckMate-459). 
6. **Validation & Downstream (`03_`):** Trains LASSO & Random forest response predictors (patient-level). Uncovers transcriptional programs via NMF. Analyzes TCR clonality and bulk signature correlations.
7. **Complementary Analyses (`04_`):** GSEA pathway enrichment, Ligand-Receptor interactions, CD8 state trajectories, Immune cell ratio biomarkers (e.g., CD8/Treg), and patient immune subtyping.

## Datasets

### HCC Immunotherapy scRNA-seq
| Study | GEO | Treatment | Patients | Cells |
|-------|-----|-----------|----------|-------|
| Liu et al. 2022 | GSE206325 | Anti-PD1 | 30 | ~680k |
| Lu et al. 2022 | GSE149614 | Primary/Met | - | ~70k |
| Ma et al. 2019 | GSE125449 | Landscape | 8 | ~10k |
| Zhang et al. 2019 | GSE140228 | Immune | - | CD45+ |
| Pfister et al. 2021 | E-MTAB-10553 | NASH-HCC + anti-PD1 | 12 | ~60k |

### Bulk RNA-seq with IO
| Dataset | Source | Treatment | Features |
|---------|--------|-----------|----------|
| GSE202069 | Zhu et al. 2023 | Anti-PDL1 | Bulk Validation |
| GSE173905 | - | - | Pseudo-bulk |

## Getting Started

### Environment

Conda environment managed via: `/data/projects/HE_ST_machine_learning/env/spatial_io_response/`
Requires `scanpy`, `harmonypy`, `sklearn`, `scipy`, `scikit-misc` (for NMF).

### Running the Pipeline

You can run the full pipeline locally or via the provided SLURM scripts.

```bash
# To run via SLURM
cd slurm
sbatch 01_convert_data.sh
sbatch 02_preprocess_integrate.sh
sbatch 03a_response_tme.sh
sbatch 03b_spatial_biomarkers.sh
sbatch 03c_bulk_deconv.sh
sbatch 04_validation.sh
sbatch 05_complementary.sh
```

## Figures and Tables

All generated outputs are grouped by their corresponding steps inside `out/`:
- `out/01_preprocessing/`
- `out/02a_response_tme/`
- `out/02c_bulk_deconvolution/`
- `out/03_validation/`
- `out/04_complementary/`

Inside each, check the `figures/` and `tables/` directories for detailed results (PNG/PDF, CSV).

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
