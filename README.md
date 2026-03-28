# Spatial Biomarkers of Immunotherapy Response in HCC

![Python](https://img.shields.io/badge/Python-3.10%2B-blue) ![R](https://img.shields.io/badge/R-4.3%2B-276DC3) ![License: MIT](https://img.shields.io/badge/License-MIT-green) ![Status](https://img.shields.io/badge/Status-Active-brightgreen)

> Defining spatial immune topology biomarkers that predict response vs. resistance to anti-PD1/anti-PD-L1 immunotherapy in hepatocellular carcinoma.

## Highlights

- Largest assembly of HCC immunotherapy scRNA-seq datasets
- Spatial TME changes upon anti-PD1 treatment
- Spatial architecture biomarkers (not just abundance)
- Composite spatial response score via LASSO
- Head-to-head comparison with TIDE, TMB, PD-L1, GEP

**Target Journal:** Nature Medicine / Nature Cancer / JCO

## Repository Structure

```
spatial_immunotherapy_response/
├── config/                          # Pipeline configuration
├── data/
│   ├── clinical_metadata/           # Patient response data
│   ├── post_treatment/              # Post-IO samples
│   └── pre_treatment/               # Pre-IO baseline samples
├── docs/                            # Documentation
├── envs/                            # Conda environment specs
├── notebooks/                       # Analysis notebooks
├── out/
│   ├── analysis/                    # Integrated results
│   ├── biomarkers/                  # Discovered biomarkers
│   └── figures/                     # Publication figures
├── prompts/                         # Analysis prompt templates
├── scripts/
│   ├── 00_data_acquisition/
│   │   ├── download_io_datasets.py
│   │   └── download_direct.sh
│   └── 02_analysis/
│       └── spatial_biomarkers.py
├── .gitignore
├── LICENSE
├── environment.yaml
├── Makefile
└── README.md
```

## Datasets

### HCC Immunotherapy scRNA-seq
| Study | GEO | Treatment | Patients | Cells |
|-------|-----|-----------|----------|-------|
| Liu et al. 2022 | GSE206325 | Anti-PD1 | 30 | ~120k |
| Zhu et al. 2023 | GSE202069 | Anti-PD1/PD-L1 | 14 | ~80k |
| Xue et al. 2023 | GSE234241 | Anti-PD1 + TKI | 20 | ~150k |
| Ma et al. 2019 | GSE125449 | Anti-PD1 | 8 | ~25k |
| Pfister et al. 2021 | E-MTAB-10553 | NASH-HCC + anti-PD1 | 12 | ~60k |

### Bulk RNA-seq with IO
| Source | Treatment | Patients |
|--------|-----------|----------|
| IMbrave150 | Atezo+Bev | 336 |
| CheckMate-459 | Nivo | 371 |
| KEYNOTE-224 | Pembro | 104 |

## Getting Started

### Prerequisites

- **Python 3.10+**: scanpy, scvi-tools, cell2location, squidpy, lifelines
- **R 4.3+**: Seurat v5, survival, glmnet, ComplexHeatmap, SPATA2

### Installation

```bash
# Clone the repository
git clone https://github.com/isadeghi87/spatial_immunotherapy_response.git
cd spatial_immunotherapy_response

# Create conda environment (if environment.yaml exists)
conda env create -f environment.yaml
# OR use envs/ specs
conda env create -f envs/environment.yaml

# Activate
conda activate spatial_immunotherapy_response
```

### Quick Start

```bash
# 1. Download data

# 2. Run analysis (if Makefile present)
make all
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Commit changes (`git commit -am 'Add feature'`)
4. Push to the branch (`git push origin feature/your-feature`)
5. Open a Pull Request

## Citation

If you use this work, please cite:

```bibtex
@article{sadeghi2026spatial_immunotherapy_response,
  title={Spatial Biomarkers of Immunotherapy Response in HCC},
  author={Sadeghi, Iman},
  year={2026},
  note={Manuscript in preparation}
}
```

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

## Contact

- **Iman Sadeghi** — [GitHub](https://github.com/isadeghi87)
