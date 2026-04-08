.PHONY: all download convert preprocess analyze_tme analyze_spatial analyze_deconv validate figures clean

all: download convert preprocess analyze_tme analyze_spatial analyze_deconv validate

# ── Step 0: Data Acquisition ──────────────────────────────────
download:
	bash scripts/00_data_acquisition/download_direct.sh
	bash scripts/00_data_acquisition/download_additional_data.sh

convert:
	Rscript scripts/00_data_acquisition/convert_rda_to_h5ad.R
	Rscript scripts/00_data_acquisition/00b_convert_emtab10553.R

# ── Step 1: Preprocessing & Integration ───────────────────────
preprocess:
	python3 scripts/01_preprocessing/01_preprocess_integrate.py

# ── Step 2a: Prompt 2 — Response TME Characterization ─────────
analyze_tme:
	python3 scripts/02_analysis/02a_response_tme_characterization.py

# ── Step 2b: Prompt 3 — Spatial Biomarker Discovery ──────────
analyze_spatial:
	python3 scripts/02_analysis/02b_spatial_biomarker_discovery.py

# ── Step 2c: Prompt 4 — Bulk IO Trial Deconvolution ──────────
analyze_deconv:
	python3 scripts/02_analysis/02c_bulk_deconvolution.py

# ── Step 3: Validation & Downstream ──────────────────────────
validate:
	python3 scripts/03_validation/03_validation_downstream.py

figures:
	@echo "Figures generated during analysis steps:"
	@echo "  out/01_preprocessing/figures/"
	@echo "  out/02a_response_tme/figures/"
	@echo "  out/02b_spatial_biomarkers/figures/"
	@echo "  out/02c_bulk_deconvolution/figures/"
	@echo "  out/03_validation/figures/"

clean:
	rm -rf data/processed/*.h5ad out/
