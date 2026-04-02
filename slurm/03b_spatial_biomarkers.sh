#!/bin/bash
#SBATCH --job-name=spatial_bio
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=slurm/logs/03b_spatial_%j.out
#SBATCH --error=slurm/logs/03b_spatial_%j.err

# ==============================================================================
# Step 3b: Spatial Biomarker Discovery (Prompt 3)
#   - Spatial immunotyping on Visium data
#   - EIES, TLS, spatial features
#   - LASSO spatial response model
# ==============================================================================
set -euo pipefail

PROJECT_DIR="/home/sadeghii/projects/spatial_immunotherapy_response"
cd "$PROJECT_DIR"

ENV_PY="/data/projects/HE_ST_machine_learning/env/spatial_io_response/bin/python"

echo "[$(date)] Step 3b: Spatial Biomarker Discovery"
$ENV_PY scripts/02_analysis/02b_spatial_biomarker_discovery.py 2>&1

echo "[$(date)] Step 3b complete."
