#!/bin/bash
#SBATCH --job-name=bulk_deconv
#SBATCH --partition=long
#SBATCH --qos=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=slurm/logs/03c_deconv_%j.out
#SBATCH --error=slurm/logs/03c_deconv_%j.err

# ==============================================================================
# Step 3c: Bulk IO Trial Deconvolution (Prompt 4)
#   - Generate CIBERSORTx reference from scRNA-seq
#   - NNLS deconvolution of bulk RNA-seq
#   - Infer spatial immunotypes
# ==============================================================================
set -euo pipefail

PROJECT_DIR="/home/sadeghii/projects/spatial_immunotherapy_response"
cd "$PROJECT_DIR"

ENV_PY="/data/projects/HE_ST_machine_learning/env/spatial_io_response/bin/python"
export PYTHONUNBUFFERED=1

echo "[$(date)] Step 3c: Bulk Deconvolution"
$ENV_PY scripts/02_analysis/02c_bulk_deconvolution.py 2>&1

echo "[$(date)] Step 3c complete."
