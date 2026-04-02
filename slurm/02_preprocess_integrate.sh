#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=slurm/logs/02_preprocess_%j.out
#SBATCH --error=slurm/logs/02_preprocess_%j.err

# ==============================================================================
# Step 2: Preprocess & Integrate all scRNA-seq datasets
#   - Loads ~800k+ cells across 5 datasets
#   - QC filtering, normalization, HVG selection
#   - Harmony batch integration
#   - UMAP, Leiden clustering
#   - Saves integrated_io_atlas.h5ad
# ==============================================================================
set -euo pipefail

PROJECT_DIR="/home/sadeghii/projects/spatial_immunotherapy_response"
cd "$PROJECT_DIR"

# Activate environment
ENV_PY="/data/projects/HE_ST_machine_learning/env/spatial_io_response/bin/python"
export PYTHONUNBUFFERED=1

echo "[$(date)] Step 2: Preprocessing & Integration"
echo "  CPUs: $SLURM_CPUS_PER_TASK"
echo "  Memory: $SLURM_MEM_PER_NODE MB"

$ENV_PY scripts/01_preprocessing/01_preprocess_integrate.py 2>&1

echo "[$(date)] Step 2 complete."
ls -lh data/processed/integrated_io_atlas.h5ad 2>/dev/null || echo "WARNING: integrated atlas not found"
