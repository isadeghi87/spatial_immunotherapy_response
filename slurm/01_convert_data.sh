#!/bin/bash
#SBATCH --job-name=convert_data
#SBATCH --partition=long
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=01:30:00
#SBATCH --output=slurm/logs/01_convert_%j.out
#SBATCH --error=slurm/logs/01_convert_%j.err

# ==============================================================================
# Step 1: Convert R objects (Rda/RDS) → MTX + CSV for Python
# ==============================================================================
set -euo pipefail

PROJECT_DIR="/home/sadeghii/projects/spatial_immunotherapy_response"
cd "$PROJECT_DIR"

echo "[$(date)] Step 1: Converting R data objects to MTX format"

R_CMD="/data/projects/HE_ST_machine_learning/env/spatial_io_response/bin/Rscript"

# Skip GSE206325 if already converted
if [[ -f "data/processed/GSE206325_matrix.mtx.gz" ]]; then
    echo "[$(date)] GSE206325 already converted — skipping."
else
    echo "[$(date)] Converting GSE206325 (SCE .Rda → MTX)..."
    $R_CMD scripts/00_data_acquisition/convert_rda_to_h5ad.R 2>&1
fi

# Skip E-MTAB-10553 if already converted
if [[ -f "data/processed/E-MTAB-10553/matrix.mtx.gz" ]]; then
    echo "[$(date)] E-MTAB-10553 already converted — skipping."
else
    echo "[$(date)] Converting E-MTAB-10553 (Seurat .rds → MTX)..."
    $R_CMD scripts/00_data_acquisition/00b_convert_emtab10553.R 2>&1
fi

echo "[$(date)] Step 1 complete."
echo "Check data/processed/ for output files."
ls -lh data/processed/
