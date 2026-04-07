#!/bin/bash
#SBATCH --job-name=regen_umap
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=04:00:00
#SBATCH --output=slurm/logs/01b_regen_umap_%j.out
#SBATCH --error=slurm/logs/01b_regen_umap_%j.err

# ==============================================================================
# Step 1b: Regenerate UMAP figures
#   - 4-panel integration UMAP (dataset, study, tissue, response)
#   - Cell-type annotation UMAP with KDE contours + abundance bar
#   - Per-cell-type marker feature plots (dark magma style)
# ==============================================================================
set -euo pipefail

PROJECT_DIR="/home/sadeghii/projects/spatial_immunotherapy_response"
cd "$PROJECT_DIR"

ENV_PY="/data/projects/HE_ST_machine_learning/env/spatial_io_response/bin/python"
export PYTHONUNBUFFERED=1

echo "[$(date)] Step 1b: Regenerate UMAP + cell-type + marker feature plots"
echo "  CPUs: $SLURM_CPUS_PER_TASK"
echo "  Memory: $SLURM_MEM_PER_NODE MB"

$ENV_PY scripts/01_preprocessing/regen_umap.py 2>&1

echo "[$(date)] Step 1b complete."
echo "Figures:"; ls -lh out/01_preprocessing/figures/*.png 2>/dev/null
