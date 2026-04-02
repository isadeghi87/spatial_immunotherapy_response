#!/bin/bash
#SBATCH --job-name=tme_analysis
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=08:00:00
#SBATCH --output=slurm/logs/03a_tme_%j.out
#SBATCH --error=slurm/logs/03a_tme_%j.err

# ==============================================================================
# Step 3a: Response-Associated TME Characterization (Prompt 2)
#   - Cell type annotation & IO signature scoring
#   - Cell type abundance R vs NR
#   - DEGs per cell type
#   - Publication figures
# ==============================================================================
set -euo pipefail

PROJECT_DIR="/home/sadeghii/projects/spatial_immunotherapy_response"
cd "$PROJECT_DIR"

ENV_PY="/data/projects/HE_ST_machine_learning/env/spatial_io_response/bin/python"
export PYTHONUNBUFFERED=1

echo "[$(date)] Step 3a: Response TME Characterization"
$ENV_PY scripts/02_analysis/02a_response_tme_characterization.py 2>&1

echo "[$(date)] Step 3a complete."
echo "Figures:"; ls out/02a_response_tme/figures/*.png 2>/dev/null | wc -l
echo "Tables:"; ls out/02a_response_tme/tables/*.csv 2>/dev/null | wc -l
