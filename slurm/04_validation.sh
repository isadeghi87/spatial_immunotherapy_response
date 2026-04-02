#!/bin/bash
#SBATCH --job-name=validation
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=08:00:00
#SBATCH --output=slurm/logs/04_validation_%j.out
#SBATCH --error=slurm/logs/04_validation_%j.err

# ==============================================================================
# Step 4: Validation & Downstream Analyses
#   - LASSO composite IO response score
#   - Random Forest comparison
#   - NMF transcriptional programs
#   - Signature correlation network
#   - Bulk RNA-seq validation
#   - TCR clonality
#   - Bootstrap stability
#   - Pathway enrichment
# ==============================================================================
set -euo pipefail

PROJECT_DIR="/home/sadeghii/projects/spatial_immunotherapy_response"
cd "$PROJECT_DIR"

ENV_PY="/data/projects/HE_ST_machine_learning/env/spatial_io_response/bin/python"
export PYTHONUNBUFFERED=1

echo "[$(date)] Step 4: Validation & Downstream"
$ENV_PY scripts/03_validation/03_validation_downstream.py 2>&1

echo "[$(date)] Step 4 complete."
echo "=== Final Output Summary ==="
for d in out/*/; do
    n_fig=$(ls "$d"figures/*.png 2>/dev/null | wc -l)
    n_tab=$(ls "$d"tables/*.csv 2>/dev/null | wc -l)
    echo "  $d: $n_fig figures, $n_tab tables"
done
