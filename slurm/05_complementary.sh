#!/bin/bash
#SBATCH --job-name=complementary
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=08:00:00
#SBATCH --output=slurm/logs/05_complementary_%j.out
#SBATCH --error=slurm/logs/05_complementary_%j.err

# ==============================================================================
# Step 5: Complementary Downstream Analyses
#   - Hallmark pathway scoring (GSEA-style)
#   - Ligand-receptor communication
#   - CD8 T cell state characterization
#   - Clonotype diversity analysis
#   - Immune cell ratio biomarkers
#   - Patient immune subtyping
#   - Cross-dataset reproducibility
#   - Multi-panel summary figure
# ==============================================================================
set -euo pipefail

PROJECT_DIR="/home/sadeghii/projects/spatial_immunotherapy_response"
cd "$PROJECT_DIR"

ENV_PY="/data/projects/HE_ST_machine_learning/env/spatial_io_response/bin/python"
export PYTHONUNBUFFERED=1

echo "[$(date)] Step 5: Complementary Analyses"
$ENV_PY scripts/04_complementary/04_complementary_analyses.py 2>&1

echo "[$(date)] Step 5 complete."
echo "Figures:"; ls out/04_complementary/figures/*.png 2>/dev/null | wc -l
echo "Tables:"; ls out/04_complementary/tables/*.csv 2>/dev/null | wc -l
