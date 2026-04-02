#!/bin/bash
#SBATCH --job-name=setup_env
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=slurm/logs/00_setup_env_%j.out
#SBATCH --error=slurm/logs/00_setup_env_%j.err

# ==============================================================================
# Step 0: Create conda environment with all required Python packages
# ==============================================================================
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"
mkdir -p slurm/logs

echo "[$(date)] Creating conda environment: spatial_io_response"

# Use mamba if available, else conda
CONDA_CMD="conda"
command -v mamba &>/dev/null && CONDA_CMD="mamba"

$CONDA_CMD env create -f environment.yaml -y 2>&1 || {
    echo "Env may already exist, trying update..."
    $CONDA_CMD env update -f environment.yaml --prune 2>&1
}

echo "[$(date)] Environment setup complete."
echo "Activate with:  conda activate spatial_io_response"
