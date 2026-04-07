#!/bin/bash
#SBATCH --job-name=03d_infercnv
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm/logs/03d_infercnv_%j.out

source slurm/00_setup_env.sh
python scripts/02_analysis/02g_infercnv.py
