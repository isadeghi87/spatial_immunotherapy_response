#!/bin/bash
# ==============================================================================
# Master script: Submit the full pipeline as a SLURM dependency chain
#
# Usage:  bash slurm/run_full_pipeline.sh [--skip-env] [--skip-convert]
#
# Pipeline order:
#   00_setup_env → 01_convert_data → 02_preprocess → 03a_tme + 03b_spatial (parallel)
#   → 03c_deconv (after 02) → 04_validation (after 03a, 03b, 03c)
# ==============================================================================
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"
mkdir -p slurm/logs

SKIP_ENV=false
SKIP_CONVERT=false
for arg in "$@"; do
    case "$arg" in
        --skip-env)     SKIP_ENV=true ;;
        --skip-convert) SKIP_CONVERT=true ;;
    esac
done

echo "=========================================="
echo " Spatial IO Response — Full SLURM Pipeline"
echo "=========================================="
echo "Project: $PROJECT_DIR"
echo "Date:    $(date)"
echo ""

# ── Step 0: Conda environment ────────────────────────────────────
if [[ "$SKIP_ENV" == false ]]; then
    JOB0=$(sbatch --parsable slurm/00_setup_env.sh)
    echo "[SUBMITTED] Step 0 — Setup env:           Job $JOB0"
    DEP_ENV="--dependency=afterok:$JOB0"
else
    echo "[SKIPPED]  Step 0 — Setup env"
    DEP_ENV=""
fi

# ── Step 1: Convert R data → MTX ─────────────────────────────────
if [[ "$SKIP_CONVERT" == false ]]; then
    JOB1=$(sbatch --parsable $DEP_ENV slurm/01_convert_data.sh)
    echo "[SUBMITTED] Step 1 — Convert data:         Job $JOB1"
    DEP_CONVERT="--dependency=afterok:$JOB1"
else
    echo "[SKIPPED]  Step 1 — Convert data"
    DEP_CONVERT="$DEP_ENV"
fi

# ── Step 2: Preprocess & Integrate ────────────────────────────────
JOB2=$(sbatch --parsable $DEP_CONVERT slurm/02_preprocess_integrate.sh)
echo "[SUBMITTED] Step 2 — Preprocess+Integrate: Job $JOB2"

# ── Step 3a: TME Characterization (depends on Step 2) ────────────
JOB3A=$(sbatch --parsable --dependency=afterok:$JOB2 slurm/03a_response_tme.sh)
echo "[SUBMITTED] Step 3a — TME Characterization: Job $JOB3A"

# ── Step 3b: Spatial Biomarkers (independent of Step 2) ───────────
JOB3B=$(sbatch --parsable $DEP_CONVERT slurm/03b_spatial_biomarkers.sh)
echo "[SUBMITTED] Step 3b — Spatial Biomarkers:   Job $JOB3B"

# ── Step 3c: Bulk Deconvolution (depends on Step 2) ──────────────
JOB3C=$(sbatch --parsable --dependency=afterok:$JOB2 slurm/03c_bulk_deconv.sh)
echo "[SUBMITTED] Step 3c — Bulk Deconvolution:   Job $JOB3C"

# ── Step 4: Validation (depends on all Step 3) ───────────────────
JOB4=$(sbatch --parsable --dependency=afterok:$JOB3A:$JOB3C slurm/04_validation.sh)
echo "[SUBMITTED] Step 4 — Validation:            Job $JOB4"

echo ""
echo "=========================================="
echo " Pipeline submitted! Monitor with:"
echo "   squeue -u $USER"
echo "   tail -f slurm/logs/*.out"
echo "=========================================="
echo ""
echo "Dependency chain:"
echo "  Step 0 (env) → Step 1 (convert) → Step 2 (preprocess)"
echo "                                   ↓"
echo "                     Step 3a (TME) + Step 3b (spatial) + Step 3c (deconv)"
echo "                                   ↓"
echo "                              Step 4 (validation)"
