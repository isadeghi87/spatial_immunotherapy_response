#!/bin/bash
set -e

# ==============================================================================
# Download data for Cappuyns et al. 2025 (J Hepatol) / Cappuyns et al. 2023 (Nat Commun)
#
# Paper: Single-cell RNA sequencing-derived signatures define response patterns
#        to atezolizumab + bevacizumab in advanced hepatocellular carcinoma
# DOI:   10.1016/j.jhep.2024.12.016 | PMID 39709141 | PMC12086051
#
# Data sources:
#   1. 21 HCC gene signatures (Table S2, mmc2.xlsx) — from J Hepatol 2025 supplementary
#   2. scRNAseq discovery cohort (91,347 cells, 31 advanced HCC tumours)
#      - From Cappuyns et al. 2023 Nat Commun (PMID 38030622)
#      - Count data available from: https://lambrechtslab.sites.vib.be/en/aHCC
#      - Raw reads in EGA: EGAS00001007547 (controlled access)
#   3. In-house bulk RNAseq (71 patients, atezo+bev) — EGA: EGAD50000001227 (controlled)
#
# Run from the project root directory.
# ==============================================================================

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$PROJECT_ROOT"

# ── 1. Download 21 HCC gene signatures (Table S2) ─────────────────────────────
echo "=============================================================================="
echo "Step 1: Downloading Cappuyns et al. 2025 gene signatures (mmc2.xlsx)"
echo "=============================================================================="
mkdir -p data/references

SIGNATURES_FILE="data/references/cappuyns2025_mmc2_signatures.xlsx"
if [ -f "$SIGNATURES_FILE" ] && python3 -c "import openpyxl; openpyxl.load_workbook('$SIGNATURES_FILE')" 2>/dev/null; then
    echo "  Already downloaded: $SIGNATURES_FILE"
else
    echo "  Attempting download from PMC..."
    # PMC supplementary (may require browser due to JavaScript redirect)
    curl -L -f -o "$SIGNATURES_FILE" \
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC12086051/bin/mmc2.xlsx" 2>/dev/null \
    && python3 -c "import openpyxl; openpyxl.load_workbook('$SIGNATURES_FILE')" 2>/dev/null \
    || {
        echo ""
        echo "  ┌─────────────────────────────────────────────────────────────────┐"
        echo "  │  MANUAL DOWNLOAD REQUIRED                                       │"
        echo "  │                                                                 │"
        echo "  │  The PMC download requires a browser.                           │"
        echo "  │  Please download mmc2.xlsx manually:                            │"
        echo "  │                                                                 │"
        echo "  │  1. Open in browser:                                            │"
        echo "  │     https://doi.org/10.1016/j.jhep.2024.12.016                 │"
        echo "  │  2. Scroll to 'Supplementary data'                              │"
        echo "  │  3. Download 'Multimedia component 2' (mmc2.xlsx, ~28 KB)       │"
        echo "  │  4. Save to: data/references/cappuyns2025_mmc2_signatures.xlsx  │"
        echo "  │                                                                 │"
        echo "  │  Alternative PMC link:                                           │"
        echo "  │     https://pmc.ncbi.nlm.nih.gov/articles/PMC12086051/          │"
        echo "  │     -> 'Supplementary Materials' -> mmc2.xlsx                   │"
        echo "  └─────────────────────────────────────────────────────────────────┘"
        echo ""
        rm -f "$SIGNATURES_FILE"
    }
fi

# ── 2. Download scRNAseq discovery cohort (Lambrechts lab) ─────────────────────
echo ""
echo "=============================================================================="
echo "Step 2: scRNAseq discovery cohort (Cappuyns et al. 2023 Nat Commun)"
echo "=============================================================================="
CAPPUYNS_DIR="data/pre_treatment/Cappuyns2023_NatCommun_atezobev"
mkdir -p "$CAPPUYNS_DIR"

if ls "$CAPPUYNS_DIR"/*.h5ad "$CAPPUYNS_DIR"/*.rds "$CAPPUYNS_DIR"/*.mtx* 2>/dev/null | head -1 > /dev/null 2>&1; then
    echo "  Data already present in: $CAPPUYNS_DIR"
else
    echo ""
    echo "  ┌─────────────────────────────────────────────────────────────────┐"
    echo "  │  MANUAL DOWNLOAD REQUIRED — Lambrechts Lab Website              │"
    echo "  │                                                                 │"
    echo "  │  The scRNAseq count data (91,347 cells, 31 advanced HCC        │"
    echo "  │  tumours treated with atezo+bev) is available from:             │"
    echo "  │                                                                 │"
    echo "  │  https://lambrechtslab.sites.vib.be/en/aHCC                    │"
    echo "  │                                                                 │"
    echo "  │  Metadata includes:                                             │"
    echo "  │    - Timepoint: week 0 (pre-treatment), week 3, week 6          │"
    echo "  │    - Treatment: atezo+bev / anti-PD(L)1 (ICP) / TKI            │"
    echo "  │    - Response: Responders vs Non-responders                     │"
    echo "  │    - Note: sample IDs differ from manuscript (GDPR)             │"
    echo "  │                                                                 │"
    echo "  │  Download all count data files and save to:                     │"
    echo "  │    $CAPPUYNS_DIR/                                               │"
    echo "  │                                                                 │"
    echo "  │  Alternative (raw reads, controlled access):                    │"
    echo "  │    EGA: EGAS00001007547                                         │"
    echo "  └─────────────────────────────────────────────────────────────────┘"
    echo ""
fi

# ── 3. In-house bulk RNAseq (controlled access) ───────────────────────────────
echo ""
echo "=============================================================================="
echo "Step 3: In-house bulk RNAseq (EGA controlled access)"
echo "=============================================================================="
echo "  The in-house bulk RNAseq cohort (71 patients, atezo+bev-treated) is"
echo "  available under controlled access at EGA: EGAD50000001227"
echo "  Apply via the UZLeuven-VIB data access committee if needed."
echo "  This data is used for bulk-level validation only, not scRNA-seq integration."
echo ""

echo "=============================================================================="
echo "Cappuyns data download script completed."
echo "=============================================================================="
