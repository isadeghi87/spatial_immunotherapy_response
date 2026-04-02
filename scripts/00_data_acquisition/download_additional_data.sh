#!/bin/bash
set -e

# ==============================================================================
# Script to download missing datasets and additional Spatial Transcriptomics data
# Run this from the root of the project directory.
# ==============================================================================

# 1. Setup Data Directories
echo "Setting up directories..."
mkdir -p data/pre_treatment/E-MTAB-10553_Pfister2021
mkdir -p data/spatial_data/GSE173905_Wu2021_Visium
mkdir -p data/clinical_trials/IMbrave150
mkdir -p data/clinical_trials/CheckMate459

# 2. Download E-MTAB-10553 (NASH-HCC + anti-PD1) from ArrayExpress
echo "Downloading E-MTAB-10553 (Pfister et al. 2021)..."
cd data/pre_treatment/E-MTAB-10553_Pfister2021
# Download metadata and processed counts matrices
wget -q -c https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10553/E-MTAB-10553.sdrf.txt
wget -q -c https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10553/E-MTAB-10553.idf.txt
# Check if processed zip exists, download if it does
wget -q -c https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10553/E-MTAB-10553.processed.1.zip || echo "Processed zip not found or requires manual extraction"
if [ -f "E-MTAB-10553.processed.1.zip" ]; then
    unzip -o E-MTAB-10553.processed.1.zip
fi
cd ../../../

# 3. Download additional Spatial Transcriptomics (ST) data for HCC
# Using a widely referenced Visium HCC dataset (e.g., Wu et al. 2021 - GSE173905)
echo "Downloading Spatial Transcriptomics data (GSE173905)..."
cd data/spatial_data/GSE173905_Wu2021_Visium
wget -q -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173905/suppl/GSE173905_RAW.tar
if [ -f "GSE173905_RAW.tar" ]; then
    tar -xf GSE173905_RAW.tar
    echo "Extracted Visium spatial data."
fi
cd ../../../

# 4. Instructions/Placeholders for Clinical Trials Processed Data
# Since you don't have EGA/dbGaP access, you will need to rely on the supplementary data 
# published alongside the trial papers (e.g., Zhu et al. Nature Medicine 2022 for IMbrave150).
echo "Setting up placeholders for Clinical Trial data..."
echo "To acquire IMbrave150 and CheckMate-459 data without controlled access:" > data/clinical_trials/README_TRIAL_DATA.txt
echo "1. Go to the Nature Medicine publication for IMbrave150 (Zhu et al. 2022)." >> data/clinical_trials/README_TRIAL_DATA.txt
echo "2. Download Supplementary Data tables containing the normalized RNA-seq count matrices." >> data/clinical_trials/README_TRIAL_DATA.txt
echo "3. Save them as CSV/TSV in data/clinical_trials/IMbrave150/ and data/clinical_trials/CheckMate459/." >> data/clinical_trials/README_TRIAL_DATA.txt

echo "=============================================================================="
echo "Download script completed successfully."
echo "Please check data/clinical_trials/README_TRIAL_DATA.txt for final manual steps."
echo "=============================================================================="
