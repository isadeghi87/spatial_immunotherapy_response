#!/bin/bash
# Direct download of IO-HCC scRNA-seq datasets
set -euo pipefail

DATA_DIR="$HOME/projects/spatial_immunotherapy_response/data/pre_treatment"
LOG="$HOME/projects/spatial_immunotherapy_response/data/download_log.txt"
echo "$(date): Starting IO-HCC downloads" > "$LOG"

geo_fetch() {
    local geo_id="$1"
    local study="$2"
    local out_dir="$DATA_DIR/${geo_id}_${study}"
    mkdir -p "$out_dir"

    if [ "$(find "$out_dir" -type f \( -name '*.gz' -o -name '*.h5' \) 2>/dev/null | wc -l)" -gt 0 ]; then
        echo "  [SKIP] $geo_id" | tee -a "$LOG"
        return
    fi

    echo "  [FETCH] $geo_id ($study)" | tee -a "$LOG"
    local prefix="${geo_id:0:$((${#geo_id}-3))}nnn"
    local base_url="https://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}/${geo_id}/suppl"
    local tmp_idx="/tmp/geo_idx_${geo_id}.html"
    
    wget -q --no-check-certificate "$base_url/" -O "$tmp_idx" 2>/dev/null || { echo "    [FAIL]"; return; }
    local files=$(grep -oP 'href="\K[^"]+' "$tmp_idx" | grep -v 'Parent\|geo/series\|hhs.gov' | grep -v '/$')
    
    for fname in $files; do
        [ -f "$out_dir/$fname" ] && continue
        echo "    [GET] $fname"
        wget -q --no-check-certificate -P "$out_dir" "${base_url}/${fname}" 2>&1 || true
    done
    rm -f "$tmp_idx"
    echo "    [OK] $(find "$out_dir" -type f | wc -l) files" | tee -a "$LOG"
}

echo "================================================================="
echo "Spatial Immunotherapy Response — Data Download"
echo "================================================================="

geo_fetch "GSE206325" "Liu2022_antiPD1"
geo_fetch "GSE202069" "Zhu2023_antiPDL1"
geo_fetch "GSE234241" "Xue2023_combo"
geo_fetch "GSE125449" "Ma2019_antiPD1"

echo "$(date): IO-HCC downloads complete" | tee -a "$LOG"
