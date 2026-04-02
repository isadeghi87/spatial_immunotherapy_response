#!/bin/bash
set -e

cd /mnt/projects/spatial_immunotherapy_response

echo "Downloading IMbrave150 supplementary data..."
cd data/clinical_trials/IMbrave150
for i in {1..7}; do
  url="https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-022-01872-6/MediaObjects/41591_2022_1872_MOESM${i}_ESM.xlsx"
  echo "Attempting: $url"
  wget -q --no-check-certificate "$url" || true
done
cd ../../../

echo "Downloading CheckMate-459 supplementary data..."
cd data/clinical_trials/CheckMate459
for i in {1..7}; do
  url="https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-021-00513-7/MediaObjects/41591_2021_513_MOESM${i}_ESM.xlsx"
  echo "Attempting: $url"
  wget -q --no-check-certificate "$url" || true
done
cd ../../../

echo "Download complete!"
