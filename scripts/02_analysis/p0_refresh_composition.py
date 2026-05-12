#!/usr/bin/env python3
"""P0.2 helper — regenerate per-spot composition tables from existing annotations.

Reads:
  out/02b_spatial_biomarkers/tables/gse238264_celltype_annotations.csv
  out/02b_spatial_biomarkers/tables/gse281759_celltype_annotations.csv

Writes:
  out/02b_spatial_biomarkers/tables/gse238264_celltype_composition.csv  (overwrites zero-spot file)
  out/02b_spatial_biomarkers/tables/gse281759_celltype_composition.csv
  out/02b_spatial_biomarkers/tables/spatial_unified_spot_metadata.csv  (joined cohort)
"""
from pathlib import Path
import pandas as pd

PROJECT = Path(__file__).resolve().parents[2]
TAB = PROJECT / "out" / "02b_spatial_biomarkers" / "tables"

# ---- GSE238264 ----
df1 = pd.read_csv(TAB / "gse238264_celltype_annotations.csv")
df1 = df1.rename(columns={"Unnamed: 0": "spot_id"})
comp1 = (df1.groupby(["sample", "response", "celltype"])
            .size().rename("n_spots").reset_index())
comp1.to_csv(TAB / "gse238264_celltype_composition.csv", index=False)
print(f"[GSE238264] composition rows: {len(comp1)}; non-zero: {(comp1['n_spots']>0).sum()}")

# ---- GSE281759 ----
df2 = pd.read_csv(TAB / "gse281759_celltype_annotations.csv")
comp2 = (df2.groupby(["sample_id", "group", "dominant_cell_type"])
            .size().rename("n_spots").reset_index())
comp2.to_csv(TAB / "gse281759_celltype_composition.csv", index=False)
print(f"[GSE281759] composition rows: {len(comp2)}; non-zero: {(comp2['n_spots']>0).sum()}")

# ---- Unified per-spot metadata for downstream Phase 1 / Phase 2 ----
# Map both schemas to the same minimal columns.
u1 = pd.DataFrame({
    "spot_id":           df1["spot_id"],
    "sample_id":         df1["sample"],
    "dataset":           "GSE238264",
    "response":          df1["response"],
    "tissue_type":       "Tumor",  # GSE238264 are all tumor sections
    "dominant_celltype": df1["celltype"],
    "spatial_immunotype": pd.NA,
})
# pull immunotype where available from existing tables
imm_path = TAB / "gse238264_immunotype_composition.csv"
if imm_path.exists():
    imm = pd.read_csv(imm_path)
    if {"sample","spot","immunotype"}.issubset(imm.columns):
        u1 = u1.merge(
            imm[["sample","spot","immunotype"]].rename(
                columns={"sample":"sample_id","spot":"spot_id","immunotype":"spatial_immunotype_imm"}),
            on=["sample_id","spot_id"], how="left"
        )
        u1["spatial_immunotype"] = u1["spatial_immunotype_imm"].fillna(u1["spatial_immunotype"])
        u1 = u1.drop(columns=["spatial_immunotype_imm"])

u2 = pd.DataFrame({
    "spot_id":           df2["sample_id"].astype(str) + "_" + df2["barcodes"].astype(str),
    "sample_id":         df2["sample_id"],
    "dataset":           "GSE281759",
    "response":          pd.NA,  # multi-region cohort has no IO response label
    "tissue_type":       df2["group"],
    "dominant_celltype": df2["dominant_cell_type"],
    "spatial_immunotype": df2["spatial_immunotype"],
})

unified = pd.concat([u1, u2], ignore_index=True)
unified.to_csv(TAB / "spatial_unified_spot_metadata.csv", index=False)
print(f"[UNIFIED] spots: {len(unified)} | columns: {list(unified.columns)}")
print(unified["dataset"].value_counts())
print(unified["dominant_celltype"].value_counts())
