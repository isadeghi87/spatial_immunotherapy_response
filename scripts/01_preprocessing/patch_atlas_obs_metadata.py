"""
Patch integrated_io_atlas.h5ad: restore obs metadata columns (response,
patient_ID, tissue, disease_etiology, treatment) that were lost during
sc.concat(join="inner").

Run with:
  /data/apps/mambaforge/25.3.1/bin/python scripts/01_preprocessing/patch_atlas_obs_metadata.py
"""
from pathlib import Path
import gzip, pandas as pd, numpy as np, anndata as ad

PROJECT  = Path(__file__).resolve().parents[2]
PROC     = PROJECT / "data" / "processed"
DATA     = PROJECT / "data"
ATLAS    = PROC / "integrated_io_atlas.h5ad"

# ── Response vocabulary normalization ────────────────────────────────────────
_RESP_MAP_206325 = {
    "anti-PD1_R":  "Responder",
    "anti-PD1_NR": "Non-Responder",
}
_RESP_MAP_CAPP = {
    "Responder":          "Responder",
    "NonResponder":       "Non-Responder",
    "Non-Responder":      "Non-Responder",
    "DeathBeforeImaging": "Non-Responder",
}

# ── Build per-dataset obs DataFrames (metadata only, no count matrix) ─────────
print("=" * 70)
print("Atlas obs metadata patch")
print("=" * 70)

# 1. GSE206325 ────────────────────────────────────────────────────────────────
print("\n[GSE206325] Building metadata...")
cell_meta = pd.read_csv(PROC / "GSE206325_metadata.csv", index_col=0)
sample_ann = pd.read_csv(PROC / "GSE206325_sample_annotations.csv", index_col=0)
sample_ann.index = sample_ann.index.astype(str)

gse206325_meta = pd.DataFrame(index=cell_meta.index)
gse206325_meta["sample_id"] = cell_meta["cell_to_sample_ID"].astype(str)

for col in ["patient_ID", "tissue", "treatment_Resp"]:
    mapping = sample_ann[col].to_dict()
    gse206325_meta[col] = gse206325_meta["sample_id"].map(mapping)

gse206325_meta["response"] = gse206325_meta["treatment_Resp"].map(_RESP_MAP_206325)
gse206325_meta["patient_ID"] = gse206325_meta["patient_ID"].astype(str)
gse206325_meta["disease_etiology"] = "HCC"
print(f"  {len(gse206325_meta):,} cells  response: {gse206325_meta['response'].value_counts().to_dict()}")

# 2. GSE140228 ────────────────────────────────────────────────────────────────
print("\n[GSE140228] Building metadata...")
ci = pd.read_csv(
    DATA / "pre_treatment" / "GSE140228_Zhang2019_immune" / "GSE140228_UMI_counts_Droplet_cellinfo.tsv.gz",
    sep="\t", index_col=0)
gse140228_meta = pd.DataFrame(index=ci.index)
gse140228_meta["patient_ID"]        = ci["Donor"].astype(str)
gse140228_meta["disease_etiology"]  = ci["Histology"].fillna("Unknown").astype(str)
gse140228_meta["tissue"]            = ci["Tissue"].astype(str)
gse140228_meta["response"]          = "Unknown"   # no IO treatment in this dataset
print(f"  {len(gse140228_meta):,} cells")

# 3. GSE149614 ────────────────────────────────────────────────────────────────
print("\n[GSE149614] Building metadata...")
lu_path = DATA / "pre_treatment" / "GSE149614_Lu2022_landscape" / "GSE149614_HCC.metadata.updated.txt.gz"
lu = pd.read_csv(lu_path, sep="\t", index_col=0)
gse149614_meta = pd.DataFrame(index=lu.index)
gse149614_meta["patient_ID"]       = lu["patient"].astype(str)
gse149614_meta["disease_etiology"] = lu["virus"].fillna("Unknown").astype(str)
gse149614_meta["tissue"]           = lu.get("site", pd.Series("Unknown", index=lu.index)).astype(str)
gse149614_meta["response"]         = "Unknown"   # landscape study, no IO treatment
print(f"  {len(gse149614_meta):,} cells")

# 4. Cappuyns2023 ─────────────────────────────────────────────────────────────
print("\n[Cappuyns2023] Building metadata...")
capp_path = DATA / "pre_treatment" / "Cappuyns2023_NatCommun_atezobev" / "HCC_data_VIB_CCB_portal" / "meta_data.txt"
capp = pd.read_csv(capp_path, sep="\t", index_col=0)
capp_meta = pd.DataFrame(index=capp.index)
capp_meta["patient_ID"]       = capp["Patient_ID"].astype(str)
capp_meta["response"]         = capp["Response"].map(_RESP_MAP_CAPP).fillna("Unknown").astype(str)
capp_meta["treatment"]        = capp["Treatment"].astype(str)
capp_meta["tissue"]           = "Tumor"
capp_meta["disease_etiology"] = "HCC"
print(f"  {len(capp_meta):,} cells  response: {capp_meta['response'].value_counts().to_dict()}")

# ── Load atlas obs ────────────────────────────────────────────────────────────
print(f"\n[ATLAS] Loading obs from {ATLAS}...")
adata = ad.read_h5ad(ATLAS, backed="r")
obs = adata.obs.copy()
print(f"  {len(obs):,} cells  columns before: {list(obs.columns)}")

# ── Merge metadata ────────────────────────────────────────────────────────────
print("\n[MERGE] Attaching metadata per dataset...")

# Initialize columns with defaults
for col in ["response", "patient_ID", "tissue", "disease_etiology", "treatment"]:
    obs[col] = "Unknown"

_datasets = {
    "GSE206325":   (gse206325_meta,  ["response", "patient_ID", "tissue", "disease_etiology"]),
    "GSE140228":   (gse140228_meta,  ["response", "patient_ID", "disease_etiology", "tissue"]),
    "GSE149614":   (gse149614_meta,  ["response", "patient_ID", "disease_etiology", "tissue"]),
    "Cappuyns2023":(capp_meta,       ["response", "patient_ID", "tissue", "disease_etiology", "treatment"]),
}

for ds_name, (meta_df, cols) in _datasets.items():
    mask = obs["dataset"] == ds_name
    n = mask.sum()
    if n == 0:
        print(f"  {ds_name}: 0 cells in atlas -- skipping")
        continue
    ds_obs = obs.loc[mask]
    matched = ds_obs.index.intersection(meta_df.index)
    print(f"  {ds_name}: {n:,} cells, {len(matched):,} matched obs_names")
    for col in cols:
        if col in meta_df.columns:
            obs.loc[matched, col] = meta_df.loc[matched, col].values

    # Report match quality
    unmatched = n - len(matched)
    if unmatched > 0:
        print(f"    WARNING: {unmatched:,} cells had no metadata match")

# Make sure no NaN remains in string columns
for col in ["response", "patient_ID", "tissue", "disease_etiology", "treatment"]:
    obs[col] = obs[col].fillna("Unknown").astype(str)

print("\n[SUMMARY] After patch:")
print(f"  response:          {obs['response'].value_counts().to_dict()}")
print(f"  patient_ID unique: {obs['patient_ID'].nunique()}")
print(f"  tissue:            {obs['tissue'].value_counts().head(6).to_dict()}")
print(f"  disease_etiology:  {obs['disease_etiology'].value_counts().to_dict()}")

# ── Write patched atlas ───────────────────────────────────────────────────────
print(f"\n[SAVE] Loading full atlas into memory (49 GB, ~5 min)...")
adata.file.close()  # close backed mode
adata = ad.read_h5ad(ATLAS)
adata.obs = obs

# Sanitize all object columns
for col in adata.obs.select_dtypes(include="object").columns:
    adata.obs[col] = adata.obs[col].astype(str)

print(f"[SAVE] Writing patched atlas to {ATLAS}...")
adata.write_h5ad(ATLAS)
print(f"[DONE] Atlas patched: {adata.n_obs:,} cells, {len(adata.obs.columns)} obs columns")
print(f"       obs columns: {list(adata.obs.columns)}")
