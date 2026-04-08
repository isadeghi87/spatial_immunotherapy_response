import scanpy as sc
import numpy as np

print("Loading h5ad...")
adata = sc.read_h5ad("data/processed/integrated_io_atlas.h5ad")

def fix_ct(x):
    if not isinstance(x, str):
        return x
    if x.lower() == "cd8_6-memory":
        return "CD8_Memory"
    if x.startswith("NK_") or x == "NK":
        return "NK"
    return x

for col in ["cell_type", "cell_type_original", "leiden"]:
    if col in adata.obs.columns:
        print(f"Fixing {col}...")
        adata.obs[col] = adata.obs[col].apply(fix_ct).astype("category")

print("Saving h5ad...")
adata.write_h5ad("data/processed/integrated_io_atlas.h5ad")
print("Done.")
