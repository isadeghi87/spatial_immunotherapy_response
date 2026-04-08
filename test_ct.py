import pandas as pd
import anndata as ad

adata = ad.read_h5ad("data/processed/integrated_io_atlas.h5ad", backed="r")
clust_ann = pd.read_csv("data/processed/GSE206325_cluster_annotations.csv", index_col=0)
clust_ann["cluster_ID"] = clust_ann["cluster_ID"].astype(str)

major_map = {}
sub_map = {}
for _, row in clust_ann.iterrows():
    cid = row["cluster_ID"]
    group = str(row["group"]) if pd.notna(row["group"]) else "Unknown"
    ct = str(row["type"]) if pd.notna(row["type"]) else ""
    sub = str(row["subgroup"]) if pd.notna(row["subgroup"]) else ""
    
    major_map[cid] = group
    
    if ct and sub:
        sub_map[cid] = f"{ct}_{sub}"
    elif ct:
        sub_map[cid] = ct
    elif sub:
        sub_map[cid] = sub
    else:
        sub_map[cid] = group

print("Unique majors:", set(major_map.values()))
print("Unique subs:", sorted(set(sub_map.values()))[:10])
