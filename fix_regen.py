import pandas as pd

clust = pd.read_csv("data/processed/GSE206325_cluster_annotations.csv")
clust["cluster_ID"] = clust["cluster_ID"].astype(str)

major_map = {}
sub_map = {}

for _, row in clust.iterrows():
    cid = row["cluster_ID"]
    comp = row["compartment"]
    group = row["group"]
    type_ = row["type"] if pd.notna(row["type"]) else ""
    subg = row["subgroup"] if pd.notna(row["subgroup"]) else ""
    
    # Major immune
    if comp in ["Lymphoid", "Myeloid"]:
        major_map[cid] = group
    else:
        major_map[cid] = "Non-Immune"
        
    # Subtype
    if comp in ["Lymphoid", "Myeloid"]:
        if pd.notna(type_) and pd.notna(subg) and type_ != "" and subg != "":
            sub_map[cid] = f"{type_}_{subg}"
        elif pd.notna(type_) and type_ != "":
            sub_map[cid] = type_
        elif pd.notna(subg) and subg != "":
            sub_map[cid] = subg
        else:
            sub_map[cid] = group
    else:
        sub_map[cid] = "Non-Immune"

import json
print("major_map = " + json.dumps(major_map, indent=4))
