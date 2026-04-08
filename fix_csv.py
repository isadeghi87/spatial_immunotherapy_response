import pandas as pd
df = pd.read_csv("data/processed/GSE206325_cluster_annotations.csv")
df.loc[df["group"] == "NK", "type"] = "NK"
df.loc[df["group"] == "NK", "subgroup"] = float("nan")
df.to_csv("data/processed/GSE206325_cluster_annotations.csv", index=False)
