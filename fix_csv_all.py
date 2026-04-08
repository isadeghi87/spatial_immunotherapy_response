import pandas as pd
import re

df = pd.read_csv("data/processed/GSE206325_cluster_annotations.csv")

def clean_name(s):
    if pd.isna(s):
        return s
    # Remove leading digits followed by dash (e.g. "7-Cytotoxic" -> "Cytotoxic")
    s = re.sub(r'^\d+-', '', str(s))
    return s

df["subgroup"] = df["subgroup"].apply(clean_name)

# Ensure all NK are completely simple:
df.loc[df["group"] == "NK", "type"] = float("nan")
df.loc[df["group"] == "NK", "subgroup"] = float("nan")

df.to_csv("data/processed/GSE206325_cluster_annotations.csv", index=False)

print(df["subgroup"].unique())
