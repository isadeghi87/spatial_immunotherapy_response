import re
with open("scripts/01_preprocessing/regen_umap.py", "r") as f:
    code = f.read()

colors_major_inject = """
if "Other" in adata_pl.obs["immune_major"].cat.categories:
    major_colors = dict(zip(adata_pl.obs["immune_major"].cat.categories, sc.pl.palettes.default_20))
    major_colors["Other"] = "#D3D3D3"
    adata_pl.uns["immune_major_colors"] = [major_colors[c] for c in adata_pl.obs["immune_major"].cat.categories]
"""
colors_sub_inject = """
if "Other" in adata_pl.obs["immune_subtype"].cat.categories:
    sub_colors = dict(zip(adata_pl.obs["immune_subtype"].cat.categories, sc.pl.palettes.default_102))
    sub_colors["Other"] = "#D3D3D3"
    adata_pl.uns["immune_subtype_colors"] = [sub_colors[c] for c in adata_pl.obs["immune_subtype"].cat.categories]
"""

# I really just want to replace standard sc.pl.umap calls with a customized palette or just grey out 'Other'. 
