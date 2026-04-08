with open("scripts/01_preprocessing/regen_umap.py", "r") as f:
    text = f.read()

import re
old = 'sc.pl.umap(adata_pl, color="immune_major", ax=axes[0], show=False, title="A. Major Immune Cell Types")\nsc.pl.umap(adata_pl, color="immune_subtype", ax=axes[1], show=False, title="B. Immune Subtypes")'

new = '''
adata_pl.obs["immune_major"] = adata_pl.obs["immune_major"].astype("category")
adata_pl.obs["immune_subtype"] = adata_pl.obs["immune_subtype"].astype("category")

if "Other" in adata_pl.obs["immune_major"].cat.categories:
    import matplotlib.pyplot as plt
    colors = plt.cm.tab20.colors
    c_list = [colors[i % len(colors)] for i in range(len(adata_pl.obs["immune_major"].cat.categories))]
    c_list[list(adata_pl.obs["immune_major"].cat.categories).index("Other")] = "#E0E0E0"
    adata_pl.uns["immune_major_colors"] = c_list

if "Other" in adata_pl.obs["immune_subtype"].cat.categories:
    import matplotlib.pyplot as plt
    colors = plt.cm.tab20.colors
    c_list = [colors[i % len(colors)] for i in range(len(adata_pl.obs["immune_subtype"].cat.categories))]
    c_list[list(adata_pl.obs["immune_subtype"].cat.categories).index("Other")] = "#E0E0E0"
    adata_pl.uns["immune_subtype_colors"] = c_list

sc.pl.umap(adata_pl, color="immune_major", ax=axes[0], show=False, title="A. Major Immune Cell Types")
sc.pl.umap(adata_pl, color="immune_subtype", ax=axes[1], show=False, title="B. Immune Subtypes")
'''

text = text.replace(old, new)
with open("scripts/01_preprocessing/regen_umap.py", "w") as f:
    f.write(text)
