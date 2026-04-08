import re

with open("scripts/01_preprocessing/regen_umap.py", "r") as f:
    text = f.read()

# I will use a simple regex replacement to remove the "digit-" pattern in keys inside the CELLTYPE_PALETTE
def repl(match):
    return f'"{match.group(1)}_{match.group(2)}":'

# "CD8_1-Dysfunctional-terminal":
text = re.sub(r'"(CD8)_\d+-(.*?)":', repl, text)
text = re.sub(r'"(CD4)_\d+-(.*?)":', repl, text)

# NK cells just NK
text = re.sub(r'"NK_NK":\s*".*?",\n', '', text)
text = re.sub(r'"NK_CD16":\s*".*?",\n', '', text)
text = re.sub(r'"NK_Proliferating":\s*".*?",\n', '', text)

with open("scripts/01_preprocessing/regen_umap.py", "w") as f:
    f.write(text)
