import re

with open("scripts/02_analysis/02a_response_tme_characterization.py", "r") as f:
    text = f.read()

def repl(match):
    return f'"{match.group(1)}_{match.group(2)}":'

text = re.sub(r'"(CD8)_\d+-(.*?)":', repl, text)
text = re.sub(r'"(CD4)_\d+-(.*?)":', repl, text)

# NK cells
text = re.sub(r'"NK_NK":\s*".*?",\s*', '', text)
text = re.sub(r'"NK_CD16":\s*".*?",\s*', '', text)
text = re.sub(r'"NK_Proliferating":\s*".*?",\s*', '', text)

with open("scripts/02_analysis/02a_response_tme_characterization.py", "w") as f:
    f.write(text)

with open("scripts/02_analysis/02g_infercnv.py", "r") as f:
    text = f.read()
text = text.replace('"CD4_4-Effector-Memory"', '"CD4_Effector-Memory"')
with open("scripts/02_analysis/02g_infercnv.py", "w") as f:
    f.write(text)

