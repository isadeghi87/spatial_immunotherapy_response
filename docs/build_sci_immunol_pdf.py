#!/usr/bin/env python3
"""Render docs/manuscript.md (+ docs/supplementary.md) to a
Science Immunology-styled PDF using pandoc + WeasyPrint.

Outputs:
    docs/manuscript_sci_immunol.pdf
    docs/supplementary_sci_immunol.pdf
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

from weasyprint import HTML, CSS

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
FIG_DIR = ROOT / "out" / "06_manuscript_figures"
PANDOC = "/data/projects/spatial_multiomics/tools/envs/r-omics/bin/pandoc"
CSS_FILE = DOCS / "science_immunology.css"

TITLE_BLOCK_MAIN = """
<h1 class='title'>Genomic instability shapes spatial immune architecture and predicts checkpoint-immunotherapy response in hepatocellular carcinoma</h1>

<p class='authors'>
Iman Sadeghi<sup>1</sup>, Bertram Bengsch<sup>1</sup>
</p>

<p class='affiliations'>
<sup>1</sup> Clinic for Internal Medicine II, Freiburg University Medical Center, Freiburg im Breisgau, Germany.
</p>

<p class='one-sentence-summary'><strong>One-sentence summary:</strong> An integrated ~970 000-cell single-cell and spatial atlas of hepatocellular carcinoma immunotherapy nominates clonal heterogeneity, a CCL5/NKG7 cytotoxic T-cell program, and composite spatial-immune topology as complementary predictors of ICB benefit.</p>
"""

TITLE_BLOCK_SUPP = """
<h1 class='title'>Supplementary Materials — Genomic instability shapes spatial immune architecture and predicts checkpoint-immunotherapy response in hepatocellular carcinoma</h1>

<p class='authors'>Iman Sadeghi<sup>1</sup>, Bertram Bengsch<sup>1</sup></p>
<p class='affiliations'>
<sup>1</sup> Clinic for Internal Medicine II, Freiburg University Medical Center, Freiburg im Breisgau, Germany.
</p>
"""

# Main-manuscript figures to embed (in order, before References).
MAIN_FIGURES = [
    ("Figure1_atlas_landscape.png",
     "Fig. 1. Integrated single-cell atlas of HCC immunotherapy.",
     "UMAP of 970,632 cells coloured by dataset, cell type, tissue and response; cohort composition; batch-effect silhouette before / after Harmony."),
    ("Figure2_tumour_intrinsic.png",
     "Fig. 2. Tumour-intrinsic axis of response (CNV + CD8 states).",
     "Per-patient CNV mean vs. standard deviation by response; CD8 state composition split by R/NR; state-level effect sizes with BH-padj."),
    ("Figure3_immune_engagement.png",
     "Fig. 3. Immune engagement and hallmark pathways.",
     "20-pathway R-vs-NR effect heatmap; top-10 ligand\u2013receptor dotplot across sender\u2013receiver pairs; CD8 activation\u2013exhaustion scatter."),
    ("Figure4_spatial_axis.png",
     "Fig. 4. Spatial axis: NMF9, TLS and Visium projection.",
     "Visium sections (GSE238264) coloured by NMF9 and TLS score; cross-cohort signature robustness across Visium cohorts."),
    ("Figure5_nmf_integrator.png",
     "Fig. 5. NMF integrator across scRNA / bulk / Visium.",
     "Program \u00d7 patient activity heatmap with response annotation; NMF9 top-gene loadings and volcano of response stats; projected activity in bulk and Visium."),
    ("Figure6_multilayer_predictor.png",
     "Fig. 6. Five-layer predictor.",
     "Per-layer LOOCV ROC curves (L1\u2013L5 + ensemble); LASSO-selected features and variance partition; balanced cumulative AUC as layers are added."),
]


def _figures_block(figs, prefix=""):
    """Build an HTML block of <figure> elements spanning both columns."""
    items = []
    for i, (fname, title, caption) in enumerate(figs, start=1):
        fpath = FIG_DIR / fname
        if not fpath.exists():
            continue
        items.append(
            f"<figure class='manuscript-figure'>"
            f"<img src='file://{fpath}' alt='{title}' />"
            f"<figcaption><strong>{title}</strong> {caption}</figcaption>"
            f"</figure>"
        )
    return "<section class='figures-section'>" + "\n".join(items) + "</section>"


def md_to_html(md_path: Path) -> str:
    """Pandoc markdown → fragment HTML (body only)."""
    out = subprocess.run(
        [PANDOC, str(md_path), "-f", "markdown+pipe_tables+smart",
         "-t", "html5", "--no-highlight"],
        capture_output=True, text=True, check=True)
    return out.stdout


def wrap_html(body_html: str, title_block: str, page_title: str,
              figures_html: str = "") -> str:
    """Wrap a pandoc fragment into a full HTML document."""
    # Strip the first <h1> (markdown title) and the following "Running title:" paragraph
    body_html = re.sub(r"<h1[^>]*>.*?</h1>", "", body_html, count=1, flags=re.S)
    body_html = re.sub(
        r"<p[^>]*>\s*<strong>Running title:</strong>[^<]*</p>", "",
        body_html, count=1, flags=re.S)
    # Also strip any leading <hr /> that immediately follows
    body_html = re.sub(r"^\s*<hr\s*/?>", "", body_html, count=1)

    # Extract the abstract section into a shaded box
    abstract_html = ""
    if "<h2 id=\"abstract\">Abstract</h2>" in body_html:
        before, rest = body_html.split("<h2 id=\"abstract\">Abstract</h2>", 1)
        after_split = rest.split("<h2", 1)
        abstract_paragraphs = after_split[0]
        body_html = before + "<h2" + after_split[1] if len(after_split) == 2 else before
        abstract_html = f"<div class='abstract'>{abstract_paragraphs}</div>"

    # Strip the markdown "Figure legends" section — figures are embedded
    # directly below their referring content via figures_html.
    body_html = re.sub(
        r"<h2 id=\"figure-legends\">.*?(?=<h2|\Z)", "",
        body_html, count=1, flags=re.S)

    # Split body at the References heading so we can render the figures
    # section BETWEEN the main text (2-col) and the references (2-col).
    ref_marker = "<h2 id=\"references\">"
    if ref_marker in body_html:
        main_body, refs_body = body_html.split(ref_marker, 1)
        refs_html = ref_marker + refs_body
    else:
        main_body, refs_html = body_html, ""

    main_col = f"<div class='two-col'>{main_body}</div>"
    refs_col = (f"<div class='references two-col'>{refs_html}</div>"
                if refs_html else "")

    return f"""<!DOCTYPE html>
<html lang='en'><head><meta charset='utf-8'>
<title>{page_title}</title>
</head><body>
{title_block}
{abstract_html}
{main_col}
{figures_html}
{refs_col}
</body></html>
"""


def build(md_name: str, out_name: str, title_block: str, page_title: str,
          figures_html: str = ""):
    md_path = DOCS / md_name
    body = md_to_html(md_path)
    html_str = wrap_html(body, title_block, page_title, figures_html)
    html_out = DOCS / (out_name.replace(".pdf", ".html"))
    html_out.write_text(html_str)
    HTML(string=html_str, base_url=str(DOCS)).write_pdf(
        DOCS / out_name, stylesheets=[CSS(filename=str(CSS_FILE))])
    print(f"wrote {DOCS / out_name}  ({(DOCS / out_name).stat().st_size // 1024} kB)")


if __name__ == "__main__":
    build("manuscript.md", "manuscript_sci_immunol.pdf",
          TITLE_BLOCK_MAIN, "Sadeghi & Bengsch — Main",
          figures_html=_figures_block(MAIN_FIGURES))
    build("supplementary.md", "supplementary_sci_immunol.pdf",
          TITLE_BLOCK_SUPP, "Sadeghi & Bengsch — Supplementary")
