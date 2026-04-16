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
PANDOC = "/data/projects/spatial_multiomics/tools/envs/r-omics/bin/pandoc"
CSS_FILE = DOCS / "science_immunology.css"

BANNER = """
<div class='journal-banner'>
  <div class='left'>Science Immunology&nbsp;|&nbsp;Research Article</div>
  <div class='right'>HCC Immunotherapy · Integrated Atlas · 2026</div>
</div>
"""

TITLE_BLOCK_MAIN = """
<h1 class='title'>Genomic instability shapes spatial immune architecture and predicts checkpoint-immunotherapy response in hepatocellular carcinoma</h1>

<p class='authors'>
Iman Sadeghi<sup>1,*</sup>
</p>

<p class='affiliations'>
<sup>1</sup> Institute for Immunology, Spatial Multi-omics Group, Bengsch Laboratory, Germany.<br>
<sup>*</sup> Corresponding author. Email: iman.sadeghi@example.org
</p>

<p class='one-sentence-summary'><strong>One-sentence summary:</strong> An integrated ~970 000-cell single-cell and spatial atlas of hepatocellular carcinoma immunotherapy nominates clonal heterogeneity, a CCL5/NKG7 cytotoxic T-cell program, and composite spatial-immune topology as complementary predictors of ICB benefit.</p>
"""

TITLE_BLOCK_SUPP = """
<h1 class='title'>Supplementary Materials — Genomic instability shapes spatial immune architecture and predicts checkpoint-immunotherapy response in hepatocellular carcinoma</h1>

<p class='authors'>Iman Sadeghi<sup>1,*</sup></p>
<p class='affiliations'>
<sup>1</sup> Institute for Immunology, Spatial Multi-omics Group, Bengsch Laboratory, Germany.<br>
<sup>*</sup> Corresponding author. Email: iman.sadeghi@example.org
</p>
"""


def md_to_html(md_path: Path) -> str:
    """Pandoc markdown → fragment HTML (body only)."""
    out = subprocess.run(
        [PANDOC, str(md_path), "-f", "markdown+pipe_tables+smart",
         "-t", "html5", "--no-highlight"],
        capture_output=True, text=True, check=True)
    return out.stdout


def wrap_html(body_html: str, title_block: str, page_title: str) -> str:
    """Wrap a pandoc fragment into a full HTML document."""
    # Strip the first <h1> (markdown title) and the following "Running title:" paragraph
    body_html = re.sub(r"<h1[^>]*>.*?</h1>", "", body_html, count=1, flags=re.S)
    body_html = re.sub(
        r"<p[^>]*>\s*<strong>Running title:</strong>[^<]*</p>", "",
        body_html, count=1, flags=re.S)
    # Also strip any leading <hr /> that immediately follows
    body_html = re.sub(r"^\s*<hr\s*/?>", "", body_html, count=1)

    # strip out the abstract block from body into a separate shaded box
    # A naive but safe approach: pandoc produces <h2>Abstract</h2><p>...</p> pattern
    # We transform: replace first "<h2 id=\"abstract\">Abstract</h2>" with <div class='abstract'>
    abstract_html = ""
    if "<h2 id=\"abstract\">Abstract</h2>" in body_html:
        before, rest = body_html.split("<h2 id=\"abstract\">Abstract</h2>", 1)
        # Grab everything up to the next <h2
        after_split = rest.split("<h2", 1)
        abstract_paragraphs = after_split[0]
        body_html = before + "<h2" + after_split[1] if len(after_split) == 2 else before
        abstract_html = f"<div class='abstract'>{abstract_paragraphs}</div>"

    # Wrap references (after "## References" pandoc will emit <h2 id="references">)
    if "<h2 id=\"references\">" in body_html:
        body_html = body_html.replace(
            "<h2 id=\"references\">", "<div class='references'><h2 id=\"references\">")
        body_html = body_html + "</div>"
        # if another top-level heading comes after (e.g. Tables), move closer to it
        # Simple approach above keeps references as inline sub-section with smaller font.

    return f"""<!DOCTYPE html>
<html lang='en'><head><meta charset='utf-8'>
<title>{page_title}</title>
</head><body>
{BANNER}
{title_block}
{abstract_html}
{body_html}
</body></html>
"""


def build(md_name: str, out_name: str, title_block: str, page_title: str):
    md_path = DOCS / md_name
    body = md_to_html(md_path)
    html_str = wrap_html(body, title_block, page_title)
    html_out = DOCS / (out_name.replace(".pdf", ".html"))
    html_out.write_text(html_str)
    HTML(string=html_str, base_url=str(DOCS)).write_pdf(
        DOCS / out_name, stylesheets=[CSS(filename=str(CSS_FILE))])
    print(f"wrote {DOCS / out_name}  ({(DOCS / out_name).stat().st_size // 1024} kB)")


if __name__ == "__main__":
    build("manuscript.md", "manuscript_sci_immunol.pdf",
          TITLE_BLOCK_MAIN, "Sadeghi et al., Sci. Immunol. 2026 — Main")
    build("supplementary.md", "supplementary_sci_immunol.pdf",
          TITLE_BLOCK_SUPP, "Sadeghi et al., Sci. Immunol. 2026 — Supplementary")
