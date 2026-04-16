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
# Legends are panel-resolved and quote effect sizes from the manuscript.
MAIN_FIGURES = [
    ("Figure1_atlas_landscape.png",
     "Fig. 1. Integrated single-cell atlas of hepatocellular carcinoma immunotherapy.",
     "(A) UMAP of 970,632 single cells from five harmonised public cohorts "
     "(GSE206325, GSE125449, GSE140228, GSE149614, E-MTAB-10553) coloured by "
     "manual cell-type annotation; 16 major populations spanning malignant "
     "hepatocytes, CD4\u2009/\u2009CD8 T, NK, B\u2009/\u2009plasma, myeloid, endothelial and "
     "stromal compartments are resolved. "
     "(B) Same embedding coloured by dataset of origin after Harmony "
     "integration (donor and study as batch covariates); overlapping "
     "distributions indicate that biological structure, rather than cohort, "
     "dominates the embedding. "
     "(C) Dot plot of canonical lineage markers (rows) across annotated "
     "cell types (columns); dot size encodes percent of cells expressing the "
     "gene, colour encodes scaled mean expression. "
     "(D) Stacked cell-type composition per patient, grouped by "
     "anti-PD-1 response label (7 Responders, 17 Non-Responders); bar "
     "heights are normalised to 100 % so compositional shifts between R "
     "and NR are directly comparable."),

    ("Figure2_tumour_intrinsic.png",
     "Fig. 2. Tumour-intrinsic axis of response: CNV burden, clonal heterogeneity and TME context.",
     "(A) UMAP coloured by per-cell InferCNV score (malignancy proxy); "
     "hepatocyte clusters display the expected elevated CNV score, "
     "non-malignant compartments remain near baseline. "
     "(B) Per-patient chromosome-level CNV heatmap (rows, patients; "
     "columns, genomic bins) with recurrent gains (1q, 8q) and losses "
     "(4q, 8p, 17p) visible in the malignant compartment, alongside "
     "sample-level CNV-score violins. "
     "(C) InferCNV dashboard: per-patient CNV-score distributions, "
     "malignant-fraction and clone-diversity (Shannon) summaries. "
     "(D) Patient-level violins comparing Responders vs Non-Responders "
     "for CNV mean (n.s., p\u2009=\u20090.42, d\u2009=\u2009\u22120.48), CNV standard "
     "deviation (clonal heterogeneity, *p\u2009=\u20090.028, d\u2009=\u2009\u22120.62, "
     "Mann\u2013Whitney U), malignant fraction (n.s.) and CNV clone diversity "
     "(Shannon, n.s.); dots, patients; box, median and IQR. "
     "(E) Tumour-microenvironment context at high CNV clonality: "
     "immune-ecosystem composition and bulk-RNA-seq concordance panels "
     "showing that high clonal heterogeneity co-occurs with depleted "
     "lymphoid and enriched myeloid signatures."),

    ("Figure3_immune_engagement.png",
     "Fig. 3. Immune engagement axis: composition, signatures and hallmark pathways.",
     "(A) UMAP of immune cells coloured by response label (R, blue; "
     "NR, red) showing Responder enrichment in CD8 effector / cytotoxic "
     "islands. "
     "(B) Per-cell-type abundance bar plot (mean \u00b1 SE of per-patient "
     "proportion) with BH-adjusted Mann\u2013Whitney p-values annotated; "
     "Responders show significantly higher B/Plasma and dysfunctional-"
     "effector CD8 populations. "
     "(C) Radar plot of mean scaled scores for 10 IO signatures "
     "(cytotoxicity, exhaustion, interferon-\u03b3, TLS, antigen presentation, "
     "TGF-\u03b2, angiogenesis, EMT, stemness, hepatocyte-differentiation) "
     "per response group. "
     "(D) Violin plots of key signatures (cytotoxicity, exhaustion, "
     "IFN-\u03b3, TLS) by response; Responders show co-elevated cytotoxicity "
     "and exhaustion\u2014the \u201cexhaustion\u2013cytotoxicity paradox\u201d. "
     "(E) Tissue-level heatmap of 20 hallmark and IO signatures across "
     "R and NR samples (z-scored); interferon-\u03b3 and allograft-rejection "
     "programmes are up in R, oxidative-phosphorylation and fatty-acid "
     "metabolism are up in NR. "
     "(F) Volcano plot of differentially expressed genes between "
     "Responders and Non-Responders (pseudo-bulk DESeq2); red, "
     "significant R-up; blue, significant NR-up (|log2FC|\u2009\u2265\u20091, "
     "padj\u2009<\u20090.05)."),

    ("Figure4_spatial_axis.png",
     "Fig. 4. Spatial axis: NMF9 and TLS hubs on Visium sections.",
     "(A) Visium (GSE238264) Responder sections with spot-level scores "
     "for key signatures (cytotoxic T, TLS, hepatocyte-differentiation); "
     "R sections show spatially clustered cytotoxic\u2013TLS hubs. "
     "(B) Same signatures on Non-Responder sections, showing diffuse "
     "cytotoxicity and weaker TLS hubs. "
     "(C) Moran\u2019s I per signature and section; Responder signatures "
     "show higher spatial autocorrelation (organised niches), Non-"
     "Responders are closer to random. "
     "(D) Violin plots of per-spot signature scores aggregated across "
     "sections, split by R vs NR, with permutation-based p-values. "
     "(E) Neighbourhood-enrichment heatmap (z-score) between cell-type "
     "deconvolved neighbourhoods; Responders show stronger CD8\u2013B-cell "
     "and CD8\u2013DC co-localisation. "
     "(F) Top TLS-score Visium tiles (spots and H&E overlays) from R and "
     "NR sections, illustrating the spatial architecture of tertiary "
     "lymphoid structure-like hubs in Responders."),

    ("Figure5_nmf_integrator.png",
     "Fig. 5. Cross-modality non-negative matrix factorisation integrator.",
     "(A) Heatmap of the top-8 gene loadings per NMF program "
     "(k\u2009=\u200910; rows, genes; columns, programs sorted by response "
     "effect); NMF9 is dominated by cytotoxic\u2013effector genes "
     "(TMSB4X, CCL5, NKG7, IL32, CCL4) while NMF3 captures a "
     "hepatocyte-differentiation program (APOA2, ALB, APOC3, APOC1). "
     "(B) Forest-style plot of standardised mean difference in per-"
     "patient program activity (Responder \u2212 Non-Responder) with "
     "annotated top genes; NMF9 is significantly elevated in Responders "
     "(*p\u2009=\u20090.028, Mann\u2013Whitney U) and NMF3 trends higher in "
     "Non-Responders (p\u2009=\u20090.099). "
     "(C) Cross-modality concordance: per-program effect size in the "
     "single-cell atlas (x-axis) versus mean activity in an independent "
     "bulk-RNA-seq cohort (GSE202069, y-axis); NMF9 is the sole program "
     "with concordant elevation in Responders across modalities. "
     "(D) Projected NMF program activity on Visium sections "
     "(rows, programs sorted by \u0394(R\u2009\u2212\u2009NR); columns, sections "
     "HCC1R\u2013HCC7NR) confirming that the sc-atlas-defined NMF9 "
     "cytotoxic program and NMF3 hepatocyte program recapitulate on "
     "spatial tissue."),

    ("Figure6_multilayer_predictor.png",
     "Fig. 6. Five-layer predictive framework and clinical translation.",
     "(A) Leave-one-out cross-validated AUC bars for each predictive "
     "layer (L1 tumour-intrinsic CNV, L2 transcription-factor regulons, "
     "L3 cell-type composition, L4 IO signatures, L5 NMF programs) and "
     "for the balanced ensemble; L3 AUC\u2009=\u20090.94, L4 AUC\u2009=\u20090.82, "
     "L1 AUC\u2009=\u20090.74, L5 AUC\u2009=\u20090.71, L2 AUC\u2009=\u20090.53, ensemble "
     "AUC\u2009=\u20090.81. "
     "(B) ROC curves per layer and for the ensemble, from the same "
     "leave-one-patient-out LASSO logistic regression. "
     "(C) Variance-partition stacked bar showing unique versus shared "
     "predictive variance across layers; tumour-intrinsic CNV "
     "contributes an independent AUC\u2009=\u20090.74 beyond cell composition. "
     "(D) Coefficient (effect) plot for the top LASSO-selected features "
     "per layer, colour-coded by direction. "
     "(E) Bulk-RNA-seq clustermap (IMbrave150 and CheckMate-459 "
     "registration-level samples) scored with the top signatures "
     "derived from the sc-atlas; rows, signatures; columns, patients. "
     "(F) Kaplan\u2013Meier survival curves stratified by composite "
     "predictor tertile on the clinical-trial bulk cohorts, showing "
     "separation in progression-free survival consistent with the "
     "sc-atlas-derived predictor."),
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
