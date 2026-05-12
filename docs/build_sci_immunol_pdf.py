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
<h1 class='title'>Genomic instability and immune-state programs define a multi-layer spatial immune signature of checkpoint-immunotherapy benefit in hepatocellular carcinoma</h1>

<p class='authors'>
Iman Sadeghi<sup>1</sup>, Bertram Bengsch<sup>1</sup>
</p>

<p class='affiliations'>
<sup>1</sup> Clinic for Internal Medicine II, Freiburg University Medical Center, Freiburg im Breisgau, Germany.
</p>

<p class='one-sentence-summary'><strong>One-sentence summary:</strong> An integrated harmonised single-cell and spatial atlas of five HCC immunotherapy cohorts and one non-tumour liver reference nominates CNV clone diversity and responder-associated immune-state programmes as candidate biomarkers, while showing that discovery-leading composition signals are not reliably portable across regimens.</p>
"""

TITLE_BLOCK_SUPP = """
<h1 class='title'>Supplementary Materials — Genomic instability and immune-state programs define a multi-layer spatial immune signature of checkpoint-immunotherapy benefit in hepatocellular carcinoma</h1>

<p class='authors'>Iman Sadeghi<sup>1</sup>, Bertram Bengsch<sup>1</sup></p>
<p class='affiliations'>
<sup>1</sup> Clinic for Internal Medicine II, Freiburg University Medical Center, Freiburg im Breisgau, Germany.
</p>
"""

SUPP_FIGURES = [
    ("SupplementaryFig_S1_atlas_QC.png",
     "Fig. S1. Atlas quality control and integration diagnostics.",
     "Per-dataset cell counts before and after QC; per-cell distributions of n_genes, n_counts, mitochondrial fraction and doublet score with thresholds; PCA variance-explained; UMAP before/after Harmony coloured by dataset; silhouette and entropy-of-batch-mixing diagnostics."),
    ("SupplementaryFig_S2_immune_UMAPs.png",
     "Fig. S2. Immune-lineage UMAPs.",
     "Major- and subtype-specific UMAPs for T-cell, B/plasma, NK, dendritic and macrophage compartments coloured by sub-cluster and by response label."),
    ("SupplementaryFig_S3_marker_annotation.png",
     "Fig. S3. Marker-based cell-type annotation.",
     "Dot plot of canonical markers across the 16 annotated populations; heatmap of per-cluster Wilcoxon markers; concordance matrix versus original GSE206325 annotations."),
    ("SupplementaryFig_S4_spatial_diagnostics.png",
     "Fig. S4. Spatial transcriptomics diagnostics.",
     "Per-section Visium QC; tissue-mask overlays; permutation null for Moran\u2019s I of TLS signature; cell-type deconvolution quality per spot."),
    ("SupplementaryFig_S5_crosscohort_visium.png",
     "Fig. S5. Cross-cohort Visium reproducibility.",
        "Side-by-side responder-associated NMF programme, cytotoxicity and TLS signature spot maps on GSE173905 and GSE238264 showing cross-cohort reproducibility of the spatial axis."),
    ("SupplementaryFig_S6_signature_stability.png",
     "Fig. S6. Signature stability under bootstrap resampling.",
     "2,000-iteration bootstrap distributions of R\u2009\u2212\u2009NR effect sizes for key signatures; NMF cophenetic-correlation curves and reconstruction-MSE elbow across k\u2009\u2208\u2009{6, 8, 10, 12, 15} supporting k\u2009=\u200910."),
    ("SupplementaryFig_S7_reproducibility.png",
     "Fig. S7. Reproducibility diagnostics.",
     "Cross-cohort Pearson correlation heatmap of mean signature scores; leave-one-cohort-out AUC per layer; Harmony vs no-integration delta in patient-level signature scores."),
    ("SupplementaryFig_S8_patient_subtyping.png",
     "Fig. S8. Patient subtyping details.",
     "Hierarchical dendrogram of the 24\u2009\u00d7\u200916 patient\u2013signature matrix; heatmap with response and cohort side-bars; patient-level pseudotime histogram by response (per-patient mean DPT; Mann\u2013Whitney U\u2009=\u200982, p\u2009=\u20090.166, Cliff\u2019s \u03b4\u2009=\u2009+0.378)."),
    ("SupplementaryFig_S9_patient_demographics.png",
     "Fig. S9. Patient and cohort demographics (primary cohort GSE206325).",
        "Patient-level metadata, sample distribution by tissue, samples per patient, cells per scRNA-seq cohort, patients per cohort across modalities, and immune-infiltration\u2009\u00d7\u2009response crosstabs for the analytical subsets."),
]


# Main-manuscript figures to embed (in order, before References).
# Legends are panel-resolved and quote effect sizes from the manuscript.
MAIN_FIGURES = [
    ("Figure1_atlas_landscape.png",
     "Fig. 1. Integrated single-cell atlas of hepatocellular carcinoma immunotherapy.",
        "(A) UMAP of 1,058,964 single cells from five harmonised HCC cohorts "
     "and one non-tumour liver reference "
        "(GSE206325, Cappuyns 2023, GSE125449, GSE140228, GSE149614 and E-MTAB-10553 / Wang 2021 + Ramachandran 2019) coloured by "
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
        "Fig. 2. Tumour-intrinsic axis of response: CNV burden, clone diversity and TME context.",
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
    "for CNV mean (n.s., p\u2009=\u20090.38), CNV standard deviation "
    "(p\u2009=\u20090.065), malignant fraction (n.s.) and CNV clone diversity "
    "(Shannon, *p\u2009=\u20090.026, Cliff\u2019s \u03b4\u2009=\u2009\u22121.03); dots, patients; box, median and IQR. "
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
        "Fig. 4. Spatial axis: responder-associated programmes and TLS hubs on Visium sections.",
     "(A) Visium (GSE238264) Responder sections with spot-level scores "
        "for responder-associated immune / TLS signatures and NMF programmes; "
        "R sections show spatially clustered immune\u2013TLS hubs. "
     "(B) Same signatures on Non-Responder sections, showing diffuse "
        "immune structure and stronger hepatocyte-differentiation dominance. "
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
        "effect); NMF10 is dominated by immunoglobulin / antigen-presentation genes "
        "(IGKC, IGHG1, SSR4, CD74, MZB1) while NMF3 captures a "
     "hepatocyte-differentiation program (APOA2, ALB, APOC3, APOC1). "
     "(B) Forest-style plot of standardised mean difference in per-"
     "patient program activity (Responder \u2212 Non-Responder) with "
        "annotated top genes; NMF10 is BH-significant in Responders, while "
        "NMF8 and NMF3 remain nominal responder- and non-responder-associated "
        "programmes. "
        "(C) Cross-modality concordance view for the sc-derived NMF programmes "
        "projected onto bulk RNA-seq and related transcriptomic datasets. "
     "(D) Projected NMF program activity on Visium sections "
     "(rows, programs sorted by \u0394(R\u2009\u2212\u2009NR); columns, sections "
        "HCC1R\u2013HCC7NR) confirming that the sc-atlas-defined NMF10 "
        "responder programme and NMF3 hepatocyte programme recapitulate on "
     "spatial tissue."),

    ("Figure6_multilayer_predictor.png",
     "Fig. 6. Five-layer predictive framework and clinical translation.",
     "(A) Leave-one-out cross-validated ROC-AUC, PR-AUC and balanced-accuracy bars for each predictive "
     "layer (L1 tumour-intrinsic CNV, L2 transcription-factor regulons, "
     "L3 cell-type composition, L4 IO signatures, L5 NMF programs) and "
        "for the balanced ensemble; L3 (no technical features) ROC-AUC\u2009=\u20090.95 (95% CI 0.83\u20131.00), "
        "L3 ROC-AUC\u2009=\u20090.96 (0.85\u20131.00), L5\u2009=\u20090.95 (0.84\u20131.00), L4\u2009=\u20090.90 (0.75\u20131.00), "
        "L1\u2009=\u20090.78 (0.48\u20131.00), L2\u2009=\u20090.55 (0.31\u20130.79), ensemble\u2009=\u20090.96 (0.87\u20131.00). "
     "Whiskers show 2,000-iteration patient-level bootstrap 95% CIs; red markers show the 1,000-iteration "
        "label-shuffle permutation 95% null; L3 (no-tech) permutation p\u2009=\u20090.022, L3 p\u2009=\u20090.014, "
        "ensemble p\u2009=\u20090.015. PR-AUC prevalence baseline 0.29. Updated panel is saved as "
     "`Fig6A_balanced_layer_metrics.png`. "
     "(B) ROC curves per layer and for the ensemble, from the same "
     "leave-one-patient-out logistic regression. "
     "(C) Variance-partition stacked bar showing unique versus shared "
        "predictive variance across layers, with the largest drop after removing L5. "
     "(D) Coefficient plot for the top-selected features "
     "per layer, colour-coded by direction. "
     "(E) IMbrave150 real per-patient Kaplan\u2013Meier stratification by intratumoural CD8 IHC "
     "density (median split), within arm: atezolizumab\u2009+\u2009bevacizumab arm PFS log-rank "
     "p\u2009=\u20090.009 (median 8.4 vs 5.6 months); sorafenib arm PFS p\u2009=\u20090.90. Cox arm\u2009\u00d7\u2009CD8 "
     "interaction HR\u2009=\u20090.57 for PFS (p\u2009=\u20090.18) and 0.36 for OS (p\u2009=\u20090.065); updated panel "
     "is `ExtendedData_Fig4_imbrave150_real.png`. "
     "(F) Bulk-RNA-seq clustermap scored with sc-atlas-derived "
     "signatures (qualitative reference)."),
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
