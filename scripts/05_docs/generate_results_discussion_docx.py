from __future__ import annotations

from pathlib import Path

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt
from PIL import Image, ImageDraw, ImageFont


REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_DIR = REPO_ROOT / "docs"
FIGURE_DIR = DOCS_DIR / "figures" / "results_discussion"
DOCX_PATH = DOCS_DIR / "manuscript_results_discussion_nature_cancer.docx"

BODY_FONT = "Arial"
HEADER_FONT = "Arial"
TITLE_FONT_SIZE = 48
PANEL_FONT_SIZE = 34
LABEL_FONT_SIZE = 30


RESULTS_SECTIONS = [
    (
        "An integrated hepatocellular carcinoma immunotherapy atlas resolves response-associated immune organization",
        [
            "We synthesized the repository’s harmonized pretreatment hepatocellular carcinoma atlas into a unified response-focused narrative spanning approximately 970,000 cells, 15,853 genes and five single-cell studies, anchored by the GSE206325 anti-PD-1 cohort with 30 responder or non-responder labels. The integrated manifold and downstream annotation outputs support the central premise that immunotherapy response is encoded at the ecosystem level rather than by a single lineage-restricted marker, because response differences emerge across T cells, myeloid compartments and tissue-context-dependent programs.",
            "Across the atlas, responders separate from non-responders by coordinated enrichment of cytotoxic, activated and late dysfunctional T-cell programs, whereas non-responders retain a comparatively progenitor-like CD8 compartment together with macrophage-polarization features that are more compatible with ineffective immune engagement. The abundance, violin, radar and tissue-stratified heatmap panels collectively indicate that the response phenotype is not simply a gain of effector cells, but a broader remodeling of the tumor microenvironment in which active immune pressure coexists with compensatory exhaustion signals.",
        ],
    ),
    (
        "Responder tumors display a coupled activation and exhaustion state rather than immune quiescence",
        [
            "A notable result emerging from the committed analyses is the paradoxical co-elevation of T-cell cytotoxicity, activation and terminal exhaustion in responders. Rather than contradicting therapeutic benefit, this pattern is consistent with a pre-existing antitumor response that has already progressed into checkpoint-constrained dysfunction before treatment and is therefore poised for reinvigoration by PD-1 or PD-L1 blockade.",
            "This interpretation is reinforced by the CD8 state analysis, which organizes the compartment into six phenotypic states from naive-like and progenitor-exhausted cells to effector, transitional, terminally exhausted and tissue-resident populations. In this framework, responders appear enriched for the downstream effector-to-terminal branches, while non-responders retain upstream progenitor-associated structure, implying that successful therapy depends on a tumor bed that has already undergone productive immune priming rather than one stalled in an early, less engaged state.",
        ],
    ),
    (
        "Multivariable modeling identifies a compact, high-performing response signature",
        [
            "The validation module nominates a compact predictive feature set at the patient level. LASSO logistic regression achieved a cross-validated area under the receiver operating characteristic curve of 0.966 and the corresponding random-forest model reached an AUC of 0.941, indicating that the response-associated programs are both biologically coherent and statistically discriminative in the available cohort.",
            "The highest-weighted predictors converge on the same biological story visible in the descriptive analyses: elevated cytotoxicity, terminal exhaustion and overall exhaustion are associated with response, whereas progenitor CD8 programs trend in the opposite direction. Bootstrap analyses further support robustness by showing non-overlapping 95% confidence intervals for the leading signatures across 50 resampling iterations, which argues that the dominant signal is stable rather than driven by a small number of influential patients.",
        ],
    ),
    (
        "Spatial and intercellular analyses localize the response program to coordinated immune neighborhoods",
        [
            "The repository’s spatial, topology and cell-cell communication outputs extend the atlas beyond cell-intrinsic states. Pseudobulk spatial scoring, high-contrast topology maps, tertiary lymphoid structure visualizations and response-prediction maps consistently suggest that beneficial immune programs are spatially organized, with exhausted-cytotoxic CD8 hubs, structured macrophage gradients and lymphoid aggregates occupying reproducible niches rather than appearing as diffuse signals.",
            "Mechanistically, ligand-receptor and pathway panels argue that these niches are coupled by conserved immunoregulatory signaling. Curated checkpoint and cytokine interaction analyses, together with hallmark pathway scoring, support a model in which responders maintain stronger interferon-licensed, antigen-engaged communication networks while also expressing inhibitory circuits expected in chronically stimulated tissues. This balance between productive immune pressure and adaptive restraint provides a plausible systems-level explanation for why checkpoint blockade succeeds in a subset of patients.",
        ],
    ),
    (
        "Cross-cohort and bulk validation support generalizability of the response architecture",
        [
            "External support comes from the independent bulk RNA-seq and cross-dataset modules. The repository includes bulk validation in GSE202069, a CIBERSORTx-compatible reference matrix spanning 15,853 genes by 60 cell types, and downstream trial-oriented deconvolution and survival plots. These outputs indicate that the response-associated immune architecture inferred from single cells is detectable at lower-resolution profiling scales relevant to translational studies.",
            "Similarly, the cross-dataset signature heatmap shows that the same broad immune axes recur across treatment-linked, treatment-naive, immune-enriched and NASH-associated cohorts. Although cohort-specific heterogeneity remains biologically meaningful, the recurrence of cytotoxic-exhausted T-cell programs, antigen-presentation-associated states and myeloid restructuring suggests that the central response biology is portable across hepatocellular carcinoma contexts.",
        ],
    ),
]

DISCUSSION_PARAGRAPHS = [
    "Taken together, the assembled results support a clinically relevant reinterpretation of pretreatment immunobiology in hepatocellular carcinoma. Response is not best understood as the presence of a uniformly inflamed microenvironment, but as the presence of a specifically organized, spatially coordinated and checkpoint-engaged immune ecosystem. Tumors that respond appear to have already crossed the threshold of immune recognition, even if that recognition is accompanied by exhaustion and compensatory inhibitory signaling.",
    "This distinction is important because superficially similar T-cell infiltration can represent very different biological states. In the current atlas, progenitor-like CD8 enrichment without corresponding cytotoxic or terminal programs aligns more closely with non-response, whereas the combination of activation, cytotoxicity and exhaustion aligns with benefit. The data therefore argue against treating exhaustion markers as simple liabilities and instead position them as context-dependent indicators of prior antitumor engagement.",
    "The spatial modules further elevate this conclusion by showing that response-associated programs cluster within structured neighborhoods that include tertiary lymphoid architectures, exhausted CD8 hubs and polarized myeloid territories. A translational implication is that clinically useful biomarkers may need to integrate composition, state and topology rather than rely on single-gene immunohistochemistry or bulk abundance estimates alone.",
    "The work also has clear limitations. The strongest quantitative claims are centered on the repository’s integrated discovery and validation outputs rather than a newly rerun analysis in this session; prospective validation in uniformly treated cohorts and orthogonal spatial assays will be required before clinical deployment. Even so, the convergence of single-cell, spatial, ligand-receptor and bulk deconvolution evidence makes a strong case that pretreatment ecosystem architecture can be developed into a practical biomarker framework for checkpoint therapy in hepatocellular carcinoma.",
]

FIGURES = [
    {
        "number": 1,
        "title": "Integrated atlas and response-associated tumor microenvironment states",
        "filename": "figure1_integrated_atlas_response_landscape.png",
        "layout": (3, 2),
        "panels": [
            ("A", "Atlas integration UMAP", "out/01_preprocessing/figures/01_integration_umap.png"),
            ("B", "Cell types and response map", "out/02a_response_tme/figures/02a_celltype_response_umap.png"),
            ("C", "Cell-type abundance by response", "out/02a_response_tme/figures/02a_abundance_barplot.png"),
            ("D", "Top immune-oncology signature violins", "out/02a_response_tme/figures/02a_signature_violins.png"),
            ("E", "Tissue-stratified signature heatmap", "out/02a_response_tme/figures/02a_tissue_signature_heatmap.png"),
        ],
        "caption": (
            "Figure 1 | Integrated atlas and response-associated tumor microenvironment states. "
            "A, Harmony-integrated single-cell manifold spanning the five hepatocellular carcinoma cohorts included in the repository, illustrating preservation of shared biological structure after batch correction. "
            "B, response-focused cell-type UMAP showing how responder and non-responder labels project onto the annotated immune and stromal landscape. "
            "C, patient-level abundance comparison highlighting the cell compartments that most clearly distinguish response classes. "
            "D, violin summaries of the highest-priority immune-oncology signatures, emphasizing the coupled increase in cytotoxic, activation and exhaustion-associated programs in responders. "
            "E, heatmap of tissue-context-dependent signature behavior, showing that response-associated programs are modulated by tumor versus normal-adjacent sampling context. "
            "Together these panels frame the central result that therapeutic benefit is encoded by ecosystem-wide immune organization rather than a single cell type."
        ),
    },
    {
        "number": 2,
        "title": "Predictive immune programs and CD8 state structure",
        "filename": "figure2_predictive_programs_cd8_states.png",
        "layout": (3, 2),
        "panels": [
            ("A", "LASSO ROC and coefficients", "out/03_validation/figures/03_roc_lasso.png"),
            ("B", "Random-forest feature importance", "out/03_validation/figures/03_rf_importance.png"),
            ("C", "Bootstrap confidence intervals", "out/03_validation/figures/03_bootstrap_ci.png"),
            ("D", "Activation versus exhaustion coupling", "out/02a_response_tme/figures/02a_activation_exhaustion_scatter.png"),
            ("E", "CD8 T-cell state distributions", "out/04_complementary/figures/04_cd8_states.png"),
            ("F", "Signature radar summary", "out/02a_response_tme/figures/02a_signature_radar.png"),
        ],
        "caption": (
            "Figure 2 | Predictive immune programs and CD8 state structure. "
            "A, patient-level LASSO performance summary, including the cross-validated receiver operating characteristic analysis with an AUC of 0.966 and the corresponding directionality of selected coefficients. "
            "B, random-forest variable importance plot showing concordant prioritization of cytotoxicity, exhaustion and related immune-state features, with model AUC reported in the repository as 0.941. "
            "C, bootstrap confidence intervals demonstrating stability of the dominant predictive signatures across 50 resampling iterations. "
            "D, scatter plot of activation versus exhaustion scores illustrating that responder tumors occupy a high-activation, high-exhaustion regime rather than an immune-quiescent state. "
            "E, six-state CD8 phenotypic framework resolving naive-like, progenitor, effector, transitional, terminally exhausted and tissue-resident programs. "
            "F, radar-chart synthesis of the principal signature shifts between response classes. "
            "These panels support the interpretation that response is associated with an antigen-experienced, checkpoint-engaged T-cell ecosystem."
        ),
    },
    {
        "number": 3,
        "title": "Spatial organization of immune response niches",
        "filename": "figure3_spatial_response_niches.png",
        "layout": (3, 2),
        "panels": [
            ("A", "Spatial topology RGB map", "out/02b_spatial_biomarkers/figures/02e_spatial_topology_rgb.png"),
            ("B", "Pseudobulk spatial signatures", "out/02b_spatial_biomarkers/figures/02b_pseudobulk_io_signatures.png"),
            ("C", "CD8 exhaustion hubs", "out/04_complementary/figures/spatial_2_cd8_exhaustion_hubs.png"),
            ("D", "Macrophage polarization map", "out/04_complementary/figures/spatial_3_macrophage_polarization.png"),
            ("E", "Tertiary lymphoid structure hubs", "out/04_complementary/figures/04g_tls_spatial_hubs.png"),
            ("F", "Spatial response prediction", "out/04_complementary/figures/spatial_4_io_response_prediction.png"),
        ],
        "caption": (
            "Figure 3 | Spatial organization of immune response niches. "
            "A, composite topology map summarizing the major spatial gradients extracted from the spatial-biomarker module. "
            "B, regional pseudobulk immune-oncology signatures, linking spatial compartments to response-relevant transcriptional programs. "
            "C, localization of exhausted CD8-rich hubs. "
            "D, mapped macrophage-polarization territories. "
            "E, tertiary lymphoid structure-focused spatial summary identifying organized immune aggregates. "
            "F, spatially resolved response-prediction overlay integrating these features into a tissue-level biomarker readout. "
            "The composite emphasizes that response-associated immune programs are organized into reproducible neighborhoods instead of being randomly distributed through the tumor bed."
        ),
    },
    {
        "number": 4,
        "title": "Mechanistic signaling, cross-cohort consistency and translational validation",
        "filename": "figure4_mechanistic_and_validation.png",
        "layout": (3, 2),
        "panels": [
            ("A", "Ligand-receptor differential scores", "out/04_complementary/figures/04_ligand_receptor.png"),
            ("B", "Ligand-receptor dot plot", "out/04_complementary/figures/04_lr_dotplot.png"),
            ("C", "Cell-cell communication network", "out/04_complementary/figures/04b_cell_communication_network.png"),
            ("D", "Hallmark pathway enrichment", "out/04_complementary/figures/04_hallmark_pathways.png"),
            ("E", "Cross-dataset signature reproducibility", "out/04_complementary/figures/04_cross_dataset_signatures.png"),
            ("F", "Bulk and trial-oriented validation", "out/04_complementary/figures/04b_bulk_trial_survival_validation.png"),
        ],
        "caption": (
            "Figure 4 | Mechanistic signaling, cross-cohort consistency and translational validation. "
            "A and B, curated ligand-receptor summaries comparing response groups across immune-checkpoint and cytokine interaction axes. "
            "C, cell-cell communication network connecting the principal cellular compartments implicated in response. "
            "D, hallmark pathway view showing how interferon, inflammatory and tissue-remodeling pathways partition across response states. "
            "E, heatmap of signature reproducibility across independent hepatocellular carcinoma cohorts, including treatment-linked and etiology-diverse datasets. "
            "F, bulk and trial-oriented validation panel illustrating how the single-cell-derived biomarker framework can be carried into lower-resolution translational cohorts. "
            "Collectively, these analyses argue that the response architecture is mechanistically coherent and generalizable beyond the discovery cohort."
        ),
    },
]

TABLES = [
    {
        "number": 1,
        "title": "Study cohorts and analytical scope represented in the manuscript-ready results",
        "caption": (
            "Table 1 | Study cohorts and analytical scope represented in the manuscript-ready results. "
            "The table summarizes the major single-cell and bulk cohorts cited throughout the Results and Discussion section and clarifies the role each dataset plays in discovery, spatial interpretation or validation."
        ),
        "columns": ["Dataset or output", "Modality", "Scale", "Primary use in this manuscript section"],
        "rows": [
            ["Integrated atlas", "scRNA-seq", "~970,000 cells × 15,853 genes", "Unified discovery framework spanning five HCC cohorts"],
            ["GSE206325", "scRNA-seq + response labels", "30 patients; ~682,000 cells", "Core responder versus non-responder discovery cohort"],
            ["GSE149614", "scRNA-seq", "~70,000 cells", "Cross-cohort consistency and ecosystem context"],
            ["GSE140228", "scRNA-seq", "~67,000 cells", "Immune-enriched validation of response-related programs"],
            ["GSE125449", "scRNA-seq", "~10,000 cells", "Baseline HCC landscape reference"],
            ["E-MTAB-10553", "scRNA-seq", "~152,000 cells", "Etiology-diverse validation, including NASH-associated disease"],
            ["GSE173905", "Spatial or pseudobulk", "Regional transcriptomic profiles", "Spatial biomarker discovery and regional program scoring"],
            ["GSE202069", "Bulk RNA-seq", "Independent validation cohort", "External confirmation of response-associated immune programs"],
        ],
    },
    {
        "number": 2,
        "title": "Principal response-associated immune features and predictive performance",
        "caption": (
            "Table 2 | Principal response-associated immune features and predictive performance. "
            "Quantitative values are taken directly from the repository documentation and linked analysis outputs. "
            "Directional terms indicate the response class in which the feature is enriched or prioritized."
        ),
        "columns": ["Analysis layer", "Key finding", "Quantitative evidence available in repository", "Interpretation"],
        "rows": [
            ["Patient-level prediction", "LASSO classifier performance", "Cross-validated AUC = 0.966", "A compact feature set strongly discriminates responders from non-responders"],
            ["Patient-level prediction", "Random-forest performance", "Cross-validated AUC = 0.941", "Independent model family reproduces the same separation"],
            ["Signature biology", "T-cell cytotoxicity", "Elevated in responders", "Responders harbor active effector programs before treatment"],
            ["Signature biology", "T-cell activation", "Elevated in responders", "Therapeutic benefit arises in already engaged immune ecosystems"],
            ["Signature biology", "Terminal exhaustion", "Elevated in responders", "Exhaustion acts as a marker of prior antigen exposure rather than simple failure"],
            ["Signature biology", "Progenitor CD8 program", "Reduced in responders", "Non-responders retain a less mature antitumor T-cell state"],
            ["Myeloid context", "Macrophage polarization signatures", "Repository summaries note macrophage-associated programs among non-responder features", "Myeloid remodeling accompanies ineffective immune control and should be interpreted with polarization-state context"],
            ["Robustness", "Bootstrap stability", "Non-overlapping 95% CIs across 50 iterations for top signatures", "Leading signals are not driven by a few outlier patients"],
        ],
    },
    {
        "number": 3,
        "title": "Mechanistic and spatial modules that extend the core response model",
        "caption": (
            "Table 3 | Mechanistic and spatial modules that extend the core response model. "
            "These analyses do not merely replicate the discovery signal; they provide biological context explaining how the response-associated programs may be organized and transmitted through tissue."
        ),
        "columns": ["Module", "Repository evidence", "Main biological message", "Translational implication"],
        "rows": [
            ["CD8 state analysis", "Six-state framework from naive-like to terminally exhausted", "Responders are enriched for downstream antigen-experienced states", "State-resolved biomarkers may outperform simple CD8 counts"],
            ["Ligand-receptor analysis", "24 curated immune-oncology interaction pairs", "Checkpoint and cytokine signaling differ systematically by response class", "Combination therapy design can target dominant inhibitory or licensing circuits"],
            ["Hallmark pathways", "20 pathway programs scored across cohorts", "Interferon, inflammatory and tissue-remodeling axes stratify response states", "Pathway-level biomarkers can complement cell-state measures"],
            ["Spatial topology", "High-contrast maps, TLS hubs and regional pseudobulk signatures", "Response-associated programs are clustered in organized immune neighborhoods", "Spatial assays may improve pretreatment patient stratification"],
            ["Bulk deconvolution", "15,853-gene by 60-cell-type reference matrix", "Single-cell signals can be propagated into translational bulk cohorts", "Biomarkers can be tested in archival and trial-scale datasets"],
            ["Cross-dataset validation", "Consistent signal across treatment-linked and etiology-diverse cohorts", "Core immune architecture is reproducible across HCC contexts", "Generalizable biomarker development is plausible"],
        ],
    },
]


def set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:fill"), fill)
    tc_pr.append(shd)


def set_default_font(document: Document) -> None:
    styles = document.styles
    for style_name in ["Normal", "Title", "Heading 1", "Heading 2"]:
        style = styles[style_name]
        style.font.name = BODY_FONT
        style._element.rPr.rFonts.set(qn("w:eastAsia"), BODY_FONT)
    styles["Normal"].font.size = Pt(10.5)
    styles["Title"].font.size = Pt(18)
    styles["Heading 1"].font.size = Pt(13)
    styles["Heading 2"].font.size = Pt(11.5)


def fit_image(image: Image.Image, box_width: int, box_height: int) -> Image.Image:
    copy = image.convert("RGB")
    copy.thumbnail((box_width, box_height))
    canvas = Image.new("RGB", (box_width, box_height), "white")
    x = (box_width - copy.width) // 2
    y = (box_height - copy.height) // 2
    canvas.paste(copy, (x, y))
    return canvas


def make_composite(spec: dict) -> Path:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    cols, rows = spec["layout"]
    cell_w, cell_h = 1700, 1180
    margin = 80
    gap = 40
    title_h = 120
    label_h = 90

    width = margin * 2 + cols * cell_w + (cols - 1) * gap
    height = margin * 2 + title_h + rows * (cell_h + label_h) + (rows - 1) * gap
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    title_font = ImageFont.load_default(size=TITLE_FONT_SIZE)
    panel_font = ImageFont.load_default(size=PANEL_FONT_SIZE)
    label_font = ImageFont.load_default(size=LABEL_FONT_SIZE)

    draw.text((margin, 20), f"Figure {spec['number']}. {spec['title']}", fill="black", font=title_font)

    for idx, panel in enumerate(spec["panels"]):
        letter, label, rel_path = panel
        row = idx // cols
        col = idx % cols
        x0 = margin + col * (cell_w + gap)
        y0 = margin + title_h + row * (cell_h + label_h + gap)
        image_path = REPO_ROOT / rel_path
        image = Image.open(image_path)
        fitted = fit_image(image, cell_w, cell_h)
        canvas.paste(fitted, (x0, y0))
        draw.rectangle((x0, y0, x0 + 86, y0 + 60), fill="black")
        draw.text((x0 + 24, y0 + 10), letter, fill="white", font=panel_font)
        draw.text((x0, y0 + cell_h + 16), label, fill="black", font=label_font)

    output_path = FIGURE_DIR / spec["filename"]
    canvas.save(output_path, dpi=(300, 300))
    return output_path


def add_paragraph(document: Document, text: str, bold_prefix: str | None = None) -> None:
    paragraph = document.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    if bold_prefix and text.startswith(bold_prefix):
        prefix_run = paragraph.add_run(bold_prefix)
        prefix_run.bold = True
        paragraph.add_run(text[len(bold_prefix):])
    else:
        paragraph.add_run(text)
    paragraph.paragraph_format.space_after = Pt(6)


def add_table(document: Document, spec: dict) -> None:
    add_paragraph(document, spec["caption"], bold_prefix=f"Table {spec['number']} | ")
    table = document.add_table(rows=1, cols=len(spec["columns"]))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.style = "Table Grid"
    hdr = table.rows[0].cells
    for idx, heading in enumerate(spec["columns"]):
        hdr[idx].text = heading
        set_cell_shading(hdr[idx], "D9EAF7")
    for row in spec["rows"]:
        cells = table.add_row().cells
        for idx, value in enumerate(row):
            cells[idx].text = value
    document.add_paragraph()


def build_docx(composites: list[tuple[dict, Path]]) -> None:
    document = Document()
    set_default_font(document)

    section = document.sections[0]
    section.top_margin = Inches(0.7)
    section.bottom_margin = Inches(0.7)
    section.left_margin = Inches(0.75)
    section.right_margin = Inches(0.75)

    title = document.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = title.add_run("Results and Discussion")
    run.bold = True
    run.font.size = Pt(18)

    subtitle = document.add_paragraph()
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER
    subtitle.add_run(
        "Nature Cancer-style manuscript section prepared from committed repository outputs"
    ).italic = True

    add_paragraph(
        document,
        "This document synthesizes the committed analysis outputs, manuscript introduction, materials and methods, and repository-level quantitative summaries into a publication-style Results and Discussion section. Composite figures were assembled from existing PNG outputs under out/ and are intended for editorial drafting rather than replacement of the original high-resolution source panels.",
    )

    document.add_heading("Results", level=1)
    for heading, paragraphs in RESULTS_SECTIONS:
        document.add_heading(heading, level=2)
        for text in paragraphs:
            add_paragraph(document, text)

    document.add_heading("Discussion", level=1)
    for text in DISCUSSION_PARAGRAPHS:
        add_paragraph(document, text)

    document.add_page_break()
    document.add_heading("Composite Figures", level=1)
    for spec, path in composites:
        document.add_picture(str(path), width=Inches(7.0))
        add_paragraph(document, spec["caption"], bold_prefix=f"Figure {spec['number']} | ")
        document.add_paragraph()

    document.add_section(WD_SECTION.NEW_PAGE)
    document.add_heading("Summary Tables", level=1)
    for spec in TABLES:
        add_table(document, spec)

    document.save(DOCX_PATH)


def main() -> None:
    composites = [(spec, make_composite(spec)) for spec in FIGURES]
    build_docx(composites)
    print(f"Wrote {DOCX_PATH}")


if __name__ == "__main__":
    main()
