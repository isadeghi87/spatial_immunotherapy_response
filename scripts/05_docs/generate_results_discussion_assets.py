#!/usr/bin/env python3
"""Generate a cautious Results/Discussion manuscript draft and composite SVG figures.

This script only uses repository-committed outputs and an audit of the analysis
scripts. It intentionally separates empirically supported figures from simulated
or placeholder analyses so the generated prose does not overstate the evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path
import json
import textwrap


PROJECT = Path(__file__).resolve().parents[2]
DOCS = PROJECT / "docs"
FIGURES = DOCS / "figures"


@dataclass(frozen=True)
class ImagePanel:
    label: str
    title: str
    rel_path: str


@dataclass(frozen=True)
class TextPanel:
    label: str
    title: str
    body: list[str]


def ensure_dirs() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)


def png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if header[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"{path} is not a PNG file")
    width = int.from_bytes(header[16:20], "big")
    height = int.from_bytes(header[20:24], "big")
    return width, height


def panel_geometry(
    index: int,
    total_cols: int,
    panel_w: int,
    panel_h: int,
    gap_x: int,
    gap_y: int,
    left: int,
    top: int,
) -> tuple[int, int]:
    row = index // total_cols
    col = index % total_cols
    x = left + col * (panel_w + gap_x)
    y = top + row * (panel_h + gap_y)
    return x, y


def add_wrapped_text(
    svg: list[str],
    text: str,
    x: int,
    y: int,
    width_chars: int,
    font_size: int = 22,
    line_height: int = 30,
    fill: str = "#203040",
    weight: str = "400",
) -> int:
    lines = textwrap.wrap(text, width=width_chars) or [""]
    svg.append(
        f'<text x="{x}" y="{y}" font-size="{font_size}" font-weight="{weight}" '
        f'fill="{fill}" font-family="Arial, Helvetica, sans-serif">'
    )
    for idx, line in enumerate(lines):
        dy = 0 if idx == 0 else line_height
        svg.append(f'<tspan x="{x}" dy="{dy}">{escape(line)}</tspan>')
    svg.append("</text>")
    return y + line_height * max(len(lines), 1)


def render_image_panel(svg: list[str], panel: ImagePanel, x: int, y: int, w: int, h: int) -> None:
    img_path = PROJECT / panel.rel_path
    img_w, img_h = png_dimensions(img_path)

    svg.append(f'<g transform="translate({x},{y})">')
    svg.append(
        f'<rect x="0" y="0" width="{w}" height="{h}" rx="18" ry="18" '
        'fill="#ffffff" stroke="#d6dee8" stroke-width="2"/>'
    )
    svg.append(
        f'<text x="24" y="42" font-size="30" font-weight="700" fill="#102a43" '
        'font-family="Arial, Helvetica, sans-serif">'
        f'{escape(panel.label)}</text>'
    )
    svg.append(
        f'<text x="78" y="42" font-size="23" font-weight="600" fill="#102a43" '
        'font-family="Arial, Helvetica, sans-serif">'
        f'{escape(panel.title)}</text>'
    )

    inner_x = 18
    inner_y = 60
    inner_w = w - 36
    inner_h = h - 78
    scale = min(inner_w / img_w, inner_h / img_h)
    draw_w = img_w * scale
    draw_h = img_h * scale
    draw_x = inner_x + (inner_w - draw_w) / 2
    draw_y = inner_y + (inner_h - draw_h) / 2

    href = escape(str(Path("../../") / Path(panel.rel_path)))
    svg.append(
        f'<image href="{href}" x="{draw_x:.1f}" y="{draw_y:.1f}" '
        f'width="{draw_w:.1f}" height="{draw_h:.1f}" preserveAspectRatio="xMidYMid meet"/>'
    )
    svg.append("</g>")


def render_text_panel(svg: list[str], panel: TextPanel, x: int, y: int, w: int, h: int) -> None:
    svg.append(f'<g transform="translate({x},{y})">')
    svg.append(
        f'<rect x="0" y="0" width="{w}" height="{h}" rx="18" ry="18" '
        'fill="#f8fbff" stroke="#d6dee8" stroke-width="2"/>'
    )
    svg.append(
        f'<text x="24" y="42" font-size="30" font-weight="700" fill="#102a43" '
        'font-family="Arial, Helvetica, sans-serif">'
        f'{escape(panel.label)}</text>'
    )
    svg.append(
        f'<text x="78" y="42" font-size="23" font-weight="600" fill="#102a43" '
        'font-family="Arial, Helvetica, sans-serif">'
        f'{escape(panel.title)}</text>'
    )
    current_y = 88
    for bullet in panel.body:
        current_y = add_wrapped_text(
            svg,
            f"• {bullet}",
            x=28,
            y=current_y,
            width_chars=44,
            font_size=20,
            line_height=28,
            fill="#334e68",
        )
        current_y += 8
    svg.append("</g>")


def write_svg(
    output_path: Path,
    title: str,
    panels: list[ImagePanel | TextPanel],
    footer: str,
    ncols: int = 3,
    panel_w: int = 500,
    panel_h: int = 410,
    gap_x: int = 24,
    gap_y: int = 24,
) -> None:
    rows = (len(panels) + ncols - 1) // ncols
    width = 80 + ncols * panel_w + (ncols - 1) * gap_x
    height = 150 + rows * panel_h + (rows - 1) * gap_y + 80

    svg: list[str] = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#eef4fb"/>',
        f'<text x="40" y="54" font-size="34" font-weight="700" fill="#102a43" '
        'font-family="Arial, Helvetica, sans-serif">'
        f'{escape(title)}</text>',
        '<text x="40" y="88" font-size="18" fill="#486581" '
        'font-family="Arial, Helvetica, sans-serif">'
        'Composite figure assembled from repository-committed outputs.</text>',
    ]

    for idx, panel in enumerate(panels):
        x, y = panel_geometry(
            index=idx,
            total_cols=ncols,
            panel_w=panel_w,
            panel_h=panel_h,
            gap_x=gap_x,
            gap_y=gap_y,
            left=40,
            top=110,
        )
        if isinstance(panel, ImagePanel):
            render_image_panel(svg, panel, x, y, panel_w, panel_h)
        else:
            render_text_panel(svg, panel, x, y, panel_w, panel_h)

    footer_y = height - 28
    svg.append(
        f'<text x="40" y="{footer_y}" font-size="17" fill="#486581" '
        'font-family="Arial, Helvetica, sans-serif">'
        f'{escape(footer)}</text>'
    )
    svg.append("</svg>")
    output_path.write_text("\n".join(svg), encoding="utf-8")


def write_provenance_json(output_path: Path) -> None:
    payload = {
        "empirical_modules": [
            {
                "module": "02a_response_tme_characterization.py",
                "status": "empirical_when_source_data_are_present",
                "notes": "Uses integrated_io_atlas.h5ad and GSE206325 metadata to derive abundance, DEG, and signature comparisons.",
            },
            {
                "module": "02b_spatial_biomarker_discovery.py",
                "status": "mixed_empirical",
                "notes": "Pseudo-bulk GSE173905 logic is real-data driven; optional HCC_omics Visium branch depends on external files not committed here.",
            },
            {
                "module": "02c_bulk_deconvolution.py",
                "status": "partially_empirical",
                "notes": "Reference matrix generation is empirical; actual trial deconvolution requires missing bulk matrices.",
            },
            {
                "module": "03_validation_downstream.py",
                "status": "empirical_when_source_data_are_present",
                "notes": "Internal validation, bootstrap, and pathway analyses are built from the integrated atlas and derived tables.",
            },
        ],
        "illustrative_modules": [
            "02e_spatial_topology.py",
            "02f_regulon_tf_networks.py",
            "02g_infercnv.py",
            "02h_cnv_tme_integration.py",
            "03b_clinical_survival.py",
            "04_complementary/04b_comprehensive_discoveries.py",
            "04_complementary/04c_tcell_trajectories.py",
            "04_complementary/04d_cell_communication.py",
            "04_complementary/04e_tcr_dynamics.py",
            "04_complementary/04f_metabolic_reprogramming.py",
            "04_complementary/04g_tls_spatial_hubs.py",
            "04_complementary/04h_spatial_high_contrast_maps.py",
        ],
    }
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_markdown(output_path: Path) -> None:
    text = """# Results and Discussion (Nature Immunology-style draft)

## Scope and provenance note

This draft is based on the figures committed in `out/` and on a script-level audit of the repository. It is intentionally conservative: only result streams that are empirically grounded in the single-cell atlas, pseudo-bulk spatial cohort, internal validation workflow, or reference-matrix construction are interpreted as evidence. Scripts that generate simulated, mock, or placeholder outputs are treated as hypothesis-generating illustrations rather than manuscript-grade results.

## Results

### An integrated hepatocellular carcinoma immunotherapy atlas resolves a responder-associated inflamed but regulated microenvironment

The preprocessing and response-characterization outputs support a coherent view of the primary atlas as a biologically structured reference rather than a response-labelled artefact. The integration, cluster, and annotation panels indicate that the atlas preserves major immune and stromal lineages, while the response projection suggests that responder and non-responder cells are distributed across shared cellular continua rather than occupying wholly distinct manifolds. This is important, because it frames immunotherapy outcome as a shift in state composition and activation depth within a common tumour ecosystem rather than as the appearance of a responder-exclusive lineage.

Within that shared atlas, the committed response panels converge on a specific biological axis. Responders appear enriched for activation, cytotoxicity and antigen-presentation programmes, yet these same samples also retain expression features associated with exhaustion and checkpoint engagement. Rather than undermining the response signal, this co-occurrence is biologically plausible: in successful anti-PD-1 settings, exhausted-like transcriptional features often mark antigen-experienced T cells that remain functionally engaged and pharmacologically targetable. The activation-versus-exhaustion scatter therefore reads less as a contradiction than as evidence that productive immunity in HCC is accompanied by sustained antigen pressure and compensatory inhibitory signalling. By contrast, the abundance and composition panels suggest that non-response is not simply an absence of T cells, but a redistribution toward less productive immune states and myeloid or stromal contexts that are less compatible with cytotoxic execution.

The tissue-stratified heatmap sharpens this interpretation by implying that the tumour compartment, rather than blood or uninvolved tissue, is where the response-associated transcriptional differences are most legible. That observation is consistent with an HCC model in which therapeutic efficacy depends on local immune organisation, tumour-proximal antigen presentation and the balance between intratumoural activation and exclusion.

![Composite Figure 1](figures/results_figure_1_atlas_response.svg)

**Figure 1 legend.** Integrated atlas and response-state synthesis. Panels assemble the committed preprocessing and `02a_response_tme` outputs to show atlas integration, response projection, signature shifts, activation–exhaustion coupling, abundance differences and tissue-specific structure.

### Spatial and pseudo-bulk analyses extend the same response axis into recurrence-associated tissue states

The pseudo-bulk spatial outputs from the GSE173905 workflow are directionally concordant with the single-cell atlas. Samples grouped as recurrence-sensitive retain stronger immune-signature structure, whereas early-recurrence samples trend toward a less favourable balance of effector and suppressive programmes. Because these data are analysed as pseudo-bulk rather than true spot-resolved histology, the strongest defensible conclusion is that recurrence state tracks with the same broad immune axis seen in the single-cell response atlas: immune-active, interferon-competent and TLS-like programmes align with the more favourable state, whereas TGF-β-rich or immune-excluded programmes align with poorer disease control.

The deconvolution module contributes a second translational layer, but only in a restricted sense. The committed outputs support the construction of a diverse scRNA-seq-derived reference matrix and show that the reference retains structured correlation among annotated cell states, which is necessary for downstream CIBERSORTx-style analyses. However, the repository does not currently version trial-level bulk count matrices or deconvolution tables, so any claim of external clinical validation should stop at methodological readiness rather than completed proof. The responsible interpretation is therefore that the project has built a plausible translational bridge from the single-cell atlas to bulk trial material, but that the bridge is only partially traversed in the committed state of the repository.

![Composite Figure 2](figures/results_figure_2_spatial_validation.svg)

**Figure 2 legend.** Spatial, pseudo-bulk and translational bridge. Panels summarize pseudo-bulk spatial signatures and embeddings from `02b_spatial_biomarkers`, together with the scRNA-seq-derived deconvolution reference and the internal prediction framework from `03_validation`.

### Internal validation supports directional stability, but not yet full external generalization

The most convincing robustness evidence in the repository comes from internal triangulation rather than from independent external cohorts. The LASSO and random-forest panels suggest that the same families of features drive classification irrespective of modelling framework, which argues against a wholly model-specific result. The signature-correlation map further indicates that responder-associated programmes cohere into a network rather than behaving as isolated markers. Finally, the bootstrap confidence interval panel suggests that the sign of the dominant effects is stable to resampling, which is especially important in a setting where class sizes and patient numbers are likely modest.

The pathway enrichment panel provides mechanistic reinforcement for this interpretation. Responder-upregulated genes map onto antigen processing, cytokine signalling, co-stimulation and related immune programmes rather than onto a purely proliferative or stromal signature. Taken together, these outputs support a consistent mechanistic model in which response emerges from an inflamed, antigen-aware microenvironment that is simultaneously under inhibitory pressure, while non-response reflects either failure to achieve that inflamed state or diversion into suppressive and excluding niches.

At the same time, the current validation should still be considered internal. The repository state does not provide versioned numerical tables for the committed ROC, feature-importance, or pathway-enrichment figures, and the clinical survival workflow present in `03b_clinical_survival.py` is simulated rather than data-derived. The direction of the findings is therefore more robust than any exact effect size, threshold or survival estimate.

## Consistency and robustness assessment

Across the empirically supported modules, the central message is consistent: favourable outcome aligns with a tumour-localized immune-active state characterized by cytotoxic and interferon-linked programmes, whereas poor outcome aligns with exclusionary, suppressive or developmentally incomplete immune states. This conclusion is reinforced independently by cell-state projection, abundance analysis, signature-level contrasts, pseudo-bulk recurrence grouping, pathway enrichment, model-based prediction and bootstrap resampling.

The main robustness limitation is provenance rather than biological contradiction. Several visually polished figures in the repository are generated from random, simulated or placeholder data, including the CNV dashboards, clinical survival curves, trajectory plots, cell–cell communication panels, metabolic reprogramming panels, TLS spatial hubs and high-contrast spatial maps. These outputs may still be useful for figure design or future hypothesis framing, but they should not be cited as evidence in a manuscript until they are replaced by analyses that consume committed source data.

![Composite Figure 3](figures/results_figure_3_consistency_robustness.svg)

**Figure 3 legend.** Robustness and provenance audit. The upper row summarizes the internally supported stability analyses; the lower row records which modules are empirically grounded versus illustrative so that manuscript claims remain matched to evidence.

## Discussion

Viewed together, the supported outputs argue for a model of HCC immunotherapy response in which successful treatment is associated with an already primed but incompletely restrained immune microenvironment. The key biological point is that responder tumours are not simply “less exhausted”; instead, they appear to host immune programmes in which effector activation, interferon exposure, antigen presentation and checkpoint engagement rise together. In this framework, checkpoint expression is best interpreted as evidence of immune engagement under chronic antigenic stress, and therapeutic benefit arises when that stressed immune compartment remains recoverable.

This interpretation helps reconcile an apparently paradoxical observation that often complicates IO biomarker studies: responder samples can show both stronger cytotoxicity and stronger exhaustion-like signatures. In liver cancer, where persistent antigen exposure and tolerogenic myeloid circuits are common, an activated-plus-exhausted state may represent the precise population most likely to benefit from checkpoint blockade. Conversely, non-responder tumours may fail for at least two non-exclusive reasons: they may never establish a sufficiently inflamed, antigen-presenting niche, or they may channel infiltrating lymphocytes into progenitor-like, excluded or suppressive states that prevent terminal effector function.

The translational implication is that abundance-only biomarkers are unlikely to be sufficient. The repository’s strongest signal is not merely that responders have more immune cells, but that they occupy a different immune state-space, with coordinated changes across T-cell, myeloid and antigen-presentation programmes. This supports a composite biomarker strategy in which cell-state composition, effector-to-suppressive balance and tissue context are integrated rather than considered separately.

The discussion must also remain candid about current limits. First, the raw and intermediate data required to recompute all numerical results are not committed in the repository snapshot used here, so the present document is a figure-driven interpretation rather than a full reanalysis. Second, external clinical validation is not yet complete in version-controlled form. Third, multiple attractive downstream figures are explicitly simulated and should be excluded from any final manuscript evidence chain until they are replaced by data-backed analyses. The next high-value step is therefore not to broaden the claim set, but to harden the existing one by versioning the missing trial-level matrices, committing derived tables, and re-running the complementary modules on real data.

## Recommended manuscript boundaries

| Safe to use as evidence now | Use only after re-analysis on real data |
| --- | --- |
| Atlas integration and annotation (`out/01_preprocessing/figures/*`) | CNV dashboards (`02g`, `02h`) |
| Response-state and signature figures in `out/02a_response_tme/figures/*` | Simulated clinical survival (`03b_clinical_survival.py`) |
| Pseudo-bulk GSE173905 figures in `out/02b_spatial_biomarkers/figures/02b_*` | Complementary trajectory, communication, TCR, metabolic and TLS figures in `out/04_complementary/figures/*` |
| Reference-matrix figures in `out/02c_bulk_deconvolution/figures/02c_reference_*` | High-contrast spatial maps in `04h_spatial_high_contrast_maps.py` |
| Internal validation figures in `out/03_validation/figures/03_*` | Any unversioned external-trial claim |

## Deliverables generated by script

- `docs/results_discussion_nature_immunology.md`
- `docs/figures/results_figure_1_atlas_response.svg`
- `docs/figures/results_figure_2_spatial_validation.svg`
- `docs/figures/results_figure_3_consistency_robustness.svg`
- `docs/results_provenance_audit.json`
"""
    output_path.write_text(text, encoding="utf-8")


def main() -> None:
    ensure_dirs()

    write_svg(
        output_path=FIGURES / "results_figure_1_atlas_response.svg",
        title="Figure 1 | Integrated atlas and response-state organization",
        panels=[
            ImagePanel("A", "Atlas integration", "out/01_preprocessing/figures/01_integration_umap.png"),
            ImagePanel("B", "Response projection on the atlas", "out/02a_response_tme/figures/02a_celltype_response_umap.png"),
            ImagePanel("C", "Responder vs non-responder signature shifts", "out/02a_response_tme/figures/02a_signature_violins.png"),
            ImagePanel("D", "Activation–exhaustion coupling", "out/02a_response_tme/figures/02a_activation_exhaustion_scatter.png"),
            ImagePanel("E", "Cell-type abundance differences", "out/02a_response_tme/figures/02a_abundance_barplot.png"),
            ImagePanel("F", "Tumour-localized signature structure", "out/02a_response_tme/figures/02a_tissue_signature_heatmap.png"),
        ],
        footer="Supported primarily by 01_preprocessing and 02a_response_tme outputs.",
    )

    write_svg(
        output_path=FIGURES / "results_figure_2_spatial_validation.svg",
        title="Figure 2 | Spatial, pseudo-bulk and translational bridge",
        panels=[
            ImagePanel("A", "Pseudo-bulk IO signature landscape", "out/02b_spatial_biomarkers/figures/02b_pseudobulk_io_signatures.png"),
            ImagePanel("B", "Pseudo-bulk sample embeddings", "out/02b_spatial_biomarkers/figures/02b_pseudobulk_embeddings.png"),
            ImagePanel("C", "scRNA-seq reference heatmap", "out/02c_bulk_deconvolution/figures/02c_reference_heatmap.png"),
            ImagePanel("D", "Reference-state correlation structure", "out/02c_bulk_deconvolution/figures/02c_reference_correlation.png"),
            ImagePanel("E", "Internal LASSO response prediction", "out/03_validation/figures/03_roc_lasso.png"),
            ImagePanel("F", "Random-forest feature ranking", "out/03_validation/figures/03_rf_importance.png"),
        ],
        footer="Combines empirical pseudo-bulk figures, reference construction, and internal validation outputs.",
    )

    write_svg(
        output_path=FIGURES / "results_figure_3_consistency_robustness.svg",
        title="Figure 3 | Consistency, robustness and provenance boundaries",
        panels=[
            ImagePanel("A", "Bootstrap stability of key effects", "out/03_validation/figures/03_bootstrap_ci.png"),
            ImagePanel("B", "Signature correlation network", "out/03_validation/figures/03_signature_correlation.png"),
            ImagePanel("C", "Responder-upregulated pathway structure", "out/03_validation/figures/03_pathway_enrichment.png"),
            TextPanel(
                "D",
                "Cross-module consistency",
                [
                    "Atlas, pseudo-bulk and internal modelling all support an immune-active responder state.",
                    "Activation, IFN-related signalling and antigen presentation recur across multiple panels.",
                    "Bootstrap resampling preserves the direction of the dominant signature effects.",
                ],
            ),
            TextPanel(
                "E",
                "Evidence that remains incomplete",
                [
                    "Trial-level bulk count matrices and deconvolution outputs are not versioned in the repository snapshot.",
                    "The committed figures do not include the numerical tables needed to report exact effect sizes in this draft.",
                    "Clinical survival validation is present only as a simulated workflow and should not be cited as evidence.",
                ],
            ),
            TextPanel(
                "F",
                "Illustrative modules excluded from claims",
                [
                    "Synthetic or placeholder analyses include 02e-02h, 03b and the 04_complementary suite.",
                    "These outputs are useful for hypothesis framing and figure design, but not for manuscript-grade inference.",
                    "Final publication claims should be limited to modules with committed, data-backed provenance.",
                ],
            ),
        ],
        footer="Lower row is text-only by design so the evidence boundary is explicit in the figure itself.",
    )

    write_markdown(DOCS / "results_discussion_nature_immunology.md")
    write_provenance_json(DOCS / "results_provenance_audit.json")

    print("Generated:")
    print(f"  - {DOCS / 'results_discussion_nature_immunology.md'}")
    print(f"  - {FIGURES / 'results_figure_1_atlas_response.svg'}")
    print(f"  - {FIGURES / 'results_figure_2_spatial_validation.svg'}")
    print(f"  - {FIGURES / 'results_figure_3_consistency_robustness.svg'}")
    print(f"  - {DOCS / 'results_provenance_audit.json'}")


if __name__ == "__main__":
    main()
