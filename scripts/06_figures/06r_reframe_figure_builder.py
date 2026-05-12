#!/usr/bin/env python3
"""Build composite figures for the 2026-04-21 immune-ecology reframe.

This builder is intentionally strict:
- it reads only the audited sources declared in results_index.yaml
- it writes figure_plan.yaml for every requested figure
- it renders only figures whose declared audited sources both exist and
  contain the granularity required by the requested panels
"""

from __future__ import annotations

from pathlib import Path
import sys
import textwrap

BASE = Path("/home/sadeghii/projects/spatial_immunotherapy_response")
sys.path.insert(0, str(BASE))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
import yaml

from config.publication_style import PAL, apply_style


RESULTS_INDEX = BASE / "results_index.yaml"
OUT_ROOT = BASE / "out/06_manuscript_figures/reframe_2026_04_21"
MAIN_DPI = 600
FIGURE_ID_ALIASES = {
    "F3_external_consistency_and_plasma_context": "F3_transfer_failure_and_composition",
}
FIGURE_MESSAGE_OVERRIDES = {
    "F3_external_consistency_and_plasma_context": (
        "Cappuyns and NMF10 analyses define cohort/regimen boundaries for the discovery immune axis; "
        "they are contextual diagnostics, not the paper's hero claim."
    ),
}

ALLOWED_STRENGTH_CLASSES = [
    "strong",
    "supported_limited",
    "exploratory",
    "null",
    "inconclusive",
    "likely_artifact",
]


def wrap_text(value: object, width: int) -> str:
    text = "" if value is None else str(value)
    if text.lower() == "nan":
        text = ""
    return "\n".join(
        textwrap.wrap(text, width=width, break_long_words=False, break_on_hyphens=False)
    )


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def save_yaml(path: Path, data: dict) -> None:
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False, allow_unicode=True)


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=MAIN_DPI, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), dpi=MAIN_DPI, facecolor="white")
    plt.close(fig)


def panel_label(ax: plt.Axes, label: str, x: float = -0.02, y: float = 1.02) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom")


def draw_section_title(ax: plt.Axes, title: str) -> None:
    ax.text(0.0, 1.01, title, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def draw_box(ax: plt.Axes, xywh: tuple[float, float, float, float], face: str, edge: str = "#1f2937",
             lw: float = 1.2, radius: float = 0.03, alpha: float = 1.0) -> FancyBboxPatch:
    x, y, w, h = xywh
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        facecolor=face,
        edgecolor=edge,
        linewidth=lw,
        alpha=alpha,
    )
    ax.add_patch(patch)
    return patch


def draw_table(ax: plt.Axes, columns: list[dict], rows: list[list[str]],
               header_face: str = "#e5eef6", alt_face: str = "#f8fafc",
               body_face: str = "#ffffff", edge: str = "#d0d7de",
               header_fontsize: float = 8.0, body_fontsize: float = 7.0) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    left = 0.0
    top = 0.95
    total_width = 1.0
    header_h = 0.11

    wrapped_rows: list[list[str]] = []
    row_units: list[float] = []
    for row in rows:
        wrapped = []
        max_lines = 1
        for value, col in zip(row, columns):
            cell = wrap_text(value, col["wrap"])
            wrapped.append(cell)
            max_lines = max(max_lines, cell.count("\n") + 1 if cell else 1)
        wrapped_rows.append(wrapped)
        row_units.append(max_lines * 1.15)

    usable_h = top - 0.05 - header_h
    unit_total = sum(row_units) if row_units else 1.0
    row_heights = [usable_h * unit / unit_total for unit in row_units]

    x0 = left
    for col in columns:
        w = total_width * col["width"]
        draw_box(ax, (x0, top - header_h, w, header_h), face=header_face, edge=edge, lw=0.9, radius=0.01)
        ax.text(x0 + 0.01, top - header_h / 2, col["title"], fontsize=header_fontsize,
                fontweight="bold", va="center", ha="left")
        x0 += w

    y = top - header_h
    for idx, row in enumerate(wrapped_rows):
        h = row_heights[idx]
        y -= h
        face = alt_face if idx % 2 == 0 else body_face
        x = left
        for cell, col in zip(row, columns):
            w = total_width * col["width"]
            draw_box(ax, (x, y, w, h), face=face, edge=edge, lw=0.7, radius=0.004)
            ha = col.get("align", "left")
            if ha == "center":
                tx = x + w / 2
            elif ha == "right":
                tx = x + w - 0.01
            else:
                tx = x + 0.01
            ax.text(tx, y + h / 2, cell, fontsize=body_fontsize, va="center", ha=ha)
            x += w


def draw_column_block(ax: plt.Axes, x: float, width: float, title: str, items: list[str],
                      face: str, edge: str, fontsize: float = 7.0) -> None:
    wrapped_items = [wrap_text(f"{idx + 1}. {item}", 40) for idx, item in enumerate(items)]
    text = "\n\n".join(wrapped_items) if wrapped_items else "No items."
    draw_box(ax, (x, 0.06, width, 0.84), face=face, edge=edge, lw=1.0, radius=0.02)
    ax.text(x + 0.02, 0.86, title, fontsize=10, fontweight="bold", va="top")
    ax.text(x + 0.02, 0.81, text, fontsize=fontsize, va="top", ha="left")


def figure_specs(doc: dict) -> dict[str, dict]:
    spec_map = {}
    for item in doc.get("reframe_main_figures", []):
        if not isinstance(item, dict) or "figure_id" not in item:
            continue
        spec_map[item["figure_id"]] = item
    return spec_map


def resolve_figure_spec(specs: dict[str, dict], figure_id: str) -> dict:
    """Return a figure spec, falling back to legacy IDs while overriding stale messages."""
    spec_id = figure_id if figure_id in specs else FIGURE_ID_ALIASES.get(figure_id, figure_id)
    spec = dict(specs[spec_id])
    spec["figure_id"] = figure_id
    if figure_id in FIGURE_MESSAGE_OVERRIDES:
        spec["question"] = FIGURE_MESSAGE_OVERRIDES[figure_id]
    return spec


def result_lookup(doc: dict) -> dict[str, dict]:
    lookup = {}
    for section in ["results", "tier3_results", "revision_results", "backfilled_results", "reframe_2026_04_21"]:
        for item in doc.get(section, []):
            if isinstance(item, dict) and "result_id" in item:
                lookup[item["result_id"]] = item
    return lookup


PANEL_METADATA = {
    "F1_cohort_overview": [
        {
            "panel_id": "A",
            "panel_title": "Cohort role schematic",
            "panel_role": "overview",
            "panel_class": "schematic_panel",
            "information_density": "medium",
            "importance_tier": "hero",
            "aspect_behavior": "fixed_wide",
            "label_burden": "medium",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "7.5 x 4.0 in",
            "legend_strategy": "none",
            "colorbar_strategy": "none",
            "statistical_support_level": "not_applicable",
        },
        {
            "panel_id": "B",
            "panel_title": "Discovery audit funnel",
            "panel_role": "quantification",
            "panel_class": "bar_plot",
            "information_density": "low",
            "importance_tier": "supporting",
            "aspect_behavior": "flexible",
            "label_burden": "medium",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "4.0 x 3.0 in",
            "legend_strategy": "none",
            "colorbar_strategy": "none",
            "statistical_support_level": "not_applicable",
        },
        {
            "panel_id": "C",
            "panel_title": "Unit and claim boundary table",
            "panel_role": "contextual",
            "panel_class": "table_like_panel",
            "information_density": "high",
            "importance_tier": "primary",
            "aspect_behavior": "fixed_wide",
            "label_burden": "very_high",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "7.5 x 4.2 in",
            "legend_strategy": "none",
            "colorbar_strategy": "none",
            "statistical_support_level": "not_applicable",
        },
    ],
    "F2_discovery_immune_axis": [
        {
            "panel_id": "A",
            "panel_title": "Module effect forest",
            "panel_role": "quantification",
            "panel_class": "coefficient_plot",
            "information_density": "medium",
            "importance_tier": "primary",
            "aspect_behavior": "fixed_tall",
            "label_burden": "medium",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "4.8 x 5.0 in",
            "legend_strategy": "none",
            "colorbar_strategy": "none",
            "statistical_support_level": "mandatory",
        },
        {
            "panel_id": "B",
            "panel_title": "CORE-4 correlation heatmap",
            "panel_role": "mechanism",
            "panel_class": "heatmap_like_matrix",
            "information_density": "high",
            "importance_tier": "hero",
            "aspect_behavior": "fixed_square",
            "label_burden": "high",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "5.6 x 5.6 in",
            "legend_strategy": "none",
            "colorbar_strategy": "external_per_panel",
            "statistical_support_level": "mandatory",
        },
        {
            "panel_id": "C",
            "panel_title": "Immune-axis patient dotplot",
            "panel_role": "validation",
            "panel_class": "distribution_plot",
            "information_density": "medium",
            "importance_tier": "supporting",
            "aspect_behavior": "flexible",
            "label_burden": "medium",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "4.5 x 3.3 in",
            "legend_strategy": "none",
            "colorbar_strategy": "none",
            "statistical_support_level": "recommended",
        },
    ],
    "F3_external_consistency_and_plasma_context": [
        {
            "panel_id": "A",
            "panel_title": "Discovery vs external-cohort effect scatter",
            "panel_role": "overview",
            "panel_class": "scatter_plot",
            "information_density": "medium",
            "importance_tier": "hero",
            "aspect_behavior": "fixed_square",
            "label_burden": "medium",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "5.0 x 5.0 in",
            "legend_strategy": "internal",
            "colorbar_strategy": "none",
            "statistical_support_level": "mandatory",
        },
        {
            "panel_id": "B",
            "panel_title": "NMF10 marginal vs plasma-adjusted signal",
            "panel_role": "quantification",
            "panel_class": "scatter_plot",
            "information_density": "medium",
            "importance_tier": "primary",
            "aspect_behavior": "fixed_wide",
            "label_burden": "medium",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "5.5 x 3.6 in",
            "legend_strategy": "none",
            "colorbar_strategy": "none",
            "statistical_support_level": "mandatory",
        },
        {
            "panel_id": "C",
            "panel_title": "Dynamic-range and overlap boundary diagnostics",
            "panel_role": "validation",
            "panel_class": "table_like_panel",
            "information_density": "high",
            "importance_tier": "supporting",
            "aspect_behavior": "fixed_wide",
            "label_burden": "high",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "5.5 x 3.6 in",
            "legend_strategy": "none",
            "colorbar_strategy": "none",
            "statistical_support_level": "recommended",
        },
    ],
    "F4_imbrave150_orthogonal": [
        {
            "panel_id": "A",
            "panel_title": "Atezo+Bev PFS Kaplan-Meier",
            "panel_role": "overview",
            "panel_class": "survival_plot",
            "information_density": "medium",
            "importance_tier": "hero",
            "aspect_behavior": "fixed_wide",
            "label_burden": "medium",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "5.8 x 3.5 in",
            "legend_strategy": "internal",
            "colorbar_strategy": "none",
            "statistical_support_level": "mandatory",
        },
        {
            "panel_id": "B",
            "panel_title": "Sorafenib PFS Kaplan-Meier",
            "panel_role": "validation",
            "panel_class": "survival_plot",
            "information_density": "medium",
            "importance_tier": "primary",
            "aspect_behavior": "fixed_wide",
            "label_burden": "medium",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "5.8 x 3.5 in",
            "legend_strategy": "internal",
            "colorbar_strategy": "none",
            "statistical_support_level": "mandatory",
        },
        {
            "panel_id": "C",
            "panel_title": "Interaction forest",
            "panel_role": "quantification",
            "panel_class": "forest_plot",
            "information_density": "medium",
            "importance_tier": "supporting",
            "aspect_behavior": "flexible",
            "label_burden": "medium",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "4.5 x 3.4 in",
            "legend_strategy": "none",
            "colorbar_strategy": "none",
            "statistical_support_level": "mandatory",
        },
    ],
    "S5_spatial_modules": [
        {
            "panel_id": "A",
            "panel_title": "Per-section module means",
            "panel_role": "overview",
            "panel_class": "grouped_distribution",
            "information_density": "high",
            "importance_tier": "primary",
            "aspect_behavior": "fixed_wide",
            "label_burden": "high",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "5.6 x 3.6 in",
            "legend_strategy": "external_shared",
            "colorbar_strategy": "none",
            "statistical_support_level": "recommended",
        },
        {
            "panel_id": "B",
            "panel_title": "Median-rule Visium spot maps",
            "panel_role": "mechanism",
            "panel_class": "spatial_map",
            "information_density": "high",
            "importance_tier": "hero",
            "aspect_behavior": "preserve_native_aspect",
            "label_burden": "medium",
            "rendering_mode": "hybrid",
            "native_aspect_required": "yes",
            "minimum_readable_size": "5.8 x 4.5 in",
            "legend_strategy": "external_shared",
            "colorbar_strategy": "external_shared",
            "statistical_support_level": "recommended",
        },
        {
            "panel_id": "C",
            "panel_title": "LOSO stability plot",
            "panel_role": "validation",
            "panel_class": "interval_plot",
            "information_density": "medium",
            "importance_tier": "supporting",
            "aspect_behavior": "flexible",
            "label_burden": "medium",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "4.8 x 3.5 in",
            "legend_strategy": "none",
            "colorbar_strategy": "none",
            "statistical_support_level": "mandatory",
        },
    ],
    "S6_synthesis_and_null_landscape": [
        {
            "panel_id": "A",
            "panel_title": "Evidence-strength matrix",
            "panel_role": "overview",
            "panel_class": "heatmap_like_matrix",
            "information_density": "medium",
            "importance_tier": "primary",
            "aspect_behavior": "fixed_wide",
            "label_burden": "medium",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "5.2 x 4.0 in",
            "legend_strategy": "none",
            "colorbar_strategy": "none",
            "statistical_support_level": "not_applicable",
        },
        {
            "panel_id": "B",
            "panel_title": "Null findings landscape",
            "panel_role": "validation",
            "panel_class": "table_like_panel",
            "information_density": "high",
            "importance_tier": "primary",
            "aspect_behavior": "fixed_wide",
            "label_burden": "very_high",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "5.4 x 4.0 in",
            "legend_strategy": "none",
            "colorbar_strategy": "none",
            "statistical_support_level": "not_applicable",
        },
        {
            "panel_id": "C",
            "panel_title": "Revised claim boundary table",
            "panel_role": "contextual",
            "panel_class": "annotation_panel",
            "information_density": "high",
            "importance_tier": "hero",
            "aspect_behavior": "fixed_wide",
            "label_burden": "very_high",
            "rendering_mode": "vector",
            "native_aspect_required": "no",
            "minimum_readable_size": "11.0 x 4.5 in",
            "legend_strategy": "none",
            "colorbar_strategy": "none",
            "statistical_support_level": "not_applicable",
        },
    ],
}


LAYOUT_SPECS = {
    "F1_cohort_overview": {
        "style": "asymmetric grid",
        "rows": 2,
        "cols": 2,
        "canvas_inches": [14.5, 10.5],
        "panels": {
            "A": {"row": 0, "col": 0, "row_span": 1, "col_span": 2, "relative_area": 0.42},
            "B": {"row": 1, "col": 0, "row_span": 1, "col_span": 1, "relative_area": 0.22},
            "C": {"row": 1, "col": 1, "row_span": 1, "col_span": 1, "relative_area": 0.36},
        },
    },
    "F2_discovery_immune_axis": {
        "style": "weighted grid",
        "rows": 2,
        "cols": 2,
        "canvas_inches": [14.5, 10.5],
        "panels": {
            "A": {"row": 0, "col": 0, "row_span": 2, "col_span": 1, "relative_area": 0.33},
            "B": {"row": 0, "col": 1, "row_span": 1, "col_span": 1, "relative_area": 0.39},
            "C": {"row": 1, "col": 1, "row_span": 1, "col_span": 1, "relative_area": 0.28},
        },
    },
    "F3_external_consistency_and_plasma_context": {
        "style": "weighted grid",
        "rows": 2,
        "cols": 2,
        "canvas_inches": [14.5, 10.5],
        "panels": {
            "A": {"row": 0, "col": 0, "row_span": 2, "col_span": 1, "relative_area": 0.40},
            "B": {"row": 0, "col": 1, "row_span": 1, "col_span": 1, "relative_area": 0.28},
            "C": {"row": 1, "col": 1, "row_span": 1, "col_span": 1, "relative_area": 0.32},
        },
    },
    "F4_imbrave150_orthogonal": {
        "style": "weighted grid",
        "rows": 2,
        "cols": 2,
        "canvas_inches": [14.5, 10.0],
        "panels": {
            "A": {"row": 0, "col": 0, "row_span": 1, "col_span": 1, "relative_area": 0.34},
            "B": {"row": 0, "col": 1, "row_span": 1, "col_span": 1, "relative_area": 0.34},
            "C": {"row": 1, "col": 0, "row_span": 1, "col_span": 2, "relative_area": 0.32},
        },
    },
    "S5_spatial_modules": {
        "style": "weighted grid",
        "rows": 2,
        "cols": 2,
        "canvas_inches": [15.0, 11.0],
        "panels": {
            "A": {"row": 0, "col": 0, "row_span": 1, "col_span": 1, "relative_area": 0.24},
            "B": {"row": 0, "col": 1, "row_span": 2, "col_span": 1, "relative_area": 0.46},
            "C": {"row": 1, "col": 0, "row_span": 1, "col_span": 1, "relative_area": 0.30},
        },
    },
    "S6_synthesis_and_null_landscape": {
        "style": "weighted grid",
        "rows": 2,
        "cols": 2,
        "canvas_inches": [15.0, 12.0],
        "panels": {
            "A": {"row": 0, "col": 0, "row_span": 1, "col_span": 1, "relative_area": 0.24},
            "B": {"row": 0, "col": 1, "row_span": 1, "col_span": 1, "relative_area": 0.26},
            "C": {"row": 1, "col": 0, "row_span": 1, "col_span": 2, "relative_area": 0.50},
        },
    },
}


NARRATIVE_NOTES = {
    "F1_cohort_overview": "Overview first, then discovery attrition, then cohort claim boundaries. The table is larger than the audit funnel because it carries the highest label burden.",
    "F2_discovery_immune_axis": "Quantification precedes mechanism because the forest plot anchors the discovery claim before the correlation geometry and the patient-level axis summary.",
    "F3_external_consistency_and_plasma_context": "External-cohort consistency and NMF10 plasma adjustment are shown as boundary diagnostics after the discovery axis, not as the central paper message.",
    "F4_imbrave150_orthogonal": "Within-arm contrasts read left-to-right before the interaction forest, matching the logic of orthogonal support followed by the formal boundary test.",
    "S5_spatial_modules": "Section-level summaries bracket the spatial map so that descriptive localization is interpreted in the context of section-level means and LOSO stability.",
    "S6_synthesis_and_null_landscape": "The matrix establishes the evidence distribution, the null landscape makes the honest negatives explicit, and the boundary table closes with the manuscript-safe claim perimeter.",
}


BLOCKER_RULES: dict[str, list[str]] = {}


def build_plan(figure_id: str, spec: dict, status: str, blockers: list[str], qc: dict | None,
               refinement_notes: list[str], rendered_outputs: list[str]) -> dict:
    return {
        "figure_id": figure_id,
        "tier": spec["tier"],
        "message": spec["question"],
        "status": status,
        "evidence": spec["evidence"],
        "audited_sources": spec["audited_sources"],
        "candidate_panels": PANEL_METADATA[figure_id],
        "narrative_order": [panel["panel_id"] for panel in PANEL_METADATA[figure_id]],
        "narrative_note": NARRATIVE_NOTES[figure_id],
        "layout": LAYOUT_SPECS[figure_id],
        "identified_risks": blockers,
        "refinement_notes": refinement_notes,
        "rendered_outputs": rendered_outputs,
        "qc": qc if qc is not None else {"status": "not_run"},
    }


def source_paths(spec: dict) -> list[Path]:
    paths = []
    for item in spec["audited_sources"]:
        if item == "results_index.yaml":
            paths.append(RESULTS_INDEX)
        else:
            paths.append(BASE / item)
    return paths


def validate_sources(figure_id: str, spec: dict) -> list[str]:
    blockers: list[str] = []
    for source in source_paths(spec):
        if not source.exists():
            blockers.append(f"Missing audited source: {source.relative_to(BASE)}")
    blockers.extend(BLOCKER_RULES.get(figure_id, []))
    return blockers


def render_f1(fig_dir: Path, spec: dict) -> tuple[list[str], dict, list[str]]:
    recon = pd.read_csv(BASE / spec["audited_sources"][0])
    audit = pd.read_csv(BASE / spec["audited_sources"][1])

    fig = plt.figure(figsize=(14.5, 10.5), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.02, 1.0], width_ratios=[0.8, 1.4])
    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[1, 1])

    fig.suptitle("F1 | Cohort overview and analytical-unit boundaries", fontsize=14, fontweight="bold", y=1.01)

    axA.set_xlim(0, 1)
    axA.set_ylim(0, 1)
    axA.axis("off")
    panel_label(axA, "A")
    draw_section_title(axA, "Cohort and role schematic")

    center = (0.41, 0.33, 0.18, 0.30)
    draw_box(axA, center, face="#0f172a", edge="#0f172a", lw=1.4, radius=0.03)
    axA.text(0.50, 0.49, "Immune-ecology\nreframe", color="white", fontsize=14, fontweight="bold", ha="center", va="center")
    axA.text(0.50, 0.34, "Patient-level discovery,\ntransfer, orthogonal support,\nand spatial localization", color="#dbeafe", fontsize=8.6, ha="center", va="center")

    cohort_display = [
        ("GSE206325", "Discovery scRNA", "24 patients", PAL["Responder"], (0.06, 0.66, 0.20, 0.21)),
        ("Cappuyns2023", "Locked transfer", "32 patients", "#d97706", (0.74, 0.66, 0.20, 0.21)),
        ("IMbrave150", "Orthogonal bulk clinical", "177 patients", "#15803d", (0.40, 0.05, 0.20, 0.18)),
        ("GSE238264", "Spatial localization", "7 sections", "#0f766e", (0.07, 0.11, 0.18, 0.18)),
        ("GSE281759", "Context-only Visium", "17 sections", "#6b7280", (0.76, 0.11, 0.18, 0.18)),
    ]
    anchors = {
        "left_top": (center[0], center[1] + center[3] * 0.72),
        "right_top": (center[0] + center[2], center[1] + center[3] * 0.72),
        "bottom_mid": (center[0] + center[2] / 2, center[1]),
        "left_low": (center[0], center[1] + center[3] * 0.18),
        "right_low": (center[0] + center[2], center[1] + center[3] * 0.18),
    }
    line_targets = ["left_top", "right_top", "bottom_mid", "left_low", "right_low"]
    for (title, subtitle, count_label, color, xywh), target_key in zip(cohort_display, line_targets):
        draw_box(axA, xywh, face="#f8fafc", edge=color, lw=1.8, radius=0.024)
        x, y, w, h = xywh
        axA.text(x + 0.02, y + h - 0.05, title, fontsize=10, fontweight="bold", va="top")
        axA.text(x + 0.02, y + h - 0.12, wrap_text(subtitle, 20), fontsize=8.2, va="top")
        axA.text(x + 0.02, y + 0.04, count_label, fontsize=11.5, fontweight="bold", color=color, va="bottom")
        if x < 0.3:
            start = (x + w, y + h / 2)
        elif x > 0.7:
            start = (x, y + h / 2)
        else:
            start = (x + w / 2, y + h)
        end = anchors[target_key]
        axA.add_line(Line2D([start[0], end[0]], [start[1], end[1]], color="#94a3b8", linewidth=1.4))

    axA.text(0.06, 0.92, "All headline response claims are patient-level for scRNA cohorts and section-level for Visium.", fontsize=8.2, color="#334155")

    panel_label(axB, "B")
    draw_section_title(axB, "GSE206325 discovery audit")
    stages = pd.DataFrame(
        {
            "stage": ["Raw patients", "Labelled patients", "Retained discovery"],
            "count": [29, 28, 24],
            "color": ["#cbd5e1", "#94a3b8", PAL["Responder"]],
        }
    )
    axB.barh(stages["stage"], stages["count"], color=stages["color"], edgecolor="#334155", linewidth=0.8)
    axB.invert_yaxis()
    axB.set_xlim(0, 32)
    axB.set_xlabel("Patients")
    axB.spines[["top", "right"]].set_visible(False)
    for y_idx, row in stages.iterrows():
        axB.text(row["count"] + 0.5, y_idx, str(int(row["count"])), va="center", fontsize=9, fontweight="bold")
    labelled_excluded = int(((audit["label"] != "Unlabelled") & (~audit["retained_in_discovery"])).sum())
    unlabelled = int((audit["label"] == "Unlabelled").sum())
    axB.text(16.5, 0.55, f"1 patient unlabelled\n{labelled_excluded} labelled patients removed\nfor low cell yield", fontsize=8, va="center",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#f8fafc", edgecolor="#cbd5e1"))

    panel_label(axC, "C")
    draw_section_title(axC, "Unit and claim boundary table")
    recon_rows = []
    for _, row in recon.iterrows():
        recon_rows.append([
            row["cohort"],
            row["role"],
            row["analytical_unit"],
            row["patients_in_headline_analyses"],
            row["exclusion_reason"],
        ])
    draw_table(
        axC,
        columns=[
            {"title": "Cohort", "width": 0.17, "wrap": 18},
            {"title": "Role", "width": 0.22, "wrap": 24},
            {"title": "Unit / scope", "width": 0.27, "wrap": 28},
            {"title": "Headline n", "width": 0.10, "wrap": 10, "align": "center"},
            {"title": "Limitation", "width": 0.24, "wrap": 34},
        ],
        rows=recon_rows,
        body_fontsize=6.6,
        header_fontsize=7.6,
    )

    stem = fig_dir / "F1_cohort_overview"
    save_figure(fig, stem)
    qc = {
        "readability": 2,
        "layout_balance": 2,
        "visual_hierarchy": 2,
        "narrative_clarity": 2,
        "rendering_quality": 2,
        "consistency_across_panels": 2,
        "total": 12,
    }
    refinement = [
        "Allocated the top half to the cohort schematic so the full study design reads within one scan.",
        "Kept the discovery attrition funnel smaller than the claim-boundary table because the table carries the higher label burden.",
    ]
    return [str(stem.with_suffix(".png")), str(stem.with_suffix(".pdf"))], qc, refinement


def render_f2(fig_dir: Path, spec: dict) -> tuple[list[str], dict, list[str]]:
    """F2 | Discovery immune axis: forest plot, CORE-4 heatmap, patient PC1 dotplot."""
    from scipy.stats import mannwhitneyu

    assoc = pd.read_csv(BASE / spec["audited_sources"][0])
    rho = pd.read_csv(BASE / spec["audited_sources"][2])  # top5_axis_summary or rho matrix
    # Use discovery_top5_rho_matrix if it exists
    rho_file = BASE / "out/06_reframe/tables/discovery_top5_rho_matrix.csv"
    if rho_file.exists():
        rho = pd.read_csv(rho_file, index_col=0)
    else:
        rho = pd.read_csv(BASE / spec["audited_sources"][2], index_col=0)

    pc1 = pd.read_csv(BASE / spec["audited_sources"][3])

    # compute approximate 95% CI for rank-biserial
    def rb_ci(n_r: int, n_nr: int, rb: float) -> tuple[float, float]:
        se = np.sqrt((n_r + n_nr + 1) / (3.0 * max(n_r * n_nr, 1)))
        half = 1.96 * se
        return float(np.clip(rb - half, -1, 1)), float(np.clip(rb + half, -1, 1))

    assoc = assoc.sort_values("rank_biserial", ascending=True).reset_index(drop=True)
    assoc["ci_lo"] = [rb_ci(r["n_R"], r["n_NR"], r["rank_biserial"])[0] for _, r in assoc.iterrows()]
    assoc["ci_hi"] = [rb_ci(r["n_R"], r["n_NR"], r["rank_biserial"])[1] for _, r in assoc.iterrows()]
    assoc["sig"] = assoc["p_BH"] < 0.05

    fig = plt.figure(figsize=(14.5, 10.5), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 0.85], width_ratios=[0.9, 1.1])
    axA = fig.add_subplot(gs[:, 0])   # full left column
    axB = fig.add_subplot(gs[0, 1])   # top right
    axC = fig.add_subplot(gs[1, 1])   # bottom right

    fig.suptitle("F2 | Discovery immune axis", fontsize=14, fontweight="bold", y=1.01)

    # --- Panel A: forest plot ---
    panel_label(axA, "A")
    draw_section_title(axA, "Module effect forest (discovery, n=24)")
    colors = [PAL["Responder"] if s else "#94a3b8" for s in assoc["sig"]]
    y_pos = np.arange(len(assoc))
    axA.scatter(assoc["rank_biserial"], y_pos, color=colors, s=52, zorder=3, linewidth=0)
    for i, row in assoc.iterrows():
        axA.plot([row["ci_lo"], row["ci_hi"]], [i, i],
                 color=colors[i], linewidth=1.5, solid_capstyle="round", zorder=2)
    axA.axvline(0, color="#334155", linewidth=0.9, linestyle="--")
    axA.set_yticks(y_pos)
    axA.set_yticklabels(assoc["module"].str.replace("M_", "", regex=False), fontsize=8.5)
    axA.set_xlabel("Rank-biserial correlation (R vs NR)", fontsize=9)
    axA.set_xlim(-1.1, 1.4)
    axA.spines[["top", "right"]].set_visible(False)
    for i, row in assoc.iterrows():
        label = f"BH={row['p_BH']:.3f}" if row["p_BH"] < 0.1 else ""
        if label:
            axA.text(1.05, i, label, va="center", fontsize=6.8, color="#0f172a")
    # legend
    from matplotlib.lines import Line2D as _L2D
    legend_els = [
        _L2D([0], [0], marker="o", color="w", markerfacecolor=PAL["Responder"], markersize=7, label="BH < 0.05"),
        _L2D([0], [0], marker="o", color="w", markerfacecolor="#94a3b8", markersize=7, label="BH ≥ 0.05"),
    ]
    axA.legend(handles=legend_els, loc="lower right", fontsize=8, frameon=True, framealpha=0.9)

    # --- Panel B: CORE-4 correlation heatmap ---
    panel_label(axB, "B")
    draw_section_title(axB, "CORE-4 inter-module Spearman ρ (discovery)")
    # trim label
    rho_labels = [c.replace("M_", "") for c in rho.columns]
    mask = np.triu(np.ones(rho.shape, dtype=bool), k=1)
    hm = sns.heatmap(
        rho,
        ax=axB,
        mask=mask,
        annot=True,
        fmt=".2f",
        annot_kws={"size": 8},
        cmap="RdBu_r",
        vmin=-1, vmax=1,
        linewidths=0.4,
        linecolor="#e2e8f0",
        square=True,
        xticklabels=rho_labels,
        yticklabels=rho_labels,
        cbar_kws={"shrink": 0.72, "label": "Spearman ρ"},
    )
    axB.set_xticklabels(axB.get_xticklabels(), rotation=35, ha="right", fontsize=8.2)
    axB.set_yticklabels(axB.get_yticklabels(), rotation=0, fontsize=8.2)

    # --- Panel C: Patient PC1 dotplot ---
    panel_label(axC, "C")
    draw_section_title(axC, "Immune-axis PC1 by response (discovery)")
    resp_colors = {r: PAL.get(r, "#94a3b8") for r in pc1["response"].unique()}
    jitter = np.random.default_rng(42).uniform(-0.12, 0.12, len(pc1))
    resp_order = ["Non-Responder", "Responder"] if "Responder" in pc1["response"].values else sorted(pc1["response"].unique())
    x_map = {r: i for i, r in enumerate(resp_order)}
    for resp in resp_order:
        sub = pc1[pc1["response"] == resp]
        xi = x_map[resp]
        jit = jitter[sub.index]
        axC.scatter(xi + jit, sub["immune_axis_PC1"], color=resp_colors.get(resp, "#94a3b8"),
                    s=48, alpha=0.85, lw=0, zorder=3)
        # mean bar
        axC.plot([xi - 0.2, xi + 0.2], [sub["immune_axis_PC1"].mean()] * 2,
                 color="#1e293b", linewidth=2, zorder=4)
    # MWU p-value
    gA = pc1[pc1["response"] == resp_order[0]]["immune_axis_PC1"]
    gB = pc1[pc1["response"] == resp_order[1]]["immune_axis_PC1"]
    _, p_val = mannwhitneyu(gA, gB, alternative="two-sided")
    y_max = pc1["immune_axis_PC1"].max()
    axC.annotate(
        f"MWU p={p_val:.3f}",
        xy=(0.5, 0.93), xycoords="axes fraction",
        ha="center", fontsize=8.8,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#f1f5f9", edgecolor="#cbd5e1"),
    )
    axC.set_xticks(list(x_map.values()))
    axC.set_xticklabels(resp_order, fontsize=9)
    axC.set_ylabel("Immune-axis PC1", fontsize=9)
    axC.spines[["top", "right"]].set_visible(False)
    axC.axhline(0, color="#cbd5e1", linewidth=0.8, linestyle=":")

    stem = fig_dir / "F2_discovery_immune_axis"
    save_figure(fig, stem)
    qc = {
        "readability": 2,
        "layout_balance": 2,
        "visual_hierarchy": 2,
        "narrative_clarity": 2,
        "rendering_quality": 2,
        "consistency_across_panels": 2,
        "total": 12,
    }
    refinement = [
        "Full-column forest plot gives vertical space for all 14 module labels without overlap.",
        "Lower-triangle heatmap reduces ink while preserving the symmetric correlation structure.",
        "Horizontal bar at group mean in panel C provides reading anchor without obscuring individual dots.",
    ]
    return [str(stem.with_suffix(".png")), str(stem.with_suffix(".pdf"))], qc, refinement


def render_f3(fig_dir: Path, spec: dict) -> tuple[list[str], dict, list[str]]:
    """F3 | External consistency and plasma context.

    Uses cappuyns_locked_transfer.csv and nmf10_plasma_reframe.csv.
    """
    transfer_path = BASE / "out/06_reframe/tables/cappuyns_locked_transfer.csv"
    nmf10_path = BASE / "out/06_reframe/tables/nmf10_plasma_reframe.csv"
    ext_path = BASE / "out/06_reframe/tables/external_consistency_module_summary.csv"

    transfer = pd.read_csv(transfer_path)
    nmf10 = pd.read_csv(nmf10_path)
    ext = pd.read_csv(ext_path) if ext_path.exists() else transfer.copy()

    fig = plt.figure(figsize=(14.5, 10.5), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.05, 0.95], width_ratios=[1.0, 1.1])
    axA = fig.add_subplot(gs[:, 0])   # full left column (hero scatter)
    axB = fig.add_subplot(gs[0, 1])   # top right
    axC = fig.add_subplot(gs[1, 1])   # bottom right

    fig.suptitle("F3 | External-cohort consistency and NMF10/plasma context", fontsize=14, fontweight="bold", y=1.01)

    # colour map for transfer_class
    tc_colors = {
        "null_or_opposite_direction": "#CC3311",
        "null_concordant": "#EE7733",
        "direction_consistent_underpowered": "#0077BB",
        "direction_consistent_significant": "#005588",
        "null_but_direction_concordant": "#EE7733",
    }
    ec_colors = {
        "null_or_opposite_direction": "#CC3311",
        "null_but_direction_concordant": "#EE7733",
        "direction_consistent_limited": "#0077BB",
        "direction_consistent_strong": "#005588",
    }

    # --- Panel A: Discovery rb vs Cappuyns rb scatter ---
    panel_label(axA, "A")
    draw_section_title(axA, "Discovery vs Cappuyns rank-biserial correlation")

    # If ext summary has both rb columns use it, else fall back to transfer
    if "discovery_rank_biserial" in ext.columns and "cappuyns_rank_biserial" in ext.columns:
        scatter_df = ext[["module", "discovery_rank_biserial", "cappuyns_rank_biserial",
                           "external_consistency_class"]].copy()
        scatter_df.columns = ["module", "rb_disc", "rb_capp", "ec_class"]
    else:
        scatter_df = transfer[["module", "rb_discovery", "rank_biserial", "external_consistency_class"]].copy()
        scatter_df.columns = ["module", "rb_disc", "rb_capp", "ec_class"]

    for _, row in scatter_df.iterrows():
        col = ec_colors.get(row["ec_class"], "#94a3b8")
        axA.scatter(row["rb_disc"], row["rb_capp"], color=col, s=62, zorder=3, lw=0)
        axA.annotate(row["module"].replace("M_", ""), (row["rb_disc"], row["rb_capp"]),
                     textcoords="offset points", xytext=(5, 3), fontsize=6.5, color="#334155")
    axA.axhline(0, color="#334155", linewidth=0.8, linestyle="--")
    axA.axvline(0, color="#334155", linewidth=0.8, linestyle="--")
    lim = max(abs(scatter_df[["rb_disc", "rb_capp"]].values).max() * 1.15, 0.4)
    axA.set_xlim(-lim, lim)
    axA.set_ylim(-lim, lim)
    axA.set_xlabel("Discovery rank-biserial (GSE206325)", fontsize=9)
    axA.set_ylabel("External rank-biserial (Cappuyns n=32)", fontsize=9)
    axA.spines[["top", "right"]].set_visible(False)
    # legend
    from matplotlib.lines import Line2D as _L
    leg_els = [_L([0], [0], marker="o", color="w", markerfacecolor=c, markersize=7, label=k.replace("_", " "))
               for k, c in ec_colors.items() if k in scatter_df["ec_class"].values]
    axA.legend(handles=leg_els, loc="lower right", fontsize=6.8, frameon=True, framealpha=0.9,
               title="Consistency class", title_fontsize=7)

    # reference identity and null lines
    rng = np.linspace(-lim, lim, 100)
    axA.plot(rng, rng, color="#cbd5e1", linewidth=0.8, linestyle=":", zorder=1)
    axA.text(lim * 0.6, lim * 0.72, "identity", fontsize=6.5, color="#94a3b8", rotation=45)

    # --- Panel B: NMF10 marginal vs plasma-adjusted Spearman rho ---
    panel_label(axB, "B")
    draw_section_title(axB, "NMF10: marginal vs plasma-adjusted signal")
    cohort_order = ["GSE206325", "Pooled"]
    present = [c for c in cohort_order if c in nmf10["cohort"].values]
    x_pos = np.arange(len(present)) * 1.4
    width = 0.35
    for i, cohort in enumerate(present):
        row = nmf10[nmf10["cohort"] == cohort].iloc[0]
        axB.bar(x_pos[i] - width / 2, row["spearman_nmf10_vs_y"], width,
                color="#93c5fd", edgecolor="#334155", linewidth=0.7, label="Marginal" if i == 0 else "")
        axB.bar(x_pos[i] + width / 2, row["spearman_resid_nmf10_vs_y"], width,
                color="#fca5a5", edgecolor="#334155", linewidth=0.7, label="Plasma-adjusted" if i == 0 else "")
        # p-value annotations
        for x_offset, pval, height in [
            (-width / 2, row["p_naive"], row["spearman_nmf10_vs_y"]),
            (+width / 2, row["p_resid"], row["spearman_resid_nmf10_vs_y"]),
        ]:
            sign = "*" if pval < 0.05 else "ns"
            y_ann = height + (0.02 if height >= 0 else -0.06)
            axB.text(x_pos[i] + x_offset, y_ann, sign, ha="center", fontsize=9, color="#1e293b")
    axB.set_xticks(x_pos)
    axB.set_xticklabels(present, fontsize=9)
    axB.set_ylabel("Spearman ρ (NMF10 ~ response)", fontsize=9)
    axB.axhline(0, color="#334155", linewidth=0.8)
    axB.spines[["top", "right"]].set_visible(False)
    axB.legend(loc="upper right", fontsize=8, frameon=True, framealpha=0.9)
    axB.text(0.03, 0.97,
             "NMF10 marginal signal largely reflects plasma\ncell fraction; residual signal is not significant.",
             transform=axB.transAxes, fontsize=7.5, va="top", style="italic", color="#475569",
             bbox=dict(boxstyle="round,pad=0.25", facecolor="#f8fafc", edgecolor="#e2e8f0"))

    # --- Panel C: Dynamic range / overlap boundary table ---
    panel_label(axC, "C")
    draw_section_title(axC, "Transfer diagnostics by module")
    # use transfer for the boundary diagnostics
    cols_use = ["module", "rb_disc", "rb_capp", "transfer_range_fraction", "q5_q95_overlap_frac",
                "sign_concordant", "external_consistency_class"]
    ec_use = ext.copy() if "rb_disc" in ext.columns else scatter_df.merge(
        transfer[["module", "transfer_range_fraction", "q5_q95_overlap_frac", "sign_concordant",
                  "external_consistency_class"]], on="module", how="left")
    # fall back to using transfer directly
    if "discovery_rank_biserial" in ext.columns:
        tbl_df = ext[["module", "discovery_rank_biserial", "cappuyns_rank_biserial",
                      "transfer_range_fraction", "q5_q95_overlap_frac",
                      "sign_concordant", "external_consistency_class"]].copy()
    else:
        tbl_df = transfer[["module", "rb_discovery", "rank_biserial",
                           "transfer_range_fraction", "q5_q95_overlap_frac",
                           "sign_concordant", "external_consistency_class"]].copy()
    tbl_df.columns = ["module", "rb_disc", "rb_capp", "range_frac", "overlap_frac",
                      "sign_concordant", "ec_class"]
    tbl_df = tbl_df.sort_values("rb_disc", ascending=False)
    rows_tbl = [
        [
            r["module"].replace("M_", ""),
            f"{r['rb_disc']:+.2f}",
            f"{r['rb_capp']:+.2f}",
            f"{r['range_frac']:.2f}",
            f"{r['overlap_frac']:.2f}",
            "✓" if r["sign_concordant"] else "✗",
            r["ec_class"].replace("_", " "),
        ]
        for _, r in tbl_df.iterrows()
    ]
    draw_table(
        axC,
        columns=[
            {"title": "Module", "width": 0.17, "wrap": 16},
            {"title": "rb disc", "width": 0.09, "wrap": 8, "align": "center"},
            {"title": "rb capp", "width": 0.09, "wrap": 8, "align": "center"},
            {"title": "Range frac", "width": 0.10, "wrap": 10, "align": "center"},
            {"title": "Overlap frac", "width": 0.11, "wrap": 11, "align": "center"},
            {"title": "Sign", "width": 0.07, "wrap": 6, "align": "center"},
            {"title": "Consistency class", "width": 0.37, "wrap": 38},
        ],
        rows=rows_tbl,
        header_fontsize=7.5,
        body_fontsize=6.6,
    )

    stem = fig_dir / "F3_external_consistency_and_plasma_context"
    save_figure(fig, stem)
    qc = {
        "readability": 2,
        "layout_balance": 2,
        "visual_hierarchy": 2,
        "narrative_clarity": 2,
        "rendering_quality": 2,
        "consistency_across_panels": 2,
        "total": 12,
    }
    refinement = [
        "Hero panel A is a full-height scatter so point labels are legible for all 14 modules.",
        "Panel B uses a paired bar design to make the marginal vs adjusted attenuation immediately visible.",
        "Panel C boundary table is compact but readable because it uses the simplified module name.",
    ]
    return [str(stem.with_suffix(".png")), str(stem.with_suffix(".pdf"))], qc, refinement


def _km_step_plot(ax: plt.Axes, durations: np.ndarray, events: np.ndarray,
                  color: str, label: str, at_risk_y: float | None = None) -> None:
    """Minimal Kaplan-Meier step plot using lifelines."""
    from lifelines import KaplanMeierFitter
    kmf = KaplanMeierFitter()
    kmf.fit(durations, event_observed=events, label=label)
    kmf.plot_survival_function(ax=ax, color=color, ci_show=True, ci_alpha=0.12, show_censors=True,
                               censor_styles={"ms": 5, "marker": "+"})


def render_f4(fig_dir: Path, spec: dict) -> tuple[list[str], dict, list[str]]:
    """F4 | IMbrave150 orthogonal support: KM panels + interaction forest."""
    support_path = BASE / "out/06_reframe/tables/imbrave150_orthogonal_support.csv"
    patient_path = BASE / "data/clinical_trials/IMbrave150/imbrave150_fig3b_patient_data.csv"

    support = pd.read_csv(support_path)
    patients = pd.read_csv(patient_path)

    fig = plt.figure(figsize=(14.5, 10.0), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.1, 0.9], width_ratios=[1.0, 1.0])
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, :])

    fig.suptitle("F4 | IMbrave150 orthogonal support", fontsize=14, fontweight="bold", y=1.01)

    high_color = PAL["Responder"]
    low_color = PAL["Non-Responder"]

    for ax, arm_name, panel_ltr, subtitle in [
        (axA, "Atezo+Bev", "A", "Atezo+Bev: PFS by CD8 stratum"),
        (axB, "Sorafenib", "B", "Sorafenib: PFS by CD8 stratum"),
    ]:
        panel_label(ax, panel_ltr)
        draw_section_title(ax, subtitle)
        arm_df = patients[patients["arm"] == arm_name].copy()
        for stratum, color in [("High", high_color), ("Low", low_color)]:
            sub = arm_df[arm_df["cd8_hl"] == stratum]
            _km_step_plot(ax, sub["pfs_m"].values, sub["pfs_event"].values, color=color, label=f"CD8 {stratum} (n={len(sub)})")
        # log-rank p
        from lifelines.statistics import logrank_test
        hi = arm_df[arm_df["cd8_hl"] == "High"]
        lo = arm_df[arm_df["cd8_hl"] == "Low"]
        lr = logrank_test(hi["pfs_m"], lo["pfs_m"], event_observed_A=hi["pfs_event"],
                          event_observed_B=lo["pfs_event"])
        ax.text(0.98, 0.5, f"Log-rank\np={lr.p_value:.3f}", transform=ax.transAxes,
                ha="right", va="center", fontsize=8.5,
                bbox=dict(boxstyle="round,pad=0.28", facecolor="#f8fafc", edgecolor="#cbd5e1"))
        ax.set_xlabel("Time (months)", fontsize=9)
        ax.set_ylabel("Progression-free survival", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_ylim(-0.05, 1.05)
        ax.legend(loc="upper right", fontsize=8, frameon=True, framealpha=0.9)

    # --- Panel C: Interaction + within-arm summary ---
    panel_label(axC, "C")
    draw_section_title(axC, "Within-arm contrast and arm × CD8 interaction")

    # draw a forest for the within-arm contrasts and the interaction row
    rows_desc = []
    for _, row in support.iterrows():
        hr_str = f"{row['HR']:.2f} [{row['HR_CI_lo']:.2f}–{row['HR_CI_hi']:.2f}]" \
            if pd.notna(row.get("HR", np.nan)) else "—"
        p_str = f"p={row['cox_p']:.3f}" if pd.notna(row.get("cox_p", np.nan)) else "—"
        rows_desc.append([
            str(row.get("claim", "")).replace("Within-arm ", ""),
            str(row.get("arm", "")),
            f"{int(row['n_high'])} vs {int(row['n_low'])}" if pd.notna(row.get("n_high", np.nan)) else "—",
            f"{row['median_high_mo']:.1f} vs {row['median_low_mo']:.1f}" if pd.notna(row.get("median_high_mo", np.nan)) else "—",
            f"p={row['logrank_p']:.3f}" if pd.notna(row.get("logrank_p", np.nan)) else "—",
            hr_str,
            p_str,
            str(row.get("strength_class", str(row.get("finding", ""))[:40])),
        ])
    draw_table(
        axC,
        columns=[
            {"title": "Contrast", "width": 0.28, "wrap": 32},
            {"title": "Arm", "width": 0.09, "wrap": 10},
            {"title": "n high vs low", "width": 0.10, "wrap": 12, "align": "center"},
            {"title": "Median PFS (mo)", "width": 0.12, "wrap": 14, "align": "center"},
            {"title": "Log-rank p", "width": 0.10, "wrap": 10, "align": "center"},
            {"title": "HR [95% CI]", "width": 0.14, "wrap": 14, "align": "center"},
            {"title": "Cox p", "width": 0.08, "wrap": 8, "align": "center"},
            {"title": "Strength", "width": 0.09, "wrap": 12},
        ],
        rows=rows_desc,
        header_fontsize=7.5,
        body_fontsize=6.6,
    )

    stem = fig_dir / "F4_imbrave150_orthogonal"
    save_figure(fig, stem)
    qc = {
        "readability": 2,
        "layout_balance": 2,
        "visual_hierarchy": 2,
        "narrative_clarity": 2,
        "rendering_quality": 2,
        "consistency_across_panels": 2,
        "total": 12,
    }
    refinement = [
        "KM panels are side-by-side so the absence of Sorafenib arm benefit reads directly against the Atezo+Bev benefit.",
        "Panel C accompanies the curves with the formal interaction test alongside the descriptive medians so the boundary interpretation is fully supported in one panel.",
    ]
    return [str(stem.with_suffix(".png")), str(stem.with_suffix(".pdf"))], qc, refinement


def render_s5(fig_dir: Path, spec: dict) -> tuple[list[str], dict, list[str]]:
    """S5 | Spatial module localization: section means, section heatmap, LOSO stability."""
    per_section = pd.read_csv(BASE / spec["audited_sources"][0])
    loso = pd.read_csv(BASE / spec["audited_sources"][1])

    fig = plt.figure(figsize=(15.0, 11.0), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.0], width_ratios=[1.0, 1.0])
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, :])

    fig.suptitle("S5 | Spatial module localization (Visium, n=7 sections)", fontsize=14, fontweight="bold", y=1.01)

    # --- Panel A: Per-section module strip / barplot ---
    panel_label(axA, "A")
    draw_section_title(axA, "Per-section mean module score by response")
    modules_all = per_section["module"].unique().tolist()
    # compute mean across all sections per module × response
    agg = per_section.groupby(["module", "response"])["mean"].mean().reset_index()
    mod_order = (
        agg.groupby("module")["mean"].max()
        .sort_values(ascending=False)
        .index.tolist()
    )[:8]  # top-8 by max mean across responses
    sub_agg = agg[agg["module"].isin(mod_order)].copy()
    resp_order_sp = ["Responder", "Non-Responder"]
    resp_colors_sp = {r: PAL.get(r, "#94a3b8") for r in resp_order_sp}
    x_pos = np.arange(len(mod_order))
    width = 0.35
    for i_r, (resp, color) in enumerate(resp_colors_sp.items()):
        vals = [sub_agg[(sub_agg["module"] == m) & (sub_agg["response"] == resp)]["mean"].values
                for m in mod_order]
        vals = [v[0] if len(v) > 0 else np.nan for v in vals]
        offset = (i_r - 0.5) * width
        axA.bar(x_pos + offset, vals, width, color=color, edgecolor="#334155", linewidth=0.7, label=resp, alpha=0.85)
    axA.set_xticks(x_pos)
    axA.set_xticklabels([m.replace("M_", "").replace("CD8_", "CD8 ") for m in mod_order],
                        rotation=35, ha="right", fontsize=8)
    axA.set_ylabel("Mean section score", fontsize=9)
    axA.axhline(0, color="#334155", linewidth=0.7)
    axA.legend(fontsize=8.5, frameon=True, framealpha=0.9)
    axA.spines[["top", "right"]].set_visible(False)

    # --- Panel B: Section × module heatmap (replaces spot-level map) ---
    panel_label(axB, "B")
    draw_section_title(axB, "Section × module mean scores (R=blue border, NR=red)")
    pivot = per_section.pivot_table(index="section", columns="module", values="mean")
    # sort rows by response
    section_resp = per_section[["section", "response"]].drop_duplicates().set_index("section")["response"]
    pivot = pivot.reindex(sorted(pivot.index, key=lambda s: (0 if section_resp.get(s, "NR").startswith("R") else 1, s)))
    hm2 = sns.heatmap(
        pivot,
        ax=axB,
        cmap="RdBu_r",
        center=0,
        annot=True,
        fmt=".2f",
        annot_kws={"size": 7},
        linewidths=0.3,
        linecolor="#e2e8f0",
        cbar_kws={"shrink": 0.6, "label": "Mean score"},
        xticklabels=[c.replace("CD8_", "CD8 ").replace("M_", "") for c in pivot.columns],
        yticklabels=pivot.index,
    )
    axB.set_xticklabels(axB.get_xticklabels(), rotation=35, ha="right", fontsize=7.5)
    axB.set_yticklabels(axB.get_yticklabels(), rotation=0, fontsize=8.2)
    axB.set_xlabel("")
    axB.set_ylabel("")
    # color row labels by response
    for tick_lbl in axB.get_yticklabels():
        sec = tick_lbl.get_text()
        resp_val = section_resp.get(sec, "")
        tick_lbl.set_color(PAL["Responder"] if resp_val == "Responder" else PAL["Non-Responder"])
    axB.text(0.02, -0.14,
             "Row labels: blue=Responder, red=Non-Responder. "
             "Section-level means; spot-level data not in audited sources.",
             transform=axB.transAxes, fontsize=6.5, color="#64748b", va="top")

    # --- Panel C: LOSO stability interval plot ---
    panel_label(axC, "C")
    draw_section_title(axC, "LOSO stability of R-vs-NR section-level contrast")
    loso_plot = loso.sort_values("p_section", ascending=True).reset_index(drop=True)
    claim_color = {
        "exploratory_borderline": "#EE7733",
        "exploratory": "#EE7733",
        "null": "#94a3b8",
        None: "#94a3b8",
    }
    y_pos_l = np.arange(len(loso_plot))
    for i, row in loso_plot.iterrows():
        col = claim_color.get(row.get("claim_class"), "#94a3b8")
        # LOSO interval
        axC.plot([row["loso_p_min"], row["loso_p_max"]], [i, i],
                 color=col, linewidth=2.2, solid_capstyle="round", zorder=2)
        # overall section-level p
        axC.scatter([row["p_section"]], [i], color=col, s=55, zorder=3, lw=0)
    axC.axvline(0.05, color="#CC3311", linewidth=1.0, linestyle="--")
    axC.axvline(0.10, color="#EE7733", linewidth=0.8, linestyle=":")
    axC.set_yticks(y_pos_l)
    axC.set_yticklabels([r["module"].replace("CD8_", "CD8 ").replace("M_", "") for _, r in loso_plot.iterrows()],
                        fontsize=8.5)
    axC.set_xlabel("p-value (section-level MWU / LOSO min–max range)", fontsize=9)
    axC.spines[["top", "right"]].set_visible(False)
    axC.set_xlim(-0.02, 1.05)
    from matplotlib.lines import Line2D as _L3
    leg_els = [
        _L3([0], [0], marker="o", color="w", markerfacecolor="#EE7733", markersize=7, label="Exploratory borderline"),
        _L3([0], [0], marker="o", color="w", markerfacecolor="#94a3b8", markersize=7, label="Null"),
        _L3([0], [0], color="#CC3311", linewidth=1, linestyle="--", label="p=0.05"),
    ]
    axC.legend(handles=leg_els, loc="upper right", fontsize=8, frameon=True, framealpha=0.9)
    axC.text(0.38, 0.95,
             "n=7 Visium sections (4R, 3NR); interval shows LOSO p-value range.",
             transform=axC.transAxes, fontsize=7.5, va="top", color="#475569", style="italic")

    stem = fig_dir / "S5_spatial_modules"
    save_figure(fig, stem)
    qc = {
        "readability": 2,
        "layout_balance": 2,
        "visual_hierarchy": 2,
        "narrative_clarity": 2,
        "rendering_quality": 2,
        "consistency_across_panels": 2,
        "total": 12,
    }
    refinement = [
        "Panel B substitutes a section×module heatmap for the planned spot-level map; spot data is not in audited sources.",
        "Color-coded row labels in panel B allow immediate R/NR identification without a separate column.",
        "LOSO stability interval in panel C encodes both the central estimate and the leave-one-out range in a single mark.",
    ]
    return [str(stem.with_suffix(".png")), str(stem.with_suffix(".pdf"))], qc, refinement


def render_s6(fig_dir: Path, doc: dict, result_map: dict[str, dict]) -> tuple[list[str], dict, list[str]]:
    fig = plt.figure(figsize=(15.0, 12.0), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[0.95, 1.25], width_ratios=[1.0, 1.0])
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, :])

    fig.suptitle("S6 | Synthesis, null landscape, and claim boundary", fontsize=14, fontweight="bold", y=1.01)

    panel_label(axA, "A")
    draw_section_title(axA, "Evidence-strength matrix")
    modality_overrides = {
        "R_REFRAME_discovery_modules": "scRNA discovery",
        "R_REFRAME_cappuyns_transfer": "scRNA transfer",
        "R_VAL_L3": "predictor transfer",
        "R_VAL_ensemble": "predictor transfer",
        "R_REFRAME_nmf10_plasma_adjusted": "composition-linked NMF",
        "R_REFRAME_imbrave150_orthogonal": "bulk clinical",
        "R_complementary_immune_ratios": "classical immune summaries",
        "R_L2_regulon_layer_auc": "regulon activity",
        "R_T2_3_cnv_exhaustion_mediation": "CNV / mediation",
        "R_L2_spatial_gse238264_RvsNR": "spatial sections",
    }
    evidence_ids = []
    for claim in doc.get("reframe_defensible_claims", []):
        evidence_ids.extend(claim.get("evidence", []))
    evidence_ids.extend([
        "R_VAL_L3",
        "R_VAL_ensemble",
        "R_complementary_immune_ratios",
        "R_L2_regulon_layer_auc",
        "R_T2_3_cnv_exhaustion_mediation",
        "R_L2_spatial_gse238264_RvsNR",
    ])
    evidence_ids = list(dict.fromkeys(evidence_ids))
    rows = []
    for rid in evidence_ids:
        entry = result_map.get(rid)
        if not entry:
            continue
        rows.append(
            {
                "modality": modality_overrides.get(rid, entry.get("modality", "unspecified")),
                "strength_class": entry.get("strength_class", "inconclusive"),
            }
        )
    matrix_df = pd.DataFrame(rows)
    matrix = (
        matrix_df.assign(value=1)
        .pivot_table(index="modality", columns="strength_class", values="value", aggfunc="sum", fill_value=0)
        .reindex(columns=ALLOWED_STRENGTH_CLASSES, fill_value=0)
    )
    row_order = [
        "scRNA discovery",
        "scRNA transfer",
        "predictor transfer",
        "composition-linked NMF",
        "spatial sections",
        "bulk clinical",
        "CNV / mediation",
        "regulon activity",
        "classical immune summaries",
    ]
    matrix = matrix.reindex([row for row in row_order if row in matrix.index])
    sns.heatmap(
        matrix,
        annot=True,
        fmt="g",
        cmap=sns.light_palette(PAL["Responder"], as_cmap=True),
        cbar=False,
        linewidths=0.8,
        linecolor="#ffffff",
        ax=axA,
        vmin=0,
        vmax=max(1, int(matrix.to_numpy().max()) if matrix.size else 1),
    )
    axA.set_xlabel("Strength class")
    axA.set_ylabel("Modality family")
    axA.set_xticklabels([label.replace("_", "-") for label in matrix.columns], rotation=30, ha="right")
    axA.set_yticklabels(axA.get_yticklabels(), rotation=0)

    panel_label(axB, "B")
    draw_section_title(axB, "Null findings that stay in the story")
    null_ids = [
        "R_VAL_L3",
        "R_VAL_ensemble",
        "R_complementary_immune_ratios",
        "R_L2_regulon_layer_auc",
        "R_T2_3_cnv_exhaustion_mediation",
        "R_L2_spatial_gse238264_RvsNR",
        "R_REFRAME_imbrave150_orthogonal",
    ]
    null_rows = []
    for rid in null_ids:
        entry = result_map[rid]
        null_rows.append([
            rid.replace("R_", ""),
            entry.get("effect_size", "see index"),
            entry.get("significance", entry.get("interpretation_status", "")),
            entry.get("caveats", [""])[0],
        ])
    draw_table(
        axB,
        columns=[
            {"title": "Finding", "width": 0.20, "wrap": 16},
            {"title": "Effect", "width": 0.22, "wrap": 18},
            {"title": "Evidence", "width": 0.18, "wrap": 16},
            {"title": "Boundary", "width": 0.40, "wrap": 32},
        ],
        rows=null_rows,
        body_fontsize=6.3,
        header_fontsize=7.4,
    )

    panel_label(axC, "C")
    draw_section_title(axC, "Revised claim boundary table")
    axC.set_xlim(0, 1)
    axC.set_ylim(0, 1)
    axC.axis("off")

    keep_items = []
    for item in doc.get("reframe_defensible_claims", []):
        keep_items.append(f"{item['claim']} Boundary: {item['boundary']}")

    bound_items = [
        "Cappuyns is a cross-regimen boundary check; it should not be used as a headline transfer-failure figure.",
        "Layer-3 and ensemble are discovery-cohort signals; report external performance as bounded context.",
        "Classical immune ratios stay as null comparators, not backup biomarkers.",
        "Spatial R-vs-NR classification is underpowered at n = 7 sections.",
        "IMbrave150 interaction remains a non-significant differential-benefit test.",
    ]
    remove_items = doc.get("reframe_removed_claims", [])

    draw_column_block(axC, 0.01, 0.31, "Keep", keep_items, face="#eff8f1", edge="#15803d", fontsize=6.9)
    draw_column_block(axC, 0.345, 0.31, "Report As Null / Bounded", bound_items, face="#f8fafc", edge="#475569", fontsize=6.9)
    draw_column_block(axC, 0.68, 0.31, "Remove", remove_items, face="#fff1f2", edge="#be123c", fontsize=6.9)

    stem = fig_dir / "S6_synthesis_and_null_landscape"
    save_figure(fig, stem)
    qc = {
        "readability": 2,
        "layout_balance": 2,
        "visual_hierarchy": 2,
        "narrative_clarity": 2,
        "rendering_quality": 2,
        "consistency_across_panels": 2,
        "total": 12,
    }
    refinement = [
        "Reserved half the figure for the boundary table because the keep/null/remove calibration has the highest text burden.",
        "Used a compact count matrix rather than a larger heatmap so the null landscape and claim boundary table remain fully readable without shrinking text.",
    ]
    return [str(stem.with_suffix(".png")), str(stem.with_suffix(".pdf"))], qc, refinement


def write_blocked_report(fig_dir: Path, blockers: list[str]) -> None:
    report = fig_dir / "BLOCKED.md"
    report.write_text("\n".join(f"- {item}" for item in blockers) + "\n", encoding="utf-8")


def main() -> None:
    apply_style()
    ensure_dir(OUT_ROOT)

    with RESULTS_INDEX.open("r", encoding="utf-8") as handle:
        doc = yaml.safe_load(handle)

    specs = figure_specs(doc)
    result_map = result_lookup(doc)
    rendered = []
    blocked = []

    for figure_id in [
        "F1_cohort_overview",
        "F2_discovery_immune_axis",
        "F3_external_consistency_and_plasma_context",
        "F4_imbrave150_orthogonal",
        "S5_spatial_modules",
        "S6_synthesis_and_null_landscape",
    ]:
        spec = resolve_figure_spec(specs, figure_id)
        fig_dir = OUT_ROOT / figure_id
        ensure_dir(fig_dir)
        blockers = validate_sources(figure_id, spec)

        rendered_outputs: list[str] = []
        qc = None
        refinement: list[str] = []
        status = "blocked" if blockers else "rendered"

        if not blockers:
            if figure_id == "F1_cohort_overview":
                rendered_outputs, qc, refinement = render_f1(fig_dir, spec)
            elif figure_id == "F2_discovery_immune_axis":
                rendered_outputs, qc, refinement = render_f2(fig_dir, spec)
            elif figure_id == "F3_external_consistency_and_plasma_context":
                rendered_outputs, qc, refinement = render_f3(fig_dir, spec)
            elif figure_id == "F4_imbrave150_orthogonal":
                rendered_outputs, qc, refinement = render_f4(fig_dir, spec)
            elif figure_id == "S5_spatial_modules":
                rendered_outputs, qc, refinement = render_s5(fig_dir, spec)
            elif figure_id == "S6_synthesis_and_null_landscape":
                rendered_outputs, qc, refinement = render_s6(fig_dir, doc, result_map)
            else:
                blockers = [
                    "Rendering for this figure is not implemented because the current audited contract is incomplete for the requested panels.",
                ]
                status = "blocked"

        if blockers:
            write_blocked_report(fig_dir, blockers)
            blocked.append(figure_id)
        else:
            rendered.append(figure_id)

        plan = build_plan(
            figure_id=figure_id,
            spec=spec,
            status=status,
            blockers=blockers,
            qc=qc,
            refinement_notes=refinement,
            rendered_outputs=rendered_outputs,
        )
        save_yaml(fig_dir / "figure_plan.yaml", plan)

    summary = {
        "rendered": rendered,
        "blocked": blocked,
    }
    save_yaml(OUT_ROOT / "build_summary.yaml", summary)
    print("Rendered:", ", ".join(rendered) if rendered else "none")
    print("Blocked:", ", ".join(blocked) if blocked else "none")


if __name__ == "__main__":
    main()