#!/usr/bin/env python3
"""Composite figure QC for the HCC immune-spatial reframe.

This script performs lightweight, reproducible checks that can run after any
figure rebuild. It does not certify biological claims; it records technical
rendering/source checks and banned narrative phrases so the figure set can be
reviewed by a human/auditor.

Outputs
-------
out/06_manuscript_figures/figure_qc_report.yaml
out/06_manuscript_figures/panel_metadata.csv
out/06_manuscript_figures/banned_phrases_audit.csv
"""
from __future__ import annotations

from pathlib import Path
import re
import sys
from typing import Any

import pandas as pd
import yaml

BASE = Path(__file__).resolve().parents[2]
FIG_ROOT = BASE / "out" / "06_manuscript_figures"
REPORT = FIG_ROOT / "figure_qc_report.yaml"
PANEL_META = FIG_ROOT / "panel_metadata.csv"
BANNED = FIG_ROOT / "banned_phrases_audit.csv"

try:
    from PIL import Image
except Exception:  # pragma: no cover - optional dependency
    Image = None

BANNED_PHRASES = {
    "transfer failure": "bounded external consistency / cohort-regimen boundary",
    "universally portable": "bounded to the tested cohort/regimen",
    "transfer-validated": "externally contextualized",
    "robust multi-cohort signature": "discovery signature with external boundary checks",
    "best section": "representative or median section",
    "treatment-predictive biomarker": "response-associated feature or within-arm stratifier",
    "NMF10 responder-programme": "NMF10 plasma/TLS programme",
}

TEXT_PATTERNS = [
    "results_index.yaml",
    "README.md",
    "docs/**/*.md",
    "docs/**/*.py",
    "scripts/**/*.py",
]


def load_yaml(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def relative(path: Path) -> str:
    try:
        return str(path.relative_to(BASE))
    except ValueError:
        return str(path)


def collect_figure_plans() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    figure_rows: list[dict[str, Any]] = []
    panel_rows: list[dict[str, Any]] = []
    for plan_path in sorted(FIG_ROOT.glob("**/figure_plan.yaml")):
        try:
            plan = load_yaml(plan_path)
        except Exception as exc:
            figure_rows.append({
                "figure_id": plan_path.parent.name,
                "plan_path": relative(plan_path),
                "status": "plan_parse_failed",
                "error": str(exc),
            })
            continue

        figure_id = plan.get("figure_id", plan_path.parent.name)
        qc = plan.get("qc", {}) or {}
        outputs = plan.get("rendered_outputs", []) or []
        blockers = plan.get("identified_risks", []) or []
        score = qc.get("total")
        figure_rows.append({
            "figure_id": figure_id,
            "plan_path": relative(plan_path),
            "status": plan.get("status", "unknown"),
            "qc_total": score,
            "n_panels": len(plan.get("candidate_panels", []) or []),
            "n_outputs": len(outputs),
            "n_blockers": len(blockers),
            "message": plan.get("message", ""),
        })
        for panel in plan.get("candidate_panels", []) or []:
            row = dict(panel)
            row["figure_id"] = figure_id
            row["plan_path"] = relative(plan_path)
            panel_rows.append(row)
    return figure_rows, panel_rows


def image_diagnostics() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for image_path in sorted(FIG_ROOT.glob("**/*.png")):
        if Image is None:
            rows.append({"file": relative(image_path), "status": "unchecked_no_pillow"})
            continue
        try:
            with Image.open(image_path) as im:
                dpi = im.info.get("dpi", (None, None))
                width, height = im.size
            min_dpi = min([d for d in dpi if isinstance(d, (int, float))] or [0])
            rows.append({
                "file": relative(image_path),
                "status": "pass" if min_dpi >= 600 else "warn_low_or_missing_dpi",
                "width_px": int(width),
                "height_px": int(height),
                "dpi_x": dpi[0],
                "dpi_y": dpi[1],
                "megapixels": round(width * height / 1e6, 2),
            })
        except Exception as exc:
            rows.append({"file": relative(image_path), "status": "image_read_failed", "error": str(exc)})
    return rows


def iter_text_files() -> list[Path]:
    files: list[Path] = []
    for pattern in TEXT_PATTERNS:
        files.extend(BASE.glob(pattern))
    # Exclude generated binary-ish directories and archives if present.
    unique = []
    seen = set()
    for path in files:
        if not path.is_file() or path in seen:
            continue
        if path.name == "06_qc_composite_figures.py":
            continue
        if "/out/_archive/" in str(path):
            continue
        seen.add(path)
        unique.append(path)
    return sorted(unique)


def banned_phrase_audit() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    files = iter_text_files()
    for phrase, replacement in BANNED_PHRASES.items():
        pattern = re.compile(re.escape(phrase), flags=re.IGNORECASE)
        hits = []
        for path in files:
            try:
                text = path.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                continue
            count = len(pattern.findall(text))
            if count:
                hits.append(f"{relative(path)}:{count}")
        rows.append({
            "phrase": phrase,
            "replacement": replacement,
            "count": sum(int(item.rsplit(":", 1)[1]) for item in hits),
            "locations": "; ".join(hits),
            "approved": "no" if hits else "yes",
        })
    return rows


def score_overall(figure_rows: list[dict[str, Any]], image_rows: list[dict[str, Any]], banned_rows: list[dict[str, Any]]) -> dict[str, Any]:
    blocked = [r for r in figure_rows if r.get("status") == "blocked"]
    low_qc = [r for r in figure_rows if isinstance(r.get("qc_total"), (int, float)) and r["qc_total"] < 9]
    image_warnings = [r for r in image_rows if not str(r.get("status", "")).startswith("pass")]
    banned_hits = [r for r in banned_rows if r.get("count", 0) > 0]
    status = "pass"
    if blocked or low_qc or banned_hits:
        status = "needs_review"
    return {
        "status": status,
        "n_figure_plans": len(figure_rows),
        "n_blocked_figures": len(blocked),
        "n_low_qc_figures": len(low_qc),
        "n_png_images": len(image_rows),
        "n_image_warnings": len(image_warnings),
        "n_banned_phrase_classes": len(banned_hits),
    }


def main() -> int:
    FIG_ROOT.mkdir(parents=True, exist_ok=True)
    figure_rows, panel_rows = collect_figure_plans()
    image_rows = image_diagnostics()
    banned_rows = banned_phrase_audit()

    pd.DataFrame(panel_rows).to_csv(PANEL_META, index=False)
    pd.DataFrame(banned_rows).to_csv(BANNED, index=False)

    report = {
        "qc_scope": relative(FIG_ROOT),
        "overall": score_overall(figure_rows, image_rows, banned_rows),
        "figures": figure_rows,
        "images": image_rows,
        "outputs": {
            "panel_metadata": relative(PANEL_META),
            "banned_phrases_audit": relative(BANNED),
            "report": relative(REPORT),
        },
        "notes": [
            "This is a technical QC report; biological claim validation remains governed by results_index.yaml and the analysis audit.",
            "Blocked figures or banned phrase hits should be resolved before final figure handoff.",
        ],
    }
    REPORT.write_text(yaml.safe_dump(report, sort_keys=False, allow_unicode=True), encoding="utf-8")
    print(f"Figure QC report written: {REPORT}")
    print(yaml.safe_dump(report["overall"], sort_keys=False))
    return 0 if report["overall"]["status"] == "pass" else 2


if __name__ == "__main__":
    sys.exit(main())
