#!/usr/bin/env python3
"""Post-run audit for the HCC immune-spatial biology reframe pipeline.

This script is intentionally read-only for analysis outputs. It summarizes SLURM
job states, expected tables/figures, figure QC, banned-phrase status, spatial
inferCNV confidence fields, and results_index.yaml parseability.
"""
from __future__ import annotations

import csv
import datetime as dt
import os
import re
import subprocess
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[2]
AUDIT_DIR = ROOT / "out" / "00_pipeline_audit"
AUDIT_DIR.mkdir(parents=True, exist_ok=True)

JOB_IDS = [j for j in os.environ.get("PIPELINE_JOB_IDS", "").replace(";", ",").split(",") if j.strip()]
NOW = dt.datetime.now().isoformat(timespec="seconds")

KEY_OUTPUTS = {
    "integrated_atlas": ROOT / "data/processed/integrated_io_atlas.h5ad",
    "integrated_cnv_atlas": ROOT / "data/processed/integrated_io_atlas_cnv.h5ad",
    "response_tme_summary": ROOT / "out/02a_response_tme/tables/celltype_response_stats.csv",
    "spatial_cnv_run_log": ROOT / "out/02b_spatial_biomarkers/tables/spatial_infercnvpy_run_log.csv",
    "spatial_cnv_spot_table": ROOT / "out/02b_spatial_biomarkers/tables/spatial_infercnvpy_spot.csv",
    "multilayer_features": ROOT / "out/05_multilayer/tables/patient_multilayer_features.csv",
    "nmf_program_stats": ROOT / "out/05_multilayer/tables/nmf_program_response_stats.csv",
    "nmf10_plasma_diag": ROOT / "out/05_multilayer/figures/05i_nmf10_plasma_diagnostic.png",
    "external_consistency_summary": ROOT / "out/06_reframe/external_consistency_module_summary.csv",
    "figure_qc_report": ROOT / "out/06_manuscript_figures/figure_qc_report.yaml",
    "banned_phrases_audit": ROOT / "out/06_manuscript_figures/banned_phrases_audit.csv",
    "panel_metadata": ROOT / "out/06_manuscript_figures/panel_metadata.csv",
}

EXPECTED_SPATIAL_LOG_COLUMNS = {
    "sample_id",
    "dataset",
    "status",
    "reason",
    "n_ref",
    "n_obs",
    "figure_status",
    "confidence_class",
    "confidence_reason",
}

ERROR_PATTERNS = re.compile(r"Traceback|\bERROR\b|Exception|Killed|CANCELLED|FAILED|DependencyNeverSatisfied", re.IGNORECASE)


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def run_cmd(cmd: list[str]) -> dict[str, Any]:
    try:
        proc = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=False, timeout=120)
        return {"ok": proc.returncode == 0, "returncode": proc.returncode, "stdout": proc.stdout.strip(), "stderr": proc.stderr.strip()}
    except Exception as exc:  # pragma: no cover - defensive audit path
        return {"ok": False, "returncode": None, "stdout": "", "stderr": str(exc)}


def file_record(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"exists": False, "path": rel(path)}
    stat = path.stat()
    return {
        "exists": True,
        "path": rel(path),
        "size_bytes": stat.st_size,
        "modified": dt.datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
    }


def csv_profile(path: Path, max_rows: int | None = None) -> dict[str, Any]:
    rec = file_record(path)
    if not rec.get("exists"):
        return rec
    try:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            columns = reader.fieldnames or []
            n_rows = 0
            counts: dict[str, dict[str, int]] = {}
            for row in reader:
                n_rows += 1
                if max_rows and n_rows > max_rows:
                    break
                for col in ("status", "confidence_class", "figure_status", "external_consistency_class", "strength_class"):
                    if col in row:
                        counts.setdefault(col, {})[row[col]] = counts.setdefault(col, {}).get(row[col], 0) + 1
        rec.update({"n_rows": n_rows, "columns": columns, "value_counts": counts})
    except Exception as exc:
        rec.update({"read_error": str(exc)})
    return rec


def read_yaml(path: Path) -> dict[str, Any]:
    rec = file_record(path)
    if not rec.get("exists"):
        return rec
    try:
        rec["content"] = yaml.safe_load(path.read_text())
        rec["parse_ok"] = True
    except Exception as exc:
        rec["parse_ok"] = False
        rec["parse_error"] = str(exc)
    return rec


def audit_jobs() -> dict[str, Any]:
    if not JOB_IDS:
        return {"job_ids": [], "sacct": {"ok": False, "stderr": "PIPELINE_JOB_IDS was not set"}}
    sacct = run_cmd([
        "sacct",
        "-j",
        ",".join(JOB_IDS),
        "--format=JobID,JobName%28,State,ExitCode,Elapsed,MaxRSS,ReqMem",
        "--parsable2",
        "--noheader",
    ])
    return {"job_ids": JOB_IDS, "sacct": sacct}


def audit_logs() -> dict[str, Any]:
    records = []
    for job_id in JOB_IDS:
        for path in sorted((ROOT / "slurm/logs").glob(f"*_{job_id}.*")):
            text = path.read_text(errors="replace")[-20000:]
            hits = [line.strip() for line in text.splitlines() if ERROR_PATTERNS.search(line)]
            records.append({
                "job_id": job_id,
                "path": rel(path),
                "exists": True,
                "n_error_like_lines_tail": len(hits),
                "error_like_lines_tail": hits[-10:],
            })
    return {"logs": records}


def audit_results_index() -> dict[str, Any]:
    path = ROOT / "results_index.yaml"
    rec = read_yaml(path)
    rec.pop("content", None)
    return rec


def audit_spatial_log() -> dict[str, Any]:
    path = KEY_OUTPUTS["spatial_cnv_run_log"]
    rec = csv_profile(path)
    columns = set(rec.get("columns", []))
    rec["missing_expected_columns"] = sorted(EXPECTED_SPATIAL_LOG_COLUMNS - columns)
    return rec


def audit_figure_qc() -> dict[str, Any]:
    qc_path = KEY_OUTPUTS["figure_qc_report"]
    qc = read_yaml(qc_path)
    content = qc.pop("content", None)
    qc["summary"] = content if isinstance(content, dict) else {}
    banned = csv_profile(KEY_OUTPUTS["banned_phrases_audit"])
    return {"figure_qc_report": qc, "banned_phrases_audit": banned}


def classify_status(audit: dict[str, Any]) -> str:
    missing = [name for name, rec in audit["key_outputs"].items() if not rec.get("exists")]
    index_ok = audit["results_index"].get("parse_ok") is True
    qc_summary = audit.get("figure_qc", {}).get("figure_qc_report", {}).get("summary", {})
    blocked = int(qc_summary.get("n_blocked_figures", 0) or 0) if isinstance(qc_summary, dict) else 0
    banned = int(qc_summary.get("n_banned_phrase_classes", 0) or 0) if isinstance(qc_summary, dict) else 0
    if missing or not index_ok:
        return "incomplete_or_failed"
    if blocked or banned:
        return "complete_needs_figure_review"
    return "complete_passed_basic_audit"


def write_markdown(audit: dict[str, Any], path: Path) -> None:
    lines = [
        "# HCC reframe pipeline post-run audit",
        "",
        f"Generated: {audit['generated']}",
        f"Overall status: **{audit['overall_status']}**",
        "",
        "## Job IDs",
        "",
        ", ".join(audit.get("jobs", {}).get("job_ids", [])) or "Not recorded",
        "",
        "## Key outputs",
        "",
        "| output | exists | path | rows/status |",
        "|---|---:|---|---|",
    ]
    for name, rec in audit["key_outputs"].items():
        details = ""
        if "n_rows" in rec:
            details = f"rows={rec['n_rows']}"
        lines.append(f"| {name} | {rec.get('exists')} | {rec.get('path', '')} | {details} |")
    lines += ["", "## Spatial inferCNV run log", ""]
    spatial = audit.get("spatial_cnv_run_log", {})
    lines.append(f"Missing expected columns: {spatial.get('missing_expected_columns', [])}")
    lines.append(f"Value counts: {spatial.get('value_counts', {})}")
    lines += ["", "## Figure QC", ""]
    fig_qc = audit.get("figure_qc", {}).get("figure_qc_report", {}).get("summary", {})
    if isinstance(fig_qc, dict):
        for key in ["status", "n_figure_plans", "n_blocked_figures", "n_low_qc_figures", "n_png_images", "n_image_warnings", "n_banned_phrase_classes"]:
            if key in fig_qc:
                lines.append(f"- {key}: {fig_qc[key]}")
    else:
        lines.append("No figure QC summary found.")
    lines += ["", "## Results index", ""]
    lines.append(f"Parse OK: {audit.get('results_index', {}).get('parse_ok')}")
    lines += ["", "## Log warnings", ""]
    for rec in audit.get("logs", {}).get("logs", []):
        if rec.get("n_error_like_lines_tail", 0):
            lines.append(f"- {rec['path']}: {rec['n_error_like_lines_tail']} error-like lines in tail")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    key_outputs: dict[str, dict[str, Any]] = {}
    for name, path in KEY_OUTPUTS.items():
        key_outputs[name] = csv_profile(path) if path.suffix.lower() == ".csv" else file_record(path)

    audit: dict[str, Any] = {
        "generated": NOW,
        "root": str(ROOT),
        "jobs": audit_jobs(),
        "logs": audit_logs(),
        "key_outputs": key_outputs,
        "spatial_cnv_run_log": audit_spatial_log(),
        "figure_qc": audit_figure_qc(),
        "results_index": audit_results_index(),
    }
    audit["overall_status"] = classify_status(audit)

    yaml_path = AUDIT_DIR / "hcc_reframe_postrun_audit.yaml"
    md_path = AUDIT_DIR / "hcc_reframe_postrun_audit.md"
    yaml_path.write_text(yaml.safe_dump(audit, sort_keys=False, allow_unicode=True))
    write_markdown(audit, md_path)

    print(f"Post-run audit written: {yaml_path}")
    print(f"Post-run audit written: {md_path}")
    print(f"Overall status: {audit['overall_status']}")


if __name__ == "__main__":
    main()
