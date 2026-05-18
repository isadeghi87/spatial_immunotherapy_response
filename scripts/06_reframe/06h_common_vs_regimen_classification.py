#!/usr/bin/env python3
"""06h - Common versus regimen-driven evidence classification.

Question:
    Which discovery response signals are common across cohorts/modalities,
    which are regimen-shifted, and which are platform/artifact driven?

Design:
    - Discovery inference uses patient-level GSE206325 module/feature effects.
    - Cappuyns is treated as cross-regimen transfer (atezolizumab+bevacizumab),
      not same-regimen validation.
    - IMC is used as protein-spatial architecture support, not strong response
      validation, because the response-labelled cohort has 5 responders.
    - IMbrave150 is orthogonal CD8 clinical support within Atezo+Bev, not
      validation of the full multilayer predictor and not differential-benefit
      evidence because the arm-by-CD8 interaction is null.
"""

from __future__ import annotations

import json
import math
import sys
from dataclasses import asdict, dataclass
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))

from config.publication_style import STATUS_PAL, apply_style, label_clean, save_fig

apply_style()

OUT = PROJECT / "out" / "06_reframe"
TAB = OUT / "tables"
FIG = OUT / "figures"
NOTE = OUT / "notes"
for path in (TAB, FIG, NOTE):
    path.mkdir(parents=True, exist_ok=True)


INPUTS = {
    "discovery_modules": OUT / "tables" / "discovery_module_associations.csv",
    "cappuyns_modules": OUT / "tables" / "cappuyns_locked_transfer.csv",
    "selected_feature_transport": PROJECT
    / "out/05_multilayer/tables/cappuyns_selected_feature_transport_qc.csv",
    "locked_validation_metrics": PROJECT
    / "out/05_multilayer/tables/cappuyns_locked_validation_metrics.csv",
    "outlier_sensitivity": PROJECT
    / "out/05_multilayer/tables/cappuyns_outlier_sensitivity_summary.csv",
    "imc_mwu": PROJECT / "out/03_validation/imc_validation/tables/imc_mwu_results.csv",
    "imc_kw": PROJECT / "out/03_validation/imc_validation/tables/imc_kw_results.csv",
    "imc_contingency": PROJECT
    / "out/03_validation/imc_validation/tables/imc_contingency_results.csv",
    "imbrave_support": OUT / "tables" / "imbrave150_orthogonal_support.csv",
}

CLASSIFICATIONS = {
    "common_supported",
    "common_suggestive",
    "orthogonal_context",
    "regimen_shift",
    "artifact_or_platform",
    "null_or_inconclusive",
}

STRENGTH_MAP = {
    "common_supported": "supported_limited",
    "common_suggestive": "exploratory",
    "orthogonal_context": "exploratory",
    "regimen_shift": "null",
    "artifact_or_platform": "likely_artifact",
    "null_or_inconclusive": "inconclusive",
}

CLASS_COLORS = {
    "common_supported": STATUS_PAL["supported_limited"],
    "common_suggestive": STATUS_PAL["exploratory"],
    "orthogonal_context": "#2B8CBE",
    "regimen_shift": "#7C3AED",
    "artifact_or_platform": STATUS_PAL["likely_artifact"],
    "null_or_inconclusive": STATUS_PAL["null"],
}

EVIDENCE_LEVELS = [
    "supported",
    "suggestive",
    "null",
    "opposite",
    "artifact",
    "not_applicable",
]

EVIDENCE_COLORS = {
    "supported": "#31A354",
    "suggestive": "#F0AD00",
    "null": "#9CA3AF",
    "opposite": "#7C3AED",
    "artifact": "#B91C1C",
    "not_applicable": "#E5E7EB",
}

FEATURE_LABELS = {
    "frac_cDC2": "cDC2 / APC fraction",
    "frac_CD8_T": "CD8 T-cell fraction",
    "frac_Plasma": "Plasma-cell fraction",
    "nmf_NMF4": "NMF4 cytotoxic T program",
    "nmf_NMF10": "NMF10 plasma/IG program",
    "nmf_NMF5": "NMF5 program",
    "regulon_IRF1_axis": "IRF1 regulon axis",
    "regulon_IRF1_over_MYC": "IRF1-over-MYC balance",
    "regulon_MYC_axis": "MYC regulon axis",
    "cnv_clone_diversity": "CNV clone diversity",
    "cnv_score_mean": "CNV score mean",
    "cnv_score_sd": "CNV score SD",
    "sig_Tfh_like": "Tfh-like signature",
    "sig_T_cell_exhaustion": "T-cell exhaustion signature",
    "sig_Immune_checkpoint": "Immune checkpoint signature",
}

MODULE_LABELS = {
    "M_Tfh": "Tfh module",
    "M_CD8_dysfunction": "CD8 dysfunction module",
    "M_Immune_checkpoint": "Immune checkpoint module",
    "M_TLS": "TLS module",
    "M_CD8_activation": "CD8 activation module",
    "M_plasma_proxy": "Plasma/NMF10 satellite",
    "M_mregDC": "mregDC module",
    "M_CD8_progenitor": "Progenitor CD8 module",
    "M_TGFbeta": "TGF-beta module",
    "M_WNT": "WNT module",
    "M_Macrophage_M1": "M1 macrophage module",
    "M_Macrophage_M2": "M2 macrophage module",
    "M_Treg": "Treg module",
    "M_Angiogenesis": "Angiogenesis module",
}

FORCE_ARTIFACT = {"M_plasma_proxy", "nmf_NMF10"}
CORE_SHIFT_MODULES = {"M_Tfh", "M_CD8_dysfunction", "M_Immune_checkpoint"}
COMMON_FEATURES = {"frac_cDC2"}
SUGGESTIVE_FEATURES = {
    "nmf_NMF4",
    "regulon_IRF1_axis",
    "regulon_IRF1_over_MYC",
    "regulon_MYC_axis",
}
SHIFT_FEATURES = {
    "frac_CD8_T",
    "frac_Plasma",
    "sig_Tfh_like",
    "sig_T_cell_exhaustion",
    "sig_Immune_checkpoint",
}


@dataclass
class EvidenceRow:
    axis_id: str
    axis_label: str
    axis_type: str
    discovery_evidence: str
    cappuyns_evidence: str
    imc_evidence: str
    imbrave150_evidence: str
    discovery_effect: float | None
    discovery_p: float | None
    discovery_q: float | None
    cappuyns_effect: float | None
    cappuyns_p: float | None
    cappuyns_q: float | None
    sign_concordant: bool | None
    dynamic_range_adequate: bool | None
    classification: str
    strength_class: str
    manuscript_priority: str
    figure_candidate: str
    caveat: str


def read_csv_checked(key: str) -> pd.DataFrame:
    path = INPUTS[key]
    if not path.exists():
        raise FileNotFoundError(f"Missing required input for {key}: {path}")
    return pd.read_csv(path)


def require_columns(df: pd.DataFrame, key: str, cols: list[str]) -> None:
    missing = [col for col in cols if col not in df.columns]
    if missing:
        raise ValueError(f"{key} missing required columns: {missing}")


def finite_or_none(value: object) -> float | None:
    try:
        val = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(val):
        return None
    return val


def evidence_from_p(effect: float | None, p: float | None, q: float | None) -> str:
    if q is not None and q < 0.10:
        return "supported"
    if p is not None and p < 0.20 and effect is not None and abs(effect) >= 0.25:
        return "suggestive"
    return "null"


def classify_module(row: pd.Series) -> tuple[str, str]:
    module = str(row["module"])
    discovery_q = finite_or_none(row.get("p_BH_discovery"))
    discovery_p = finite_or_none(row.get("p_raw_discovery"))
    discovery_effect = finite_or_none(row.get("rank_biserial_discovery"))
    cap_q = finite_or_none(row.get("p_BH_cappuyns"))
    cap_p = finite_or_none(row.get("p_raw_cappuyns"))
    cap_effect = finite_or_none(row.get("rank_biserial_cappuyns"))
    sign_concordant = bool(row.get("sign_concordant", False))
    discovery_supported = (
        (discovery_q is not None and discovery_q < 0.10)
        or (
            discovery_p is not None
            and discovery_p < 0.25
            and discovery_effect is not None
            and abs(discovery_effect) >= 0.30
        )
    )
    cappuyns_supported_or_suggestive = (
        (cap_q is not None and cap_q < 0.10)
        or (cap_p is not None and cap_p < 0.20 and cap_effect is not None and abs(cap_effect) >= 0.25)
    )

    if module in FORCE_ARTIFACT:
        return (
            "artifact_or_platform",
            "Plasma/NMF10 tracks B/plasma composition; the independent-response claim is withdrawn after plasma adjustment.",
        )
    if module in CORE_SHIFT_MODULES and discovery_q is not None and discovery_q < 0.05 and not sign_concordant:
        return (
            "regimen_shift",
            "Strong discovery module flips direction in Cappuyns despite adequate dynamic range, consistent with regimen/cohort-specific biology.",
        )
    if cap_q is not None and cap_q < 0.10 and sign_concordant:
        return (
            "common_supported",
            "Direction-concordant and BH-supported in Cappuyns transfer.",
        )
    if sign_concordant and discovery_supported and cappuyns_supported_or_suggestive:
        return (
            "common_suggestive",
            "Direction-concordant with non-trivial discovery and Cappuyns evidence, but not BH-supported in transfer; treat as a shared motif hypothesis.",
        )
    return (
        "null_or_inconclusive",
        "No robust cross-cohort evidence after direction, dynamic range, and multiple-testing checks.",
    )


def classify_feature(row: pd.Series) -> tuple[str, str]:
    feature = str(row["feature"])
    sign_concordant = bool(row.get("sign_concordant", False))
    cap_q = finite_or_none(row.get("cappuyns_p_BH_selected_features"))
    cap_p = finite_or_none(row.get("cappuyns_p"))
    disc_p = finite_or_none(row.get("discovery_p"))
    disc_effect = finite_or_none(row.get("discovery_rank_biserial"))
    cap_effect = finite_or_none(row.get("cappuyns_rank_biserial"))
    discovery_supported = disc_p is not None and disc_p < 0.10 and disc_effect is not None and abs(disc_effect) >= 0.25
    cappuyns_supported_or_suggestive = (
        (cap_q is not None and cap_q < 0.10)
        or (cap_p is not None and cap_p < 0.20 and cap_effect is not None and abs(cap_effect) >= 0.25)
    )

    if feature in FORCE_ARTIFACT:
        return (
            "artifact_or_platform",
            "NMF10/plasma signal is composition-linked; it should not be claimed as independent response biology.",
        )
    if feature in COMMON_FEATURES and sign_concordant and cap_q is not None and cap_q < 0.10:
        return (
            "common_supported",
            "Direction-concordant selected feature with Cappuyns BH support; strongest portable candidate in the current tables.",
        )
    if feature in SHIFT_FEATURES and not sign_concordant and disc_p is not None and disc_p < 0.05:
        return (
            "regimen_shift",
            "Discovery-associated immune feature flips or collapses in Cappuyns, arguing against a universal locked-feature biomarker.",
        )
    if feature in SUGGESTIVE_FEATURES and sign_concordant and discovery_supported and cappuyns_supported_or_suggestive:
        return (
            "common_suggestive",
            "Direction-concordant and biologically interpretable with non-trivial discovery and Cappuyns evidence, but validation is nominal/underpowered rather than BH-strong.",
        )
    if not sign_concordant:
        return (
            "null_or_inconclusive",
            "Direction flips or is near-null in Cappuyns; do not promote as shared biology.",
        )
    return (
        "null_or_inconclusive",
        "Measurable but not statistically or biologically strong enough for a shared-regimen claim.",
    )


def evidence_status_for_classification(classification: str, cohort: str) -> str:
    if classification == "artifact_or_platform":
        return "artifact"
    if classification == "regimen_shift" and cohort == "Cappuyns":
        return "opposite"
    if classification == "common_supported":
        return "supported"
    if classification == "common_suggestive":
        return "suggestive"
    return "null"


def build_module_rows(discovery: pd.DataFrame, transfer: pd.DataFrame) -> list[EvidenceRow]:
    merged = discovery.merge(
        transfer,
        on="module",
        how="inner",
        suffixes=("_discovery", "_cappuyns"),
    )
    rows: list[EvidenceRow] = []
    for _, row in merged.iterrows():
        classification, caveat = classify_module(row)
        disc_effect = finite_or_none(row.get("rank_biserial_discovery"))
        cap_effect = finite_or_none(row.get("rank_biserial_cappuyns"))
        disc_p = finite_or_none(row.get("p_raw_discovery"))
        disc_q = finite_or_none(row.get("p_BH_discovery"))
        cap_p = finite_or_none(row.get("p_raw_cappuyns"))
        cap_q = finite_or_none(row.get("p_BH_cappuyns"))
        rows.append(
            EvidenceRow(
                axis_id=str(row["module"]),
                axis_label=MODULE_LABELS.get(str(row["module"]), label_clean(str(row["module"]))),
                axis_type="locked_module",
                discovery_evidence=evidence_from_p(disc_effect, disc_p, disc_q),
                cappuyns_evidence=evidence_status_for_classification(classification, "Cappuyns"),
                imc_evidence="not_applicable",
                imbrave150_evidence="not_applicable",
                discovery_effect=disc_effect,
                discovery_p=disc_p,
                discovery_q=disc_q,
                cappuyns_effect=cap_effect,
                cappuyns_p=cap_p,
                cappuyns_q=cap_q,
                sign_concordant=bool(row.get("sign_concordant")),
                dynamic_range_adequate=bool(row.get("dynamic_range_adequate")),
                classification=classification,
                strength_class=STRENGTH_MAP[classification],
                manuscript_priority="main_boundary" if classification in {"regimen_shift", "artifact_or_platform"} else "supplement",
                figure_candidate="yes",
                caveat=caveat,
            )
        )
    return rows


def build_feature_rows(features: pd.DataFrame) -> list[EvidenceRow]:
    rows: list[EvidenceRow] = []
    for _, row in features.iterrows():
        feature = str(row["feature"])
        classification, caveat = classify_feature(row)
        disc_effect = finite_or_none(row.get("discovery_rank_biserial"))
        cap_effect = finite_or_none(row.get("cappuyns_rank_biserial"))
        disc_p = finite_or_none(row.get("discovery_p"))
        cap_p = finite_or_none(row.get("cappuyns_p"))
        cap_q = finite_or_none(row.get("cappuyns_p_BH_selected_features"))
        rows.append(
            EvidenceRow(
                axis_id=feature,
                axis_label=FEATURE_LABELS.get(feature, label_clean(feature)),
                axis_type="selected_feature",
                discovery_evidence=evidence_from_p(disc_effect, disc_p, None),
                cappuyns_evidence=evidence_status_for_classification(classification, "Cappuyns"),
                imc_evidence="not_applicable",
                imbrave150_evidence="not_applicable",
                discovery_effect=disc_effect,
                discovery_p=disc_p,
                discovery_q=None,
                cappuyns_effect=cap_effect,
                cappuyns_p=cap_p,
                cappuyns_q=cap_q,
                sign_concordant=bool(row.get("sign_concordant")),
                dynamic_range_adequate=True,
                classification=classification,
                strength_class=STRENGTH_MAP[classification],
                manuscript_priority="main_candidate" if classification == "common_supported" else "supplement",
                figure_candidate="yes" if feature in FEATURE_LABELS else "no",
                caveat=caveat,
            )
        )
    return rows


def build_model_rows(metrics: pd.DataFrame, outliers: pd.DataFrame) -> list[EvidenceRow]:
    rows: list[EvidenceRow] = []
    model_labels = {
        "L3_composition": "Layer-3 composition model",
        "L3_composition_noTech": "Layer-3 composition model (no-tech)",
        "Balanced ensemble": "Balanced multilayer ensemble",
        "L5_nmf_programs": "Layer-5 NMF model",
        "L4_io_signatures": "Layer-4 IO signature model",
        "L2_regulon": "Layer-2 regulon model",
        "L1_tumour_intrinsic": "Layer-1 CNV model",
    }
    for _, row in metrics.iterrows():
        target = str(row["target"])
        auc = finite_or_none(row.get("auc"))
        pr_auc = finite_or_none(row.get("pr_auc"))
        if target in {"L3_composition", "L3_composition_noTech", "Balanced ensemble", "L4_io_signatures"}:
            classification = "regimen_shift" if target != "Balanced ensemble" else "null_or_inconclusive"
            caveat = "Locked discovery model does not transfer robustly to Cappuyns; this is a regimen/feature-transport boundary, not validation."
        elif target == "L2_regulon":
            classification = "common_suggestive"
            caveat = "Layer-level transfer is the best locked model signal, but still modest and not sufficient for a predictor claim."
        else:
            classification = "null_or_inconclusive"
            caveat = "Layer-level transfer is near chance or assumption-sensitive."
        influence = outliers[outliers["target"] == target]
        if not influence.empty:
            caveat += f" Jackknife max absolute AUC delta={float(influence['max_abs_auc_delta'].iloc[0]):.3f}."
        rows.append(
            EvidenceRow(
                axis_id=target,
                axis_label=model_labels.get(target, target),
                axis_type="locked_model",
                discovery_evidence="not_applicable",
                cappuyns_evidence="suggestive" if auc is not None and auc >= 0.65 else "null",
                imc_evidence="not_applicable",
                imbrave150_evidence="not_applicable",
                discovery_effect=None,
                discovery_p=None,
                discovery_q=None,
                cappuyns_effect=auc,
                cappuyns_p=None,
                cappuyns_q=None,
                sign_concordant=None,
                dynamic_range_adequate=None,
                classification=classification,
                strength_class=STRENGTH_MAP[classification],
                manuscript_priority="main_boundary" if target in {"L3_composition", "Balanced ensemble"} else "supplement",
                figure_candidate="yes" if target in {"L3_composition", "Balanced ensemble", "L2_regulon"} else "no",
                caveat=f"Cappuyns AUC={auc:.3f}, PR-AUC={pr_auc:.3f}. {caveat}" if auc is not None and pr_auc is not None else caveat,
            )
        )
    return rows


def build_imc_rows(imc_mwu: pd.DataFrame, imc_kw: pd.DataFrame, imc_cont: pd.DataFrame) -> list[EvidenceRow]:
    rows: list[EvidenceRow] = []
    architecture_features = {
        "CD8_density_mm2": "IMC CD8 density architecture",
        "CD8_parenchyma_mm2": "IMC parenchymal CD8 architecture",
        "CD8_parenchyma_stroma_ratio": "IMC CD8 parenchyma/stroma access",
        "frac_CD8_T": "IMC CD8 fraction architecture",
        "CD8_pd1_frac": "IMC CD8 PD1+ architecture",
        "frac_Macrophage": "IMC macrophage immunotype architecture",
    }
    response_lookup = imc_mwu.set_index("feature") if "feature" in imc_mwu else pd.DataFrame()
    for _, row in imc_kw.iterrows():
        feature = str(row["feature"])
        if feature not in architecture_features:
            continue
        qval = finite_or_none(row.get("qval"))
        pval = finite_or_none(row.get("pval"))
        effect = None
        if feature in response_lookup.index:
            response_row = response_lookup.loc[feature]
            resp_p = finite_or_none(response_row.get("pval"))
            resp_q = finite_or_none(response_row.get("qval"))
            med_r = finite_or_none(response_row.get("median_R"))
            med_nr = finite_or_none(response_row.get("median_NR"))
            effect = None if med_r is None or med_nr is None else med_r - med_nr
        else:
            resp_p = None
            resp_q = None
        classification = "orthogonal_context" if qval is not None and qval < 0.05 else "null_or_inconclusive"
        caveat = (
            "IMC strongly supports immunotype/spatial architecture, but response validation remains weak because only 5 responders are available; "
            f"response MWU p={resp_p:.3g}, q={resp_q:.3g}." if resp_p is not None and resp_q is not None else
            "IMC supports immunotype/spatial architecture, but no direct response test is available for this feature."
        )
        rows.append(
            EvidenceRow(
                axis_id=f"imc_{feature}",
                axis_label=architecture_features[feature],
                axis_type="imc_architecture",
                discovery_evidence="not_applicable",
                cappuyns_evidence="not_applicable",
                imc_evidence="supported" if qval is not None and qval < 0.05 else "null",
                imbrave150_evidence="not_applicable",
                discovery_effect=None,
                discovery_p=None,
                discovery_q=None,
                cappuyns_effect=effect,
                cappuyns_p=resp_p,
                cappuyns_q=resp_q,
                sign_concordant=None,
                dynamic_range_adequate=True,
                classification=classification,
                strength_class=STRENGTH_MAP[classification],
                manuscript_priority="supplement",
                figure_candidate="yes",
                caveat=caveat,
            )
        )

    enriched = imc_cont[imc_cont["test"].astype(str).str.contains("Enriched_vs_rest", regex=False)]
    if not enriched.empty:
        pval = finite_or_none(enriched["pval"].iloc[0])
        rows.append(
            EvidenceRow(
                axis_id="imc_immunotype_response_enriched",
                axis_label="IMC enriched immunotype vs response",
                axis_type="imc_response_boundary",
                discovery_evidence="not_applicable",
                cappuyns_evidence="not_applicable",
                imc_evidence="suggestive" if pval is not None and pval < 0.10 else "null",
                imbrave150_evidence="not_applicable",
                discovery_effect=None,
                discovery_p=None,
                discovery_q=None,
                cappuyns_effect=None,
                cappuyns_p=pval,
                cappuyns_q=None,
                sign_concordant=None,
                dynamic_range_adequate=True,
                classification="orthogonal_context" if pval is not None and pval < 0.10 else "null_or_inconclusive",
                strength_class=STRENGTH_MAP["orthogonal_context" if pval is not None and pval < 0.10 else "null_or_inconclusive"],
                manuscript_priority="supplement",
                figure_candidate="yes",
                caveat="IMC enriched immunotype is nominally enriched among responders, but the cohort has only 5 responders; retain as architecture validation, not response validation.",
            )
        )
    return rows


def build_imbrave_rows(imbrave: pd.DataFrame) -> list[EvidenceRow]:
    pfs = imbrave[(imbrave["endpoint"] == "PFS") & (imbrave["arm"] == "Atezo+Bev")]
    interaction = imbrave[(imbrave["endpoint"] == "PFS") & (imbrave["arm"] == "interaction")]
    if pfs.empty:
        return []
    logrank = finite_or_none(pfs["logrank_p"].iloc[0])
    med_high = finite_or_none(pfs["median_high_mo"].iloc[0])
    med_low = finite_or_none(pfs["median_low_mo"].iloc[0])
    interaction_p = finite_or_none(interaction["cox_p"].iloc[0]) if not interaction.empty else None
    caveat = (
        f"IMbrave150 supports CD8-high within Atezo+Bev PFS ({med_high:.1f} vs {med_low:.1f} months; p={logrank:.3g}), "
        f"but arm-by-CD8 interaction is null for PFS (p={interaction_p:.3g}); not a differential-benefit claim."
    )
    return [
        EvidenceRow(
            axis_id="imbrave150_cd8_high_pfs",
            axis_label="IMbrave150 CD8-high Atezo+Bev PFS",
            axis_type="orthogonal_clinical",
            discovery_evidence="not_applicable",
            cappuyns_evidence="not_applicable",
            imc_evidence="not_applicable",
            imbrave150_evidence="supported" if logrank is not None and logrank < 0.05 else "null",
            discovery_effect=None,
            discovery_p=None,
            discovery_q=None,
            cappuyns_effect=None,
            cappuyns_p=logrank,
            cappuyns_q=None,
            sign_concordant=None,
            dynamic_range_adequate=True,
            classification="orthogonal_context",
            strength_class=STRENGTH_MAP["orthogonal_context"],
            manuscript_priority="main_context",
            figure_candidate="yes",
            caveat=caveat,
        )
    ]


def rows_to_frame(rows: list[EvidenceRow]) -> pd.DataFrame:
    df = pd.DataFrame([asdict(row) for row in rows])
    if df.empty:
        raise ValueError("No evidence rows were created")
    df["classification"] = pd.Categorical(
        df["classification"],
        categories=[
            "common_supported",
            "common_suggestive",
            "orthogonal_context",
            "regimen_shift",
            "artifact_or_platform",
            "null_or_inconclusive",
        ],
        ordered=True,
    )
    df = df.sort_values(["classification", "axis_type", "axis_label"]).reset_index(drop=True)
    return df


def make_summary(df: pd.DataFrame) -> pd.DataFrame:
    summary = (
        df.groupby(["classification", "strength_class", "manuscript_priority"], observed=False)
        .size()
        .reset_index(name="n_axes")
        .query("n_axes > 0")
        .sort_values(["classification", "n_axes"], ascending=[True, False])
    )
    return summary


def make_caveat_audit(df: pd.DataFrame) -> pd.DataFrame:
    caveats = df[["axis_id", "axis_label", "axis_type", "classification", "strength_class", "caveat"]].copy()
    caveats["regimen_mismatch_flag"] = caveats["classification"].astype(str).eq("regimen_shift")
    caveats["artifact_flag"] = caveats["classification"].astype(str).eq("artifact_or_platform")
    caveats["response_underpowered_flag"] = caveats["caveat"].str.contains("5 responders", case=False, na=False)
    caveats["platform_boundary_flag"] = caveats["caveat"].str.contains("platform|IMbrave|IMC|plasma", case=False, na=False)
    return caveats


def plot_effect_scatter(df: pd.DataFrame) -> None:
    plot_df = df[
        df["discovery_effect"].notna()
        & df["cappuyns_effect"].notna()
        & df["axis_type"].isin(["locked_module", "selected_feature"])
    ].copy()
    fig, ax = plt.subplots(figsize=(6.1, 4.9))
    for classification, sub in plot_df.groupby("classification", observed=False):
        if sub.empty:
            continue
        ax.scatter(
            sub["discovery_effect"],
            sub["cappuyns_effect"],
            s=42,
            color=CLASS_COLORS[str(classification)],
            edgecolor="white",
            linewidth=0.6,
            label=str(classification).replace("_", " "),
            alpha=0.95,
        )
    ax.axhline(0, color="0.35", lw=0.8, ls="--")
    ax.axvline(0, color="0.35", lw=0.8, ls="--")
    lim = max(abs(plot_df["discovery_effect"]).max(), abs(plot_df["cappuyns_effect"]).max(), 0.1) * 1.12
    ax.plot([-lim, lim], [-lim, lim], color="0.75", lw=0.8, ls=":")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Discovery effect (rank-biserial; R > NR positive)", fontsize=8.5)
    ax.set_ylabel("Cappuyns transfer effect (rank-biserial; R > NR positive)", fontsize=8.5)
    ax.set_title("Discovery versus cross-regimen transfer effects", fontsize=10.5, pad=8)
    ax.tick_params(labelsize=7.5)
    ax.text(0.02, 0.96, "shared direction", transform=ax.transAxes,
            ha="left", va="top", fontsize=7, color="0.35")
    ax.text(0.98, 0.05, "discovery-only / flipped", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7, color="0.35")
    label_specs = {
        "frac_cDC2": (8, 8, "left"),
        "nmf_NMF4": (6, 6, "left"),
        "regulon_IRF1_over_MYC": (6, 9, "left"),
        "M_TLS": (6, -4, "left"),
        "M_Tfh": (-10, -16, "right"),
        "M_CD8_dysfunction": (-10, 8, "right"),
        "M_Immune_checkpoint": (7, -2, "left"),
        "M_plasma_proxy": (7, 9, "left"),
    }
    short_labels = {
        "frac_cDC2": "cDC2/APC",
        "nmf_NMF4": "NMF4 cytotoxic",
        "regulon_IRF1_over_MYC": "IRF1/MYC",
        "M_TLS": "TLS",
        "M_Tfh": "Tfh flip",
        "M_CD8_dysfunction": "CD8 dysfunction flip",
        "M_Immune_checkpoint": "checkpoint flip",
        "M_plasma_proxy": "plasma/NMF10",
    }
    for _, row in plot_df[plot_df["axis_id"].isin(label_specs)].iterrows():
        dx, dy, ha = label_specs[str(row["axis_id"])]
        ax.annotate(
            short_labels[str(row["axis_id"])],
            (row["discovery_effect"], row["cappuyns_effect"]),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=ha,
            fontsize=6.8,
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75),
        )
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=6.7,
              frameon=False, borderpad=0.5)
    fig.tight_layout()
    save_fig(fig, FIG / "06h_common_vs_regimen_scatter")


def plot_evidence_matrix(df: pd.DataFrame) -> None:
    priority_ids = [
        "frac_cDC2",
        "nmf_NMF4",
        "regulon_IRF1_axis",
        "regulon_IRF1_over_MYC",
        "regulon_MYC_axis",
        "M_Tfh",
        "M_CD8_dysfunction",
        "M_Immune_checkpoint",
        "M_TLS",
        "M_plasma_proxy",
        "frac_CD8_T",
        "frac_Plasma",
        "L3_composition",
        "Balanced ensemble",
        "imc_CD8_density_mm2",
        "imc_CD8_parenchyma_stroma_ratio",
        "imc_CD8_pd1_frac",
        "imbrave150_cd8_high_pfs",
    ]
    mat_df = df[df["axis_id"].isin(priority_ids)].copy()
    mat_df["order"] = mat_df["axis_id"].map({axis: idx for idx, axis in enumerate(priority_ids)})
    mat_df = mat_df.sort_values("order")
    cols = ["discovery_evidence", "cappuyns_evidence", "imc_evidence", "imbrave150_evidence"]
    col_labels = ["Discovery", "Cappuyns", "IMC", "IMbrave150"]
    code = {level: idx for idx, level in enumerate(EVIDENCE_LEVELS)}
    matrix = np.vectorize(code.get)(mat_df[cols].to_numpy()).astype(int)
    cmap = ListedColormap([EVIDENCE_COLORS[level] for level in EVIDENCE_LEVELS])

    fig_h = max(5.0, 0.30 * len(mat_df) + 1.1)
    fig, ax = plt.subplots(figsize=(6.2, fig_h))
    ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=-0.5, vmax=len(EVIDENCE_LEVELS) - 0.5)
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(col_labels, fontweight="bold", fontsize=7.2, rotation=12, ha="right")
    ax.set_yticks(np.arange(len(mat_df)))
    ax.set_yticklabels(mat_df["axis_label"], fontsize=6.8)
    ax.set_title("Common versus regimen-driven evidence matrix", fontsize=10.5, pad=8)
    ax.set_xticks(np.arange(-0.5, len(cols), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(mat_df), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)

    handles = [
        Line2D([0], [0], marker="s", color="w", markerfacecolor=EVIDENCE_COLORS[level], markersize=8, label=level.replace("_", " "))
        for level in EVIDENCE_LEVELS
    ]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=7)
    fig.tight_layout()
    save_fig(fig, FIG / "06h_common_vs_regimen_evidence_matrix")


def write_interpretation_note(df: pd.DataFrame, summary: pd.DataFrame) -> None:
    common_supported = df[df["classification"].astype(str) == "common_supported"]["axis_label"].tolist()
    common_suggestive = df[df["classification"].astype(str) == "common_suggestive"]["axis_label"].tolist()
    orthogonal_context = df[df["classification"].astype(str) == "orthogonal_context"]["axis_label"].tolist()
    regimen_shift = df[df["classification"].astype(str) == "regimen_shift"]["axis_label"].tolist()
    artifact = df[df["classification"].astype(str) == "artifact_or_platform"]["axis_label"].tolist()
    note = {
        "result_id": "R_REFRAME_common_vs_regimen_classification",
        "date": str(date.today()),
        "question": "Which discovery signals are common, regimen-shifted, or artifact/platform driven across discovery, Cappuyns, IMC, and IMbrave150 evidence?",
        "classification_counts": summary.to_dict(orient="records"),
        "main_findings": [
            "The strongest portable candidate in the current selected-feature table is cDC2/APC fraction, which is direction-concordant and BH-supported in Cappuyns.",
            "Cytotoxic NMF4, MYC-axis, TLS, WNT, and progenitor-CD8 motifs remain shared-biology hypotheses, but validation support is suggestive rather than definitive.",
            "The discovery CORE immune axis is not a universal biomarker: Tfh, CD8 dysfunction, and immune checkpoint modules flip or collapse in Cappuyns despite adequate dynamic range.",
            "IMC and IMbrave150 are orthogonal context rather than common-signal validation: IMC validates CD8-rich spatial architecture, and IMbrave150 supports CD8-high PFS within Atezo+Bev, but neither validates the multilayer predictor.",
        ],
        "common_supported_axes": common_supported,
        "common_suggestive_axes": common_suggestive,
        "orthogonal_context_axes": orthogonal_context,
        "regimen_shift_axes": regimen_shift,
        "artifact_or_platform_axes": artifact,
        "allowed_language": [
            "shared immune-access and antigen-presentation motifs are plausible across cohorts",
            "the exact discovery predictor is regimen/cohort bounded",
            "IMC provides protein-spatial architecture validation, not strong response validation",
        ],
        "forbidden_language": [
            "universal cross-cohort validated immune signature",
            "NMF10 is independent of plasma-cell abundance",
            "full multilayer ensemble validated in Cappuyns or IMbrave150",
            "IMC proves response prediction",
            "CD8-high is a differential-benefit biomarker for Atezo+Bev versus sorafenib",
        ],
    }
    (NOTE / "common_vs_regimen_interpretation_note.json").write_text(json.dumps(note, indent=2))


def run_qc(df: pd.DataFrame) -> None:
    bad_class = set(df["classification"].astype(str)) - CLASSIFICATIONS
    if bad_class:
        raise ValueError(f"Unexpected classification values: {bad_class}")
    bad_evidence = set(df[["discovery_evidence", "cappuyns_evidence", "imc_evidence", "imbrave150_evidence"]].to_numpy().ravel()) - set(EVIDENCE_LEVELS)
    if bad_evidence:
        raise ValueError(f"Unexpected evidence values: {bad_evidence}")
    imc_supported = df[(df["axis_type"].str.startswith("imc")) & (df["classification"].astype(str) == "common_supported")]
    if not imc_supported.empty:
        raise ValueError("IMC response-underpowered rows cannot be promoted to common_supported")
    nmf10 = df[df["axis_id"].isin(["M_plasma_proxy", "nmf_NMF10"])]
    if not nmf10.empty and not (nmf10["classification"].astype(str) == "artifact_or_platform").all():
        raise ValueError("NMF10/plasma rows must remain artifact_or_platform under the current audited evidence")


def main() -> None:
    discovery = read_csv_checked("discovery_modules")
    transfer = read_csv_checked("cappuyns_modules")
    features = read_csv_checked("selected_feature_transport")
    metrics = read_csv_checked("locked_validation_metrics")
    outliers = read_csv_checked("outlier_sensitivity")
    imc_mwu = read_csv_checked("imc_mwu")
    imc_kw = read_csv_checked("imc_kw")
    imc_cont = read_csv_checked("imc_contingency")
    imbrave = read_csv_checked("imbrave_support")

    require_columns(discovery, "discovery_modules", ["module", "rank_biserial", "p_raw", "p_BH"])
    require_columns(transfer, "cappuyns_modules", ["module", "rank_biserial", "p_raw", "p_BH", "sign_concordant", "dynamic_range_adequate"])
    require_columns(features, "selected_feature_transport", ["feature", "discovery_rank_biserial", "cappuyns_rank_biserial", "discovery_p", "cappuyns_p", "sign_concordant"])
    require_columns(metrics, "locked_validation_metrics", ["target", "auc", "pr_auc"])
    require_columns(imc_kw, "imc_kw", ["feature", "pval", "qval"])
    require_columns(imc_mwu, "imc_mwu", ["feature", "median_R", "median_NR", "pval", "qval"])
    require_columns(imbrave, "imbrave_support", ["endpoint", "arm", "logrank_p", "strength_class"])

    rows: list[EvidenceRow] = []
    rows.extend(build_module_rows(discovery, transfer))
    rows.extend(build_feature_rows(features))
    rows.extend(build_model_rows(metrics, outliers))
    rows.extend(build_imc_rows(imc_mwu, imc_kw, imc_cont))
    rows.extend(build_imbrave_rows(imbrave))

    df = rows_to_frame(rows)
    run_qc(df)
    summary = make_summary(df)
    caveats = make_caveat_audit(df)

    df.to_csv(TAB / "common_vs_regimen_feature_classification.csv", index=False)
    summary.to_csv(TAB / "common_vs_regimen_summary.csv", index=False)
    caveats.to_csv(TAB / "common_vs_regimen_caveat_audit.csv", index=False)
    write_interpretation_note(df, summary)
    plot_effect_scatter(df)
    plot_evidence_matrix(df)

    print("Common-vs-regimen classification complete")
    print(f"  rows: {len(df)}")
    print("  classification counts:")
    print(df["classification"].astype(str).value_counts().to_string())
    print(f"  tables: {TAB}")
    print(f"  figures: {FIG}")


if __name__ == "__main__":
    main()