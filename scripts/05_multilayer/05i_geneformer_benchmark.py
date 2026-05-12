#!/usr/bin/env python3
"""
05i_geneformer_benchmark.py - Benchmark Geneformer as an L6 representation layer.

Question
--------
Does a frozen Geneformer patient embedding improve discovery discrimination or,
more importantly, external transfer to the Cappuyns atezo+bev cohort?

Outputs
-------
  out/05_multilayer/tables/geneformer_selected_features_<model>.csv
  out/05_multilayer/tables/geneformer_layer_metrics_<model>.csv
  out/05_multilayer/tables/geneformer_validation_predictions_<model>.csv
  out/05_multilayer/figures/05i_geneformer_layer_benchmark_<model>.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

from config.publication_style import apply_style
from geneformer_utils import DEFAULT_MODEL_VARIANT, sanitize_model_name

apply_style()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark Geneformer as an L6 layer")
    parser.add_argument(
        "--feature-file",
        type=Path,
        default=None,
        help="Patient-level Geneformer feature file. Defaults to the standard output of 05h.",
    )
    parser.add_argument(
        "--model-variant",
        type=str,
        default=DEFAULT_MODEL_VARIANT,
        help="Model variant tag used to derive the default feature filename.",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=3,
        help="Number of top Geneformer dimensions to keep for the L6 layer.",
    )
    parser.add_argument(
        "--output-tag",
        type=str,
        default=None,
        help="Optional tag appended to benchmark outputs. Defaults to the model tag, or derives from the feature filename.",
    )
    return parser.parse_args()


def valid_feature(data: pd.DataFrame, column: str) -> bool:
    values = data[column]
    return values.notna().sum() >= len(data) * 0.5 and values.nunique(dropna=True) > 1


def pick_top_features(data: pd.DataFrame, response: np.ndarray, prefix: str, k: int) -> pd.DataFrame:
    rows = []
    for column in data.columns:
        if not column.startswith(prefix) or not valid_feature(data, column):
            continue
        responder = data.loc[response == 1, column].dropna()
        non_responder = data.loc[response == 0, column].dropna()
        if len(responder) < 2 or len(non_responder) < 2:
            continue
        t_stat, p_value = ttest_ind(responder, non_responder, equal_var=False)
        rows.append({"feature": column, "t_stat": t_stat, "abs_t": abs(t_stat), "p_value": p_value})
    if not rows:
        return pd.DataFrame(columns=["feature", "t_stat", "abs_t", "p_value"])
    return pd.DataFrame(rows).sort_values("abs_t", ascending=False).head(k).reset_index(drop=True)


def impute_scale_train_val(train_df: pd.DataFrame, val_df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    medians = train_df.median()
    train_imputed = train_df.fillna(medians).fillna(0)
    val_imputed = val_df.fillna(medians).fillna(0)
    scaler = StandardScaler().fit(train_imputed.values)
    return scaler.transform(train_imputed.values), scaler.transform(val_imputed.values)


def loocv_predictions(features: pd.DataFrame, labels: np.ndarray) -> np.ndarray:
    loo = LeaveOneOut()
    predictions = np.zeros(len(labels), dtype=float)
    for train_idx, test_idx in loo.split(features):
        x_train, x_test = impute_scale_train_val(features.iloc[train_idx], features.iloc[test_idx])
        model = LogisticRegression(C=1.0, max_iter=3000, class_weight="balanced")
        model.fit(x_train, labels[train_idx])
        predictions[test_idx] = model.predict_proba(x_test)[:, 1]
    return predictions


def fit_and_score(train_df: pd.DataFrame, train_y: np.ndarray, val_df: pd.DataFrame, val_y: np.ndarray) -> tuple[np.ndarray, float, float]:
    x_train, x_val = impute_scale_train_val(train_df, val_df)
    model = LogisticRegression(C=1.0, max_iter=3000, class_weight="balanced")
    model.fit(x_train, train_y)
    val_pred = model.predict_proba(x_val)[:, 1]
    val_auc = roc_auc_score(val_y, val_pred)
    val_bal_acc = balanced_accuracy_score(val_y, (val_pred >= 0.5).astype(int))
    return val_pred, val_auc, val_bal_acc


def main() -> None:
    args = parse_args()
    model_tag = sanitize_model_name(args.model_variant)

    tables_dir = PROJECT / "out" / "05_multilayer" / "tables"
    figures_dir = PROJECT / "out" / "05_multilayer" / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    multilayer_path = tables_dir / "patient_multilayer_features.csv"
    nmf_path = tables_dir / "nmf_programs_patient_activity.csv"
    selected_path = tables_dir / "balanced_selected_features.csv"
    default_feature_path = tables_dir / f"geneformer_patient_features_{model_tag}.csv"
    feature_path = args.feature_file or default_feature_path
    if args.output_tag is not None:
        output_tag = args.output_tag
    elif feature_path.resolve() == default_feature_path.resolve():
        output_tag = model_tag
    else:
        output_tag = feature_path.stem.replace("geneformer_patient_features_", "")

    output_metrics = tables_dir / f"geneformer_layer_metrics_{output_tag}.csv"
    output_selected = tables_dir / f"geneformer_selected_features_{output_tag}.csv"
    output_predictions = tables_dir / f"geneformer_validation_predictions_{output_tag}.csv"
    figure_path = figures_dir / f"05i_geneformer_layer_benchmark_{output_tag}.png"

    print("=" * 70)
    print("[05i] Geneformer L6 benchmark")
    print("=" * 70)
    print(f"Feature file    : {feature_path}")

    base = pd.read_csv(multilayer_path)
    if nmf_path.exists():
        nmf = pd.read_csv(nmf_path)
        nmf_cols = [column for column in nmf.columns if column.startswith("NMF")]
        nmf = nmf.rename(columns={column: f"nmf_{column}" for column in nmf_cols})
        base = base.merge(nmf[["patient_id"] + [f"nmf_{column}" for column in nmf_cols]], on="patient_id", how="left")

    base_has_geneformer = any(column.startswith("gf_cls_") for column in base.columns)
    if base_has_geneformer:
        print("Feature source  : patient_multilayer_features.csv (canonical L6 merge)")
        df = base.copy()
    else:
        if not feature_path.exists():
            raise FileNotFoundError(
                "Geneformer features were not found in the canonical patient table and the requested feature file does not exist."
            )
        print("Feature source  : standalone Geneformer feature file")
        geneformer = pd.read_csv(feature_path)
        df = base.merge(geneformer, on=["patient_id", "response", "dataset"], how="inner")

    discovery = df[df["dataset"].astype(str) != "Cappuyns2023"].copy()
    validation = df[df["dataset"].astype(str) == "Cappuyns2023"].copy()
    # Exclude DeathBeforeImaging patients — not treatment-response evaluable.
    # Matches the n=32 filter in 03c_cappuyns_validation.py and 06b.
    _DEATH_BEFORE_IMAGING = {"pat5", "pat28", "pat29", "pat38"}
    validation = validation[~validation["patient_id"].astype(str).isin(_DEATH_BEFORE_IMAGING)].copy()
    assert len(validation) == 32, f"Expected 32 Cappuyns patients (19R/13NR), got {len(validation)}"
    y_disc = (discovery["response"] == "Responder").astype(int).values
    y_val = (validation["response"] == "Responder").astype(int).values

    if len(discovery) == 0 or len(validation) == 0:
        raise RuntimeError("Expected both discovery and Cappuyns patients after merging Geneformer features.")

    existing_selected = pd.read_csv(selected_path) if selected_path.exists() else pd.DataFrame(columns=["layer", "feature"])
    selected_by_layer: dict[str, list[str]] = {}
    for layer in existing_selected.get("layer", pd.Series(dtype=str)).dropna().unique():
        features = existing_selected.loc[existing_selected["layer"] == layer, "feature"].tolist()
        features = [feature for feature in features if feature in discovery.columns and feature in validation.columns]
        if features:
            selected_by_layer[layer] = features

    geneformer_selected = pick_top_features(discovery, y_disc, prefix="gf_cls_", k=args.k)
    if geneformer_selected.empty:
        raise RuntimeError("No valid Geneformer dimensions were available for selection.")
    geneformer_selected.insert(0, "layer", "L6_geneformer")
    geneformer_selected.insert(1, "analysis_cohort", "GSE206325_discovery")
    geneformer_selected.to_csv(output_selected, index=False)
    selected_by_layer["L6_geneformer"] = geneformer_selected["feature"].tolist()

    print(f"Discovery cohort : {len(discovery)} patients")
    print(f"Validation cohort: {len(validation)} patients")
    print(f"L6 features      : {selected_by_layer['L6_geneformer']}")

    metric_rows = []
    prediction_rows = []
    layer_order = [
        "L1_tumour_intrinsic",
        "L2_regulon",
        "L3_composition",
        "L4_io_signatures",
        "L5_nmf_programs",
        "L6_geneformer",
    ]

    for layer in layer_order:
        if layer not in selected_by_layer:
            continue
        features = selected_by_layer[layer]
        disc_pred = loocv_predictions(discovery[features], y_disc)
        disc_auc = roc_auc_score(y_disc, disc_pred)
        val_pred, val_auc, val_bal_acc = fit_and_score(discovery[features], y_disc, validation[features], y_val)
        val_prauc = average_precision_score(y_val, val_pred)
        metric_rows.append(
            {
                "layer": layer,
                "n_features": len(features),
                "disc_loocv_auc": disc_auc,
                "val_auc": val_auc,
                "val_prauc": val_prauc,
                "val_balanced_acc": val_bal_acc,
            }
        )
        for patient_id, response, score in zip(validation["patient_id"], validation["response"], val_pred):
            prediction_rows.append(
                {
                    "layer": layer,
                    "patient_id": patient_id,
                    "response": response,
                    "score": score,
                }
            )

    ensemble_features = sum([selected_by_layer[layer] for layer in layer_order if layer in selected_by_layer], [])
    disc_pred = loocv_predictions(discovery[ensemble_features], y_disc)
    disc_auc = roc_auc_score(y_disc, disc_pred)
    val_pred, val_auc, val_bal_acc = fit_and_score(discovery[ensemble_features], y_disc, validation[ensemble_features], y_val)
    val_prauc = average_precision_score(y_val, val_pred)
    metric_rows.append(
        {
            "layer": "Balanced_plus_Geneformer",
            "n_features": len(ensemble_features),
            "disc_loocv_auc": disc_auc,
            "val_auc": val_auc,
            "val_prauc": val_prauc,
            "val_balanced_acc": val_bal_acc,
        }
    )
    for patient_id, response, score in zip(validation["patient_id"], validation["response"], val_pred):
        prediction_rows.append(
            {
                "layer": "Balanced_plus_Geneformer",
                "patient_id": patient_id,
                "response": response,
                "score": score,
            }
        )

    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(output_metrics, index=False)
    pd.DataFrame(prediction_rows).to_csv(output_predictions, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)
    plot_metrics = metrics.copy()
    plot_metrics["label"] = plot_metrics["layer"].str.replace("_", " ", regex=False)

    axes[0].barh(plot_metrics["label"], plot_metrics["disc_loocv_auc"], color="#4477AA")
    axes[0].set_xlim(0.0, 1.0)
    axes[0].set_xlabel("Discovery LOOCV AUC")
    axes[0].set_title("GSE206325", fontsize=11, fontweight="bold")

    axes[1].barh(plot_metrics["label"], plot_metrics["val_auc"], color="#228833")
    axes[1].set_xlim(0.0, 1.0)
    axes[1].set_xlabel("Cappuyns validation AUC")
    axes[1].set_title("External transfer", fontsize=11, fontweight="bold")
    axes[1].tick_params(axis="y", labelleft=False)

    fig.suptitle("Geneformer as L6 patient representation", fontsize=12, fontweight="bold")
    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(metrics.to_string(index=False))
    print(f"Saved metrics   : {output_metrics}")
    print(f"Saved features  : {output_selected}")
    print(f"Saved figure    : {figure_path}")
    print("[DONE]")


if __name__ == "__main__":
    main()