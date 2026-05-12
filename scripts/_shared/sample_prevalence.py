#!/usr/bin/env python3
"""
sample_prevalence.py — generic per-feature sample/patient prevalence helper.

Purpose
-------
Many results in this project are reported as cohort-level effects
(e.g. "CNV clone diversity is higher in NR"). The reviewer-relevant question
is then: *what proportion of patients actually exhibit the effect, and how
do those carrying patients differ from the rest of the cohort?*

This module provides one entry point — :func:`prevalence_panel` — that:

1. For each feature in ``features``, classifies every patient as
   "carrier" (above the cohort threshold, default Q3) vs "other".
2. Computes the carrier fraction, both overall and within ``group_col``
   (typically ``response``).
3. Runs Mann–Whitney U tests of carrier vs other on every column in
   ``compare_cols`` — those are the molecular / clinical / compositional
   properties to compare across the carrier subset.
4. Persists a tidy CSV (``<tag>_prevalence.csv``) and renders a 3-panel
   figure (``<tag>_prevalence.png/pdf``):

     A. carrier fraction by group (proportion of patients showing the result)
     B. per-feature carrier-vs-other distributions of a chosen property
     C. heatmap of effect sizes / p-values across features × properties

The helper is intentionally tolerant: it skips features that are missing,
constant, or have <3 carriers.

Usage
-----
>>> from scripts._shared.sample_prevalence import prevalence_panel
>>> prevalence_panel(
...     patient_df,
...     features=["cnv_clone_diversity", "NMF10", "CD8_score"],
...     compare_cols=["malig_frac", "tls_score", "treg_frac"],
...     group_col="response",
...     out_dir=TAB,
...     fig_dir=FIG,
...     tag="05a_multilayer",
... )
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu


def _carrier_mask(values: pd.Series, threshold: str | float = "q75",
                  direction: str = "high") -> pd.Series:
    """Return a boolean carrier mask using a quantile or absolute threshold."""
    s = pd.to_numeric(values, errors="coerce")
    if isinstance(threshold, str) and threshold.startswith("q"):
        q = float(threshold[1:]) / 100.0
        cut = s.quantile(q)
    else:
        cut = float(threshold)
    if direction == "low":
        return s <= cut
    return s >= cut


def prevalence_panel(
    patient_df: pd.DataFrame,
    features: Sequence[str],
    compare_cols: Sequence[str] = (),
    group_col: str | None = "response",
    threshold: str | float = "q75",
    direction: str = "high",
    out_dir: Path | str | None = None,
    fig_dir: Path | str | None = None,
    tag: str = "prevalence",
    palette: dict | None = None,
    title: str | None = None,
) -> pd.DataFrame:
    """Compute per-feature carrier prevalence + carrier-vs-other comparisons.

    Returns a tidy DataFrame keyed by ``(feature, property)`` with columns
    ``n_carrier, n_other, frac_carrier, frac_carrier_<group>, mean_carrier,
    mean_other, effect_size, U, p_value``. Also writes a CSV + figure when
    ``out_dir`` / ``fig_dir`` are provided.
    """
    df = patient_df.copy()
    feats = [f for f in features if f in df.columns]
    if not feats:
        return pd.DataFrame()

    if palette is None:
        palette = {"carrier": "#117733", "other": "#BBBBBB"}

    rows: list[dict] = []
    carrier_masks: dict[str, pd.Series] = {}
    for feat in feats:
        mask = _carrier_mask(df[feat], threshold=threshold, direction=direction)
        if mask.sum() < 3 or (~mask).sum() < 3:
            continue
        carrier_masks[feat] = mask

        base = {
            "feature": feat,
            "threshold_kind": str(threshold),
            "direction": direction,
            "n_carrier": int(mask.sum()),
            "n_other": int((~mask).sum()),
            "frac_carrier": float(mask.mean()),
        }
        if group_col and group_col in df.columns:
            grp = df.groupby(group_col, observed=True)
            for level, sub_idx in grp.indices.items():
                sub_mask = mask.iloc[sub_idx] if isinstance(sub_idx, np.ndarray) else mask.loc[sub_idx]
                base[f"frac_carrier_{level}"] = float(np.asarray(sub_mask).mean())

        # Self-stats: carrier feature value vs others
        carrier_vals = pd.to_numeric(df.loc[mask, feat], errors="coerce").dropna()
        other_vals = pd.to_numeric(df.loc[~mask, feat], errors="coerce").dropna()
        base["mean_carrier_self"] = float(carrier_vals.mean()) if len(carrier_vals) else np.nan
        base["mean_other_self"] = float(other_vals.mean()) if len(other_vals) else np.nan

        # Compare against each property
        cmp_cols = [c for c in compare_cols if c in df.columns and c != feat]
        if not cmp_cols:
            rows.append({**base, "property": None})
        else:
            for prop in cmp_cols:
                vals = pd.to_numeric(df[prop], errors="coerce")
                a = vals[mask].dropna()
                b = vals[~mask].dropna()
                if len(a) < 3 or len(b) < 3:
                    continue
                pooled_std = np.std(np.concatenate([a, b]), ddof=1)
                eff = (a.mean() - b.mean()) / (pooled_std + 1e-12)
                try:
                    U, p = mannwhitneyu(a, b, alternative="two-sided")
                except ValueError:
                    U, p = np.nan, np.nan
                rows.append({
                    **base,
                    "property": prop,
                    "mean_carrier": float(a.mean()),
                    "mean_other": float(b.mean()),
                    "effect_size": float(eff),
                    "U": float(U) if np.isfinite(U) else np.nan,
                    "p_value": float(p) if np.isfinite(p) else np.nan,
                })

    if not rows:
        return pd.DataFrame()

    summary = pd.DataFrame(rows)
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        summary.to_csv(out_dir / f"{tag}_prevalence.csv", index=False)

    if fig_dir is not None and carrier_masks:
        fig_dir = Path(fig_dir)
        fig_dir.mkdir(parents=True, exist_ok=True)
        _render_panel(df, carrier_masks, summary, group_col, compare_cols,
                      fig_dir / tag, palette, title)
    return summary


def _render_panel(df, carrier_masks, summary, group_col, compare_cols,
                  fig_stem: Path, palette: dict, title: str | None) -> None:
    feats = list(carrier_masks)
    n_feats = len(feats)

    fig, axes = plt.subplots(1, 3, figsize=(min(22, 6 + 1.0 * n_feats), 5.5))

    # ── Panel A: carrier fraction (overall + per-group) ──────────────
    ax = axes[0]
    overall = pd.Series({f: carrier_masks[f].mean() for f in feats}, name="overall")
    bars = pd.DataFrame({"overall": overall})
    if group_col and group_col in df.columns:
        for level, sub in df.groupby(group_col, observed=True):
            bars[str(level)] = pd.Series({f: carrier_masks[f].loc[sub.index].mean()
                                          for f in feats})
    bars = bars.sort_values("overall", ascending=True)
    bars.plot(kind="barh", ax=ax, width=0.8,
              color=["#444444"] + sns.color_palette("Set2", max(0, bars.shape[1] - 1)))
    ax.set_xlim(0, 1)
    ax.set_xlabel("Proportion of patients (carrier)")
    ax.set_title("A. Fraction of patients showing each result", fontweight="bold")
    ax.legend(fontsize=7, frameon=False)

    # ── Panel B: feature value distribution carrier vs other for top feature ─
    ax = axes[1]
    top_feat = bars.index[-1]  # largest fraction
    plot_df = df[[top_feat]].copy()
    plot_df["group"] = np.where(carrier_masks[top_feat], "carrier", "other")
    sns.boxplot(data=plot_df, x="group", y=top_feat, ax=ax, hue="group",
                palette={"carrier": palette["carrier"], "other": palette["other"]},
                width=0.55, fliersize=0, linewidth=1, legend=False)
    sns.stripplot(data=plot_df, x="group", y=top_feat, ax=ax,
                  color="black", size=3.2, alpha=0.7, jitter=0.18)
    ax.set_title(f"B. {top_feat}: carrier vs other", fontweight="bold")
    ax.set_xlabel("")

    # ── Panel C: effect-size heatmap features × compare_cols ─────────
    ax = axes[2]
    cmp_cols = [c for c in compare_cols if c in df.columns]
    if cmp_cols and "property" in summary.columns:
        eff_mat = (summary.dropna(subset=["property"])
                          .pivot_table(index="feature", columns="property",
                                       values="effect_size", aggfunc="first"))
        if eff_mat.size > 0:
            sns.heatmap(eff_mat, cmap="RdBu_r", center=0, ax=ax,
                        annot=True, fmt=".2f", linewidths=0.4,
                        cbar_kws=dict(label="Cohen-style d", shrink=0.6))
            ax.set_title("C. Carrier vs other effect across properties",
                         fontweight="bold")
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.tick_params(axis="x", rotation=30, labelsize=8)
            ax.tick_params(axis="y", labelsize=8)
        else:
            ax.axis("off")
    else:
        ax.axis("off")

    if title:
        fig.suptitle(title, fontweight="bold", y=1.02)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{fig_stem}_prevalence.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


__all__ = ["prevalence_panel"]
