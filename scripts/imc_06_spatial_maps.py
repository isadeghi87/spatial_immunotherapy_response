#!/usr/bin/env python3
"""
imc_06_spatial_maps.py

Spatial cell-map figure for IMC validation.

Layout: 3 columns (Enriched / Compartmentalized / Depleted representative ROI)
        × 3 rows:
  Row 1 — Cell-type map       : all cells coloured by phenotype class
  Row 2 — Compartment overlay : stroma/parenchyma background + CD8 highlighted
  Row 3 — CD8 functional state: CD8 T cells only, coloured by PD1/GzmB/TCF1

Representative ROIs are chosen as the largest ROI from the patient whose
patient-level CD8 density is closest to the immunotype median (rank-based).

Outputs → out/03_validation/imc_validation/figures/
  imc_spatial_maps.{png,pdf}
"""

from __future__ import annotations
import os, sys, warnings
warnings.filterwarnings("ignore")

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR))
from config.publication_style import apply_style, PAL

apply_style()

TABLE_DIR = PROJECT_DIR / "out" / "03_validation" / "imc_validation" / "tables"
FIG_DIR   = PROJECT_DIR / "out" / "03_validation" / "imc_validation" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Palettes ───────────────────────────────────────────────────────────────
IMMUNOTYPE_ORDER  = ["Enriched", "Compartmentalized", "Depleted"]
IMMUNOTYPE_COLORS = {
    "Enriched":          "#0077BB",
    "Compartmentalized": "#EE7733",
    "Depleted":          "#CC3311",
}

CELLTYPE_PALETTE = {
    "CD8_T":        "#0077BB",
    "CD4_T":        "#33BBEE",
    "Treg":         "#AA3377",
    "NK":           "#009988",
    "B_cell":       "#7B4F9E",
    "Macrophage":   "#E69F00",
    "APC_DC":       "#F0E442",
    "Tumor":        "#CC3311",
    "Stroma_other": "#BBBBBB",
    "Other_immune": "#44BB99",
}

COMPARTMENT_COLORS = {
    "parenchyma": "#FADADD",   # pale rose
    "stroma":     "#D9EAD3",   # pale green
}

CD8_FUNC_STATES = ["PD1+", "GzmB+", "TCF1+", "Other"]
CD8_FUNC_COLORS = {
    "PD1+":  "#CC3311",   # exhausted/activated
    "GzmB+": "#0077BB",   # cytotoxic
    "TCF1+": "#009988",   # progenitor-like
    "Other": "#CCCCCC",
}

ROW_LABELS = [
    "Cell-type map",
    "Compartment  +  CD8",
    "CD8 functional state",
]


CROP_SIZE_PX = 1000   # crop all ROIs to this square (≈ 1 mm × 1 mm at 1µm/px)


# ── Helper: pick representative patient+ROI ────────────────────────────────
def pick_representative_roi(
    pf: pd.DataFrame,
    ct: pd.DataFrame,
    immunotype: str,
) -> tuple[str, str]:
    """Return (patient_id, roi_id) for the patient closest to median CD8
    density in the given immunotype, using the ROI with the most cells."""
    sub = pf[pf["immunotype_paper"] == immunotype].copy()
    sub = sub.sort_values("CD8_density_mm2")
    med_idx = len(sub) // 2
    pid = sub.iloc[med_idx]["patient_id"]
    roi_sizes = (
        ct[ct["patient_id"] == pid]
        .groupby("roi_id")["cell_id"]
        .count()
    )
    roi_id = roi_sizes.idxmax()
    return pid, roi_id


def crop_to_square(cells: pd.DataFrame, size: int = CROP_SIZE_PX) -> pd.DataFrame:
    """Crop cells to a centred square of `size` px.
    Uses the centroid of all cells as the crop centre so tissue edges are avoided."""
    cx = cells["centroid_x"].mean()
    cy = cells["centroid_y"].mean()
    half = size / 2
    mask = (
        (cells["centroid_x"] >= cx - half) & (cells["centroid_x"] <= cx + half) &
        (cells["centroid_y"] >= cy - half) & (cells["centroid_y"] <= cy + half)
    )
    return cells.loc[mask].copy()


# ── Row 1: Cell-type map ───────────────────────────────────────────────────
def plot_celltype_map(ax: plt.Axes, cells: pd.DataFrame) -> None:
    order = list(CELLTYPE_PALETTE.keys())
    # Plot non-stroma first (background), then immune on top
    bg_types = ["Stroma_other", "Tumor"]
    fg_types  = [t for t in order if t not in bg_types]

    for ct_group, pt_size, alpha in [(bg_types, 3, 0.25), (fg_types, 5, 0.80)]:
        for ctype in ct_group:
            mask = cells["cell_type"] == ctype
            if not mask.any():
                continue
            ax.scatter(
                cells.loc[mask, "centroid_x"],
                cells.loc[mask, "centroid_y"],
                c=CELLTYPE_PALETTE.get(ctype, "#888888"),
                s=pt_size, alpha=alpha,
                linewidths=0, rasterized=True,
            )

    ax.axis("off")


# ── Row 2: Compartment overlay + CD8 ──────────────────────────────────────
def plot_compartment_map(ax: plt.Axes, cells: pd.DataFrame) -> None:
    # Background: stroma / parenchyma
    for comp, color in COMPARTMENT_COLORS.items():
        mask = cells["compartment"] == comp
        if not mask.any():
            continue
        ax.scatter(
            cells.loc[mask, "centroid_x"],
            cells.loc[mask, "centroid_y"],
            c=color, s=4, alpha=0.55, linewidths=0, rasterized=True,
        )

    # CD8 T cells on top, coloured darker
    cd8 = cells[cells["cell_type"] == "CD8_T"]
    if len(cd8):
        ax.scatter(
            cd8["centroid_x"], cd8["centroid_y"],
            c=CELLTYPE_PALETTE["CD8_T"], s=12, alpha=0.90,
            linewidths=0.3, edgecolors="white", rasterized=True, zorder=3,
        )

    ax.axis("off")


# ── Row 3: CD8 functional state ────────────────────────────────────────────
def assign_cd8_func_state(cells: pd.DataFrame) -> pd.Series:
    """Priority: PD1+ > GzmB+ > TCF1+ > Other."""
    cd8 = cells[cells["cell_type"] == "CD8_T"].copy()
    state = pd.Series("Other", index=cd8.index)
    state[cd8["CD8_TCF1pos"]]  = "TCF1+"
    state[cd8["CD8_GzmBpos"]]  = "GzmB+"
    state[cd8["CD8_PD1pos"]]   = "PD1+"
    return state


def plot_cd8_functional_map(ax: plt.Axes, cells: pd.DataFrame) -> None:
    # Greyed background: non-CD8 cells
    non_cd8 = cells[cells["cell_type"] != "CD8_T"]
    ax.scatter(
        non_cd8["centroid_x"], non_cd8["centroid_y"],
        c="#E8E8E8", s=2, alpha=0.15, linewidths=0, rasterized=True,
    )

    cd8 = cells[cells["cell_type"] == "CD8_T"].copy()
    if len(cd8) == 0:
        ax.axis("off")
        return

    func_state = assign_cd8_func_state(cells)
    for state in ["Other", "TCF1+", "GzmB+", "PD1+"]:
        idx = func_state[func_state == state].index
        if len(idx) == 0:
            continue
        sub = cd8.loc[cd8.index.intersection(idx)]
        ax.scatter(
            sub["centroid_x"], sub["centroid_y"],
            c=CD8_FUNC_COLORS[state],
            s=9 if state != "Other" else 4,
            alpha=0.90 if state != "Other" else 0.35,
            linewidths=0.3 if state != "Other" else 0,
            edgecolors="white" if state != "Other" else "none",
            rasterized=True, zorder=3 if state != "Other" else 1,
        )

    ax.axis("off")


# ── Shared scale-bar ───────────────────────────────────────────────────────
def add_scale_bar(ax: plt.Axes, cells: pd.DataFrame,
                  bar_um: float = 200.0, px_per_um: float = 1.0) -> None:
    """200 µm scale bar, bottom-left, using axes-fraction coordinates so
    it is independent of data-range and axis-inversion direction."""
    bar_frac = (bar_um * px_per_um) / CROP_SIZE_PX  # fraction of axes width
    # Position: 3% from left, 5% from bottom (in axes coords)
    ax.annotate(
        "", xy=(0.03 + bar_frac, 0.05), xytext=(0.03, 0.05),
        xycoords="axes fraction", textcoords="axes fraction",
        arrowprops=dict(arrowstyle="-", color="black", lw=1.8),
    )
    ax.text(0.03 + bar_frac / 2, 0.08, f"{int(bar_um)} µm",
            ha="center", va="bottom", fontsize=6.5, color="black",
            transform=ax.transAxes,
            path_effects=[pe.withStroke(linewidth=1.5, foreground="white")])


# ── Legend builders ────────────────────────────────────────────────────────
def celltype_legend_handles() -> list:
    handles = []
    for ct, col in CELLTYPE_PALETTE.items():
        handles.append(
            Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=col, markersize=6, label=ct.replace("_", " "))
        )
    return handles


def compartment_legend_handles() -> list:
    handles = [
        Patch(facecolor=COMPARTMENT_COLORS["parenchyma"], label="Parenchyma"),
        Patch(facecolor=COMPARTMENT_COLORS["stroma"],     label="Stroma"),
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=CELLTYPE_PALETTE["CD8_T"], markersize=6, label="CD8 T"),
    ]
    return handles


def cd8func_legend_handles() -> list:
    return [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=CD8_FUNC_COLORS[s], markersize=6, label=s)
        for s in CD8_FUNC_STATES
    ]


# ── Main figure assembly ───────────────────────────────────────────────────
def build_spatial_figure(
    ct_all: pd.DataFrame,
    pf: pd.DataFrame,
) -> plt.Figure:

    n_rows, n_cols = 3, 3
    # Each panel: aspect-free axes, we'll let invert+equal handle it
    fig_w = n_cols * 3.8 + 1.5   # +1.5 for row-label margin
    fig_h = n_rows * 3.6 + 1.2   # +1.2 for col headers

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")

    # Grid: leave left margin for row labels, top margin for col headers
    gs = fig.add_gridspec(
        n_rows, n_cols,
        left=0.10, right=0.88, top=0.90, bottom=0.04,
        hspace=0.06, wspace=0.04,
    )

    row_legend_ax = [None] * n_rows
    axes_grid = [[None] * n_cols for _ in range(n_rows)]

    for r in range(n_rows):
        for c in range(n_cols):
            axes_grid[r][c] = fig.add_subplot(gs[r, c])

    # ── Pick representatives ──────────────────────────────────────────
    reps = {}   # immunotype → (pid, roi_id, cells_df)
    for it in IMMUNOTYPE_ORDER:
        pid, roi_id = pick_representative_roi(pf, ct_all, it)
        cells = ct_all[
            (ct_all["patient_id"] == pid) & (ct_all["roi_id"] == roi_id)
        ].copy()
        cells = crop_to_square(cells, CROP_SIZE_PX)
        reps[it] = (pid, roi_id, cells)
        print(f"  {it:<20s}: {pid} / {roi_id}  ({len(cells):,} cells)")

    # ── Per-column metadata for header labels ─────────────────────────
    meta = {}
    for it in IMMUNOTYPE_ORDER:
        pid, roi_id, cells = reps[it]
        prow = pf[pf["patient_id"] == pid].iloc[0]
        n_cd8  = (cells["cell_type"] == "CD8_T").sum()
        cd8_d  = prow["CD8_density_mm2"]
        resp   = "R" if prow["response_binary"] == 1 else "NR"
        meta[it] = dict(pid=pid, roi=roi_id, n_cells=len(cells),
                        n_cd8=n_cd8, cd8_density=cd8_d, response=resp)

    # ── Draw panels ───────────────────────────────────────────────────
    plot_fns = [plot_celltype_map, plot_compartment_map, plot_cd8_functional_map]

    for col_idx, it in enumerate(IMMUNOTYPE_ORDER):
        pid, roi_id, cells = reps[it]
        m = meta[it]

        # Shared axis limits for this column (the cropped square)
        xmid = cells["centroid_x"].mean()
        ymid = cells["centroid_y"].mean()
        half = CROP_SIZE_PX / 2
        xlim = (xmid - half, xmid + half)
        ylim = (ymid - half, ymid + half)

        for row_idx, plot_fn in enumerate(plot_fns):
            ax = axes_grid[row_idx][col_idx]
            plot_fn(ax, cells)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim[1], ylim[0])   # invert y: smaller y = top
            add_scale_bar(ax, cells)

            # Immunotype-coloured border on all rows for visual grouping
            for sp in ["top", "bottom", "left", "right"]:
                ax.spines[sp].set_visible(True)
                ax.spines[sp].set_linewidth(1.5 if row_idx > 0 else 2.5)
                ax.spines[sp].set_edgecolor(IMMUNOTYPE_COLORS[it])

        # Column header via set_title on top-row axis (respects aspect=equal shifts)
        ax_top = axes_grid[0][col_idx]
        ax_top.set_title(
            f"{it}\n"
            f"{m['pid']} / {m['roi']}  ·  {m['response']}\n"
            f"CD8: {m['cd8_density']:.0f} cells/mm²  ·  n={m['n_cells']:,}",
            fontsize=9, fontweight="bold",
            color=IMMUNOTYPE_COLORS[it],
            pad=6, loc="center",
            linespacing=1.5,
        )

    # ── Row labels (left side via ylabel on leftmost axes) ────────────
    for row_idx, rlabel in enumerate(ROW_LABELS):
        ax_left = axes_grid[row_idx][0]
        ax_left.set_ylabel(rlabel, fontsize=9, fontweight="bold",
                           labelpad=8, rotation=90, va="center")
        # Actually ylabel won't show with axis('off'); use fig.text instead
        ax_left.set_ylabel("")  # clear it

    # Re-add row labels using fig.text after layout
    fig.canvas.draw()  # force layout so get_position() is accurate
    for row_idx, rlabel in enumerate(ROW_LABELS):
        ax_left = axes_grid[row_idx][0]
        pos = ax_left.get_position()
        ry = pos.y0 + pos.height / 2
        rx = pos.x0 - 0.006
        fig.text(rx, ry, rlabel,
                 ha="right", va="center", fontsize=8.5, fontweight="bold",
                 rotation=90, transform=fig.transFigure)

    # ── Legends (right side) ──────────────────────────────────────────
    legend_specs = [
        (celltype_legend_handles(),    axes_grid[0][2], "Cell type"),
        (compartment_legend_handles(), axes_grid[1][2], "Compartment"),
        (cd8func_legend_handles(),     axes_grid[2][2], "CD8 state"),
    ]
    for handles, ref_ax, title in legend_specs:
        ref_pos = ref_ax.get_position()
        lx = 0.895
        ly = ref_pos.y0 + ref_pos.height / 2
        leg = fig.legend(
            handles=handles, loc="center left",
            bbox_to_anchor=(lx, ly),
            bbox_transform=fig.transFigure,
            fontsize=6.5, frameon=True, framealpha=0.85,
            edgecolor="0.80", title=title,
            title_fontsize=7,
        )

    # ── Suptitle ─────────────────────────────────────────────────────
    n_total = pf["immunotype_paper"].value_counts()
    ns_str = "  ".join(
        [f"{it[0]}:n={n_total.get(it,0)}" for it in IMMUNOTYPE_ORDER]
    )
    fig.suptitle(
        f"IMC spatial maps by immunotype  ·  representative ROIs\n"
        f"{ns_str}   |   centroid_x/y from segmentation masks",
        fontsize=10, fontweight="bold", y=0.97, linespacing=1.6,
    )

    return fig


# ── Entry point ────────────────────────────────────────────────────────────
def main():
    print("Loading data ...")
    ct_all = pd.read_parquet(TABLE_DIR / "imc_cell_typed.parquet")
    pf     = pd.read_csv(TABLE_DIR / "imc_patient_features.csv")
    pf     = pf[pf["immunotype_paper"].notna()].copy()

    print(f"  {len(ct_all):,} cells  |  {ct_all['patient_id'].nunique()} patients")
    print(f"  Immunotype counts: {pf['immunotype_paper'].value_counts().to_dict()}")

    print("\nBuilding spatial figure ...")
    fig = build_spatial_figure(ct_all, pf)

    out_stem = str(FIG_DIR / "imc_spatial_maps")
    for ext in ("png", "pdf"):
        fig.savefig(f"{out_stem}.{ext}", dpi=200, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print(f"\nSaved: {out_stem}.{{png,pdf}}")


if __name__ == "__main__":
    main()
