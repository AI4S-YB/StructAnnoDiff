#!/usr/bin/env python3
"""Create a vertical single-species curation summary figure."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd

from analysis_config import ANALYSIS, FIGURES_DIR, RESULTS_DIR, find_annotation


MPLCONFIGDIR = ANALYSIS / ".cache" / "matplotlib"
MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REQUIRED_COLUMNS = [
    "species_id",
    "Species",
    "total_before_genes",
    "total_after_genes",
    "new_loci_no_overlap",
    "deleted_loci_no_overlap",
    "split_events",
    "merge_events",
    "rep_exon_changed",
    "rep_exon_changed_before_pct",
    "rep_exon_changed_after_pct",
]

COLORS = {
    "deleted": "#B65A5A",
    "new": "#4F8DBA",
    "split": "#7A6BB0",
    "merge": "#C78C3D",
    "file1": "#4D4D4D",
    "file2": "#2F6C99",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot a vertical single-species curation summary figure."
    )
    parser.add_argument(
        "--input",
        default=str(RESULTS_DIR / "curation_core_metrics.csv"),
        help="Input curation core metrics CSV.",
    )
    parser.add_argument("--species-id", default="Pineapple", help="Species ID to plot.")
    parser.add_argument(
        "--output-prefix",
        default=str(FIGURES_DIR / "Pineapple_curation_summary_vertical"),
        help="Output path prefix without extension.",
    )
    parser.add_argument("--dpi", type=int, default=600, help="PNG output resolution.")
    return parser.parse_args()


def load_row(path: Path, species_id: str) -> pd.Series:
    df = pd.read_csv(path)
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        raise SystemExit(f"{path}: missing required columns: {', '.join(missing)}")
    rows = df[df["species_id"] == species_id]
    if rows.empty:
        raise SystemExit(f"{path}: species_id not found: {species_id}")
    if len(rows) > 1:
        raise SystemExit(f"{path}: species_id is not unique: {species_id}")
    row = rows.iloc[0]
    for col in ("rep_exon_changed_before_pct", "rep_exon_changed_after_pct"):
        if row[col] < 0 or row[col] > 100:
            raise SystemExit(f"{path}: {col} must be within 0-100")
    return row


def fmt_count(value: float) -> str:
    return f"{int(round(value)):,}"


def fmt_percent(pct: float, count: int, denom: int) -> str:
    return f"{pct:.2f}% ({fmt_count(count)}/{fmt_count(denom)})"


def style_panel(ax, x_max: float, xlabel: str):
    ax.set_xlim(0, x_max)
    ax.set_ylim(-0.65, 1.65)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.7)
    ax.tick_params(axis="y", length=0, labelsize=8.8)
    ax.tick_params(axis="x", labelsize=8.2, width=0.6, length=3)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.65)
    ax.set_axisbelow(True)
    ax.set_xlabel(xlabel, fontsize=8.8)


def panel_header(ax, label: str, title: str):
    ax.text(
        0.0,
        1.10,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11.0,
        fontweight="bold",
    )
    ax.text(
        0.07,
        1.10,
        title,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10.0,
        fontweight="bold",
    )


def annotate(ax, x: float, y: float, text: str, x_max: float):
    pad = x_max * 0.018
    if x > x_max * 0.72:
        ax.text(x - pad, y, text, ha="right", va="center", fontsize=8.2, color="white")
    else:
        ax.text(x + pad, y, text, ha="left", va="center", fontsize=8.2, color="#222222")


def draw_pair_panel(ax, labels, values, colors, x_max, xlabel, header_label, header_title):
    y = [1, 0]
    ax.barh(y, values, height=0.44, color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    style_panel(ax, x_max, xlabel)
    panel_header(ax, header_label, header_title)
    for yi, value in zip(y, values):
        annotate(ax, value, yi, fmt_count(value), x_max)


def plot(row: pd.Series, species_id: str, output_prefix: Path, dpi: int):
    label = str(row["Species"])
    before = find_annotation(species_id, "before")
    after = find_annotation(species_id, "after")
    before_name = before.name if before else f"{species_id}.before"
    after_name = after.name if after else f"{species_id}.after"

    deleted = int(row["deleted_loci_no_overlap"])
    new = int(row["new_loci_no_overlap"])
    split = int(row["split_events"])
    merge = int(row["merge_events"])
    exon_changed = int(row["rep_exon_changed"])
    before_total = int(row["total_before_genes"])
    after_total = int(row["total_after_genes"])
    before_pct = float(row["rep_exon_changed_before_pct"])
    after_pct = float(row["rep_exon_changed_after_pct"])

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(6.2, 6.9),
        constrained_layout=False,
        gridspec_kw={"height_ratios": [1, 1, 1.08], "hspace": 0.78},
    )
    fig.patch.set_facecolor("white")
    fig.suptitle(
        f"{label} annotation curation summary",
        x=0.08,
        y=0.985,
        ha="left",
        va="top",
        fontsize=13.0,
        fontweight="bold",
    )
    fig.text(
        0.08,
        0.945,
        f"{before_name} -> {after_name}",
        ha="left",
        va="top",
        fontsize=8.6,
        color="#444444",
    )

    draw_pair_panel(
        axes[0],
        ["Deleted loci", "New loci"],
        [deleted, new],
        [COLORS["deleted"], COLORS["new"]],
        max(deleted, new, 1) * 1.24,
        "No-overlap loci",
        "A",
        "Locus gain/loss",
    )
    draw_pair_panel(
        axes[1],
        ["Split events", "Merge events"],
        [split, merge],
        [COLORS["split"], COLORS["merge"]],
        max(split, merge, 1) * 1.24,
        "Events",
        "B",
        "Split and merge events",
    )

    y = [1, 0]
    pct_values = [before_pct, after_pct]
    axes[2].barh(y, pct_values, height=0.44, color=[COLORS["file1"], COLORS["file2"]])
    axes[2].set_yticks(y)
    axes[2].set_yticklabels(["Input file 1", "Input file 2"])
    style_panel(axes[2], 100, "Representative exon changes (%)")
    panel_header(axes[2], "C", "Representative transcript exon changes")
    for yi, pct, denom in zip(y, pct_values, [before_total, after_total]):
        annotate(axes[2], pct, yi, fmt_percent(pct, exon_changed, denom), 100)

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(left=0.22, right=0.96, top=0.87, bottom=0.08)
    fig.savefig(output_prefix.with_suffix(".png"), dpi=dpi)
    fig.savefig(output_prefix.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    args = parse_args()
    row = load_row(Path(args.input), args.species_id)
    plot(row, args.species_id, Path(args.output_prefix), args.dpi)


if __name__ == "__main__":
    main()
