#!/usr/bin/env python3
"""Create a publication-style three-panel figure for core curation metrics."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd

from analysis_config import ANALYSIS, FIGURES_DIR, RESULTS_DIR, SPECIES


MPLCONFIGDIR = ANALYSIS / ".cache" / "matplotlib"
MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


REQUIRED_COLUMNS = [
    "species_id",
    "Species",
    "new_loci_no_overlap",
    "deleted_loci_no_overlap",
    "split_events",
    "merge_events",
    "total_before_genes",
    "total_after_genes",
    "rep_exon_changed",
    "rep_exon_changed_before_pct",
    "rep_exon_changed_after_pct",
]

COLORS = {
    "before_ref": "#4D4D4D",
    "after_ref": "#2F6C99",
    "deleted": "#B65A5A",
    "new": "#4F8DBA",
    "split": "#7A6BB0",
    "merge": "#C78C3D",
    "representative_change": "#2F6C99",
    "exon_count": "#3B7F5E",
    "exon_boundary": "#8FBF9F",
    "cds_count": "#B06A3C",
    "cds_boundary": "#D9A66A",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot a publication-style horizontal-bar figure from curation_core_metrics.csv."
    )
    parser.add_argument(
        "--input",
        default=str(RESULTS_DIR / "curation_core_metrics.csv"),
        help="Input curation core metrics CSV.",
    )
    parser.add_argument(
        "--output-prefix",
        default=str(FIGURES_DIR / "curation_core_metrics_publication"),
        help="Output path prefix without extension.",
    )
    parser.add_argument("--dpi", type=int, default=600, help="PNG output resolution.")
    return parser.parse_args()


def load_metrics(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        raise SystemExit(f"{path}: missing required columns: {', '.join(missing)}")

    order = {sp.id: i for i, sp in enumerate(SPECIES)}
    observed = set(df["species_id"])
    expected = set(order)
    if observed != expected:
        raise SystemExit(
            f"{path}: species mismatch; missing={sorted(expected - observed)}, "
            f"extra={sorted(observed - expected)}"
        )

    df = df.copy()
    df["_order"] = df["species_id"].map(order)
    df = df.sort_values("_order").drop(columns="_order")

    for col in ("rep_exon_changed_before_pct", "rep_exon_changed_after_pct"):
        if ((df[col] < 0) | (df[col] > 100)).any():
            raise SystemExit(f"{path}: {col} must be within 0-100")
    if (df["rep_exon_changed"] > df["total_before_genes"]).any():
        raise SystemExit(f"{path}: rep_exon_changed exceeds total_before_genes")
    if (df["rep_exon_changed"] > df["total_after_genes"]).any():
        raise SystemExit(f"{path}: rep_exon_changed exceeds total_after_genes")

    return df


def _format_count(value: float) -> str:
    return f"{int(round(value)):,}"


def _count_formatter(value, _pos):
    value = abs(value)
    if value >= 1000:
        return f"{value / 1000:.0f}K"
    return f"{value:.0f}"


def _style_axis(ax, show_y=False):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6 if show_y else 0.0)
    ax.spines["bottom"].set_linewidth(0.6)
    if not show_y:
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0, labelleft=False)
    ax.tick_params(axis="both", labelsize=6.8, width=0.6, length=2.5)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.55)
    ax.set_axisbelow(True)


def _panel_header(ax, label, title):
    ax.text(
        0.0,
        1.055,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10.5,
        fontweight="bold",
    )
    ax.text(
        0.075,
        1.055,
        title,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.8,
        fontweight="bold",
    )


def _annotate_barh(ax, x, y, text, color="#222222", pad_frac=0.012, ha=None):
    xmin, xmax = ax.get_xlim()
    span = xmax - xmin
    pad = span * pad_frac
    if ha is None:
        ha = "left" if x >= 0 else "right"
    xpos = x + pad if x >= 0 else x - pad
    ax.text(
        xpos,
        y,
        text,
        ha=ha,
        va="center",
        fontsize=5.2,
        color=color,
        clip_on=False,
    )


def _annotate_percent_bar(ax, pct, y, text):
    if pct >= 52:
        ax.text(
            pct - 1.4,
            y,
            text,
            ha="right",
            va="center",
            fontsize=5.2,
            color="white",
            clip_on=True,
        )
    else:
        _annotate_barh(ax, pct, y, text, pad_frac=0.008)


def _set_shared_y(ax, y, labels, show_labels):
    ax.set_yticks(y)
    if show_labels:
        ax.set_yticklabels(labels, fontsize=6.8)
    else:
        ax.set_yticklabels([])
    ax.set_ylim(len(y) - 0.45, -1.25)


def _compact_legend(ax, ncol, fontsize=5.6):
    ax.legend(
        frameon=False,
        fontsize=fontsize,
        loc="upper left",
        bbox_to_anchor=(0.0, 0.985),
        ncol=ncol,
        borderaxespad=0.0,
        handlelength=0.9,
        handletextpad=0.35,
        columnspacing=0.7,
        labelspacing=0.35,
    )


def plot_locus_gain_loss(ax, df, y, labels, panel_label="A", show_labels=False):
    offset = 0.16
    height = 0.27
    deleted = df["deleted_loci_no_overlap"].to_numpy()
    new = df["new_loci_no_overlap"].to_numpy()
    max_value = max(float(deleted.max()), float(new.max()), 1.0)

    ax.barh(y - offset, deleted, height=height, color=COLORS["deleted"], label="Deleted loci")
    ax.barh(y + offset, new, height=height, color=COLORS["new"], label="New loci")
    ax.set_xlim(0, max_value * 1.2)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(_count_formatter))
    ax.set_xlabel("No-overlap loci", fontsize=7.8)
    _set_shared_y(ax, y, labels, show_labels=show_labels)
    _style_axis(ax, show_y=show_labels)
    _panel_header(ax, panel_label, "Locus gain/loss")
    _compact_legend(ax, ncol=2)

    for yi, value in zip(y - offset, deleted):
        if value:
            _annotate_barh(ax, value, yi, _format_count(value), pad_frac=0.01)
    for yi, value in zip(y + offset, new):
        if value:
            _annotate_barh(ax, value, yi, _format_count(value), pad_frac=0.01)


def plot_split_merge(ax, df, y, labels, panel_label="B", show_labels=False):
    offset = 0.16
    height = 0.27
    split = df["split_events"].to_numpy()
    merge = df["merge_events"].to_numpy()
    max_value = max(float(split.max()), float(merge.max()), 1.0)
    ax.set_xlim(0, max_value * 1.2)

    ax.barh(y - offset, split, height=height, color=COLORS["split"], label="Split")
    ax.barh(y + offset, merge, height=height, color=COLORS["merge"], label="Merge")
    ax.set_xlabel("Events", fontsize=7.8)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(_count_formatter))
    _set_shared_y(ax, y, labels, show_labels=show_labels)
    _style_axis(ax, show_y=show_labels)
    _panel_header(ax, panel_label, "Split and merge events")
    _compact_legend(ax, ncol=2)

    zero_pad = max_value * 0.012
    for yi, value in zip(y - offset, split):
        _annotate_barh(ax, value if value else zero_pad, yi, _format_count(value), pad_frac=0.01)
    for yi, value in zip(y + offset, merge):
        _annotate_barh(ax, value if value else zero_pad, yi, _format_count(value), pad_frac=0.01)


def plot_representative_transcript_changes(ax, df, y, labels, panel_label="C", show_labels=False):
    offset = 0.16
    height = 0.27
    changed = df["rep_exon_changed"].to_numpy()
    before_total = df["total_before_genes"].to_numpy()
    after_total = df["total_after_genes"].to_numpy()
    before_pct = df["rep_exon_changed_before_pct"].to_numpy()
    after_pct = df["rep_exon_changed_after_pct"].to_numpy()

    ax.barh(
        y - offset,
        before_pct,
        height=height,
        color=COLORS["before_ref"],
        label="Before ref.",
    )
    ax.barh(
        y + offset,
        after_pct,
        height=height,
        color=COLORS["after_ref"],
        label="After ref.",
    )
    ax.set_xlim(0, 100)
    ax.set_xlabel("Genes with representative exon changes (%)", fontsize=7.8)
    ax.xaxis.set_major_formatter(mticker.PercentFormatter(xmax=100, decimals=0))
    _set_shared_y(ax, y, labels, show_labels=show_labels)
    _style_axis(ax, show_y=show_labels)
    _panel_header(ax, panel_label, "Representative exon changes")
    _compact_legend(ax, ncol=2)

    for yi, value, count, denom in zip(y - offset, before_pct, changed, before_total):
        _annotate_percent_bar(
            ax,
            value,
            yi,
            f"{value:.1f}% ({_format_count(count)}/{_format_count(denom)})",
        )
    for yi, value, count, denom in zip(y + offset, after_pct, changed, after_total):
        _annotate_percent_bar(
            ax,
            value,
            yi,
            f"{value:.1f}% ({_format_count(count)}/{_format_count(denom)})",
        )


def plot_figure(df: pd.DataFrame, output_prefix: Path, dpi: int) -> None:
    labels = df["Species"].tolist()
    y = np.arange(len(df))

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.6,
    })

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(10.1, 4.35),
        gridspec_kw={"width_ratios": [1.34, 1.02, 1.35], "wspace": 0.22},
    )

    plot_locus_gain_loss(axes[0], df, y, labels, panel_label="A", show_labels=True)
    plot_split_merge(axes[1], df, y, labels, panel_label="B")
    plot_representative_transcript_changes(axes[2], df, y, labels, panel_label="C")

    fig.subplots_adjust(left=0.095, right=0.99, top=0.84, bottom=0.16)

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = output_prefix.with_suffix(".pdf")
    png_path = output_prefix.with_suffix(".png")
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")


def main() -> None:
    args = parse_args()
    df = load_metrics(Path(args.input))
    plot_figure(df, Path(args.output_prefix), args.dpi)


if __name__ == "__main__":
    main()
