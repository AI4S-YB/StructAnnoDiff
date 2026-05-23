"""Generate the two primary figures for the curated species set."""

from __future__ import annotations

import os

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
import seaborn as sns


FIGURES_DIR.mkdir(exist_ok=True)
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.15)

SPECIES_IDS = [sp.id for sp in SPECIES]
SPECIES_LABELS = [sp.label for sp in SPECIES]

BEFORE_COLOR = "#C44E52"
AFTER_COLOR = "#4C72B0"


def _ordered_rows(df: pd.DataFrame, species_col: str = "species") -> pd.DataFrame:
    order = {species_id: i for i, species_id in enumerate(SPECIES_IDS)}
    out = df[df[species_col].isin(order)].copy()
    out["_order"] = out[species_col].map(order)
    return out.sort_values("_order").drop(columns="_order")


def _metric_pair(row: pd.Series, metric: str) -> tuple[float, float]:
    before = row.get(f"{metric}_before", 0)
    after = row.get(f"{metric}_after", 0)
    return (
        float(before) if pd.notna(before) else 0.0,
        float(after) if pd.notna(after) else 0.0,
    )


def _grouped_before_after(ax, labels, before, after, ylabel, title, value_fmt=None):
    x = np.arange(len(labels))
    width = 0.36
    ax.bar(x - width / 2, before, width, color=BEFORE_COLOR, label="Before", alpha=0.9)
    ax.bar(x + width / 2, after, width, color=AFTER_COLOR, label="After", alpha=0.9)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=28, ha="right")
    if value_fmt:
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(value_fmt))


def plot_quantity_changes() -> None:
    stats = _ordered_rows(pd.read_csv(RESULTS_DIR / "summary_stats.csv"))

    metrics = [
        ("Number of gene", "Gene count", "A. Total genes", lambda v, _: f"{v/1000:.0f}K"),
        ("Number of single exon gene", "Single-exon genes", "B. Single-exon genes", lambda v, _: f"{v/1000:.0f}K"),
        ("mean exons per mrna", "Mean exons per mRNA", "C. Exons per transcript", None),
        ("mean cds length (bp)", "Mean CDS length (bp)", "D. Mean CDS length", lambda v, _: f"{v/1000:.1f}K"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    labels = [sp.label for sp in SPECIES]

    for ax, (metric, ylabel, title, formatter) in zip(axes.flatten(), metrics):
        before = []
        after = []
        for _, row in stats.iterrows():
            b, a = _metric_pair(row, metric)
            before.append(b)
            after.append(a)
        _grouped_before_after(ax, labels, before, after, ylabel, title, formatter)

    handles, legend_labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="upper center", ncol=2, frameon=False)
    fig.suptitle("Annotation quantity changes after manual curation", y=0.99, fontweight="bold")
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    out = FIGURES_DIR / "figure1_quantity_changes.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.name}")


def plot_syntenic_structure_changes() -> None:
    summary = pd.read_csv(RESULTS_DIR / "locus_comparison_summary.csv")
    multilabel = pd.read_csv(RESULTS_DIR / "locus_comparison_multilabel.csv")
    label_to_order = {sp.label: i for i, sp in enumerate(SPECIES)}
    summary = summary[summary["Species"].isin(label_to_order)].copy()
    summary["_order"] = summary["Species"].map(label_to_order)
    summary = summary.sort_values("_order").drop(columns="_order")
    multilabel = multilabel[multilabel["Species"].isin(label_to_order)].copy()
    multilabel["_order"] = multilabel["Species"].map(label_to_order)
    multilabel = multilabel.sort_values("_order").drop(columns="_order")

    plot_df = multilabel.merge(
        summary[["Species", "Before"]],
        on="Species",
        how="left",
    )

    categories = [
        ("Exact", "Exact"),
        ("Any_gene_boundary_changed", "Gene boundary"),
        ("Any_UTR_added", "UTR added"),
        ("Any_UTR_lost", "UTR lost"),
        ("Any_UTR_exon_gained", "UTR exon gained"),
        ("Any_UTR_exon_removed", "UTR exon removed"),
        ("Any_UTR_refined", "UTR refined"),
        ("Any_coding_exon_gain", "Coding exon gain"),
        ("Any_coding_exon_loss", "Coding exon loss"),
        ("Any_exon_boundary_refined", "Exon boundary"),
        ("Any_CDS_change", "CDS changed"),
        ("Any_CDS_boundary_refined", "CDS boundary"),
        ("Any_isoform_change", "Isoform changed"),
    ]
    categories = [
        (col, label)
        for col, label in categories
        if col == "Exact" or (col in plot_df and plot_df[col].fillna(0).sum() > 0)
    ]

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(13, 9),
        gridspec_kw={"height_ratios": [3.2, 1.15]},
    )
    labels = plot_df["Species"].tolist()

    ax = axes[0]
    heatmap = []
    for col, _ in categories:
        values = plot_df[col] if col in plot_df else pd.Series(0, index=plot_df.index)
        heatmap.append((values / plot_df["Syntenic"] * 100).to_numpy())
    heatmap = np.array(heatmap)
    vmax = min(100, np.nanpercentile(heatmap, 96))
    vmax = vmax if np.isfinite(vmax) and vmax > 0 else 1
    im = ax.imshow(heatmap, aspect="auto", cmap="YlGnBu", vmin=0, vmax=vmax)
    ax.set_title("A. Non-exclusive structural attributes among confirmed one-to-one genes", loc="left", fontweight="bold")
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=28, ha="right")
    ax.set_yticks(np.arange(len(categories)))
    ax.set_yticklabels([label for _, label in categories])
    for i in range(heatmap.shape[0]):
        for j in range(heatmap.shape[1]):
            value = heatmap[i, j]
            if value >= 0.5:
                color = "white" if value >= vmax * 0.58 else "#2b2b2b"
                ax.text(j, i, f"{value:.1f}", ha="center", va="center", fontsize=7, color=color)
    fig.colorbar(im, ax=ax, label="% of one-to-one genes", fraction=0.035, pad=0.015)

    ax = axes[1]
    x = np.arange(len(plot_df))
    syntenic_rate = plot_df["Syntenic"] / plot_df["Before"] * 100
    ax.bar(x, syntenic_rate, color="#4C72B0", alpha=0.85)
    ax.set_ylabel("1:1 genes (%)")
    ax.set_title("B. One-to-one genes as a fraction of before genes", loc="left", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=28, ha="right")
    ax.set_ylim(0, 100)
    for i, value in enumerate(syntenic_rate):
        ax.text(i, value + 1.2, f"{value:.1f}", ha="center", va="bottom", fontsize=8)

    fig.suptitle("Structure changes in confidently matched gene pairs", y=0.99, fontweight="bold")
    plt.tight_layout(rect=(0, 0, 1, 0.96))
    out = FIGURES_DIR / "figure2_syntenic_structure_changes.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.name}")


def main() -> None:
    plot_quantity_changes()
    plot_syntenic_structure_changes()


if __name__ == "__main__":
    main()
