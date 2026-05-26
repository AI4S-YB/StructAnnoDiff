#!/usr/bin/env python3
"""Build A/B/C/D single-species tables and a compact summary figure."""

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create A/B/C/D tables and a figure for one species."
    )
    parser.add_argument(
        "--species",
        default="Pineapple",
        help="Species ID or display label. Defaults to Pineapple.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(RESULTS_DIR / "single_species"),
        help="Directory for generated CSV/Markdown tables.",
    )
    parser.add_argument(
        "--figure-dir",
        default=str(FIGURES_DIR),
        help="Directory for generated PNG figure.",
    )
    return parser.parse_args()


def _species_from_arg(value: str):
    for sp in SPECIES:
        if value in {sp.id, sp.label, sp.short_label}:
            return sp
    known = ", ".join(sp.id for sp in SPECIES)
    raise SystemExit(f"Unknown species {value!r}. Known species IDs: {known}")


def _as_float(value) -> float:
    if pd.isna(value):
        return 0.0
    return float(value)


def _fmt_num(value) -> str:
    if pd.isna(value):
        return ""
    value = float(value)
    if abs(value - round(value)) < 1e-9:
        return f"{int(round(value)):,}"
    return f"{value:,.2f}"


def _pct(after: float, before: float) -> float:
    return (after - before) / before * 100 if before else 0.0


def build_figure_a_table(sp_id: str) -> pd.DataFrame:
    stats = pd.read_csv(RESULTS_DIR / "summary_stats.csv")
    row = stats.loc[stats["species"] == sp_id]
    if row.empty:
        raise SystemExit(f"{sp_id}: missing in results/summary_stats.csv")
    row = row.iloc[0]

    metrics = [
        ("Number of gene", "Genes", "count"),
        ("Number of mrna", "mRNAs", "count"),
        ("Number of exon", "Exons", "count"),
        ("Number of cds", "CDS features", "count"),
        ("Number of five_prime_utr", "5' UTR features", "count"),
        ("Number of three_prime_utr", "3' UTR features", "count"),
        ("Number of single exon gene", "Single-exon genes", "count"),
        ("mean exons per mrna", "Mean exons per mRNA", "mean"),
        ("mean cds length (bp)", "Mean CDS length (bp)", "mean"),
        ("mean gene length (bp)", "Mean gene length (bp)", "mean"),
        ("median gene length (bp)", "Median gene length (bp)", "median"),
    ]

    rows = []
    for metric, label, metric_type in metrics:
        before = _as_float(row.get(f"{metric}_before", 0))
        after = _as_float(row.get(f"{metric}_after", 0))
        rows.append({
            "figure": "A",
            "layer": "Global annotation quantity",
            "metric": label,
            "metric_type": metric_type,
            "before": before,
            "after": after,
            "delta": after - before,
            "pct_change": _pct(after, before),
        })
    return pd.DataFrame(rows)


def build_figure_b_table(sp_id: str) -> pd.DataFrame:
    path = RESULTS_DIR / "locus" / f"{sp_id}_change_summary.csv"
    row = pd.read_csv(path).iloc[0]
    before_total = _as_float(row["total_before_genes"])
    after_total = _as_float(row["total_after_genes"])
    specs = [
        ("Syntenic 1:1", "syntenic_total", "syntenic_total"),
        ("Split", "before_in_splits", "after_in_splits"),
        ("Merge", "before_in_merges", "after_in_merges"),
        ("Complex", "before_in_complex", "after_in_complex"),
        (
            "Unresolved weak-overlap",
            "unresolved_overlap_before_genes",
            "unresolved_overlap_after_genes",
        ),
        ("Strict unmatched", "deleted_genes", "novel_genes"),
    ]
    rows = []
    for label, before_col, after_col in specs:
        before = _as_float(row.get(before_col, 0))
        after = _as_float(row.get(after_col, 0))
        rows.append({
            "figure": "B",
            "layer": "Locus fate",
            "category": label,
            "before_genes": before,
            "before_pct": before / before_total * 100 if before_total else 0.0,
            "after_genes": after,
            "after_pct": after / after_total * 100 if after_total else 0.0,
        })
    return pd.DataFrame(rows)


def build_figure_c_table(sp_label: str) -> pd.DataFrame:
    multilabel = pd.read_csv(RESULTS_DIR / "locus_comparison_multilabel.csv")
    row = multilabel.loc[multilabel["Species"] == sp_label]
    if row.empty:
        raise SystemExit(f"{sp_label}: missing in locus_comparison_multilabel.csv")
    row = row.iloc[0]
    syntenic = _as_float(row["Syntenic"])
    specs = [
        ("Exact 1:1", "Exact"),
        ("Gene boundary changed", "Any_gene_boundary_changed"),
        ("UTR added", "Any_UTR_added"),
        ("UTR lost", "Any_UTR_lost"),
        ("UTR exon gained", "Any_UTR_exon_gained"),
        ("UTR exon removed", "Any_UTR_exon_removed"),
        ("UTR refined", "Any_UTR_refined"),
        ("Coding exon gain", "Any_coding_exon_gain"),
        ("Coding exon loss", "Any_coding_exon_loss"),
        ("Exon boundary refined", "Any_exon_boundary_refined"),
        ("CDS changed", "Any_CDS_change"),
        ("CDS boundary refined", "Any_CDS_boundary_refined"),
        ("Isoform changed", "Any_isoform_change"),
    ]
    rows = []
    for label, col in specs:
        count = _as_float(row.get(col, 0))
        rows.append({
            "figure": "C",
            "layer": "Confirmed 1:1 structural attributes",
            "attribute": label,
            "count": count,
            "pct_of_syntenic": count / syntenic * 100 if syntenic else 0.0,
        })
    return pd.DataFrame(rows)


def _metric_delta_summary(df: pd.DataFrame, label: str, before_col: str, after_col: str) -> dict:
    before = pd.to_numeric(df[before_col], errors="coerce").fillna(0)
    after = pd.to_numeric(df[after_col], errors="coerce").fillna(0)
    delta = after - before
    changed = delta != 0
    valid_rel = before > 0
    rel = pd.Series(np.nan, index=df.index, dtype=float)
    rel.loc[valid_rel] = delta.loc[valid_rel] / before.loc[valid_rel] * 100
    total = len(df)
    return {
        "figure": "D",
        "layer": "Confirmed 1:1 change magnitude",
        "metric": label,
        "n_pairs": total,
        "mean_before": before.mean(),
        "mean_after": after.mean(),
        "mean_delta": delta.mean(),
        "median_delta": delta.median(),
        "p05_delta": delta.quantile(0.05),
        "p95_delta": delta.quantile(0.95),
        "mean_pct_delta": rel.mean(skipna=True),
        "decreased": int((delta < 0).sum()),
        "unchanged": int((delta == 0).sum()),
        "increased": int((delta > 0).sum()),
        "changed": int(changed.sum()),
        "changed_pct": changed.sum() / total * 100 if total else 0.0,
    }


def build_figure_d_tables(sp_id: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = RESULTS_DIR / "locus" / f"{sp_id}_change_log.csv"
    log = pd.read_csv(path)
    syntenic = log.loc[log["match_type"] == "syntenic"].copy()
    metrics = [
        ("Gene span length", "before_gene_length", "after_gene_length"),
        ("Model span length", "before_length", "after_length"),
        ("CDS length", "before_cds", "after_cds"),
        ("Exon count", "before_exons", "after_exons"),
        ("mRNA count", "before_mrnas", "after_mrnas"),
    ]
    summary = pd.DataFrame([
        _metric_delta_summary(syntenic, label, before_col, after_col)
        for label, before_col, after_col in metrics
    ])

    pair_rows = syntenic[[
        "before_gene",
        "after_gene",
        "seqid",
        "strand",
        "change_subtype",
        "before_gene_length",
        "after_gene_length",
        "before_length",
        "after_length",
        "before_cds",
        "after_cds",
        "before_exons",
        "after_exons",
        "before_mrnas",
        "after_mrnas",
    ]].copy()
    pair_rows["gene_length_delta"] = pair_rows["after_gene_length"] - pair_rows["before_gene_length"]
    pair_rows["model_length_delta"] = pair_rows["after_length"] - pair_rows["before_length"]
    pair_rows["cds_delta"] = pair_rows["after_cds"] - pair_rows["before_cds"]
    pair_rows["exon_count_delta"] = pair_rows["after_exons"] - pair_rows["before_exons"]
    pair_rows["mrna_count_delta"] = pair_rows["after_mrnas"] - pair_rows["before_mrnas"]
    return summary, pair_rows


def _write_markdown(path: Path, title: str, tables: list[tuple[str, pd.DataFrame]]) -> None:
    lines = [f"# {title}", ""]
    for heading, df in tables:
        lines.extend([f"## {heading}", ""])
        out = df.copy()
        for col in out.columns:
            if pd.api.types.is_numeric_dtype(out[col]):
                out[col] = out[col].map(_fmt_num)
        lines.append(out.to_csv(index=False))
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def _bar_label(ax, bars, fmt="{:.0f}", pad=3):
    for bar in bars:
        value = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + pad,
            fmt.format(value),
            ha="center",
            va="bottom",
            fontsize=8,
        )


def plot_abcd(sp_label: str, table_a: pd.DataFrame, table_b: pd.DataFrame,
              table_c: pd.DataFrame, table_d: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(17, 12))
    fig.suptitle(f"{sp_label}: annotation curation impact by analysis layer", y=0.995,
                 fontsize=16, fontweight="bold")

    # A. Global quantities
    ax = axes[0, 0]
    a_counts = table_a.loc[table_a["metric_type"] == "count"].copy()
    labels = a_counts["metric"].tolist()
    x = np.arange(len(labels))
    width = 0.36
    before = a_counts["before"].to_numpy()
    after = a_counts["after"].to_numpy()
    ax.bar(x - width / 2, before, width, label="Before", color="#C44E52", alpha=0.88)
    ax.bar(x + width / 2, after, width, label="After", color="#4C72B0", alpha=0.88)
    ax.set_title("A. Global annotation quantities", loc="left", fontweight="bold")
    ax.set_ylabel("Count")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=32, ha="right")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v/1000:.0f}K" if v >= 1000 else f"{v:.0f}"))
    ax.legend(frameon=False, ncol=2)

    # B. Locus fate accounting
    ax = axes[0, 1]
    fate_colors = {
        "Syntenic 1:1": "#4C72B0",
        "Split": "#F28E2B",
        "Merge": "#8A60B0",
        "Complex": "#8C8C8C",
        "Unresolved weak-overlap": "#EDC948",
        "Strict unmatched": "#E15759",
    }
    before_pct = table_b["before_pct"].to_numpy()
    after_pct = table_b["after_pct"].to_numpy()
    left_before = 0.0
    left_after = 0.0
    for _, row in table_b.iterrows():
        color = fate_colors[row["category"]]
        ax.barh(1, row["before_pct"], left=left_before, color=color, edgecolor="white", height=0.36)
        ax.barh(0, row["after_pct"], left=left_after, color=color, edgecolor="white", height=0.36,
                label=row["category"])
        if row["before_pct"] >= 3:
            ax.text(left_before + row["before_pct"] / 2, 1, f"{row['before_pct']:.1f}%",
                    ha="center", va="center", fontsize=8, color="white" if row["category"] == "Syntenic 1:1" else "#222")
        if row["after_pct"] >= 3:
            ax.text(left_after + row["after_pct"] / 2, 0, f"{row['after_pct']:.1f}%",
                    ha="center", va="center", fontsize=8, color="white" if row["category"] == "Syntenic 1:1" else "#222")
        left_before += row["before_pct"]
        left_after += row["after_pct"]
    ax.set_title("B. Locus fate accounting", loc="left", fontweight="bold")
    ax.set_xlim(0, 100)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["After genes", "Before genes"])
    ax.set_xlabel("% of genes")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.36), ncol=2, frameon=False, fontsize=8)

    # C. 1:1 structural attributes
    ax = axes[1, 0]
    c_plot = table_c.loc[table_c["count"] > 0].copy()
    c_plot = c_plot.sort_values("pct_of_syntenic")
    colors = ["#59A14F" if attr == "Exact 1:1" else "#76B7B2" for attr in c_plot["attribute"]]
    ax.barh(c_plot["attribute"], c_plot["pct_of_syntenic"], color=colors, alpha=0.9)
    ax.set_title("C. Structural attributes in confirmed 1:1 genes", loc="left", fontweight="bold")
    ax.set_xlabel("% of syntenic 1:1 genes")
    ax.xaxis.set_major_formatter(mticker.PercentFormatter(xmax=100))
    for y, value in enumerate(c_plot["pct_of_syntenic"]):
        ax.text(value + 0.7, y, f"{value:.1f}%", va="center", fontsize=8)
    ax.set_xlim(0, max(5, min(100, c_plot["pct_of_syntenic"].max() * 1.18)))

    # D. Paired magnitude directions
    ax = axes[1, 1]
    d_plot = table_d.copy()
    d_plot["decreased_pct"] = d_plot["decreased"] / d_plot["n_pairs"] * 100
    d_plot["unchanged_pct"] = d_plot["unchanged"] / d_plot["n_pairs"] * 100
    d_plot["increased_pct"] = d_plot["increased"] / d_plot["n_pairs"] * 100
    y = np.arange(len(d_plot))
    left = np.zeros(len(d_plot))
    for col, label, color in [
        ("decreased_pct", "Decreased", "#E15759"),
        ("unchanged_pct", "Unchanged", "#BAB0AC"),
        ("increased_pct", "Increased", "#59A14F"),
    ]:
        vals = d_plot[col].to_numpy()
        ax.barh(y, vals, left=left, color=color, edgecolor="white", height=0.62, label=label)
        for i, val in enumerate(vals):
            if val >= 8:
                ax.text(left[i] + val / 2, i, f"{val:.0f}%", ha="center", va="center", fontsize=8)
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels(d_plot["metric"])
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of syntenic 1:1 gene pairs")
    ax.set_title("D. Direction of paired changes", loc="left", fontweight="bold")
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.22))

    for ax in axes.flatten():
        ax.grid(axis="y", alpha=0.18)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout(rect=(0, 0.025, 1, 0.965))
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    sp = _species_from_arg(args.species)
    output_dir = Path(args.output_dir)
    figure_dir = Path(args.figure_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    table_a = build_figure_a_table(sp.id)
    table_b = build_figure_b_table(sp.id)
    table_c = build_figure_c_table(sp.label)
    table_d, pair_deltas = build_figure_d_tables(sp.id)

    prefix = sp.id
    paths = {
        "A": output_dir / f"{prefix}_figureA_quantity_table.csv",
        "B": output_dir / f"{prefix}_figureB_locus_fate_table.csv",
        "C": output_dir / f"{prefix}_figureC_syntenic_structure_table.csv",
        "D": output_dir / f"{prefix}_figureD_pair_delta_summary.csv",
        "pairs": output_dir / f"{prefix}_syntenic_pair_deltas.csv",
        "md": output_dir / f"{prefix}_ABCD_tables.md",
        "fig": figure_dir / f"{prefix}_ABCD_single_species.png",
    }
    table_a.to_csv(paths["A"], index=False)
    table_b.to_csv(paths["B"], index=False)
    table_c.to_csv(paths["C"], index=False)
    table_d.to_csv(paths["D"], index=False)
    pair_deltas.to_csv(paths["pairs"], index=False)
    _write_markdown(
        paths["md"],
        f"{sp.label} A/B/C/D analysis tables",
        [
            ("Figure A table: global quantities", table_a),
            ("Figure B table: locus fate", table_b),
            ("Figure C table: 1:1 structural attributes", table_c),
            ("Figure D table: paired change magnitude", table_d),
        ],
    )
    plot_abcd(sp.label, table_a, table_b, table_c, table_d, paths["fig"])

    print(f"Saved: {paths['A']}")
    print(f"Saved: {paths['B']}")
    print(f"Saved: {paths['C']}")
    print(f"Saved: {paths['D']}")
    print(f"Saved: {paths['pairs']}")
    print(f"Saved: {paths['md']}")
    print(f"Saved: {paths['fig']}")


if __name__ == "__main__":
    main()
