#!/usr/bin/env python3
"""Validate cross-table consistency for annotation comparison outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from analysis_config import ANALYSIS, RESULTS_DIR, SPECIES, SPECIES_LABELS

EXCLUSIVE_LOCUS_COLS = [
    "Exact",
    "Boundary_refined_only",
    "UTR_added",
    "UTR_lost",
    "UTR_exon_gained",
    "UTR_exon_removed",
    "UTR_refined",
    "Coding_exon_gain",
    "Coding_exon_loss",
    "Exon_boundary_refined",
    "CDS_change_only",
    "CDS_boundary_refined",
    "Isoform_change",
]

MULTILABEL_DIRECT_MAP = {
    "Exact": "one_to_one_exact",
    "Any_gene_boundary_changed": "one_to_one_gene_boundary_changed",
    "Any_UTR_added": "one_to_one_utr_added",
    "Any_UTR_lost": "one_to_one_utr_lost",
    "Any_UTR_exon_gained": "one_to_one_utr_exon_added",
    "Any_UTR_exon_removed": "one_to_one_utr_exon_removed",
    "Any_UTR_refined": "one_to_one_utr_refined",
    "Any_coding_exon_gain": "one_to_one_coding_exon_gain",
    "Any_coding_exon_loss": "one_to_one_coding_exon_loss",
    "Any_exon_boundary_refined": "one_to_one_exon_boundary_refined",
    "Any_CDS_change": "one_to_one_cds_change",
    "Any_CDS_boundary_refined": "one_to_one_cds_boundary_refined",
    "Any_isoform_change": "one_to_one_isoform_change",
}


def _read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def _as_int(value) -> int:
    return int(float(value)) if pd.notna(value) else 0


def _locus_with_species_ids(locus_df: pd.DataFrame) -> pd.DataFrame:
    if locus_df.empty:
        return locus_df
    label_to_id = {sp.label: sp.id for sp in SPECIES}
    locus_df = locus_df.copy()
    locus_df["species"] = locus_df["Species"].map(label_to_id)
    return locus_df


def build_validation_table(results_dir: Path = RESULTS_DIR) -> pd.DataFrame:
    stats = _read_csv(results_dir / "summary_stats.csv")
    compare = _read_csv(results_dir / "comparison_matrix.csv")
    compare_all = _read_csv(results_dir / "comparison_matrix_all_gene_types.csv")
    locus = _locus_with_species_ids(_read_csv(results_dir / "locus_comparison_summary.csv"))

    rows = []
    for sp in SPECIES:
        row = {
            "species": sp.id,
            "label": SPECIES_LABELS[sp.id],
        }

        stats_row = stats[stats["species"] == sp.id] if not stats.empty else pd.DataFrame()
        if not stats_row.empty:
            stats_row = stats_row.iloc[0]
            row["stats_coding_section_before"] = stats_row.get("coding_section_before")
            row["stats_coding_section_after"] = stats_row.get("coding_section_after")
            row["stats_gene_before"] = stats_row.get("Number of gene_before")
            row["stats_gene_after"] = stats_row.get("Number of gene_after")

        compare_row = compare[compare["species"] == sp.id] if not compare.empty else pd.DataFrame()
        if not compare_row.empty:
            compare_row = compare_row.iloc[0]
            row["compare_coding_before"] = compare_row.get("total_before")
            row["compare_coding_after"] = compare_row.get("total_after")

        compare_all_row = compare_all[compare_all["species"] == sp.id] if not compare_all.empty else pd.DataFrame()
        if not compare_all_row.empty:
            compare_all_row = compare_all_row.iloc[0]
            row["compare_all_before"] = compare_all_row.get("total_before")
            row["compare_all_after"] = compare_all_row.get("total_after")

        locus_row = locus[locus["species"] == sp.id] if not locus.empty else pd.DataFrame()
        if not locus_row.empty:
            locus_row = locus_row.iloc[0]
            row["locus_gene_before"] = locus_row.get("Before")
            row["locus_gene_after"] = locus_row.get("After")

        for lhs, rhs, name in [
            ("stats_gene_before", "compare_coding_before", "stats_vs_compare_coding_before_delta"),
            ("stats_gene_after", "compare_coding_after", "stats_vs_compare_coding_after_delta"),
            ("stats_gene_before", "locus_gene_before", "stats_vs_locus_before_delta"),
            ("stats_gene_after", "locus_gene_after", "stats_vs_locus_after_delta"),
        ]:
            if lhs in row and rhs in row and pd.notna(row[lhs]) and pd.notna(row[rhs]):
                row[name] = row[lhs] - row[rhs]

        rows.append(row)

    return pd.DataFrame(rows)


def run_integrity_checks(results_dir: Path = RESULTS_DIR) -> list[str]:
    """Return result-level consistency issues that should fail validation."""
    issues = []
    expected_labels = [sp.label for sp in SPECIES]

    for table_name in [
        "summary_stats.csv",
        "comparison_matrix.csv",
        "comparison_matrix_all_gene_types.csv",
        "accuracy_metrics.csv",
    ]:
        table = _read_csv(results_dir / table_name)
        if table.empty or "species" not in table:
            issues.append(f"{table_name}: missing or has no species column")
            continue
        observed = set(table["species"].dropna())
        expected = {sp.id for sp in SPECIES}
        extra = sorted(observed - expected)
        missing = sorted(expected - observed)
        if extra or missing:
            issues.append(f"{table_name}: extra={extra}, missing={missing}")

    summary = _read_csv(results_dir / "locus_comparison_summary.csv")
    if summary.empty:
        issues.append("locus_comparison_summary.csv: missing or empty")
    else:
        labels = summary["Species"].tolist()
        if labels != expected_labels:
            issues.append(
                "locus_comparison_summary.csv: species order/content mismatch "
                f"{labels} != {expected_labels}"
            )
        for _, row in summary.iterrows():
            subtotal = sum(_as_int(row[col]) for col in EXCLUSIVE_LOCUS_COLS)
            syntenic = _as_int(row["Syntenic"])
            if subtotal != syntenic:
                issues.append(
                    f"{row['Species']}: exclusive locus categories {subtotal} "
                    f"!= Syntenic {syntenic}"
                )

    multilabel = _read_csv(results_dir / "locus_comparison_multilabel.csv")
    if multilabel.empty:
        issues.append("locus_comparison_multilabel.csv: missing or empty")
    else:
        labels = multilabel["Species"].tolist()
        if labels != expected_labels:
            issues.append(
                "locus_comparison_multilabel.csv: species order/content mismatch "
                f"{labels} != {expected_labels}"
            )

    label_to_id = {sp.label: sp.id for sp in SPECIES}
    for sp in SPECIES:
        locus_summary_path = results_dir / "locus" / f"{sp.id}_change_summary.csv"
        if not locus_summary_path.exists():
            issues.append(f"{locus_summary_path}: missing")
            continue

        row = pd.read_csv(locus_summary_path).iloc[0]
        syntenic = _as_int(row["syntenic_total"])
        subtype_sum = sum(
            _as_int(value)
            for col, value in row.items()
            if col.startswith("syntenic_") and col != "syntenic_total"
        )
        if subtype_sum != syntenic:
            issues.append(f"{sp.id}: syntenic subtypes {subtype_sum} != {syntenic}")

        before_accounted = (
            syntenic
            + _as_int(row["before_in_splits"])
            + _as_int(row["before_in_merges"])
            + _as_int(row["before_in_complex"])
            + _as_int(row.get("unresolved_overlap_before_genes", 0))
            + _as_int(row["deleted_genes"])
        )
        after_accounted = (
            syntenic
            + _as_int(row["after_in_splits"])
            + _as_int(row["after_in_merges"])
            + _as_int(row["after_in_complex"])
            + _as_int(row.get("unresolved_overlap_after_genes", 0))
            + _as_int(row["novel_genes"])
        )
        if before_accounted != _as_int(row["total_before_genes"]):
            issues.append(f"{sp.id}: before locus accounting does not close")
        if after_accounted != _as_int(row["total_after_genes"]):
            issues.append(f"{sp.id}: after locus accounting does not close")

    if not multilabel.empty:
        for _, ml_row in multilabel.iterrows():
            species_id = label_to_id.get(ml_row["Species"])
            if species_id is None:
                continue
            locus_summary_path = results_dir / "locus" / f"{species_id}_change_summary.csv"
            if not locus_summary_path.exists():
                continue
            summary_row = pd.read_csv(locus_summary_path).iloc[0]
            syntenic = _as_int(ml_row["Syntenic"])
            for output_col, source_col in MULTILABEL_DIRECT_MAP.items():
                output_value = _as_int(ml_row[output_col])
                source_value = _as_int(summary_row[source_col])
                if output_value != source_value:
                    issues.append(
                        f"{ml_row['Species']}: {output_col}={output_value} "
                        f"!= {source_col}={source_value}"
                    )
                if output_value > syntenic:
                    issues.append(
                        f"{ml_row['Species']}: {output_col}={output_value} "
                        f"> Syntenic={syntenic}"
                    )

    return issues


def write_markdown_report(df: pd.DataFrame, path: Path, integrity_issues=None) -> None:
    delta_cols = [c for c in df.columns if c.endswith("_delta")]
    flagged = df[df[delta_cols].fillna(0).abs().sum(axis=1) > 0] if delta_cols else pd.DataFrame()
    integrity_issues = integrity_issues or []

    lines = [
        "# Validation Report",
        "",
        "This report compares mRNA/transcript-scope AGAT stats, AGAT compare, and locus totals.",
        "Non-zero deltas require either a documented feature-type explanation or a parser fix.",
        "",
        f"Species checked: {len(df)}",
        f"Species with non-zero deltas: {len(flagged)}",
        f"Integrity issues: {len(integrity_issues)}",
        "",
    ]

    if not flagged.empty:
        cols = ["species"] + delta_cols
        lines.append("| " + " | ".join(cols) + " |")
        lines.append("| " + " | ".join(["---"] * len(cols)) + " |")
        for _, row in flagged[cols].iterrows():
            lines.append("| " + " | ".join(str(row.get(col, "")) for col in cols) + " |")
        lines.append("")

    if integrity_issues:
        lines.append("## Integrity Issues")
        lines.append("")
        for issue in integrity_issues:
            lines.append(f"- {issue}")
        lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate generated result tables.")
    parser.add_argument(
        "--analysis-dir",
        default=str(ANALYSIS),
        help="Analysis directory containing results/. Defaults to ANALYSIS_DIR or this directory.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    results_dir = Path(args.analysis_dir).resolve() / "results"
    results_dir.mkdir(exist_ok=True)

    df = build_validation_table(results_dir)
    integrity_issues = run_integrity_checks(results_dir)
    csv_path = results_dir / "validation_report.csv"
    md_path = results_dir / "validation_report.md"
    df.to_csv(csv_path, index=False)
    write_markdown_report(df, md_path, integrity_issues)
    print(f"Saved: {csv_path}")
    print(f"Saved: {md_path}")
    print(df.to_string(index=False))
    if integrity_issues:
        print("\nIntegrity issues:")
        for issue in integrity_issues:
            print(f"  - {issue}")
        raise SystemExit(1)
    print("\nIntegrity checks: OK")


if __name__ == "__main__":
    main()
