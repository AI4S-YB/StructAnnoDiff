#!/usr/bin/env python3
"""Validate cross-table consistency for annotation comparison outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from analysis_config import ANALYSIS, RESULTS_DIR, SPECIES, SPECIES_LABELS, find_annotation
from gff_utils import normalize_attribute_id, normalize_attribute_ids, parse_gff3

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


def _as_float(value) -> float:
    return float(value) if pd.notna(value) else 0.0


def _locus_with_species_ids(locus_df: pd.DataFrame) -> pd.DataFrame:
    if locus_df.empty:
        return locus_df
    label_to_id = {sp.label: sp.id for sp in SPECIES}
    locus_df = locus_df.copy()
    locus_df["species"] = locus_df["Species"].map(label_to_id)
    return locus_df


def _attr(attrs, key, default=""):
    return attrs.get(key, attrs.get(key.lower(), default))


def count_genes_filtered_by_locus_evidence(annotation_path: Path | None) -> int:
    """Count mRNA-scope genes excluded because no mRNA has exon/CDS evidence."""
    if annotation_path is None or not annotation_path.exists():
        return 0

    gene_ids = set()
    gene_mrnas = {}
    mrna_has_locus_evidence = {}

    for feat in parse_gff3(annotation_path):
        ftype = feat["type"]
        attrs = feat["attributes"]
        if ftype == "gene":
            gid = normalize_attribute_id(_attr(attrs, "ID"), prefixes=("gene:",))
            if gid:
                gene_ids.add(gid)
                gene_mrnas.setdefault(gid, set())
        elif ftype in ("mRNA", "transcript"):
            mid = normalize_attribute_id(_attr(attrs, "ID"), prefixes=("transcript:",))
            if not mid:
                continue
            parents = normalize_attribute_ids(_attr(attrs, "Parent"), prefixes=("gene:",))
            mrna_has_locus_evidence.setdefault(mid, False)
            for parent in parents:
                gene_mrnas.setdefault(parent, set()).add(mid)
        elif ftype in ("exon", "CDS"):
            parent_ids = normalize_attribute_ids(_attr(attrs, "Parent"), prefixes=("transcript:",))
            for parent in parent_ids:
                mrna_has_locus_evidence[parent] = True

    mrna_scope_genes = {
        gid for gid in gene_ids
        if gene_mrnas.get(gid)
    }
    locus_supported_genes = {
        gid for gid in mrna_scope_genes
        if any(mrna_has_locus_evidence.get(mid, False) for mid in gene_mrnas.get(gid, ()))
    }
    return len(mrna_scope_genes - locus_supported_genes)


def build_validation_table(results_dir: Path = RESULTS_DIR) -> pd.DataFrame:
    analysis_dir = results_dir.parent
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
        row["locus_evidence_filtered_before"] = count_genes_filtered_by_locus_evidence(
            find_annotation(sp.id, "before", analysis_dir)
        )
        row["locus_evidence_filtered_after"] = count_genes_filtered_by_locus_evidence(
            find_annotation(sp.id, "after", analysis_dir)
        )

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

        for state in ("before", "after"):
            delta_name = f"stats_vs_locus_{state}_delta"
            filtered_name = f"locus_evidence_filtered_{state}"
            unexplained_name = f"stats_vs_locus_{state}_unexplained_delta"
            if delta_name in row:
                row[unexplained_name] = row[delta_name] - row.get(filtered_name, 0)

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

    curation_core = _read_csv(results_dir / "curation_core_metrics.csv")
    if curation_core.empty:
        issues.append("curation_core_metrics.csv: missing or empty")
    else:
        observed_ids = curation_core.get("species_id", pd.Series(dtype=str)).tolist()
        expected_ids = [sp.id for sp in SPECIES]
        if observed_ids != expected_ids:
            issues.append(
                "curation_core_metrics.csv: species order/content mismatch "
                f"{observed_ids} != {expected_ids}"
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

        expected_changed_before = (
            _as_int(row["total_before_genes"]) - _as_int(row["one_to_one_exact"])
        )
        expected_changed_after = (
            _as_int(row["total_after_genes"]) - _as_int(row["one_to_one_exact"])
        )
        if _as_int(row.get("changed_before_genes", 0)) != expected_changed_before:
            issues.append(f"{sp.id}: changed_before_genes does not match exact 1:1 complement")
        if _as_int(row.get("changed_after_genes", 0)) != expected_changed_after:
            issues.append(f"{sp.id}: changed_after_genes does not match exact 1:1 complement")
        if _as_int(row.get("no_overlap_after_loci", 0)) > _as_int(row["total_after_genes"]):
            issues.append(f"{sp.id}: no_overlap_after_loci exceeds total_after_genes")
        if _as_int(row.get("no_overlap_before_loci", 0)) > _as_int(row["total_before_genes"]):
            issues.append(f"{sp.id}: no_overlap_before_loci exceeds total_before_genes")
        rep_pairs = _as_int(row.get("rep_transcript_pairs", 0))
        rep_structural_changed = _as_int(row.get("rep_structural_changed", 0))
        if rep_structural_changed > rep_pairs:
            issues.append(f"{sp.id}: representative structural changes exceed pair count")
        expected_rep_structural_pct = (
            rep_structural_changed / rep_pairs * 100 if rep_pairs else 0.0
        )
        if abs(
            _as_float(row.get("rep_structural_changed_pct", 0.0))
            - expected_rep_structural_pct
        ) > 1e-9:
            issues.append(f"{sp.id}: representative structural change percent mismatch")
        if (
            _as_int(row.get("rep_exon_count_changed", 0))
            + _as_int(row.get("rep_exon_boundary_changed_same_count", 0))
            > rep_pairs
        ):
            issues.append(f"{sp.id}: representative exon change categories exceed pair count")
        rep_exon_changed = (
            _as_int(row.get("rep_exon_count_changed", 0))
            + _as_int(row.get("rep_exon_boundary_changed_same_count", 0))
        )
        rep_cds_changed = (
            _as_int(row.get("rep_cds_count_changed", 0))
            + _as_int(row.get("rep_cds_boundary_changed_same_count", 0))
        )
        if (
            rep_cds_changed > rep_pairs
        ):
            issues.append(f"{sp.id}: representative CDS change categories exceed pair count")
        if rep_structural_changed < max(rep_exon_changed, rep_cds_changed):
            issues.append(f"{sp.id}: representative structural changes undercount exon/CDS changes")
        if rep_structural_changed > rep_exon_changed + rep_cds_changed:
            issues.append(f"{sp.id}: representative structural changes exceed exon/CDS change union")

        if not curation_core.empty:
            core_row = curation_core[curation_core["species_id"] == sp.id]
            if not core_row.empty:
                core_row = core_row.iloc[0]
                core_to_summary = {
                    "total_before_genes": "total_before_genes",
                    "total_after_genes": "total_after_genes",
                    "changed_before_genes": "changed_before_genes",
                    "changed_after_genes": "changed_after_genes",
                    "new_loci_no_overlap": "no_overlap_after_loci",
                    "deleted_loci_no_overlap": "no_overlap_before_loci",
                    "split_events": "split_events",
                    "merge_events": "merge_events",
                    "rep_transcript_pairs": "rep_transcript_pairs",
                    "rep_structural_changed": "rep_structural_changed",
                    "rep_exon_count_changed": "rep_exon_count_changed",
                    "rep_exon_boundary_changed_same_count": "rep_exon_boundary_changed_same_count",
                    "rep_cds_count_changed": "rep_cds_count_changed",
                    "rep_cds_boundary_changed_same_count": "rep_cds_boundary_changed_same_count",
                }
                for core_col, summary_col in core_to_summary.items():
                    if core_col not in core_row.index:
                        issues.append(f"{sp.id}: curation_core_metrics missing {core_col}")
                        continue
                    if _as_int(core_row[core_col]) != _as_int(row[summary_col]):
                        issues.append(
                            f"{sp.id}: curation_core_metrics {core_col} "
                            f"!= {summary_col}"
                        )
                if "rep_exon_changed" not in core_row.index:
                    issues.append(f"{sp.id}: curation_core_metrics missing rep_exon_changed")
                elif _as_int(core_row["rep_exon_changed"]) != rep_exon_changed:
                    issues.append(f"{sp.id}: curation_core_metrics rep_exon_changed mismatch")
                pct_checks = {
                    "rep_exon_changed_before_pct": (
                        rep_exon_changed / _as_int(row["total_before_genes"]) * 100
                        if _as_int(row["total_before_genes"]) else 0.0
                    ),
                    "rep_exon_changed_after_pct": (
                        rep_exon_changed / _as_int(row["total_after_genes"]) * 100
                        if _as_int(row["total_after_genes"]) else 0.0
                    ),
                }
                for core_col, expected_value in pct_checks.items():
                    if core_col not in core_row.index:
                        issues.append(f"{sp.id}: curation_core_metrics missing {core_col}")
                    elif abs(_as_float(core_row[core_col]) - expected_value) > 1e-9:
                        issues.append(f"{sp.id}: curation_core_metrics {core_col} mismatch")
                if "rep_structural_changed_pct" not in core_row.index:
                    issues.append(f"{sp.id}: curation_core_metrics missing rep_structural_changed_pct")
                elif abs(
                    _as_float(core_row["rep_structural_changed_pct"])
                    - _as_float(row["rep_structural_changed_pct"])
                ) > 1e-9:
                    issues.append(
                        f"{sp.id}: curation_core_metrics rep_structural_changed_pct "
                        "!= rep_structural_changed_pct"
                    )

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
    unexplained_cols = [c for c in df.columns if c.endswith("_unexplained_delta")]
    flagged = (
        df[df[unexplained_cols].fillna(0).abs().sum(axis=1) > 0]
        if unexplained_cols else pd.DataFrame()
    )
    integrity_issues = integrity_issues or []

    lines = [
        "# Validation Report",
        "",
        "This report compares mRNA/transcript-scope AGAT stats, AGAT compare, and locus totals.",
        "Locus totals exclude genes whose mRNA/transcript records have no explicit exon or CDS evidence.",
        "Non-zero unexplained deltas require either a documented feature-type explanation or a parser fix.",
        "",
        f"Species checked: {len(df)}",
        f"Species with non-zero unexplained deltas: {len(flagged)}",
        f"Integrity issues: {len(integrity_issues)}",
        "",
    ]

    if not flagged.empty:
        cols = [
            "species",
            "locus_evidence_filtered_before",
            "locus_evidence_filtered_after",
        ] + unexplained_cols
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
