"""Generate tab-separated summary tables from locus_compare results."""

import pandas as pd

from analysis_config import ANALYSIS, SPECIES, SPECIES_LABELS

LOCUS_DIR = ANALYSIS / "results" / "locus"
RESULTS_DIR = ANALYSIS / "results"

SPECIES_ORDER = [sp.id for sp in SPECIES]
SPECIES_SHORT = SPECIES_LABELS


EXCLUSIVE_SUBTYPE_RULES = [
    ('exact', lambda subtype: subtype == 'exact'),
    ('boundary_refined', lambda subtype: subtype == 'boundary_refined'),
    # Keep isoform-first precedence for the mutually exclusive summary so that
    # isoform restructuring is not hidden by accompanying primary-transcript edits.
    ('isoform', lambda subtype: 'isoform_' in subtype),
    ('utr_exon_added', lambda subtype: subtype == 'utr_exon_added' or subtype.startswith('utr_exon_added_')),
    ('utr_exon_removed', lambda subtype: subtype == 'utr_exon_removed' or subtype.startswith('utr_exon_removed_')),
    ('utr_added', lambda subtype: subtype == 'utr_added' or subtype.startswith('utr_added_')),
    ('utr_lost', lambda subtype: subtype == 'utr_lost' or subtype.startswith('utr_lost_')),
    ('utr_refined', lambda subtype: subtype == 'utr_refined' or subtype.startswith('utr_refined_')),
    ('exon_gain_cds_change', lambda subtype: subtype.startswith('exon_gain_')),
    ('exon_loss_cds_change', lambda subtype: subtype.startswith('exon_loss_')),
    ('cds_boundary_refined', lambda subtype: 'cds_boundary_refined' in subtype),
    ('exon_boundary_refined', lambda subtype: 'exon_boundary_refined' in subtype),
    ('cds_change_only', lambda subtype: 'cds_extended' in subtype or 'cds_truncated' in subtype),
]


def classify_exclusive_subtype(subtype):
    """Map one detailed syntenic subtype to one broad exclusive category."""
    if subtype == 'total':
        return None
    for category, predicate in EXCLUSIVE_SUBTYPE_RULES:
        if predicate(subtype):
            return category
    raise ValueError(f"Unclassified syntenic subtype: {subtype}")


def aggregate_subtypes(row):
    """Aggregate detailed syntenic subtypes into broad categories."""
    cats = {
        'exact': 0, 'boundary_refined': 0,
        'utr_added': 0, 'utr_lost': 0,
        'utr_exon_added': 0, 'utr_exon_removed': 0,
        'utr_refined': 0,
        'exon_gain_cds_change': 0, 'exon_loss_cds_change': 0,
        'exon_boundary_refined': 0,
        'cds_change_only': 0, 'cds_boundary_refined': 0, 'isoform': 0,
    }
    for col, val in row.items():
        if pd.isna(val) or val == 0:
            continue
        if not col.startswith('syntenic_'):
            continue
        subtype = col.replace('syntenic_', '')
        category = classify_exclusive_subtype(subtype)
        if category is None:
            continue
        cats[category] += int(val)
    return cats


def aggregate_multilabel_subtypes(row):
    """Count syntenic subtype attributes independently instead of exclusively."""
    cats = {
        'Exact': 0,
        'Boundary_refined': 0,
        'Any_gene_boundary_changed': 0,
        'Any_UTR_added': 0,
        'Any_UTR_lost': 0,
        'Any_UTR_exon_gained': 0,
        'Any_UTR_exon_removed': 0,
        'Any_UTR_refined': 0,
        'Any_coding_exon_gain': 0,
        'Any_coding_exon_loss': 0,
        'Any_exon_boundary_refined': 0,
        'Any_CDS_change': 0,
        'Any_CDS_boundary_refined': 0,
        'Any_isoform_change': 0,
    }
    for col, val in row.items():
        if pd.isna(val) or val == 0 or not col.startswith('syntenic_'):
            continue
        subtype = col.replace('syntenic_', '')
        if subtype == 'total':
            continue
        classify_exclusive_subtype(subtype)
        v = int(val)
        if subtype == 'exact':
            cats['Exact'] += v
        if subtype == 'boundary_refined':
            cats['Boundary_refined'] += v
            cats['Any_gene_boundary_changed'] += v
        if subtype == 'utr_added' or subtype.startswith('utr_added_') or '_utr_added' in subtype:
            cats['Any_UTR_added'] += v
        if subtype == 'utr_lost' or subtype.startswith('utr_lost_') or '_utr_lost' in subtype:
            cats['Any_UTR_lost'] += v
        if subtype == 'utr_exon_added' or subtype.startswith('utr_exon_added_'):
            cats['Any_UTR_exon_gained'] += v
        if subtype == 'utr_exon_removed' or subtype.startswith('utr_exon_removed_'):
            cats['Any_UTR_exon_removed'] += v
        if subtype == 'utr_refined' or subtype.startswith('utr_refined_') or '_utr_refined' in subtype:
            cats['Any_UTR_refined'] += v
        if subtype.startswith('exon_gain_'):
            cats['Any_coding_exon_gain'] += v
        if subtype.startswith('exon_loss_'):
            cats['Any_coding_exon_loss'] += v
        if 'exon_boundary_refined' in subtype:
            cats['Any_exon_boundary_refined'] += v
        if 'cds_extended' in subtype or 'cds_truncated' in subtype or 'cds_boundary_refined' in subtype:
            cats['Any_CDS_change'] += v
        if 'cds_boundary_refined' in subtype:
            cats['Any_CDS_boundary_refined'] += v
        if 'isoform_' in subtype:
            cats['Any_isoform_change'] += v
    return cats


def build_master_table():
    """Build the master comparison table."""
    rows = []
    for sp in SPECIES_ORDER:
        fpath = LOCUS_DIR / f"{sp}_change_summary.csv"
        if not fpath.exists():
            continue
        df = pd.read_csv(fpath)
        row = df.iloc[0]
        cats = aggregate_subtypes(row)

        total_synt = int(row['syntenic_total'])
        subtype_sum = sum(cats.values())
        if subtype_sum != total_synt:
            raise ValueError(
                f"{sp}: exclusive syntenic subtype total {subtype_sum} "
                f"does not match syntenic_total {total_synt}"
            )

        rows.append({
            'Species': SPECIES_SHORT[sp],
            'Before': int(row['total_before_genes']),
            'After': int(row['total_after_genes']),
            'Syntenic': total_synt,
            'Exact': cats['exact'],
            'Boundary_refined_only': cats['boundary_refined'],
            'UTR_added': cats['utr_added'],
            'UTR_lost': cats['utr_lost'],
            'UTR_exon_gained': cats['utr_exon_added'],
            'UTR_exon_removed': cats['utr_exon_removed'],
            'UTR_refined': cats['utr_refined'],
            'Coding_exon_gain': cats['exon_gain_cds_change'],
            'Coding_exon_loss': cats['exon_loss_cds_change'],
            'Exon_boundary_refined': cats['exon_boundary_refined'],
            'CDS_change_only': cats['cds_change_only'],
            'CDS_boundary_refined': cats['cds_boundary_refined'],
            'Isoform_change': cats['isoform'],
            'Split_events': int(row['split_events']),
            'Merge_events': int(row['merge_events']),
            'Complex_events': int(row['complex_events']),
            'Unresolved_overlap_after_genes': int(row.get('unresolved_overlap_after_genes', 0)),
            'Unresolved_overlap_before_genes': int(row.get('unresolved_overlap_before_genes', 0)),
            'Novel_genes': int(row['novel_genes']),
            'Deleted_genes': int(row['deleted_genes']),
        })
    return pd.DataFrame(rows)


def build_multilabel_table():
    """Build a non-exclusive syntenic subtype attribute table."""
    rows = []
    direct_counts = {
        'Exact': 'one_to_one_exact',
        'Any_gene_boundary_changed': 'one_to_one_gene_boundary_changed',
        'Any_UTR_added': 'one_to_one_utr_added',
        'Any_UTR_lost': 'one_to_one_utr_lost',
        'Any_UTR_exon_gained': 'one_to_one_utr_exon_added',
        'Any_UTR_exon_removed': 'one_to_one_utr_exon_removed',
        'Any_UTR_refined': 'one_to_one_utr_refined',
        'Any_coding_exon_gain': 'one_to_one_coding_exon_gain',
        'Any_coding_exon_loss': 'one_to_one_coding_exon_loss',
        'Any_exon_boundary_refined': 'one_to_one_exon_boundary_refined',
        'Any_CDS_change': 'one_to_one_cds_change',
        'Any_CDS_boundary_refined': 'one_to_one_cds_boundary_refined',
        'Any_isoform_change': 'one_to_one_isoform_change',
    }
    for sp in SPECIES_ORDER:
        fpath = LOCUS_DIR / f"{sp}_change_summary.csv"
        if not fpath.exists():
            continue
        row = pd.read_csv(fpath).iloc[0]
        cats = aggregate_multilabel_subtypes(row)
        missing_direct_cols = [
            source_col for source_col in direct_counts.values() if source_col not in row
        ]
        if missing_direct_cols:
            raise ValueError(
                f"{sp}: missing direct multilabel count columns: "
                f"{', '.join(missing_direct_cols)}. Re-run run_locus_comparisons.py."
            )
        for output_col, source_col in direct_counts.items():
            cats[output_col] = int(row[source_col]) if pd.notna(row[source_col]) else 0
        rows.append({
            'Species': SPECIES_SHORT[sp],
            'Syntenic': int(row['syntenic_total']),
            **cats,
        })
    return pd.DataFrame(rows)


def build_locus_diagnostics_table():
    """Build overlap-mode diagnostics from per-species locus summaries."""
    rows = []
    for sp in SPECIES_ORDER:
        fpath = LOCUS_DIR / f"{sp}_change_summary.csv"
        if not fpath.exists():
            continue
        row = pd.read_csv(fpath).iloc[0]
        rows.append({
            'Species': SPECIES_SHORT[sp],
            'Overlap_mode': row.get('overlap_mode', ''),
            'Overlap_threshold': row.get('overlap_threshold', ''),
            'Candidate_pairs': row.get('candidate_pairs', ''),
            'Same_strand_overlaps': row.get('same_strand_overlaps', ''),
            'Containment_pairs_filtered_by_reciprocal': row.get(
                'containment_pairs_filtered_by_reciprocal', ''
            ),
        })
    return pd.DataFrame(rows)


def write_tsv(df, filepath):
    """Write DataFrame as formatted TSV with aligned numbers."""
    cols = df.columns.tolist()
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write('\t'.join(cols) + '\n')
        for _, row in df.iterrows():
            vals = []
            for c in cols:
                v = row[c]
                if isinstance(v, int):
                    vals.append(f'{v:,}')
                else:
                    vals.append(str(v))
            f.write('\t'.join(vals) + '\n')


def write_detailed_subtypes():
    """Write detailed subtype breakdown per species."""
    lines = []
    for sp in SPECIES_ORDER:
        fpath = LOCUS_DIR / f"{sp}_change_summary.csv"
        if not fpath.exists():
            continue
        df = pd.read_csv(fpath)
        row = df.iloc[0]

        lines.append(f"\n{'=' * 60}")
        lines.append(f"Species: {SPECIES_SHORT[sp]}")
        lines.append(f"  Total before: {int(row['total_before_genes']):,}  "
                     f"Total after: {int(row['total_after_genes']):,}  "
                     f"Syntenic: {int(row['syntenic_total']):,}")
        lines.append(f"  Splits: {int(row['split_events'])}  "
                     f"Merges: {int(row['merge_events'])}  "
                     f"Complex: {int(row['complex_events'])}  "
                     f"Unresolved after: {int(row.get('unresolved_overlap_after_genes', 0))}  "
                     f"Unresolved before: {int(row.get('unresolved_overlap_before_genes', 0))}  "
                     f"Novel: {int(row['novel_genes'])}  "
                     f"Deleted: {int(row['deleted_genes'])}")
        lines.append(f"  Syntenic subtypes:")

        subtypes = []
        for col, val in row.items():
            if col.startswith('syntenic_') and col != 'syntenic_total' and pd.notna(val) and val > 0:
                subtypes.append((col.replace('syntenic_', ''), int(val)))
        for k, v in sorted(subtypes, key=lambda x: -x[1]):
            lines.append(f"    {k}: {v:,}")

    return '\n'.join(lines)


def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    print("Building master comparison table...")
    df = build_master_table()

    # Write TSV
    tsv_path = RESULTS_DIR / "locus_comparison_summary.tsv"
    write_tsv(df, tsv_path)
    print(f"  Saved: {tsv_path}")

    # Also save CSV version
    csv_path = RESULTS_DIR / "locus_comparison_summary.csv"
    df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    multilabel_df = build_multilabel_table()
    multilabel_csv = RESULTS_DIR / "locus_comparison_multilabel.csv"
    multilabel_tsv = RESULTS_DIR / "locus_comparison_multilabel.tsv"
    multilabel_df.to_csv(multilabel_csv, index=False)
    write_tsv(multilabel_df, multilabel_tsv)
    print(f"  Saved: {multilabel_csv}")
    print(f"  Saved: {multilabel_tsv}")

    diagnostics_df = build_locus_diagnostics_table()
    diagnostics_csv = RESULTS_DIR / "locus_diagnostics.csv"
    diagnostics_tsv = RESULTS_DIR / "locus_diagnostics.tsv"
    diagnostics_df.to_csv(diagnostics_csv, index=False)
    write_tsv(diagnostics_df, diagnostics_tsv)
    print(f"  Saved: {diagnostics_csv}")
    print(f"  Saved: {diagnostics_tsv}")

    # Write detailed subtypes
    detail_text = write_detailed_subtypes()
    detail_path = RESULTS_DIR / "locus_subtype_details.txt"
    with open(detail_path, 'w', encoding='utf-8') as f:
        f.write(detail_text)
    print(f"  Saved: {detail_path}")

    # Print table for verification
    print(f"\n{'=' * 130}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
