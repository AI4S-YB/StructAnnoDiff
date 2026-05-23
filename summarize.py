"""Aggregate AGAT and gffcompare outputs into summary CSV tables."""

import argparse
import re
from pathlib import Path
from collections import Counter, defaultdict

import pandas as pd
import numpy as np

from analysis_config import ANALYSIS as DEFAULT_ANALYSIS
from analysis_config import SPECIES_IDS, find_annotation
from gff_utils import parse_gff3

ANALYSIS = DEFAULT_ANALYSIS
STATS_DIR = ANALYSIS / "stats"
COMPARE_DIR = ANALYSIS / "compare"
TCOMPARE_DIR = ANALYSIS / "tcompare"
RESULTS_DIR = ANALYSIS / "results"
SPECIES_LIST = SPECIES_IDS
CODING_STATS_SECTIONS = ("mrna", "transcript")
CODING_COMPARE_PATHS = (
    "gene@mrna@cds",
    "gene@mrna@exon",
    "gene@transcript@cds",
    "gene@transcript@exon",
)


def configure_paths(analysis_dir):
    """Configure module-level paths from a chosen analysis directory."""
    global ANALYSIS, STATS_DIR, COMPARE_DIR, TCOMPARE_DIR, RESULTS_DIR
    ANALYSIS = Path(analysis_dir).resolve()
    STATS_DIR = ANALYSIS / "stats"
    COMPARE_DIR = ANALYSIS / "compare"
    TCOMPARE_DIR = ANALYSIS / "tcompare"
    RESULTS_DIR = ANALYSIS / "results"
    RESULTS_DIR.mkdir(exist_ok=True)


def _parse_metric_value(value):
    try:
        return float(value)
    except ValueError:
        return value


def parse_agat_statistics_sections(filepath):
    """Parse AGAT statistics output into section -> metric -> value.

    AGAT repeats metric names across feature sections (for example mrna and
    trna both have "Number of gene"). Keeping the section avoids silently
    overwriting coding-gene metrics with later non-coding sections.
    """
    sections = defaultdict(dict)
    current_section = "global"

    with open(filepath, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            section_match = re.match(r'^-+\s+(.+?)\s+-+$', line)
            if section_match:
                current_section = section_match.group(1).strip().lower().replace(' ', '_')
                continue
            if line.startswith('-'):
                continue
            parts = re.split(r'\s{2,}', line)
            if len(parts) >= 2:
                key = parts[0].strip()
                sections[current_section][key] = _parse_metric_value(parts[-1].strip())
    return dict(sections)


def select_coding_stats_section(sections):
    """Return the preferred coding section name and metrics."""
    for section in CODING_STATS_SECTIONS:
        if section in sections:
            return section, sections[section]
    for section, metrics in sections.items():
        if "Number of gene" in metrics:
            return section, metrics
    return "", {}


def build_stats_by_section_table():
    """Build long-form AGAT statistics by species, state, section, and metric."""
    rows = []
    for sp in SPECIES_LIST:
        for state in ['before', 'after']:
            stat_file = STATS_DIR / f"{sp}.{state}"
            if not stat_file.exists():
                continue
            sections = parse_agat_statistics_sections(stat_file)
            for section, metrics in sections.items():
                for metric, value in metrics.items():
                    rows.append({
                        'species': sp,
                        'state': state,
                        'section': section,
                        'metric': metric,
                        'value': value,
                    })
    return pd.DataFrame(rows)


def build_stats_table():
    """Build a master coding statistics table: rows=species, cols=metric_before/after."""
    all_data = []
    for sp in SPECIES_LIST:
        row = {'species': sp}
        for state in ['before', 'after']:
            stat_file = STATS_DIR / f"{sp}.{state}"
            if not stat_file.exists():
                print(f"  WARNING: {stat_file} not found")
                continue
            sections = parse_agat_statistics_sections(stat_file)
            section, metrics = select_coding_stats_section(sections)
            row[f'coding_section_{state}'] = section
            for key, val in metrics.items():
                col_name = f"{key}_{state}"
                row[col_name] = val
        all_data.append(row)
    return pd.DataFrame(all_data)


def parse_agat_report_tables(filepath):
    """Parse AGAT report.txt into feature_path -> (n1,n2) -> count."""
    tables = {}
    current_table = None

    with open(filepath, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()

            parts = [p.strip() for p in line.split('|') if p.strip()] if '|' in line else []
            if len(parts) == 1 and '@' in parts[0] and 'Number of cases' not in parts[0]:
                current_table = parts[0].lower()
                tables.setdefault(current_table, {})
                continue

            if current_table and len(parts) >= 3:
                try:
                    n1 = int(parts[0])
                    n2 = int(parts[1])
                    count = int(parts[2])
                except ValueError:
                    continue
                tables[current_table][(n1, n2)] = count
    return tables


def summarize_cases(cases):
    """Summarize AGAT case counts into broad before/after categories."""
    return {
        'removed': cases.get((1, 0), 0),
        'added': cases.get((0, 1), 0),
        'concordant': cases.get((1, 1), 0),
        'split': sum(v for (n1, n2), v in cases.items() if n1 == 1 and n2 > 1),
        'fusion': sum(v for (n1, n2), v in cases.items() if n1 > 1 and n2 == 1),
        'complex': sum(v for (n1, n2), v in cases.items() if n1 > 1 and n2 > 1),
        'total_before': sum(n1 * v for (n1, n2), v in cases.items() if n1 > 0),
        'total_after': sum(n2 * v for (n1, n2), v in cases.items() if n2 > 0),
    }


def combine_cases(tables, feature_paths):
    combined = defaultdict(int)
    for feature_path in feature_paths:
        for key, value in tables.get(feature_path, {}).items():
            combined[key] += value
    return dict(combined)


def build_comparison_by_feature_path_table():
    """Build long-form AGAT comparison summaries for each feature path table."""
    rows = []
    for sp in SPECIES_LIST:
        report_file = COMPARE_DIR / sp / "report.txt"
        if not report_file.exists():
            continue
        tables = parse_agat_report_tables(report_file)
        for feature_path, cases in sorted(tables.items()):
            row = {'species': sp, 'feature_path': feature_path}
            row.update(summarize_cases(cases))
            rows.append(row)
    return pd.DataFrame(rows)


def build_comparison_all_gene_types_table():
    """Build broad AGAT comparison across all feature-path tables."""
    rows = []
    for sp in SPECIES_LIST:
        report_file = COMPARE_DIR / sp / "report.txt"
        if not report_file.exists():
            continue
        tables = parse_agat_report_tables(report_file)
        combined = defaultdict(int)
        for cases in tables.values():
            for key, value in cases.items():
                combined[key] += value
        row = {'species': sp}
        row.update(summarize_cases(combined))
        rows.append(row)
    return pd.DataFrame(rows)


def build_comparison_table():
    """Build AGAT comparison counts for mRNA/transcript-bearing gene paths."""
    rows = []
    for sp in SPECIES_LIST:
        report_file = COMPARE_DIR / sp / "report.txt"
        if not report_file.exists():
            print(f"  WARNING: {report_file} not found")
            continue

        tables = parse_agat_report_tables(report_file)
        cases = combine_cases(tables, CODING_COMPARE_PATHS)
        row = {'species': sp}
        row.update(summarize_cases(cases))
        rows.append(row)

    return pd.DataFrame(rows)


def parse_gffcompare_stats(filepath):
    """Parse gffcompare .stats file."""
    metrics = {}
    with open(filepath) as f:
        content = f.read()

    # Lines like "        Base level:    94.6     |    87.8    |"
    for match in re.finditer(r'(\w[\w\s]+level):\s+([\d.]+)\s*\|\s*([\d.]+)', content):
        level = match.group(1).strip().lower().replace(' ', '_')
        metrics[f'{level}_sensitivity'] = float(match.group(2))
        metrics[f'{level}_precision'] = float(match.group(3))

    # Count novel/missed features
    for match in re.finditer(r'(Missed|Novel)\s+(\w+):\s+(\d+)/(\d+)\s+([\d.]+)%', content):
        key = f'{match.group(1).lower()}_{match.group(2).lower()}'
        metrics[key] = int(match.group(3))
        metrics[f'{key}_total'] = int(match.group(4))
        metrics[f'{key}_pct'] = float(match.group(5))

    # Matching counts
    for match in re.finditer(r'Matching\s+(\w[\w\s]+):\s+(\d+)', content):
        key = f'matching_{match.group(1).strip().lower().replace(" ", "_")}'
        metrics[key] = int(match.group(2))

    return metrics


def parse_gffcompare_class_codes(filepath):
    """Extract class code distribution from annotated GTF."""
    codes = defaultdict(int)
    if not filepath.exists():
        return dict(codes)
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            m = re.search(r'class_code\s+"([^"]+)"', line)
            if m:
                codes[m.group(1)] += 1
    return dict(codes)


def build_gffcompare_table():
    """Build table of gffcompare accuracy metrics per species."""
    rows = []
    for sp in SPECIES_LIST:
        stats_file = TCOMPARE_DIR / f"{sp}.stats"
        gtf_file = TCOMPARE_DIR / f"{sp}.annotated.gtf"
        if stats_file.exists():
            metrics = parse_gffcompare_stats(stats_file)
            metrics['species'] = sp
            # Add class codes
            if gtf_file.exists():
                codes = parse_gffcompare_class_codes(gtf_file)
                for code, count in codes.items():
                    metrics[f'class_{code}'] = count
            rows.append(metrics)
        else:
            print(f"  WARNING: {stats_file} not found")
    return pd.DataFrame(rows)


def _attribute_values(attrs, name):
    value = attrs.get(name)
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _top_counts_text(counter, limit=10):
    return ' | '.join(f'{count} {value}' for value, count in counter.most_common(limit))


def extract_quality_from_gff(filepath, aed_limit=5000):
    """Extract quality attribute summaries from one annotation file."""
    manual_scores = Counter()
    maninfo = Counter()
    aed_values = []

    for feat in parse_gff3(filepath):
        attrs = feat['attributes']

        for value in _attribute_values(attrs, 'manual_quality_score'):
            manual_scores[value] += 1
        for value in _attribute_values(attrs, 'ManInfo'):
            maninfo[value] += 1
        if len(aed_values) < aed_limit:
            for value in _attribute_values(attrs, '_AED'):
                if len(aed_values) >= aed_limit:
                    break
                try:
                    aed_values.append(float(value))
                except ValueError:
                    continue

    return manual_scores, maninfo, aed_values


def extract_quality_attributes():
    """Extract quality attributes from GFF3 files."""
    rows = []
    for sp in SPECIES_LIST:
        row = {'species': sp}
        for state in ['before', 'after']:
            gff = find_annotation(sp, state, ANALYSIS)
            if gff is None:
                continue

            manual_scores, maninfo, aed_vals = extract_quality_from_gff(gff)
            if manual_scores:
                row[f'quality_scores_{state}'] = _top_counts_text(manual_scores)
            if maninfo:
                row[f'maninfo_{state}'] = _top_counts_text(maninfo)
            if aed_vals:
                row[f'aed_mean_{state}'] = round(np.mean(aed_vals), 4)
                row[f'aed_median_{state}'] = round(np.median(aed_vals), 4)

        rows.append(row)
    return pd.DataFrame(rows)


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate AGAT and gffcompare output tables.")
    parser.add_argument(
        "--analysis-dir",
        default=str(DEFAULT_ANALYSIS),
        help="Analysis directory containing stats/, compare/, tcompare/, and annotation files.",
    )
    return parser.parse_args()


def main(analysis_dir=None):
    configure_paths(analysis_dir or DEFAULT_ANALYSIS)

    print("=" * 60)
    print("1. Building statistics summary table...")
    print("=" * 60)
    stats_df = build_stats_table()
    stats_path = RESULTS_DIR / "summary_stats.csv"
    stats_df.to_csv(stats_path, index=False)
    print(f"  Saved: {stats_path} ({len(stats_df)} rows, {len(stats_df.columns)} cols)")

    stats_section_df = build_stats_by_section_table()
    stats_section_path = RESULTS_DIR / "summary_stats_by_section.csv"
    stats_section_df.to_csv(stats_section_path, index=False)
    print(f"  Saved: {stats_section_path} ({len(stats_section_df)} rows)")

    # Print key metrics comparison
    key_metrics = [
        ('Number of gene', 'Gene count'),
        ('Number of mrna', 'mRNA count'),
        ('mean gene length (bp)', 'Mean gene len'),
        ('mean cds length (bp)', 'Mean CDS len'),
        ('mean exons per mrna', 'Mean exons/mRNA'),
        ('Number of single exon gene', 'Single-exon genes'),
        ('Number gene overlapping', 'Overlapping genes'),
    ]
    print(f"\n{'Species':<25}", end='')
    for _, label in key_metrics:
        print(f'  {label:<18}', end='')
    print()
    for sp in SPECIES_LIST:
        row = stats_df[stats_df['species'] == sp]
        if row.empty:
            continue
        row = row.iloc[0]
        print(f'{sp:<25}', end='')
        for metric, _ in key_metrics:
            b = row.get(f'{metric}_before', None)
            a = row.get(f'{metric}_after', None)
            if b is not None and a is not None and isinstance(b, (int, float)) and isinstance(a, (int, float)):
                if isinstance(b, float) or isinstance(a, float):
                    print(f'  {b:.1f}→{a:.1f}', end='')
                else:
                    delta = a - b
                    sign = '+' if delta >= 0 else ''
                    print(f'  {b}→{a} ({sign}{delta})', end='')
            else:
                print(f'  N/A', end='')
        print()

    print(f"\n{'=' * 60}")
    print("2. Building AGAT comparison tables...")
    print(f"{'=' * 60}")
    comp_df = build_comparison_table()
    comp_path = RESULTS_DIR / "comparison_matrix.csv"
    comp_df.to_csv(comp_path, index=False)
    print(f"  Saved: {comp_path} ({len(comp_df)} mRNA/transcript-scope rows)")

    comp_feature_df = build_comparison_by_feature_path_table()
    comp_feature_path = RESULTS_DIR / "comparison_by_feature_path.csv"
    comp_feature_df.to_csv(comp_feature_path, index=False)
    print(f"  Saved: {comp_feature_path} ({len(comp_feature_df)} rows)")

    comp_all_df = build_comparison_all_gene_types_table()
    comp_all_path = RESULTS_DIR / "comparison_matrix_all_gene_types.csv"
    comp_all_df.to_csv(comp_all_path, index=False)
    print(f"  Saved: {comp_all_path} ({len(comp_all_df)} rows)")

    for _, row in comp_df.iterrows():
        print(f"\n  {row['species']}:")
        print(f"    Removed (1→0): {row['removed']:<8} Added (0→1): {row['added']}")
        print(f"    Concordant (1→1): {row['concordant']:<8} Split (1→many): {row['split']}")
        print(f"    Fusion (many→1): {row['fusion']:<8} Complex: {row['complex']}")

    print(f"\n{'=' * 60}")
    print("3. Building gffcompare accuracy table...")
    print(f"{'=' * 60}")
    gc_df = build_gffcompare_table()
    gc_path = RESULTS_DIR / "accuracy_metrics.csv"
    gc_df.to_csv(gc_path, index=False)
    print(f"  Saved: {gc_path} ({len(gc_df)} rows, {len(gc_df.columns)} cols)")

    # Print key gffcompare metrics
    accuracy_cols = [c for c in gc_df.columns if 'sensitivity' in c or 'precision' in c]
    if accuracy_cols:
        print(f"\n{'Species':<25}", end='')
        for c in accuracy_cols:
            print(f'  {c:<25}', end='')
        print()
        for _, row in gc_df.iterrows():
            print(f'{row.get("species", ""):<25}', end='')
            for c in accuracy_cols:
                v = row.get(c, None)
                if v is not None:
                    print(f'  {v:<25.1f}%', end='')
                else:
                    print(f'  {"N/A":<25}', end='')
            print()

    print(f"\n{'=' * 60}")
    print("4. Extracting quality attributes...")
    print(f"{'=' * 60}")
    qual_df = extract_quality_attributes()
    qual_path = RESULTS_DIR / "quality_attributes.csv"
    qual_df.to_csv(qual_path, index=False)
    print(f"  Saved: {qual_path} ({len(qual_df)} rows)")

    # Print quality highlights
    for _, row in qual_df.iterrows():
        print(f"\n  {row['species']}:")
        if 'quality_scores_after' in row and pd.notna(row['quality_scores_after']):
            print(f"    Manual quality scores (after): {row['quality_scores_after'][:120]}")
        if 'maninfo_after' in row and pd.notna(row['maninfo_after']):
            print(f"    ManInfo (after): {row['maninfo_after'][:120]}")
        if 'aed_mean_before' in row and pd.notna(row['aed_mean_before']):
            print(f"    AED before: mean={row['aed_mean_before']}, median={row['aed_median_before']}")
        if 'aed_mean_after' in row and pd.notna(row['aed_mean_after']):
            print(f"    AED after: mean={row['aed_mean_after']}, median={row['aed_median_after']}")

    print(f"\n{'=' * 60}")
    print("All tables saved to results/")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    args = parse_args()
    main(args.analysis_dir)
