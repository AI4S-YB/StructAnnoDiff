"""Generate publication-quality figures from locus_compare results."""

import os

import pandas as pd
import numpy as np

from analysis_config import ANALYSIS, SPECIES

MPLCONFIGDIR = ANALYSIS / ".cache" / "matplotlib"
MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns

RESULTS_DIR = ANALYSIS / "results"
FIGURES_DIR = ANALYSIS / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.15)

SPECIES_LABELS = [sp.label for sp in SPECIES]
N = len(SPECIES_LABELS)

# Color palettes
SYNTENIC_COLORS = {
    'Exact':       '#2ECC71',
    'Boundary':    '#27AE60',
    'UTR added':   '#3498DB',
    'UTR exon +/-':'#2980B9',
    'UTR refined': '#1ABC9C',
    'Coding exon +/-': '#E74C3C',
    'Exon boundary': '#F1948A',
    'CDS only':    '#E67E22',
    'CDS boundary': '#F5B041',
    'Isoform':     '#9B59B6',
}
EVENT_COLORS = {
    'Split':   '#F39C12',
    'Merge':   '#E67E22',
    'Complex': '#95A5A6',
    'Unresolved': '#7F8C8D',
    'Novel':   '#3498DB',
    'Deleted': '#E74C3C',
}


def load_data():
    """Load verified broad locus categories from generate_tables.py output."""
    summary_path = RESULTS_DIR / "locus_comparison_summary.csv"
    if not summary_path.exists():
        raise FileNotFoundError(
            f"{summary_path} not found. Run generate_tables.py before visualize_locus.py."
        )

    summary = pd.read_csv(summary_path)
    label_order = {sp.label: i for i, sp in enumerate(SPECIES)}
    summary = summary[summary["Species"].isin(label_order)].copy()
    summary["_order"] = summary["Species"].map(label_order)
    summary = summary.sort_values("_order").drop(columns="_order")

    rows = []
    exclusive_cols = [
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

    for _, row in summary.iterrows():
        total_synt = int(row["Syntenic"])
        accounted = sum(int(row[col]) for col in exclusive_cols)
        if accounted != total_synt:
            raise ValueError(
                f"{row['Species']}: exclusive syntenic categories sum to "
                f"{accounted}, expected {total_synt}"
            )

        rows.append({
            "species": row["Species"],
            "total_before": int(row["Before"]),
            "total_after": int(row["After"]),
            "syntenic_total": total_synt,
            "exact": int(row["Exact"]),
            "boundary": int(row["Boundary_refined_only"]),
            "utr_added": int(row["UTR_added"]),
            "utr_lost": int(row["UTR_lost"]),
            "utr_exon_added": int(row["UTR_exon_gained"]),
            "utr_exon_removed": int(row["UTR_exon_removed"]),
            "utr_refined": int(row["UTR_refined"]),
            "exon_gain_cds": int(row["Coding_exon_gain"]),
            "exon_loss_cds": int(row["Coding_exon_loss"]),
            "exon_boundary": int(row["Exon_boundary_refined"]),
            "cds_only": int(row["CDS_change_only"]),
            "cds_boundary": int(row["CDS_boundary_refined"]),
            "isoform": int(row["Isoform_change"]),
            "split_events": int(row["Split_events"]),
            "merge_events": int(row["Merge_events"]),
            "complex_events": int(row["Complex_events"]),
            "unresolved_after": int(row.get("Unresolved_overlap_after_genes", 0)),
            "unresolved_before": int(row.get("Unresolved_overlap_before_genes", 0)),
            "novel_genes": int(row["Novel_genes"]),
            "deleted_genes": int(row["Deleted_genes"]),
        })
    return pd.DataFrame(rows)


def plot_syntenic_subtypes(df):
    """Stacked bar chart: syntenic change subtypes as % of total genes."""
    fig, ax = plt.subplots(figsize=(16, 7))
    x = np.arange(N)
    width = 0.65

    categories = [
        ('exact', 'Exact match', SYNTENIC_COLORS['Exact']),
        ('boundary', 'Boundary refined', SYNTENIC_COLORS['Boundary']),
        ('utr_added', 'UTR added', SYNTENIC_COLORS['UTR added']),
        ('utr_lost', 'UTR lost', '#E74C3C'),
        ('utr_exon_added', 'UTR exon added', '#5DADE2'),
        ('utr_exon_removed', 'UTR exon removed', '#85C1E9'),
        ('utr_refined', 'UTR refined', SYNTENIC_COLORS['UTR refined']),
        ('exon_gain_cds', 'Coding exon gain + CDS', '#E74C3C'),
        ('exon_loss_cds', 'Coding exon loss + CDS', '#C0392B'),
        ('exon_boundary', 'Exon boundary refined', SYNTENIC_COLORS['Exon boundary']),
        ('cds_only', 'CDS change only', SYNTENIC_COLORS['CDS only']),
        ('cds_boundary', 'CDS boundary refined', SYNTENIC_COLORS['CDS boundary']),
        ('isoform', 'Isoform change', SYNTENIC_COLORS['Isoform']),
    ]

    bottom = np.zeros(N)
    for key, label, color in categories:
        vals = np.array([df[key].iloc[i] for i in range(N)])
        # Convert to percentage of total genes
        totals = np.array([df['total_before'].iloc[i] for i in range(N)])
        pcts = vals / totals * 100
        bars = ax.bar(x, pcts, width, bottom=bottom, label=label, color=color, alpha=0.88, edgecolor='white', linewidth=0.3)
        bottom += pcts

    ax.set_ylabel('% of before genes')
    ax.set_title('Syntenic gene change subtypes (% of total before-annotation genes)')
    ax.set_xticks(x)
    ax.set_xticklabels(SPECIES_LABELS, rotation=25, ha='right', fontsize=10)
    ax.legend(loc='upper left', fontsize=7.5, ncol=2, framealpha=0.9)
    ax.set_ylim(0, 105)

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "locus_syntenic_subtypes.png", dpi=200, bbox_inches='tight')
    plt.close()
    print("  Saved: locus_syntenic_subtypes.png")


def plot_event_summary(df):
    """Grouped bar chart: split, merge, complex, novel, deleted events."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))

    x = np.arange(N)
    width = 0.35

    # Panel A: Split / Merge / Complex events
    ax = axes[0]
    b1 = ax.bar(x - width, df['split_events'], width, label='Split (1→N)',
                color=EVENT_COLORS['Split'], alpha=0.85)
    b2 = ax.bar(x, df['merge_events'], width, label='Merge (N→1)',
                color=EVENT_COLORS['Merge'], alpha=0.85)
    b3 = ax.bar(x + width, df['complex_events'], width, label='Complex (M→N)',
                color=EVENT_COLORS['Complex'], alpha=0.85)
    ax.set_ylabel('Event count')
    ax.set_title('A. Gene locus restructuring events')
    ax.set_xticks(x)
    ax.set_xticklabels(SPECIES_LABELS, rotation=25, ha='right', fontsize=9)
    ax.legend(fontsize=8)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{v/1000:.1f}K' if v >= 1000 else f'{v:.0f}'))

    # Panel B: unresolved weak-overlap / novel / deleted genes
    ax = axes[1]
    narrow = 0.18
    ax.bar(x - 1.5 * narrow, df['unresolved_after'], narrow, label='Unresolved after',
           color=EVENT_COLORS['Unresolved'], alpha=0.75)
    ax.bar(x - 0.5 * narrow, df['unresolved_before'], narrow, label='Unresolved before',
           color='#AAB2B6', alpha=0.75)
    ax.bar(x + 0.5 * narrow, df['novel_genes'], narrow, label='Novel genes',
                color=EVENT_COLORS['Novel'], alpha=0.85)
    ax.bar(x + 1.5 * narrow, df['deleted_genes'], narrow, label='Deleted genes',
                color=EVENT_COLORS['Deleted'], alpha=0.85)
    ax.set_ylabel('Gene count')
    ax.set_title('B. Unresolved weak-overlap, novel, and deleted genes')
    ax.set_xticks(x)
    ax.set_xticklabels(SPECIES_LABELS, rotation=25, ha='right', fontsize=9)
    ax.legend(fontsize=8)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{v/1000:.1f}K' if v >= 1000 else f'{v:.0f}'))

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "locus_events.png", dpi=200, bbox_inches='tight')
    plt.close()
    print("  Saved: locus_events.png")


def plot_heatmap(df):
    """Heatmap showing % of genes in each category across species."""
    # Categories to show
    cats = ['exact', 'boundary', 'utr_added', 'utr_lost',
            'utr_exon_added', 'utr_exon_removed', 'utr_refined',
            'exon_gain_cds', 'exon_loss_cds', 'exon_boundary',
            'cds_only', 'cds_boundary', 'isoform',
            'split_events', 'merge_events', 'unresolved_after', 'unresolved_before',
            'novel_genes', 'deleted_genes']
    labels = ['Exact', 'Boundary\nrefined', 'UTR\nadded', 'UTR\nlost',
              'UTR exon\ngained', 'UTR exon\nremoved', 'UTR\nrefined',
              'Coding exon\ngain+CDS', 'Coding exon\nloss+CDS',
              'Exon boundary\nrefined', 'CDS change\nonly',
              'CDS boundary\nrefined', 'Isoform\nchange', 'Split\nevents', 'Merge\nevents',
              'Unresolved\nafter', 'Unresolved\nbefore', 'Novel\ngenes', 'Deleted\ngenes']

    data = np.zeros((N, len(cats)))
    for i in range(N):
        row = df.iloc[i]
        for j, cat in enumerate(cats):
            val = row[cat]
            total = row['total_before']
            data[i, j] = val / total * 100 if total > 0 else 0

    fig, ax = plt.subplots(figsize=(18, 7))
    cmap = sns.color_palette("YlOrRd", as_cmap=True)
    im = ax.imshow(data, cmap=cmap, aspect='auto', vmin=0, vmax=np.percentile(data[data > 0], 95))

    ax.set_xticks(np.arange(len(cats)))
    ax.set_xticklabels(labels, rotation=35, ha='right', fontsize=8.5)
    ax.set_yticks(np.arange(N))
    ax.set_yticklabels(SPECIES_LABELS, fontsize=10)

    # Annotate with values
    for i in range(N):
        for j in range(len(cats)):
            val = data[i, j]
            if val > 0.05:
                text = f'{val:.1f}'
                color = 'white' if val > np.percentile(data[data > 0], 60) else 'black'
                size = 7 if val < 10 else 7.5
                ax.text(j, i, text, ha='center', va='center', fontsize=size, color=color)

    ax.set_title('Gene annotation changes across species (% of before-annotation genes)', fontsize=13)
    fig.colorbar(im, ax=ax, label='% of before genes', shrink=0.78)
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "locus_heatmap.png", dpi=200, bbox_inches='tight')
    plt.close()
    print("  Saved: locus_heatmap.png")


def plot_conservation_index(df):
    """Bar chart: conservation index = (exact + boundary) / total_before * 100."""
    fig, ax = plt.subplots(figsize=(12, 6))

    x = np.arange(N)
    width = 0.55

    conserved = (df['exact'] + df['boundary']) / df['total_before'] * 100
    colors = ['#2ECC71' if v > 50 else '#F39C12' if v > 20 else '#E74C3C' for v in conserved]

    bars = ax.bar(x, conserved, width, color=colors, alpha=0.85, edgecolor='white', linewidth=0.5)
    ax.set_ylabel('% of before genes')
    ax.set_title('Annotation conservation index (exact + boundary refined genes)')
    ax.set_xticks(x)
    ax.set_xticklabels(SPECIES_LABELS, rotation=25, ha='right', fontsize=10)
    ax.axhline(y=50, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
    ax.axhline(y=90, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)

    # Annotate bars
    for i, (bar, val) in enumerate(zip(bars, conserved)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.8,
                f'{val:.1f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax.set_ylim(0, 105)
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "locus_conservation.png", dpi=200, bbox_inches='tight')
    plt.close()
    print("  Saved: locus_conservation.png")


def plot_gene_count_comparison(df):
    """Before vs after gene counts with delta annotation."""
    fig, ax = plt.subplots(figsize=(12, 6))

    x = np.arange(N)
    width = 0.3

    before = df['total_before'].values
    after = df['total_after'].values

    ax.bar(x - width/2, before, width, label='Before curation', color='#E74C3C', alpha=0.85)
    ax.bar(x + width/2, after, width, label='After curation', color='#2ECC71', alpha=0.85)

    # Delta annotations
    for i in range(N):
        delta = after[i] - before[i]
        pct = delta / before[i] * 100
        sign = '+' if delta >= 0 else ''
        y = max(before[i], after[i])
        color = '#27AE60' if delta >= 0 else '#C0392B'
        ax.annotate(f'{sign}{delta:,} ({sign}{pct:.1f}%)', (x[i], y),
                    textcoords="offset points", xytext=(0, 6),
                    ha='center', fontsize=8, color=color, fontweight='bold')

    ax.set_ylabel('Gene count')
    ax.set_title('Total gene count before vs after manual curation')
    ax.set_xticks(x)
    ax.set_xticklabels(SPECIES_LABELS, rotation=25, ha='right', fontsize=10)
    ax.legend(fontsize=9)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{v/1000:.0f}K'))

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "locus_gene_counts.png", dpi=200, bbox_inches='tight')
    plt.close()
    print("  Saved: locus_gene_counts.png")


def main():
    global SPECIES_LABELS, N

    print("Loading locus comparison data...")
    df = load_data()
    print(f"  Loaded {len(df)} species")
    if df.empty:
        print("ERROR: no locus comparison data found. Run run_locus_comparisons.py first.")
        return

    SPECIES_LABELS = df['species'].tolist()
    N = len(df)

    print("\nGenerating figures...")
    plot_syntenic_subtypes(df)
    plot_event_summary(df)
    plot_heatmap(df)
    plot_conservation_index(df)
    plot_gene_count_comparison(df)

    print(f"\nAll figures saved to {FIGURES_DIR}/")


if __name__ == "__main__":
    main()
