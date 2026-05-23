"""Generate publication-quality figures comparing before vs after annotations."""

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
sns.set_context("paper", font_scale=1.2)

SPECIES_IDS = [sp.id for sp in SPECIES]
SPECIES_LABELS = [sp.label.replace(" ", "\n") for sp in SPECIES]
SPECIES_SHORT = [sp.short_label for sp in SPECIES]


def plot_gene_counts(df_stats, df_comp):
    """Grouped bar chart: gene count before/after + added/removed."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Panel A: Gene counts
    ax = axes[0]
    x = np.arange(len(SPECIES_SHORT))
    width = 0.35
    before_genes = []
    after_genes = []
    for sp_label, sp_file in zip(SPECIES_SHORT, SPECIES_IDS):
        row = df_stats[df_stats['species'] == sp_file]
        if not row.empty:
            row = row.iloc[0]
            before_genes.append(row.get('Number of gene_before', 0))
            after_genes.append(row.get('Number of gene_after', 0))
        else:
            before_genes.append(0)
            after_genes.append(0)

    bars1 = ax.bar(x - width / 2, before_genes, width, label='Before', color='#E74C3C', alpha=0.85)
    bars2 = ax.bar(x + width / 2, after_genes, width, label='After', color='#2ECC71', alpha=0.85)
    ax.set_ylabel('Gene count')
    ax.set_title('A. Gene count before vs after manual correction')
    ax.set_xticks(x)
    ax.set_xticklabels(SPECIES_SHORT, rotation=30, ha='right', fontsize=9)
    ax.legend()
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{v/1000:.0f}K'))

    # Panel B: Added/removed genes from AGAT comparison
    ax = axes[1]
    if not df_comp.empty:
        added = []
        removed = []
        for sp_file in SPECIES_IDS:
            row = df_comp[df_comp['species'] == sp_file]
            if not row.empty:
                row = row.iloc[0]
                added.append(row.get('added', 0))
                removed.append(row.get('removed', 0))
            else:
                added.append(0)
                removed.append(0)

        bars3 = ax.bar(x - width / 2, removed, width, label='Removed (1→0)', color='#E74C3C', alpha=0.85)
        bars4 = ax.bar(x + width / 2, added, width, label='Added (0→1)', color='#3498DB', alpha=0.85)
        ax.set_ylabel('Gene count')
        ax.set_title('B. Genes removed/added during curation')
        ax.set_xticks(x)
        ax.set_xticklabels(SPECIES_SHORT, rotation=30, ha='right', fontsize=9)
        ax.legend()

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "gene_counts.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: gene_counts.png")


def plot_change_categories(df_comp):
    """Stacked bar chart of annotation change categories."""
    if df_comp.empty:
        print("  WARNING: No comparison data, skipping change_categories plot")
        return

    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(SPECIES_SHORT))
    width = 0.6

    categories = ['concordant', 'removed', 'added', 'split', 'fusion', 'complex']
    labels = ['Concordant (1→1)', 'Removed (1→0)', 'Added (0→1)',
              'Split (1→many)', 'Fusion (many→1)', 'Complex (many→many)']
    colors = ['#2ECC71', '#E74C3C', '#3498DB', '#F39C12', '#9B59B6', '#95A5A6']

    bottom = np.zeros(len(SPECIES_SHORT))
    for cat, label, color in zip(categories, labels, colors):
        vals = []
        for sp_file in SPECIES_IDS:
            row = df_comp[df_comp['species'] == sp_file]
            if not row.empty:
                vals.append(row.iloc[0].get(cat, 0))
            else:
                vals.append(0)
        ax.bar(x, vals, width, bottom=bottom, label=label, color=color, alpha=0.85)
        bottom += np.array(vals)

    ax.set_ylabel('Number of genes')
    ax.set_title('Annotation change categories after manual correction')
    ax.set_xticks(x)
    ax.set_xticklabels(SPECIES_SHORT, rotation=30, ha='right', fontsize=9)
    ax.legend(loc='upper right', fontsize=9)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{v/1000:.0f}K'))

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "change_categories.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: change_categories.png")


def plot_gene_length_distributions(df_stats):
    """Overlayed histograms of gene length for each species."""
    n_species = len(SPECIES_SHORT)
    n_cols = 3
    n_rows = (n_species + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 5 * n_rows))
    axes = axes.flatten()

    for i, (sp_label, sp_file) in enumerate(zip(SPECIES_SHORT, SPECIES_IDS)):
        ax = axes[i]
        row = df_stats[df_stats['species'] == sp_file]
        if row.empty:
            continue
        row = row.iloc[0]

        # Use median and 90th percentile as summary
        metrics_plot = [
            ('mean gene\nlength', 'mean gene length (bp)_'),
            ('median gene\nlength', 'median gene length (bp)_'),
            ('90%ile gene\nlength', '90 percentile gene length (bp)_'),
        ]
        before_vals = [row.get(m + 'before', 0) for _, m in metrics_plot]
        after_vals = [row.get(m + 'after', 0) for _, m in metrics_plot]

        x_pos = np.arange(len(metrics_plot))
        width = 0.3
        ax.bar(x_pos - width / 2, before_vals, width, label='Before', color='#E74C3C', alpha=0.85)
        ax.bar(x_pos + width / 2, after_vals, width, label='After', color='#2ECC71', alpha=0.85)
        ax.set_xticks(x_pos)
        ax.set_xticklabels([m[0] for m in metrics_plot], fontsize=8)
        ax.set_title(sp_label, fontsize=11)
        if i == 0:
            ax.legend(fontsize=8)
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{v/1000:.0f}k'))

    # Hide unused axes
    for j in range(n_species, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle('Gene length distribution metrics (before vs after)', fontsize=14, y=1.01)
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "gene_lengths.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: gene_lengths.png")


def plot_exon_counts(df_stats):
    """Exon count and monoexonic gene comparison."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    x = np.arange(len(SPECIES_SHORT))
    width = 0.35

    # Panel A: mean exons per mRNA
    ax = axes[0]
    before_exons = []
    after_exons = []
    for sp_file in SPECIES_IDS:
        row = df_stats[df_stats['species'] == sp_file]
        if not row.empty:
            row = row.iloc[0]
            before_exons.append(row.get('mean exons per mrna_before', 0))
            after_exons.append(row.get('mean exons per mrna_after', 0))
        else:
            before_exons.append(0)
            after_exons.append(0)

    ax.bar(x - width / 2, before_exons, width, label='Before', color='#E74C3C', alpha=0.85)
    ax.bar(x + width / 2, after_exons, width, label='After', color='#2ECC71', alpha=0.85)
    ax.set_ylabel('Mean exons per mRNA')
    ax.set_title('A. Mean exon count per transcript')
    ax.set_xticks(x)
    ax.set_xticklabels(SPECIES_SHORT, rotation=30, ha='right', fontsize=9)
    ax.legend()

    # Panel B: Single exon gene ratio
    ax = axes[1]
    before_ratio = []
    after_ratio = []
    for sp_file in SPECIES_IDS:
        row = df_stats[df_stats['species'] == sp_file]
        if not row.empty:
            row = row.iloc[0]
            b_single = row.get('Number of single exon gene_before', 0)
            b_total = row.get('Number of gene_before', 1)
            a_single = row.get('Number of single exon gene_after', 0)
            a_total = row.get('Number of gene_after', 1)
            before_ratio.append(b_single / b_total * 100 if b_total else 0)
            after_ratio.append(a_single / a_total * 100 if a_total else 0)
        else:
            before_ratio.append(0)
            after_ratio.append(0)

    ax.bar(x - width / 2, before_ratio, width, label='Before', color='#E74C3C', alpha=0.85)
    ax.bar(x + width / 2, after_ratio, width, label='After', color='#2ECC71', alpha=0.85)
    ax.set_ylabel('Single-exon genes (%)')
    ax.set_title('B. Single-exon gene proportion')
    ax.set_xticks(x)
    ax.set_xticklabels(SPECIES_SHORT, rotation=30, ha='right', fontsize=9)
    ax.legend()

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "exon_stats.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: exon_stats.png")


def plot_gffcompare_class_codes(df_gc):
    """Stacked bar of gffcompare transcript class codes."""
    if df_gc.empty:
        print("  WARNING: No gffcompare data, skipping class_codes plot")
        return

    class_cols = [c for c in df_gc.columns if c.startswith('class_')]
    if not class_cols:
        print("  WARNING: No class code columns found")
        return

    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(SPECIES_SHORT))
    width = 0.7

    # Color scheme for class codes
    code_colors = {
        'class_=': '#2ECC71',     # complete match - green
        'class_c': '#27AE60',     # contained
        'class_j': '#3498DB',     # multi-exon with junction match
        'class_k': '#2980B9',     # contains reference
        'class_m': '#9B59B6',     # retained intron
        'class_n': '#E67E22',     # novel intron
        'class_o': '#F39C12',     # other overlap
        'class_e': '#D35400',     # single exon covering intron
        'class_i': '#95A5A6',     # fully in intron
        'class_u': '#E74C3C',     # intergenic (novel)
        'class_p': '#C0392B',     # polymerase run-on
        'class_r': '#7F8C8D',     # repeat
        'class_x': '#34495E',     # opposite strand
        'class_s': '#BDC3C7',     # opposite strand
        'class_y': '#1ABC9C',     # contains reference in intron
    }

    # Group into major categories
    major_groups = {
        'Exact match (=)': ['class_='],
        'Contained/overlap (c,k,m,n,j,o,e)': ['class_c', 'class_k', 'class_m', 'class_n', 'class_j', 'class_o', 'class_e'],
        'Novel/Intergenic (u,i,y,p)': ['class_u', 'class_i', 'class_y', 'class_p'],
        'Other (r,s,x)': ['class_r', 'class_s', 'class_x'],
    }
    group_colors = ['#2ECC71', '#3498DB', '#E74C3C', '#95A5A6']

    bottom = np.zeros(len(SPECIES_SHORT))
    for (group_name, cols), color in zip(major_groups.items(), group_colors):
        vals = []
        for sp_file in SPECIES_IDS:
            row = df_gc[df_gc['species'] == sp_file]
            if not row.empty:
                total = sum(row.iloc[0].get(c, 0) for c in cols if c in df_gc.columns)
                vals.append(total)
            else:
                vals.append(0)
        ax.bar(x, vals, width, bottom=bottom, label=group_name, color=color, alpha=0.85)
        bottom += np.array(vals)

    ax.set_ylabel('Number of transcripts')
    ax.set_title('gffcompare transcript classification (after vs before)')
    ax.set_xticks(x)
    ax.set_xticklabels(SPECIES_SHORT, rotation=30, ha='right', fontsize=9)
    ax.legend(loc='upper right', fontsize=9)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{v/1000:.0f}K'))

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "class_codes.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: class_codes.png")


def plot_summary_heatmap(df_stats):
    """Heatmap showing key changes across all species."""
    metrics = [
        ('Number of gene', 'Gene count'),
        ('Number of mrna', 'mRNA count'),
        ('mean gene length (bp)', 'Mean gene len'),
        ('mean cds length (bp)', 'Mean CDS len'),
        ('mean exons per mrna', 'Mean exons/mRNA'),
        ('Number of single exon gene', 'Single-exon genes'),
    ]

    data = []
    row_labels = []
    for sp_file, sp_label in zip(SPECIES_IDS, SPECIES_SHORT):
        row = df_stats[df_stats['species'] == sp_file]
        if row.empty:
            continue
        row = row.iloc[0]
        row_data = []
        for metric, _ in metrics:
            b = row.get(f'{metric}_before', 0)
            a = row.get(f'{metric}_after', 0)
            if b and b != 0:
                pct_change = (a - b) / b * 100
            else:
                pct_change = 0
            row_data.append(pct_change)
        data.append(row_data)
        row_labels.append(sp_label)

    data = np.array(data)

    fig, ax = plt.subplots(figsize=(12, 8))
    cmap = sns.diverging_palette(10, 130, as_cmap=True)
    im = ax.imshow(data, cmap=cmap, aspect='auto', vmin=-50, vmax=50)

    ax.set_xticks(np.arange(len(metrics)))
    ax.set_xticklabels([m[1] for m in metrics], rotation=30, ha='right', fontsize=9)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=10)

    # Annotate cells
    for i in range(len(row_labels)):
        for j in range(len(metrics)):
            val = data[i, j]
            text = f'{val:+.1f}%'
            color = 'white' if abs(val) > 25 else 'black'
            ax.text(j, i, text, ha='center', va='center', fontsize=8, color=color)

    ax.set_title('Percentage change after manual correction (% change)', fontsize=13)
    fig.colorbar(im, ax=ax, label='% Change', shrink=0.8)
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "summary_heatmap.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: summary_heatmap.png")


def main():
    print("=" * 60)
    print("Loading aggregated data...")
    print("=" * 60)

    # Load data
    stats_df = pd.read_csv(RESULTS_DIR / "summary_stats.csv") if (RESULTS_DIR / "summary_stats.csv").exists() else pd.DataFrame()
    comp_df = pd.read_csv(RESULTS_DIR / "comparison_matrix.csv") if (RESULTS_DIR / "comparison_matrix.csv").exists() else pd.DataFrame()
    gc_df = pd.read_csv(RESULTS_DIR / "accuracy_metrics.csv") if (RESULTS_DIR / "accuracy_metrics.csv").exists() else pd.DataFrame()

    if stats_df.empty:
        print("ERROR: summary_stats.csv not found. Run summarize.py first.")
        return

    print(f"  Stats: {len(stats_df)} species, {len(stats_df.columns)} columns")
    print(f"  Comparison: {len(comp_df)} species" if not comp_df.empty else "  Comparison: not available")
    print(f"  gffcompare: {len(gc_df)} species" if not gc_df.empty else "  gffcompare: not available")

    print(f"\n{'=' * 60}")
    print("Generating figures...")
    print(f"{'=' * 60}")

    plot_gene_counts(stats_df, comp_df)
    plot_change_categories(comp_df)
    plot_gene_length_distributions(stats_df)
    plot_exon_counts(stats_df)
    plot_gffcompare_class_codes(gc_df)
    plot_summary_heatmap(stats_df)

    print(f"\n{'=' * 60}")
    print(f"All figures saved to {FIGURES_DIR}/")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
