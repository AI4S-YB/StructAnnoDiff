# Annotation Curation Comparison Analysis

This directory compares gene structure annotations before and after manual
curation for the configured species in `species.json`.

The GitHub repository is intended to contain the analysis code, tests,
documentation, lightweight summary tables, and publication figures. Large raw
FASTA/GFF/GTF inputs and large per-gene logs are intentionally excluded by
`.gitignore`.

## Inputs

Each species is expected to have one `before` and one `after` annotation file in
the analysis directory:

```text
<species_id>.before.gff
<species_id>.before.gff3
<species_id>.after.gff
<species_id>.after.gff3
```

Compressed `.gff.gz` and `.gff3.gz` files are also supported by the Python
parsers. Derived files such as `.tmap` and `.refmap` are ignored when resolving
primary annotations.

Raw input files are not committed to GitHub. To reproduce the analysis, either
place the input annotations in the project directory using the names above, or
set `ANALYSIS_DIR=/path/to/data` before running the commands below.

## Environment

Install the locked environment with:

```bash
pixi install
```

The main dependencies are Python, pandas, numpy, matplotlib, seaborn, AGAT,
gffcompare, and bedtools.

## Reproducible Workflow

Run the external comparison tools:

```bash
pixi run analyze
```

Build summary tables:

```bash
pixi run summarize
```

Run coordinate-based locus comparisons:

```bash
pixi run locus
```

The default locus scope is `mrna`, which excludes gene features without an
`mRNA` or `transcript` child. Use `python run_locus_comparisons.py --gene-scope all`
when non-mRNA gene features such as tRNA or ncRNA should be included.
The default overlap mode is strict `reciprocal`, requiring overlap to cover at
least the threshold fraction of both before and after gene spans. Use
`--overlap-mode containment` only to reproduce the legacy containment-style
matching behavior.

Generate final locus summary tables and figures:

```bash
pixi run tables
pixi run figures
pixi run target-figures
```

Generate A/B/C/D single-species summary tables and a four-panel figure:

```bash
python plot_single_species_abcd.py --species Pineapple
```

For table and figure regeneration from existing AGAT/gffcompare/locus outputs:

```bash
pixi run report
```

Run the complete workflow, including external tools and locus comparisons:

```bash
pixi run full
```

Validate cross-table consistency:

```bash
pixi run validate
```

## Outputs

- `stats/`: AGAT per-annotation statistics, generated locally and not tracked.
- `compare/`: AGAT before/after comparison reports, generated locally and not tracked.
- `tcompare/`: gffcompare outputs, generated locally and not tracked.
- `results/`: aggregate CSV/TSV tables and locus comparison logs.
- `figures/`: generated publication figures.
- `logs/`: full command logs from `run_analyses.sh`.

Important result tables:

- `summary_stats.csv`: coding-section AGAT statistics used by plots.
- `summary_stats_by_section.csv`: long-form AGAT statistics preserving all sections.
- `comparison_matrix.csv`: AGAT mRNA/transcript-scope comparison using `gene@mrna@cds`, `gene@mrna@exon`, `gene@transcript@cds`, and `gene@transcript@exon`.
- `comparison_by_feature_path.csv`: AGAT comparison split by every feature path table.
- `comparison_matrix_all_gene_types.csv`: broad AGAT comparison across all feature paths.
- `locus_comparison_summary.csv`: mutually exclusive locus subtype summary; one syntenic gene contributes to one broad category.
- `locus_comparison_multilabel.csv`: non-exclusive locus subtype attributes; one syntenic gene can count in multiple columns.
- `locus_diagnostics.csv`: overlap-mode diagnostics, including candidate pairs and containment-style pairs filtered by strict reciprocal overlap.
- `validation_report.csv` / `.md`: consistency checks across summary, compare, and locus outputs.

Large per-gene files such as `results/locus/*_change_log.csv` and
`results/single_species/*_syntenic_pair_deltas.csv` are generated locally but
excluded from GitHub. The tracked `*_change_summary.csv`, A/B/C/D summary
tables, and figures are sufficient for quick review.

Primary figures:

- `figure1_quantity_changes.png`: before/after quantity changes, including gene counts.
- `figure2_syntenic_structure_changes.png`: non-exclusive structural attributes for confirmed one-to-one gene pairs.
- `{Species}_ABCD_single_species.png`: per-species four-panel summary covering global quantities, locus fate, 1:1 structural attributes, and paired change magnitude.

## Java implementation

A JDK 11 / NetBeans Ant reimplementation lives in `java-locus-compare/`.
Generated Java outputs are excluded from GitHub; regenerate and validate them
locally with:

```bash
cd java-locus-compare
ant test
ant run-all
ant validate-python-parity
```

## Configuration

- Edit `species.json` to add, remove, or reorder species.
- Set `ANALYSIS_DIR=/path/to/analysis` or pass `--analysis-dir DIR` to run
  against another directory.
- Pass one or more species IDs to process a subset:

```bash
bash run_analyses.sh Rice Peach
python run_locus_comparisons.py Rice Peach
```

## Validation

Run the lightweight parser tests:

```bash
pixi run test
```

The key consistency expectation is that all configured species appear in
`summary_stats.csv`, `comparison_matrix.csv`, `accuracy_metrics.csv`, locus
summary tables, and the generated figures.
