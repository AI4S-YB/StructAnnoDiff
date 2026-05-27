# Java Locus Compare

JDK 11 / NetBeans Ant implementation of the current locus-based annotation comparison algorithm.
The runtime path is pure Java: it parses GFF/GFF3 files, performs locus matching, writes CSV outputs, and does not call Python or plotting libraries.

## Commands

```bash
ant clean test
ant run-all
ant validate-python-parity
ant plot-pineapple
```

`ant run-all` reads `../species.json`, processes only the configured seven species, and writes Java outputs under `java-locus-compare/results/`.

Main outputs:

- `results/locus/<species>_change_summary.csv`: per-species summary counters.
- `results/locus/<species>_change_log.csv`: per matched or unmatched locus record.
- `results/locus_comparison_summary.csv`: publication-style broad subtype summary.
- `results/locus_comparison_multilabel.csv`: non-exclusive subtype counts.
- `results/locus_diagnostics.csv`: overlap candidate and pruning diagnostics.
- `results/curation_core_metrics.csv`: compact per-species metrics used by the final figure.
- `results/curation_core_panel_metrics.csv`: figure-equivalent table with one row per species/bar. Panel A is locus gain/loss, Panel B is split/merge events, and Panel C is representative-transcript exon-change percentage against before and after references.
- `results/curation_summary_table.csv`: wide, human-readable summary table with input file names, deleted/new loci, split/merge counts, and exon-change ratios formatted as `percent (count/total)`.
- `results/curation_summary_table_columns.csv`: column descriptions for `curation_summary_table.csv`.

`ant validate-python-parity` compares Java outputs against the Python baseline under `../results/`:

- per-species `*_change_summary.csv` fields
- per-species `change_log` match-type/subtype counts
- `locus_comparison_summary.csv`
- `locus_comparison_multilabel.csv`
- `locus_diagnostics.csv`
- `curation_core_metrics.csv`

The Java implementation does not overwrite the Python result directory and does not generate figures. Use `curation_summary_table.csv` when a compact pure-Java table equivalent to the figure is needed.

`ant plot-pineapple` creates `results/figures/Pineapple_curation_summary_jigplot.{png,pdf,svg}` from `curation_core_metrics.csv`. It loads JIGplot at runtime and does not modify TBtools. Override the JIGplot path when needed:

```bash
ant -Djigplot.jar=/path/to/JIGplot.jar plot-pineapple
```
