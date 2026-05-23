# Java Locus Compare

JDK 11 / NetBeans Ant implementation of the current locus-based annotation comparison algorithm.

## Commands

```bash
ant clean test
ant run-all
ant validate-python-parity
```

`ant run-all` reads `../species.json`, processes only the configured seven species, and writes Java outputs under `java-locus-compare/results/`.

`ant validate-python-parity` compares Java outputs against the Python baseline under `../results/`:

- per-species `*_change_summary.csv` fields
- per-species `change_log` match-type/subtype counts
- `locus_comparison_summary.csv`
- `locus_comparison_multilabel.csv`
- `locus_diagnostics.csv`

The Java implementation does not overwrite the Python result directory and does not generate figures; use the existing Python plotting workflow for PNG figures.
