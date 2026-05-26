#!/usr/bin/env python3
"""Batch-run locus_compare.py for configured species."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from analysis_config import ANALYSIS, ROOT_DIR, SPECIES_IDS, find_annotation


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run locus-based comparisons for all configured species.")
    parser.add_argument(
        "species",
        nargs="*",
        help="Optional species IDs to process. Defaults to all species in species.json.",
    )
    parser.add_argument(
        "--analysis-dir",
        default=str(ANALYSIS),
        help="Analysis directory containing annotation files and results/. Defaults to ANALYSIS_DIR or this directory.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output directory. Defaults to <analysis-dir>/results/locus.",
    )
    parser.add_argument(
        "--gene-scope",
        choices=("mrna", "coding", "all"),
        default="mrna",
        help="Gene models to compare. Defaults to mrna; use all to include non-mRNA gene features.",
    )
    parser.add_argument(
        "--overlap-mode",
        choices=("reciprocal", "containment", "hybrid"),
        default="hybrid",
        help="Overlap score for candidate matching. Defaults to hybrid: reciprocal for strict 1:1 and containment for split/merge events.",
    )
    parser.add_argument(
        "--reciprocal-overlap",
        type=float,
        default=0.5,
        help="Minimum reciprocal overlap for candidate matching. Defaults to 0.5.",
    )
    parser.add_argument(
        "--boundary-tol",
        type=int,
        default=10,
        help="Boundary tolerance in bp for exact gene-boundary matching. Defaults to 10.",
    )
    parser.add_argument(
        "--cds-change-pct",
        type=float,
        default=0.1,
        help="CDS length change threshold for refinement. Defaults to 0.1.",
    )
    parser.add_argument(
        "--utr-change-pct",
        type=float,
        default=0.1,
        help="UTR length change threshold for refinement. Defaults to 0.1.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    analysis_dir = Path(args.analysis_dir).resolve()
    output_dir = Path(args.output).resolve() if args.output else analysis_dir / "results" / "locus"
    output_dir.mkdir(parents=True, exist_ok=True)

    species_ids = args.species or SPECIES_IDS
    unknown = sorted(set(species_ids) - set(SPECIES_IDS))
    if unknown:
        raise SystemExit(f"Unknown species ID(s): {', '.join(unknown)}")

    for species_id in species_ids:
        before = find_annotation(species_id, "before", analysis_dir)
        after = find_annotation(species_id, "after", analysis_dir)
        if before is None or after is None:
            raise SystemExit(f"Missing before/after annotation for {species_id}")

        print(f"=== {species_id} ===", flush=True)
        cmd = [
            sys.executable,
            str(ROOT_DIR / "locus_compare.py"),
            "--before",
            str(before),
            "--after",
            str(after),
            "--output",
            str(output_dir),
            "--name",
            species_id,
            "--gene-scope",
            args.gene_scope,
            "--overlap-mode",
            args.overlap_mode,
            "--reciprocal-overlap",
            str(args.reciprocal_overlap),
            "--boundary-tol",
            str(args.boundary_tol),
            "--cds-change-pct",
            str(args.cds_change_pct),
            "--utr-change-pct",
            str(args.utr_change_pct),
        ]
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
