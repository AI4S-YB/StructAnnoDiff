"""Shared project configuration for annotation comparison analyses."""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parent
ANALYSIS = Path(os.environ.get("ANALYSIS_DIR", ROOT_DIR)).resolve()
STATS_DIR = ANALYSIS / "stats"
COMPARE_DIR = ANALYSIS / "compare"
TCOMPARE_DIR = ANALYSIS / "tcompare"
RESULTS_DIR = ANALYSIS / "results"
FIGURES_DIR = ANALYSIS / "figures"

ANNOTATION_SUFFIXES = (".gff", ".gff3", ".gff.gz", ".gff3.gz")


@dataclass(frozen=True)
class Species:
    """Species metadata used by tables, figures, and batch runners."""

    id: str
    label: str
    short_label: str


def _load_species() -> list[Species]:
    with open(ROOT_DIR / "species.json", encoding="utf-8") as fh:
        rows = json.load(fh)
    return [Species(**row) for row in rows]


SPECIES = _load_species()
SPECIES_IDS = [sp.id for sp in SPECIES]
SPECIES_LABELS = {sp.id: sp.label for sp in SPECIES}
SPECIES_SHORT_LABELS = {sp.id: sp.short_label for sp in SPECIES}


def is_annotation_path(path: Path) -> bool:
    """Return True for primary GFF/GFF3 annotation inputs."""
    return path.name.endswith(ANNOTATION_SUFFIXES)


def find_annotation(species_id: str, state: str, analysis_dir: Path = ANALYSIS) -> Path | None:
    """Find the before/after annotation file for a species.

    This intentionally ignores derived files such as .tmap and .refmap.
    """
    analysis_dir = Path(analysis_dir)
    candidates = [
        path
        for path in analysis_dir.glob(f"{species_id}.{state}.*")
        if path.is_file() and is_annotation_path(path)
    ]
    return sorted(candidates)[0] if candidates else None
