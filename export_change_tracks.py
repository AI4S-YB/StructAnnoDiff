#!/usr/bin/env python3
"""Export locus change logs as IGV-friendly BED and GFF3 tracks."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any
from urllib.parse import quote

import pandas as pd

from analysis_config import RESULTS_DIR, SPECIES_IDS


MATCH_TYPE_COLORS = {
    "syntenic": "65,105,225",
    "split": "245,133,24",
    "merge": "128,85,170",
    "complex": "214,39,40",
    "novel": "44,160,44",
    "deleted": "120,120,120",
    "unresolved_overlap_before": "150,95,50",
    "unresolved_overlap_after": "150,95,50",
}

DEFAULT_COLOR = "80,80,80"

DERIVED_TSV_COLUMNS = [
    "event_id",
    "event_group_id",
    "species_id",
    "seqid",
    "event_start",
    "event_end",
    "bed_start",
    "bed_end",
    "strand",
    "color_rgb",
    "before_region",
    "after_region",
    "bed_name",
]


def clean_value(value: Any, missing: str = "") -> str:
    """Return a stable string representation for CSV/GFF/BED output."""
    if value is None:
        return missing
    try:
        if pd.isna(value):
            return missing
    except (TypeError, ValueError):
        pass
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def safe_token(value: Any, missing: str = "NA") -> str:
    """Return a compact token suitable for event IDs and display names."""
    text = clean_value(value, missing=missing) or missing
    return re.sub(r"[^A-Za-z0-9_.:-]+", "_", text)


def parse_int(value: Any) -> int | None:
    text = clean_value(value)
    if not text:
        return None
    try:
        number = int(float(text))
    except ValueError:
        return None
    return number


def interval_from_row(row: pd.Series, prefix: str) -> tuple[int, int] | None:
    start = parse_int(row.get(f"{prefix}_start"))
    end = parse_int(row.get(f"{prefix}_end"))
    if start is None or end is None or start <= 0 or end < start:
        return None
    return start, end


def event_interval(row: pd.Series) -> tuple[int, int]:
    """Return the union of valid before/after event coordinates."""
    intervals = [
        interval
        for interval in (
            interval_from_row(row, "before"),
            interval_from_row(row, "after"),
        )
        if interval is not None
    ]
    if not intervals:
        raise ValueError(f"Row has no valid before/after interval: {row.to_dict()}")
    return min(start for start, _end in intervals), max(end for _start, end in intervals)


def region_text(seqid: str, interval: tuple[int, int] | None) -> str:
    if interval is None:
        return "NA"
    return f"{seqid}:{interval[0]}-{interval[1]}"


def gff3_escape(value: Any) -> str:
    text = clean_value(value, missing="NA") or "NA"
    return quote(text, safe="A-Za-z0-9_.:|+-")


def gff3_attributes(attrs: dict[str, Any]) -> str:
    return ";".join(f"{key}={gff3_escape(value)}" for key, value in attrs.items())


def bed_name(row: pd.Series) -> str:
    match_type = safe_token(row.get("match_type"))
    subtype = safe_token(row.get("change_subtype"))
    before_gene = safe_token(row.get("before_gene"))
    after_gene = safe_token(row.get("after_gene"))
    return f"{match_type}|{subtype}|{before_gene}->{after_gene}"


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def find(self, item: str) -> str:
        if item not in self.parent:
            self.parent[item] = item
        if self.parent[item] != item:
            self.parent[item] = self.find(self.parent[item])
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def assign_event_group_ids(df: pd.DataFrame, species_id: str) -> list[str]:
    """Assign stable group IDs for event rows.

    Split rows are grouped by the before gene, merge rows by the after gene,
    and complex rows by connected components of the before/after overlap graph.
    """
    group_ids: dict[int, str] = {}
    complex_rows = df[df["match_type"].astype(str).eq("complex")]

    if not complex_rows.empty:
        uf = UnionFind()
        row_nodes: dict[int, tuple[str, str]] = {}
        for index, row in complex_rows.iterrows():
            before = f"b:{clean_value(row.get('seqid'))}:{clean_value(row.get('strand'))}:{safe_token(row.get('before_gene'))}"
            after = f"a:{clean_value(row.get('seqid'))}:{clean_value(row.get('strand'))}:{safe_token(row.get('after_gene'))}"
            uf.union(before, after)
            row_nodes[index] = (before, after)

        component_rows: dict[str, list[int]] = {}
        for index, (before, _after) in row_nodes.items():
            root = uf.find(before)
            component_rows.setdefault(root, []).append(index)

        sortable_components = []
        for root, indices in component_rows.items():
            rows = df.loc[indices]
            starts = [event_interval(row)[0] for _idx, row in rows.iterrows()]
            sortable_components.append((
                clean_value(rows.iloc[0].get("seqid")),
                min(starts),
                sorted(indices),
                root,
            ))

        for number, (_seqid, _start, indices, _root) in enumerate(sorted(sortable_components), start=1):
            group_id = f"{species_id}:complex:{number:05d}"
            for index in indices:
                group_ids[index] = group_id

    assigned = []
    for index, row in df.iterrows():
        match_type = clean_value(row.get("match_type"))
        before_gene = safe_token(row.get("before_gene"))
        after_gene = safe_token(row.get("after_gene"))
        if match_type == "split":
            group_id = f"{species_id}:split:{before_gene}"
        elif match_type == "merge":
            group_id = f"{species_id}:merge:{after_gene}"
        elif match_type == "complex":
            group_id = group_ids[index]
        elif match_type == "novel":
            group_id = f"{species_id}:novel:{after_gene}"
        elif match_type == "deleted":
            group_id = f"{species_id}:deleted:{before_gene}"
        elif match_type == "unresolved_overlap_before":
            group_id = f"{species_id}:unresolved_before:{before_gene}"
        elif match_type == "unresolved_overlap_after":
            group_id = f"{species_id}:unresolved_after:{after_gene}"
        else:
            group_id = f"{species_id}:{safe_token(match_type)}:{before_gene}:{after_gene}"
        assigned.append(group_id)
    return assigned


def build_track_record(row: pd.Series, species_id: str, event_id: str, event_group_id: str) -> dict[str, Any]:
    seqid = clean_value(row.get("seqid"), missing="NA") or "NA"
    strand = clean_value(row.get("strand"), missing=".") or "."
    start, end = event_interval(row)
    bed_start = start - 1
    bed_end = end
    before_interval = interval_from_row(row, "before")
    after_interval = interval_from_row(row, "after")
    match_type = clean_value(row.get("match_type"), missing="unknown") or "unknown"
    color = MATCH_TYPE_COLORS.get(match_type, DEFAULT_COLOR)

    return {
        "event_id": event_id,
        "event_group_id": event_group_id,
        "species_id": species_id,
        "seqid": seqid,
        "event_start": start,
        "event_end": end,
        "bed_start": bed_start,
        "bed_end": bed_end,
        "strand": strand,
        "color_rgb": color,
        "before_region": region_text(seqid, before_interval),
        "after_region": region_text(seqid, after_interval),
        "bed_name": bed_name(row),
    }


def bed_line(record: dict[str, Any]) -> str:
    return "\t".join([
        record["seqid"],
        str(record["bed_start"]),
        str(record["bed_end"]),
        record["bed_name"],
        "0",
        record["strand"],
        str(record["bed_start"]),
        str(record["bed_end"]),
        record["color_rgb"],
    ])


def gff3_line(row: pd.Series, record: dict[str, Any]) -> str:
    attrs = {
        "ID": record["event_id"],
        "Name": record["bed_name"],
        "event_group_id": record["event_group_id"],
        "species_id": record["species_id"],
        "before_gene": clean_value(row.get("before_gene"), missing="NA") or "NA",
        "after_gene": clean_value(row.get("after_gene"), missing="NA") or "NA",
        "match_type": clean_value(row.get("match_type"), missing="NA") or "NA",
        "change_subtype": clean_value(row.get("change_subtype"), missing="NA") or "NA",
        "before_region": record["before_region"],
        "after_region": record["after_region"],
        "before_exons": clean_value(row.get("before_exons"), missing="0") or "0",
        "after_exons": clean_value(row.get("after_exons"), missing="0") or "0",
        "before_cds": clean_value(row.get("before_cds"), missing="0") or "0",
        "after_cds": clean_value(row.get("after_cds"), missing="0") or "0",
        "before_mrnas": clean_value(row.get("before_mrnas"), missing="0") or "0",
        "after_mrnas": clean_value(row.get("after_mrnas"), missing="0") or "0",
        "color_rgb": record["color_rgb"],
    }
    return "\t".join([
        record["seqid"],
        "StructAnnoDiff",
        "annotation_change",
        str(record["event_start"]),
        str(record["event_end"]),
        ".",
        record["strand"],
        ".",
        gff3_attributes(attrs),
    ])


def export_species(species_id: str, results_dir: Path = RESULTS_DIR, output_dir: Path | None = None) -> dict[str, Any]:
    results_dir = Path(results_dir)
    output_dir = Path(output_dir) if output_dir is not None else results_dir / "tracks"
    input_path = results_dir / "locus" / f"{species_id}_change_log.csv"
    if not input_path.exists():
        raise FileNotFoundError(f"Missing change log: {input_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    bed_path = output_dir / f"{species_id}_annotation_changes.bed"
    gff3_path = output_dir / f"{species_id}_annotation_changes.gff3"
    tsv_path = output_dir / f"{species_id}_annotation_changes.tsv"

    df = pd.read_csv(input_path)
    if df.empty:
        raise ValueError(f"Change log is empty: {input_path}")
    group_ids = assign_event_group_ids(df, species_id)

    tsv_records = []
    with bed_path.open("w", encoding="utf-8") as bed_fh, gff3_path.open("w", encoding="utf-8") as gff_fh:
        gff_fh.write("##gff-version 3\n")
        for number, ((index, row), group_id) in enumerate(zip(df.iterrows(), group_ids), start=1):
            event_id = f"{species_id}_change_{number:06d}"
            record = build_track_record(row, species_id, event_id, group_id)
            bed_fh.write(bed_line(record) + "\n")
            gff_fh.write(gff3_line(row, record) + "\n")
            tsv_records.append({
                **record,
                **{column: clean_value(row.get(column)) for column in df.columns},
                "source_row_index": index,
            })

    tsv_columns = [
        *DERIVED_TSV_COLUMNS,
        "source_row_index",
        *[column for column in df.columns if column not in DERIVED_TSV_COLUMNS],
    ]
    pd.DataFrame(tsv_records).reindex(columns=tsv_columns).to_csv(tsv_path, sep="\t", index=False)

    return {
        "species_id": species_id,
        "rows": len(df),
        "bed": bed_path,
        "gff3": gff3_path,
        "tsv": tsv_path,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export locus change logs as BED/GFF3/TSV tracks for IGV review."
    )
    parser.add_argument(
        "--species",
        nargs="+",
        help="Species ID(s) to export. Defaults to all configured species when omitted.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Export all species configured in species.json.",
    )
    parser.add_argument(
        "--results-dir",
        default=str(RESULTS_DIR),
        help="Directory containing results/locus. Defaults to the project results directory.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output directory. Defaults to <results-dir>/tracks.",
    )
    return parser.parse_args()


def selected_species(args: argparse.Namespace) -> list[str]:
    if args.all or not args.species:
        return list(SPECIES_IDS)
    unknown = sorted(set(args.species) - set(SPECIES_IDS))
    if unknown:
        raise SystemExit(f"Unknown species ID(s): {', '.join(unknown)}")
    return args.species


def main() -> None:
    args = parse_args()
    results_dir = Path(args.results_dir).resolve()
    output_dir = Path(args.output).resolve() if args.output else results_dir / "tracks"
    for species_id in selected_species(args):
        result = export_species(species_id, results_dir=results_dir, output_dir=output_dir)
        print(
            f"{result['species_id']}: exported {result['rows']} events -> "
            f"{result['bed']}, {result['gff3']}, {result['tsv']}"
        )


if __name__ == "__main__":
    main()
