import csv
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from export_change_tracks import (
    assign_event_group_ids,
    build_track_record,
    event_interval,
    export_species,
    gff3_attributes,
)


def base_row(**overrides):
    row = {
        "before_gene": "before1",
        "after_gene": "after1",
        "seqid": "chr1",
        "before_start": 10,
        "before_end": 50,
        "after_start": 12,
        "after_end": 60,
        "strand": "+",
        "match_type": "syntenic",
        "change_subtype": "exon_boundary_refined",
        "before_exons": 2,
        "after_exons": 2,
        "before_cds": 1,
        "after_cds": 1,
        "before_mrnas": 1,
        "after_mrnas": 1,
    }
    row.update(overrides)
    return row


class ExportChangeTracksTests(unittest.TestCase):
    def test_event_interval_uses_available_side_for_novel_and_deleted(self):
        novel = pd.Series(base_row(
            before_gene="",
            before_start=0,
            before_end=0,
            after_start=100,
            after_end=200,
            match_type="novel",
        ))
        deleted = pd.Series(base_row(
            after_gene="",
            before_start=300,
            before_end=400,
            after_start=0,
            after_end=0,
            match_type="deleted",
        ))

        self.assertEqual(event_interval(novel), (100, 200))
        self.assertEqual(event_interval(deleted), (300, 400))

    def test_build_track_record_converts_to_bed_coordinates(self):
        row = pd.Series(base_row(before_start=1, before_end=10, after_start=3, after_end=20))

        record = build_track_record(row, "Test", "Test_change_000001", "Test:syntenic:before1:after1")

        self.assertEqual(record["event_start"], 1)
        self.assertEqual(record["event_end"], 20)
        self.assertEqual(record["bed_start"], 0)
        self.assertEqual(record["bed_end"], 20)
        self.assertEqual(record["before_region"], "chr1:1-10")
        self.assertEqual(record["after_region"], "chr1:3-20")

    def test_split_merge_and_complex_group_ids_are_stable(self):
        df = pd.DataFrame([
            base_row(before_gene="b1", after_gene="a1", match_type="split"),
            base_row(before_gene="b1", after_gene="a2", match_type="split"),
            base_row(before_gene="b2", after_gene="a3", match_type="merge"),
            base_row(before_gene="b3", after_gene="a3", match_type="merge"),
            base_row(before_gene="b4", after_gene="a4", match_type="complex", before_start=500, before_end=700),
            base_row(before_gene="b4", after_gene="a5", match_type="complex", before_start=500, before_end=700),
            base_row(before_gene="b5", after_gene="a5", match_type="complex", before_start=800, before_end=900),
        ])

        groups = assign_event_group_ids(df, "Test")

        self.assertEqual(groups[0], groups[1])
        self.assertEqual(groups[2], groups[3])
        self.assertEqual(groups[4], groups[5])
        self.assertEqual(groups[5], groups[6])
        self.assertTrue(groups[4].startswith("Test:complex:"))

    def test_gff3_attributes_escape_special_characters(self):
        attributes = gff3_attributes({
            "ID": "event1",
            "change_subtype": "boundary;with space",
            "before_gene": "gene=1",
        })

        self.assertIn("change_subtype=boundary%3Bwith%20space", attributes)
        self.assertIn("before_gene=gene%3D1", attributes)
        self.assertEqual(attributes.count(";"), 2)

    def test_export_species_writes_matching_bed_gff3_and_tsv_rows(self):
        with TemporaryDirectory() as tmpdir:
            results_dir = Path(tmpdir) / "results"
            locus_dir = results_dir / "locus"
            locus_dir.mkdir(parents=True)
            pd.DataFrame([
                base_row(match_type="syntenic"),
                base_row(
                    before_gene="",
                    before_start=0,
                    before_end=0,
                    after_gene="new1",
                    after_start=100,
                    after_end=200,
                    match_type="novel",
                    change_subtype="new_gene",
                ),
                base_row(
                    before_gene="old1",
                    before_start=300,
                    before_end=400,
                    after_gene="",
                    after_start=0,
                    after_end=0,
                    match_type="deleted",
                    change_subtype="lost_gene",
                ),
            ]).to_csv(locus_dir / "Test_change_log.csv", index=False)

            result = export_species("Test", results_dir=results_dir)

            bed_rows = result["bed"].read_text(encoding="utf-8").strip().splitlines()
            gff_rows = [
                line
                for line in result["gff3"].read_text(encoding="utf-8").splitlines()
                if not line.startswith("#")
            ]
            with result["tsv"].open(encoding="utf-8", newline="") as fh:
                tsv_rows = list(csv.DictReader(fh, delimiter="\t"))

        self.assertEqual(len(bed_rows), 3)
        self.assertEqual(len(gff_rows), 3)
        self.assertEqual(len(tsv_rows), 3)
        self.assertEqual(bed_rows[0].split("\t")[1:3], ["9", "60"])
        self.assertIn("match_type=novel", gff_rows[1])
        self.assertEqual(tsv_rows[2]["before_region"], "chr1:300-400")
        self.assertEqual(tsv_rows[2]["after_region"], "NA")


if __name__ == "__main__":
    unittest.main()
