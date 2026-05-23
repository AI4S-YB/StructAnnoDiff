import unittest

import pandas as pd

from generate_tables import (
    aggregate_multilabel_subtypes,
    aggregate_subtypes,
    build_multilabel_table,
    classify_exclusive_subtype,
)


class GenerateTablesTests(unittest.TestCase):
    def test_exclusive_subtype_priority_is_explicit(self):
        self.assertEqual(classify_exclusive_subtype("exact"), "exact")
        self.assertEqual(
            classify_exclusive_subtype("exon_gain_cds_extended_isoform_mrna_1x2_restructured"),
            "isoform",
        )
        self.assertEqual(
            classify_exclusive_subtype("exon_boundary_refined_cds_boundary_refined"),
            "cds_boundary_refined",
        )

    def test_unknown_exclusive_subtype_raises(self):
        row = pd.Series({
            "syntenic_total": 1,
            "syntenic_unexpected_change": 1,
        })

        with self.assertRaisesRegex(ValueError, "Unclassified syntenic subtype"):
            aggregate_subtypes(row)

    def test_unknown_multilabel_subtype_raises(self):
        row = pd.Series({
            "syntenic_total": 1,
            "syntenic_unexpected_change": 1,
        })

        with self.assertRaisesRegex(ValueError, "Unclassified syntenic subtype"):
            aggregate_multilabel_subtypes(row)

    def test_multilabel_boundary_has_explicit_any_bucket(self):
        row = pd.Series({
            "syntenic_total": 1,
            "syntenic_boundary_refined": 1,
        })

        cats = aggregate_multilabel_subtypes(row)

        self.assertEqual(cats["Boundary_refined"], 1)
        self.assertEqual(cats["Any_gene_boundary_changed"], 1)

    def test_multilabel_table_requires_direct_counts(self):
        import generate_tables

        original_species_order = generate_tables.SPECIES_ORDER
        original_species_short = generate_tables.SPECIES_SHORT
        original_locus_dir = generate_tables.LOCUS_DIR
        try:
            from tempfile import TemporaryDirectory
            from pathlib import Path

            with TemporaryDirectory() as tmpdir:
                locus_dir = Path(tmpdir)
                pd.DataFrame([{
                    "syntenic_total": 1,
                    "syntenic_exact": 1,
                }]).to_csv(locus_dir / "Test_change_summary.csv", index=False)
                generate_tables.SPECIES_ORDER = ["Test"]
                generate_tables.SPECIES_SHORT = {"Test": "Test"}
                generate_tables.LOCUS_DIR = locus_dir

                with self.assertRaisesRegex(ValueError, "missing direct multilabel"):
                    build_multilabel_table()
        finally:
            generate_tables.SPECIES_ORDER = original_species_order
            generate_tables.SPECIES_SHORT = original_species_short
            generate_tables.LOCUS_DIR = original_locus_dir


if __name__ == "__main__":
    unittest.main()
