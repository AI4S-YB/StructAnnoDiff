import unittest
from collections import defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory

from locus_compare import (
    CDSFeature,
    ExonFeature,
    GeneModel,
    MRNAModel,
    UTRFeature,
    classify_syntenic_change,
    compare_annotations,
    compute_syntenic_attributes,
    containment_overlap,
    find_overlapping_pairs,
    parse_gff3_to_models,
    reciprocal_overlap,
    resolve_matches,
)


class LocusResolutionTests(unittest.TestCase):
    def _gene(self, gene_id, start, end, raw_start=None, raw_end=None):
        return GeneModel(
            gene_id=gene_id,
            seqid="chr1",
            start=start,
            end=end,
            strand="+",
            source="test",
            mrnas=[],
            raw_start=raw_start,
            raw_end=raw_end,
        )

    def _gene_with_mrnas(self, gene_id, mrnas):
        gene = self._gene(gene_id, 1, 300)
        gene.mrnas = mrnas
        return gene

    def _mrna(self, mrna_id, exons, cds, utrs=None):
        return MRNAModel(
            mrna_id=mrna_id,
            exons=[ExonFeature(start, end) for start, end in exons],
            cds=[CDSFeature(start, end, phase) for start, end, phase in cds],
            utrs=[UTRFeature(start, end, utr_type) for start, end, utr_type in (utrs or [])],
        )

    def test_before_after_same_ids_are_distinct_graph_nodes(self):
        before = {
            "A": self._gene("A", 1, 100),
            "B": self._gene("B", 201, 300),
        }
        after = {
            "A": self._gene("A", 201, 300),
            "B": self._gene("B", 1, 100),
        }
        pairs = [
            (before["A"], after["B"], 1.0, 1.0),
            (before["B"], after["A"], 1.0, 1.0),
        ]

        result = resolve_matches(pairs, before, after)

        self.assertEqual(len(result["syntenic"]), 2)
        self.assertEqual(len(result["complex"]), 0)
        self.assertEqual(len(result["novel"]), 0)
        self.assertEqual(len(result["deleted"]), 0)

    def test_reciprocal_overlap_rejects_simple_containment(self):
        before = self._gene("small", 101, 200)
        after = self._gene("large", 1, 1000)

        self.assertAlmostEqual(reciprocal_overlap(before, after), 0.1)
        self.assertAlmostEqual(containment_overlap(before, after), 1.0)

        strict_pairs = find_overlapping_pairs(
            [before],
            [after],
            min_reciprocal=0.5,
            overlap_mode="reciprocal",
        )
        legacy_pairs = find_overlapping_pairs(
            [before],
            [after],
            min_reciprocal=0.5,
            overlap_mode="containment",
        )

        self.assertEqual(strict_pairs, [])
        self.assertEqual(len(legacy_pairs), 1)

    def test_overlap_diagnostics_count_only_real_overlaps(self):
        before = [
            self._gene("non_overlap", 1, 10),
            self._gene("overlap", 1005, 1010),
        ]
        after = [self._gene("after", 1000, 1100)]
        diagnostics = defaultdict(int)

        pairs = find_overlapping_pairs(
            before,
            after,
            min_reciprocal=0.01,
            overlap_mode="reciprocal",
            diagnostics=diagnostics,
        )

        self.assertEqual(len(pairs), 1)
        self.assertEqual(diagnostics["same_strand_overlaps"], 1)

    def test_weak_rejected_does_not_depend_on_diagnostics(self):
        before = self._gene("small", 101, 200)
        after = self._gene("large", 1, 1000)
        weak_rejected = {"before": set(), "after": set()}

        pairs = find_overlapping_pairs(
            [before],
            [after],
            min_reciprocal=0.5,
            overlap_mode="reciprocal",
            weak_rejected=weak_rejected,
        )

        self.assertEqual(pairs, [])
        self.assertEqual(weak_rejected["before"], {"small"})
        self.assertEqual(weak_rejected["after"], {"large"})

    def test_primary_mrna_tie_uses_file_order(self):
        first = self._mrna("tx_a", [(1, 100)], [(1, 90, "0")])
        second = self._mrna("tx_z", [(200, 299)], [(200, 289, "0")])
        gene = self._gene_with_mrnas("gene", [first, second])

        self.assertIs(gene.primary_mrna, first)

    def test_exact_requires_cds_coordinate_identity(self):
        before = self._gene_with_mrnas(
            "before",
            [self._mrna("b1", [(1, 100)], [(10, 90, "0")])],
        )
        after = self._gene_with_mrnas(
            "after",
            [self._mrna("a1", [(1, 100)], [(11, 91, "0")])],
        )

        subtype = classify_syntenic_change(before, after)

        self.assertNotEqual(subtype, "exact")
        self.assertIn("cds_boundary_refined", subtype)

    def test_boundary_refinement_uses_raw_gene_boundaries(self):
        before = self._gene_with_mrnas(
            "before",
            [self._mrna("b1", [(1, 100)], [(10, 90, "0")])],
        )
        after = self._gene_with_mrnas(
            "after",
            [self._mrna("a1", [(1, 100)], [(10, 90, "0")])],
        )
        before.raw_start = 1
        before.raw_end = 100
        after.raw_start = 20
        after.raw_end = 100

        subtype = classify_syntenic_change(before, after, boundary_tol=10)

        self.assertEqual(subtype, "boundary_refined")

    def test_exact_requires_non_primary_isoform_coordinate_identity(self):
        before = self._gene_with_mrnas(
            "before",
            [
                self._mrna("b1", [(1, 100)], [(10, 90, "0")]),
                self._mrna("b2", [(1, 60), (80, 100)], [(10, 50, "0")]),
            ],
        )
        after = self._gene_with_mrnas(
            "after",
            [
                self._mrna("a1", [(1, 100)], [(10, 90, "0")]),
                self._mrna("a2", [(1, 55), (75, 100)], [(10, 50, "0")]),
            ],
        )

        subtype = classify_syntenic_change(before, after)

        self.assertNotEqual(subtype, "exact")
        self.assertIn("isoform_restructured", subtype)

    def test_utr_order_does_not_trigger_refinement(self):
        before = self._gene_with_mrnas(
            "before",
            [self._mrna(
                "b1",
                [(1, 100)],
                [(20, 80, "0")],
                [(1, 19, "five_prime_UTR"), (81, 100, "three_prime_UTR")],
            )],
        )
        after = self._gene_with_mrnas(
            "after",
            [self._mrna(
                "a1",
                [(1, 100)],
                [(20, 80, "0")],
                [(81, 100, "three_prime_UTR"), (1, 19, "five_prime_UTR")],
            )],
        )

        subtype = classify_syntenic_change(before, after)

        self.assertEqual(subtype, "exact")

    def test_after_transcript_is_matched_to_before_primary(self):
        before = self._gene_with_mrnas(
            "before",
            [self._mrna("b1", [(1, 100)], [(10, 90, "0")])],
        )
        after = self._gene_with_mrnas(
            "after",
            [
                self._mrna("a_extra", [(1, 100)], [(11, 91, "0")]),
                self._mrna("a_best", [(1, 100)], [(10, 90, "0")]),
            ],
        )

        subtype = classify_syntenic_change(before, after)

        self.assertIn("isoform_mrna_1x2_restructured", subtype)
        self.assertNotIn("cds_boundary_refined", subtype)

    def test_syntenic_attributes_keep_utr_refined_with_utr_gain(self):
        before = self._gene_with_mrnas(
            "before",
            [self._mrna(
                "b1",
                [(1, 100)],
                [(20, 80, "0")],
                [(1, 19, "five_prime_UTR")],
            )],
        )
        after = self._gene_with_mrnas(
            "after",
            [self._mrna(
                "a1",
                [(1, 120)],
                [(20, 80, "0")],
                [(1, 30, "five_prime_UTR"), (81, 120, "three_prime_UTR")],
            )],
        )

        attrs = compute_syntenic_attributes(before, after)

        self.assertTrue(attrs["utr_added"])
        self.assertTrue(attrs["utr_refined"])
        self.assertFalse(attrs["exact"])

    def test_syntenic_attributes_keep_utr_gain_with_utr_exon_gain(self):
        before = self._gene_with_mrnas(
            "before",
            [self._mrna("b1", [(20, 80)], [(20, 80, "0")])],
        )
        after = self._gene_with_mrnas(
            "after",
            [self._mrna(
                "a1",
                [(1, 10), (20, 80)],
                [(20, 80, "0")],
                [(1, 10, "five_prime_UTR")],
            )],
        )

        attrs = compute_syntenic_attributes(before, after)

        self.assertTrue(attrs["utr_exon_added"])
        self.assertTrue(attrs["utr_added"])
        self.assertFalse(attrs["coding_exon_gain"])

    def test_utr_exon_gain_is_not_coding_gain_when_cds_boundary_changes(self):
        before = self._gene_with_mrnas(
            "before",
            [self._mrna("b1", [(20, 80)], [(20, 80, "0")])],
        )
        after = self._gene_with_mrnas(
            "after",
            [self._mrna(
                "a1",
                [(1, 10), (20, 81)],
                [(20, 81, "0")],
                [(1, 10, "five_prime_UTR")],
            )],
        )

        attrs = compute_syntenic_attributes(before, after)
        subtype = classify_syntenic_change(before, after)

        self.assertTrue(attrs["utr_exon_added"])
        self.assertFalse(attrs["coding_exon_gain"])
        self.assertIn("utr_exon_added", subtype)
        self.assertIn("cds_extended", subtype)
        self.assertNotIn("exon_gain", subtype)

    def test_missing_utrs_are_derived_from_transcript_bounds_without_exons(self):
        with TemporaryDirectory() as tmpdir:
            gff = Path(tmpdir) / "input.gff3"
            gff.write_text(
                "\n".join([
                    "chr1\ttest\tgene\t1\t300\t.\t+\t.\tID=g1",
                    "chr1\ttest\tmRNA\t1\t300\t.\t+\t.\tID=t1;Parent=g1",
                    "chr1\ttest\tCDS\t101\t200\t.\t+\t0\tID=cds1;Parent=t1",
                ]),
                encoding="utf-8",
            )

            genes = parse_gff3_to_models(gff)

        mrna = genes["g1"].mrnas[0]
        self.assertEqual(
            [(u.start, u.end, u.utr_type) for u in mrna.utrs],
            [(1, 100, "five_prime_UTR"), (201, 300, "three_prime_UTR")],
        )
        self.assertEqual([(e.start, e.end) for e in mrna.exons], [(1, 300)])

    def test_missing_one_sided_utr_is_derived_from_exons(self):
        with TemporaryDirectory() as tmpdir:
            gff = Path(tmpdir) / "input.gff3"
            gff.write_text(
                "\n".join([
                    "chr1\ttest\tgene\t1\t300\t.\t+\t.\tID=g1",
                    "chr1\ttest\tmRNA\t1\t300\t.\t+\t.\tID=t1;Parent=g1",
                    "chr1\ttest\texon\t1\t300\t.\t+\t.\tID=e1;Parent=t1",
                    "chr1\ttest\tfive_prime_UTR\t1\t100\t.\t+\t.\tID=u1;Parent=t1",
                    "chr1\ttest\tCDS\t101\t200\t.\t+\t0\tID=cds1;Parent=t1",
                ]),
                encoding="utf-8",
            )

            genes = parse_gff3_to_models(gff)

        mrna = genes["g1"].mrnas[0]
        self.assertEqual(
            [(u.start, u.end, u.utr_type) for u in mrna.utrs],
            [(1, 100, "five_prime_UTR"), (201, 300, "three_prime_UTR")],
        )

    def test_no_common_seqids_returns_empty_result_sentinel(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            before = tmp / "before.gff3"
            after = tmp / "after.gff3"
            before.write_text(
                "\n".join([
                    "chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=g1",
                    "chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1",
                    "chr1\ttest\texon\t1\t100\t.\t+\t.\tParent=t1",
                    "chr1\ttest\tCDS\t10\t90\t.\t+\t0\tParent=t1",
                ]),
                encoding="utf-8",
            )
            after.write_text(
                "\n".join([
                    "chr2\ttest\tgene\t1\t100\t.\t+\t.\tID=g2",
                    "chr2\ttest\tmRNA\t1\t100\t.\t+\t.\tID=t2;Parent=g2",
                    "chr2\ttest\texon\t1\t100\t.\t+\t.\tParent=t2",
                    "chr2\ttest\tCDS\t10\t90\t.\t+\t0\tParent=t2",
                ]),
                encoding="utf-8",
            )

            summary, change_log = compare_annotations(before, after)

        self.assertIsNone(summary)
        self.assertIsNone(change_log)

    def test_weak_containment_overlap_is_not_counted_as_novel_deleted(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            before = tmp / "before.gff3"
            after = tmp / "after.gff3"
            before.write_text(
                "\n".join([
                    "chr1\ttest\tgene\t101\t200\t.\t+\t.\tID=g1",
                    "chr1\ttest\tmRNA\t101\t200\t.\t+\t.\tID=t1;Parent=g1",
                    "chr1\ttest\texon\t101\t200\t.\t+\t.\tParent=t1",
                    "chr1\ttest\tCDS\t101\t200\t.\t+\t0\tParent=t1",
                ]),
                encoding="utf-8",
            )
            after.write_text(
                "\n".join([
                    "chr1\ttest\tgene\t1\t1000\t.\t+\t.\tID=g2",
                    "chr1\ttest\tmRNA\t1\t1000\t.\t+\t.\tID=t2;Parent=g2",
                    "chr1\ttest\texon\t1\t1000\t.\t+\t.\tParent=t2",
                    "chr1\ttest\tCDS\t1\t1000\t.\t+\t0\tParent=t2",
                ]),
                encoding="utf-8",
            )

            summary, change_log = compare_annotations(before, after)

        self.assertEqual(summary["syntenic_total"], 0)
        self.assertEqual(summary["novel_genes"], 0)
        self.assertEqual(summary["deleted_genes"], 0)
        self.assertEqual(summary["unresolved_overlap_after_genes"], 1)
        self.assertEqual(summary["unresolved_overlap_before_genes"], 1)
        self.assertEqual(
            {row["match_type"] for row in change_log},
            {"unresolved_overlap_after", "unresolved_overlap_before"},
        )


if __name__ == "__main__":
    unittest.main()
