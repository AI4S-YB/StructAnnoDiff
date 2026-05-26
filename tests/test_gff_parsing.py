import tempfile
import unittest
from pathlib import Path

from gff_utils import build_gene_index, normalize_attribute_id, normalize_attribute_ids
from locus_compare import parse_gff3_to_models


class GffParsingTests(unittest.TestCase):
    def test_normalize_attribute_id(self):
        self.assertEqual(normalize_attribute_id("gene:Gene1", prefixes=("gene:",)), "Gene1")
        self.assertEqual(normalize_attribute_id("transcript:Tx1", prefixes=("transcript:",)), "Tx1")
        self.assertEqual(normalize_attribute_id(["gene:Gene1"], prefixes=("gene:",)), "Gene1")
        self.assertEqual(
            normalize_attribute_ids("transcript:Tx1,transcript:Tx2", prefixes=("transcript:",)),
            ["Tx1", "Tx2"],
        )

    def test_build_gene_index_cleans_parent_prefixes(self):
        gff = "\n".join([
            "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=Gene1",
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=Tx1;Parent=gene:Gene1",
            "chr1\tsrc\texon\t1\t50\t.\t+\t.\tParent=transcript:Tx1",
            "chr1\tsrc\tCDS\t10\t40\t.\t+\t0\tParent=transcript:Tx1",
            "",
        ])
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test.gff3"
            path.write_text(gff, encoding="utf-8")
            genes = build_gene_index(path)

        self.assertIn("Gene1", genes)
        self.assertIn("Tx1", genes["Gene1"]["mrnas"])
        self.assertEqual(len(genes["Gene1"]["mrnas"]["Tx1"]["exons"]), 1)
        self.assertEqual(len(genes["Gene1"]["mrnas"]["Tx1"]["cdss"]), 1)

    def test_locus_models_use_shared_parser_and_prefix_cleanup(self):
        gff = "\n".join([
            "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene:Gene1",
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=transcript:Tx1;Parent=gene:Gene1",
            "chr1\tsrc\texon\t1\t100\t.\t+\t.\tParent=transcript:Tx1",
            "chr1\tsrc\tCDS\t10\t90\t.\t+\t0\tParent=transcript:Tx1",
            "",
        ])
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test.gff3"
            path.write_text(gff, encoding="utf-8")
            genes = parse_gff3_to_models(path)

        self.assertIn("Gene1", genes)
        self.assertEqual(genes["Gene1"].mrnas[0].mrna_id, "Tx1")
        self.assertEqual(genes["Gene1"].total_exon_count, 1)
        self.assertEqual(genes["Gene1"].total_cds_length, 81)

    def test_locus_models_exclude_genes_without_exon_or_cds_evidence(self):
        gff = "\n".join([
            "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=Gene1",
            "chr1\tsrc\tgene\t200\t300\t.\t+\t.\tID=Gene2",
            "chr1\tsrc\tgene\t400\t500\t.\t+\t.\tID=Gene3",
            "chr1\tsrc\tgene\t600\t700\t.\t+\t.\tID=Gene4",
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=Tx1;Parent=Gene1",
            "chr1\tsrc\tmRNA\t200\t300\t.\t+\t.\tID=Tx2;Parent=Gene2",
            "chr1\tsrc\tmRNA\t400\t500\t.\t+\t.\tID=Tx3;Parent=Gene3",
            "chr1\tsrc\tmRNA\t600\t700\t.\t+\t.\tID=Tx4;Parent=Gene4",
            "chr1\tsrc\texon\t1\t100\t.\t+\t.\tParent=Tx1",
            "chr1\tsrc\tfive_prime_UTR\t400\t450\t.\t+\t.\tParent=Tx3",
            "chr1\tsrc\tCDS\t620\t680\t.\t+\t0\tParent=Tx4",
            "",
        ])
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test.gff3"
            path.write_text(gff, encoding="utf-8")
            mrna_scope = parse_gff3_to_models(path)
            all_scope = parse_gff3_to_models(path, gene_scope="all")

        self.assertEqual(set(mrna_scope), {"Gene1", "Gene4"})
        self.assertEqual(set(all_scope), {"Gene1", "Gene4"})
        self.assertEqual(mrna_scope["Gene4"].mrnas[0].exon_count, 1)
        self.assertTrue(mrna_scope["Gene4"].mrnas[0].exons_are_inferred)

    def test_invalid_transcripts_are_removed_before_gene_span_is_recomputed(self):
        gff = "\n".join([
            "chr1\tsrc\tgene\t1\t1000\t.\t+\t.\tID=Gene1",
            "chr1\tsrc\tmRNA\t1\t1000\t.\t+\t.\tID=BadTx;Parent=Gene1",
            "chr1\tsrc\tmRNA\t500\t700\t.\t+\t.\tID=GoodTx;Parent=Gene1",
            "chr1\tsrc\texon\t500\t700\t.\t+\t.\tParent=GoodTx",
            "",
        ])
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test.gff3"
            path.write_text(gff, encoding="utf-8")
            genes = parse_gff3_to_models(path)

        self.assertEqual(set(genes), {"Gene1"})
        self.assertEqual([m.mrna_id for m in genes["Gene1"].mrnas], ["GoodTx"])
        self.assertEqual((genes["Gene1"].start, genes["Gene1"].end), (500, 700))

    def test_multi_parent_children_are_assigned_to_all_transcripts(self):
        gff = "\n".join([
            "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=Gene1",
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=Tx1;Parent=Gene1",
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=Tx2;Parent=Gene1",
            "chr1\tsrc\texon\t1\t100\t.\t+\t.\tParent=Tx1,Tx2",
            "chr1\tsrc\tCDS\t10\t90\t.\t+\t0\tParent=Tx1,Tx2",
            "",
        ])
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test.gff3"
            path.write_text(gff, encoding="utf-8")
            index = build_gene_index(path)
            models = parse_gff3_to_models(path)

        self.assertEqual(len(index["Gene1"]["mrnas"]["Tx1"]["exons"]), 1)
        self.assertEqual(len(index["Gene1"]["mrnas"]["Tx2"]["exons"]), 1)
        self.assertEqual(len(index["Gene1"]["mrnas"]["Tx1"]["cdss"]), 1)
        self.assertEqual(len(index["Gene1"]["mrnas"]["Tx2"]["cdss"]), 1)

        by_id = {m.mrna_id: m for m in models["Gene1"].mrnas}
        self.assertEqual(by_id["Tx1"].exon_count, 1)
        self.assertEqual(by_id["Tx2"].exon_count, 1)
        self.assertEqual(by_id["Tx1"].cds_length, 81)
        self.assertEqual(by_id["Tx2"].cds_length, 81)


if __name__ == "__main__":
    unittest.main()
