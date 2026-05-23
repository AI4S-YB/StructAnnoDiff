import tempfile
import unittest
from pathlib import Path

from summarize import (
    build_comparison_by_feature_path_table,
    build_comparison_table,
    configure_paths,
    parse_agat_report_tables,
    parse_agat_statistics_sections,
    select_coding_stats_section,
)


class SummarizeParsingTests(unittest.TestCase):
    def test_agat_stats_sections_do_not_overwrite_metric_names(self):
        content = "\n".join([
            "------------------------------------- mrna -------------------------------------",
            "Number of gene                                              54347",
            "mean exons per mrna                                         4.2",
            "------------------------------------- trna -------------------------------------",
            "Number of gene                                              1092",
            "",
        ])
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "stats.txt"
            path.write_text(content, encoding="utf-8")
            sections = parse_agat_statistics_sections(path)

        self.assertEqual(sections["mrna"]["Number of gene"], 54347.0)
        self.assertEqual(sections["trna"]["Number of gene"], 1092.0)
        section, metrics = select_coding_stats_section(sections)
        self.assertEqual(section, "mrna")
        self.assertEqual(metrics["Number of gene"], 54347.0)

    def test_agat_compare_tables_are_kept_separate(self):
        report = "\n".join([
            "|                                        gene@mrna@cds                                       |",
            "|              0               |              1               |            42726             |",
            "|                                        gene@mrna@exon                                      |",
            "|              0               |              1               |              24              |",
            "|                                    gene@transcript@cds                                     |",
            "|              1               |              0               |            64295             |",
            "",
        ])
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "report.txt"
            path.write_text(report, encoding="utf-8")
            tables = parse_agat_report_tables(path)

        self.assertEqual(tables["gene@mrna@cds"][(0, 1)], 42726)
        self.assertEqual(tables["gene@mrna@exon"][(0, 1)], 24)
        self.assertEqual(tables["gene@transcript@cds"][(1, 0)], 64295)

    def test_comparison_matrix_matches_mrna_transcript_gene_scope(self):
        report = "\n".join([
            "|                                        gene@mrna@cds                                       |",
            "|              0               |              1               |              10              |",
            "|                                        gene@mrna@exon                                      |",
            "|              0               |              1               |               5              |",
            "|                                    gene@transcript@cds                                     |",
            "|              1               |              0               |               7              |",
            "|                                   gene@transcript@exon                                    |",
            "|              1               |              0               |               3              |",
            "|                                      gene@ncrna@exon                                      |",
            "|              0               |              1               |              99              |",
            "",
        ])
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            compare_dir = root / "compare" / "Cucumber"
            compare_dir.mkdir(parents=True)
            (compare_dir / "report.txt").write_text(report, encoding="utf-8")
            configure_paths(root)

            coding = build_comparison_table()
            by_path = build_comparison_by_feature_path_table()

        row = coding[coding["species"] == "Cucumber"].iloc[0]
        self.assertEqual(row["added"], 15)
        self.assertEqual(row["removed"], 10)
        self.assertEqual(row["total_after"], 15)
        self.assertIn("gene@mrna@exon", set(by_path["feature_path"]))
        self.assertIn("gene@ncrna@exon", set(by_path["feature_path"]))


if __name__ == "__main__":
    unittest.main()
