package locuscompare;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public final class CoreTests {
    private CoreTests() {
    }

    public static void main(String[] args) throws Exception {
        testNormalizeAttributeIds();
        testDefaultScopeExcludesGenesWithoutMrna();
        testMultiParentChildren();
        testReciprocalOverlapAcceptsContainment();
        testOverlapScoreUsesBestTranscriptPairExonLength();
        testHybridGraphPrunesWeakBridgeBetweenStrongPairs();
        testSameIdsAreDistinctGraphNodes();
        testPrimaryMrnaTieUsesFileOrder();
        testUtrDerivationWithoutExons();
        testPhaseOnlyCdsDifferenceIsIgnored();
        testUtrExonGainIsNotCodingGainWhenCdsBoundaryChanges();
        testNoOverlapLociIgnoreStrand();
        testRepresentativeTranscriptChangeSummary();
        testCurationPanelMetricsTable();
        testCurationSummaryTable();
        testNoCommonSeqidsReturnsNull();
        testWeakContainmentOverlapIsNotCountedAsNovelDeleted();
        System.out.println("CoreTests: OK");
    }

    private static void testNormalizeAttributeIds() {
        assertEquals("Gene1", Core.normalizeAttributeId(List.of("gene:Gene1"), "gene:"), "gene prefix");
        assertEquals("Tx1", Core.normalizeAttributeId(List.of("transcript:Tx1"), "transcript:"), "transcript prefix");
        assertEquals(Arrays.asList("Tx1", "Tx2"),
                Core.normalizeAttributeIds(List.of("transcript:Tx1,transcript:Tx2"), "transcript:"),
                "multi parent");
    }

    private static void testDefaultScopeExcludesGenesWithoutMrna() throws Exception {
        Path dir = Files.createTempDirectory("jlc-test-");
        Path gff = dir.resolve("test.gff3");
        Files.writeString(gff, String.join("\n",
                "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=Gene1",
                "chr1\tsrc\tgene\t200\t300\t.\t+\t.\tID=Gene2",
                "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=Tx1;Parent=Gene1",
                "chr1\tsrc\texon\t1\t100\t.\t+\t.\tParent=Tx1",
                ""));
        Map<String, Core.GeneModel> mrnaScope = Core.parseGff3ToModels(gff, "mrna");
        Map<String, Core.GeneModel> allScope = Core.parseGff3ToModels(gff, "all");
        assertTrue(mrnaScope.containsKey("Gene1") && !mrnaScope.containsKey("Gene2"), "mrna scope");
        assertTrue(allScope.containsKey("Gene1") && !allScope.containsKey("Gene2"), "all scope");
    }

    private static void testMultiParentChildren() throws Exception {
        Path dir = Files.createTempDirectory("jlc-test-");
        Path gff = dir.resolve("test.gff3");
        Files.writeString(gff, String.join("\n",
                "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=Gene1",
                "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=Tx1;Parent=Gene1",
                "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=Tx2;Parent=Gene1",
                "chr1\tsrc\texon\t1\t100\t.\t+\t.\tParent=Tx1,Tx2",
                "chr1\tsrc\tCDS\t10\t90\t.\t+\t0\tParent=Tx1,Tx2",
                ""));
        Core.GeneModel gene = Core.parseGff3ToModels(gff, "mrna").get("Gene1");
        assertEquals(2, gene.mrnaCount(), "mrna count");
        assertEquals(2, gene.totalExonCount(), "shared exons");
        assertEquals(162, gene.totalCdsLength(), "shared cds");
    }

    private static void testReciprocalOverlapAcceptsContainment() {
        Core.GeneModel small = gene("small", 101, 200);
        Core.GeneModel large = gene("large", 1, 1000);
        assertEquals(1.0, Core.reciprocalOverlap(small, large), "reciprocal overlap");
        assertEquals(1.0, Core.containmentOverlap(small, large), "containment overlap");
        assertEquals(1, Core.findOverlappingPairs(List.of(small), List.of(large), 0.5,
                "reciprocal", new LinkedHashMap<>()).size(), "strict pairs");
        assertEquals(1, Core.findOverlappingPairs(List.of(small), List.of(large), 0.5,
                "containment", new LinkedHashMap<>()).size(), "containment pairs");
    }

    private static void testOverlapScoreUsesBestTranscriptPairExonLength() {
        Core.GeneModel before = geneWithMrnas("before", List.of(
                mrna("b.match", new int[][]{{100, 199}}, new Object[][]{{100, 199, "0"}}, null),
                mrna("b.alt1", new int[][]{{1000, 1099}}, new Object[][]{{1000, 1099, "0"}}, null),
                mrna("b.alt2", new int[][]{{2000, 2099}}, new Object[][]{{2000, 2099, "0"}}, null)
        ));
        Core.GeneModel after = geneWithMrnas("after", List.of(
                mrna("a.match", new int[][]{{100, 199}}, new Object[][]{{100, 199, "0"}}, null),
                mrna("a.alt1", new int[][]{{3000, 3099}}, new Object[][]{{3000, 3099, "0"}}, null),
                mrna("a.alt2", new int[][]{{4000, 4099}}, new Object[][]{{4000, 4099, "0"}}, null)
        ));

        Core.TranscriptOverlapMetrics best = Core.bestTranscriptOverlapMetrics(before, after);

        assertEquals(100, Core.transcriptFeatureLength(before.mrnas.get(0)), "transcript length");
        assertEquals("b.match", best.leftMrna.mrnaId, "best before transcript");
        assertEquals("a.match", best.rightMrna.mrnaId, "best after transcript");
        assertEquals(100, Core.featureOverlapLen(before, after), "best overlap");
        assertEquals(1.0, Core.reciprocalOverlap(before, after), "best reciprocal");
        assertEquals(1.0, best.jaccard, "best jaccard");
    }

    private static void testHybridGraphPrunesWeakBridgeBetweenStrongPairs() throws Exception {
        Path dir = Files.createTempDirectory("jlc-test-");
        Path before = dir.resolve("before.gff3");
        Path after = dir.resolve("after.gff3");
        Files.writeString(before, String.join("\n",
                "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=g1",
                "chr1\ttest\tmRNA\t100\t200\t.\t+\t.\tID=t1;Parent=g1",
                "chr1\ttest\texon\t100\t200\t.\t+\t.\tParent=t1",
                "chr1\ttest\tCDS\t100\t200\t.\t+\t0\tParent=t1",
                "chr1\ttest\tgene\t220\t240\t.\t+\t.\tID=g2",
                "chr1\ttest\tmRNA\t220\t240\t.\t+\t.\tID=t2;Parent=g2",
                "chr1\ttest\texon\t220\t240\t.\t+\t.\tParent=t2",
                "chr1\ttest\tCDS\t220\t240\t.\t+\t0\tParent=t2",
                ""));
        Files.writeString(after, String.join("\n",
                "chr1\ttest\tgene\t100\t230\t.\t+\t.\tID=a1",
                "chr1\ttest\tmRNA\t100\t230\t.\t+\t.\tID=ta1;Parent=a1",
                "chr1\ttest\texon\t100\t230\t.\t+\t.\tParent=ta1",
                "chr1\ttest\tCDS\t100\t200\t.\t+\t0\tParent=ta1",
                "chr1\ttest\tgene\t220\t240\t.\t+\t.\tID=a2",
                "chr1\ttest\tmRNA\t220\t240\t.\t+\t.\tID=ta2;Parent=a2",
                "chr1\ttest\texon\t220\t240\t.\t+\t.\tParent=ta2",
                "chr1\ttest\tCDS\t220\t240\t.\t+\t0\tParent=ta2",
                ""));

        Core.ComparisonResult result = Core.compareAnnotations(before, after, 0.5, 10, 0.1, 0.1, "mrna", "hybrid");

        assertEquals(2, Core.intObject(result.summary.get("syntenic_total")), "bridge syntenic");
        assertEquals(0, Core.intObject(result.summary.get("complex_events")), "bridge complex");
        assertEquals(1, Core.intObject(result.summary.get("containment_bridge_edges_pruned")), "bridge pruned");
    }

    private static void testSameIdsAreDistinctGraphNodes() {
        LinkedHashMap<String, Core.GeneModel> before = new LinkedHashMap<>();
        before.put("A", gene("A", 1, 100));
        before.put("B", gene("B", 201, 300));
        LinkedHashMap<String, Core.GeneModel> after = new LinkedHashMap<>();
        after.put("A", gene("A", 201, 300));
        after.put("B", gene("B", 1, 100));
        List<Core.Pair> pairs = List.of(
                new Core.Pair(before.get("A"), after.get("B"), 1.0, 1.0),
                new Core.Pair(before.get("B"), after.get("A"), 1.0, 1.0)
        );
        Core.ResolvedMatches result = Core.resolveMatches(pairs, before, after);
        assertEquals(2, result.syntenic.size(), "same IDs syntenic");
        assertEquals(0, result.complex.size(), "same IDs complex");
    }

    private static void testPrimaryMrnaTieUsesFileOrder() {
        Core.MRNAModel first = mrna("tx_a", new int[][]{{1, 100}}, new Object[][]{{1, 90, "0"}}, null);
        Core.MRNAModel second = mrna("tx_z", new int[][]{{200, 299}}, new Object[][]{{200, 289, "0"}}, null);
        Core.GeneModel gene = gene("gene", 1, 300);
        gene.mrnas.add(first);
        gene.mrnas.add(second);
        assertSame(first, gene.primaryMrna(), "primary tie");
    }

    private static void testUtrDerivationWithoutExons() throws Exception {
        Path dir = Files.createTempDirectory("jlc-test-");
        Path gff = dir.resolve("input.gff3");
        Files.writeString(gff, String.join("\n",
                "chr1\ttest\tgene\t1\t300\t.\t+\t.\tID=g1",
                "chr1\ttest\tmRNA\t1\t300\t.\t+\t.\tID=t1;Parent=g1",
                "chr1\ttest\tCDS\t101\t200\t.\t+\t0\tID=cds1;Parent=t1",
                ""));
        Core.MRNAModel mrna = Core.parseGff3ToModels(gff, "mrna").get("g1").mrnas.get(0);
        assertEquals(0, mrna.utrs.size(), "no raw-span derived utr");
        assertEquals("", Core.utrSignature(mrna), "no raw-span utr signature");
        assertEquals("101:200", Core.exonSignature(mrna), "cds-derived exon");
    }

    private static void testPhaseOnlyCdsDifferenceIsIgnored() {
        Core.GeneModel before = geneWithMrnas("before", List.of(
                mrna("b1", new int[][]{{1, 100}}, new Object[][]{{10, 90, "0"}}, null)
        ));
        Core.GeneModel after = geneWithMrnas("after", List.of(
                mrna("a1", new int[][]{{1, 100}}, new Object[][]{{10, 90, "2"}}, null)
        ));

        String subtype = Core.classifySyntenicChange(before, after, 10, 0.1, 0.1);
        Map<String, Boolean> attrs = Core.computeSyntenicAttributes(before, after, 10, 0.1, 0.1);
        LinkedHashMap<String, Integer> summary = Core.summarizeRepresentativeTranscriptChanges(
                List.of(new Core.GenePair(before, after))
        );

        assertEquals("exact", subtype, "phase-only subtype");
        assertTrue(attrs.get("exact"), "phase-only exact");
        assertTrue(!attrs.get("cds_change"), "phase-only not cds change");
        assertEquals(0, summary.get("rep_structural_changed").intValue(), "phase-only not structural");
        assertEquals(0, summary.get("rep_cds_boundary_changed_same_count").intValue(),
                "phase-only not cds boundary");
    }

    private static void testUtrExonGainIsNotCodingGainWhenCdsBoundaryChanges() {
        Core.GeneModel before = geneWithMrnas("before", List.of(
                mrna("b1", new int[][]{{20, 80}}, new Object[][]{{20, 80, "0"}}, null)
        ));
        Core.GeneModel after = geneWithMrnas("after", List.of(
                mrna("a1", new int[][]{{1, 10}, {20, 81}}, new Object[][]{{20, 81, "0"}},
                        new Object[][]{{1, 10, "five_prime_UTR"}})
        ));
        Map<String, Boolean> attrs = Core.computeSyntenicAttributes(before, after, 10, 0.1, 0.1);
        String subtype = Core.classifySyntenicChange(before, after, 10, 0.1, 0.1);
        assertTrue(attrs.get("utr_exon_added"), "utr exon added");
        assertTrue(!attrs.get("coding_exon_gain"), "not coding gain");
        assertTrue(subtype.contains("utr_exon_added"), "subtype has utr exon");
        assertTrue(subtype.contains("cds_extended"), "subtype has cds extended");
        assertTrue(!subtype.contains("exon_gain"), "subtype lacks exon gain");
    }

    private static void testNoOverlapLociIgnoreStrand() {
        Core.GeneModel beforeOverlap = gene("before_overlap", 100, 200);
        Core.GeneModel beforeDeleted = gene("before_deleted", 300, 400);
        Core.GeneModel afterOverlap = new Core.GeneModel(
                "after_overlap", "chr1", 150, 250, "-", "test", 150, 250);
        Core.GeneModel afterNew = gene("after_new", 500, 600);
        LinkedHashMap<String, Core.GeneModel> before = new LinkedHashMap<>();
        before.put(beforeOverlap.geneId, beforeOverlap);
        before.put(beforeDeleted.geneId, beforeDeleted);
        LinkedHashMap<String, Core.GeneModel> after = new LinkedHashMap<>();
        after.put(afterOverlap.geneId, afterOverlap);
        after.put(afterNew.geneId, afterNew);

        LinkedHashMap<String, Integer> counts = Core.countNoOverlapLoci(before, after);

        assertEquals(1, counts.get("no_overlap_after_loci").intValue(), "no-overlap after");
        assertEquals(1, counts.get("no_overlap_before_loci").intValue(), "no-overlap before");
    }

    private static void testRepresentativeTranscriptChangeSummary() {
        Core.GeneModel pairExonCountBefore = geneWithMrnas("b_exon_count", List.of(
                mrna("bt1", new int[][]{{1, 50}, {101, 150}},
                        new Object[][]{{10, 40, "0"}, {110, 140, "0"}}, null)
        ));
        Core.GeneModel pairExonCountAfter = geneWithMrnas("a_exon_count", List.of(
                mrna("at1", new int[][]{{1, 50}, {101, 150}, {201, 250}},
                        new Object[][]{{10, 40, "0"}, {110, 140, "0"}}, null)
        ));
        Core.GeneModel pairBoundaryBefore = geneWithMrnas("b_boundary", List.of(
                mrna("bt2", new int[][]{{1, 50}, {101, 150}},
                        new Object[][]{{10, 40, "0"}, {110, 140, "0"}}, null)
        ));
        Core.GeneModel pairBoundaryAfter = geneWithMrnas("a_boundary", List.of(
                mrna("at2", new int[][]{{1, 60}, {101, 150}},
                        new Object[][]{{10, 40, "0"}, {110, 145, "0"}}, null)
        ));
        Core.GeneModel pairCdsCountBefore = geneWithMrnas("b_cds_count", List.of(
                mrna("bt3", new int[][]{{1, 200}}, new Object[][]{{10, 40, "0"}}, null)
        ));
        Core.GeneModel pairCdsCountAfter = geneWithMrnas("a_cds_count", List.of(
                mrna("at3", new int[][]{{1, 200}},
                        new Object[][]{{10, 40, "0"}, {110, 140, "0"}}, null)
        ));

        LinkedHashMap<String, Integer> summary = Core.summarizeRepresentativeTranscriptChanges(List.of(
                new Core.GenePair(pairExonCountBefore, pairExonCountAfter),
                new Core.GenePair(pairBoundaryBefore, pairBoundaryAfter),
                new Core.GenePair(pairCdsCountBefore, pairCdsCountAfter)
        ));

        assertEquals(3, summary.get("rep_transcript_pairs").intValue(), "rep pairs");
        assertEquals(3, summary.get("rep_structural_changed").intValue(), "rep structural changed");
        assertEquals(1, summary.get("rep_exon_count_changed").intValue(), "rep exon count");
        assertEquals(1, summary.get("rep_exon_boundary_changed_same_count").intValue(), "rep exon boundary");
        assertEquals(1, summary.get("rep_cds_count_changed").intValue(), "rep cds count");
        assertEquals(1, summary.get("rep_cds_boundary_changed_same_count").intValue(), "rep cds boundary");
    }

    private static void testCurationPanelMetricsTable() {
        LinkedHashMap<String, Object> core = new LinkedHashMap<>();
        core.put("species_id", "Demo_species");
        core.put("Species", "Demo");
        core.put("total_before_genes", 100);
        core.put("total_after_genes", 80);
        core.put("deleted_loci_no_overlap", 5);
        core.put("new_loci_no_overlap", 7);
        core.put("split_events", 2);
        core.put("merge_events", 3);
        core.put("rep_exon_changed", 25);

        List<LinkedHashMap<String, Object>> rows = GenerateTablesMain.buildCurationPanelMetricsTable(List.of(core));

        assertEquals(6, rows.size(), "panel row count");
        assertEquals("A", rows.get(0).get("panel"), "panel A deleted");
        assertEquals("deleted_loci_no_overlap", rows.get(0).get("metric"), "deleted metric");
        assertEquals(5, rows.get(0).get("count"), "deleted count");
        assertEquals(100, rows.get(0).get("denominator"), "deleted denominator");
        assertEquals(5.0, ((Number) rows.get(0).get("percent")).doubleValue(), "deleted percent");
        assertEquals("", rows.get(2).get("percent"), "event percent empty");
        assertEquals(25.0, ((Number) rows.get(4).get("axis_value")).doubleValue(), "before axis pct");
        assertEquals(31.25, ((Number) rows.get(5).get("axis_value")).doubleValue(), "after axis pct");
    }

    private static void testCurationSummaryTable() throws Exception {
        Path dir = Files.createTempDirectory("jlc-test-");
        Path locus = dir.resolve("results").resolve("locus");
        Files.createDirectories(locus);
        Files.writeString(dir.resolve("Demo_species.before.gff3"), "");
        Files.writeString(dir.resolve("Demo_species.after.gff3"), "");

        LinkedHashMap<String, Object> summary = new LinkedHashMap<>();
        summary.put("total_before_genes", 100);
        summary.put("total_after_genes", 80);
        summary.put("no_overlap_before_loci", 5);
        summary.put("no_overlap_after_loci", 7);
        summary.put("split_events", 2);
        summary.put("merge_events", 3);
        summary.put("rep_exon_count_changed", 10);
        summary.put("rep_exon_boundary_changed_same_count", 15);
        Core.writeCsv(locus.resolve("Demo_species_change_summary.csv"),
                List.copyOf(summary.keySet()), List.of(summary));

        List<LinkedHashMap<String, Object>> rows = GenerateTablesMain.buildCurationSummaryTable(
                List.of(new Core.Species("Demo_species", "Demo", "Demo")), dir, locus);
        LinkedHashMap<String, Object> row = rows.get(0);

        assertEquals(1, rows.size(), "summary table row count");
        assertEquals(List.of(
                "input_file_1", "input_file_2", "deleted_loci", "new_loci", "split_events", "merge_events",
                "exon_changes_in_file_1", "exon_changes_in_file_2"
        ), List.copyOf(row.keySet()), "summary table columns");
        assertEquals("Demo_species.before.gff3", row.get("input_file_1"), "input file 1");
        assertEquals("Demo_species.after.gff3", row.get("input_file_2"), "input file 2");
        assertEquals(5, row.get("deleted_loci"), "deleted loci");
        assertEquals(7, row.get("new_loci"), "new loci");
        assertEquals(2, row.get("split_events"), "split events");
        assertEquals(3, row.get("merge_events"), "merge events");
        assertEquals("25.00% (25/100)", row.get("exon_changes_in_file_1"), "before ratio");
        assertEquals("31.25% (25/80)", row.get("exon_changes_in_file_2"), "after ratio");

        List<LinkedHashMap<String, Object>> infoRows = GenerateTablesMain.buildCurationSummaryColumnInfoTable();
        assertEquals(8, infoRows.size(), "summary column info row count");
        assertEquals("input_file_1", infoRows.get(0).get("column"), "summary column info first column");
        assertTrue(String.valueOf(infoRows.get(0).get("description")).contains("original annotation"),
                "summary column info description");
    }

    private static void testNoCommonSeqidsReturnsNull() throws Exception {
        Path dir = Files.createTempDirectory("jlc-test-");
        Path before = dir.resolve("before.gff3");
        Path after = dir.resolve("after.gff3");
        Files.writeString(before, String.join("\n",
                "chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=g1",
                "chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1",
                "chr1\ttest\texon\t1\t100\t.\t+\t.\tParent=t1",
                "chr1\ttest\tCDS\t10\t90\t.\t+\t0\tParent=t1",
                ""));
        Files.writeString(after, String.join("\n",
                "chr2\ttest\tgene\t1\t100\t.\t+\t.\tID=g2",
                "chr2\ttest\tmRNA\t1\t100\t.\t+\t.\tID=t2;Parent=g2",
                "chr2\ttest\texon\t1\t100\t.\t+\t.\tParent=t2",
                "chr2\ttest\tCDS\t10\t90\t.\t+\t0\tParent=t2",
                ""));
        Core.ComparisonResult result = Core.compareAnnotations(before, after, 0.5, 10, 0.1, 0.1, "mrna", "reciprocal");
        assertTrue(result == null, "no common seqids");
    }

    private static void testWeakContainmentOverlapIsNotCountedAsNovelDeleted() throws Exception {
        Path dir = Files.createTempDirectory("jlc-test-");
        Path before = dir.resolve("before.gff3");
        Path after = dir.resolve("after.gff3");
        Files.writeString(before, String.join("\n",
                "chr1\ttest\tgene\t101\t200\t.\t+\t.\tID=g1",
                "chr1\ttest\tmRNA\t101\t200\t.\t+\t.\tID=t1;Parent=g1",
                "chr1\ttest\texon\t101\t200\t.\t+\t.\tParent=t1",
                "chr1\ttest\tCDS\t101\t200\t.\t+\t0\tParent=t1",
                ""));
        Files.writeString(after, String.join("\n",
                "chr1\ttest\tgene\t1\t1000\t.\t+\t.\tID=g2",
                "chr1\ttest\tmRNA\t1\t1000\t.\t+\t.\tID=t2;Parent=g2",
                "chr1\ttest\texon\t1\t1000\t.\t+\t.\tParent=t2",
                "chr1\ttest\tCDS\t1\t1000\t.\t+\t0\tParent=t2",
                ""));
        Core.ComparisonResult result = Core.compareAnnotations(before, after, 0.5, 10, 0.1, 0.1, "mrna", "reciprocal");
        assertEquals(1, Core.intObject(result.summary.get("syntenic_total")), "containment syntenic");
        assertEquals(0, Core.intObject(result.summary.get("novel_genes")), "strict novel");
        assertEquals(0, Core.intObject(result.summary.get("deleted_genes")), "strict deleted");
        assertEquals(0, Core.intObject(result.summary.get("unresolved_overlap_after_genes")), "unresolved after");
        assertEquals(0, Core.intObject(result.summary.get("unresolved_overlap_before_genes")), "unresolved before");
    }

    private static Core.GeneModel gene(String id, int start, int end) {
        return new Core.GeneModel(id, "chr1", start, end, "+", "test", start, end);
    }

    private static Core.GeneModel geneWithMrnas(String id, List<Core.MRNAModel> mrnas) {
        Core.GeneModel gene = gene(id, 1, 300);
        gene.mrnas.addAll(mrnas);
        return gene;
    }

    private static Core.MRNAModel mrna(String id, int[][] exons, Object[][] cds, Object[][] utrs) {
        Core.MRNAModel mrna = new Core.MRNAModel(id, null, null, 0);
        for (int[] exon : exons) {
            mrna.exons.add(new Core.ExonFeature(exon[0], exon[1]));
        }
        for (Object[] c : cds) {
            mrna.cds.add(new Core.CDSFeature((Integer) c[0], (Integer) c[1], (String) c[2]));
        }
        if (utrs != null) {
            for (Object[] u : utrs) {
                mrna.utrs.add(new Core.UTRFeature((Integer) u[0], (Integer) u[1], (String) u[2]));
            }
        }
        return mrna;
    }

    private static void assertEquals(Object expected, Object actual, String message) {
        if (!expected.equals(actual)) {
            throw new AssertionError(message + ": expected " + expected + ", got " + actual);
        }
    }

    private static void assertEquals(double expected, double actual, String message) {
        if (Math.abs(expected - actual) > 1e-9) {
            throw new AssertionError(message + ": expected " + expected + ", got " + actual);
        }
    }

    private static void assertTrue(boolean condition, String message) {
        if (!condition) {
            throw new AssertionError(message);
        }
    }

    private static void assertSame(Object expected, Object actual, String message) {
        if (expected != actual) {
            throw new AssertionError(message + ": objects differ");
        }
    }
}
