package locuscompare;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

public final class CoreTests {
    private CoreTests() {
    }

    public static void main(String[] args) throws Exception {
        testNormalizeAttributeIds();
        testDefaultScopeExcludesGenesWithoutMrna();
        testMultiParentChildren();
        testReciprocalOverlapRejectsContainment();
        testWeakRejectedDoesNotDependOnDiagnostics();
        testSameIdsAreDistinctGraphNodes();
        testPrimaryMrnaTieUsesFileOrder();
        testUtrDerivationWithoutExons();
        testUtrExonGainIsNotCodingGainWhenCdsBoundaryChanges();
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
        assertTrue(allScope.containsKey("Gene1") && allScope.containsKey("Gene2"), "all scope");
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

    private static void testReciprocalOverlapRejectsContainment() {
        Core.GeneModel small = gene("small", 101, 200);
        Core.GeneModel large = gene("large", 1, 1000);
        assertEquals(0.1, Core.reciprocalOverlap(small, large), "reciprocal overlap");
        assertEquals(1.0, Core.containmentOverlap(small, large), "containment overlap");
        assertEquals(0, Core.findOverlappingPairs(List.of(small), List.of(large), 0.5,
                "reciprocal", new LinkedHashMap<>()).size(), "strict pairs");
        assertEquals(1, Core.findOverlappingPairs(List.of(small), List.of(large), 0.5,
                "containment", new LinkedHashMap<>()).size(), "containment pairs");
    }

    private static void testWeakRejectedDoesNotDependOnDiagnostics() {
        Core.GeneModel small = gene("small", 101, 200);
        Core.GeneModel large = gene("large", 1, 1000);
        Map<String, Set<String>> weakRejected = new HashMap<>();
        weakRejected.put("before", new HashSet<>());
        weakRejected.put("after", new HashSet<>());

        List<Core.Pair> pairs = Core.findOverlappingPairs(
                List.of(small), List.of(large), 0.5, "reciprocal", null, weakRejected);

        assertEquals(0, pairs.size(), "strict weak pairs");
        assertEquals(Set.of("small"), weakRejected.get("before"), "weak before");
        assertEquals(Set.of("large"), weakRejected.get("after"), "weak after");
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
        assertEquals(2, mrna.utrs.size(), "derived utr count");
        assertEquals("1:100:five_prime_UTR|201:300:three_prime_UTR", Core.utrSignature(mrna), "derived utr signature");
        assertEquals("1:300", Core.exonSignature(mrna), "derived exon");
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
        assertEquals(0, Core.intObject(result.summary.get("novel_genes")), "strict novel");
        assertEquals(0, Core.intObject(result.summary.get("deleted_genes")), "strict deleted");
        assertEquals(1, Core.intObject(result.summary.get("unresolved_overlap_after_genes")), "unresolved after");
        assertEquals(1, Core.intObject(result.summary.get("unresolved_overlap_before_genes")), "unresolved before");
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
