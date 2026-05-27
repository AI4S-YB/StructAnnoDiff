package locuscompare;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Predicate;
import java.util.zip.GZIPInputStream;

public final class Core {
    private Core() {
    }

    public static final List<String> SYNTENIC_ATTRIBUTE_KEYS = Collections.unmodifiableList(Arrays.asList(
            "exact",
            "gene_boundary_changed",
            "utr_added",
            "utr_lost",
            "utr_exon_added",
            "utr_exon_removed",
            "utr_refined",
            "coding_exon_gain",
            "coding_exon_loss",
            "exon_boundary_refined",
            "cds_change",
            "cds_boundary_refined",
            "isoform_change"
    ));

    public static final List<String> CHANGE_LOG_COLUMNS = Collections.unmodifiableList(Arrays.asList(
            "before_gene",
            "after_gene",
            "seqid",
            "before_start",
            "before_end",
            "after_start",
            "after_end",
            "before_gene_start",
            "before_gene_end",
            "after_gene_start",
            "after_gene_end",
            "strand",
            "match_type",
            "change_subtype",
            "before_length",
            "after_length",
            "before_gene_length",
            "after_gene_length",
            "before_exons",
            "after_exons",
            "before_cds",
            "after_cds",
            "before_mrnas",
            "after_mrnas"
    ));

    public static final class ExonFeature {
        public final int start;
        public final int end;

        public ExonFeature(int start, int end) {
            this.start = start;
            this.end = end;
        }
    }

    public static final class CDSFeature {
        public final int start;
        public final int end;
        public final String phase;

        public CDSFeature(int start, int end, String phase) {
            this.start = start;
            this.end = end;
            this.phase = phase;
        }
    }

    public static final class UTRFeature {
        public final int start;
        public final int end;
        public final String utrType;

        public UTRFeature(int start, int end, String utrType) {
            this.start = start;
            this.end = end;
            this.utrType = utrType;
        }
    }

    public static final class MRNAModel {
        public final String mrnaId;
        public final List<ExonFeature> exons = new ArrayList<>();
        public final List<CDSFeature> cds = new ArrayList<>();
        public final List<UTRFeature> utrs = new ArrayList<>();
        public Integer rawStart;
        public Integer rawEnd;
        public boolean hasExplicitExon = false;
        public boolean exonsAreInferred = false;
        private final int order;

        public MRNAModel(String mrnaId, Integer rawStart, Integer rawEnd, int order) {
            this.mrnaId = mrnaId;
            this.rawStart = rawStart;
            this.rawEnd = rawEnd;
            this.order = order;
        }

        public int exonCount() {
            return exons.size();
        }

        public int cdsLength() {
            int total = 0;
            for (CDSFeature c : cds) {
                total += c.end - c.start + 1;
            }
            return total;
        }

        public boolean hasFivePrimeUtr() {
            return utrs.stream().anyMatch(u -> "five_prime_UTR".equals(u.utrType));
        }

        public boolean hasThreePrimeUtr() {
            return utrs.stream().anyMatch(u -> "three_prime_UTR".equals(u.utrType));
        }

        public int fivePrimeUtrLength() {
            int total = 0;
            for (UTRFeature u : utrs) {
                if ("five_prime_UTR".equals(u.utrType)) {
                    total += u.end - u.start + 1;
                }
            }
            return total;
        }

        public int threePrimeUtrLength() {
            int total = 0;
            for (UTRFeature u : utrs) {
                if ("three_prime_UTR".equals(u.utrType)) {
                    total += u.end - u.start + 1;
                }
            }
            return total;
        }
    }

    public static final class GeneModel {
        public final String geneId;
        public final String seqid;
        public int start;
        public int end;
        public final String strand;
        public final String source;
        public final List<MRNAModel> mrnas = new ArrayList<>();
        public Integer rawStart;
        public Integer rawEnd;

        public GeneModel(String geneId, String seqid, int start, int end, String strand,
                         String source, Integer rawStart, Integer rawEnd) {
            this.geneId = geneId;
            this.seqid = seqid;
            this.start = start;
            this.end = end;
            this.strand = strand;
            this.source = source;
            this.rawStart = rawStart;
            this.rawEnd = rawEnd;
        }

        public int length() {
            return end - start + 1;
        }

        public int mrnaCount() {
            return mrnas.size();
        }

        public int totalExonCount() {
            int total = 0;
            for (MRNAModel m : mrnas) {
                total += m.exonCount();
            }
            return total;
        }

        public int totalCdsLength() {
            int total = 0;
            for (MRNAModel m : mrnas) {
                total += m.cdsLength();
            }
            return total;
        }

        public MRNAModel primaryMrna() {
            MRNAModel best = null;
            for (MRNAModel m : mrnas) {
                if (best == null
                        || m.cdsLength() > best.cdsLength()
                        || (m.cdsLength() == best.cdsLength() && m.exonCount() > best.exonCount())) {
                    best = m;
                }
            }
            return best;
        }
    }

    public static final class GffFeature {
        public final String seqid;
        public final String source;
        public final String type;
        public final int start;
        public final int end;
        public final String score;
        public final String strand;
        public final String phase;
        public final Map<String, List<String>> attributes;

        GffFeature(String seqid, String source, String type, int start, int end,
                   String score, String strand, String phase, Map<String, List<String>> attributes) {
            this.seqid = seqid;
            this.source = source;
            this.type = type;
            this.start = start;
            this.end = end;
            this.score = score;
            this.strand = strand;
            this.phase = phase;
            this.attributes = attributes;
        }
    }

    public static final class Pair {
        public final GeneModel before;
        public final GeneModel after;
        public final double score;
        public final double jaccard;

        Pair(GeneModel before, GeneModel after, double score, double jaccard) {
            this.before = before;
            this.after = after;
            this.score = score;
            this.jaccard = jaccard;
        }
    }

    public static final class TranscriptOverlapMetrics {
        public final double score;
        public final double jaccard;
        public final int overlap;
        public final int leftLength;
        public final int rightLength;
        public final MRNAModel leftMrna;
        public final MRNAModel rightMrna;

        TranscriptOverlapMetrics(double score, double jaccard, int overlap,
                                 int leftLength, int rightLength,
                                 MRNAModel leftMrna, MRNAModel rightMrna) {
            this.score = score;
            this.jaccard = jaccard;
            this.overlap = overlap;
            this.leftLength = leftLength;
            this.rightLength = rightLength;
            this.leftMrna = leftMrna;
            this.rightMrna = rightMrna;
        }
    }

    private static final class TranscriptIntervals {
        final MRNAModel mrna;
        final List<int[]> intervals;

        TranscriptIntervals(MRNAModel mrna, List<int[]> intervals) {
            this.mrna = mrna;
            this.intervals = intervals;
        }
    }

    private static final class PruneResult {
        final List<Pair> pairs;
        final int prunedCount;

        PruneResult(List<Pair> pairs, int prunedCount) {
            this.pairs = pairs;
            this.prunedCount = prunedCount;
        }
    }

    public static final class ResolvedMatches {
        public final List<GenePair> syntenic = new ArrayList<>();
        public final List<SplitEvent> split = new ArrayList<>();
        public final List<MergeEvent> merge = new ArrayList<>();
        public final List<ComplexEvent> complex = new ArrayList<>();
        public final List<GeneModel> novel = new ArrayList<>();
        public final List<GeneModel> deleted = new ArrayList<>();
    }

    public static final class GenePair {
        public final GeneModel before;
        public final GeneModel after;

        GenePair(GeneModel before, GeneModel after) {
            this.before = before;
            this.after = after;
        }
    }

    public static final class SplitEvent {
        public final GeneModel before;
        public final List<GeneModel> afters;

        SplitEvent(GeneModel before, List<GeneModel> afters) {
            this.before = before;
            this.afters = afters;
        }
    }

    public static final class MergeEvent {
        public final List<GeneModel> befores;
        public final GeneModel after;

        MergeEvent(List<GeneModel> befores, GeneModel after) {
            this.befores = befores;
            this.after = after;
        }
    }

    public static final class ComplexEvent {
        public final List<GeneModel> befores;
        public final List<GeneModel> afters;
        public final List<GenePair> edgePairs;

        ComplexEvent(List<GeneModel> befores, List<GeneModel> afters, List<GenePair> edgePairs) {
            this.befores = befores;
            this.afters = afters;
            this.edgePairs = edgePairs;
        }
    }

    public static final class ComparisonResult {
        public final LinkedHashMap<String, Object> summary;
        public final List<LinkedHashMap<String, Object>> changeLog;

        ComparisonResult(LinkedHashMap<String, Object> summary, List<LinkedHashMap<String, Object>> changeLog) {
            this.summary = summary;
            this.changeLog = changeLog;
        }
    }

    private static final class Node {
        final String side;
        final String id;

        Node(String side, String id) {
            this.side = side;
            this.id = id;
        }

        @Override
        public boolean equals(Object o) {
            if (!(o instanceof Node)) {
                return false;
            }
            Node other = (Node) o;
            return side.equals(other.side) && id.equals(other.id);
        }

        @Override
        public int hashCode() {
            return Objects.hash(side, id);
        }
    }

    public static List<GffFeature> parseGff3(Path path) throws IOException {
        List<GffFeature> features = new ArrayList<>();
        try (InputStream raw = Files.newInputStream(path);
             InputStream in = path.getFileName().toString().endsWith(".gz") ? new GZIPInputStream(raw) : raw;
             BufferedReader reader = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.replaceAll("[\\r\\n]+$", "");
                if (line.isEmpty() || line.startsWith("#")) {
                    continue;
                }
                String[] fields = line.split("\t", -1);
                if (fields.length != 9) {
                    continue;
                }
                features.add(new GffFeature(
                        fields[0],
                        fields[1],
                        fields[2],
                        Integer.parseInt(fields[3]),
                        Integer.parseInt(fields[4]),
                        fields[5],
                        fields[6],
                        fields[7],
                        parseAttributes(fields[8])
                ));
            }
        }
        return features;
    }

    public static Map<String, List<String>> parseAttributes(String attrString) {
        Map<String, List<String>> attrs = new LinkedHashMap<>();
        for (String rawPart : attrString.split(";")) {
            String part = rawPart.trim();
            if (part.isEmpty() || !part.contains("=")) {
                continue;
            }
            String[] kv = part.split("=", 2);
            attrs.computeIfAbsent(kv[0].trim(), k -> new ArrayList<>()).add(kv[1].trim());
        }
        return attrs;
    }

    private static List<String> attr(Map<String, List<String>> attrs, String key) {
        List<String> direct = attrs.get(key);
        if (direct != null) {
            return direct;
        }
        return attrs.getOrDefault(key.toLowerCase(Locale.ROOT), Collections.emptyList());
    }

    public static String normalizeAttributeId(List<String> values, String... prefixes) {
        List<String> ids = normalizeAttributeIds(values, prefixes);
        return ids.isEmpty() ? "" : ids.get(0);
    }

    public static List<String> normalizeAttributeIds(List<String> values, String... prefixes) {
        List<String> normalized = new ArrayList<>();
        Set<String> seen = new HashSet<>();
        if (values == null) {
            values = Collections.emptyList();
        }
        for (String raw : values) {
            for (String part : raw.split(",")) {
                String item = part.trim();
                if (item.isEmpty()) {
                    continue;
                }
                for (String prefix : prefixes) {
                    if (item.startsWith(prefix)) {
                        item = item.substring(prefix.length());
                        break;
                    }
                }
                if (!item.isEmpty() && seen.add(item)) {
                    normalized.add(item);
                }
            }
        }
        return normalized;
    }

    public static LinkedHashMap<String, GeneModel> parseGff3ToModels(Path path, String geneScope) throws IOException {
        LinkedHashMap<String, GeneModel> genes = new LinkedHashMap<>();
        LinkedHashMap<String, MRNAModel> mrnas = new LinkedHashMap<>();
        Map<String, List<GffFeature>> mrnaChildren = new LinkedHashMap<>();

        List<GffFeature> features = new ArrayList<>();
        Set<String> keepTypes = new HashSet<>(Arrays.asList(
                "gene", "mRNA", "transcript", "exon", "CDS", "five_prime_UTR", "three_prime_UTR"
        ));
        for (GffFeature feature : parseGff3(path)) {
            if (keepTypes.contains(feature.type)) {
                features.add(feature);
            }
        }

        for (GffFeature feature : features) {
            if ("gene".equals(feature.type)) {
                String gid = normalizeAttributeId(attr(feature.attributes, "ID"), "gene:");
                if (gid.isEmpty()) {
                    continue;
                }
                genes.put(gid, new GeneModel(
                        gid,
                        feature.seqid,
                        feature.start,
                        feature.end,
                        feature.strand,
                        feature.source,
                        feature.start,
                        feature.end
                ));
            }
        }

        int mrnaOrder = 0;
        for (GffFeature feature : features) {
            if ("mRNA".equals(feature.type) || "transcript".equals(feature.type)) {
                String mid = normalizeAttributeId(attr(feature.attributes, "ID"), "transcript:");
                if (mid.isEmpty()) {
                    continue;
                }
                MRNAModel model = new MRNAModel(mid, feature.start, feature.end, mrnaOrder++);
                mrnas.put(mid, model);
                for (String parent : normalizeAttributeIds(attr(feature.attributes, "Parent"), "gene:")) {
                    mrnaChildren.computeIfAbsent(parent, k -> new ArrayList<>()).add(feature);
                }
            }
        }

        Map<String, String> mrnaToGene = new HashMap<>();
        for (Map.Entry<String, List<GffFeature>> entry : mrnaChildren.entrySet()) {
            String gid = entry.getKey();
            GeneModel gene = genes.get(gid);
            if (gene == null) {
                continue;
            }
            for (GffFeature mf : entry.getValue()) {
                String mid = normalizeAttributeId(attr(mf.attributes, "ID"), "transcript:");
                MRNAModel mrna = mrnas.get(mid);
                if (mrna != null && !containsMrnaIdentity(gene.mrnas, mrna)) {
                    gene.mrnas.add(mrna);
                    mrnaToGene.put(mid, gid);
                }
            }
        }

        if (genes.values().stream().noneMatch(g -> !g.mrnas.isEmpty())) {
            for (String gid : genes.keySet()) {
                List<GffFeature> children = mrnaChildren.get(gid);
                if (children == null) {
                    continue;
                }
                GeneModel gene = genes.get(gid);
                for (GffFeature mf : children) {
                    String mid = normalizeAttributeId(attr(mf.attributes, "ID"), "transcript:");
                    MRNAModel mrna = mrnas.get(mid);
                    if (mrna != null && !containsMrnaIdentity(gene.mrnas, mrna)) {
                        gene.mrnas.add(mrna);
                        mrnaToGene.put(mid, gid);
                    }
                }
            }
        }

        for (GffFeature feature : features) {
            if (!("exon".equals(feature.type) || "CDS".equals(feature.type)
                    || "five_prime_UTR".equals(feature.type) || "three_prime_UTR".equals(feature.type))) {
                continue;
            }
            for (String mrnaId : normalizeAttributeIds(attr(feature.attributes, "Parent"), "transcript:")) {
                MRNAModel mrna = mrnas.get(mrnaId);
                if (mrna == null) {
                    continue;
                }
                if ("exon".equals(feature.type)) {
                    mrna.exons.add(new ExonFeature(feature.start, feature.end));
                    mrna.hasExplicitExon = true;
                } else if ("CDS".equals(feature.type)) {
                    mrna.cds.add(new CDSFeature(feature.start, feature.end, feature.phase));
                } else {
                    mrna.utrs.add(new UTRFeature(feature.start, feature.end, feature.type));
                }
            }
        }

        for (GeneModel gene : genes.values()) {
            for (MRNAModel mrna : gene.mrnas) {
                sortMrnaFeatures(mrna);
            }

            for (MRNAModel mrna : gene.mrnas) {
                deriveMissingUtrs(gene, mrna);
            }

            for (MRNAModel mrna : gene.mrnas) {
                deriveMissingExons(mrna);
            }

            for (MRNAModel mrna : gene.mrnas) {
                sortMrnaFeatures(mrna);
            }

            gene.mrnas.removeIf(mrna -> !(mrna.hasExplicitExon || !mrna.cds.isEmpty()));

            List<Integer> starts = new ArrayList<>();
            List<Integer> ends = new ArrayList<>();
            for (MRNAModel mrna : gene.mrnas) {
                if (!mrna.exons.isEmpty()) {
                    starts.add(mrna.exons.stream().mapToInt(e -> e.start).min().orElse(gene.start));
                    ends.add(mrna.exons.stream().mapToInt(e -> e.end).max().orElse(gene.end));
                } else if (!mrna.cds.isEmpty()) {
                    starts.add(mrna.cds.stream().mapToInt(c -> c.start).min().orElse(gene.start));
                    ends.add(mrna.cds.stream().mapToInt(c -> c.end).max().orElse(gene.end));
                } else if (!mrna.utrs.isEmpty()) {
                    starts.add(mrna.utrs.stream().mapToInt(u -> u.start).min().orElse(gene.start));
                    ends.add(mrna.utrs.stream().mapToInt(u -> u.end).max().orElse(gene.end));
                }
            }
            if (!starts.isEmpty()) {
                gene.start = starts.stream().mapToInt(Integer::intValue).min().orElse(gene.start);
                gene.end = ends.stream().mapToInt(Integer::intValue).max().orElse(gene.end);
            }
        }

        LinkedHashMap<String, GeneModel> filtered = new LinkedHashMap<>();
        for (Map.Entry<String, GeneModel> entry : genes.entrySet()) {
            GeneModel gene = entry.getValue();
            if (gene.mrnaCount() == 0) {
                continue;
            }
            boolean keep;
            if ("mrna".equals(geneScope)) {
                keep = gene.mrnaCount() > 0;
            } else if ("coding".equals(geneScope)) {
                keep = gene.totalCdsLength() > 0;
            } else if ("all".equals(geneScope)) {
                keep = true;
            } else {
                throw new IllegalArgumentException("Unsupported gene_scope: " + geneScope);
            }
            if (keep) {
                filtered.put(entry.getKey(), gene);
            }
        }
        return filtered;
    }

    private static boolean containsMrnaIdentity(List<MRNAModel> mrnas, MRNAModel target) {
        for (MRNAModel mrna : mrnas) {
            if (mrna == target) {
                return true;
            }
        }
        return false;
    }

    private static void sortMrnaFeatures(MRNAModel mrna) {
        mrna.exons.sort(Comparator.comparingInt((ExonFeature e) -> e.start).thenComparingInt(e -> e.end));
        mrna.cds.sort(Comparator.comparingInt((CDSFeature c) -> c.start).thenComparingInt(c -> c.end));
        mrna.utrs.sort(Comparator.comparingInt((UTRFeature u) -> u.start)
                .thenComparingInt(u -> u.end)
                .thenComparing(u -> u.utrType));
    }

    private static void deriveMissingUtrs(GeneModel gene, MRNAModel mrna) {
        if (mrna.cds.isEmpty()) {
            return;
        }
        int cdsStart = mrna.cds.stream().mapToInt(c -> c.start).min().orElse(0);
        int cdsEnd = mrna.cds.stream().mapToInt(c -> c.end).max().orElse(0);
        String beforeCdsType = "-".equals(gene.strand) ? "three_prime_UTR" : "five_prime_UTR";
        String afterCdsType = "-".equals(gene.strand) ? "five_prime_UTR" : "three_prime_UTR";

        Set<String> presentUtrTypes = new HashSet<>();
        for (UTRFeature u : mrna.utrs) {
            presentUtrTypes.add(u.utrType);
        }
        List<UTRFeature> missing = new ArrayList<>();
        if (mrna.exons.isEmpty()) {
            return;
        }
        for (ExonFeature e : mrna.exons) {
            if (e.end < cdsStart) {
                addMissingUtr(missing, presentUtrTypes, e.start, e.end, beforeCdsType);
            } else if (e.start > cdsEnd) {
                addMissingUtr(missing, presentUtrTypes, e.start, e.end, afterCdsType);
            } else {
                addMissingUtr(missing, presentUtrTypes, e.start, cdsStart - 1, beforeCdsType);
                addMissingUtr(missing, presentUtrTypes, cdsEnd + 1, e.end, afterCdsType);
            }
        }
        mrna.utrs.addAll(missing);
        sortMrnaFeatures(mrna);
    }

    private static void addMissingUtr(List<UTRFeature> missing, Set<String> presentUtrTypes,
                                      int start, int end, String type) {
        if (start <= end && !presentUtrTypes.contains(type)) {
            missing.add(new UTRFeature(start, end, type));
        }
    }

    private static void deriveMissingExons(MRNAModel mrna) {
        if (!mrna.exons.isEmpty()) {
            return;
        }
        List<int[]> intervals = new ArrayList<>();
        for (CDSFeature c : mrna.cds) {
            intervals.add(new int[]{c.start, c.end});
        }
        for (UTRFeature u : mrna.utrs) {
            intervals.add(new int[]{u.start, u.end});
        }
        if (intervals.isEmpty()) {
            return;
        }
        intervals.sort(Comparator.<int[]>comparingInt(a -> a[0]).thenComparingInt(a -> a[1]));
        List<int[]> merged = new ArrayList<>();
        merged.add(Arrays.copyOf(intervals.get(0), 2));
        for (int i = 1; i < intervals.size(); i++) {
            int[] current = intervals.get(i);
            int[] last = merged.get(merged.size() - 1);
            if (current[0] <= last[1] + 1) {
                last[1] = Math.max(last[1], current[1]);
            } else {
                merged.add(Arrays.copyOf(current, 2));
            }
        }
        for (int[] interval : merged) {
            mrna.exons.add(new ExonFeature(interval[0], interval[1]));
        }
        mrna.exonsAreInferred = true;
    }

    public static int geneBoundaryStart(GeneModel gene) {
        return gene.start;
    }

    public static int geneBoundaryEnd(GeneModel gene) {
        return gene.end;
    }

    public static int geneBoundaryLength(GeneModel gene) {
        return geneBoundaryEnd(gene) - geneBoundaryStart(gene) + 1;
    }

    public static int rawGeneBoundaryStart(GeneModel gene) {
        return gene.rawStart != null ? gene.rawStart : gene.start;
    }

    public static int rawGeneBoundaryEnd(GeneModel gene) {
        return gene.rawEnd != null ? gene.rawEnd : gene.end;
    }

    public static int rawGeneBoundaryLength(GeneModel gene) {
        return rawGeneBoundaryEnd(gene) - rawGeneBoundaryStart(gene) + 1;
    }

    public static boolean geneBoundaryChanged(GeneModel before, GeneModel after, int boundaryTol) {
        return Math.abs(geneBoundaryStart(before) - geneBoundaryStart(after)) > boundaryTol
                || Math.abs(geneBoundaryEnd(before) - geneBoundaryEnd(after)) > boundaryTol;
    }

    public static int overlapLen(int aStart, int aEnd, int bStart, int bEnd) {
        return Math.max(0, Math.min(aEnd, bEnd) - Math.max(aStart, bStart) + 1);
    }

    public static List<int[]> mergeIntervals(List<int[]> intervals) {
        List<int[]> valid = new ArrayList<>();
        for (int[] interval : intervals) {
            if (interval[0] > 0 && interval[1] > 0 && interval[0] <= interval[1]) {
                valid.add(Arrays.copyOf(interval, 2));
            }
        }
        valid.sort(Comparator.<int[]>comparingInt(a -> a[0]).thenComparingInt(a -> a[1]));
        if (valid.isEmpty()) {
            return valid;
        }
        List<int[]> merged = new ArrayList<>();
        merged.add(valid.get(0));
        for (int i = 1; i < valid.size(); i++) {
            int[] current = valid.get(i);
            int[] last = merged.get(merged.size() - 1);
            if (current[0] <= last[1] + 1) {
                last[1] = Math.max(last[1], current[1]);
            } else {
                merged.add(current);
            }
        }
        return merged;
    }

    public static int intervalTotalLength(List<int[]> intervals) {
        int total = 0;
        for (int[] interval : mergeIntervals(intervals)) {
            total += interval[1] - interval[0] + 1;
        }
        return total;
    }

    public static int intervalOverlapLength(List<int[]> leftIntervals, List<int[]> rightIntervals) {
        List<int[]> left = mergeIntervals(leftIntervals);
        List<int[]> right = mergeIntervals(rightIntervals);
        int total = 0;
        int i = 0;
        int j = 0;
        while (i < left.size() && j < right.size()) {
            int[] l = left.get(i);
            int[] r = right.get(j);
            total += overlapLen(l[0], l[1], r[0], r[1]);
            if (l[1] < r[1]) {
                i++;
            } else {
                j++;
            }
        }
        return total;
    }

    public static List<int[]> mrnaFeatureIntervals(MRNAModel mrna) {
        List<int[]> intervals = new ArrayList<>();
        if (!mrna.exons.isEmpty()) {
            for (ExonFeature exon : mrna.exons) {
                intervals.add(new int[]{exon.start, exon.end});
            }
            return intervals;
        }
        for (CDSFeature cds : mrna.cds) {
            intervals.add(new int[]{cds.start, cds.end});
        }
        for (UTRFeature utr : mrna.utrs) {
            intervals.add(new int[]{utr.start, utr.end});
        }
        return intervals;
    }

    public static int transcriptFeatureLength(MRNAModel mrna) {
        return intervalTotalLength(mrnaFeatureIntervals(mrna));
    }

    public static List<int[]> geneFeatureIntervals(GeneModel gene) {
        List<int[]> intervals = new ArrayList<>();
        for (MRNAModel mrna : gene.mrnas) {
            intervals.addAll(mrnaFeatureIntervals(mrna));
        }
        if (!intervals.isEmpty()) {
            return mergeIntervals(intervals);
        }
        return Collections.singletonList(new int[]{gene.start, gene.end});
    }

    public static List<int[]> mrnaCdsIntervals(MRNAModel mrna) {
        List<int[]> intervals = new ArrayList<>();
        for (CDSFeature cds : mrna.cds) {
            intervals.add(new int[]{cds.start, cds.end});
        }
        return intervals;
    }

    public static boolean geneLacksExplicitExons(GeneModel gene) {
        for (MRNAModel mrna : gene.mrnas) {
            if (!mrna.exons.isEmpty() && !mrna.exonsAreInferred) {
                return false;
            }
        }
        return true;
    }

    private static List<TranscriptIntervals> geneTranscriptIntervalSets(
            GeneModel gene, boolean cdsOnly, boolean fallbackToGeneSpan) {
        List<TranscriptIntervals> sets = new ArrayList<>();
        for (MRNAModel mrna : gene.mrnas) {
            List<int[]> intervals = cdsOnly ? mrnaCdsIntervals(mrna) : mrnaFeatureIntervals(mrna);
            intervals = mergeIntervals(intervals);
            if (!intervals.isEmpty()) {
                sets.add(new TranscriptIntervals(mrna, intervals));
            }
        }
        if (!sets.isEmpty() || !fallbackToGeneSpan) {
            return sets;
        }
        sets.add(new TranscriptIntervals(null, Collections.singletonList(new int[]{gene.start, gene.end})));
        return sets;
    }

    public static TranscriptOverlapMetrics bestTranscriptOverlapMetrics(GeneModel a, GeneModel b) {
        return bestTranscriptOverlapMetrics(a, b, false, true);
    }

    public static TranscriptOverlapMetrics bestTranscriptOverlapMetrics(
            GeneModel a, GeneModel b, boolean cdsOnly, boolean fallbackToGeneSpan) {
        TranscriptOverlapMetrics best = new TranscriptOverlapMetrics(0.0, 0.0, 0, 0, 0, null, null);
        for (TranscriptIntervals left : geneTranscriptIntervalSets(a, cdsOnly, fallbackToGeneSpan)) {
            int leftLength = intervalTotalLength(left.intervals);
            if (leftLength == 0) {
                continue;
            }
            for (TranscriptIntervals right : geneTranscriptIntervalSets(b, cdsOnly, fallbackToGeneSpan)) {
                int rightLength = intervalTotalLength(right.intervals);
                if (rightLength == 0) {
                    continue;
                }
                int overlap = intervalOverlapLength(left.intervals, right.intervals);
                if (overlap == 0) {
                    continue;
                }
                double score = overlap / (double) Math.min(leftLength, rightLength);
                int union = leftLength + rightLength - overlap;
                double jaccard = union == 0 ? 0.0 : overlap / (double) union;
                TranscriptOverlapMetrics candidate = new TranscriptOverlapMetrics(
                        score, jaccard, overlap, leftLength, rightLength, left.mrna, right.mrna);
                if (compareTranscriptOverlapMetrics(candidate, best) > 0) {
                    best = candidate;
                }
            }
        }
        return best;
    }

    private static int compareTranscriptOverlapMetrics(TranscriptOverlapMetrics a, TranscriptOverlapMetrics b) {
        int c;
        c = Double.compare(a.score, b.score);
        if (c != 0) return c;
        c = Double.compare(a.jaccard, b.jaccard);
        if (c != 0) return c;
        c = Integer.compare(a.overlap, b.overlap);
        if (c != 0) return c;
        c = Integer.compare(-Math.abs(a.leftLength - a.rightLength), -Math.abs(b.leftLength - b.rightLength));
        if (c != 0) return c;
        c = mrnaId(a.leftMrna).compareTo(mrnaId(b.leftMrna));
        if (c != 0) return c;
        return mrnaId(a.rightMrna).compareTo(mrnaId(b.rightMrna));
    }

    private static String mrnaId(MRNAModel mrna) {
        return mrna == null ? "" : mrna.mrnaId;
    }

    public static double cdsCompatibleReciprocalOverlap(GeneModel a, GeneModel b) {
        if (!(geneLacksExplicitExons(a) || geneLacksExplicitExons(b))) {
            return 0.0;
        }
        return bestTranscriptOverlapMetrics(a, b, true, false).score;
    }

    public static double cdsOverlapRatio(GeneModel a, GeneModel b) {
        return bestTranscriptOverlapMetrics(a, b, true, false).score;
    }

    public static int anyFeatureOverlapLen(GeneModel a, GeneModel b) {
        return intervalOverlapLength(geneFeatureIntervals(a), geneFeatureIntervals(b));
    }

    public static int featureOverlapLen(GeneModel a, GeneModel b) {
        return bestTranscriptOverlapMetrics(a, b).overlap;
    }

    public static double reciprocalOverlap(GeneModel a, GeneModel b) {
        return Math.max(bestTranscriptOverlapMetrics(a, b).score, cdsCompatibleReciprocalOverlap(a, b));
    }

    public static double containmentOverlap(GeneModel a, GeneModel b) {
        return bestTranscriptOverlapMetrics(a, b).score;
    }

    public static double jaccardOverlap(GeneModel a, GeneModel b) {
        return bestTranscriptOverlapMetrics(a, b).jaccard;
    }

    public static double containmentRatio(GeneModel a, GeneModel b) {
        double best = 0.0;
        for (TranscriptIntervals left : geneTranscriptIntervalSets(a, false, true)) {
            int leftLength = intervalTotalLength(left.intervals);
            if (leftLength == 0) {
                continue;
            }
            for (TranscriptIntervals right : geneTranscriptIntervalSets(b, false, true)) {
                int overlap = intervalOverlapLength(left.intervals, right.intervals);
                if (overlap > 0) {
                    best = Math.max(best, overlap / (double) leftLength);
                }
            }
        }
        return best;
    }

    public static List<Pair> findOverlappingPairs(List<GeneModel> beforeGenes, List<GeneModel> afterGenes,
                                                  double minReciprocal, String overlapMode,
                                                  Map<String, Integer> diagnostics) {
        return findOverlappingPairs(beforeGenes, afterGenes, minReciprocal, overlapMode, diagnostics, null);
    }

    public static List<Pair> findOverlappingPairs(List<GeneModel> beforeGenes, List<GeneModel> afterGenes,
                                                  double minReciprocal, String overlapMode,
                                                  Map<String, Integer> diagnostics,
                                                  Map<String, Set<String>> weakRejected) {
        if (!"reciprocal".equals(overlapMode) && !"containment".equals(overlapMode)) {
            throw new IllegalArgumentException("Unsupported overlap_mode: " + overlapMode);
        }
        List<Pair> pairs = new ArrayList<>();
        int i = 0;
        int j = 0;
        List<GeneModel> activeBefore = new ArrayList<>();
        while (j < afterGenes.size()) {
            GeneModel after = afterGenes.get(j);
            activeBefore.removeIf(bg -> bg.end < after.start);
            while (i < beforeGenes.size() && beforeGenes.get(i).start <= after.end) {
                activeBefore.add(beforeGenes.get(i));
                i++;
            }
            for (GeneModel before : activeBefore) {
                if (!before.strand.equals(after.strand)) {
                    continue;
                }
                if (featureOverlapLen(before, after) == 0) {
                    continue;
                }
                double ro = reciprocalOverlap(before, after);
                double co = containmentOverlap(before, after);
                if (diagnostics != null) {
                    diagnostics.merge("same_strand_overlaps", 1, Integer::sum);
                }
                if (co >= minReciprocal && ro < minReciprocal) {
                    if (diagnostics != null) {
                        diagnostics.merge("containment_pairs_filtered_by_reciprocal", 1, Integer::sum);
                    }
                    if (weakRejected != null) {
                        weakRejected.computeIfAbsent("before", k -> new HashSet<>()).add(before.geneId);
                        weakRejected.computeIfAbsent("after", k -> new HashSet<>()).add(after.geneId);
                    }
                }
                double score = "reciprocal".equals(overlapMode) ? ro : co;
                if (score >= minReciprocal) {
                    pairs.add(new Pair(before, after, score, jaccardOverlap(before, after)));
                }
            }
            j++;
        }
        return pairs;
    }

    public static LinkedHashMap<String, Integer> countNoOverlapLoci(LinkedHashMap<String, GeneModel> beforeGenes,
                                                                    LinkedHashMap<String, GeneModel> afterGenes) {
        Map<String, List<GeneModel>> beforeBySeq = groupBySeqid(beforeGenes);
        Map<String, List<GeneModel>> afterBySeq = groupBySeqid(afterGenes);
        LinkedHashMap<String, Integer> counts = new LinkedHashMap<>();
        counts.put("no_overlap_after_loci", countWithoutOverlap(afterBySeq, beforeBySeq));
        counts.put("no_overlap_before_loci", countWithoutOverlap(beforeBySeq, afterBySeq));
        return counts;
    }

    private static PruneResult pruneContainmentBridgeEdges(List<Pair> containmentPairs, List<Pair> reciprocalPairs) {
        List<Pair> allPairs = new ArrayList<>();
        Set<String> seen = new HashSet<>();
        for (Pair pair : concatPairs(containmentPairs, reciprocalPairs)) {
            String key = pair.before.geneId + "\u0000" + pair.after.geneId;
            if (seen.add(key)) {
                allPairs.add(pair);
            }
        }

        Map<String, Pair> bestAfterByBefore = new HashMap<>();
        Map<String, Pair> bestBeforeByAfter = new HashMap<>();
        for (Pair pair : allPairs) {
            Pair bestAfter = bestAfterByBefore.get(pair.before.geneId);
            if (bestAfter == null || comparePairQuality(pair, bestAfter) > 0) {
                bestAfterByBefore.put(pair.before.geneId, pair);
            }
            Pair bestBefore = bestBeforeByAfter.get(pair.after.geneId);
            if (bestBefore == null || comparePairQuality(pair, bestBefore) > 0) {
                bestBeforeByAfter.put(pair.after.geneId, pair);
            }
        }

        Map<String, String> mutualAnchorAfterByBefore = new HashMap<>();
        Map<String, String> mutualAnchorBeforeByAfter = new HashMap<>();
        for (Pair pair : allPairs) {
            if (bestAfterByBefore.get(pair.before.geneId) == pair
                    && bestBeforeByAfter.get(pair.after.geneId) == pair
                    && (cdsOverlapRatio(pair.before, pair.after) >= 0.5
                    || pair.jaccard >= 0.5
                    || pair.score >= 0.95)) {
                mutualAnchorAfterByBefore.put(pair.before.geneId, pair.after.geneId);
                mutualAnchorBeforeByAfter.put(pair.after.geneId, pair.before.geneId);
            }
        }

        List<Pair> pruned = new ArrayList<>();
        int prunedCount = 0;
        for (Pair pair : containmentPairs) {
            String beforeAnchor = mutualAnchorAfterByBefore.get(pair.before.geneId);
            String afterAnchor = mutualAnchorBeforeByAfter.get(pair.after.geneId);
            boolean bridgesTwoIndependentAnchors = beforeAnchor != null
                    && afterAnchor != null
                    && !beforeAnchor.equals(pair.after.geneId)
                    && !afterAnchor.equals(pair.before.geneId);
            boolean lowSupportBridge = cdsOverlapRatio(pair.before, pair.after) < 0.5
                    && pair.jaccard < 0.25;
            if (bridgesTwoIndependentAnchors && lowSupportBridge) {
                prunedCount++;
                continue;
            }
            pruned.add(pair);
        }
        return new PruneResult(pruned, prunedCount);
    }

    private static List<Pair> concatPairs(List<Pair> left, List<Pair> right) {
        List<Pair> out = new ArrayList<>(left.size() + right.size());
        out.addAll(left);
        out.addAll(right);
        return out;
    }

    private static int comparePairQuality(Pair a, Pair b) {
        int c;
        double aCds = cdsOverlapRatio(a.before, a.after);
        double bCds = cdsOverlapRatio(b.before, b.after);
        c = Boolean.compare(aCds >= 0.5, bCds >= 0.5);
        if (c != 0) return c;
        c = Double.compare(aCds, bCds);
        if (c != 0) return c;
        c = Double.compare(a.jaccard, b.jaccard);
        if (c != 0) return c;
        return Double.compare(a.score, b.score);
    }

    private static Map<String, List<GeneModel>> groupBySeqid(LinkedHashMap<String, GeneModel> genes) {
        Map<String, List<GeneModel>> groups = new LinkedHashMap<>();
        for (GeneModel gene : genes.values()) {
            groups.computeIfAbsent(gene.seqid, k -> new ArrayList<>()).add(gene);
        }
        for (List<GeneModel> group : groups.values()) {
            group.sort(Comparator.comparingInt((GeneModel g) -> g.start)
                    .thenComparingInt(g -> g.end)
                    .thenComparing(g -> g.geneId));
        }
        return groups;
    }

    private static int countWithoutOverlap(Map<String, List<GeneModel>> queryBySeq,
                                           Map<String, List<GeneModel>> targetBySeq) {
        int count = 0;
        for (Map.Entry<String, List<GeneModel>> entry : queryBySeq.entrySet()) {
            List<GeneModel> targetGenes = targetBySeq.getOrDefault(entry.getKey(), Collections.emptyList());
            List<GeneModel> active = new ArrayList<>();
            int i = 0;
            for (GeneModel query : entry.getValue()) {
                active.removeIf(target -> target.end < query.start);
                while (i < targetGenes.size() && targetGenes.get(i).start <= query.end) {
                    active.add(targetGenes.get(i));
                    i++;
                }
                boolean hasOverlap = false;
                for (GeneModel target : active) {
                    if (anyFeatureOverlapLen(query, target) > 0) {
                        hasOverlap = true;
                        break;
                    }
                }
                if (!hasOverlap) {
                    count++;
                }
            }
        }
        return count;
    }

    public static ResolvedMatches resolveMatches(List<Pair> pairs,
                                                 LinkedHashMap<String, GeneModel> beforeGenes,
                                                 LinkedHashMap<String, GeneModel> afterGenes) {
        Map<Node, Set<Node>> adj = new HashMap<>();
        Set<Node> beforeNodes = new LinkedHashSet<>();
        Set<Node> afterNodes = new LinkedHashSet<>();
        Set<String> directPairIds = new HashSet<>();
        for (Pair pair : pairs) {
            Node before = new Node("before", pair.before.geneId);
            Node after = new Node("after", pair.after.geneId);
            adj.computeIfAbsent(before, k -> new LinkedHashSet<>()).add(after);
            adj.computeIfAbsent(after, k -> new LinkedHashSet<>()).add(before);
            beforeNodes.add(before);
            afterNodes.add(after);
            directPairIds.add(pair.before.geneId + "\u0000" + pair.after.geneId);
        }

        Set<Node> allNodes = new LinkedHashSet<>();
        allNodes.addAll(beforeNodes);
        allNodes.addAll(afterNodes);
        Set<Node> visited = new HashSet<>();
        List<List<List<String>>> loci = new ArrayList<>();

        for (Node node : allNodes) {
            if (visited.contains(node)) {
                continue;
            }
            Set<Node> component = new HashSet<>();
            ArrayDeque<Node> stack = new ArrayDeque<>();
            stack.push(node);
            while (!stack.isEmpty()) {
                Node current = stack.pop();
                if (!visited.add(current)) {
                    continue;
                }
                component.add(current);
                for (Node neighbor : adj.getOrDefault(current, Collections.emptySet())) {
                    if (!visited.contains(neighbor)) {
                        stack.push(neighbor);
                    }
                }
            }
            List<String> compBefore = new ArrayList<>();
            List<String> compAfter = new ArrayList<>();
            for (Node current : component) {
                if ("before".equals(current.side)) {
                    compBefore.add(current.id);
                } else {
                    compAfter.add(current.id);
                }
            }
            Collections.sort(compBefore);
            Collections.sort(compAfter);
            if (!compBefore.isEmpty() || !compAfter.isEmpty()) {
                loci.add(Arrays.asList(compBefore, compAfter));
            }
        }

        ResolvedMatches resolved = new ResolvedMatches();
        Set<String> matchedBefore = new HashSet<>();
        Set<String> matchedAfter = new HashSet<>();
        for (List<List<String>> locus : loci) {
            List<String> compBefore = locus.get(0);
            List<String> compAfter = locus.get(1);
            int nb = compBefore.size();
            int na = compAfter.size();
            if (nb == 1 && na == 1) {
                GeneModel before = beforeGenes.get(compBefore.get(0));
                GeneModel after = afterGenes.get(compAfter.get(0));
                resolved.syntenic.add(new GenePair(before, after));
                matchedBefore.add(compBefore.get(0));
                matchedAfter.add(compAfter.get(0));
            } else if (nb == 1 && na >= 2) {
                GeneModel before = beforeGenes.get(compBefore.get(0));
                List<GeneModel> afters = new ArrayList<>();
                for (String id : compAfter) {
                    afters.add(afterGenes.get(id));
                    matchedAfter.add(id);
                }
                resolved.split.add(new SplitEvent(before, afters));
                matchedBefore.add(compBefore.get(0));
            } else if (nb >= 2 && na == 1) {
                List<GeneModel> befores = new ArrayList<>();
                for (String id : compBefore) {
                    befores.add(beforeGenes.get(id));
                    matchedBefore.add(id);
                }
                GeneModel after = afterGenes.get(compAfter.get(0));
                resolved.merge.add(new MergeEvent(befores, after));
                matchedAfter.add(compAfter.get(0));
            } else if (nb >= 2 && na >= 2) {
                List<GeneModel> befores = new ArrayList<>();
                for (String id : compBefore) {
                    befores.add(beforeGenes.get(id));
                    matchedBefore.add(id);
                }
                List<GeneModel> afters = new ArrayList<>();
                for (String id : compAfter) {
                    afters.add(afterGenes.get(id));
                    matchedAfter.add(id);
                }
                List<GenePair> edgePairs = new ArrayList<>();
                for (String beforeId : compBefore) {
                    for (String afterId : compAfter) {
                        if (directPairIds.contains(beforeId + "\u0000" + afterId)) {
                            edgePairs.add(new GenePair(beforeGenes.get(beforeId), afterGenes.get(afterId)));
                        }
                    }
                }
                resolved.complex.add(new ComplexEvent(befores, afters, edgePairs));
            }
        }

        for (Map.Entry<String, GeneModel> entry : afterGenes.entrySet()) {
            if (!matchedAfter.contains(entry.getKey())) {
                resolved.novel.add(entry.getValue());
            }
        }
        for (Map.Entry<String, GeneModel> entry : beforeGenes.entrySet()) {
            if (!matchedBefore.contains(entry.getKey())) {
                resolved.deleted.add(entry.getValue());
            }
        }
        return resolved;
    }

    public static String exonSignature(MRNAModel mrna) {
        List<String> parts = new ArrayList<>();
        for (ExonFeature e : mrna.exons) {
            parts.add(e.start + ":" + e.end);
        }
        Collections.sort(parts);
        return String.join("|", parts);
    }

    public static String cdsSignature(MRNAModel mrna) {
        List<String> parts = new ArrayList<>();
        for (CDSFeature c : mrna.cds) {
            // Phase-only differences are ignored for structural comparison.
            parts.add(c.start + ":" + c.end);
        }
        Collections.sort(parts);
        return String.join("|", parts);
    }

    public static String utrSignature(MRNAModel mrna) {
        List<String> parts = new ArrayList<>();
        for (UTRFeature u : mrna.utrs) {
            parts.add(u.start + ":" + u.end + ":" + u.utrType);
        }
        Collections.sort(parts);
        return String.join("|", parts);
    }

    public static String transcriptSignature(MRNAModel mrna) {
        return "(" + exonSignature(mrna) + ")(" + cdsSignature(mrna) + ")(" + utrSignature(mrna) + ")";
    }

    public static String geneStructureSignature(GeneModel gene) {
        List<String> parts = new ArrayList<>();
        for (MRNAModel mrna : gene.mrnas) {
            parts.add(transcriptSignature(mrna));
        }
        Collections.sort(parts);
        return String.join("||", parts);
    }

    public static boolean exonOverlapsCds(ExonFeature exon, List<CDSFeature> cdsFeatures) {
        for (CDSFeature cds : cdsFeatures) {
            if (overlapLen(exon.start, exon.end, cds.start, cds.end) > 0) {
                return true;
            }
        }
        return false;
    }

    public static int codingExonCount(MRNAModel mrna) {
        if (!mrna.exons.isEmpty()) {
            int count = 0;
            for (ExonFeature exon : mrna.exons) {
                if (exonOverlapsCds(exon, mrna.cds)) {
                    count++;
                }
            }
            return count;
        }
        return mrna.cds.size();
    }

    public static int utrOnlyExonCount(MRNAModel mrna) {
        if (!mrna.exons.isEmpty()) {
            int count = 0;
            for (ExonFeature exon : mrna.exons) {
                if (!exonOverlapsCds(exon, mrna.cds)) {
                    count++;
                }
            }
            return count;
        }
        return !mrna.utrs.isEmpty() ? mrna.utrs.size() : 0;
    }

    public static MRNAModel[] selectRepresentativeMrnaPair(GeneModel beforeGene, GeneModel afterGene) {
        if (beforeGene.mrnas.isEmpty() || afterGene.mrnas.isEmpty()) {
            return new MRNAModel[]{beforeGene.primaryMrna(), afterGene.primaryMrna()};
        }
        MRNAModel before = beforeGene.primaryMrna();
        MRNAModel best = null;
        TranscriptScore bestScore = null;
        for (MRNAModel candidate : afterGene.mrnas) {
            TranscriptScore score = transcriptMatchScore(before, candidate);
            if (best == null || score.compareTo(bestScore) > 0) {
                best = candidate;
                bestScore = score;
            }
        }
        return new MRNAModel[]{before, best};
    }

    public static LinkedHashMap<String, Integer> summarizeRepresentativeTranscriptChanges(List<GenePair> syntenicPairs) {
        LinkedHashMap<String, Integer> summary = new LinkedHashMap<>();
        summary.put("rep_transcript_pairs", 0);
        summary.put("rep_structural_changed", 0);
        summary.put("rep_exon_count_changed", 0);
        summary.put("rep_exon_boundary_changed_same_count", 0);
        summary.put("rep_cds_count_changed", 0);
        summary.put("rep_cds_boundary_changed_same_count", 0);
        for (GenePair pair : syntenicPairs) {
            MRNAModel[] mrnas = selectRepresentativeMrnaPair(pair.before, pair.after);
            MRNAModel before = mrnas[0];
            MRNAModel after = mrnas[1];
            if (before == null || after == null) {
                continue;
            }
            summary.merge("rep_transcript_pairs", 1, Integer::sum);
            boolean exonChanged = false;
            boolean cdsChanged = false;
            if (before.exonCount() != after.exonCount()) {
                summary.merge("rep_exon_count_changed", 1, Integer::sum);
                exonChanged = true;
            } else if (!exonSignature(before).equals(exonSignature(after))) {
                summary.merge("rep_exon_boundary_changed_same_count", 1, Integer::sum);
                exonChanged = true;
            }
            if (before.cds.size() != after.cds.size()) {
                summary.merge("rep_cds_count_changed", 1, Integer::sum);
                cdsChanged = true;
            } else if (!cdsSignature(before).equals(cdsSignature(after))) {
                summary.merge("rep_cds_boundary_changed_same_count", 1, Integer::sum);
                cdsChanged = true;
            }
            if (exonChanged || cdsChanged) {
                summary.merge("rep_structural_changed", 1, Integer::sum);
            }
        }
        return summary;
    }

    private static TranscriptScore transcriptMatchScore(MRNAModel before, MRNAModel after) {
        String bExons = exonSignature(before);
        String aExons = exonSignature(after);
        String bCds = cdsSignature(before);
        String aCds = cdsSignature(after);
        String bUtrs = utrSignature(before);
        String aUtrs = utrSignature(after);
        return new TranscriptScore(
                bExons.equals(aExons) && bCds.equals(aCds) && bUtrs.equals(aUtrs),
                bCds.equals(aCds),
                bExons.equals(aExons),
                bUtrs.equals(aUtrs),
                intersectionSize(signatureSet(bCds), signatureSet(aCds)),
                intersectionSize(signatureSet(bExons), signatureSet(aExons)),
                intersectionSize(signatureSet(bUtrs), signatureSet(aUtrs)),
                -Math.abs(before.cdsLength() - after.cdsLength()),
                -Math.abs(before.exonCount() - after.exonCount()),
                after.cdsLength(),
                after.exonCount(),
                after.mrnaId
        );
    }

    private static Set<String> signatureSet(String signature) {
        if (signature.isEmpty()) {
            return Collections.emptySet();
        }
        return new HashSet<>(Arrays.asList(signature.split("\\|", -1)));
    }

    private static int intersectionSize(Set<String> a, Set<String> b) {
        int count = 0;
        for (String item : a) {
            if (b.contains(item)) {
                count++;
            }
        }
        return count;
    }

    private static final class TranscriptScore implements Comparable<TranscriptScore> {
        final boolean exact;
        final boolean cdsEqual;
        final boolean exonsEqual;
        final boolean utrsEqual;
        final int commonCds;
        final int commonExons;
        final int commonUtrs;
        final int cdsLengthDeltaNeg;
        final int exonCountDeltaNeg;
        final int afterCdsLength;
        final int afterExonCount;
        final String afterMrnaId;

        TranscriptScore(boolean exact, boolean cdsEqual, boolean exonsEqual, boolean utrsEqual,
                        int commonCds, int commonExons, int commonUtrs,
                        int cdsLengthDeltaNeg, int exonCountDeltaNeg,
                        int afterCdsLength, int afterExonCount, String afterMrnaId) {
            this.exact = exact;
            this.cdsEqual = cdsEqual;
            this.exonsEqual = exonsEqual;
            this.utrsEqual = utrsEqual;
            this.commonCds = commonCds;
            this.commonExons = commonExons;
            this.commonUtrs = commonUtrs;
            this.cdsLengthDeltaNeg = cdsLengthDeltaNeg;
            this.exonCountDeltaNeg = exonCountDeltaNeg;
            this.afterCdsLength = afterCdsLength;
            this.afterExonCount = afterExonCount;
            this.afterMrnaId = afterMrnaId;
        }

        @Override
        public int compareTo(TranscriptScore other) {
            int c;
            c = Boolean.compare(exact, other.exact);
            if (c != 0) return c;
            c = Boolean.compare(cdsEqual, other.cdsEqual);
            if (c != 0) return c;
            c = Boolean.compare(exonsEqual, other.exonsEqual);
            if (c != 0) return c;
            c = Boolean.compare(utrsEqual, other.utrsEqual);
            if (c != 0) return c;
            c = Integer.compare(commonCds, other.commonCds);
            if (c != 0) return c;
            c = Integer.compare(commonExons, other.commonExons);
            if (c != 0) return c;
            c = Integer.compare(commonUtrs, other.commonUtrs);
            if (c != 0) return c;
            c = Integer.compare(cdsLengthDeltaNeg, other.cdsLengthDeltaNeg);
            if (c != 0) return c;
            c = Integer.compare(exonCountDeltaNeg, other.exonCountDeltaNeg);
            if (c != 0) return c;
            c = Integer.compare(afterCdsLength, other.afterCdsLength);
            if (c != 0) return c;
            c = Integer.compare(afterExonCount, other.afterExonCount);
            if (c != 0) return c;
            return afterMrnaId.compareTo(other.afterMrnaId);
        }
    }

    public static LinkedHashMap<String, Boolean> computeSyntenicAttributes(GeneModel before, GeneModel after,
                                                                           int boundaryTol, double cdsChangePct,
                                                                           double utrChangePct) {
        LinkedHashMap<String, Boolean> attrs = new LinkedHashMap<>();
        for (String key : SYNTENIC_ATTRIBUTE_KEYS) {
            attrs.put(key, false);
        }
        attrs.put("gene_boundary_changed", geneBoundaryChanged(before, after, boundaryTol));
        MRNAModel[] pair = selectRepresentativeMrnaPair(before, after);
        MRNAModel bm = pair[0];
        MRNAModel am = pair[1];
        if (bm == null || am == null) {
            attrs.put("exact", !attrs.get("gene_boundary_changed"));
            return attrs;
        }

        boolean exactBoundaries = !attrs.get("gene_boundary_changed");
        int exonDiff = am.exonCount() - bm.exonCount();
        int codingExonDiff = codingExonCount(am) - codingExonCount(bm);
        int utrOnlyExonDiff = utrOnlyExonCount(am) - utrOnlyExonCount(bm);
        boolean primaryExonCoordsChanged = !exonSignature(bm).equals(exonSignature(am));

        int bCds = bm.cdsLength();
        int aCds = am.cdsLength();
        boolean primaryCdsCoordsChanged = !cdsSignature(bm).equals(cdsSignature(am));
        double cdsChange = bCds > 0 ? Math.abs(aCds - bCds) / (double) bCds : (aCds > 0 ? 1.0 : 0.0);
        boolean cdsChanged = cdsChange >= cdsChangePct || primaryCdsCoordsChanged;

        boolean bUtr5 = bm.hasFivePrimeUtr();
        boolean bUtr3 = bm.hasThreePrimeUtr();
        boolean aUtr5 = am.hasFivePrimeUtr();
        boolean aUtr3 = am.hasThreePrimeUtr();
        attrs.put("utr_added", (!bUtr5 && aUtr5) || (!bUtr3 && aUtr3));
        attrs.put("utr_lost", (bUtr5 && !aUtr5) || (bUtr3 && !aUtr3));
        boolean primaryUtrCoordsChanged = !utrSignature(bm).equals(utrSignature(am));

        boolean utr5Refined = false;
        boolean utr3Refined = false;
        if (bUtr5 && aUtr5 && bm.fivePrimeUtrLength() > 0) {
            int delta = Math.abs(am.fivePrimeUtrLength() - bm.fivePrimeUtrLength());
            utr5Refined = delta / (double) bm.fivePrimeUtrLength() >= utrChangePct;
        }
        if (bUtr3 && aUtr3 && bm.threePrimeUtrLength() > 0) {
            int delta = Math.abs(am.threePrimeUtrLength() - bm.threePrimeUtrLength());
            utr3Refined = delta / (double) bm.threePrimeUtrLength() >= utrChangePct;
        }
        attrs.put("utr_refined", utr5Refined || utr3Refined
                || (primaryUtrCoordsChanged && !(attrs.get("utr_added") || attrs.get("utr_lost"))));
        attrs.put("coding_exon_gain", codingExonDiff > 0);
        attrs.put("coding_exon_loss", codingExonDiff < 0);
        attrs.put("utr_exon_added", utrOnlyExonDiff > 0);
        attrs.put("utr_exon_removed", utrOnlyExonDiff < 0);
        attrs.put("exon_boundary_refined",
                codingExonDiff == 0 && utrOnlyExonDiff == 0 && exonDiff == 0 && primaryExonCoordsChanged);
        attrs.put("cds_change", cdsChanged);
        attrs.put("cds_boundary_refined", cdsChanged && aCds == bCds);

        boolean mrnaCountChanged = before.mrnaCount() != after.mrnaCount();
        int totalExonDiff = Math.abs(after.totalExonCount() - before.totalExonCount());
        int bTotalCds = before.totalCdsLength();
        int aTotalCds = after.totalCdsLength();
        double totalCdsChange = bTotalCds > 0
                ? Math.abs(aTotalCds - bTotalCds) / (double) bTotalCds
                : (aTotalCds > 0 ? 1.0 : 0.0);
        boolean isoformRestructured =
                (totalExonDiff > 0 && totalExonDiff != Math.abs(exonDiff))
                        || (totalCdsChange >= cdsChangePct && Math.abs(totalCdsChange - cdsChange) > 0.001)
                        || (!geneStructureSignature(before).equals(geneStructureSignature(after))
                        && !(primaryExonCoordsChanged || primaryCdsCoordsChanged || primaryUtrCoordsChanged));
        attrs.put("isoform_change", mrnaCountChanged || isoformRestructured);
        attrs.put("exact", exactBoundaries
                && geneStructureSignature(before).equals(geneStructureSignature(after))
                && !attrs.get("utr_added")
                && !attrs.get("utr_lost")
                && !attrs.get("utr_refined")
                && !attrs.get("isoform_change"));
        return attrs;
    }

    public static String classifySyntenicChange(GeneModel before, GeneModel after,
                                                int boundaryTol, double cdsChangePct, double utrChangePct) {
        MRNAModel[] pair = selectRepresentativeMrnaPair(before, after);
        MRNAModel bm = pair[0];
        MRNAModel am = pair[1];
        if (bm == null || am == null) {
            return !geneBoundaryChanged(before, after, boundaryTol) ? "exact" : "boundary_refined";
        }

        boolean exactBoundaries = !geneBoundaryChanged(before, after, boundaryTol);
        int exonDiff = am.exonCount() - bm.exonCount();
        int codingExonDiff = codingExonCount(am) - codingExonCount(bm);
        int utrOnlyExonDiff = utrOnlyExonCount(am) - utrOnlyExonCount(bm);
        boolean primaryExonCoordsChanged = !exonSignature(bm).equals(exonSignature(am));

        int bCds = bm.cdsLength();
        int aCds = am.cdsLength();
        boolean primaryCdsCoordsChanged = !cdsSignature(bm).equals(cdsSignature(am));
        double cdsChange = bCds > 0 ? Math.abs(aCds - bCds) / (double) bCds : (aCds > 0 ? 1.0 : 0.0);

        boolean bUtr5 = bm.hasFivePrimeUtr();
        boolean bUtr3 = bm.hasThreePrimeUtr();
        boolean aUtr5 = am.hasFivePrimeUtr();
        boolean aUtr3 = am.hasThreePrimeUtr();
        boolean utr5Gained = !bUtr5 && aUtr5;
        boolean utr5Lost = bUtr5 && !aUtr5;
        boolean utr3Gained = !bUtr3 && aUtr3;
        boolean utr3Lost = bUtr3 && !aUtr3;
        boolean anyUtrGained = utr5Gained || utr3Gained;
        boolean anyUtrLost = utr5Lost || utr3Lost;
        boolean anyUtrTypeChanged = anyUtrGained || anyUtrLost;
        boolean primaryUtrCoordsChanged = !utrSignature(bm).equals(utrSignature(am));

        boolean utr5Refined = false;
        boolean utr3Refined = false;
        if (bUtr5 && aUtr5 && bm.fivePrimeUtrLength() > 0) {
            int delta = Math.abs(am.fivePrimeUtrLength() - bm.fivePrimeUtrLength());
            utr5Refined = delta / (double) bm.fivePrimeUtrLength() >= utrChangePct;
        }
        if (bUtr3 && aUtr3 && bm.threePrimeUtrLength() > 0) {
            int delta = Math.abs(am.threePrimeUtrLength() - bm.threePrimeUtrLength());
            utr3Refined = delta / (double) bm.threePrimeUtrLength() >= utrChangePct;
        }
        boolean utrBoundaryRefined = utr5Refined || utr3Refined
                || (primaryUtrCoordsChanged && !anyUtrTypeChanged);

        boolean mrnaCountChanged = before.mrnaCount() != after.mrnaCount();
        int totalExonDiff = Math.abs(after.totalExonCount() - before.totalExonCount());
        int bTotalCds = before.totalCdsLength();
        int aTotalCds = after.totalCdsLength();
        double totalCdsChange = bTotalCds > 0
                ? Math.abs(aTotalCds - bTotalCds) / (double) bTotalCds
                : (aTotalCds > 0 ? 1.0 : 0.0);
        boolean isoformRestructured =
                (totalExonDiff > 0 && totalExonDiff != Math.abs(exonDiff))
                        || (totalCdsChange >= cdsChangePct && Math.abs(totalCdsChange - cdsChange) > 0.001)
                        || (!geneStructureSignature(before).equals(geneStructureSignature(after))
                        && !(primaryExonCoordsChanged || primaryCdsCoordsChanged || primaryUtrCoordsChanged));
        boolean isoformChanged = mrnaCountChanged || isoformRestructured;
        boolean cdsChanged = cdsChange >= cdsChangePct || primaryCdsCoordsChanged;

        if (exactBoundaries
                && geneStructureSignature(before).equals(geneStructureSignature(after))
                && !anyUtrTypeChanged
                && !utrBoundaryRefined
                && !isoformChanged) {
            return "exact";
        }

        List<String> changes = new ArrayList<>();
        if (codingExonDiff > 0) {
            changes.add("exon_gain");
        } else if (codingExonDiff < 0) {
            changes.add("exon_loss");
        }
        if (utrOnlyExonDiff > 0) {
            changes.add("utr_exon_added");
        } else if (utrOnlyExonDiff < 0) {
            changes.add("utr_exon_removed");
        }
        if (codingExonDiff == 0 && utrOnlyExonDiff == 0 && exonDiff == 0 && primaryExonCoordsChanged) {
            changes.add("exon_boundary_refined");
        }
        if (cdsChanged) {
            if (aCds > bCds) {
                changes.add("cds_extended");
            } else if (aCds < bCds) {
                changes.add("cds_truncated");
            } else {
                changes.add("cds_boundary_refined");
            }
        }
        boolean hasExonUtrChange = changes.contains("utr_exon_added") || changes.contains("utr_exon_removed");
        if (!hasExonUtrChange) {
            if (anyUtrGained) {
                changes.add("utr_added");
            }
            if (anyUtrLost) {
                changes.add("utr_lost");
            }
            if (utrBoundaryRefined && !anyUtrTypeChanged) {
                changes.add("utr_refined");
            }
        }
        if (isoformChanged) {
            List<String> isoTags = new ArrayList<>();
            if (mrnaCountChanged) {
                isoTags.add("mrna_" + before.mrnaCount() + "x" + after.mrnaCount());
            }
            if (isoformRestructured) {
                isoTags.add("restructured");
            }
            changes.add("isoform_" + String.join("_", isoTags));
        }
        if (changes.isEmpty()) {
            return !exactBoundaries ? "boundary_refined" : "exact";
        }
        return String.join("_", changes);
    }

    public static ComparisonResult compareAnnotations(Path beforeGff, Path afterGff,
                                                      double minReciprocal, int boundaryTol,
                                                      double cdsChangePct, double utrChangePct,
                                                      String geneScope, String overlapMode) throws IOException {
        if (!"reciprocal".equals(overlapMode)
                && !"containment".equals(overlapMode)
                && !"hybrid".equals(overlapMode)) {
            throw new IllegalArgumentException("Unsupported overlap_mode: " + overlapMode);
        }
        System.out.println("Parsing before: " + beforeGff);
        LinkedHashMap<String, GeneModel> beforeGenes = parseGff3ToModels(beforeGff, geneScope);
        System.out.println("  Found " + beforeGenes.size() + " gene models");
        System.out.println("Parsing after: " + afterGff);
        LinkedHashMap<String, GeneModel> afterGenes = parseGff3ToModels(afterGff, geneScope);
        System.out.println("  Found " + afterGenes.size() + " gene models");

        Map<String, List<GeneModel>> beforeGroups = groupByLocus(beforeGenes);
        Map<String, List<GeneModel>> afterGroups = groupByLocus(afterGenes);
        Set<String> beforeSeqs = seqsFromGroups(beforeGroups);
        Set<String> afterSeqs = seqsFromGroups(afterGroups);
        Set<String> commonSeqs = new HashSet<>(beforeSeqs);
        commonSeqs.retainAll(afterSeqs);
        if (commonSeqs.isEmpty()) {
            System.out.println("  FATAL: No common seqids between before and after!");
            return null;
        }
        System.out.println("  Common seqids: " + commonSeqs.size());

        List<Pair> reciprocalPairs = new ArrayList<>();
        List<Pair> containmentPairs = new ArrayList<>();
        Map<String, Integer> diagnostics = new HashMap<>();
        Map<String, Set<String>> weakRejected = new HashMap<>();
        weakRejected.put("before", new HashSet<>());
        weakRejected.put("after", new HashSet<>());
        for (Map.Entry<String, List<GeneModel>> entry : beforeGroups.entrySet()) {
            String key = entry.getKey();
            String seqid = key.substring(0, key.indexOf('\u0000'));
            if (!afterSeqs.contains(seqid)) {
                continue;
            }
            List<GeneModel> afterList = afterGroups.getOrDefault(key, Collections.emptyList());
            if (entry.getValue().isEmpty() || afterList.isEmpty()) {
                continue;
            }
            List<Pair> pairs = findOverlappingPairs(
                    entry.getValue(), afterList, minReciprocal,
                    "hybrid".equals(overlapMode) ? "reciprocal" : overlapMode,
                    diagnostics, weakRejected);
            if ("hybrid".equals(overlapMode)) {
                reciprocalPairs.addAll(pairs);
                containmentPairs.addAll(findOverlappingPairs(
                        entry.getValue(), afterList, minReciprocal, "containment", null, null));
            } else {
                reciprocalPairs.addAll(pairs);
            }
        }
        List<Pair> graphPairs;
        Set<String> reciprocalPairKeys = new HashSet<>();
        int prunedContainmentBridgeEdges = 0;
        if ("hybrid".equals(overlapMode)) {
            PruneResult pruneResult = pruneContainmentBridgeEdges(containmentPairs, reciprocalPairs);
            graphPairs = pruneResult.pairs;
            prunedContainmentBridgeEdges = pruneResult.prunedCount;
            for (Pair pair : graphPairs) {
                reciprocalPairKeys.add(pair.before.geneId + "\u0000" + pair.after.geneId);
            }
            System.out.println("  Found " + reciprocalPairs.size() + " strict reciprocal candidate pairs");
            System.out.println("  Found " + containmentPairs.size() + " containment candidate pairs");
            System.out.println("  Pruned weak bridge edges between strong 1:1 anchors: "
                    + prunedContainmentBridgeEdges);
        } else {
            graphPairs = reciprocalPairs;
            System.out.println("  Found " + graphPairs.size() + " candidate matching pairs");
        }

        ResolvedMatches graphResolved = resolveMatches(graphPairs, beforeGenes, afterGenes);
        ResolvedMatches resolved = graphResolved;
        List<GenePair> weakOneToOne = new ArrayList<>();
        if ("hybrid".equals(overlapMode)) {
            resolved = new ResolvedMatches();
            for (GenePair pair : graphResolved.syntenic) {
                String key = pair.before.geneId + "\u0000" + pair.after.geneId;
                if (reciprocalPairKeys.contains(key)) {
                    resolved.syntenic.add(pair);
                } else {
                    weakOneToOne.add(pair);
                }
            }
            resolved.split.addAll(graphResolved.split);
            resolved.merge.addAll(graphResolved.merge);
            resolved.complex.addAll(graphResolved.complex);
            resolved.novel.addAll(graphResolved.novel);
            resolved.deleted.addAll(graphResolved.deleted);
        }
        LinkedHashMap<String, Integer> noOverlapCounts = countNoOverlapLoci(beforeGenes, afterGenes);
        List<GeneModel> strictNovel = new ArrayList<>();
        List<GeneModel> strictDeleted = new ArrayList<>();
        List<GeneModel> unresolvedNovel = new ArrayList<>();
        List<GeneModel> unresolvedDeleted = new ArrayList<>();
        if ("hybrid".equals(overlapMode)) {
            strictNovel.addAll(resolved.novel);
            strictDeleted.addAll(resolved.deleted);
            for (GenePair pair : weakOneToOne) {
                unresolvedDeleted.add(pair.before);
                unresolvedNovel.add(pair.after);
            }
        } else {
            Set<String> weakBefore = weakRejected.getOrDefault("before", Collections.emptySet());
            Set<String> weakAfter = weakRejected.getOrDefault("after", Collections.emptySet());
            for (GeneModel after : resolved.novel) {
                if (weakAfter.contains(after.geneId)) {
                    unresolvedNovel.add(after);
                } else {
                    strictNovel.add(after);
                }
            }
            for (GeneModel before : resolved.deleted) {
                if (weakBefore.contains(before.geneId)) {
                    unresolvedDeleted.add(before);
                } else {
                    strictDeleted.add(before);
                }
            }
        }
        List<LinkedHashMap<String, Object>> changeLog = new ArrayList<>();
        LinkedHashMap<String, Integer> syntenicBySubtype = new LinkedHashMap<>();
        LinkedHashMap<String, Integer> attributeCounts = new LinkedHashMap<>();
        for (String key : SYNTENIC_ATTRIBUTE_KEYS) {
            attributeCounts.put(key, 0);
        }
        LinkedHashMap<String, Integer> representativeCounts =
                summarizeRepresentativeTranscriptChanges(resolved.syntenic);

        for (GenePair pair : resolved.syntenic) {
            String subtype = classifySyntenicChange(pair.before, pair.after, boundaryTol, cdsChangePct, utrChangePct);
            syntenicBySubtype.merge(subtype, 1, Integer::sum);
            LinkedHashMap<String, Boolean> attrs = computeSyntenicAttributes(
                    pair.before, pair.after, boundaryTol, cdsChangePct, utrChangePct);
            for (Map.Entry<String, Boolean> attr : attrs.entrySet()) {
                if (Boolean.TRUE.equals(attr.getValue())) {
                    attributeCounts.merge(attr.getKey(), 1, Integer::sum);
                }
            }
            changeLog.add(makeChangeLogRow(pair.before, pair.after, "syntenic", subtype));
        }
        for (SplitEvent event : resolved.split) {
            for (GeneModel after : event.afters) {
                changeLog.add(makeChangeLogRow(event.before, after, "split", "split_into_" + event.afters.size()));
            }
        }
        for (MergeEvent event : resolved.merge) {
            for (GeneModel before : event.befores) {
                changeLog.add(makeChangeLogRow(before, event.after, "merge", "merge_from_" + event.befores.size()));
            }
        }
        for (ComplexEvent event : resolved.complex) {
            for (GenePair pair : event.edgePairs) {
                changeLog.add(makeChangeLogRow(pair.before, pair.after, "complex",
                        "complex_" + event.befores.size() + "x" + event.afters.size()));
            }
        }
        for (GeneModel after : strictNovel) {
            changeLog.add(makeChangeLogRow(null, after, "novel", "new_gene"));
        }
        for (GeneModel before : strictDeleted) {
            changeLog.add(makeChangeLogRow(before, null, "deleted", "lost_gene"));
        }
        for (GeneModel after : unresolvedNovel) {
            changeLog.add(makeChangeLogRow(null, after, "unresolved_overlap_after", "weak_overlap_new_gene"));
        }
        for (GeneModel before : unresolvedDeleted) {
            changeLog.add(makeChangeLogRow(before, null, "unresolved_overlap_before", "weak_overlap_lost_gene"));
        }

        LinkedHashMap<String, Object> summary = new LinkedHashMap<>();
        summary.put("gene_scope", geneScope);
        summary.put("overlap_mode", overlapMode);
        summary.put("overlap_threshold", minReciprocal);
        summary.put("boundary_tolerance_bp", boundaryTol);
        summary.put("cds_change_threshold", cdsChangePct);
        summary.put("utr_change_threshold", utrChangePct);
        summary.put("candidate_pairs", "hybrid".equals(overlapMode) ? reciprocalPairs.size() : graphPairs.size());
        summary.put("same_strand_overlaps", diagnostics.getOrDefault("same_strand_overlaps", 0));
        summary.put("containment_pairs_filtered_by_reciprocal",
                diagnostics.getOrDefault("containment_pairs_filtered_by_reciprocal", 0));
        summary.put("containment_candidate_pairs", "hybrid".equals(overlapMode) ? containmentPairs.size() : "");
        summary.put("containment_bridge_edges_pruned", prunedContainmentBridgeEdges);
        summary.put("total_before_genes", beforeGenes.size());
        summary.put("total_after_genes", afterGenes.size());
        summary.put("syntenic_total", resolved.syntenic.size());
        int exact = attributeCounts.getOrDefault("exact", 0);
        summary.put("changed_before_genes", beforeGenes.size() - exact);
        summary.put("changed_before_pct",
                beforeGenes.isEmpty() ? 0.0 : (beforeGenes.size() - exact) / (double) beforeGenes.size() * 100.0);
        summary.put("changed_after_genes", afterGenes.size() - exact);
        summary.put("changed_after_pct",
                afterGenes.isEmpty() ? 0.0 : (afterGenes.size() - exact) / (double) afterGenes.size() * 100.0);
        for (Map.Entry<String, Integer> entry : noOverlapCounts.entrySet()) {
            summary.put(entry.getKey(), entry.getValue());
        }
        for (Map.Entry<String, Integer> entry : representativeCounts.entrySet()) {
            summary.put(entry.getKey(), entry.getValue());
        }
        int repPairs = representativeCounts.getOrDefault("rep_transcript_pairs", 0);
        int repStructuralChanged = representativeCounts.getOrDefault("rep_structural_changed", 0);
        summary.put("rep_structural_changed_pct",
                repPairs == 0 ? 0.0 : repStructuralChanged / (double) repPairs * 100.0);
        summary.put("split_events", resolved.split.size());
        summary.put("merge_events", resolved.merge.size());
        summary.put("complex_events", resolved.complex.size());
        summary.put("unresolved_overlap_after_genes", unresolvedNovel.size());
        summary.put("unresolved_overlap_before_genes", unresolvedDeleted.size());
        summary.put("novel_genes", strictNovel.size());
        summary.put("deleted_genes", strictDeleted.size());
        for (String key : SYNTENIC_ATTRIBUTE_KEYS) {
            summary.put("one_to_one_" + key, attributeCounts.getOrDefault(key, 0));
        }
        for (Map.Entry<String, Integer> entry : syntenicBySubtype.entrySet()) {
            summary.put("syntenic_" + entry.getKey(), entry.getValue());
        }
        summary.put("before_in_splits", resolved.split.size());
        summary.put("before_in_merges", resolved.merge.stream().mapToInt(e -> e.befores.size()).sum());
        summary.put("before_in_complex", resolved.complex.stream().mapToInt(e -> e.befores.size()).sum());
        summary.put("after_in_splits", resolved.split.stream().mapToInt(e -> e.afters.size()).sum());
        summary.put("after_in_merges", resolved.merge.size());
        summary.put("after_in_complex", resolved.complex.stream().mapToInt(e -> e.afters.size()).sum());
        return new ComparisonResult(summary, changeLog);
    }

    private static Map<String, List<GeneModel>> groupByLocus(LinkedHashMap<String, GeneModel> genes) {
        Map<String, List<GeneModel>> groups = new LinkedHashMap<>();
        for (GeneModel gene : genes.values()) {
            groups.computeIfAbsent(groupKey(gene.seqid, gene.strand), k -> new ArrayList<>()).add(gene);
        }
        for (List<GeneModel> group : groups.values()) {
            group.sort(Comparator.comparingInt((GeneModel g) -> g.start).thenComparing(g -> g.geneId));
        }
        return groups;
    }

    private static String groupKey(String seqid, String strand) {
        return seqid + '\u0000' + strand;
    }

    private static Set<String> seqsFromGroups(Map<String, List<GeneModel>> groups) {
        Set<String> seqs = new HashSet<>();
        for (String key : groups.keySet()) {
            seqs.add(key.substring(0, key.indexOf('\u0000')));
        }
        return seqs;
    }

    public static LinkedHashMap<String, Object> makeChangeLogRow(GeneModel before, GeneModel after,
                                                                  String matchType, String changeSubtype) {
        GeneModel anchor = before != null ? before : after;
        LinkedHashMap<String, Object> row = new LinkedHashMap<>();
        row.put("before_gene", before != null ? before.geneId : "");
        row.put("after_gene", after != null ? after.geneId : "");
        row.put("seqid", anchor != null ? anchor.seqid : "");
        row.put("before_start", before != null ? before.start : 0);
        row.put("before_end", before != null ? before.end : 0);
        row.put("after_start", after != null ? after.start : 0);
        row.put("after_end", after != null ? after.end : 0);
        row.put("before_gene_start", before != null ? rawGeneBoundaryStart(before) : 0);
        row.put("before_gene_end", before != null ? rawGeneBoundaryEnd(before) : 0);
        row.put("after_gene_start", after != null ? rawGeneBoundaryStart(after) : 0);
        row.put("after_gene_end", after != null ? rawGeneBoundaryEnd(after) : 0);
        row.put("strand", anchor != null ? anchor.strand : "");
        row.put("match_type", matchType);
        row.put("change_subtype", changeSubtype);
        row.put("before_length", before != null ? before.length() : 0);
        row.put("after_length", after != null ? after.length() : 0);
        row.put("before_gene_length", before != null ? rawGeneBoundaryLength(before) : 0);
        row.put("after_gene_length", after != null ? rawGeneBoundaryLength(after) : 0);
        row.put("before_exons", before != null ? before.totalExonCount() : 0);
        row.put("after_exons", after != null ? after.totalExonCount() : 0);
        row.put("before_cds", before != null ? before.totalCdsLength() : 0);
        row.put("after_cds", after != null ? after.totalCdsLength() : 0);
        row.put("before_mrnas", before != null ? before.mrnaCount() : 0);
        row.put("after_mrnas", after != null ? after.mrnaCount() : 0);
        return row;
    }

    public static void writeCsv(Path path, List<String> columns, List<? extends Map<String, ?>> rows) throws IOException {
        Files.createDirectories(path.getParent());
        try (BufferedWriter writer = Files.newBufferedWriter(path, StandardCharsets.UTF_8)) {
            writer.write(String.join(",", columns));
            writer.newLine();
            for (Map<String, ?> row : rows) {
                for (int i = 0; i < columns.size(); i++) {
                    if (i > 0) {
                        writer.write(',');
                    }
                    Object value = row.get(columns.get(i));
                    writer.write(csvEscape(value == null ? "" : String.valueOf(value)));
                }
                writer.newLine();
            }
        }
    }

    public static void writeSummaryCsv(Path path, LinkedHashMap<String, Object> summary) throws IOException {
        writeCsv(path, new ArrayList<>(summary.keySet()), Collections.singletonList(summary));
    }

    public static void writeChangeLogCsv(Path path, List<LinkedHashMap<String, Object>> changeLog) throws IOException {
        writeCsv(path, CHANGE_LOG_COLUMNS, changeLog);
    }

    private static String csvEscape(String value) {
        if (value.contains(",") || value.contains("\"") || value.contains("\n") || value.contains("\r")) {
            return "\"" + value.replace("\"", "\"\"") + "\"";
        }
        return value;
    }

    public static List<LinkedHashMap<String, String>> readCsv(Path path) throws IOException {
        List<LinkedHashMap<String, String>> rows = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            String headerLine = reader.readLine();
            if (headerLine == null) {
                return rows;
            }
            List<String> headers = parseCsvLine(headerLine);
            String line;
            while ((line = reader.readLine()) != null) {
                List<String> values = parseCsvLine(line);
                LinkedHashMap<String, String> row = new LinkedHashMap<>();
                for (int i = 0; i < headers.size(); i++) {
                    row.put(headers.get(i), i < values.size() ? values.get(i) : "");
                }
                rows.add(row);
            }
        }
        return rows;
    }

    public static List<String> readCsvHeader(Path path) throws IOException {
        try (BufferedReader reader = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            String header = reader.readLine();
            return header == null ? Collections.emptyList() : parseCsvLine(header);
        }
    }

    private static List<String> parseCsvLine(String line) {
        List<String> values = new ArrayList<>();
        StringBuilder current = new StringBuilder();
        boolean quoted = false;
        for (int i = 0; i < line.length(); i++) {
            char ch = line.charAt(i);
            if (quoted) {
                if (ch == '"') {
                    if (i + 1 < line.length() && line.charAt(i + 1) == '"') {
                        current.append('"');
                        i++;
                    } else {
                        quoted = false;
                    }
                } else {
                    current.append(ch);
                }
            } else {
                if (ch == ',') {
                    values.add(current.toString());
                    current.setLength(0);
                } else if (ch == '"') {
                    quoted = true;
                } else {
                    current.append(ch);
                }
            }
        }
        values.add(current.toString());
        return values;
    }

    public static List<Species> readSpecies(Path speciesJson) throws IOException {
        String text = Files.readString(speciesJson, StandardCharsets.UTF_8);
        List<Species> species = new ArrayList<>();
        int pos = 0;
        while (true) {
            int idKey = text.indexOf("\"id\"", pos);
            if (idKey < 0) {
                break;
            }
            String id = jsonStringValue(text, idKey);
            int labelKey = text.indexOf("\"label\"", idKey);
            String label = jsonStringValue(text, labelKey);
            int shortKey = text.indexOf("\"short_label\"", labelKey);
            String shortLabel = jsonStringValue(text, shortKey);
            species.add(new Species(id, label, shortLabel));
            pos = shortKey + 1;
        }
        return species;
    }

    private static String jsonStringValue(String text, int keyPos) {
        int colon = text.indexOf(':', keyPos);
        int firstQuote = text.indexOf('"', colon + 1);
        int secondQuote = text.indexOf('"', firstQuote + 1);
        return text.substring(firstQuote + 1, secondQuote);
    }

    public static final class Species {
        public final String id;
        public final String label;
        public final String shortLabel;

        Species(String id, String label, String shortLabel) {
            this.id = id;
            this.label = label;
            this.shortLabel = shortLabel;
        }
    }

    public static Path findAnnotation(String speciesId, String state, Path analysisDir) throws IOException {
        List<Path> candidates = new ArrayList<>();
        try (var stream = Files.list(analysisDir)) {
            stream.filter(Files::isRegularFile)
                    .filter(path -> path.getFileName().toString().startsWith(speciesId + "." + state + "."))
                    .filter(Core::isAnnotationPath)
                    .forEach(candidates::add);
        }
        Collections.sort(candidates);
        return candidates.isEmpty() ? null : candidates.get(0);
    }

    private static boolean isAnnotationPath(Path path) {
        String name = path.getFileName().toString();
        return name.endsWith(".gff") || name.endsWith(".gff3")
                || name.endsWith(".gff.gz") || name.endsWith(".gff3.gz");
    }

    public static int intValue(Map<String, String> row, String column) {
        String value = row.get(column);
        if (value == null || value.isEmpty()) {
            return 0;
        }
        return (int) Double.parseDouble(value);
    }

    public static int intObject(Object value) {
        if (value instanceof Number) {
            return ((Number) value).intValue();
        }
        if (value == null || String.valueOf(value).isEmpty()) {
            return 0;
        }
        return (int) Double.parseDouble(String.valueOf(value));
    }

    public static LinkedHashMap<String, Integer> newCategoryMap() {
        LinkedHashMap<String, Integer> cats = new LinkedHashMap<>();
        cats.put("exact", 0);
        cats.put("boundary_refined", 0);
        cats.put("utr_added", 0);
        cats.put("utr_lost", 0);
        cats.put("utr_exon_added", 0);
        cats.put("utr_exon_removed", 0);
        cats.put("utr_refined", 0);
        cats.put("exon_gain_cds_change", 0);
        cats.put("exon_loss_cds_change", 0);
        cats.put("exon_boundary_refined", 0);
        cats.put("cds_change_only", 0);
        cats.put("cds_boundary_refined", 0);
        cats.put("isoform", 0);
        return cats;
    }

    public static String classifyExclusiveSubtype(String subtype) {
        if ("total".equals(subtype)) {
            return null;
        }
        List<Map.Entry<String, Predicate<String>>> rules = Arrays.asList(
                Map.entry("exact", s -> s.equals("exact")),
                Map.entry("boundary_refined", s -> s.equals("boundary_refined")),
                Map.entry("isoform", s -> s.contains("isoform_")),
                Map.entry("utr_exon_added", s -> s.equals("utr_exon_added") || s.startsWith("utr_exon_added_")),
                Map.entry("utr_exon_removed", s -> s.equals("utr_exon_removed") || s.startsWith("utr_exon_removed_")),
                Map.entry("utr_added", s -> s.equals("utr_added") || s.startsWith("utr_added_")),
                Map.entry("utr_lost", s -> s.equals("utr_lost") || s.startsWith("utr_lost_")),
                Map.entry("utr_refined", s -> s.equals("utr_refined") || s.startsWith("utr_refined_")),
                Map.entry("exon_gain_cds_change", s -> s.startsWith("exon_gain_")),
                Map.entry("exon_loss_cds_change", s -> s.startsWith("exon_loss_")),
                Map.entry("cds_boundary_refined", s -> s.contains("cds_boundary_refined")),
                Map.entry("exon_boundary_refined", s -> s.contains("exon_boundary_refined")),
                Map.entry("cds_change_only", s -> s.contains("cds_extended") || s.contains("cds_truncated"))
        );
        for (Map.Entry<String, Predicate<String>> rule : rules) {
            if (rule.getValue().test(subtype)) {
                return rule.getKey();
            }
        }
        throw new IllegalArgumentException("Unclassified syntenic subtype: " + subtype);
    }

    public static Set<String> unionColumns(Path a, Path b) throws IOException {
        Set<String> columns = new TreeSet<>();
        columns.addAll(readCsvHeader(a));
        columns.addAll(readCsvHeader(b));
        return columns;
    }
}
