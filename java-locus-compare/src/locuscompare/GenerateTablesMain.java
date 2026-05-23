package locuscompare;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public final class GenerateTablesMain {
    private static final List<String> SUMMARY_COLUMNS = Arrays.asList(
            "Species", "Before", "After", "Syntenic", "Exact", "Boundary_refined_only",
            "UTR_added", "UTR_lost", "UTR_exon_gained", "UTR_exon_removed", "UTR_refined",
            "Coding_exon_gain", "Coding_exon_loss", "Exon_boundary_refined", "CDS_change_only",
            "CDS_boundary_refined", "Isoform_change", "Split_events", "Merge_events",
            "Complex_events", "Unresolved_overlap_after_genes", "Unresolved_overlap_before_genes",
            "Novel_genes", "Deleted_genes"
    );

    private static final LinkedHashMap<String, String> DIRECT_COUNTS = new LinkedHashMap<>();

    static {
        DIRECT_COUNTS.put("Exact", "one_to_one_exact");
        DIRECT_COUNTS.put("Any_gene_boundary_changed", "one_to_one_gene_boundary_changed");
        DIRECT_COUNTS.put("Any_UTR_added", "one_to_one_utr_added");
        DIRECT_COUNTS.put("Any_UTR_lost", "one_to_one_utr_lost");
        DIRECT_COUNTS.put("Any_UTR_exon_gained", "one_to_one_utr_exon_added");
        DIRECT_COUNTS.put("Any_UTR_exon_removed", "one_to_one_utr_exon_removed");
        DIRECT_COUNTS.put("Any_UTR_refined", "one_to_one_utr_refined");
        DIRECT_COUNTS.put("Any_coding_exon_gain", "one_to_one_coding_exon_gain");
        DIRECT_COUNTS.put("Any_coding_exon_loss", "one_to_one_coding_exon_loss");
        DIRECT_COUNTS.put("Any_exon_boundary_refined", "one_to_one_exon_boundary_refined");
        DIRECT_COUNTS.put("Any_CDS_change", "one_to_one_cds_change");
        DIRECT_COUNTS.put("Any_CDS_boundary_refined", "one_to_one_cds_boundary_refined");
        DIRECT_COUNTS.put("Any_isoform_change", "one_to_one_isoform_change");
    }

    private GenerateTablesMain() {
    }

    public static void main(String[] args) throws Exception {
        Map<String, String> opts = LocusCompareMain.parseArgs(args);
        Path analysisDir = Path.of(opts.getOrDefault("--analysis-dir", "..")).toAbsolutePath().normalize();
        Path javaResults = Path.of(opts.getOrDefault("--java-results", "results"));
        Path locusDir = javaResults.resolve("locus");
        Files.createDirectories(javaResults);

        List<Core.Species> species = Core.readSpecies(analysisDir.resolve("species.json"));
        List<LinkedHashMap<String, Object>> summaryRows = buildMasterTable(species, locusDir);
        Core.writeCsv(javaResults.resolve("locus_comparison_summary.csv"), SUMMARY_COLUMNS, summaryRows);

        List<LinkedHashMap<String, Object>> multilabelRows = buildMultilabelTable(species, locusDir);
        Core.writeCsv(javaResults.resolve("locus_comparison_multilabel.csv"), multilabelColumns(), multilabelRows);

        List<LinkedHashMap<String, Object>> diagnosticsRows = buildDiagnosticsTable(species, locusDir);
        Core.writeCsv(javaResults.resolve("locus_diagnostics.csv"), diagnosticsColumns(), diagnosticsRows);

        System.out.println("Saved: " + javaResults.resolve("locus_comparison_summary.csv"));
        System.out.println("Saved: " + javaResults.resolve("locus_comparison_multilabel.csv"));
        System.out.println("Saved: " + javaResults.resolve("locus_diagnostics.csv"));
    }

    static List<LinkedHashMap<String, Object>> buildMasterTable(List<Core.Species> species, Path locusDir) throws Exception {
        List<LinkedHashMap<String, Object>> rows = new ArrayList<>();
        for (Core.Species sp : species) {
            Path path = locusDir.resolve(sp.id + "_change_summary.csv");
            if (!Files.exists(path)) {
                continue;
            }
            Map<String, String> row = Core.readCsv(path).get(0);
            LinkedHashMap<String, Integer> cats = aggregateSubtypes(row);
            int syntenic = Core.intValue(row, "syntenic_total");
            int subtypeSum = cats.values().stream().mapToInt(Integer::intValue).sum();
            if (subtypeSum != syntenic) {
                throw new IllegalStateException(sp.id + ": exclusive syntenic subtype total "
                        + subtypeSum + " != syntenic_total " + syntenic);
            }
            LinkedHashMap<String, Object> out = new LinkedHashMap<>();
            out.put("Species", sp.label);
            out.put("Before", Core.intValue(row, "total_before_genes"));
            out.put("After", Core.intValue(row, "total_after_genes"));
            out.put("Syntenic", syntenic);
            out.put("Exact", cats.get("exact"));
            out.put("Boundary_refined_only", cats.get("boundary_refined"));
            out.put("UTR_added", cats.get("utr_added"));
            out.put("UTR_lost", cats.get("utr_lost"));
            out.put("UTR_exon_gained", cats.get("utr_exon_added"));
            out.put("UTR_exon_removed", cats.get("utr_exon_removed"));
            out.put("UTR_refined", cats.get("utr_refined"));
            out.put("Coding_exon_gain", cats.get("exon_gain_cds_change"));
            out.put("Coding_exon_loss", cats.get("exon_loss_cds_change"));
            out.put("Exon_boundary_refined", cats.get("exon_boundary_refined"));
            out.put("CDS_change_only", cats.get("cds_change_only"));
            out.put("CDS_boundary_refined", cats.get("cds_boundary_refined"));
            out.put("Isoform_change", cats.get("isoform"));
            out.put("Split_events", Core.intValue(row, "split_events"));
            out.put("Merge_events", Core.intValue(row, "merge_events"));
            out.put("Complex_events", Core.intValue(row, "complex_events"));
            out.put("Unresolved_overlap_after_genes", Core.intValue(row, "unresolved_overlap_after_genes"));
            out.put("Unresolved_overlap_before_genes", Core.intValue(row, "unresolved_overlap_before_genes"));
            out.put("Novel_genes", Core.intValue(row, "novel_genes"));
            out.put("Deleted_genes", Core.intValue(row, "deleted_genes"));
            rows.add(out);
        }
        return rows;
    }

    static LinkedHashMap<String, Integer> aggregateSubtypes(Map<String, String> row) {
        LinkedHashMap<String, Integer> cats = Core.newCategoryMap();
        for (Map.Entry<String, String> entry : row.entrySet()) {
            if (!entry.getKey().startsWith("syntenic_")) {
                continue;
            }
            int value = parseInt(entry.getValue());
            if (value == 0) {
                continue;
            }
            String subtype = entry.getKey().substring("syntenic_".length());
            String category = Core.classifyExclusiveSubtype(subtype);
            if (category != null) {
                cats.put(category, cats.get(category) + value);
            }
        }
        return cats;
    }

    static List<LinkedHashMap<String, Object>> buildMultilabelTable(List<Core.Species> species, Path locusDir) throws Exception {
        List<LinkedHashMap<String, Object>> rows = new ArrayList<>();
        for (Core.Species sp : species) {
            Path path = locusDir.resolve(sp.id + "_change_summary.csv");
            if (!Files.exists(path)) {
                continue;
            }
            Map<String, String> row = Core.readCsv(path).get(0);
            LinkedHashMap<String, Integer> cats = aggregateMultilabelSubtypes(row);
            for (String sourceCol : DIRECT_COUNTS.values()) {
                if (!row.containsKey(sourceCol)) {
                    throw new IllegalStateException(sp.id + ": missing direct multilabel count column " + sourceCol);
                }
            }
            for (Map.Entry<String, String> entry : DIRECT_COUNTS.entrySet()) {
                cats.put(entry.getKey(), Core.intValue(row, entry.getValue()));
            }
            LinkedHashMap<String, Object> out = new LinkedHashMap<>();
            out.put("Species", sp.label);
            out.put("Syntenic", Core.intValue(row, "syntenic_total"));
            for (String col : multilabelColumns()) {
                if (!"Species".equals(col) && !"Syntenic".equals(col)) {
                    out.put(col, cats.getOrDefault(col, 0));
                }
            }
            rows.add(out);
        }
        return rows;
    }

    static LinkedHashMap<String, Integer> aggregateMultilabelSubtypes(Map<String, String> row) {
        LinkedHashMap<String, Integer> cats = new LinkedHashMap<>();
        for (String col : multilabelColumns()) {
            if (!"Species".equals(col) && !"Syntenic".equals(col)) {
                cats.put(col, 0);
            }
        }
        for (Map.Entry<String, String> entry : row.entrySet()) {
            if (!entry.getKey().startsWith("syntenic_")) {
                continue;
            }
            int value = parseInt(entry.getValue());
            if (value == 0) {
                continue;
            }
            String subtype = entry.getKey().substring("syntenic_".length());
            if ("total".equals(subtype)) {
                continue;
            }
            Core.classifyExclusiveSubtype(subtype);
            if ("exact".equals(subtype)) {
                cats.merge("Exact", value, Integer::sum);
            }
            if ("boundary_refined".equals(subtype)) {
                cats.merge("Boundary_refined", value, Integer::sum);
                cats.merge("Any_gene_boundary_changed", value, Integer::sum);
            }
            if ("utr_added".equals(subtype) || subtype.startsWith("utr_added_") || subtype.contains("_utr_added")) {
                cats.merge("Any_UTR_added", value, Integer::sum);
            }
            if ("utr_lost".equals(subtype) || subtype.startsWith("utr_lost_") || subtype.contains("_utr_lost")) {
                cats.merge("Any_UTR_lost", value, Integer::sum);
            }
            if ("utr_exon_added".equals(subtype) || subtype.startsWith("utr_exon_added_")) {
                cats.merge("Any_UTR_exon_gained", value, Integer::sum);
            }
            if ("utr_exon_removed".equals(subtype) || subtype.startsWith("utr_exon_removed_")) {
                cats.merge("Any_UTR_exon_removed", value, Integer::sum);
            }
            if ("utr_refined".equals(subtype) || subtype.startsWith("utr_refined_") || subtype.contains("_utr_refined")) {
                cats.merge("Any_UTR_refined", value, Integer::sum);
            }
            if (subtype.startsWith("exon_gain_")) {
                cats.merge("Any_coding_exon_gain", value, Integer::sum);
            }
            if (subtype.startsWith("exon_loss_")) {
                cats.merge("Any_coding_exon_loss", value, Integer::sum);
            }
            if (subtype.contains("exon_boundary_refined")) {
                cats.merge("Any_exon_boundary_refined", value, Integer::sum);
            }
            if (subtype.contains("cds_extended") || subtype.contains("cds_truncated")
                    || subtype.contains("cds_boundary_refined")) {
                cats.merge("Any_CDS_change", value, Integer::sum);
            }
            if (subtype.contains("cds_boundary_refined")) {
                cats.merge("Any_CDS_boundary_refined", value, Integer::sum);
            }
            if (subtype.contains("isoform_")) {
                cats.merge("Any_isoform_change", value, Integer::sum);
            }
        }
        return cats;
    }

    static List<LinkedHashMap<String, Object>> buildDiagnosticsTable(List<Core.Species> species, Path locusDir) throws Exception {
        List<LinkedHashMap<String, Object>> rows = new ArrayList<>();
        for (Core.Species sp : species) {
            Path path = locusDir.resolve(sp.id + "_change_summary.csv");
            if (!Files.exists(path)) {
                continue;
            }
            Map<String, String> row = Core.readCsv(path).get(0);
            LinkedHashMap<String, Object> out = new LinkedHashMap<>();
            out.put("Species", sp.label);
            out.put("Overlap_mode", row.getOrDefault("overlap_mode", ""));
            out.put("Overlap_threshold", row.getOrDefault("overlap_threshold", ""));
            out.put("Candidate_pairs", row.getOrDefault("candidate_pairs", ""));
            out.put("Same_strand_overlaps", row.getOrDefault("same_strand_overlaps", ""));
            out.put("Containment_pairs_filtered_by_reciprocal",
                    row.getOrDefault("containment_pairs_filtered_by_reciprocal", ""));
            rows.add(out);
        }
        return rows;
    }

    private static List<String> multilabelColumns() {
        List<String> cols = new ArrayList<>();
        cols.add("Species");
        cols.add("Syntenic");
        cols.add("Exact");
        cols.add("Boundary_refined");
        cols.add("Any_gene_boundary_changed");
        cols.add("Any_UTR_added");
        cols.add("Any_UTR_lost");
        cols.add("Any_UTR_exon_gained");
        cols.add("Any_UTR_exon_removed");
        cols.add("Any_UTR_refined");
        cols.add("Any_coding_exon_gain");
        cols.add("Any_coding_exon_loss");
        cols.add("Any_exon_boundary_refined");
        cols.add("Any_CDS_change");
        cols.add("Any_CDS_boundary_refined");
        cols.add("Any_isoform_change");
        return cols;
    }

    private static List<String> diagnosticsColumns() {
        return Arrays.asList(
                "Species", "Overlap_mode", "Overlap_threshold", "Candidate_pairs",
                "Same_strand_overlaps", "Containment_pairs_filtered_by_reciprocal"
        );
    }

    private static int parseInt(String value) {
        if (value == null || value.isEmpty()) {
            return 0;
        }
        return (int) Double.parseDouble(value);
    }
}
