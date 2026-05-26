package locuscompare;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public final class RunLocusComparisonsMain {
    private RunLocusComparisonsMain() {
    }

    public static void main(String[] args) throws Exception {
        Map<String, String> opts = LocusCompareMain.parseArgs(args);
        Path analysisDir = Path.of(opts.getOrDefault("--analysis-dir", "..")).toAbsolutePath().normalize();
        Path output = Path.of(opts.getOrDefault("--output", "results/locus"));
        double reciprocal = Double.parseDouble(opts.getOrDefault("--reciprocal-overlap", "0.5"));
        int boundaryTol = Integer.parseInt(opts.getOrDefault("--boundary-tol", "10"));
        double cdsPct = Double.parseDouble(opts.getOrDefault("--cds-change-pct", "0.1"));
        double utrPct = Double.parseDouble(opts.getOrDefault("--utr-change-pct", "0.1"));
        String geneScope = opts.getOrDefault("--gene-scope", "mrna");
        String overlapMode = opts.getOrDefault("--overlap-mode", "hybrid");

        List<Core.Species> species = Core.readSpecies(analysisDir.resolve("species.json"));
        List<String> requested = parseSpeciesArg(opts.get("--species"));
        if (requested.isEmpty()) {
            for (Core.Species sp : species) {
                requested.add(sp.id);
            }
        }
        Set<String> known = new HashSet<>();
        for (Core.Species sp : species) {
            known.add(sp.id);
        }
        for (String id : requested) {
            if (!known.contains(id)) {
                throw new IllegalArgumentException("Unknown species ID: " + id);
            }
        }

        Files.createDirectories(output);
        for (String id : requested) {
            Path before = Core.findAnnotation(id, "before", analysisDir);
            Path after = Core.findAnnotation(id, "after", analysisDir);
            if (before == null || after == null) {
                throw new IllegalStateException("Missing before/after annotation for " + id);
            }
            System.out.println("=== " + id + " ===");
            Core.ComparisonResult result = Core.compareAnnotations(
                    before, after, reciprocal, boundaryTol, cdsPct, utrPct, geneScope, overlapMode);
            if (result == null) {
                throw new IllegalStateException("Comparison failed for " + id);
            }
            Core.writeSummaryCsv(output.resolve(id + "_change_summary.csv"), result.summary);
            Core.writeChangeLogCsv(output.resolve(id + "_change_log.csv"), result.changeLog);
        }
    }

    private static List<String> parseSpeciesArg(String value) {
        List<String> out = new ArrayList<>();
        if (value == null || value.isBlank()) {
            return out;
        }
        for (String part : value.split(",")) {
            String trimmed = part.trim();
            if (!trimmed.isEmpty()) {
                out.add(trimmed);
            }
        }
        return out;
    }
}
