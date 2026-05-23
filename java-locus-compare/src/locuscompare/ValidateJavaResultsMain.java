package locuscompare;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

public final class ValidateJavaResultsMain {
    private ValidateJavaResultsMain() {
    }

    public static void main(String[] args) throws Exception {
        Map<String, String> opts = LocusCompareMain.parseArgs(args);
        Path pythonResults = Path.of(opts.getOrDefault("--python-results", "../results")).toAbsolutePath().normalize();
        Path javaResults = Path.of(opts.getOrDefault("--java-results", "results")).toAbsolutePath().normalize();

        Path analysisDir = pythonResults.getParent();
        List<Core.Species> species = Core.readSpecies(analysisDir.resolve("species.json"));
        List<String> issues = new ArrayList<>();

        for (Core.Species sp : species) {
            compareCsvSingleRow(
                    pythonResults.resolve("locus").resolve(sp.id + "_change_summary.csv"),
                    javaResults.resolve("locus").resolve(sp.id + "_change_summary.csv"),
                    sp.id + " summary",
                    issues
            );
            compareChangeLogCounts(
                    pythonResults.resolve("locus").resolve(sp.id + "_change_log.csv"),
                    javaResults.resolve("locus").resolve(sp.id + "_change_log.csv"),
                    sp.id + " change_log",
                    issues
            );
        }

        compareCsvTable(
                pythonResults.resolve("locus_comparison_summary.csv"),
                javaResults.resolve("locus_comparison_summary.csv"),
                "locus_comparison_summary.csv",
                issues
        );
        compareCsvTable(
                pythonResults.resolve("locus_comparison_multilabel.csv"),
                javaResults.resolve("locus_comparison_multilabel.csv"),
                "locus_comparison_multilabel.csv",
                issues
        );
        compareCsvTable(
                pythonResults.resolve("locus_diagnostics.csv"),
                javaResults.resolve("locus_diagnostics.csv"),
                "locus_diagnostics.csv",
                issues
        );

        assertNoExcludedOutputs(javaResults, issues);

        if (!issues.isEmpty()) {
            System.err.println("Java/Python parity failed:");
            for (String issue : issues) {
                System.err.println("  - " + issue);
            }
            System.exit(1);
        }
        System.out.println("Java/Python parity checks: OK");
    }

    private static void compareCsvSingleRow(Path expectedPath, Path actualPath, String label,
                                            List<String> issues) throws Exception {
        if (!Files.exists(expectedPath) || !Files.exists(actualPath)) {
            issues.add(label + ": missing expected or actual file");
            return;
        }
        List<LinkedHashMap<String, String>> expectedRows = Core.readCsv(expectedPath);
        List<LinkedHashMap<String, String>> actualRows = Core.readCsv(actualPath);
        if (expectedRows.size() != 1 || actualRows.size() != 1) {
            issues.add(label + ": expected one row, got " + expectedRows.size() + " and " + actualRows.size());
            return;
        }
        compareRows(expectedRows.get(0), actualRows.get(0), label, issues);
    }

    private static void compareCsvTable(Path expectedPath, Path actualPath, String label,
                                        List<String> issues) throws Exception {
        if (!Files.exists(expectedPath) || !Files.exists(actualPath)) {
            issues.add(label + ": missing expected or actual file");
            return;
        }
        List<LinkedHashMap<String, String>> expectedRows = Core.readCsv(expectedPath);
        List<LinkedHashMap<String, String>> actualRows = Core.readCsv(actualPath);
        if (expectedRows.size() != actualRows.size()) {
            issues.add(label + ": row count " + actualRows.size() + " != " + expectedRows.size());
            return;
        }
        for (int i = 0; i < expectedRows.size(); i++) {
            compareRows(expectedRows.get(i), actualRows.get(i), label + " row " + i, issues);
        }
    }

    private static void compareRows(Map<String, String> expected, Map<String, String> actual,
                                    String label, List<String> issues) {
        Set<String> cols = new TreeSet<>();
        cols.addAll(expected.keySet());
        cols.addAll(actual.keySet());
        for (String col : cols) {
            if (!expected.containsKey(col)) {
                issues.add(label + ": unexpected column " + col);
                continue;
            }
            if (!actual.containsKey(col)) {
                issues.add(label + ": missing column " + col);
                continue;
            }
            String e = expected.getOrDefault(col, "");
            String a = actual.getOrDefault(col, "");
            if (isNumeric(e) && isNumeric(a)) {
                double ed = Double.parseDouble(e);
                double ad = Double.parseDouble(a);
                if (Math.abs(ed - ad) > 1e-9) {
                    issues.add(label + ": " + col + " " + a + " != " + e);
                }
            } else if (!e.equals(a)) {
                issues.add(label + ": " + col + " " + a + " != " + e);
            }
        }
    }

    private static void compareChangeLogCounts(Path expectedPath, Path actualPath, String label,
                                               List<String> issues) throws Exception {
        if (!Files.exists(expectedPath) || !Files.exists(actualPath)) {
            issues.add(label + ": missing expected or actual file");
            return;
        }
        Map<String, Integer> expected = changeLogCounts(expectedPath);
        Map<String, Integer> actual = changeLogCounts(actualPath);
        Set<String> keys = new TreeSet<>();
        keys.addAll(expected.keySet());
        keys.addAll(actual.keySet());
        for (String key : keys) {
            int e = expected.getOrDefault(key, 0);
            int a = actual.getOrDefault(key, 0);
            if (e != a) {
                issues.add(label + ": " + key + " count " + a + " != " + e);
            }
        }
    }

    private static Map<String, Integer> changeLogCounts(Path path) throws Exception {
        Map<String, Integer> counts = new HashMap<>();
        for (Map<String, String> row : Core.readCsv(path)) {
            String key = row.getOrDefault("match_type", "") + "|" + row.getOrDefault("change_subtype", "");
            counts.merge(key, 1, Integer::sum);
        }
        return counts;
    }

    private static boolean isNumeric(String value) {
        if (value == null || value.isEmpty()) {
            return false;
        }
        try {
            Double.parseDouble(value);
            return true;
        } catch (NumberFormatException ex) {
            return false;
        }
    }

    private static void assertNoExcludedOutputs(Path javaResults, List<String> issues) throws Exception {
        Path locus = javaResults.resolve("locus");
        if (!Files.exists(locus)) {
            issues.add("Java locus results directory missing");
            return;
        }
        try (var stream = Files.list(locus)) {
            stream.map(path -> path.getFileName().toString()).forEach(name -> {
                if (name.contains("Actinidia_chinensis") || name.contains("Ipomoea_batatas")) {
                    issues.add("Excluded species output exists: " + name);
                }
            });
        }
    }
}
