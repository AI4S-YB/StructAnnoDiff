package locuscompare;

import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;

public final class LocusCompareMain {
    private LocusCompareMain() {
    }

    public static void main(String[] args) throws Exception {
        Map<String, String> opts = parseArgs(args);
        Path before = requiredPath(opts, "--before");
        Path after = requiredPath(opts, "--after");
        Path output = Path.of(opts.getOrDefault("--output", "results/locus"));
        String name = opts.getOrDefault("--name", inferName(before));
        double reciprocal = Double.parseDouble(opts.getOrDefault("--reciprocal-overlap", "0.5"));
        int boundaryTol = Integer.parseInt(opts.getOrDefault("--boundary-tol", "10"));
        double cdsPct = Double.parseDouble(opts.getOrDefault("--cds-change-pct", "0.1"));
        double utrPct = Double.parseDouble(opts.getOrDefault("--utr-change-pct", "0.1"));
        String geneScope = opts.getOrDefault("--gene-scope", "mrna");
        String overlapMode = opts.getOrDefault("--overlap-mode", "reciprocal");

        System.out.println("============================================================");
        System.out.println("Species: " + name);
        System.out.println("  Before: " + before);
        System.out.println("  After:  " + after);
        System.out.println("  Gene scope: " + geneScope);
        System.out.println("  Overlap mode: " + overlapMode);
        System.out.println("  Reciprocal overlap threshold: " + reciprocal);
        System.out.println("============================================================");

        Core.ComparisonResult result = Core.compareAnnotations(
                before, after, reciprocal, boundaryTol, cdsPct, utrPct, geneScope, overlapMode);
        if (result == null) {
            System.exit(1);
        }

        Core.writeSummaryCsv(output.resolve(name + "_change_summary.csv"), result.summary);
        Core.writeChangeLogCsv(output.resolve(name + "_change_log.csv"), result.changeLog);

        System.out.println("Summary saved: " + output.resolve(name + "_change_summary.csv"));
        System.out.println("Change log saved: " + output.resolve(name + "_change_log.csv")
                + " (" + result.changeLog.size() + " entries)");
    }

    static Map<String, String> parseArgs(String[] args) {
        LinkedHashMap<String, String> opts = new LinkedHashMap<>();
        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            if (!arg.startsWith("--")) {
                throw new IllegalArgumentException("Unexpected positional argument: " + arg);
            }
            if (i + 1 >= args.length || args[i + 1].startsWith("--")) {
                throw new IllegalArgumentException("Missing value for " + arg);
            }
            opts.put(arg, args[++i]);
        }
        return opts;
    }

    private static Path requiredPath(Map<String, String> opts, String key) {
        String value = opts.get(key);
        if (value == null) {
            throw new IllegalArgumentException("Missing required option " + key);
        }
        return Path.of(value);
    }

    private static String inferName(Path before) {
        String name = before.getFileName().toString();
        int dot = name.indexOf(".before");
        if (dot >= 0) {
            return name.substring(0, dot);
        }
        int lastDot = name.lastIndexOf('.');
        return lastDot > 0 ? name.substring(0, lastDot) : name;
    }
}
