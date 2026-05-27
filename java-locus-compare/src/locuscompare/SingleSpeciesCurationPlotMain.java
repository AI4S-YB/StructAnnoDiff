package locuscompare;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

public final class SingleSpeciesCurationPlotMain {
    private static final Color DELETED = new Color(182, 90, 90);
    private static final Color NEW = new Color(79, 141, 186);
    private static final Color SPLIT = new Color(122, 107, 176);
    private static final Color MERGE = new Color(199, 140, 61);
    private static final Color FILE1 = new Color(77, 77, 77);
    private static final Color FILE2 = new Color(47, 108, 153);
    private static final Color AXIS = new Color(55, 55, 55);
    private static final Color GRID = new Color(230, 230, 230);
    private static final Color TEXT = new Color(34, 34, 34);

    private SingleSpeciesCurationPlotMain() {
    }

    public static void main(String[] args) throws Exception {
        Map<String, String> opts = LocusCompareMain.parseArgs(args);
        Path input = Path.of(opts.getOrDefault("--input", "results/curation_core_metrics.csv"));
        String speciesId = opts.getOrDefault("--species-id", "Pineapple");
        Path analysisDir = Path.of(opts.getOrDefault("--analysis-dir", "..")).toAbsolutePath().normalize();
        Path outputPrefix = resolveOutputPrefix(opts, speciesId);
        Path jigplotJar = Path.of(opts.getOrDefault(
                "--jigplot-jar", defaultJigplotJar())).toAbsolutePath().normalize();

        if (!Files.exists(jigplotJar)) {
            throw new IllegalArgumentException("JIGplot jar not found: " + jigplotJar
                    + ". Use --jigplot-jar or ant -Djigplot.jar=<path-to-JIGplot.jar>.");
        }

        SpeciesMetrics metrics = readMetrics(input, speciesId, analysisDir);
        Files.createDirectories(outputPrefix.toAbsolutePath().getParent());
        try (JigplotRenderer renderer = new JigplotRenderer(jigplotJar)) {
            renderer.render(metrics, outputPrefix);
        }
    }

    private static Path resolveOutputPrefix(Map<String, String> opts, String speciesId) {
        if (opts.containsKey("--output-prefix")) {
            return Path.of(opts.get("--output-prefix"));
        }
        if (opts.containsKey("--output")) {
            return stripKnownImageSuffix(Path.of(opts.get("--output")));
        }
        return Path.of("results/figures/" + speciesId + "_curation_summary_jigplot");
    }

    private static Path stripKnownImageSuffix(Path path) {
        String name = path.getFileName().toString();
        String lower = name.toLowerCase(Locale.ROOT);
        for (String suffix : List.of(".png", ".pdf", ".svg")) {
            if (lower.endsWith(suffix)) {
                String stripped = name.substring(0, name.length() - suffix.length());
                Path parent = path.getParent();
                return parent == null ? Path.of(stripped) : parent.resolve(stripped);
            }
        }
        return path;
    }

    private static String defaultJigplotJar() {
        String env = System.getenv("JIGPLOT_JAR");
        if (env != null && !env.isBlank()) {
            return env;
        }
        return "/mnt/d/NetbeansProject2/TBtools/dist/lib/JIGplot.jar";
    }

    private static SpeciesMetrics readMetrics(Path input, String speciesId, Path analysisDir) throws Exception {
        List<LinkedHashMap<String, String>> rows = Core.readCsv(input);
        LinkedHashMap<String, String> selected = null;
        for (LinkedHashMap<String, String> row : rows) {
            if (speciesId.equals(row.get("species_id"))) {
                if (selected != null) {
                    throw new IllegalArgumentException(input + ": duplicated species_id: " + speciesId);
                }
                selected = row;
            }
        }
        if (selected == null) {
            throw new IllegalArgumentException(input + ": species_id not found: " + speciesId);
        }
        SpeciesMetrics m = new SpeciesMetrics();
        m.speciesId = speciesId;
        m.speciesLabel = selected.getOrDefault("Species", speciesId);
        Path before = Core.findAnnotation(speciesId, "before", analysisDir);
        Path after = Core.findAnnotation(speciesId, "after", analysisDir);
        m.inputFile1 = before == null ? speciesId + ".before" : before.getFileName().toString();
        m.inputFile2 = after == null ? speciesId + ".after" : after.getFileName().toString();
        m.totalBefore = Core.intValue(selected, "total_before_genes");
        m.totalAfter = Core.intValue(selected, "total_after_genes");
        m.deletedLoci = Core.intValue(selected, "deleted_loci_no_overlap");
        m.newLoci = Core.intValue(selected, "new_loci_no_overlap");
        m.splitEvents = Core.intValue(selected, "split_events");
        m.mergeEvents = Core.intValue(selected, "merge_events");
        m.exonChanged = Core.intValue(selected, "rep_exon_changed");
        m.exonChangedFile1Pct = Double.parseDouble(selected.get("rep_exon_changed_before_pct"));
        m.exonChangedFile2Pct = Double.parseDouble(selected.get("rep_exon_changed_after_pct"));
        return m;
    }

    private static String countLabel(int value) {
        return String.format(Locale.ROOT, "%,d", value);
    }

    private static String ratioLabel(double pct, int count, int total) {
        return String.format(Locale.ROOT, "%.2f%% (%s/%s)", pct, countLabel(count), countLabel(total));
    }

    private static String tickLabel(double value) {
        if (value >= 1000.0) {
            return String.format(Locale.ROOT, "%.0fK", value / 1000.0);
        }
        return String.format(Locale.ROOT, "%.0f", value);
    }

    private static double roundUp(double value) {
        if (value <= 1.0) {
            return 1.0;
        }
        double exponent = Math.floor(Math.log10(value));
        double base = Math.pow(10.0, exponent);
        double scaled = value / base;
        double nice;
        if (scaled <= 2.0) {
            nice = 2.0;
        } else if (scaled <= 5.0) {
            nice = 5.0;
        } else {
            nice = 10.0;
        }
        return nice * base;
    }

    private static final class SpeciesMetrics {
        String speciesId;
        String speciesLabel;
        String inputFile1;
        String inputFile2;
        int totalBefore;
        int totalAfter;
        int deletedLoci;
        int newLoci;
        int splitEvents;
        int mergeEvents;
        int exonChanged;
        double exonChangedFile1Pct;
        double exonChangedFile2Pct;
    }

    private static final class JigplotRenderer implements AutoCloseable {
        private static final int WIDTH = 1200;
        private static final int HEIGHT = 1450;
        private static final int BAR_LEFT = 320;
        private static final int BAR_WIDTH = 700;
        private static final int LABEL_X = 95;
        private static final int VALUE_PAD = 14;

        private final URLClassLoader loader;
        private final Class<?> basePanelClass;
        private final Class<?> subPanelClass;
        private final Class<?> elementClass;
        private final Class<?> scalarClass;
        private final Class<?> elementTypeClass;
        private final Constructor<?> basePanelCtor;
        private final Constructor<?> subPanelCtor;
        private final Constructor<?> elementCtor;
        private final Constructor<?> scalarCtor;
        private final Method setSize;
        private final Method setPreferredSize;
        private final Method addSubPanel;
        private final Method addElement;
        private final Method setPanelScalarX;
        private final Method setPanelScalarY;
        private final Method save2Png;
        private final Method save2Pdf;
        private final Method save2Svg;
        private final Method setPoints;
        private final Method setElementType;
        private final Method setFillColor;
        private final Method setDrawColor;
        private final Method setDrawStroke;
        private final Method setText;
        private final Method setFont;
        private final Method setRelativePosArr;
        private final Object rectangleType;
        private final Object lineType;
        private final Object textType;

        JigplotRenderer(Path jigplotJar) throws Exception {
            this.loader = new URLClassLoader(classpathUrls(jigplotJar), ClassLoader.getSystemClassLoader());
            this.basePanelClass = Class.forName("jigplot.engine.JIGBasePanel", true, loader);
            this.subPanelClass = Class.forName("jigplot.engine.JIGSubPanel", true, loader);
            this.elementClass = Class.forName("jigplot.engine.JIGElement", true, loader);
            this.scalarClass = Class.forName("jigplot.OtherTools.Scalar", true, loader);
            this.elementTypeClass = Class.forName("jigplot.engine.JIGConstants$ElementType", true, loader);

            this.basePanelCtor = basePanelClass.getConstructor(int.class, int.class);
            this.subPanelCtor = subPanelClass.getConstructor();
            this.elementCtor = elementClass.getConstructor();
            this.scalarCtor = scalarClass.getConstructor(double.class, double.class, double.class, double.class);
            this.setSize = basePanelClass.getMethod("setSize", int.class, int.class);
            this.setPreferredSize = basePanelClass.getMethod("setPreferredSize", Dimension.class);
            this.addSubPanel = basePanelClass.getMethod("addSubPanel", subPanelClass);
            this.addElement = subPanelClass.getMethod("addElement", elementClass);
            this.setPanelScalarX = subPanelClass.getMethod("setPanelScalarX", scalarClass);
            this.setPanelScalarY = subPanelClass.getMethod("setPanelScalarY", scalarClass);
            this.save2Png = basePanelClass.getMethod("save2PNG", File.class);
            this.save2Pdf = basePanelClass.getMethod("save2PDF", File.class);
            this.save2Svg = basePanelClass.getMethod("save2SVG", File.class);
            this.setPoints = elementClass.getMethod("setPoints", Point2D[].class);
            this.setElementType = elementClass.getMethod("setElementType", elementTypeClass);
            this.setFillColor = elementClass.getMethod("setFillColor", Color.class);
            this.setDrawColor = elementClass.getMethod("setDrawColor", Color.class);
            this.setDrawStroke = elementClass.getMethod("setDrawStroke", Float.class);
            this.setText = elementClass.getMethod("setText", String.class);
            this.setFont = elementClass.getMethod("setFont", Font.class);
            this.setRelativePosArr = elementClass.getMethod("setRelativePosArr", Float[].class);
            this.rectangleType = enumValue("Rectangle");
            this.lineType = enumValue("Line");
            this.textType = enumValue("Text");
        }

        private static URL[] classpathUrls(Path jigplotJar) throws Exception {
            List<URL> urls = new ArrayList<>();
            urls.add(jigplotJar.toUri().toURL());
            Path libDir = jigplotJar.getParent();
            for (String name : List.of("JFreeSVG_Ant.jar", "OrsonPDF_Ant.jar")) {
                Path candidate = libDir == null ? null : libDir.resolve(name);
                if (candidate != null && Files.exists(candidate)) {
                    urls.add(candidate.toUri().toURL());
                }
            }
            return urls.toArray(new URL[0]);
        }

        @SuppressWarnings({"unchecked", "rawtypes"})
        private Object enumValue(String name) {
            return Enum.valueOf((Class<? extends Enum>) elementTypeClass.asSubclass(Enum.class), name);
        }

        void render(SpeciesMetrics m, Path outputPrefix) throws Exception {
            Object panel = buildPanel(m);
            Path png = withSuffix(outputPrefix, ".png");
            Path pdf = withSuffix(outputPrefix, ".pdf");
            Path svg = withSuffix(outputPrefix, ".svg");
            save2Png.invoke(panel, png.toFile());
            save2Pdf.invoke(panel, pdf.toFile());
            save2Svg.invoke(panel, svg.toFile());
            System.out.println("Saved: " + png);
            System.out.println("Saved: " + pdf);
            System.out.println("Saved: " + svg);
        }

        private static Path withSuffix(Path prefix, String suffix) {
            Path parent = prefix.getParent();
            String name = prefix.getFileName().toString() + suffix;
            return parent == null ? Path.of(name) : parent.resolve(name);
        }

        private Object buildPanel(SpeciesMetrics m) throws Exception {
            Object panel = basePanelCtor.newInstance(WIDTH, HEIGHT);
            setSize.invoke(panel, WIDTH + 40, HEIGHT + 40);
            setPreferredSize.invoke(panel, new Dimension(WIDTH + 40, HEIGHT + 40));
            Object subPanel = subPanelCtor.newInstance();
            setPoints.invoke(subPanel, (Object) new Point2D[]{
                    new Point2D.Double(0, 0), new Point2D.Double(WIDTH, HEIGHT)
            });
            Object xScalar = scalarCtor.newInstance(0.0, (double) WIDTH, 0.0, (double) WIDTH);
            Object yScalar = scalarCtor.newInstance(0.0, (double) HEIGHT, 0.0, (double) HEIGHT);
            setPanelScalarX.invoke(subPanel, xScalar);
            setPanelScalarY.invoke(subPanel, yScalar);

            title(subPanel, 70, m.speciesLabel + " annotation curation summary", 28, Font.BOLD);
            title(subPanel, 112, m.inputFile1 + " -> " + m.inputFile2, 17, Font.PLAIN);
            drawCountPanel(subPanel, 230, "A", "Locus gain/loss",
                    new String[]{"Deleted loci", "New loci"},
                    new int[]{m.deletedLoci, m.newLoci},
                    new Color[]{DELETED, NEW},
                    "No-overlap loci");
            drawCountPanel(subPanel, 620, "B", "Split and merge events",
                    new String[]{"Split events", "Merge events"},
                    new int[]{m.splitEvents, m.mergeEvents},
                    new Color[]{SPLIT, MERGE},
                    "Events");
            drawPercentPanel(subPanel, 1010, m);

            addSubPanel.invoke(panel, subPanel);
            return panel;
        }

        private void drawCountPanel(Object subPanel, int top, String panelLabel, String title,
                                    String[] labels, int[] values, Color[] colors, String axisLabel) throws Exception {
            double xMax = roundUp(Math.max(values[0], values[1]) * 1.18);
            drawPanelHeader(subPanel, top, panelLabel, title);
            drawAxis(subPanel, top, xMax, false, axisLabel);
            int[] ys = new int[]{top + 120, top + 195};
            for (int i = 0; i < labels.length; i++) {
                text(subPanel, LABEL_X, ys[i] + 6, labels[i], 18, Font.PLAIN, TEXT, 0.0f, 0.5f);
                double width = values[i] / xMax * BAR_WIDTH;
                rect(subPanel, BAR_LEFT, ys[i] - 22, BAR_LEFT + width, ys[i] + 22, colors[i]);
                text(subPanel, BAR_LEFT + width + VALUE_PAD, ys[i] + 5, countLabel(values[i]),
                        17, Font.PLAIN, TEXT, 0.0f, 0.5f);
            }
        }

        private void drawPercentPanel(Object subPanel, int top, SpeciesMetrics m) throws Exception {
            drawPanelHeader(subPanel, top, "C", "Representative transcript exon changes");
            drawAxis(subPanel, top, 100.0, true, "Representative exon changes (%)");
            String[] labels = new String[]{"Input file 1", "Input file 2"};
            double[] values = new double[]{m.exonChangedFile1Pct, m.exonChangedFile2Pct};
            int[] totals = new int[]{m.totalBefore, m.totalAfter};
            Color[] colors = new Color[]{FILE1, FILE2};
            int[] ys = new int[]{top + 120, top + 195};
            for (int i = 0; i < labels.length; i++) {
                text(subPanel, LABEL_X, ys[i] + 6, labels[i], 18, Font.PLAIN, TEXT, 0.0f, 0.5f);
                double width = values[i] / 100.0 * BAR_WIDTH;
                rect(subPanel, BAR_LEFT, ys[i] - 22, BAR_LEFT + width, ys[i] + 22, colors[i]);
                String label = ratioLabel(values[i], m.exonChanged, totals[i]);
                if (values[i] >= 52.0) {
                    text(subPanel, BAR_LEFT + width - VALUE_PAD, ys[i] + 5, label,
                            17, Font.PLAIN, Color.WHITE, -1.0f, 0.5f);
                } else {
                    text(subPanel, BAR_LEFT + width + VALUE_PAD, ys[i] + 5, label,
                            17, Font.PLAIN, TEXT, 0.0f, 0.5f);
                }
            }
        }

        private void drawPanelHeader(Object subPanel, int top, String label, String title) throws Exception {
            text(subPanel, 80, top - 42, label, 24, Font.BOLD, TEXT, 0.0f, 0.5f);
            text(subPanel, 128, top - 42, title, 21, Font.BOLD, TEXT, 0.0f, 0.5f);
        }

        private void drawAxis(Object subPanel, int top, double xMax, boolean percent, String axisLabel) throws Exception {
            int axisY = top + 255;
            line(subPanel, BAR_LEFT, axisY, BAR_LEFT + BAR_WIDTH, axisY, AXIS, 1.2f);
            int tickCount = percent ? 4 : 3;
            DecimalFormat df = new DecimalFormat(percent ? "0" : "0.##");
            for (int i = 0; i <= tickCount; i++) {
                double v = xMax * i / tickCount;
                double x = BAR_LEFT + BAR_WIDTH * i / (double) tickCount;
                line(subPanel, x, top + 70, x, axisY, GRID, 0.8f);
                line(subPanel, x, axisY, x, axisY + 8, AXIS, 1.0f);
                String label = percent ? df.format(v) : tickLabel(v);
                text(subPanel, x, axisY + 34, label, 15, Font.PLAIN, TEXT, -0.5f, 0.5f);
            }
            text(subPanel, BAR_LEFT, axisY + 70, axisLabel, 16, Font.PLAIN, TEXT, 0.0f, 0.5f);
        }

        private void title(Object subPanel, int y, String s, int size, int style) throws Exception {
            text(subPanel, 80, y, s, size, style, TEXT, 0.0f, 0.5f);
        }

        private Object element(String type) throws Exception {
            Object element = elementCtor.newInstance();
            if ("Rectangle".equals(type)) {
                setElementType.invoke(element, rectangleType);
            } else if ("Line".equals(type)) {
                setElementType.invoke(element, lineType);
            } else if ("Text".equals(type)) {
                setElementType.invoke(element, textType);
            }
            return element;
        }

        private void rect(Object subPanel, double x1, double y1, double x2, double y2, Color fill) throws Exception {
            Object element = element("Rectangle");
            setPoints.invoke(element, (Object) new Point2D[]{
                    new Point2D.Double(x1, y1), new Point2D.Double(x2, y2)
            });
            setFillColor.invoke(element, fill);
            setDrawColor.invoke(element, fill.darker());
            setDrawStroke.invoke(element, 0.7f);
            addElement.invoke(subPanel, element);
        }

        private void line(Object subPanel, double x1, double y1, double x2, double y2,
                          Color color, float stroke) throws Exception {
            Object element = element("Line");
            setPoints.invoke(element, (Object) new Point2D[]{
                    new Point2D.Double(x1, y1), new Point2D.Double(x2, y2)
            });
            setDrawColor.invoke(element, color);
            setDrawStroke.invoke(element, stroke);
            addElement.invoke(subPanel, element);
        }

        private void text(Object subPanel, double x, double y, String value, int size, int style,
                          Color color, float relX, float relY) throws Exception {
            Object element = element("Text");
            setPoints.invoke(element, (Object) new Point2D[]{new Point2D.Double(x, y)});
            setText.invoke(element, value);
            setFont.invoke(element, new Font("Arial", style, size));
            setDrawColor.invoke(element, color);
            setRelativePosArr.invoke(element, (Object) new Float[]{relX, relY});
            addElement.invoke(subPanel, element);
        }

        @Override
        public void close() throws Exception {
            loader.close();
        }
    }
}
