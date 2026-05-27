#!/usr/bin/env python3
"""Plot core curation metrics with rows ordered by a supplied species tree."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.ticker as mticker

from analysis_config import FIGURES_DIR, RESULTS_DIR
from plot_curation_core_figure import (
    COLORS,
    _compact_legend,
    _panel_header,
    _set_shared_y,
    _style_axis,
    load_metrics,
    plot_locus_gain_loss,
    plot_split_merge,
)


DEFAULT_NEWICK = (
    "((Artemisia_annua:125.11960000,(Cucumis_sativus:102.61000000,"
    "((Fragaria_ananassa:7.01170000,Fragaria_vesca:7.01170000):62.05377000,"
    "Prunus_persica:69.06547000):33.54453000):22.50960000):34.49316000,"
    "(Ananas_comosus:110.07878000,Oryza_sativa:110.07878000):49.53398000)"
)

TREE_NAME_TO_SPECIES_ID = {
    "Artemisia_annua": "Artemisia_annua",
    "Cucumis_sativus": "Cucumber",
    "Fragaria_ananassa": "Fragaria_ananassa",
    "Fragaria_vesca": "Fragaria_vesca",
    "Prunus_persica": "Peach",
    "Ananas_comosus": "Pineapple",
    "Oryza_sativa": "Rice",
}

TREE_DISPLAY_LABELS = {
    "Artemisia_annua": "A. annua",
    "Cucumis_sativus": "C. sativus",
    "Fragaria_ananassa": "F. ananassa",
    "Fragaria_vesca": "F. vesca",
    "Prunus_persica": "P. persica",
    "Ananas_comosus": "A. comosus",
    "Oryza_sativa": "O. sativa",
}

QUALITY_COLUMNS = [
    "species",
    "before_busco_complete_pct",
    "after_busco_complete_pct",
    "before_psauron_overall_score",
    "after_psauron_overall_score",
]


@dataclass
class TreeNode:
    name: str = ""
    length: float = 0.0
    children: list["TreeNode"] = field(default_factory=list)
    x: float = 0.0
    y: float = 0.0

    @property
    def is_leaf(self) -> bool:
        return not self.children


class NewickParser:
    def __init__(self, text: str):
        self.text = text.strip()
        self.pos = 0

    def parse(self) -> TreeNode:
        node = self._parse_subtree()
        self._skip_ws()
        if self._peek() == ";":
            self.pos += 1
            self._skip_ws()
        if self.pos != len(self.text):
            raise ValueError(f"Unexpected Newick content at position {self.pos}: {self.text[self.pos:]!r}")
        return node

    def _parse_subtree(self) -> TreeNode:
        self._skip_ws()
        if self._peek() == "(":
            self.pos += 1
            children = []
            while True:
                children.append(self._parse_subtree())
                self._skip_ws()
                char = self._peek()
                if char == ",":
                    self.pos += 1
                    continue
                if char == ")":
                    self.pos += 1
                    break
                raise ValueError(f"Expected ',' or ')' at position {self.pos}")
            name = self._parse_name()
            length = self._parse_length()
            return TreeNode(name=name, length=length, children=children)

        name = self._parse_name()
        if not name:
            raise ValueError(f"Expected leaf name at position {self.pos}")
        length = self._parse_length()
        return TreeNode(name=name, length=length)

    def _parse_name(self) -> str:
        self._skip_ws()
        start = self.pos
        while self.pos < len(self.text) and self.text[self.pos] not in ":,();":
            self.pos += 1
        return self.text[start:self.pos].strip()

    def _parse_length(self) -> float:
        self._skip_ws()
        if self._peek() != ":":
            return 0.0
        self.pos += 1
        self._skip_ws()
        start = self.pos
        while self.pos < len(self.text) and self.text[self.pos] not in ",();":
            self.pos += 1
        value = self.text[start:self.pos].strip()
        if not value:
            raise ValueError(f"Expected branch length at position {start}")
        return float(value)

    def _peek(self) -> str:
        if self.pos >= len(self.text):
            return ""
        return self.text[self.pos]

    def _skip_ws(self) -> None:
        while self.pos < len(self.text) and self.text[self.pos].isspace():
            self.pos += 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot all-species curation metrics ordered by a Newick species tree."
    )
    parser.add_argument(
        "--input",
        default=str(RESULTS_DIR / "curation_core_metrics.csv"),
        help="Input curation core metrics CSV.",
    )
    parser.add_argument(
        "--newick",
        default=DEFAULT_NEWICK,
        help="Newick tree text. Leaf names must be known in TREE_NAME_TO_SPECIES_ID.",
    )
    parser.add_argument(
        "--output-prefix",
        default=str(FIGURES_DIR / "curation_core_metrics_tree_ordered"),
        help="Output path prefix without extension.",
    )
    parser.add_argument(
        "--quality-input",
        default=str(RESULTS_DIR / "busco_psauron_metrics.tsv"),
        help="Input TSV with BUSCO and Psauron before/after metrics.",
    )
    parser.add_argument("--dpi", type=int, default=600, help="PNG output resolution.")
    parser.add_argument(
        "--width-scale",
        type=float,
        default=1.0,
        help="Scale the physical width of every panel; use 0.333 for the earlier compressed layout.",
    )
    parser.add_argument(
        "--font-scale",
        type=float,
        default=1.0,
        help="Deprecated compatibility option; use --font-size for the final text size.",
    )
    parser.add_argument(
        "--font-size",
        type=float,
        default=14.0,
        help="Final font size in points for all text in the figure.",
    )
    parser.add_argument(
        "--line-width",
        type=float,
        default=0.35,
        help="Thin line width in points for axes, grids, and the species tree.",
    )
    parser.add_argument(
        "--height",
        type=float,
        default=7.2,
        help="Figure height in inches.",
    )
    return parser.parse_args()


def load_quality_metrics(path: Path, expected_species: set[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    missing = [col for col in QUALITY_COLUMNS if col not in df.columns]
    if missing:
        raise SystemExit(f"{path}: missing required columns: {', '.join(missing)}")

    if df["species"].duplicated().any():
        duplicates = sorted(df.loc[df["species"].duplicated(), "species"].unique())
        raise SystemExit(f"{path}: duplicate species rows: {', '.join(duplicates)}")

    observed = set(df["species"])
    if observed != expected_species:
        raise SystemExit(
            f"{path}: species mismatch; missing={sorted(expected_species - observed)}, "
            f"extra={sorted(observed - expected_species)}"
        )

    numeric_cols = [col for col in QUALITY_COLUMNS if col != "species"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="raise")
        if ((df[col] < 0) | (df[col] > 100)).any():
            raise SystemExit(f"{path}: {col} must be within 0-100")

    return df


def iter_leaves(node: TreeNode) -> list[TreeNode]:
    if node.is_leaf:
        return [node]
    leaves: list[TreeNode] = []
    for child in node.children:
        leaves.extend(iter_leaves(child))
    return leaves


def assign_tree_coordinates(root: TreeNode, leaf_names: list[str]) -> float:
    leaf_y = {name: float(index) for index, name in enumerate(leaf_names)}

    def walk(node: TreeNode, parent_x: float) -> None:
        node.x = parent_x + node.length
        if node.is_leaf:
            node.y = leaf_y[node.name]
            return
        for child in node.children:
            walk(child, node.x)
        node.y = float(np.mean([child.y for child in node.children]))

    # The root length in user-supplied Newick is absent; ignore it if present so
    # the drawing starts at zero.
    root.x = 0.0
    if root.is_leaf:
        root.y = leaf_y[root.name]
    else:
        for child in root.children:
            walk(child, 0.0)
        root.y = float(np.mean([child.y for child in root.children]))
    return max(leaf.x for leaf in iter_leaves(root))


def reorder_metrics_by_tree(df: pd.DataFrame, root: TreeNode) -> tuple[pd.DataFrame, list[str], list[str]]:
    leaf_names = [leaf.name for leaf in iter_leaves(root)]
    unmapped = [name for name in leaf_names if name not in TREE_NAME_TO_SPECIES_ID]
    if unmapped:
        raise SystemExit(f"Newick contains unmapped species names: {', '.join(unmapped)}")

    species_order = [TREE_NAME_TO_SPECIES_ID[name] for name in leaf_names]
    observed = set(df["species_id"])
    tree_species = set(species_order)
    if observed != tree_species:
        raise SystemExit(
            "Tree/metrics species mismatch; "
            f"missing_in_tree={sorted(observed - tree_species)}, "
            f"missing_in_metrics={sorted(tree_species - observed)}"
        )

    ordered = df.set_index("species_id").loc[species_order].reset_index()
    labels = [TREE_DISPLAY_LABELS.get(name, name.replace("_", " ")) for name in leaf_names]
    return ordered, leaf_names, labels


def draw_species_tree(ax, root: TreeNode, leaf_labels: dict[str, str], max_depth: float) -> None:
    line_color = "#404040"

    def draw_node(node: TreeNode) -> None:
        if node.children:
            child_ys = [child.y for child in node.children]
            ax.plot(
                [node.x, node.x],
                [min(child_ys), max(child_ys)],
                color=line_color,
                linewidth=0.35,
                solid_capstyle="round",
            )
            for child in node.children:
                ax.plot(
                    [node.x, child.x],
                    [child.y, child.y],
                    color=line_color,
                    linewidth=0.35,
                    solid_capstyle="round",
                )
                draw_node(child)
            return

        ax.text(
            max_depth * 1.025,
            node.y,
            leaf_labels[node.name],
            ha="left",
            va="center",
            fontsize=6.9,
            color="#202020",
        )

    draw_node(root)
    ax.set_xlim(-max_depth * 0.035, max_depth * 1.52)
    ax.set_ylim(len(leaf_labels) - 0.45, -1.25)
    ax.axis("off")
    ax.text(
        0.0,
        1.055,
        "Species tree",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.8,
        fontweight="bold",
    )


def plot_quality_panel(
    ax,
    df: pd.DataFrame,
    y: np.ndarray,
    labels: list[str],
    before_col: str,
    after_col: str,
    panel_label: str,
    title: str,
    xlabel: str,
    suffix: str,
    x_min: float = 50.0,
    x_max: float = 100.0,
    show_labels: bool = False,
) -> None:
    offset = 0.16
    height = 0.27
    before = df[before_col].to_numpy()
    after = df[after_col].to_numpy()

    ax.barh(
        y - offset,
        before,
        height=height,
        color=COLORS["before_ref"],
        label="Before",
    )
    ax.barh(
        y + offset,
        after,
        height=height,
        color=COLORS["after_ref"],
        label="After",
    )
    ax.set_xlim(x_min, x_max)
    ax.set_xlabel(xlabel, fontsize=7.8)
    ax.xaxis.set_major_locator(mticker.FixedLocator([50, 75, 100]))
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda value, _pos: f"{value:.0f}"))
    _set_shared_y(ax, y, labels, show_labels=show_labels)
    _style_axis(ax, show_y=show_labels)
    _panel_header(ax, panel_label, title)
    _compact_legend(ax, ncol=2)

    for yi, value in zip(y - offset, before):
        ax.text(
            value - 1.2,
            yi,
            f"{value:.1f}{suffix}",
            ha="right",
            va="center",
            fontsize=5.2,
            color="white",
            clip_on=True,
        )
    for yi, value in zip(y + offset, after):
        ax.text(
            value - 1.2,
            yi,
            f"{value:.1f}{suffix}",
            ha="right",
            va="center",
            fontsize=5.2,
            color="white",
            clip_on=True,
        )


def plot_representative_exon_changes_compact(
    ax,
    df: pd.DataFrame,
    y: np.ndarray,
    labels: list[str],
    panel_label: str = "C",
    show_labels: bool = False,
) -> None:
    offset = 0.16
    height = 0.27
    before_pct = df["rep_exon_changed_before_pct"].to_numpy()
    after_pct = df["rep_exon_changed_after_pct"].to_numpy()

    ax.barh(
        y - offset,
        before_pct,
        height=height,
        color=COLORS["before_ref"],
        label="Before ref.",
    )
    ax.barh(
        y + offset,
        after_pct,
        height=height,
        color=COLORS["after_ref"],
        label="After ref.",
    )
    ax.set_xlim(0, 100)
    ax.set_xlabel("Exon change (%)", fontsize=7.8)
    ax.xaxis.set_major_locator(mticker.FixedLocator([0, 50, 100]))
    ax.xaxis.set_major_formatter(mticker.PercentFormatter(xmax=100, decimals=0))
    _set_shared_y(ax, y, labels, show_labels=show_labels)
    _style_axis(ax, show_y=show_labels)
    _panel_header(ax, panel_label, "Exon change")
    _compact_legend(ax, ncol=2)

    for yi, value in zip(y - offset, before_pct):
        color = "white" if value >= 45 else "#222222"
        xpos = value - 1.2 if value >= 45 else value + 1.2
        ha = "right" if value >= 45 else "left"
        ax.text(xpos, yi, f"{value:.1f}%", ha=ha, va="center", fontsize=5.2, color=color, clip_on=False)
    for yi, value in zip(y + offset, after_pct):
        color = "white" if value >= 45 else "#222222"
        xpos = value - 1.2 if value >= 45 else value + 1.2
        ha = "right" if value >= 45 else "left"
        ax.text(xpos, yi, f"{value:.1f}%", ha=ha, va="center", fontsize=5.2, color=color, clip_on=False)


def combine_panel_header_text(axes) -> None:
    panel_labels = {"A", "B", "C", "D", "E"}
    for ax in axes:
        label_text = None
        title_text = None
        for text in ax.texts:
            if text.get_transform() != ax.transAxes:
                continue
            x, y = text.get_position()
            if abs(y - 1.055) > 0.01:
                continue
            if text.get_text() in panel_labels:
                label_text = text
            else:
                title_text = text
        if label_text is None or title_text is None:
            continue
        label_text.set_text(f"{label_text.get_text()} {title_text.get_text()}")
        label_text.set_position((0.0, label_text.get_position()[1]))
        title_text.set_visible(False)


def shorten_panel_text(axes) -> None:
    replacements = {
        "Locus gain/loss": "Gain/loss",
        "Split and merge events": "Split/merge",
        "BUSCO complete": "BUSCO",
        "Psauron score": "Psauron",
        "No-overlap loci": "Loci",
        "Deleted loci": "Del.",
        "New loci": "New",
        "Before ref.": "Before",
        "After ref.": "After",
    }
    for ax in axes:
        xlabel = ax.get_xlabel()
        if xlabel in replacements:
            ax.set_xlabel(replacements[xlabel])
        legend = ax.get_legend()
        if legend is not None:
            for legend_text in legend.get_texts():
                replacement = replacements.get(legend_text.get_text())
                if replacement:
                    legend_text.set_text(replacement)
        for text in ax.texts:
            replacement = replacements.get(text.get_text())
            if replacement:
                text.set_text(replacement)


def apply_global_text_and_line_style(fig, axes, font_size: float, line_width: float) -> None:
    from matplotlib.text import Text

    for text in fig.findobj(match=Text):
        text.set_fontsize(font_size)

    for ax in axes:
        for spine in ax.spines.values():
            spine.set_linewidth(line_width)
        ax.tick_params(axis="both", width=line_width, length=max(1.4, line_width * 5.0))
        for line in ax.lines:
            line.set_linewidth(line_width)
        for gridline in ax.get_xgridlines() + ax.get_ygridlines():
            gridline.set_linewidth(line_width)
        for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(line_width)
            tick.tick2line.set_markeredgewidth(line_width)
        for patch in ax.patches:
            patch.set_linewidth(0.0)


def plot_tree_ordered_figure(
    df: pd.DataFrame,
    root: TreeNode,
    leaf_names: list[str],
    labels: list[str],
    output_prefix: Path,
    dpi: int,
    width_scale: float,
    font_scale: float,
    font_size: float,
    line_width: float,
    height: float,
) -> None:
    import matplotlib.pyplot as plt

    y = np.arange(len(df))
    max_depth = assign_tree_coordinates(root, leaf_names)
    leaf_labels = dict(zip(leaf_names, labels))

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": font_size,
            "axes.linewidth": line_width,
        }
    )

    width_ratios = [1.6, 0.875, 0.625, 1.0, 0.65, 0.65]
    base_unit_width = 12.15 / sum([1.38, 1.32, 1.02, 1.52])
    fig, axes = plt.subplots(
        1,
        6,
        figsize=(base_unit_width * sum(width_ratios) * width_scale, height),
        gridspec_kw={"width_ratios": width_ratios, "wspace": 0.38},
    )

    draw_species_tree(axes[0], root, leaf_labels, max_depth)
    plot_locus_gain_loss(axes[1], df, y, labels, panel_label="A", show_labels=False)
    plot_split_merge(axes[2], df, y, labels, panel_label="B", show_labels=False)
    plot_representative_exon_changes_compact(axes[3], df, y, labels, panel_label="C", show_labels=False)
    plot_quality_panel(
        axes[4],
        df,
        y,
        labels,
        before_col="before_busco_complete_pct",
        after_col="after_busco_complete_pct",
        panel_label="D",
        title="BUSCO complete",
        xlabel="BUSCO (%)",
        suffix="%",
        x_min=50.0,
        x_max=100.0,
    )
    plot_quality_panel(
        axes[5],
        df,
        y,
        labels,
        before_col="before_psauron_overall_score",
        after_col="after_psauron_overall_score",
        panel_label="E",
        title="Psauron score",
        xlabel="Score",
        suffix="",
        x_min=50.0,
        x_max=100.0,
    )

    shorten_panel_text(axes[1:])
    axes[1].xaxis.set_major_locator(mticker.FixedLocator([0, 5000, 10000, 15000]))
    axes[2].xaxis.set_major_locator(mticker.FixedLocator([0, 2000, 4000]))
    combine_panel_header_text(axes[1:])
    if font_scale != 1.0:
        font_size *= font_scale
    apply_global_text_and_line_style(fig, axes, font_size, line_width)

    fig.subplots_adjust(left=0.035, right=0.99, top=0.84, bottom=0.16)

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = output_prefix.with_suffix(".pdf")
    svg_path = output_prefix.with_suffix(".svg")
    png_path = output_prefix.with_suffix(".png")
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {pdf_path}")
    print(f"Saved: {svg_path}")
    print(f"Saved: {png_path}")


def main() -> None:
    args = parse_args()
    root = NewickParser(args.newick).parse()
    df = load_metrics(Path(args.input))
    quality_df = load_quality_metrics(Path(args.quality_input), set(df["species_id"]))
    ordered_df, leaf_names, labels = reorder_metrics_by_tree(df, root)
    ordered_df = ordered_df.merge(
        quality_df.rename(columns={"species": "species_id"}),
        on="species_id",
        how="left",
        validate="one_to_one",
    )
    plot_tree_ordered_figure(
        ordered_df,
        root,
        leaf_names,
        labels,
        Path(args.output_prefix),
        args.dpi,
        args.width_scale,
        args.font_scale,
        args.font_size,
        args.line_width,
        args.height,
    )


if __name__ == "__main__":
    main()
