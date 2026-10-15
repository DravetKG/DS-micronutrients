"""Figure of the micronutrient gene set counts and their overlap.

Panel A gives the number of genes of each set. Panel B is an UpSet plot of the
intersections, because a Venn diagram of ten sets is not readable.
"""

from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ORDER = ["Vitamin A", "Vitamin D", "Vitamin C", "Vitamin E", "Vitamin B6",
         "Vitamin B12", "Iron", "Folic Acid", "Zinc", "Selenium"]

COLORS = {
    "Vitamin A": "#c2453b",
    "Vitamin D": "#2a78d6",
    "Vitamin C": "#eb6834",
    "Vitamin E": "#1baf7a",
    "Vitamin B6": "#eda100",
    "Vitamin B12": "#e87ba4",
    "Iron": "#8c6d46",
    "Folic Acid": "#38a3a5",
    "Zinc": "#7b5ea7",
    "Selenium": "#6b8e23",
    "Keto diet": "#4a3aa7",
    "Dravet": "#0b0b0b",
    "CAKUT": "#7a7974",
}

INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#dcdbd6"
DOT_OFF = "#dcdbd6"
DOT_ON = "#0b0b0b"
GROUND = "#fcfcfb"

MAX_ROWS = 40


def intersection_sizes(sets, order=None, max_rows=None):
    """Count the genes that belong to each exact combination of the sets."""
    order = order or ORDER
    counts = Counter()
    universe = set().union(*sets.values()) if sets else set()
    for gene in universe:
        counts[tuple(name for name in order if gene in sets.get(name, ()))] += 1
    rows = sorted(counts.items(), key=lambda item: (-item[1], len(item[0])))
    return rows[:max_rows] if max_rows else rows


def style_axes(axes):
    """Remove the top and right spines and make the axes recessive."""
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axes.spines[side].set_color(GRID)
    axes.tick_params(colors=INK_MUTED, labelsize=8, length=3)


def panel_counts(axes, sets, order):
    """Draw the bar chart of the gene count of each micronutrient."""
    values = [len(sets.get(name, ())) for name in order]
    positions = np.arange(len(order))
    axes.bar(positions, values, width=0.62,
             color=[COLORS[name] for name in order], linewidth=0)
    for position, value in zip(positions, values):
        axes.text(position, value + max(values) * 0.02, str(value),
                  ha="center", va="bottom", fontsize=8, color=INK)
    axes.set_xticks(positions)
    axes.set_xticklabels(order, rotation=30, ha="right")
    axes.set_ylabel("Genes", fontsize=9, color=INK_MUTED)
    axes.set_ylim(0, max(values) * 1.15)
    axes.yaxis.grid(True, color=GRID, linewidth=0.6)
    axes.set_axisbelow(True)
    style_axes(axes)
    axes.set_title("A  Genes per gene set", loc="left", fontsize=10,
                   color=INK, pad=10)


def panel_upset(axes_bars, axes_dots, rows, order, title):
    """Draw the UpSet plot of the exact set intersections."""
    positions = np.arange(len(rows))
    values = [size for _, size in rows]
    colors = [COLORS[members[0]] if len(members) == 1 else INK_MUTED
              for members, _ in rows]
    axes_bars.bar(positions, values, width=0.6, color=colors, linewidth=0)
    for position, value in zip(positions, values):
        axes_bars.text(position, value + max(values) * 0.02, str(value),
                       ha="center", va="bottom", fontsize=7.5, color=INK)
    axes_bars.set_ylabel("Genes in intersection", fontsize=9, color=INK_MUTED)
    axes_bars.set_ylim(0, max(values) * 1.18)
    axes_bars.set_xticks([])
    axes_bars.set_xlim(-0.7, len(rows) - 0.3)
    axes_bars.yaxis.grid(True, color=GRID, linewidth=0.6)
    axes_bars.set_axisbelow(True)
    style_axes(axes_bars)
    axes_bars.set_title(title, loc="left", fontsize=10, color=INK, pad=10)

    for index in range(0, len(order), 2):
        axes_dots.axhspan(index - 0.5, index + 0.5, color="#f2f1ec",
                          linewidth=0, zorder=0)
    for position, (members, _) in enumerate(rows):
        indices = [order.index(name) for name in members]
        axes_dots.scatter([position] * len(order), range(len(order)), s=42,
                          color=DOT_OFF, zorder=1, linewidths=0)
        axes_dots.scatter([position] * len(indices), indices, s=42,
                          color=DOT_ON, zorder=2, linewidths=0)
        if len(indices) > 1:
            axes_dots.plot([position, position], [min(indices), max(indices)],
                           color=DOT_ON, linewidth=1.6, zorder=2)
    axes_dots.set_xlim(-0.7, len(rows) - 0.3)
    axes_dots.set_ylim(len(order) - 0.5, -0.5)
    axes_dots.set_yticks(range(len(order)))
    axes_dots.set_yticklabels(order, fontsize=8, color=INK_MUTED)
    axes_dots.set_xticks([])
    for side in ("top", "right", "bottom", "left"):
        axes_dots.spines[side].set_visible(False)
    axes_dots.tick_params(length=0)


def make_figure(sets, output_path=None, subtitle="", order=None,
                max_rows=MAX_ROWS, title="CTD micronutrient gene sets"):
    """Draw the two panels and write the figure, keeping the ``max_rows`` largest
    intersections in panel B.
    """
    order = order or ORDER
    rows = intersection_sizes(sets, order, max_rows)
    # The dot matrix gets 0.22 inch per set.
    dots = 0.22 * len(order)
    figure = plt.figure(figsize=(12, 5.8 + dots), facecolor=GROUND)
    grid = figure.add_gridspec(2, 1, height_ratios=[2.6, 3.4 + dots], hspace=0.5)
    axes_counts = figure.add_subplot(grid[0])

    lower = grid[1].subgridspec(2, 1, height_ratios=[3.4, dots], hspace=0.05)
    axes_upset_bars = figure.add_subplot(lower[0])
    axes_upset_dots = figure.add_subplot(lower[1])

    for axes in (axes_counts, axes_upset_bars, axes_upset_dots):
        axes.set_facecolor(GROUND)

    panel_counts(axes_counts, sets, order)
    total = len(intersection_sizes(sets, order))
    upset_title = "B  Overlap of the gene sets"
    if len(rows) < total:
        upset_title += f"  ({len(rows)} of {total} intersections)"
    panel_upset(axes_upset_bars, axes_upset_dots, rows, order, upset_title)

    handles = [Line2D([0], [0], marker="s", linestyle="", markersize=7,
                      color=COLORS[name], label=name) for name in order]
    handles.append(Line2D([0], [0], marker="s", linestyle="", markersize=7,
                          color=INK_MUTED, label="Shared by two or more"))
    figure.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
                  fontsize=8, labelcolor=INK_MUTED, bbox_to_anchor=(0.5, -0.005))

    figure.suptitle(title + subtitle, fontsize=11.5,
                    color=INK, x=0.02, ha="left", y=0.98)
    if output_path:
        for path in (Path(output_path), Path(output_path).with_suffix(".pdf")):
            figure.savefig(path, dpi=300, bbox_inches="tight",
                           facecolor=figure.get_facecolor())
            print(f"Written to {path}")
    return figure


SOURCE_COLORS = {"CTD only": "#2a78d6", "CTD and pathway": "#eb6834",
                 "Pathway only": "#1baf7a"}


def source_figure(ctd, pathway, output_path=None, order=None,
                  title="Seed genes by source"):
    """Draw one stacked bar per micronutrient: CTD only, both, pathway only.

    ``ctd`` and ``pathway`` map each micronutrient to its set of gene symbols.
    """
    order = [name for name in (order or ORDER) if name in ctd]
    parts = {"CTD only": [len(ctd[n] - pathway[n]) for n in order],
             "CTD and pathway": [len(ctd[n] & pathway[n]) for n in order],
             "Pathway only": [len(pathway[n] - ctd[n]) for n in order]}
    positions = np.arange(len(order))
    figure, axes = plt.subplots(figsize=(7.5, 0.42 * len(order) + 1.4),
                                facecolor=GROUND)
    axes.set_facecolor(GROUND)
    left = np.zeros(len(order))
    # Segments narrower than this fraction of the x range get no count label.
    min_width = 0.012 * max(sum(column) for column in zip(*parts.values()))
    for label, values in parts.items():
        axes.barh(positions, values, left=left, height=0.62,
                  color=SOURCE_COLORS[label], edgecolor=GROUND, linewidth=1,
                  label=label)
        for position, start, value in zip(positions, left, values):
            if value >= min_width:
                axes.text(start + value / 2, position, str(value), ha="center",
                          va="center", fontsize=7, color=INK)
        left += values
    for position, total in zip(positions, left):
        axes.text(total + left.max() * 0.01, position, f"{int(total)}",
                  va="center", fontsize=8, color=INK_MUTED)
    axes.set_yticks(positions)
    axes.set_yticklabels(order)
    axes.set_ylim(len(order) - 0.5, -0.5)
    axes.set_xlabel("Protein-coding genes", fontsize=9, color=INK_MUTED)
    axes.xaxis.grid(True, color=GRID, linewidth=0.6)
    axes.set_axisbelow(True)
    style_axes(axes)
    axes.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3,
                frameon=False, fontsize=8, labelcolor=INK_MUTED)
    axes.set_title(title, loc="left", fontsize=10.5, color=INK, pad=24)
    if output_path:
        for path in (Path(output_path), Path(output_path).with_suffix(".pdf")):
            figure.savefig(path, dpi=300, bbox_inches="tight",
                           facecolor=figure.get_facecolor())
            print(f"Written to {path}")
    return figure
