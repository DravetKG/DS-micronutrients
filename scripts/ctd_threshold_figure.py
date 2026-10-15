"""Figure of the per set reference threshold.

Panel A gives the set size against the reference threshold, with the chosen
threshold of each micronutrient marked. Panel B gives the size before and after
the filter.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from ctd_set_figure import (COLORS, GRID, GROUND, INK, INK_MUTED, ORDER,
                            style_axes)


def size_curve(counts, thresholds):
    """Count the genes that stay at each reference threshold."""
    counts = np.asarray(counts)
    return [int((counts >= t).sum()) for t in thresholds]


def panel_curves(axes, counts_by_set, chosen, order, min_size, max_size,
                 max_threshold):
    """Draw the set size against the reference threshold."""
    thresholds = np.arange(1, max_threshold + 1)
    axes.axhline(min_size, color=GRID, linewidth=0.9, zorder=1)
    axes.axhline(max_size, color=GRID, linewidth=0.9, zorder=1)
    axes.axhspan(min_size, max_size, color="#f2f1ec", linewidth=0, zorder=0)
    for name in order:
        sizes = size_curve(counts_by_set[name], thresholds)
        drawn = [index for index, size in enumerate(sizes) if size > 0]
        last = drawn[-1] + 1
        axes.plot(thresholds[:last], sizes[:last], color=COLORS[name],
                  linewidth=1.6, zorder=2)
        threshold = chosen[name]
        axes.scatter([threshold], [sizes[threshold - 1]],
                     s=52, color=COLORS[name], edgecolor=GROUND,
                     linewidth=1.2, zorder=3)
    axes.set_yscale("log")
    axes.set_xticks(thresholds)
    axes.set_xlabel("Reference threshold", fontsize=9, color=INK_MUTED)
    axes.set_ylabel("Genes", fontsize=9, color=INK_MUTED)
    axes.set_xlim(0.7, max_threshold + 0.3)
    axes.yaxis.grid(True, color=GRID, linewidth=0.6)
    axes.set_axisbelow(True)
    style_axes(axes)
    axes.text(max_threshold + 0.2, np.sqrt(min_size * max_size),
              f"{min_size}-{max_size}\ngenes", fontsize=7.5, color=INK_MUTED,
              ha="left", va="center")
    axes.set_title("A.  Set size against the reference threshold",
                   loc="left", fontsize=10, color=INK,
                   pad=10)


def panel_before_after(axes, counts_by_set, chosen, order, min_size, max_size):
    """Draw the set size before and after the filter."""
    positions = np.arange(len(order))
    before = [len(counts_by_set[name]) for name in order]
    after = [int((np.asarray(counts_by_set[name]) >= chosen[name]).sum())
             for name in order]

    axes.axhspan(min_size, max_size, color="#f2f1ec", linewidth=0, zorder=0)
    axes.bar(positions - 0.19, before, width=0.34, color=GRID, linewidth=0,
             zorder=2)
    axes.bar(positions + 0.19, after, width=0.34,
             color=[COLORS[name] for name in order], linewidth=0, zorder=2)
    for position, value in zip(positions - 0.19, before):
        axes.text(position, value * 1.06, str(value), ha="center", va="bottom",
                  fontsize=7.5, color=INK_MUTED)
    for position, value, name in zip(positions + 0.19, after, order):
        axes.text(position, value * 1.06, f"{value}\nn≥{chosen[name]}",
                  ha="center", va="bottom", fontsize=7.5, color=INK)
    axes.axhline(min_size, color=INK_MUTED, linewidth=0.7, linestyle=(0, (4, 3)),
                 zorder=3)
    axes.axhline(max_size, color=INK_MUTED, linewidth=0.7, linestyle=(0, (4, 3)),
                 zorder=3)
    axes.set_yscale("log")
    axes.set_xticks(positions)
    axes.set_xticklabels(order, rotation=30, ha="right")
    axes.set_ylabel("Genes", fontsize=9, color=INK_MUTED)
    axes.set_ylim(0.8, max(before) * 3.2)
    axes.yaxis.grid(True, color=GRID, linewidth=0.6)
    axes.set_axisbelow(True)
    style_axes(axes)
    axes.set_title("B.  Set size before and after the filter", loc="left",
                   fontsize=10, color=INK, pad=10)


def make_figure(counts_by_set, chosen, output_path=None, order=None,
                min_size=5, max_size=200, max_threshold=None):
    """Draw the two panels and write the figure.

    ``counts_by_set`` holds the reference counts of the genes of each set and
    ``chosen`` the threshold the filter selected.
    """
    order = order or ORDER
    if max_threshold is None:
        max_threshold = max(max(chosen.values()) + 2, 8)
    figure = plt.figure(figsize=(12, 8.4), facecolor=GROUND)
    grid = figure.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.42)
    axes_curves = figure.add_subplot(grid[0])
    axes_bars = figure.add_subplot(grid[1])
    for axes in (axes_curves, axes_bars):
        axes.set_facecolor(GROUND)

    panel_curves(axes_curves, counts_by_set, chosen, order, min_size, max_size,
                 max_threshold)
    panel_before_after(axes_bars, counts_by_set, chosen, order, min_size,
                       max_size)

    handles = [Line2D([0], [0], marker="s", linestyle="", markersize=7,
                      color=COLORS[name], label=name) for name in order]
    handles.append(Line2D([0], [0], marker="s", linestyle="", markersize=7,
                          color=GRID, label="Before the filter"))
    figure.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
                  fontsize=8, labelcolor=INK_MUTED, bbox_to_anchor=(0.5, -0.005))

    figure.suptitle("CTD micronutrient gene sets, reference filter",
                    fontsize=11.5, color=INK, x=0.02, ha="left", y=0.98)
    if output_path:
        for path in (Path(output_path), Path(output_path).with_suffix(".pdf")):
            figure.savefig(path, dpi=300, bbox_inches="tight",
                           facecolor=figure.get_facecolor())
            print(f"Written to {path}")
    return figure
