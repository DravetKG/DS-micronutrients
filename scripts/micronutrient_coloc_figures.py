"""Figures of the disease against micronutrient colocalization.

``pair_dotplot`` grids the pair statistics, ``gene_heatmap`` the genes of the
colocalized modules, and ``enrichment_dotplot`` the enriched terms.
"""

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
YELLOW, MAGENTA = "#eda100", "#e87ba4"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
SURFACE = "#fcfcfb"

def gene_heatmap(matrix, out_png=None, max_genes=60, title=None):
    """Genes x micronutrient grid of the combined z of the colocalized genes.

    An empty cell means the gene is not in that module, not a combined z of zero.
    """
    if not len(matrix):
        print("no colocalized gene")
        return None
    matrix = matrix.head(max_genes)

    fig, ax = plt.subplots(figsize=(1.15 * len(matrix.columns) + 3.4,
                                    0.21 * len(matrix) + 2.4),
                           facecolor=SURFACE)
    image = ax.imshow(matrix.values, cmap="YlGnBu", aspect="auto")
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, rotation=30, ha="right", fontsize=8.5,
                       color=MUTED)
    ax.set_yticks(range(len(matrix)))
    ax.set_yticklabels(matrix.index, fontsize=6.5, color=MUTED)
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.set_title(title or "Colocalized genes",
                 loc="left", fontsize=10.5, color=INK, pad=10)

    bar = fig.colorbar(image, ax=ax, fraction=0.03, pad=0.03)
    bar.set_label("combined z", fontsize=8, color=MUTED)
    bar.ax.tick_params(labelsize=7.5, colors=MUTED, length=2)
    bar.outline.set_visible(False)
    fig.tight_layout()
    if out_png:
        for path in (Path(out_png), Path(out_png).with_suffix(".pdf")):
            fig.savefig(path, dpi=300, bbox_inches="tight", facecolor=SURFACE)
    return fig


def annotation_table(annotations, labels, out_tsv=None):
    """One long table of the g:Profiler terms of every significant pair."""
    frames = []
    for (disease, nutrient), ann in annotations.items():
        if ann is None or not len(ann):
            continue
        frame = ann.copy()
        frame.insert(0, "disease", labels[disease])
        frame.insert(1, "micronutrient", labels[nutrient])
        frames.append(frame)
    if not frames:
        return pd.DataFrame()

    table = pd.concat(frames, ignore_index=True)
    keep = [c for c in ["disease", "micronutrient", "source", "native", "name",
                        "p_value", "term_size", "query_size",
                        "intersection_size", "intersections"]
            if c in table.columns]
    table = table[keep].sort_values(["disease", "micronutrient", "p_value"])
    if out_tsv:
        table.to_csv(out_tsv, sep="\t", index=False)
    return table.reset_index(drop=True)


def _dot_area(hits, size_max):
    """Dot area of a term, on a square root scale of the gene count."""
    return 22 + 300 * np.sqrt(np.asarray(hits, dtype=float) / size_max)


def _shorten(name, width=46):
    """Cut a term name to ``width`` characters at a word boundary."""
    if len(name) <= width:
        return name
    cut = name[:width].rsplit(" ", 1)[0]
    return (cut or name[:width]) + "\u2026"


def _wrap_label(name, width=16):
    """Wrap a column label on two lines to keep the columns close."""
    return "\n".join(textwrap.wrap(name, width)[:2])


def enrichment_dotplot(table, top_terms=8, out_png=None, title=None):
    """Terms x micronutrient dot plot of the enrichment of one disease.

    Colour is -log10 p and dot area is the query gene count. The size scale is per
    figure, so do not compare dot areas between figures.
    """
    if not len(table):
        print("no term to draw")
        return None

    top = (table.sort_values("p_value")
           .groupby("micronutrient", group_keys=False).head(top_terms))
    terms = list(dict.fromkeys(top.sort_values("p_value")["name"]))
    nutrients = sorted(top["micronutrient"].unique())
    subset = table[table["name"].isin(terms)]

    index = {name: i for i, name in enumerate(terms)}
    column = {name: j for j, name in enumerate(nutrients)}
    x = [column[n] for n in subset["micronutrient"]]
    y = [index[t] for t in subset["name"]]
    colour = -np.log10(subset["p_value"].clip(lower=1e-300))
    hits = subset["intersection_size"].astype(float)
    size_max = float(hits.max()) or 1.0
    size = _dot_area(hits, size_max)

    fig, ax = plt.subplots(figsize=(0.62 * len(nutrients) + 5.2,
                                    0.30 * len(terms) + 1.9),
                           facecolor=SURFACE)
    dots = ax.scatter(x, y, c=colour, s=size, cmap="viridis",
                      edgecolors=SURFACE, linewidths=0.6)
    ax.set_xticks(range(len(nutrients)))
    ax.set_xticklabels([_wrap_label(n) for n in nutrients], rotation=45,
                       ha="right", va="top", fontsize=8, color=INK)
    ax.set_yticks(range(len(terms)))
    ax.set_yticklabels([_shorten(t) for t in terms], fontsize=8, color=INK)
    ax.set_ylim(len(terms) - 0.5, -0.5)
    ax.set_xlim(-0.5, len(nutrients) - 0.5)
    ax.set_axisbelow(True)
    ax.grid(True, color=GRID, linewidth=0.5, alpha=0.7)
    for edge in np.arange(0.5, len(terms) - 0.5):
        ax.axhline(edge, color=GRID, linewidth=0.4, alpha=0.45)
    ax.tick_params(length=0, pad=2)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.set_title(title or "Enriched terms of the significant pairs",
                 loc="left", fontsize=10.5, color=INK, pad=10)

    bar = fig.colorbar(dots, ax=ax, fraction=0.03, pad=0.02)
    bar.set_label("-log10 p", fontsize=8, color=MUTED)
    bar.ax.tick_params(labelsize=7.5, colors=MUTED, length=2)
    bar.outline.set_visible(False)

    counts = sorted({int(v) for v in np.percentile(hits, [5, 50, 95])})
    ax.legend(handles=[Line2D([], [], marker="o", linestyle="", color=MUTED,
                              markersize=np.sqrt(_dot_area(c, size_max)),
                              label=f"{c} genes") for c in counts],
              frameon=False, fontsize=8, loc="upper left",
              bbox_to_anchor=(1.02, 1.0), labelcolor=MUTED,
              borderpad=0.8, labelspacing=1.4, handletextpad=1.0,
              title="genes in term", title_fontsize=8)
    fig.tight_layout()
    if out_png:
        for path in (Path(out_png), Path(out_png).with_suffix(".pdf")):
            fig.savefig(path, dpi=300, bbox_inches="tight", facecolor=SURFACE)
    return fig


SET_COLORS = [BLUE, ORANGE, AQUA, MAGENTA, YELLOW, "#7b5ad1", "#00a3b5"]


def _label_points(ax, positions, fontsize=6.5, step=0.16, rounds=60):
    """Write one label per point, pushed apart so the labels do not overlap."""
    genes = list(positions)
    points = np.array([positions[g] for g in genes], dtype=float)
    labels = points + np.array([0.0, step])
    for _ in range(rounds):
        moved = False
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                delta = labels[i] - labels[j]
                distance = np.linalg.norm(delta)
                if distance < step * 1.6:
                    push = (delta / max(distance, 1e-6)) * step * 0.5
                    labels[i] += push
                    labels[j] -= push
                    moved = True
        if not moved:
            break
    for gene, point, label in zip(genes, points, labels):
        ax.annotate(gene, point, xytext=label, fontsize=fontsize, color=INK,
                    ha="center", va="center", zorder=6,
                    arrowprops=dict(arrowstyle="-", color=GRID, linewidth=0.5,
                                    shrinkA=0, shrinkB=2))


def pair_dotplot(pair_stats, out_png=None, q_thresh=0.05,
                 title="Colocalization of each disease with each gene set"):
    """One ranked panel per disease of the pair statistic.

    A filled dot is BH-significant at ``q_thresh`` and dot area is the colocalized
    module size. The margin text is the seed count, a power proxy, not a result.
    """
    diseases = list(dict.fromkeys(pair_stats["disease"]))
    heights = [pair_stats["disease"].eq(d).sum() for d in diseases]
    fig, axes = plt.subplots(
        len(diseases), 1, sharex=True, squeeze=False,
        figsize=(8.6, 0.34 * sum(heights) + 1.6 * len(diseases)),
        gridspec_kw={"height_ratios": heights}, facecolor=SURFACE)
    axes = axes[:, 0]

    z_all = pair_stats["z_mean_z"].astype(float)
    left = min(0.0, float(z_all.min())) - 0.6
    right = float(z_all.max()) * 1.06 + 0.6
    size_max = float(pair_stats["coloc_size"].max()) or 1.0

    for ax, disease in zip(axes, diseases):
        rows = (pair_stats[pair_stats["disease"] == disease]
                .sort_values("z_mean_z", ascending=False).reset_index(drop=True))
        y = np.arange(len(rows))[::-1]
        z = rows["z_mean_z"].astype(float).values
        sig = rows["significant"].astype(bool).values
        area = 40 + 260 * (rows["coloc_size"].astype(float).values / size_max)

        ax.hlines(y, 0, z, color=GRID, linewidth=1.4, zorder=1)
        ax.scatter(z[sig], y[sig], s=area[sig], color=BLUE,
                   edgecolors=SURFACE, linewidths=0.8, zorder=3)
        ax.scatter(z[~sig], y[~sig], s=area[~sig], facecolors=SURFACE,
                   edgecolors=MUTED, linewidths=1.1, zorder=3)
        ax.axvline(0, color=GRID, linewidth=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{n}  ({s} seeds)" for n, s in
                            zip(rows["micronutrient"],
                                rows["n_micronutrient_seeds"])],
                           fontsize=8.5, color=MUTED)
        ax.set_ylim(-0.7, len(rows) - 0.3)
        ax.set_xlim(left, right)
        ax.grid(True, axis="x", color=GRID, linewidth=0.5, alpha=0.7)
        ax.set_axisbelow(True)
        ax.tick_params(length=0, labelsize=8, colors=MUTED)
        for side in ax.spines.values():
            side.set_visible(False)
        ax.set_title(disease, loc="left", fontsize=9.5, color=INK, pad=6)

    axes[-1].set_xlabel("network colocalization z-score (vs. permuted gene labels)", fontsize=9, color=MUTED)
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, color=INK)

    counts = sorted({int(v) for v in np.percentile(
        pair_stats["coloc_size"].astype(float), [10, 50, 90])})
    size_handles = [Line2D([], [], marker="o", linestyle="", color=BLUE,
                           markersize=np.sqrt(40 + 260 * (c / size_max)),
                           label=f"{c} genes") for c in counts]
    fill_handles = [Line2D([], [], marker="o", linestyle="", color=BLUE,
                           markersize=7, label=f"BH q < {q_thresh}"),
                    Line2D([], [], marker="o", linestyle="", color=SURFACE,
                           markeredgecolor=MUTED, markersize=7,
                           label="not significant")]
    size_legend = axes[0].legend(handles=size_handles, frameon=False,
                                 fontsize=8, loc="upper left",
                                 bbox_to_anchor=(1.02, 1.0), labelcolor=MUTED,
                                 title="module size", title_fontsize=8,
                                 labelspacing=1.5, borderpad=0.8,
                                 handletextpad=1.4)
    axes[0].add_artist(size_legend)
    axes[-1].legend(handles=fill_handles, frameon=False, fontsize=8,
                    loc="lower left", bbox_to_anchor=(1.02, 0.0),
                    labelcolor=MUTED)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    if out_png:
        for path in (Path(out_png), Path(out_png).with_suffix(".pdf")):
            fig.savefig(path, dpi=300, bbox_inches="tight", facecolor=SURFACE)
    return fig
