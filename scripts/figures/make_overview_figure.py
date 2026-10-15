#!/usr/bin/env python3
"""Make docs/project_overview.png and docs/project_overview.pdf.

One figure with four panels for Dravet: gene sets, propagation and
colocalization, enrichment and the knowledge graph. Panels 1, 2 (bottom) and 3
read the tables in results/micronutrient_coloc. Panel 2 (top) and
panel 4 are schematics.
All colours are colourblind-safe (Okabe-Ito, Tol muted, Blues).
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Wedge

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from micronutrient_coloc_figures import GRID, INK, MUTED, SURFACE

RESULTS = ROOT / "results" / "micronutrient_coloc"
OUT = ROOT / "docs" / "project_overview"

DRAVET = "#D55E00"
SET_GREY = "#999999"
SOURCES = {"GO:BP": "#882255", "GO:CC": "#44AA99", "REAC": "#DDCC77"}
SEED = "#08306B"
TOY_SETS = ["#CC6677", "#117733", "#AA4499"]

pairs = pd.read_csv(RESULTS / "tables" / "pair_stats.tsv", sep="\t")

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8, "text.color": INK,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.spines.top": False, "axes.spines.right": False,
})

fig = plt.figure(figsize=(15, 7), facecolor=SURFACE)
grid = fig.add_gridspec(2, 4, height_ratios=[0.8, 5.6], width_ratios=[0.55, 1, 1, 1.1],
                        wspace=0.35,
                        hspace=0.08, left=0.05, right=0.99, top=0.89, bottom=0.1)
fig.suptitle("Dravet syndrome and micronutrients: project overview",
             x=0.05, ha="left", fontsize=15, fontweight="bold", color=INK, y=0.975)

HEADERS = [
    ("1  Gene sets", "Dravet against 11 gene sets:\n"
     "10 micronutrients and the keto diet"),
    ("2  Propagation and colocalization", "Random walk with restart from each seed set,\n"
     "NetColoc Python package; 11 Dravet pairs"),
    ("3  Enrichment", "g:Profiler on the colocalized genes of each\n"
     "significant pair, collapsed by orsum"),
    ("4  Knowledge graph", ""),
]

# Panel titles and short descriptions
heads = []
for i, (title, text) in enumerate(HEADERS):
    head = fig.add_subplot(grid[0, i])
    heads.append(head)
    head.axis("off")
    head.text(0, 0.95, title, fontsize=11.5, fontweight="bold", color=INK,
              va="top", transform=head.transAxes)
    head.text(0, 0.45, text, fontsize=8.5, color=INK, va="top",
              transform=head.transAxes)

# Dravet pairs, sorted by colocalization z-score; this order is used in panels 1 to 3
rows = (pairs[pairs["disease"] == "Dravet"].sort_values("z_mean_z")
        .set_index("micronutrient"))
SET_ORDER = rows.index.tolist()

# Panel 1: seed genes of each set on the brain network
ax = fig.add_subplot(grid[1, 0])
counts = pd.concat([rows["n_micronutrient_seeds"],
                    pd.Series({"Dravet": rows["n_disease_seeds"].iloc[0]})])
ax.barh(counts.index, counts, height=0.7,
        color=[DRAVET if name == "Dravet" else SET_GREY for name in counts.index])
label = ax.get_yticklabels()[-1]
label.set_color(DRAVET)
label.set_fontweight("bold")
ax.tick_params(axis="y", length=0, labelsize=7.5)
ax.set_xlim(0, 210)
ax.set_xlabel("seed genes")
ax.set_title("brain network", fontsize=8.5, color=MUTED, loc="left")

# Panel 2, top: random walk with restart on a toy graph
sub = grid[1, 1].subgridspec(2, 1, height_ratios=[1, 1.1], hspace=0.3)
ax = fig.add_subplot(sub[0])
toy = nx.random_geometric_graph(42, 0.26, seed=4)
toy = toy.subgraph(max(nx.connected_components(toy), key=len)).copy()
order = list(toy.nodes)
A = nx.to_numpy_array(toy, nodelist=order)
W = A / A.sum(axis=0)
toy_seeds = [order[0], order[len(order) // 2]]
Y0 = np.isin(order, toy_seeds).astype(float)
score = 0.5 * np.linalg.solve(np.eye(len(order)) - 0.5 * W, Y0)
pos = nx.get_node_attributes(toy, "pos")
nx.draw_networkx_edges(toy, pos, ax=ax, edge_color=GRID, width=0.8)
nx.draw_networkx_nodes(toy, pos, ax=ax, node_size=70, node_color=np.log10(score),
                       cmap="Blues", edgecolors=MUTED, linewidths=0.4)
nx.draw_networkx_nodes(toy, pos, nodelist=toy_seeds, ax=ax, node_size=110,
                       node_color=SEED, edgecolors=INK, linewidths=1.4)
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=SEED, mec=INK,
                          label="seed gene"),
                   Line2D([], [], marker="o", ls="", color="#6baed6", mec=MUTED,
                          label="visit score (darker = more)")],
          frameon=False, loc="upper center", bbox_to_anchor=(0.5, 0.02),
          ncol=2, fontsize=7.5)
ax.axis("off")
ax.set_aspect("equal")

# Panel 2, bottom: colocalization z-score of each Dravet pair
size_max = float(pairs["coloc_size"].max())
area = lambda n: 15 + 160 * n / size_max
ax = fig.add_subplot(sub[1])
y = np.arange(len(rows))
z = rows["z_mean_z"].to_numpy()
sig = rows["significant"].to_numpy(bool)
size = area(rows["coloc_size"].to_numpy())
ax.hlines(y, 0, z, color=GRID, lw=1.4, zorder=1)
ax.scatter(z[sig], y[sig], s=size[sig], color=DRAVET, edgecolor=SURFACE,
           lw=0.6, zorder=3)
ax.scatter(z[~sig], y[~sig], s=size[~sig], facecolor=SURFACE,
           edgecolor=DRAVET, lw=1.1, zorder=3)
ax.axvline(0, color=GRID, lw=0.8)
ax.set_xlim(-12, 68)
ax.set_yticks(y, rows.index)
ax.tick_params(axis="y", length=0, labelsize=7, labelcolor=INK)
ax.set_xlabel("colocalization z-score")
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=MUTED,
                          markersize=np.sqrt(area(n)), label=f"{n} genes")
                   for n in [100, 500, 1000]]
          + [Line2D([], [], marker="o", ls="", color=MUTED, label="significant"),
             Line2D([], [], marker="o", ls="", mfc=SURFACE, mec=MUTED,
                    label="not significant")],
          frameon=False, loc="lower right", fontsize=7, labelspacing=0.9,
          title="module size", title_fontsize=7)

# Panel 3: the 4 orsum terms per source that represent the most terms
sub = grid[1, 2].subgridspec(2, 2, width_ratios=[0.75, 1], height_ratios=[4, 1],
                             wspace=0, hspace=0)
ax = fig.add_subplot(sub[0, 1])
# Dot area by term rank: rank 1 to 5, 6 to 20, more than 20
RANK_BINS = [(5, 90), (20, 45), (np.inf, 15)]
columns = [m for m in SET_ORDER[::-1] if m != "selenium"]
term_rows = []
for source, colour in SOURCES.items():
    table = pd.read_csv(RESULTS / "tables" / f"orsum_collapsed_dravet_{source.replace(':', '')}.tsv",
                        sep="\t")
    top = table.nlargest(4, "Represented term number")
    for _, term in top.iterrows():
        term_rows.append((term["Representing term name"], colour,
                          [term.get(f"{m} term rank", np.nan) for m in columns]))
for r, (name, colour, ranks) in enumerate(term_rows):
    for c, rank in enumerate(ranks):
        if pd.notna(rank):
            size = next(area for limit, area in RANK_BINS if rank <= limit)
            ax.scatter(c, r, s=size, color=colour, edgecolor=SURFACE, lw=0.5)
ax.set_xticks(range(len(columns)), columns, rotation=45, ha="right",
              rotation_mode="anchor")
ax.set_yticks(range(len(term_rows)), [name for name, _, _ in term_rows])
for label, (_, colour, _) in zip(ax.get_yticklabels(), term_rows):
    label.set_color(colour if colour != SOURCES["REAC"] else "#8a7a2a")
ax.set_xlim(-0.6, len(columns) - 0.4)
ax.set_ylim(len(term_rows) - 0.5, -0.5)
ax.tick_params(length=0, labelsize=8)
ax.tick_params(axis="y", labelsize=8.5)
ax.grid(color=GRID, lw=0.5)
ax.set_axisbelow(True)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=colour, label=source)
                   for source, colour in SOURCES.items()],
          frameon=False, loc="upper center", bbox_to_anchor=(0.3, -0.14),
          ncol=3, fontsize=6.5, handletextpad=0.2, columnspacing=0.8)

# Panel 4: knowledge graph schematic with Dravet at the centre
ax = fig.add_subplot(grid[1, 3])
# y is stretched to fill the tall panel
polar = lambda r, deg: (r * np.cos(np.radians(deg)), 1.4 * r * np.sin(np.radians(deg)))
node = lambda xy, text, fc, colour, box: ax.text(
    *xy, text, ha="center", va="center", fontsize=7.5, color=colour, zorder=3,
    bbox=dict(boxstyle=box, fc=fc, ec=SURFACE, lw=1))
LINKS = {"coloc": dict(color=DRAVET, lw=2), "module": dict(color=MUTED, lw=1),
         "enriched": dict(color=MUTED, lw=1, ls="--")}
edge = lambda a, b, kind: ax.plot(*zip(a, b), zorder=1, **LINKS[kind])
for deg, colour in zip((90, 210, 330), TOY_SETS):
    nutrient = polar(0.8, deg)
    edge((0, 0), nutrient, "coloc")
    angles = [deg - 30, deg + 30]
    genes = [polar(1.5, a) for a in angles]
    # Put the process and pathway on the upper gene, away from the legend
    k = int(genes[1][1] >= genes[0][1])
    for shift, term, fc, ink in ((-25, "process", SOURCES["GO:BP"], SURFACE),
                                 (25, "pathway", SOURCES["REAC"], INK)):
        end = np.add(genes[k], polar(0.7, angles[k] + shift))
        edge(genes[k], end, "enriched")
        node(end, term, fc, ink, "round,pad=0.3")
    for gene in genes:
        # The shared gene replaces the lower genes of the two bottom micronutrients
        if gene[1] < -1.5:
            continue
        edge(nutrient, gene, "module")
        node(gene, "gene", colour, SURFACE, "circle,pad=0.3")
    node(nutrient, "micro-\nnutrient", colour, SURFACE, "square,pad=0.4")
# Gene in the modules of two micronutrients, drawn as a pie
shared = polar(1.65, 270)
for deg in (210, 330):
    edge(polar(0.8, deg), shared, "module")
for (start, stop), colour in zip(((90, 270), (270, 90)), TOY_SETS[1:]):
    ax.add_patch(Wedge(shared, 0.23, start, stop, fc=colour, ec=SURFACE, lw=1, zorder=3))
ax.text(*shared, "gene", ha="center", va="center", fontsize=7.5, color=SURFACE, zorder=4)
node((0, 0), "Dravet", DRAVET, SURFACE, "round,pad=0.5")
ax.set_xlim(-2.5, 2.5)
ax.set_ylim(-3.5, 2.95)
ax.legend(handles=[Line2D([], [], label="significant colocalization", **LINKS["coloc"]),
                   Line2D([], [], label="gene in colocalized module", **LINKS["module"]),
                   Line2D([], [], label="gene in enriched term", **LINKS["enriched"])],
          frameon=False, loc="upper center", bbox_to_anchor=(0.5, 0.12), fontsize=8)
ax.set_aspect("equal")
ax.axis("off")

# Arrows between the panel titles
for a, b in zip(heads, heads[1:]):
    x0 = a.get_position().x1 + 0.006
    x1 = b.get_position().x0 - 0.012
    fig.patches.append(FancyArrowPatch(
        (x0, 0.883), (x1, 0.883), transform=fig.transFigure,
        arrowstyle="-|>", mutation_scale=12, color=MUTED, lw=1.2))

fig.text(0.05, 0.008, "CAKUT (kidney network) was used as validation of the workflow.",
         fontsize=8.5, color=MUTED, style="italic")

for ext in ("png", "pdf"):
    fig.savefig(OUT.with_suffix(f".{ext}"), dpi=300, facecolor=SURFACE)
print(f"wrote {OUT}.png and {OUT}.pdf")
