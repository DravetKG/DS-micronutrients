"""Knowledge graph of the micronutrient colocalization of one disease.

Nodes are micronutrients, the genes of the colocalized modules, and enriched
terms. A ``coloc`` edge carries the combined z and the gene role.
"""

from pathlib import Path

import networkx as nx
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
TABLES = ROOT / "results" / "micronutrient_coloc" / "tables"

TOP_TERMS = 5


def build(edges, annotation=None, top_terms=TOP_TERMS, min_nutrients=1):
    """Build the graph from the colocalization edges and the enriched terms.

    ``min_nutrients`` drops the genes that only one micronutrient reaches.
    """
    per_gene = edges.groupby("gene")["micronutrient"].nunique()
    kept = edges[edges["gene"].isin(per_gene[per_gene >= min_nutrients].index)]

    graph = nx.Graph(disease=str(edges["disease"].iloc[0]))
    for nutrient, rows in kept.groupby("micronutrient"):
        graph.add_node(nutrient, node_type="micronutrient", label=nutrient,
                       n_genes=int(rows["gene"].nunique()))

    for gene, rows in kept.groupby("gene"):
        roles = sorted(set(rows["role"]))
        best = rows.loc[rows["combined_z"].idxmax()]
        graph.add_node(gene, node_type="gene", label=gene,
                       entrez=int(rows["entrez"].iloc[0]),
                       n_micronutrients=int(rows["micronutrient"].nunique()),
                       micronutrients=";".join(sorted(set(rows["micronutrient"]))),
                       top_micronutrient=str(best["micronutrient"]),
                       max_combined_z=float(rows["combined_z"].max()),
                       role=";".join(roles),
                       is_seed=any(r.endswith("seed") for r in roles),
                       is_disease_seed="disease seed" in roles or "both seeds" in roles)

    for row in kept.itertuples():
        graph.add_edge(row.micronutrient, row.gene, edge_type="coloc",
                       combined_z=float(row.combined_z), role=row.role)

    if annotation is not None and len(annotation):
        _add_pathways(graph, annotation, kept, top_terms)
    return graph


def _add_pathways(graph, annotation, edges, top_terms):
    """Join the top terms of every micronutrient to the graph."""
    nutrients = set(edges["micronutrient"])
    top = (annotation[annotation["micronutrient"].isin(nutrients)]
           .sort_values("p_value")
           .groupby("micronutrient", group_keys=False)
           .head(top_terms))

    for row in top.itertuples():
        if not graph.has_node(row.native):
            graph.add_node(row.native, node_type="pathway", label=row.name,
                           source=row.source, term_size=int(row.term_size))
        graph.add_edge(row.micronutrient, row.native, edge_type="enriched",
                       p_value=float(row.p_value),
                       intersection_size=int(row.intersection_size))


def summary(graph):
    """Count the nodes of every type and the edges of every type."""
    nodes = pd.Series([d["node_type"] for _, d in graph.nodes(data=True)])
    edges = pd.Series([d["edge_type"] for _, _, d in graph.edges(data=True)])
    return nodes.value_counts(), edges.value_counts()


DISEASE_SEED_COLOR = "#e6002e"
SEED_RING_COLOR = "#0b0b0b"
BRIDGE_COLOR = "#b8b7b3"
PATHWAY_COLOR = "#fcfcfb"
GENE_RING_WIDTH = 0.6
SEED_RING_WIDTH = 1.6
DISEASE_SEED_RING_WIDTH = 1.6
EXTRA_SET_COLORS = ["#8c6239", "#5c8a2f", "#c2185b", "#0b7285", "#a06cd5"]


def nutrient_colors(graph):
    """Map every micronutrient of the graph to one palette colour.

    The drawing and the Cytoscape export share this map, so the two agree.
    ``EXTRA_SET_COLORS`` carries the palette past the sets of ``SET_COLORS``,
    so two micronutrients never take the same colour.
    """
    import micronutrient_coloc_figures as mcf

    nutrients = sorted(n for n, k in
                       nx.get_node_attributes(graph, "node_type").items()
                       if k == "micronutrient")
    palette = [c for c in list(mcf.SET_COLORS) + EXTRA_SET_COLORS
               if c != DISEASE_SEED_COLOR]
    if len(nutrients) > len(palette):
        raise ValueError(f"{len(nutrients)} micronutrients but only "
                         f"{len(palette)} palette colours")
    return {n: palette[i] for i, n in enumerate(nutrients)}


def export(graph, prefix, out_dir=None):
    """Write the graph as GraphML, as node and edge tables, and as a style.

    The style columns are baked in, so Cytoscape draws it in the colours of ``draw``.
    """
    out = Path(out_dir or TABLES.parent / "kg")
    out.mkdir(parents=True, exist_ok=True)

    styled = with_cytoscape_style(graph)
    write_graphml(styled, out / f"{prefix}.graphml")
    nodes = pd.DataFrame([{"node": n, **d} for n, d in styled.nodes(data=True)])
    links = pd.DataFrame([{"source": u, "target": v, **d}
                          for u, v, d in styled.edges(data=True)])
    nodes.to_csv(out / f"{prefix}_nodes.tsv", sep="\t", index=False)
    links.to_csv(out / f"{prefix}_edges.tsv", sep="\t", index=False)
    write_cytoscape_style(out / f"{prefix}_style.xml", f"{prefix} style")
    return out / f"{prefix}.graphml"


NODE_SIZE = {"micronutrient": 55.0, "gene": 22.0, "pathway": 40.0}
NODE_SHAPE = {"micronutrient": "RECTANGLE", "gene": "ELLIPSE",
              "pathway": "ROUND_RECTANGLE"}


def with_cytoscape_style(graph):
    """Copy the graph with the drawing colours baked in as node columns.

    A gene several micronutrients reach carries a ``piechart`` string for
    enhancedGraphics; ``fill_color`` holds the colour of its strongest pair.
    """
    import micronutrient_coloc_figures as mcf

    styled = nx.Graph()
    colour_of = nutrient_colors(graph)

    for node, data in graph.nodes(data=True):
        kind = data["node_type"]
        fill = BRIDGE_COLOR
        border, width = mcf.SURFACE, 1.0
        chart = ""
        if kind == "micronutrient":
            fill = colour_of[node]
        elif kind == "pathway":
            fill, border = PATHWAY_COLOR, mcf.MUTED
        else:
            members = sorted(n for n in graph.neighbors(node)
                             if graph.nodes[n]["node_type"] == "micronutrient")
            top = data.get("top_micronutrient")
            if top in colour_of:
                fill = colour_of[top]
            elif members:
                fill = colour_of[members[0]]
            border = SEED_RING_COLOR
            if _is_disease_seed(graph, node):
                fill, width = DISEASE_SEED_COLOR, 4.0
            else:
                width = 3.0 if _is_seed(graph, node) else 1.5
                if len(members) > 1:
                    chart = _piechart(members, colour_of)

        size = NODE_SIZE[kind]
        if kind == "gene":
            size *= 1.0 + 0.30 * (int(data.get("n_micronutrients", 1)) - 1)
            if _is_disease_seed(graph, node):
                size *= 1.5
        styled.add_node(node, **data, fill_color=fill, border_color=border,
                        border_width=float(width), shape=NODE_SHAPE[kind],
                        size=round(float(size), 1), chart=chart)

    for u, v, data in graph.edges(data=True):
        colour = (colour_of.get(u) or colour_of.get(v) or mcf.GRID)
        styled.add_edge(u, v, **data, line_color=colour,
                        line_width=1.0 if data["edge_type"] == "coloc" else 2.0)
    return styled


def _piechart(members, colour_of):
    """Write the enhancedGraphics string that draws one gene as a pie.

    The slices are equal: z values of two pairs come from different nulls, so
    their ratio is not a quantity to size by.
    """
    values = ",".join(["1.0"] * len(members))
    colours = ",".join(colour_of[m] for m in members)
    labels = ",".join(m.split(" (")[0] for m in members)
    return (f'piechart: valuelist="{values}" colorlist="{colours}" '
            f'labellist="{labels}" showlabels="false" '
            f'arcstart="90" arcdirection="clockwise" borderwidth="0"')


GRAPHML_TYPE = {bool: "boolean", int: "int", float: "double", str: "string"}


def write_graphml(graph, path):
    """Write GraphML with the standard library only.

    ``networkx.write_graphml`` of this release reads ``numpy.float_``, which NumPy
    2 removed, so every value is cast to a plain Python type here.
    """
    import xml.etree.ElementTree as ET

    def cast(value):
        if isinstance(value, bool):
            return bool(value), bool
        if isinstance(value, (int,)) and not isinstance(value, bool):
            return int(value), int
        if isinstance(value, float):
            return float(value), float
        try:
            import numpy as np
            if isinstance(value, np.bool_):
                return bool(value), bool
            if isinstance(value, np.integer):
                return int(value), int
            if isinstance(value, np.floating):
                return float(value), float
        except ImportError:
            pass
        return str(value), str

    ns = "http://graphml.graphdrawing.org/xmlns"
    root = ET.Element("graphml", {"xmlns": ns})
    keys, order = {}, []
    for scope, items in (("node", (d for _, d in graph.nodes(data=True))),
                         ("edge", (d for _, _, d in graph.edges(data=True)))):
        for data in items:
            for name, value in data.items():
                if (scope, name) in keys:
                    continue
                key_id = f"{scope[0]}_{name}"
                keys[(scope, name)] = key_id
                order.append((key_id, scope, name, GRAPHML_TYPE[cast(value)[1]]))

    for key_id, scope, name, kind in order:
        ET.SubElement(root, "key", {"id": key_id, "for": scope,
                                    "attr.name": name, "attr.type": kind})

    element = ET.SubElement(root, "graph", {"id": "G", "edgedefault": "undirected"})
    for node, data in graph.nodes(data=True):
        item = ET.SubElement(element, "node", {"id": str(node)})
        _graphml_data(ET, item, data, "node", keys, cast)
    for i, (u, v, data) in enumerate(graph.edges(data=True)):
        item = ET.SubElement(element, "edge",
                             {"id": f"e{i}", "source": str(u), "target": str(v)})
        _graphml_data(ET, item, data, "edge", keys, cast)

    ET.indent(root, space="  ")
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)
    return path


def _graphml_data(ET, item, data, scope, keys, cast):
    """Write one ``data`` element per attribute of a node or an edge."""
    for name, value in data.items():
        clean, kind = cast(value)
        child = ET.SubElement(item, "data", {"key": keys[(scope, name)]})
        child.text = str(clean).lower() if kind is bool else str(clean)


def write_cytoscape_style(path, name="coloc kg"):
    """Write a Cytoscape visual style that passes the baked columns through.

    The pie of a shared gene needs enhancedGraphics; without it the node keeps
    ``fill_color``.
    """
    mappings = [
        ("NODE_FILL_COLOR", "fill_color", "string"),
        ("NODE_BORDER_PAINT", "border_color", "string"),
        ("NODE_BORDER_WIDTH", "border_width", "float"),
        ("NODE_SHAPE", "shape", "string"),
        ("NODE_SIZE", "size", "float"),
        ("NODE_LABEL", "label", "string"),
        ("NODE_CUSTOMGRAPHICS_1", "chart", "string"),
        ("EDGE_STROKE_UNSELECTED_PAINT", "line_color", "string"),
        ("EDGE_WIDTH", "line_width", "float"),
    ]
    lines = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
             '<vizmap documentVersion="3.0" id="coloc_kg">',
             f'  <visualStyle name="{name}">',
             '    <network>',
             '      <visualProperty name="NETWORK_BACKGROUND_PAINT" '
             'default="#FCFCFB"/>', '    </network>', '    <node>',
             '      <visualProperty name="NODE_TRANSPARENCY" default="255"/>',
             '      <visualProperty name="NODE_LABEL_FONT_SIZE" default="9"/>',
             '      <visualProperty name="NODE_LABEL_COLOR" default="#0B0B0B"/>',
             '      <visualProperty name="NODE_CUSTOMGRAPHICS_SIZE_1" '
             'default="1.0"/>',
             '      <visualProperty name="NODE_CUSTOMGRAPHICS_POSITION_1" '
             'default="C,C,c,0.00,0.00"/>']
    for prop, column, kind in mappings:
        scope = "edge" if prop.startswith("EDGE") else "node"
        if scope == "node":
            lines.append(f'      <visualProperty name="{prop}">')
            lines.append(f'        <passthroughMapping attributeName="{column}" '
                         f'attributeType="{kind}"/>')
            lines.append('      </visualProperty>')
    lines += ['    </node>', '    <edge>',
              '      <visualProperty name="EDGE_TRANSPARENCY" default="120"/>']
    for prop, column, kind in mappings:
        if prop.startswith("EDGE"):
            lines.append(f'      <visualProperty name="{prop}">')
            lines.append(f'        <passthroughMapping attributeName="{column}" '
                         f'attributeType="{kind}"/>')
            lines.append('      </visualProperty>')
    lines += ['    </edge>', '  </visualStyle>', '</vizmap>', '']
    Path(path).write_text("\n".join(lines), encoding="utf-8")
    return path




def draw(graph, out_png=None, seed=0, label_min=3, with_pathways=True,
         title=None):
    """Draw the graph with the micronutrients as hubs and the genes around.

    A gene takes the colour of the micronutrient it colocalizes with, or a pie
    when several reach it. A disease seed is red, because it is an input.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.lines import Line2D
    from matplotlib.patches import Wedge
    import micronutrient_coloc_figures as mcf

    kinds = nx.get_node_attributes(graph, "node_type")
    nutrients = sorted(n for n, k in kinds.items() if k == "micronutrient")
    genes = sorted(n for n, k in kinds.items() if k == "gene")
    pathways = sorted(n for n, k in kinds.items() if k == "pathway")
    if not with_pathways and pathways:
        graph = graph.subgraph(nutrients + genes)
        pathways = []

    colour_of = nutrient_colors(graph)
    sets_of = {g: sorted(n for n in graph.neighbors(g) if n in nutrients)
               for g in genes}
    layout = _spring_layout(graph, nutrients, seed=seed)

    fig, ax = plt.subplots(figsize=(13.5, 11.0), facecolor=mcf.SURFACE)
    nx.draw_networkx_edges(graph, layout, ax=ax, width=0.6, alpha=0.5,
                           edge_color=mcf.GRID)

    radius = _node_radius(layout)
    for gene in genes:
        _draw_gene(ax, gene, graph, layout, sets_of, colour_of, radius, mcf,
                   Wedge)
    _draw_nutrients(ax, nutrients, graph, layout, colour_of, radius, mcf)
    if pathways:
        _draw_pathway_boxes(ax, pathways, graph, layout, mcf)

    labelled = [g for g in genes
                if len(sets_of[g]) >= label_min or graph.nodes[g]["is_seed"]]
    if labelled:
        mcf._label_points(ax, {g: np.asarray(layout[g]) for g in labelled},
                          fontsize=6.5, step=radius * 2.6)

    _draw_legend(ax, nutrients, colour_of, bool(pathways), mcf, Line2D,
                 graph.graph["disease"])
    ax.set_title(title or f"{graph.graph['disease']} micronutrient colocalization",
                 loc="left", fontsize=12, color=mcf.INK)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.margins(0.10)
    fig.tight_layout()
    if out_png:
        for path in (Path(out_png), Path(out_png).with_suffix(".pdf")):
            fig.savefig(path, dpi=300, bbox_inches="tight", facecolor=mcf.SURFACE)
    return fig


def shared_genes_table(graph):
    """Rank the genes by the number of micronutrients that reach them.

    The columns are the gene, the micronutrient count, the micronutrients, the
    largest combined z over the pairs, and the seed flags.
    """
    kinds = nx.get_node_attributes(graph, "node_type")
    nutrients = {n for n, k in kinds.items() if k == "micronutrient"}
    rows = []
    for gene, data in graph.nodes(data=True):
        if data["node_type"] != "gene":
            continue
        sets = sorted(n for n in graph.neighbors(gene) if n in nutrients)
        rows.append({"gene": gene, "n_micronutrients": len(sets),
                     "micronutrients": ";".join(sets),
                     "max_combined_z": data["max_combined_z"],
                     "is_seed": data["is_seed"],
                     "is_disease_seed": data["is_disease_seed"]})
    table = pd.DataFrame(rows)
    return table.sort_values(["n_micronutrients", "max_combined_z"],
                             ascending=False).reset_index(drop=True)


def draw_shared_genes(graph, out_png=None, top=30, min_nutrients=2,
                      title=None):
    """Bar the genes that most micronutrients reach.

    One bar per gene, split into one segment per micronutrient, in the colours
    of ``draw``. A red gene name is a disease seed.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import micronutrient_coloc_figures as mcf

    table = shared_genes_table(graph)
    table = table[table["n_micronutrients"] >= min_nutrients].head(top)
    if not len(table):
        print(f"no gene is reached by {min_nutrients} micronutrients or more")
        return None

    colour_of = nutrient_colors(graph)
    order = sorted(colour_of, key=lambda n: -graph.nodes[n]["n_genes"])
    table = table.iloc[::-1]

    fig, ax = plt.subplots(figsize=(8.4, 0.26 * len(table) + 1.6),
                           facecolor=mcf.SURFACE)
    left = pd.Series(0.0, index=table.index)
    for nutrient in order:
        width = table["micronutrients"].str.split(";").apply(
            lambda s: float(nutrient in s))
        ax.barh(range(len(table)), width, left=left, height=0.72,
                color=colour_of[nutrient], edgecolor=mcf.SURFACE,
                linewidth=0.6, label=nutrient.split(" (")[0])
        left = left + width

    ax.set_yticks(range(len(table)))
    ax.set_yticklabels(table["gene"], fontsize=7.5)
    for tick, seed in zip(ax.get_yticklabels(), table["is_disease_seed"]):
        tick.set_color(DISEASE_SEED_COLOR if seed else mcf.MUTED)
    ax.set_xlabel("micronutrient gene sets reaching the gene", fontsize=8.5,
                  color=mcf.MUTED)
    ax.set_xticks(range(0, int(table["n_micronutrients"].max()) + 1))
    ax.tick_params(length=0, labelsize=8, colors=mcf.MUTED)
    ax.set_ylim(-0.7, len(table) - 0.3)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, color=mcf.GRID, linewidth=0.6)

    handles, _ = ax.get_legend_handles_labels()
    handles.append(Line2D([], [], marker="s", linestyle="", markersize=0,
                          color=DISEASE_SEED_COLOR,
                          label=f"red name: {graph.graph['disease']} seed"))
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper left",
              bbox_to_anchor=(1.02, 1.0), labelcolor=mcf.MUTED)
    ax.set_title(title or f"{graph.graph['disease']}: the most shared colocalized "
                          "genes", loc="left", fontsize=10.5, color=mcf.INK,
                 pad=10)
    fig.tight_layout()
    if out_png:
        for path in (Path(out_png), Path(out_png).with_suffix(".pdf")):
            fig.savefig(path, dpi=300, bbox_inches="tight", facecolor=mcf.SURFACE)
    return fig


def _node_radius(layout):
    """Scale the nodes to the extent of the layout."""
    import numpy as np

    points = np.array(list(layout.values()))
    span = float(np.ptp(points, axis=0).max()) or 1.0
    return span * 0.0075


def _draw_gene(ax, gene, graph, layout, sets_of, colour_of, radius, mcf, Wedge):
    """Draw one gene as a disc, or as a pie if several micronutrients reach it.

    The ring says what the gene is: thin black for a bridge, thick black for a
    micronutrient seed, red and larger for a disease seed.
    """
    import numpy as np

    centre = np.asarray(layout[gene], dtype=float)
    members = sets_of[gene]
    disease_seed = _is_disease_seed(graph, gene)
    ring = _ring_width(graph, gene)
    ring_colour = DISEASE_SEED_COLOR if disease_seed else SEED_RING_COLOR
    size = radius * (1.0 + 0.32 * (len(members) - 1))
    if disease_seed:
        size *= 1.5
    zorder = 5 if disease_seed else 3

    if disease_seed:
        ax.add_patch(plt_circle(centre, size, DISEASE_SEED_COLOR,
                                SEED_RING_COLOR, zorder=zorder, linewidth=ring))
        return
    if not members:
        ax.add_patch(plt_circle(centre, radius, BRIDGE_COLOR, ring_colour,
                                zorder=2, linewidth=ring))
        return
    if len(members) == 1:
        ax.add_patch(plt_circle(centre, size, colour_of[members[0]],
                                ring_colour, zorder=zorder, linewidth=ring))
        return

    step = 360.0 / len(members)
    for i, nutrient in enumerate(members):
        ax.add_patch(Wedge(centre, size, i * step, (i + 1) * step,
                           facecolor=colour_of[nutrient], edgecolor=mcf.SURFACE,
                           linewidth=0.4, zorder=zorder))
    ax.add_patch(plt_circle(centre, size, "none", ring_colour,
                            zorder=zorder + 1, linewidth=ring))


def plt_circle(centre, radius, facecolor, edgecolor, zorder=3, linewidth=0.5):
    """One circle patch, which the pies and the plain nodes share."""
    from matplotlib.patches import Circle

    return Circle(centre, radius, facecolor=facecolor, edgecolor=edgecolor,
                  linewidth=linewidth, zorder=zorder)


def _is_disease_seed(graph, gene):
    """Say whether the gene is a seed of the disease."""
    return bool(graph.nodes[gene].get("is_disease_seed",
                                      "disease seed" in graph.nodes[gene]["role"]))


def _is_seed(graph, gene):
    """Say whether the gene is a seed of the disease or of a micronutrient."""
    return bool(graph.nodes[gene].get("is_seed",
                                      "seed" in graph.nodes[gene]["role"]))


def _ring_width(graph, gene):
    """Give the ring width that says what the gene is to the colocalization."""
    if _is_disease_seed(graph, gene):
        return DISEASE_SEED_RING_WIDTH
    return SEED_RING_WIDTH if _is_seed(graph, gene) else GENE_RING_WIDTH


def _draw_nutrients(ax, nutrients, graph, layout, colour_of, radius, mcf):
    """Draw every micronutrient as a square in its own colour."""
    import numpy as np
    from matplotlib.patches import Rectangle

    for nutrient in nutrients:
        x, y = layout[nutrient]
        side = radius * 5.2
        ax.add_patch(Rectangle((x - side / 2, y - side / 2), side, side,
                               facecolor=colour_of[nutrient],
                               edgecolor=SEED_RING_COLOR, linewidth=1.2,
                               zorder=10))
        ax.text(x, y - side * 0.85, nutrient.split(" (")[0],
                fontsize=9.5, fontweight="bold", color=colour_of[nutrient],
                ha="center", va="top", zorder=11,
                bbox=dict(boxstyle="round,pad=0.25", facecolor=mcf.SURFACE,
                          edgecolor="none", alpha=0.85))


def _draw_pathway_boxes(ax, pathways, graph, layout, mcf, step=0.20):
    """Draw every enriched term as a rectangle, pushed apart and joined by a leader."""
    import numpy as np

    points = np.array([layout[t] for t in pathways], dtype=float)
    boxes = points + np.array([0.0, step])
    for _ in range(120):
        moved = False
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                delta = boxes[i] - boxes[j]
                distance = np.linalg.norm(delta)
                if distance < step * 2.4:
                    push = (delta / max(distance, 1e-6)) * step * 0.6
                    boxes[i] += push
                    boxes[j] -= push
                    moved = True
        if not moved:
            break

    for term, point, box in zip(pathways, points, boxes):
        ax.annotate(mcf._wrap_label(graph.nodes[term]["label"], 20), point,
                    xytext=box, fontsize=6.0, color=mcf.INK, ha="center",
                    va="center", zorder=9,
                    bbox=dict(boxstyle="square,pad=0.42",
                              facecolor=mcf.SURFACE, edgecolor=mcf.MUTED,
                              linewidth=0.9),
                    arrowprops=dict(arrowstyle="-", color=mcf.MUTED,
                                    linewidth=0.6, shrinkA=0, shrinkB=2))


def _draw_legend(ax, nutrients, colour_of, with_pathways, mcf, Line2D, disease):
    """One handle per micronutrient, plus the node kinds."""
    handles = [Line2D([], [], marker="s", linestyle="", markersize=8,
                      color=colour_of[n], label=n.split(" (")[0])
               for n in nutrients]
    handles += [
        Line2D([], [], marker="o", linestyle="", markersize=10,
               color=DISEASE_SEED_COLOR, markeredgecolor=SEED_RING_COLOR,
               markeredgewidth=DISEASE_SEED_RING_WIDTH,
               label=f"{disease} seed"),
        Line2D([], [], marker="o", linestyle="", markersize=8,
               color=BRIDGE_COLOR, markeredgecolor=SEED_RING_COLOR,
               markeredgewidth=SEED_RING_WIDTH, label="micronutrient seed"),
    ]
    if with_pathways:
        handles.append(Line2D([], [], marker="s", linestyle="", markersize=8,
                              color=mcf.SURFACE, markeredgecolor=mcf.MUTED,
                              label="enriched term"))
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="lower center",
              ncol=5, bbox_to_anchor=(0.5, -0.06), labelcolor=mcf.MUTED)


def _spring_layout(graph, nutrients, seed=0):
    """Relax the graph with the micronutrients pinned on a ring."""
    import numpy as np

    angles = np.linspace(0, 2 * np.pi, len(nutrients), endpoint=False)
    start = {n: np.array([np.cos(a), np.sin(a)]) * 1.9
             for n, a in zip(nutrients, angles)}
    rng = np.random.default_rng(seed)
    for node in graph.nodes:
        start.setdefault(node, rng.normal(scale=0.45, size=2))
    return nx.spring_layout(graph, pos=start, fixed=list(nutrients),
                            k=0.30, iterations=400, seed=seed)
