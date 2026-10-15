#!/usr/bin/env python3
"""Make docs/micronutrient_coloc_pipeline_diagram.svg.

Follows dravet_cakut_micronutrients_netcoloc.ipynb, Part A to Part E. The
layout is computed rather than hand placed.
"""
import html
from pathlib import Path

OUT = (Path(__file__).resolve().parents[2] / "docs"
       / "micronutrient_coloc_pipeline_diagram.svg")

FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"
FS = 12          # body text
FS_SMALL = 9.5   # captions and section labels
LH = 15          # body line height
CHARW = 6.05     # mean advance width at FS, for wrapping
CHARW_SMALL = 4.9

INK = "#12171d"      # box outlines and body text
ACCENT = "#b3400c"   # colocalization, tissue filter, knowledge graph
MUTED = "#79828c"    # captions
RULE = "#d7dbe0"     # section rules

HAIR = 1.0
ACCENT_W = 1.6

W = 1000
MARGIN = 46
PAD_X = 14          # text inset inside a box
BOX_PAD_Y = 11      # text inset above and below
CAP_GAP = 13        # baseline of the first caption under a box
CAP_LH = 12.5
ROW_GAP = 34        # clear space between one row and the next
SECTION_GAP = 26    # extra space above a section label

# Three columns: the side columns carry the inputs and outputs, the centre the main chain.
CENTRE_W = 440
SIDE_W = 224
CX = W / 2
LX = MARGIN + SIDE_W / 2
RX = W - MARGIN - SIDE_W / 2


def esc(s):
    return html.escape(s)


def wrap(text, box_w, charw=CHARW, pad=PAD_X):
    """Greedy wrap to the usable width of a box."""
    maxc = max(6, int((box_w - 2 * pad) / charw))
    lines, cur = [], ""
    for word in text.split(" "):
        t = (cur + " " + word).strip()
        if len(t) <= maxc:
            cur = t
        else:
            if cur:
                lines.append(cur)
            cur = word
    if cur:
        lines.append(cur)
    return lines


class Node:
    """One box, its wrapped text and its captions. `y` is set by the layout."""

    def __init__(self, key, cx, w, text, captions=(), accent=False,
                 section=None, output=False):
        self.key = key
        self.cx = cx
        self.w = w
        self.lines = wrap(text, w)
        self.captions = list(captions)
        self.accent = accent
        self.section = section
        self.output = output
        self.h = len(self.lines) * LH + 2 * BOX_PAD_Y
        self.y = 0.0

    @property
    def top(self):
        return self.y

    @property
    def bottom(self):
        return self.y + self.h

    @property
    def mid(self):
        return self.y + self.h / 2

    @property
    def left(self):
        return self.cx - self.w / 2

    @property
    def right(self):
        return self.cx + self.w / 2

    @property
    def caption_bottom(self):
        """Bottom of the caption block, so the next row clears it."""
        if not self.captions:
            return self.bottom
        return self.bottom + CAP_GAP + (len(self.captions) - 1) * CAP_LH + 3


# The pipeline as rows. A row is a list of nodes laid out side by side; the
# row's height is the tallest node in it.
ROWS = [
    ("Part A - seed sets", [
        Node("disease", LX, SIDE_W,
             "2 disease sets: Dravet (Orphanet), CAKUT (Ozan et al. 2021)"),
        Node("micronutrient", CX, CENTRE_W,
             "10 micronutrient sets: CTD genes and pathway genes",
             ["results/ctd/ctd_micronutrient_sets.xlsx",
              "input/micronutrients/micronutrients_pathways.xlsx"]),
        Node("filesets", RX, SIDE_W,
             "1 dietary set: keto diet pathway genes",
             ["input/keto_diet/keto_diet_pathways.xlsx"]),
    ]),
    (None, [
        Node("sets", CX, CENTRE_W,
             "11 gene sets against 2 diseases: 22 pairs"),
    ]),
    (None, [
        Node("qc", CX, CENTRE_W,
             "Seed overlap, seed degree and seed sources",
             ["tables/seed_provenance.tsv",
              "tables/seed_overlap_micronutrients_{network}.tsv"]),
        Node("seedfig", RX, SIDE_W,
             "fig_seed_upset_{network}, fig_seed_sources",
             output=True),
    ]),
    ("Part B - tissue filter and propagation", [
        Node("tissue", CX, CENTRE_W,
             "Tissue filter (HPA nTPM > 1) on the interactome and the seeds: "
             "brain for Dravet, kidney for CAKUT",
             ["input/hpa/rna_brain_region_hpa.tsv",
              "input/hpa/rna_tissue_consensus.tsv"],
             accent=True),
    ]),
    (None, [
        Node("nps", CX, CENTRE_W,
             "Network propagation per seed set, degree-matched null, "
             "1000 reps -> NPS z",
             ["scripts/propagate.py", "cached in work/nps_sets/{network}"]),
    ]),
    ("Part C - colocalized genes", [
        Node("coloc", CX, CENTRE_W,
             "Colocalization of Dravet-micronutrients and "
             "CAKUT-micronutrients",
             ["scripts/run_one_pair.py", "one interpreter per pair, cached",
              "tables/gene_matrix_{disease}.tsv"],
             accent=True),
        Node("heat", RX, SIDE_W, "fig_gene_heatmap_{disease}",
             output=True),
    ]),
    ("Part D - enrichment", [
        Node("orsum", LX, SIDE_W,
             "orsum: representative terms per disease and source"),
        Node("gprofiler", CX, CENTRE_W,
             "g:Profiler (GO:BP, GO:CC, REAC) per significant pair, against "
             "the nodes of the tissue network as background",
             ["tables/coloc_annotation_all_pairs.tsv",
              "tables/coloc_annotation_{disease}.tsv"]),
        Node("enrich", RX, SIDE_W, "fig_enrichment_{disease}",
             output=True),
    ]),
    ("Part E - the knowledge graph of each disease", [
        Node("kg", CX, CENTRE_W,
             "Colocalized module genes of the significant pairs and their "
             "top 5 terms per micronutrient as one graph of micronutrient, "
             "gene and pathway nodes",
             ["scripts/dravet_coloc_kg.py"],
             accent=True),
        Node("kgfig", RX, SIDE_W, "fig_{disease}_kg_full", output=True),
    ]),
    (None, [
        Node("shared", CX, CENTRE_W,
             "Genes that the most micronutrient sets reach",
             ["tables/{disease}_shared_genes.tsv"]),
        Node("sharedfig", RX, SIDE_W, "fig_{disease}_kg_shared_genes",
             output=True),
    ]),
    (None, [
        Node("export", CX, CENTRE_W,
             "Export of the full graph and of the core, the genes reached by "
             "at least 2 micronutrients",
             ["{disease}_coloc_kg.graphml",
              "{disease}_coloc_kg_core.graphml"],
             accent=True),
    ]),
]

# Resolve the vertical position of every row.
N = {}
sections = []
y = 96
for title, row in ROWS:
    if title:
        y += SECTION_GAP
        sections.append((title, y))
        y += 16
    for node in row:
        node.y = y
        N[node.key] = node
    y = max(n.caption_bottom for n in row) + ROW_GAP
H = y + 24

svg = []


def add(s):
    svg.append(s)


def draw_box(n):
    stroke = ACCENT if n.accent else (MUTED if n.output else INK)
    dash = ' stroke-dasharray="4 3"' if n.output else ""
    add(f'<rect x="{n.left:.1f}" y="{n.top:.1f}" width="{n.w}" '
        f'height="{n.h:.1f}" rx="3.5" fill="#ffffff" stroke="{stroke}" '
        f'stroke-width="{ACCENT_W if n.accent else HAIR}"{dash}/>')
    y0 = n.top + BOX_PAD_Y + LH - 4
    for i, line in enumerate(n.lines):
        add(f'<text x="{n.cx:.1f}" y="{y0 + i * LH:.1f}" text-anchor="middle" '
            f'font-size="{FS}" fill="{stroke}">{esc(line)}</text>')
    cap_x, anchor = n.cx, "middle"
    if n.captions:
        if n.right > W - MARGIN - 60:
            cap_x, anchor = n.right, "end"
        else:
            cap_x, anchor = n.cx + 10, "start"
    for i, cap in enumerate(n.captions):
        add(f'<text x="{cap_x:.1f}" y="{n.bottom + CAP_GAP + i * CAP_LH:.1f}" '
            f'text-anchor="{anchor}" font-size="{FS_SMALL}" fill="{MUTED}">'
            f'{esc(cap)}</text>')


def section(title, y):
    add(f'<text x="{MARGIN}" y="{y:.1f}" font-size="{FS_SMALL}" '
        f'font-weight="700" letter-spacing="1.5" fill="{MUTED}">'
        f'{esc(title.upper())}</text>')
    tw = len(title) * (CHARW_SMALL * 1.34 + 1.5) + 34
    add(f'<path d="M{MARGIN + tw:.1f} {y - 4:.1f} L{W - MARGIN} {y - 4:.1f}" '
        f'stroke="{RULE}" stroke-width="{HAIR}"/>')


def down(a, b):
    """Vertical connector between two boxes in the centre column.

    The captions of a box sit to the right of this line, so it runs straight
    down the centre without striking through them.
    """
    add(f'<path d="M{a.cx:.1f} {a.bottom:.1f} L{a.cx:.1f} {b.top - 2:.1f}" '
        f'fill="none" stroke="{INK}" stroke-width="{HAIR}" '
        f'marker-end="url(#arw)"/>')


def side_in(a, b):
    """Horizontal arrow from box a to box b on the same row, at the mid of the shorter box."""
    y = a.mid if a.h < b.h else b.mid
    x1 = a.left if a.cx > b.cx else a.right
    x2 = b.right if a.cx > b.cx else b.left
    sign = -1 if a.cx > b.cx else 1
    add(f'<path d="M{x1:.1f} {y:.1f} L{x2 + sign * 2:.1f} {y:.1f}" '
        f'fill="none" stroke="{INK}" stroke-width="{HAIR}" '
        f'marker-end="url(#arw)"/>')


side_out = side_in


def elbow_down(a, b, drop=None):
    """Out of the bottom of a side box, along, then into the target's side."""
    y_mid = drop if drop is not None else b.mid
    x2 = b.right if a.cx > b.cx else b.left
    sign = -1 if a.cx > b.cx else 1
    add(f'<path d="M{a.cx:.1f} {a.caption_bottom:.1f} L{a.cx:.1f} {y_mid:.1f} '
        f'L{x2 + sign * 2:.1f} {y_mid:.1f}" fill="none" stroke="{INK}" '
        f'stroke-width="{HAIR}" marker-end="url(#arw)"/>')


add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H:.0f}" '
    f'width="{W}" height="{H:.0f}" font-family="{FONT}">')
add('<defs><marker id="arw" viewBox="0 0 10 10" refX="8.5" refY="5" '
    'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
    f'<path d="M1.5 1.5L8 5L1.5 8.5" fill="none" stroke="{INK}" '
    'stroke-width="1.4" stroke-linecap="round" stroke-linejoin="round"/>'
    '</marker></defs>')
add(f'<rect x="0" y="0" width="{W}" height="{H:.0f}" fill="#ffffff"/>')

add(f'<text x="{CX}" y="46" text-anchor="middle" font-size="16" '
    f'font-weight="700" fill="{INK}">Dravet and CAKUT against the '
    f'micronutrient gene sets</text>')
add(f'<text x="{CX}" y="66" text-anchor="middle" font-size="{FS_SMALL + 0.5}" '
    f'fill="{MUTED}">netcoloc on PCNet 2.2  ·  '
    f'dravet_cakut_micronutrients_netcoloc.ipynb  ·  '
    f'figures as PNG and PDF</text>')

for title, y in sections:
    section(title, y)

# The three seed families meet at the pair grid.
elbow_down(N["disease"], N["sets"])
down(N["micronutrient"], N["sets"])
elbow_down(N["filesets"], N["sets"])

down(N["sets"], N["qc"])
side_out(N["qc"], N["seedfig"])
down(N["qc"], N["tissue"])
down(N["tissue"], N["nps"])
down(N["nps"], N["coloc"])
side_out(N["coloc"], N["heat"])

down(N["coloc"], N["gprofiler"])
side_out(N["gprofiler"], N["enrich"])
side_out(N["gprofiler"], N["orsum"])

# The module and its enriched terms both become knowledge-graph nodes.
down(N["gprofiler"], N["kg"])
side_out(N["kg"], N["kgfig"])
down(N["kg"], N["shared"])
side_out(N["shared"], N["sharedfig"])
down(N["shared"], N["export"])

# Boxes last, so their white fill masks the connector ends.
for _, row in ROWS:
    for node in row:
        draw_box(node)

# Accent key.
key_y = H - 20
add(f'<rect x="{W - MARGIN - 250:.1f}" y="{key_y - 8:.1f}" width="14" '
    f'height="9.5" rx="2" fill="none" stroke="{ACCENT}" '
    f'stroke-width="{ACCENT_W}"/>')
add(f'<text x="{W - MARGIN - 230:.1f}" y="{key_y:.1f}" '
    f'font-size="{FS_SMALL}" fill="{MUTED}">colocalization, tissue filter '
    f'and knowledge graph</text>')

add('</svg>')
OUT.write_text("\n".join(svg))
print(f"wrote {OUT}")
