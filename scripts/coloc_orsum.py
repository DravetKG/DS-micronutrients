"""orsum over the colocalization enrichment, one ranked term list per micronutrient.

orsum collapses a ranked list of terms into representative terms, using the
GMT the enrichment was run against. One run holds the micronutrients of a single
disease, so the collapsed terms are comparable over the micronutrients and the
two diseases stay separate. The uncollapsed table stays untouched, so both views
are available.
"""

import os
import subprocess
from pathlib import Path

import pandas as pd

import geneset_coloc
import micronutrient_coloc as mc

ROOT = Path(__file__).resolve().parent.parent
GMT_DIR = ROOT / "input" / "gprofiler_gmt"

# The g:Profiler sources that are collapsed, and the GMT of each one.
SOURCES = {
    "GO:BP": "hsapiens.GO:BP.name.gmt",
    "GO:CC": "hsapiens.GO:CC.name.gmt",
    "REAC": "hsapiens.REAC.name.gmt",
}


def out_dir():
    """The orsum sub-tree of the current network."""
    return mc.OUT / "orsum"


def gmt_path(source):
    path = GMT_DIR / SOURCES[source]
    if not path.exists():
        raise FileNotFoundError(f"missing GMT for {source}: {path}")
    return path


def micronutrients_of(ann_table, disease, source):
    """The micronutrients of one disease and source that carry a term."""
    subset = ann_table[(ann_table["disease"] == geneset_coloc.LABELS[disease])
                       & (ann_table["source"] == source)]
    return list(dict.fromkeys(subset["micronutrient"]))


def ranked_terms(ann_table, disease, source, micronutrient):
    """The terms of one disease, source and micronutrient, best p-value first.

    A term reported twice under one micronutrient enters once, with its smallest
    p-value.
    """
    subset = ann_table[(ann_table["disease"] == geneset_coloc.LABELS[disease])
                       & (ann_table["source"] == source)
                       & (ann_table["micronutrient"] == micronutrient)]
    if not len(subset):
        return pd.DataFrame(columns=["native", "name", "p_value",
                                     "term_size", "intersection_size"])
    best = (subset.groupby("native")
            .agg(name=("name", "first"), p_value=("p_value", "min"),
                 term_size=("term_size", "first"),
                 intersection_size=("intersection_size", "max"))
            .reset_index()
            .sort_values("p_value"))
    return best.reset_index(drop=True)


def write_input(terms, path):
    """Write a ranked term-ID list, the input format orsum reads."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(terms["native"]) + "\n")
    return path


def _slug(text):
    """A file name safe form of a label."""
    keep = [c if c.isalnum() else "_" for c in text.lower()]
    return "".join(keep).strip("_").replace("__", "_")


def run(ann_table, sources=None, max_rep_size=2000, min_term_size=10,
        max_term_size=None, n_plot=50, diseases=None):
    """Collapse the micronutrients of each disease, one orsum run per source.

    Every micronutrient of a disease is a separate input file, so the orsum
    plots hold one column per micronutrient. The diseases run separately.
    """
    base = out_dir()
    summary = {}
    for disease in diseases or mc.DISEASES:
        for source in sources or SOURCES:
            files, aliases, ranked = [], [], {}
            tag = f"{disease}_{source.replace(':', '').lower()}"
            for micronutrient in micronutrients_of(ann_table, disease, source):
                terms = ranked_terms(ann_table, disease, source, micronutrient)
                if not len(terms):
                    continue
                path = base / tag / f"{_slug(micronutrient)}_terms.txt"
                write_input(terms, path)
                terms.to_csv(path.with_name(f"{_slug(micronutrient)}_uncollapsed.tsv"),
                             sep="\t", index=False)
                files.append(str(path))
                aliases.append(micronutrient)
                ranked[micronutrient] = terms
            if not files:
                print(f"{geneset_coloc.LABELS[disease]} / {source}: no term to collapse")
                continue

            result_dir = base / tag / "orsum"
            result_dir.mkdir(parents=True, exist_ok=True)
            command = ["orsum.py", "--gmt", str(gmt_path(source)),
                       "--files", *files, "--fileAliases", *aliases,
                       "--outputFolder", str(result_dir),
                       "--maxRepSize", str(max_rep_size),
                       "--minTermSize", str(min_term_size),
                       "--numberOfTermsToPlot", str(n_plot)]
            if max_term_size is not None:
                command += ["--maxTermSize", str(max_term_size)]
            print(f"\n{geneset_coloc.LABELS[disease]} / {source}: collapsing "
                  f"{len(files)} micronutrients, "
                  f"{sum(len(t) for t in ranked.values())} terms")
            # orsum draws its plots itself, which needs a non-interactive backend.
            env = {**os.environ, "MPLBACKEND": "Agg"}
            done = subprocess.run(command, capture_output=True, text=True, env=env)
            if done.returncode:
                print(done.stdout[-2000:])
                print(done.stderr[-2000:])
                raise RuntimeError(f"orsum failed on {disease} / {source}")
            summary[(disease, source)] = {"input": ranked, "output": result_dir}
    return summary


def collapsed_table(summary, disease, source):
    """The representative terms orsum wrote for one disease and source."""
    result_dir = Path(summary[(disease, source)]["output"])
    matches = sorted(result_dir.glob("filteredResult-Summary.tsv"))
    if not matches:
        matches = sorted(result_dir.glob("*Summary.tsv"))
    if not matches:
        raise FileNotFoundError(f"no orsum summary in {result_dir}")
    return pd.read_csv(matches[0], sep="\t")
