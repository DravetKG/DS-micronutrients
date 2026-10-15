"""Run one disease-micronutrient pair in its own process and cache the result.

    python scripts/run_one_pair.py <disease> <micronutrient> [cache_dir]

The pair runs on the tissue network of the disease. The outputs go next to
``cache_dir``.
"""
import pickle
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import geneset_coloc
import micronutrient_coloc as mc

disease, nutrient = sys.argv[1], sys.argv[2]
mc.use_disease_network(disease)
cache = Path(sys.argv[3]) if len(sys.argv) > 3 else mc.CACHE
out_dir = cache.parent
cache.mkdir(parents=True, exist_ok=True)
out = cache / f"{disease}_{nutrient}.pkl"
if out.exists():
    print("cached", out.name); sys.exit(0)

geneset_coloc.set_pair(f"{out_dir.name}/{cache.name}/{disease}_{nutrient}")
res = geneset_coloc.main_geneset_pair(disease, nutrient, verbose=False)
tsv = out_dir / "tables" / f"genes_{disease}_{nutrient}.tsv"
tsv.parent.mkdir(parents=True, exist_ok=True)
res["coloc_genes"].to_csv(tsv, sep="\t", index=False)

_, _, _, degrees = geneset_coloc.network()
stats = res["stats"].iloc[0].to_dict()
pickle.dump({
    "disease": disease, "nutrient": nutrient, "stats": stats,
    "seeds_disease": res[f"seeds_{disease}"], "seeds_nutrient": res[f"seeds_{nutrient}"],
    "coloc_genes": res["coloc_genes"],
    "median_seed_degree": mc._median_seed_degree(res[f"seeds_{nutrient}"], degrees),
}, open(out, "wb"))
print("wrote", out.name, "coloc", int(stats["coloc_size"]), "z", round(stats["z_mean_z"],2))
