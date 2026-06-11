"""Exploratory: PCA on (DE union stability-knee) genes vs equally-many random genes.

Selection = main DE genes (union of main_ligands per-ligand contrasts) UNION the
stability-knee gene set. The knee is the n_features at the elbow of the mean_jaccard
curve in results/feature_selection/feature_count_stability.csv (FeatureSelector._elbow_index),
and its genes are the top-`knee` by selection_frequency in feature_count_core_genes.csv.
One random N-gene control (seed 42) from genes QC-passing in both batches. Mirrors
scripts/de_vs_random_pca.py; main pipeline untouched.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from src.feature_engineering import _elbow_index
from de_vs_random_pca import RESULTS, load_main_de_genes, build_conditions, make_pca

SEED = 42
STABILITY_CSV = RESULTS / "feature_selection" / "feature_count_stability.csv"
CORE_CSV = RESULTS / "feature_selection" / "feature_count_core_genes.csv"


def stability_knee_genes() -> set:
    """Top-`knee` genes by selection_frequency, knee = elbow of the mean_jaccard curve."""
    stab = pd.read_csv(STABILITY_CSV)
    knee = int(stab["n_features"].iloc[_elbow_index(stab["mean_jaccard"].to_numpy()) - 1])
    core = pd.read_csv(CORE_CSV).sort_values("selection_frequency", ascending=False)
    print(f"Stability knee = {knee} genes")
    return set(core["gene"].head(knee))


def main():
    op = sys.argv[1] if len(sys.argv) > 1 else "union"
    if op not in ("union", "intersect"):
        raise SystemExit("usage: de_stability_vs_random_pca.py [union|intersect]")
    out_subdir = "explore_de_stability" if op == "union" else "explore_de_stability_intersect"
    (RESULTS / "figures" / "pca" / out_subdir).mkdir(parents=True, exist_ok=True)

    de_genes = load_main_de_genes()
    stab_genes = stability_knee_genes()
    print(f"DE = {len(de_genes)}, stability = {len(stab_genes)}, overlap = {len(de_genes & stab_genes)}")

    conditions, Xm = build_conditions()
    combined = (de_genes | stab_genes) if op == "union" else (de_genes & stab_genes)
    selected = sorted(combined & set(Xm.columns))
    n = len(selected)
    print(f"\n=== N = {n} genes (DE {op} stability-knee) ===\n")

    qc_pool = set(Xm.columns[Xm.sum(axis=0) >= 10])
    if "external_test" in conditions:
        Xe = conditions["external_test"][1]
        qc_pool &= set(Xe.columns[Xe.sum(axis=0) >= 10])
    qc_genes = sorted(qc_pool)
    print(f"Random pool (QC-passing in both batches): {len(qc_genes)} genes")
    random_set = list(np.random.default_rng(SEED).choice(qc_genes, size=n, replace=False))

    summary = []
    for cond_key, (palette_key, X, y) in conditions.items():
        used = make_pca(cond_key, palette_key, X, y, selected, "destab", out_subdir=out_subdir, equal_aspect=True)
        used_r = make_pca(cond_key, palette_key, X, y, random_set, f"random_seed{SEED}", out_subdir=out_subdir, equal_aspect=True)
        print(f"{cond_key}: DE {op} stability used = {used}, random = {used_r}")
        summary.append({"condition": cond_key, "gene_set": "destab", "n_genes": used})
        summary.append({"condition": cond_key, "gene_set": f"random_seed{SEED}", "n_genes": used_r})

    out = RESULTS / "figures" / "pca" / out_subdir / "gene_counts.csv"
    pd.DataFrame(summary).to_csv(out, index=False)
    print(f"\nN (DE {op} stability-knee) = {n}")
    print(f"Figures: {RESULTS / 'figures' / 'pca' / out_subdir}")
    print(f"Gene-count summary: {out}")


if __name__ == "__main__":
    main()
