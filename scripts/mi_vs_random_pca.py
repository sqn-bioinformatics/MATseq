"""Exploratory: PCA on the top-MI genes (MI elbow = 4063) vs equally-many random genes.

Selection = top 4063 genes by mutual information with the main_ligands labels, i.e.
the MI-elbow set from feature_count_analysis (results/feature_selection/feature_count_scores.csv,
elbow computed with FeatureSelector._elbow_index). The scores CSV has no gene names, so MI is
recomputed on the main subset (same computation: mean of mutual_info_classif over seeds 42-44)
and verified against the CSV. One random N-gene control (seed 42) is drawn from genes QC-passing
in both batches. Mirrors scripts/de_vs_random_pca.py; main pipeline untouched.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.feature_selection import mutual_info_classif

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from src import FeatureSelector
from de_vs_random_pca import RESULTS, build_conditions, make_pca

MI_TOP = 4063          # MI elbow over feature_count_scores.csv
OUT_SUBDIR = "explore_mi_random"
SEED = 42
SCORES_CSV = RESULTS / "feature_selection" / "feature_count_scores.csv"


def mi_ranked_genes(Xm, ym) -> pd.Series:
    """Genes ranked by MI with labels (mean of mutual_info_classif over seeds 42-44)."""
    X_pre = FeatureSelector.preprocessing_pipeline().set_output(transform="pandas").fit_transform(Xm)
    mi = np.mean(
        [mutual_info_classif(X_pre, ym, random_state=42 + i) for i in range(3)], axis=0
    )
    return pd.Series(mi, index=X_pre.columns).sort_values(ascending=False)


def main():
    (RESULTS / "figures" / "pca" / OUT_SUBDIR).mkdir(parents=True, exist_ok=True)

    conditions, Xm = build_conditions()
    ym = conditions["main_ligands"][2]

    mi_series = mi_ranked_genes(Xm, ym)
    csv_mi = pd.read_csv(SCORES_CSV)["mi_sorted"].dropna().to_numpy()
    diff = float(np.abs(np.sort(mi_series.to_numpy())[::-1][: len(csv_mi)] - csv_mi).max())
    print(f"MI reproduced: {len(mi_series)} genes, max|recomputed - CSV| = {diff:.2e}")

    top_genes = list(mi_series.index[:MI_TOP])
    n = len(top_genes)
    print(f"\n=== N = {n} top-MI genes (MI elbow) ===\n")

    qc_pool = set(Xm.columns[Xm.sum(axis=0) >= 10])
    if "external_test" in conditions:
        Xe = conditions["external_test"][1]
        qc_pool &= set(Xe.columns[Xe.sum(axis=0) >= 10])
    qc_genes = sorted(qc_pool)
    print(f"Random pool (QC-passing in both batches): {len(qc_genes)} genes")
    random_set = list(np.random.default_rng(SEED).choice(qc_genes, size=n, replace=False))

    summary = []
    for cond_key, (palette_key, X, y) in conditions.items():
        used = make_pca(cond_key, palette_key, X, y, top_genes, "mi", out_subdir=OUT_SUBDIR)
        used_r = make_pca(cond_key, palette_key, X, y, random_set, f"random_seed{SEED}", out_subdir=OUT_SUBDIR)
        print(f"{cond_key}: MI genes used = {used}, random = {used_r}")
        summary.append({"condition": cond_key, "gene_set": "mi", "n_genes": used})
        summary.append({"condition": cond_key, "gene_set": f"random_seed{SEED}", "n_genes": used_r})

    out = RESULTS / "figures" / "pca" / OUT_SUBDIR / "gene_counts.csv"
    pd.DataFrame(summary).to_csv(out, index=False)
    print(f"\nN (top-MI genes) = {n}")
    print(f"Figures: {RESULTS / 'figures' / 'pca' / OUT_SUBDIR}")
    print(f"Gene-count summary: {out}")


if __name__ == "__main__":
    main()
