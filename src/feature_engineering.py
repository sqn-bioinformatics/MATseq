from collections import defaultdict
from functools import partial
from typing import Dict, Sequence, Union

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
from sklearn.base import BaseEstimator, TransformerMixin, OneToOneFeatureMixin
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, SelectFromModel, mutual_info_classif
from sklearn.preprocessing import StandardScaler, FunctionTransformer
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

from .preprocessing import normalize_rpm

SEEDS = (9419, 3374, 7796)


class LibraryLengthNormalizer(OneToOneFeatureMixin, BaseEstimator, TransformerMixin):
    """Normalize gene counts to library size (reads per million)."""

    def fit(self, X, y=None):
        self.n_features_in_ = X.shape[1]
        if hasattr(X, "columns"):
            self.feature_names_in_ = np.asarray(X.columns, dtype=object)
        return self

    def transform(self, X):
        return normalize_rpm(X)


class ColumnSelector(BaseEstimator, TransformerMixin):
    """Keep a fixed set of named columns"""

    def __init__(self, genes):
        self.genes = genes

    def fit(self, X, y=None):
        self.columns_ = [g for g in self.genes if g in X.columns]
        if not self.columns_:
            raise ValueError("none of the requested genes are present in X")
        return self

    def transform(self, X):
        return X.reindex(columns=self.columns_)

    def get_feature_names_out(self, input_features=None):
        return np.asarray(self.columns_, dtype=object)


def preprocessing_pipeline() -> Pipeline:
    return Pipeline(
        [
            ("normalise_for_library_size", LibraryLengthNormalizer()),
            ("log1p", FunctionTransformer(np.log1p, feature_names_out="one-to-one")),
            ("standard_scale", StandardScaler()),
        ]
    )


def selection_pipeline(
    k_best,
    n_estimators,
    max_depth,
    max_features,
    random_state: int = 42,
) -> Pipeline:
    en = ExtraTreesClassifier(
        n_estimators=n_estimators, max_depth=max_depth,
        random_state=random_state, 
        n_jobs=-1, class_weight="balanced"
    )
    score_func = partial(mutual_info_classif, random_state=random_state)
    return Pipeline(
        [
            ("select_k_best", SelectKBest(score_func, k=k_best)),
            ("select_forest", SelectFromModel(en, max_features=max_features)),
        ]
    )


def feature_pipeline(**kwargs) -> Pipeline:
    return Pipeline(
        [
            *preprocessing_pipeline().steps,
            *selection_pipeline(**kwargs).steps,
        ]
    )


def selected_with_importance(fitted_pipeline: Pipeline) -> pd.DataFrame:
    """Extract selected genes ranked by ExtraTrees importance."""
    sfm = fitted_pipeline.named_steps["select_forest"]
    names_in = fitted_pipeline[:-1].get_feature_names_out()
    importances = sfm.estimator_.feature_importances_

    support = sfm.get_support()
    genes = names_in[support]
    gene_importances = importances[support]

    order = np.argsort(gene_importances)[::-1]
    return pd.DataFrame(
        {
            "gene": genes[order],
            "importance": gene_importances[order],
            "rank": np.arange(1, len(order) + 1),
        }
    )


def mutual_information(
    fitted_pipeline: Pipeline,
    X: pd.DataFrame,
    y: Union[np.ndarray, pd.Series],
    seeds: Sequence[int] = SEEDS,
) -> Dict:
    """MI curve over all genes plus its elbow, on the fitted preprocessing."""
    X_pre = fitted_pipeline.transform(X)
    mi_mean, per_run_elbows, curves = _mi_elbow_ranking(X_pre, y, seeds)
    mi_sorted = np.sort(mi_mean)[::-1]
    scores = {"rank": np.arange(1, len(mi_sorted) + 1), "mi_sorted": mi_sorted}
    for seed, curve in zip(seeds, curves):
        scores[f"mi_sorted_seed_{seed}"] = np.sort(curve)[::-1]
    return {
        "mi_elbow": elbow_index(mi_sorted),
        "per_run_elbows": per_run_elbows,
        "scores": pd.DataFrame(scores),
    }


def elbow_index(scores) -> int:
    y = np.asarray(scores, dtype=float)
    x = np.arange(len(y), dtype=float)
    xn = (x - x.min()) / (np.ptp(x) or 1)
    yn = (y - y.min()) / (np.ptp(y) or 1)
    num = np.abs(
        (yn[-1] - yn[0]) * xn
        - (xn[-1] - xn[0]) * yn
        + xn[-1] * yn[0]
        - yn[-1] * xn[0]
    )
    return int(np.argmax(num)) + 1


def _mi_elbow_ranking(X_pre, y, seeds: Sequence[int] = SEEDS):
    """Per-seed MI; returns the gene-wise mean MI and the per-seed elbows."""
    curves = list(
        tqdm(
            Parallel(n_jobs=len(seeds), return_as="generator")(
                delayed(mutual_info_classif)(X_pre, y, random_state=seed)
                for seed in seeds
            ),
            total=len(seeds),
            desc="Mutual information",
            unit="seed",
            dynamic_ncols=True,
        )
    )
    per_run_elbows = [elbow_index(np.sort(mi)[::-1]) for mi in curves]
    return np.mean(curves, axis=0), per_run_elbows, curves


def forest_kmeans(
    fitted_pipeline: Pipeline,
    X: pd.DataFrame,
    y: Union[np.ndarray, pd.Series],
    grid: Dict,
    k_best: int,
    seeds: Sequence[int] = SEEDS,
    random_state: int = 42,
) -> Dict:
    """Scan ExtraTrees settings by KMeans/ARI separation on the top k_best genes.
    """
    X_pre = fitted_pipeline.transform(X)
    X_k = (
        SelectKBest(
            partial(mutual_info_classif, random_state=random_state),
            k=k_best,
        )
        .set_output(transform="pandas")
        .fit_transform(X_pre, y)
    )

    triplets = [
        (n_estimators, max_depth, seed)
        for n_estimators in grid["n_estimators"]
        for max_depth in grid["max_depth"]
        for seed in seeds
    ]
    n_classes = len(np.unique(y))
    rows = defaultdict(list)

    for n_estimators, max_depth, seed in tqdm(
        triplets, desc="  Forest/k-means scan", unit="cell", dynamic_ncols=True
    ):
        et = ExtraTreesClassifier(
            n_estimators=n_estimators,
            max_depth=max_depth,
            random_state=seed,
        )
        et.fit(X_k, y)
        order = np.argsort(et.feature_importances_)[::-1]
        for n_selected in grid["n_selected"]:
            cols = X_k.columns[order[:n_selected]]
            labels = KMeans(
                n_clusters=n_classes, n_init=10, random_state=seed
            ).fit_predict(X_k[cols])
            rows[(n_estimators, max_depth, n_selected)].append(
                adjusted_rand_score(y, labels)
            )

    ari_cols = [f"ari_seed_{seed}" for seed in seeds]
    scan = pd.DataFrame(
        [
            {
                "n_estimators": ne,
                "max_depth": md,
                "n_selected": ns,
                **dict(zip(ari_cols, aris)),
                "ari_mean": float(np.mean(aris)),
                "ari_std": float(np.std(aris, ddof=1)),
            }
            for (ne, md, ns), aris in rows.items()
        ]
    ).sort_values(["n_estimators", "max_depth", "n_selected"]).reset_index(drop=True)

    best_row = scan.loc[scan["ari_mean"].idxmax()]
    best = {
        "n_estimators": int(best_row["n_estimators"]),
        "max_depth": int(best_row["max_depth"]),
        "n_selected": int(best_row["n_selected"]),
        "ari": float(best_row["ari_mean"]),
    }
    return {"best": best, "plateau": best["n_selected"], "scan": scan}
