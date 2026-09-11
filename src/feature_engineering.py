from functools import partial
from typing import Dict, Optional, Sequence, Union

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, SelectFromModel, mutual_info_classif
from sklearn.preprocessing import StandardScaler, FunctionTransformer
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

from .preprocessing import normalize_rpm

SEEDS = (9419, 3374, 7796)


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
            (
                "normalise_for_library_size",
                FunctionTransformer(normalize_rpm, feature_names_out="one-to-one"),
            ),
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
        random_state=random_state, n_jobs=-1, class_weight="balanced"
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
    X_pre: pd.DataFrame,
    y: Union[np.ndarray, pd.Series],
    seeds: Sequence[int] = SEEDS,
) -> Dict:
    """Seed-averaged MI curve over all preprocessed genes plus its elbow."""
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
    mi_sorted = np.sort(np.mean(curves, axis=0))[::-1]
    scores = {"rank": np.arange(1, len(mi_sorted) + 1), "mi_sorted": mi_sorted}
    for seed, curve in zip(seeds, curves):
        scores[f"mi_sorted_seed_{seed}"] = np.sort(curve)[::-1]
    return {"mi_elbow": elbow_index(mi_sorted), "scores": pd.DataFrame(scores)}


def elbow_index(scores) -> int:
    y = np.asarray(scores, dtype=float)
    xn = np.linspace(0.0, 1.0, len(y))
    yn = (y - y.min()) / (np.ptp(y) or 1)
    return int(np.argmax(np.abs(yn[0] + (yn[-1] - yn[0]) * xn - yn))) + 1


def forest_kmeans(
    X_pre: pd.DataFrame,
    y: Union[np.ndarray, pd.Series],
    k_best: int,
    n_estimators: int,
    max_depth: Optional[int],
    seeds: Sequence[int] = (*SEEDS, 6434, 7362),
    random_state: int = 42,
) -> pd.DataFrame:
    """k-means/ligand ARI per ExtraTrees gene rank within the top k_best MI genes."""
    X_k = (
        SelectKBest(
            partial(mutual_info_classif, random_state=random_state),
            k=k_best,
        )
        .set_output(transform="pandas")
        .fit_transform(X_pre, y)
    )
    ari_cols = [f"ari_seed_{seed}" for seed in seeds]
    n_classes = len(np.unique(y))
    scan = pd.DataFrame({"n_selected": np.arange(1, k_best + 1)})
    for col, seed in tqdm(
        zip(ari_cols, seeds),
        total=len(seeds),
        desc="  k-means per gene rank",
        unit="seed",
        dynamic_ncols=True,
    ):
        forest = ExtraTreesClassifier(
            n_estimators=n_estimators,
            max_depth=max_depth,
            random_state=seed,
            n_jobs=-1,
            class_weight="balanced",
        ).fit(X_k, y)
        X_ranked = X_k.to_numpy()[:, np.argsort(forest.feature_importances_)[::-1]]
        labels = Parallel(n_jobs=-1)(
            delayed(
                KMeans(n_clusters=n_classes, n_init=10, random_state=seed).fit_predict
            )(X_ranked[:, :n])
            for n in scan["n_selected"]
        )
        scan[col] = [adjusted_rand_score(y, lab) for lab in labels]
    scan["ari_mean"] = scan[ari_cols].mean(axis=1)
    scan["ari_std"] = scan[ari_cols].std(axis=1)
    return scan
