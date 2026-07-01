from functools import partial
from pathlib import Path
from typing import Dict, Union

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin, OneToOneFeatureMixin
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, SelectFromModel, mutual_info_classif
from sklearn.preprocessing import StandardScaler, FunctionTransformer
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.utils.validation import (
    _check_feature_names,
    _check_n_features,
)

from .preprocessing import normalize_rpm


def _record_input(estimator, X):
    _check_n_features(estimator, X, reset=True)
    _check_feature_names(estimator, X, reset=True)


class LibraryLengthNormalizer(OneToOneFeatureMixin, BaseEstimator, TransformerMixin):
    """Normalize gene counts to library size (reads per million)."""

    def fit(self, X, y=None):
        _record_input(self, X)
        return self

    def transform(self, X):
        """Normalize counts to library size (RPM)."""
        if X is None:
            raise ValueError("X cannot be None")
        return normalize_rpm(X)


class ColumnSelector(BaseEstimator, TransformerMixin):
    """Keep a fixed set of named columns; used to deploy a frozen gene list."""

    def __init__(self, genes):
        self.genes = list(genes)

    def fit(self, X, y=None):
        if not isinstance(X, pd.DataFrame):
            raise TypeError("ColumnSelector requires a DataFrame input")
        self.columns_ = [g for g in self.genes if g in X.columns]
        if not self.columns_:
            raise ValueError("none of the requested genes are present in X")
        _record_input(self, X)
        return self

    def transform(self, X):
        if not isinstance(X, pd.DataFrame):
            raise TypeError("ColumnSelector requires a DataFrame input")
        return X.reindex(columns=self.columns_)

    def get_feature_names_out(self, input_features=None):
        return np.asarray(self.columns_, dtype=object)


def _prepare_output_dir(output_dir: Path = None) -> Path:
    if output_dir is None:
        output_dir = Path.cwd() / "results" / "feature_analysis"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def _elbow_index(scores) -> int:
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


class FeatureSelector:
    """Feature-selection pipelines plus stability and count analyses."""

    def __init__(self, output_dir: Path = None):
        """Initialize selector.

        Args:
            output_dir: Directory for output files. Defaults to results/feature_analysis.
        """
        self.output_dir = _prepare_output_dir(output_dir)
        self.cv_results_ = None

    @staticmethod
    def preprocessing_pipeline() -> Pipeline:
        return Pipeline(
            [
                ("normalise_for_library_size", LibraryLengthNormalizer()),
                ("log1p", FunctionTransformer(np.log1p, feature_names_out="one-to-one")),
                ("standard_scale", StandardScaler()),
            ]
        )

    @staticmethod
    def selection_pipeline(
        k_best,
        n_estimators,
        max_depth,
        max_features,
        random_state: int = 42,
    ) -> Pipeline:
        en = ExtraTreesClassifier(
            n_estimators=n_estimators, max_depth=max_depth, random_state=random_state
        )
        score_func = partial(mutual_info_classif, random_state=random_state)
        return Pipeline(
            [
                ("select_k_best", SelectKBest(score_func, k=k_best)),
                ("select_forest", SelectFromModel(en, max_features=max_features)),
            ]
        )

    @staticmethod
    def feature_pipeline(**kwargs) -> Pipeline:
        return Pipeline(
            [
                *FeatureSelector.preprocessing_pipeline().steps,
                *FeatureSelector.selection_pipeline(**kwargs).steps,
            ]
        )

    @staticmethod
    def count_analysis(
        X: pd.DataFrame,
        y: Union[np.ndarray, pd.Series],
        random_state: int = 42,
    ) -> Dict:
        pre = FeatureSelector.preprocessing_pipeline().set_output(transform="pandas")
        X_pre = pre.fit_transform(X, y)
        mi_sorted, per_run_elbows = FeatureSelector._mi_elbow_ranking(
            X_pre, y, random_state=random_state
        )
        return {
            "mi_elbow": int(round(float(np.mean(per_run_elbows)))),
            "per_run_elbows": per_run_elbows,
            "scores": pd.DataFrame(
                {"rank": np.arange(1, len(mi_sorted) + 1), "mi_sorted": mi_sorted}
            ),
        }

    @staticmethod
    def _mi_elbow_ranking(X_pre, y, n_seeds: int = 3, random_state: int = 42):
        """Per-seed sorted MI curves; returns the averaged curve and per-seed elbows."""
        sorted_curves = []
        per_run_elbows = []
        for i in range(n_seeds):
            mi = mutual_info_classif(X_pre, y, random_state=random_state + i)
            mi_sorted = np.sort(mi)[::-1]
            sorted_curves.append(mi_sorted)
            per_run_elbows.append(_elbow_index(mi_sorted))
        averaged_sorted = np.mean(sorted_curves, axis=0)
        return averaged_sorted, per_run_elbows

    @staticmethod
    def forest_kmeans_separation(
        X: pd.DataFrame,
        y: Union[np.ndarray, pd.Series],
        k_best: int,
        grid: Dict,
        random_state: int = 42,
    ) -> Dict:
        pre = FeatureSelector.preprocessing_pipeline().set_output(transform="pandas")
        X_pre = pre.fit_transform(X, y)
        k = min(k_best, X_pre.shape[1])
        X_k = (
            SelectKBest(partial(mutual_info_classif, random_state=random_state), k=k)
            .set_output(transform="pandas")
            .fit_transform(X_pre, y)
        )
        n_classes = len(np.unique(y))

        rows = []
        for n_estimators in grid["n_estimators"]:
            for max_depth in grid["max_depth"]:
                et = ExtraTreesClassifier(
                    n_estimators=n_estimators,
                    max_depth=max_depth,
                    random_state=random_state,
                )
                et.fit(X_k, y)
                order = np.argsort(et.feature_importances_)[::-1]
                for n_selected in grid["n_selected"]:
                    n = min(n_selected, X_k.shape[1])
                    cols = X_k.columns[order[:n]]
                    labels = KMeans(
                        n_clusters=n_classes, n_init=10, random_state=random_state
                    ).fit_predict(X_k[cols])
                    rows.append(
                        {
                            "n_estimators": int(n_estimators),
                            "max_depth": int(max_depth),
                            "n_selected": int(n_selected),
                            "ari": float(adjusted_rand_score(y, labels)),
                        }
                    )

        scan = pd.DataFrame(rows)
        best_row = scan.loc[scan["ari"].idxmax()]
        best = {
            "n_estimators": int(best_row["n_estimators"]),
            "max_depth": int(best_row["max_depth"]),
            "n_selected": int(best_row["n_selected"]),
            "ari": float(best_row["ari"]),
        }
        return {"best": best, "scan": scan}
