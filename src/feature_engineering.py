"""Feature engineering and transformation pipeline for RNA-seq data."""

from functools import partial

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin, OneToOneFeatureMixin
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, SelectFromModel, mutual_info_classif
from sklearn.preprocessing import StandardScaler, FunctionTransformer
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.utils.validation import (
    _check_feature_names,
    _check_feature_names_in,
    _check_n_features,
)
from feature_engine.selection import DropDuplicateFeatures


def _record_input(estimator, X):
    _check_n_features(estimator, X, reset=True)
    _check_feature_names(estimator, X, reset=True)


class QCLowerCountRemover(BaseEstimator, TransformerMixin):
    def fit(self, X, y=None):
        values = (
            X.to_numpy(dtype=float)
            if isinstance(X, pd.DataFrame)
            else np.asarray(X, dtype=float)
        )
        self.mask_ = values.sum(axis=0) >= 10
        _record_input(self, X)
        return self

    def transform(self, X):
        if X is None:
            raise ValueError("X cannot be None")
        X_is_df = isinstance(X, pd.DataFrame)

        if X_is_df:
            values = X.to_numpy(dtype=float)
            index = X.index
            columns = X.columns
        else:
            values = np.asarray(X, dtype=float)
            index = None
            columns = None

        if values.ndim != 2:
            raise ValueError(f"Expected 2D array, got shape {values.shape}")

        filtered_values = values[:, self.mask_]
        if X_is_df:
            return pd.DataFrame(
                filtered_values,
                index=index,
                columns=columns[self.mask_],
            )

        return filtered_values

    def get_feature_names_out(self, input_features=None):
        input_features = _check_feature_names_in(self, input_features)
        return np.asarray(input_features, dtype=object)[self.mask_]


class LibraryLengthNormalizer(OneToOneFeatureMixin, BaseEstimator, TransformerMixin):
    """Normalize gene counts to library size (reads per million)."""

    def fit(self, X, y=None):
        _record_input(self, X)
        return self

    def transform(self, X):
        """Normalize counts to library size (RPM)."""
        if X is None:
            raise ValueError("X cannot be None")

        X_is_df = isinstance(X, pd.DataFrame)

        if X_is_df:
            values = X.to_numpy(dtype=float)
            index = X.index
            columns = X.columns
        else:
            values = np.asarray(X, dtype=float)
            index = None
            columns = None

        if values.ndim != 2:
            raise ValueError(f"Expected 2D array, got shape {values.shape}")

        totals = values.sum(axis=1, keepdims=True)
        totals[totals == 0.0] = 1.0
        normalized = (values / totals) * 1_000_000.0

        if X_is_df:
            return pd.DataFrame(
                normalized,
                index=index,
                columns=columns,
            )

        return normalized


def create_preprocessing_pipeline(
    drop_low_count: bool = True, scale: bool = False
) -> Pipeline:
    steps = []
    if drop_low_count:
        steps.append(("drop_low_count_features", QCLowerCountRemover()))
    steps += [
        ("normalise_for_library_size", LibraryLengthNormalizer()),
        ("log1p", FunctionTransformer(np.log1p, feature_names_out="one-to-one")),
    ]
    if scale:
        steps.append(("standard_scale", StandardScaler()))
    return Pipeline(steps)


def create_selection_pipeline(
    k_best: int = 1000,
    n_estimators: int = 250,
    max_depth: int = 5,
    max_features: int = 250,
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


def create_feature_pipeline(
    k_best: int = 1000,
    n_estimators: int = 250,
    max_depth: int = 5,
    max_features: int = 250,
    random_state: int = 42,
) -> Pipeline:
    """Create a feature selection and preprocessing pipeline.

    Args:
        k_best: Number of top features to select with chi2.
        n_estimators: Number of trees in ExtraTreesClassifier.
        max_depth: Maximum depth of trees.
        max_features: Maximum number of features after forest selection.
        random_state: Random state for reproducibility.

    Returns:
        Pipeline: Sklearn pipeline for feature engineering.
    """
    return Pipeline(
        [
            *create_preprocessing_pipeline(drop_low_count=True, scale=False).steps,
            *create_selection_pipeline(
                k_best, n_estimators, max_depth, max_features, random_state
            ).steps,
            ("standard_scale", StandardScaler()),
        ]
    )
