"""Feature engineering and transformation pipeline for RNA-seq data."""

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin, OneToOneFeatureMixin
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, SelectFromModel, mutual_info_classif
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import ExtraTreesClassifier
from feature_engine.selection import DropDuplicateFeatures


class QCLowerCountRemover(BaseEstimator, TransformerMixin):
    def fit(self, X, y=None):
        values = (
            X.to_numpy(dtype=float)
            if isinstance(X, pd.DataFrame)
            else np.asarray(X, dtype=float)
        )
        self.mask_ = values.sum(axis=0) >= 10
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
        if input_features is None:
            raise ValueError(
                "input_features must be provided for "
                "QCLowerCountRemover.get_feature_names_out"
            )
        return np.asarray(input_features, dtype=object)[self.mask_]


class LibraryLengthNormalizer(BaseEstimator, TransformerMixin):
    """Normalize gene counts to library size (reads per million)."""

    def fit(self, X, y=None):
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
    if k_best <= 0:
        raise ValueError(f"k_best must be positive, got {k_best}")
    if n_estimators <= 0:
        raise ValueError(f"n_estimators must be positive, got {n_estimators}")
    if max_depth <= 0:
        raise ValueError(f"max_depth must be positive, got {max_depth}")
    if max_features <= 0:
        raise ValueError(f"max_features must be positive, got {max_features}")

    en = ExtraTreesClassifier(
        n_estimators=n_estimators, max_depth=max_depth, random_state=random_state
    )

    pipe = Pipeline(
        [
            ("drop_low_count_features", QCLowerCountRemover()),
            ("normalise_for_library_size", LibraryLengthNormalizer()),
            ("select_k_best", SelectKBest(mutual_info_classif, k=k_best)),
            (
                "select_forest",
                SelectFromModel(en, max_features=max_features),
            ),
            ("standard_scale", StandardScaler()),
        ]
    )
    return pipe
