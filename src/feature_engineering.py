"""Feature engineering and transformation pipeline for RNA-seq data."""

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, SelectFromModel, chi2
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import ExtraTreesClassifier
from feature_engine.selection import DropDuplicateFeatures


class LibraryLengthNormalizer(BaseEstimator, TransformerMixin):
    """Normalize gene counts to library size (reads per million)."""

    def fit(self, X, y=None):
        """Fit method (no-op for this transformer).

        Args:
            X: Feature matrix.
            y: Target labels (optional).

        Returns:
            self
        """
        return self

    def transform(self, X, y=None):
        """Normalize counts to library size.

        Args:
            X: Feature matrix with gene counts.
            y: Target labels (optional).

        Returns:
            Normalized feature matrix.
        """
        # Normalize the gene counts to the library size
        X = X.apply(lambda x: (x / (x.sum() if x.sum() != 0 else 1)) * 1000000, axis=1)
        return X


def create_feature_pipeline(
    k_best: int = 1000,
    n_estimators: int = 250,
    max_depth: int = 5,
    max_features: int = 250,
    feature_threshold: float = 0.001,
    random_state: int = 42,
) -> Pipeline:
    """Create a feature selection and preprocessing pipeline.

    Args:
        k_best: Number of top features to select with chi2.
        n_estimators: Number of trees in ExtraTreesClassifier.
        max_depth: Maximum depth of trees.
        max_features: Maximum number of features after forest selection.
        feature_threshold: Threshold for feature importance in forest selection.
        random_state: Random state for reproducibility.

    Returns:
        Pipeline: Sklearn pipeline for feature engineering.
    """
    en = ExtraTreesClassifier(
        n_estimators=n_estimators, max_depth=max_depth, random_state=random_state
    )

    pipe = Pipeline(
        [
            ("drop_duplicates", DropDuplicateFeatures()),
            ("normalise_for_library_size", LibraryLengthNormalizer()),
            ("select_k_best", SelectKBest(chi2, k=k_best)),
            (
                "select_forest",
                SelectFromModel(
                    en, threshold=feature_threshold, max_features=max_features
                ),
            ),
            ("standard_scale", StandardScaler()),
        ]
    )

    return pipe
