"""Model training, evaluation, and prediction for multiclass classification."""

import warnings
import numpy as np
import pandas as pd
import pickle
from pathlib import Path
from typing import Dict, Tuple, Optional
from sklearn.linear_model import SGDClassifier
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from xgboost import XGBClassifier
from imblearn.over_sampling import SMOTE
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    ConfusionMatrixDisplay,
)
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.base import clone
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message="Liblinear failed to converge")


def multiclass_roc_auc_score(y_test, y_pred, average: str = "macro") -> float:
    """Calculate averaged AUC per class for multiclass classification.

    Args:
        y_test: True labels.
        y_pred: Predicted labels.
        average: Averaging method for AUC calculation.

    Returns:
        float: Average ROC AUC score across all classes.
    """
    unique_classes = set(y_test)
    roc_auc_dict = {}

    for tested_class in unique_classes:
        other_classes = [x for x in unique_classes if x != tested_class]

        binary_y_test = [0 if x in other_classes else 1 for x in y_test]
        binary_y_pred = [0 if x in other_classes else 1 for x in y_pred]

        roc_auc = roc_auc_score(binary_y_test, binary_y_pred, average=average)
        roc_auc_dict[tested_class] = roc_auc

    return sum(roc_auc_dict.values()) / len(roc_auc_dict.values())


def make_score(y_test, y_pred) -> dict:
    """Calculate multiple classification metrics.

    Args:
        y_test: True labels.
        y_pred: Predicted labels.

    Returns:
        dict: Dictionary of metric names to scores.
    """
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred, average="macro", zero_division=np.nan)
    recall = recall_score(y_test, y_pred, average="macro", zero_division=np.nan)
    f1 = f1_score(y_test, y_pred, average="macro", zero_division=np.nan)
    roc_auc = multiclass_roc_auc_score(y_test, y_pred)

    return {
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "roc_auc": roc_auc,
    }


def get_confusion_matrix(test_pred, y_test, name: str, output_dir: Path = None) -> None:
    """Generate and save confusion matrix visualization.

    Args:
        test_pred: Predicted labels.
        y_test: True labels.
        name: Name for the plot title and file.
        output_dir: Directory to save figure.
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    ConfusionMatrixDisplay.from_predictions(
        y_test,
        test_pred,
        ax=ax,
        xticks_rotation="vertical",
        colorbar=False,
        normalize="true",
        values_format=".0%",
    )

    ax.set_title(f"Confusion Matrix {name}")
    plt.tight_layout()

    output_dir = Path(output_dir) if output_dir else Path("results/figures")
    output_dir.mkdir(parents=True, exist_ok=True)
    save_path = output_dir / f"Confusion_Matrix_{name}.png"
    try:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Figure saved: {save_path}")
    finally:
        plt.close()


class ModelFactory:
    """Factory for creating and initializing classifier models."""

    @staticmethod
    def create_models(
        random_state: int = 42,
        calibrate: bool = True,
    ) -> Dict:
        """Create a dictionary of classifier models.

        Args:
            random_state: Random state for reproducibility.
            calibrate: Whether to use CalibratedClassifierCV for SVM.

        Returns:
            Dictionary of model names to model instances.
        """
        models = {}

        svc = LinearSVC(
            max_iter=50000, tol=1e-3, random_state=random_state, dual="auto", verbose=0
        )
        if calibrate:
            models["LinearSVC"] = CalibratedClassifierCV(svc, cv=5)
        else:
            models["LinearSVC"] = svc

        # SGD Classifier
        models["SGDClassifier"] = SGDClassifier(
            loss="modified_huber",
            max_iter=1000,
            early_stopping=True,
            n_iter_no_change=20,
            random_state=random_state,
            verbose=0,
        )

        # Random Forest
        models["RandomForest"] = RandomForestClassifier(
            max_depth=5, n_estimators=500, random_state=random_state, n_jobs=-1
        )

        # XGBoost
        models["XGBoost"] = XGBClassifier(
            objective="multi:softmax",
            max_depth=5,
            n_estimators=500,
            random_state=random_state,
            verbosity=0,
        )

        return models

    @staticmethod
    def create_smote(
        sampling_strategy: str = "not majority",
        k_neighbors: int = 1,
        random_state: int = 42,
    ) -> SMOTE:
        """Create SMOTE oversampler for class imbalance.

        Args:
            sampling_strategy: Oversampling strategy ('not majority').
            k_neighbors: Number of nearest neighbors for SMOTE.
            random_state: Random state for reproducibility.

        Returns:
            SMOTE instance configured for oversampling.
        """
        return SMOTE(
            sampling_strategy=sampling_strategy,
            k_neighbors=k_neighbors,
            random_state=random_state,
        )


class ModelTrainer:
    """Trainer for multiclass classification models."""

    def __init__(
        self,
        X: np.ndarray,
        y: np.ndarray,
        models: Optional[Dict] = None,
        apply_smote: bool = True,
        random_state: int = 42,
        smote_sampling_strategy: str = "not majority",
        smote_k_neighbors: int = 1,
    ) -> None:
        """Initialize trainer with data and models.

        Args:
            X: Feature matrix (samples x features).
            y: Target labels as 1D array or string labels.
            models: Dictionary of model names to instances. If None, creates default set.
            apply_smote: Whether to apply SMOTE oversampling.
            random_state: Random state for reproducibility.
            smote_sampling_strategy: SMOTE sampling strategy ('not majority' or 'all').
            smote_k_neighbors: Number of neighbors for SMOTE.
        """
        assert isinstance(
            X, (np.ndarray, pd.DataFrame)
        ), f"X must be ndarray or DataFrame, got {type(X)}"
        assert X.ndim == 2, f"X must be 2D, got shape {X.shape}"
        assert len(X) > 0, f"X must have samples, got {len(X)}"
        assert isinstance(
            y, (np.ndarray, pd.Series)
        ), f"y must be ndarray or Series, got {type(y)}"
        assert y.ndim == 1, f"y must be 1D, got shape {y.shape}"
        assert len(X) == len(y), f"X and y must have same length: {len(X)} vs {len(y)}"
        assert isinstance(
            apply_smote, bool
        ), f"apply_smote must be bool, got {type(apply_smote)}"
        assert isinstance(
            random_state, int
        ), f"random_state must be int, got {type(random_state)}"

        self.X = X
        self.y = y
        self.models = (
            models
            if models is not None
            else ModelFactory.create_models(random_state=random_state)
        )
        self.apply_smote = apply_smote
        self.random_state = random_state
        self.smote_sampling_strategy = smote_sampling_strategy
        self.smote_k_neighbors = smote_k_neighbors
        self.label_encoder = LabelEncoder()
        self.trained_models = {}
        self.X_train = None
        self.y_train = None

    def prepare_data(self) -> Tuple[np.ndarray, np.ndarray]:
        """Prepare data: encode labels and optionally apply SMOTE.

        Returns:
            Tuple of (X_resampled, y_resampled_encoded) or (X, y_encoded) if no SMOTE.
        """
        y_encoded = self.label_encoder.fit_transform(self.y)
        assert isinstance(
            y_encoded, np.ndarray
        ), f"Expected numpy.ndarray, got {type(y_encoded)}"

        if self.apply_smote:
            smote = ModelFactory.create_smote(
                sampling_strategy=self.smote_sampling_strategy,
                k_neighbors=self.smote_k_neighbors,
                random_state=self.random_state,
            )
            X_resampled, y_resampled = smote.fit_resample(self.X, y_encoded)
            self.X_train = X_resampled
            self.y_train = y_resampled
            assert isinstance(
                self.y_train, np.ndarray
            ), f"Expected numpy.ndarray, got {type(self.y_train)}"
            return X_resampled, y_resampled
        else:
            self.X_train = self.X.values if isinstance(self.X, pd.DataFrame) else self.X
            self.y_train = y_encoded
            assert isinstance(
                self.y_train, np.ndarray
            ), f"Expected numpy.ndarray, got {type(self.y_train)}"
            return self.X_train, y_encoded

    def train_all_models(self) -> Dict:
        """Train all models on prepared data.

        Returns:
            Dictionary of model names to trained model instances.
        """
        X_train, y_train = self.prepare_data()

        for name, model in self.models.items():
            print(f"Training {name}...")
            model.fit(X_train, y_train)
            self.trained_models[name] = model

        return self.trained_models

    def predict(self, X_test: np.ndarray, model_name: str) -> np.ndarray:
        """Predict using a trained model.

        Args:
            X_test: Test feature matrix as numpy.ndarray.
            model_name: Name of the model to use for prediction.

        Returns:
            Predicted labels (encoded).
        """
        assert isinstance(
            X_test, np.ndarray
        ), f"Expected numpy.ndarray, got {type(X_test)}"
        if model_name not in self.trained_models:
            raise ValueError(f"Model '{model_name}' not trained. Train models first.")
        return self.trained_models[model_name].predict(X_test)

    def predict_proba(self, X_test: np.ndarray, model_name: str) -> np.ndarray:
        """Get prediction probabilities using a trained model.

        Args:
            X_test: Test feature matrix as numpy.ndarray.
            model_name: Name of the model to use.

        Returns:
            Probability matrix (n_samples x n_classes).
        """
        assert isinstance(
            X_test, np.ndarray
        ), f"Expected numpy.ndarray, got {type(X_test)}"
        if model_name not in self.trained_models:
            raise ValueError(f"Model '{model_name}' not trained. Train models first.")

        model = self.trained_models[model_name]
        if not hasattr(model, "predict_proba"):
            raise ValueError(f"Model '{model_name}' does not support predict_proba.")

        return model.predict_proba(X_test)

    def get_label_encoder(self) -> LabelEncoder:
        """Get the label encoder used during training.

        Returns:
            LabelEncoder instance with fitted classes.
        """
        return self.label_encoder

    def decode_predictions(self, y_pred_encoded: np.ndarray) -> np.ndarray:
        """Convert encoded predictions back to original labels.

        Args:
            y_pred_encoded: Encoded predictions.

        Returns:
            Decoded predictions with original labels.
        """
        return self.label_encoder.inverse_transform(y_pred_encoded)

    def save_models(self, output_dir: Path) -> Path:
        """Save trained models to disk.

        Args:
            output_dir: Directory to save models.

        Returns:
            Path to the saved models directory.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        encoder_path = output_dir / "label_encoder.pkl"
        with open(encoder_path, "wb") as f:
            pickle.dump(self.label_encoder, f)

        for model_name, model in self.trained_models.items():
            model_path = output_dir / f"{model_name}.pkl"
            with open(model_path, "wb") as f:
                pickle.dump(model, f)
            print(f"Model {model_name} saved to {model_path}")

        print(f"All models saved to {output_dir}")
        return output_dir

    def evaluate(self, X: np.ndarray, y: np.ndarray, eval_dir: Path = None,
                 cv: int = 5, eval_name: str = None) -> pd.DataFrame:
        """Evaluate all trained models using cross-validation.

        Args:
            X: Feature matrix (already feature-selected).
            y: True labels.
            eval_dir: Directory to save evaluation results and figures.
            cv: Number of cross-validation splits.
            eval_name: Name prefix for output files (e.g., 'full', 'fs', 'fs_de').

        Returns:
            DataFrame with evaluation metrics for all models.
        """
        if eval_dir is not None:
            eval_dir = Path(eval_dir)
            eval_dir.mkdir(parents=True, exist_ok=True)
            fig_dir = eval_dir.parent / "figures" / "model_evaluation"
            fig_dir.mkdir(parents=True, exist_ok=True)
        else:
            fig_dir = None

        prefix = f"{eval_name}_" if eval_name else ""

        X_array = X.values if isinstance(X, pd.DataFrame) else X
        y_array = y.values if isinstance(y, pd.Series) else y

        if not hasattr(self.label_encoder, 'classes_') or self.label_encoder.classes_ is None:
            self.label_encoder.fit(y_array)
        y_encoded = self.label_encoder.transform(y_array)
        smote = ModelFactory.create_smote(random_state=self.random_state)
        X_resampled, y_resampled = smote.fit_resample(X_array, y_encoded)

        all_fold_results = []
        summary_results = []
        sss = StratifiedShuffleSplit(n_splits=cv, test_size=0.2, random_state=self.random_state)

        for model_name, model in self.trained_models.items():
            print(f"Evaluating {model_name}...")
            fold_scores = []

            for fold_idx, (train_idx, test_idx) in enumerate(sss.split(X_resampled, y_resampled)):
                X_train_fold = X_resampled[train_idx]
                X_test_fold = X_resampled[test_idx]
                y_train_fold = y_resampled[train_idx]
                y_test_fold = y_resampled[test_idx]

                model_clone = clone(model)
                model_clone.fit(X_train_fold, y_train_fold)
                y_pred_encoded = model_clone.predict(X_test_fold)

                scores = make_score(y_test_fold, y_pred_encoded)
                fold_scores.append(scores)

                all_fold_results.append({
                    "model": model_name,
                    "fold": fold_idx,
                    **scores,
                })

                if fig_dir:
                    y_pred = self.label_encoder.inverse_transform(y_pred_encoded)
                    y_test = self.label_encoder.inverse_transform(y_test_fold)
                    self._save_confusion_matrix(
                        y_pred, y_test,
                        f"{prefix}{model_name}_fold_{fold_idx}",
                        fig_dir
                    )

            avg_scores = {k: np.mean([s[k] for s in fold_scores]) for k in fold_scores[0].keys()}
            std_scores = {k: np.std([s[k] for s in fold_scores]) for k in fold_scores[0].keys()}

            summary_results.append({
                "model": model_name,
                **{f"{k}_mean": v for k, v in avg_scores.items()},
                **{f"{k}_std": v for k, v in std_scores.items()},
            })

        fold_df = pd.DataFrame(all_fold_results)
        summary_df = pd.DataFrame(summary_results)

        if eval_dir:
            fold_csv = eval_dir / f"{prefix}model_evaluation_per_fold.csv"
            fold_df.to_csv(fold_csv, index=False)
            print(f"CSV saved: {fold_csv}")

            summary_csv = eval_dir / f"{prefix}model_evaluation.csv"
            summary_df.to_csv(summary_csv, index=False)
            print(f"CSV saved: {summary_csv}")

        return summary_df

    def _save_confusion_matrix(self, y_pred, y_test, name: str, output_dir: Path) -> None:
        """Save confusion matrix visualization.

        Args:
            y_pred: Predicted labels.
            y_test: True labels.
            name: Name for the plot.
            output_dir: Directory to save figure.
        """
        fig, ax = plt.subplots(figsize=(8, 8))

        ConfusionMatrixDisplay.from_predictions(
            y_test,
            y_pred,
            ax=ax,
            xticks_rotation="vertical",
            colorbar=False,
            normalize="true",
            values_format=".0%",
        )

        ax.set_title(f"Confusion Matrix {name}")
        plt.tight_layout()

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        save_path = output_dir / f"Confusion_Matrix_{name}.png"
        try:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            print(f"Figure saved: {save_path}")
        finally:
            plt.close()

    @classmethod
    def load_models(cls, model_dir: Path) -> "ModelTrainer":
        """Load trained models from disk.

        Args:
            model_dir: Directory containing saved models.

        Returns:
            ModelTrainer instance with loaded models.
        """
        model_dir = Path(model_dir)

        with open(model_dir / "label_encoder.pkl", "rb") as f:
            label_encoder = pickle.load(f)

        trained_models = {}
        for model_file in model_dir.glob("*.pkl"):
            if model_file.name == "label_encoder.pkl":
                continue
            model_name = model_file.stem
            with open(model_file, "rb") as f:
                trained_models[model_name] = pickle.load(f)

        trainer = cls.__new__(cls)
        trainer.trained_models = trained_models
        trainer.label_encoder = label_encoder
        trainer.X = None
        trainer.y = None
        print(f"Models loaded from {model_dir}")
        return trainer
