"""Model training, evaluation, and prediction for multiclass classification."""

import json
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
from imblearn.pipeline import Pipeline as ImbPipeline
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    ConfusionMatrixDisplay,
)
from sklearn.model_selection import (
    StratifiedShuffleSplit,
    StratifiedKFold,
    GridSearchCV,
)
from sklearn.base import clone
import matplotlib.pyplot as plt

from .feature_engineering import create_feature_pipeline
from .config import FEATURE_SELECTION_CONFIG

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

        models["SGDClassifier"] = SGDClassifier(
            loss="modified_huber",
            max_iter=1000,
            early_stopping=True,
            n_iter_no_change=20,
            eta0=0.01,
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
        if not isinstance(X, (np.ndarray, pd.DataFrame)):
            raise TypeError(f"X must be ndarray or DataFrame, got {type(X)}")
        if X.ndim != 2:
            raise ValueError(f"X must be 2D, got shape {X.shape}")
        if len(X) == 0:
            raise ValueError(f"X must have samples, got {len(X)}")
        if not isinstance(y, (np.ndarray, pd.Series)):
            raise TypeError(f"y must be ndarray or Series, got {type(y)}")
        if y.ndim != 1:
            raise ValueError(f"y must be 1D, got shape {y.shape}")
        if len(X) != len(y):
            raise ValueError(f"X and y must have same length: {len(X)} vs {len(y)}")
        if not isinstance(apply_smote, bool):
            raise TypeError(f"apply_smote must be bool, got {type(apply_smote)}")
        if not isinstance(random_state, int):
            raise TypeError(f"random_state must be int, got {type(random_state)}")

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

        if self.apply_smote:
            smote = ModelFactory.create_smote(
                sampling_strategy=self.smote_sampling_strategy,
                k_neighbors=self.smote_k_neighbors,
                random_state=self.random_state,
            )
            X_resampled, y_resampled = smote.fit_resample(self.X, y_encoded)
            self.X_train = X_resampled
            self.y_train = y_resampled
            return X_resampled, y_resampled
        else:
            self.X_train = self.X.values if isinstance(self.X, pd.DataFrame) else self.X
            self.y_train = y_encoded
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
        if not isinstance(X_test, np.ndarray):
            raise TypeError(f"X_test must be ndarray, got {type(X_test)}")
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
        if not isinstance(X_test, np.ndarray):
            raise TypeError(f"X_test must be ndarray, got {type(X_test)}")
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

    def _build_tuning_pipeline(self, model, fold_seed: int) -> ImbPipeline:
        fs = create_feature_pipeline(
            **FEATURE_SELECTION_CONFIG, random_state=fold_seed
        )
        smote = SMOTE(
            sampling_strategy=self.smote_sampling_strategy,
            k_neighbors=self.smote_k_neighbors,
            random_state=fold_seed,
        )
        return ImbPipeline([*fs.steps, ("smote", smote), ("clf", clone(model))])

    def tune_nested(
        self,
        X,
        y,
        param_grids: Dict[str, Dict],
        output_dir: Path,
        outer_cv: int = 10,
        inner_cv: int = 5,
        scoring: str = "f1_macro",
        eval_name: Optional[str] = None,
    ) -> pd.DataFrame:
        """Nested CV tuning + deployment refit.

        Outer split estimates generalization of the tuning procedure.
        Inner GridSearchCV (with SMOTE inside via imblearn pipeline) jointly
        selects FS and model hyperparameters per outer fold. A final refit on
        the full data populates self.trained_models with deployed pipelines.
        """
        output_dir = Path(output_dir)
        inner_results_dir = output_dir / "inner_cv_results"
        inner_results_dir.mkdir(parents=True, exist_ok=True)
        fig_dir = output_dir.parent / "figures" / "model_evaluation"
        fig_dir.mkdir(parents=True, exist_ok=True)

        prefix = f"{eval_name}_" if eval_name else ""

        y_array = y.values if isinstance(y, pd.Series) else y
        self.label_encoder.fit(y_array)
        y_encoded = self.label_encoder.transform(y_array)

        outer = StratifiedShuffleSplit(
            n_splits=outer_cv, test_size=0.2, random_state=self.random_state
        )

        per_fold_rows = []

        for fold_idx, (train_idx, test_idx) in enumerate(
            outer.split(np.arange(len(X)), y_encoded)
        ):
            if isinstance(X, pd.DataFrame):
                X_tr, X_te = X.iloc[train_idx], X.iloc[test_idx]
            else:
                X_tr, X_te = X[train_idx], X[test_idx]
            y_tr = y_encoded[train_idx]
            y_te = y_encoded[test_idx]
            fold_seed = self.random_state + fold_idx

            for model_name, model in self.models.items():
                pipe = self._build_tuning_pipeline(model, fold_seed)
                inner = StratifiedKFold(
                    n_splits=inner_cv, shuffle=True, random_state=fold_seed
                )
                gs = GridSearchCV(
                    pipe,
                    param_grids[model_name],
                    cv=inner,
                    scoring=scoring,
                    n_jobs=-1,
                    refit=True,
                )
                print(f"  Outer fold {fold_idx} — tuning {model_name}...")
                gs.fit(X_tr, y_tr)

                y_pred_enc = gs.best_estimator_.predict(X_te)
                scores = make_score(y_te, y_pred_enc)

                per_fold_rows.append(
                    {
                        "model": model_name,
                        "outer_fold": fold_idx,
                        "best_params": json.dumps(gs.best_params_, default=str),
                        "inner_best_score": gs.best_score_,
                        **scores,
                    }
                )

                pd.DataFrame(gs.cv_results_).to_csv(
                    inner_results_dir / f"{model_name}_fold_{fold_idx}.csv",
                    index=False,
                )

                y_pred_dec = self.label_encoder.inverse_transform(y_pred_enc)
                y_te_dec = self.label_encoder.inverse_transform(y_te)
                self._save_confusion_matrix(
                    y_pred_dec,
                    y_te_dec,
                    f"{prefix}{model_name}_fold_{fold_idx}",
                    fig_dir,
                )

        per_fold_df = pd.DataFrame(per_fold_rows)
        per_fold_df.to_csv(output_dir / "nested_cv_per_fold.csv", index=False)

        metric_cols = [
            "accuracy",
            "precision",
            "recall",
            "f1",
            "roc_auc",
            "inner_best_score",
        ]
        summary = per_fold_df.groupby("model")[metric_cols].agg(["mean", "std"])
        summary.columns = [f"{m}_{s}" for m, s in summary.columns]
        summary = summary.reset_index()
        summary.to_csv(output_dir / "nested_cv_summary.csv", index=False)
        print(f"Nested CV summary saved: {output_dir / 'nested_cv_summary.csv'}")

        selected_params = {}
        for model_name, model in self.models.items():
            pipe = self._build_tuning_pipeline(model, self.random_state)
            inner = StratifiedKFold(
                n_splits=inner_cv, shuffle=True, random_state=self.random_state
            )
            gs = GridSearchCV(
                pipe,
                param_grids[model_name],
                cv=inner,
                scoring=scoring,
                n_jobs=-1,
                refit=True,
            )
            print(f"  Final refit — tuning {model_name} on full data...")
            gs.fit(X, y_encoded)
            self.trained_models[model_name] = gs.best_estimator_
            selected_params[model_name] = gs.best_params_

        with open(output_dir / "selected_params.json", "w") as f:
            json.dump(selected_params, f, indent=2, default=str)
        print(f"Selected params saved: {output_dir / 'selected_params.json'}")

        return summary

    def _save_confusion_matrix(
        self, y_pred, y_test, name: str, output_dir: Path
    ) -> None:
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
