"""Model training, evaluation, and prediction for multiclass classification."""

import json
import shutil
import tempfile
import numpy as np
import pandas as pd
import pickle
from pathlib import Path
from typing import Dict, Optional
from joblib import Memory
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from xgboost import XGBClassifier
from sklearn.pipeline import Pipeline as SkPipeline
from sklearn.utils.class_weight import compute_sample_weight
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    classification_report,
    confusion_matrix,
    ConfusionMatrixDisplay,
)
from sklearn.model_selection import (
    StratifiedKFold,
    GridSearchCV,
)
from sklearn.base import clone
import matplotlib.pyplot as plt

from .feature_engineering import create_feature_pipeline
from .config import FEATURE_SELECTION_CONFIG


def make_score(y_test, y_pred) -> dict:
    """Calculate multiple classification metrics.

    Args:
        y_test: True labels.
        y_pred: Predicted labels.

    Returns:
        dict: Dictionary of metric names to scores.
    """
    return {
        "accuracy": accuracy_score(y_test, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_test, y_pred),
        "precision": precision_score(y_test, y_pred, average="macro", zero_division=0),
        "recall": recall_score(y_test, y_pred, average="macro", zero_division=0),
        "f1": f1_score(y_test, y_pred, average="macro", zero_division=0),
        "f1_weighted": f1_score(y_test, y_pred, average="weighted", zero_division=0),
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
            max_iter=10000,
            tol=1e-3,
            random_state=random_state,
            dual="auto",
            verbose=0,
            class_weight="balanced",
        )
        if calibrate:
            models["LinearSVC"] = CalibratedClassifierCV(svc, cv=2)
        else:
            models["LinearSVC"] = svc

        models["SGDClassifier"] = SGDClassifier(
            loss="modified_huber",
            early_stopping=False,
            max_iter=1000,
            tol=1e-3,
            class_weight="balanced",
            random_state=random_state,
        )

        models["LogisticRegression"] = LogisticRegression(
            penalty="l2",
            class_weight="balanced",
            solver="liblinear",
            max_iter=10000,
            random_state=random_state,
        )

        models["RandomForest"] = RandomForestClassifier(
            max_depth=5,
            n_estimators=500,
            random_state=random_state,
            n_jobs=-1,
            class_weight="balanced",
        )

        models["XGBoost"] = XGBClassifier(
            objective="multi:softmax",
            max_depth=5,
            n_estimators=500,
            random_state=random_state,
            verbosity=0,
        )

        return models


class ModelTrainer:
    """Trainer for multiclass classification models."""

    def __init__(
        self,
        X: np.ndarray,
        y: np.ndarray,
        models: Optional[Dict] = None,
        random_state: int = 42,
    ) -> None:
        """Initialize trainer with data and models.

        Args:
            X: Feature matrix (samples x features).
            y: Target labels as 1D array or string labels.
            models: Dictionary of model names to instances. If None, creates default set.
            random_state: Random state for reproducibility.
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
        if not isinstance(random_state, int):
            raise TypeError(f"random_state must be int, got {type(random_state)}")

        self.X = X
        self.y = y
        self.models = (
            models
            if models is not None
            else ModelFactory.create_models(random_state=random_state)
        )
        self.random_state = random_state
        self.label_encoder = LabelEncoder()
        self.trained_models = {}

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

    def _build_tuning_pipeline(self, model, fold_seed: int) -> SkPipeline:
        fs = create_feature_pipeline(
            **FEATURE_SELECTION_CONFIG, random_state=fold_seed
        )
        # The FS steps are identical across all clf hyperparameter combos in a
        # fold; caching them means SelectKBest/ExtraTrees run once per fold, not
        # once per grid point.
        memory = getattr(self, "_memory", None)
        return SkPipeline([*fs.steps, ("clf", clone(model))], memory=memory)

    @staticmethod
    def _fit_kwargs(model_name: str, y) -> dict:
        if model_name == "XGBoost":
            return {"clf__sample_weight": compute_sample_weight("balanced", y)}
        return {}

    def tune_nested(
        self,
        X,
        y,
        param_grids: Dict[str, Dict],
        output_dir: Path,
        outer_cv: int = 5,
        inner_cv: int = 3,
        scoring: str = "f1_macro",
        eval_name: Optional[str] = None,
    ) -> pd.DataFrame:
        """Nested CV tuning + deployment refit.

        Leakage boundary: every step that learns from X or y (library-size
        normalization, log1p, feature selection, scaling, class-weight
        derivation, calibration) lives inside the sklearn Pipeline and is
        refit per outer fold. No selection or filtering on labels happens
        outside outer CV.

        Outer split estimates generalization of the tuning procedure. Inner
        GridSearchCV jointly selects FS and model hyperparameters per outer
        fold using f1_macro. Final deployment params are selected by majority
        vote across outer folds (tie-break: mean inner f1_macro, then pooled
        outer f1_macro); pipeline is refit on the full data.

        X and y are method arguments — not instance state.
        """
        output_dir = Path(output_dir)
        inner_results_dir = output_dir / "inner_cv_results"
        inner_results_dir.mkdir(parents=True, exist_ok=True)
        fig_dir = output_dir.parent / "figures" / "model_evaluation"
        fig_dir.mkdir(parents=True, exist_ok=True)

        memory_dir = tempfile.mkdtemp(prefix="matseq_fs_cache_")
        self._memory = Memory(location=memory_dir, verbose=0)

        prefix = f"{eval_name}_" if eval_name else ""

        y_array = y.values if isinstance(y, pd.Series) else y
        self.label_encoder.fit(y_array)
        y_encoded = self.label_encoder.transform(y_array)

        min_class_count = np.bincount(y_encoded).min()
        if outer_cv > min_class_count:
            raise ValueError(
                f"outer_cv={outer_cv} exceeds smallest class count={min_class_count}"
            )

        sample_ids = (
            X.index.to_numpy() if isinstance(X, pd.DataFrame) else np.arange(len(X))
        )

        outer = StratifiedKFold(
            n_splits=outer_cv, shuffle=True, random_state=self.random_state
        )

        per_fold_rows = []
        pooled: Dict[str, list] = {name: [] for name in self.models}

        fold_seeds = np.random.SeedSequence(self.random_state).generate_state(outer_cv)

        for fold_idx, (train_idx, test_idx) in enumerate(
            outer.split(np.arange(len(X)), y_encoded)
        ):
            if isinstance(X, pd.DataFrame):
                X_tr, X_te = X.iloc[train_idx], X.iloc[test_idx]
            else:
                X_tr, X_te = X[train_idx], X[test_idx]
            y_tr = y_encoded[train_idx]
            y_te = y_encoded[test_idx]
            fold_seed = int(fold_seeds[fold_idx])
            ids_te = sample_ids[test_idx]

            min_train_class_count = np.bincount(y_tr).min()
            if inner_cv > min_train_class_count:
                raise ValueError(
                    f"inner_cv={inner_cv} exceeds smallest training class count={min_train_class_count} in outer fold {fold_idx}"
                )

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
                gs.fit(X_tr, y_tr, **self._fit_kwargs(model_name, y_tr))

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

                pooled[model_name].extend(
                    {"sample_id": sid, "true_enc": t, "pred_enc": p, "outer_fold": fold_idx}
                    for sid, t, p in zip(ids_te, y_te.tolist(), y_pred_enc.tolist())
                )

        per_fold_df = pd.DataFrame(per_fold_rows)
        per_fold_df.to_csv(output_dir / "nested_cv_per_fold.csv", index=False)

        metric_cols = [
            "accuracy",
            "balanced_accuracy",
            "precision",
            "recall",
            "f1",
            "f1_weighted",
            "inner_best_score",
        ]
        summary = per_fold_df.groupby("model")[metric_cols].agg(["mean", "std"])
        summary.columns = [f"{m}_{s}" for m, s in summary.columns]
        summary = summary.reset_index()

        oof_frames = []
        pooled_rows = []
        class_names = list(self.label_encoder.classes_)
        for model_name in self.models:
            df = pd.DataFrame(pooled[model_name])
            y_true_dec = self.decode_predictions(df["true_enc"].to_numpy())
            y_pred_dec = self.decode_predictions(df["pred_enc"].to_numpy())

            pooled_scores = make_score(df["true_enc"], df["pred_enc"])
            pooled_rows.append(
                {"model": model_name, **{f"pooled_{k}": v for k, v in pooled_scores.items()}}
            )

            oof_frames.append(
                pd.DataFrame(
                    {
                        "sample_id": df["sample_id"],
                        "true_label": y_true_dec,
                        "pred_label": y_pred_dec,
                        "outer_fold": df["outer_fold"],
                        "model": model_name,
                        "subset": eval_name or "",
                    }
                )
            )

            report = classification_report(
                y_true_dec, y_pred_dec, labels=class_names,
                zero_division=0, output_dict=True,
            )
            pd.DataFrame(report).transpose().to_csv(
                output_dir / f"{prefix}{model_name}_classification_report.csv"
            )

            cm = confusion_matrix(y_true_dec, y_pred_dec, labels=class_names)
            cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True).clip(min=1)
            pd.DataFrame(cm, index=class_names, columns=class_names).to_csv(
                fig_dir / f"{prefix}{model_name}_confusion_matrix.csv"
            )
            pd.DataFrame(cm_norm, index=class_names, columns=class_names).to_csv(
                fig_dir / f"{prefix}{model_name}_confusion_matrix_normalized.csv"
            )
            self._save_confusion_matrix(
                cm_norm, class_names, f"{prefix}{model_name}", fig_dir
            )

        pd.concat(oof_frames, ignore_index=True).to_csv(
            output_dir / f"{prefix}oof_predictions.csv", index=False
        )

        pooled_df = pd.DataFrame(pooled_rows)
        summary = summary.merge(pooled_df, on="model")
        summary.to_csv(output_dir / "nested_cv_summary.csv", index=False)
        print(f"Nested CV summary saved: {output_dir / 'nested_cv_summary.csv'}")

        selected_params = self._select_and_refit(
            X, y_encoded, per_fold_df
        )
        with open(output_dir / "selected_params.json", "w") as f:
            json.dump(selected_params, f, indent=2, default=str)
        print(f"Selected params saved: {output_dir / 'selected_params.json'}")

        self._memory = None
        shutil.rmtree(memory_dir, ignore_errors=True)

        return summary

    def _select_and_refit(
        self, X, y_encoded, per_fold_df: pd.DataFrame
    ) -> Dict[str, dict]:
        """Select deployment params by majority vote across outer folds; tie-break by mean inner f1_macro, then mean outer f1; refit on full data."""
        selected = {}
        for model_name, model in self.models.items():
            rows = per_fold_df[per_fold_df["model"] == model_name]
            ranking = (
                rows.groupby("best_params")
                .agg(
                    count=("inner_best_score", "size"),
                    mean_inner=("inner_best_score", "mean"),
                    mean_f1=("f1", "mean"),
                )
                .sort_values(
                    ["count", "mean_inner", "mean_f1"], ascending=False
                )
            )
            chosen_json = ranking.index[0]
            chosen = json.loads(chosen_json)
            selected[model_name] = {
                "chosen_params": chosen,
                "selection_rule": "majority_vote_then_mean_inner_f1_then_outer_f1",
                "outer_fold_params": [
                    json.loads(p) for p in rows["best_params"].tolist()
                ],
            }

            pipe = self._build_tuning_pipeline(model, self.random_state)
            pipe.set_params(**chosen)
            print(f"  Deploying {model_name} with {chosen}...")
            pipe.fit(X, y_encoded, **self._fit_kwargs(model_name, y_encoded))
            self.trained_models[model_name] = pipe

        return selected

    def _save_confusion_matrix(
        self, cm_norm, class_names, name: str, output_dir: Path
    ) -> None:
        fig, ax = plt.subplots(figsize=(8, 8))
        ConfusionMatrixDisplay(cm_norm, display_labels=class_names).plot(
            ax=ax,
            xticks_rotation="vertical",
            colorbar=False,
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
