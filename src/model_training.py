"""Model training, evaluation, and prediction for multiclass classification."""

import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
)
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder
from sklearn.svm import LinearSVC
from sklearn.utils.class_weight import compute_sample_weight
from xgboost import XGBClassifier

from .config import CLASS_ORDER, FEATURE_SELECTION_CONFIG, SUBSET_DISPLAY_NAMES
from .feature_engineering import ColumnSelector, feature_pipeline, preprocessing_pipeline
from .visualization import plot_confusion_matrix


def make_score(y_true, y_pred) -> dict[str, float]:
    """Accuracy, balanced accuracy, macro precision/recall/f1, and weighted f1."""
    return {
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, average="macro", zero_division=0),
        "recall": recall_score(y_true, y_pred, average="macro", zero_division=0),
        "f1": f1_score(y_true, y_pred, average="macro", zero_division=0),
        "f1_weighted": f1_score(y_true, y_pred, average="weighted", zero_division=0),
    }


def evaluate(y_true, y_pred, model_name: str, subset: str, output_dir: Path,
             fig_dir: Path) -> dict[str, float]:
    """Write the classification report and confusion matrix (CSV + PNG); return scores."""
    present = set(y_true) | set(y_pred)
    order = CLASS_ORDER[subset]
    labels = [c for c in order if c in present] + sorted(present - set(order))
    report = classification_report(
        y_true, y_pred, labels=labels, zero_division=0, output_dict=True
    )
    pd.DataFrame(report).transpose().to_csv(
        output_dir / f"{model_name}_classification_report.csv"
    )
    cm = confusion_matrix(y_true, y_pred, labels=labels, normalize="true")
    pd.DataFrame(cm, index=labels, columns=labels).to_csv(
        output_dir / f"{model_name}_confusion_matrix.csv"
    )
    plot_confusion_matrix(
        cm, labels,
        title=f"{model_name} Confusion Matrix {SUBSET_DISPLAY_NAMES[subset]}",
        output_path=fig_dir,
        output_filename=f"Confusion_Matrix_{model_name}.png",
    )
    return make_score(y_true, y_pred)


class ModelTrainer:
    """Nested-CV tuning and per-gene-set refitting of the five classifiers."""

    def __init__(self, X: pd.DataFrame, y: pd.Series, random_state: int = 42) -> None:
        self.X = X
        self.y = y
        self.random_state = random_state
        self.label_encoder = LabelEncoder()
        self.y_enc = self.label_encoder.fit_transform(y)
        self.trained_models: dict[str, Pipeline] = {}
        self.models = {
            "LinearSVC": CalibratedClassifierCV(
                LinearSVC(
                    max_iter=10000,
                    tol=1e-3,
                    random_state=random_state,
                    dual="auto",
                    verbose=0,
                    class_weight="balanced",
                ),
                cv=2,
            ),
            "SGDClassifier": SGDClassifier(
                loss="modified_huber",
                early_stopping=False,
                max_iter=1000,
                tol=1e-3,
                class_weight="balanced",
                random_state=random_state,
            ),
            "LogisticRegression": LogisticRegression(
                penalty="l2",
                class_weight="balanced",
                solver="liblinear",
                max_iter=10000,
                random_state=random_state,
            ),
            "RandomForest": RandomForestClassifier(
                max_depth=5,
                n_estimators=500,
                random_state=random_state,
                n_jobs=-1,
                class_weight="balanced",
            ),
            "XGBoost": XGBClassifier(
                objective="multi:softmax",
                max_depth=5,
                n_estimators=500,
                random_state=random_state,
                verbosity=0,
            ),
        }

    def save_models(self, output_dir: Path) -> None:
        """Pickle the trained models and the label encoder into output_dir."""
        output_dir.mkdir(parents=True, exist_ok=True)
        with open(output_dir / "label_encoder.pkl", "wb") as f:
            pickle.dump(self.label_encoder, f)
        for model_name, model in self.trained_models.items():
            with open(output_dir / f"{model_name}.pkl", "wb") as f:
                pickle.dump(model, f)
        print(f"All models saved to {output_dir}")

    def tune_nested(
        self,
        param_grids: dict[str, dict[str, list]],
        output_dir: Path,
        fig_dir: Path,
        outer_cv: int = 5,
        inner_cv: int = 3,
        scoring: str = "f1_macro",
    ) -> pd.DataFrame:
        """Nested CV tuning; selects per-model params by majority vote across outer folds."""
        inner_results_dir = output_dir / "inner_cv_results"
        inner_results_dir.mkdir(parents=True, exist_ok=True)

        outer = StratifiedKFold(
            n_splits=outer_cv, shuffle=True, random_state=self.random_state
        )
        fold_seeds = np.random.SeedSequence(self.random_state).generate_state(outer_cv)
        per_fold_rows = []
        oof_frames = []

        for fold_idx, (train_idx, test_idx) in enumerate(outer.split(self.X, self.y_enc)):
            X_tr, X_te = self.X.iloc[train_idx], self.X.iloc[test_idx]
            y_tr, y_te = self.y_enc[train_idx], self.y.to_numpy()[test_idx]
            fold_seed = int(fold_seeds[fold_idx])

            for model_name, model in self.models.items():
                pipe = Pipeline([
                    *feature_pipeline(
                        **FEATURE_SELECTION_CONFIG, random_state=fold_seed
                    ).steps,
                    ("clf", clone(model)),
                ])
                gs = GridSearchCV(
                    pipe,
                    param_grids[model_name],
                    cv=StratifiedKFold(n_splits=inner_cv, shuffle=True, random_state=fold_seed),
                    scoring=scoring,
                    n_jobs=-1,
                    refit=True,
                )
                print(f"  Outer fold {fold_idx} — tuning {model_name}...")
                gs.fit(X_tr, y_tr, **(
                    {"clf__sample_weight": compute_sample_weight("balanced", y_tr)}
                    if model_name == "XGBoost" else {}
                ))
                y_pred = self.label_encoder.inverse_transform(gs.best_estimator_.predict(X_te))

                per_fold_rows.append({
                    "model": model_name,
                    "outer_fold": fold_idx,
                    "best_params": json.dumps(gs.best_params_, default=str),
                    "inner_best_score": gs.best_score_,
                    **make_score(y_te, y_pred),
                })
                pd.DataFrame(gs.cv_results_).to_csv(
                    inner_results_dir / f"{model_name}_fold_{fold_idx}.csv", index=False
                )
                oof_frames.append(pd.DataFrame({
                    "sample_id": X_te.index,
                    "true_label": y_te,
                    "pred_label": y_pred,
                    "outer_fold": fold_idx,
                    "model": model_name,
                }))

        per_fold = pd.DataFrame(per_fold_rows)
        per_fold.to_csv(output_dir / "nested_cv_per_fold.csv", index=False)
        oof = pd.concat(oof_frames, ignore_index=True)
        oof.to_csv(output_dir / "oof_predictions.csv", index=False)

        pooled = pd.DataFrame({
            m: evaluate(d.true_label, d.pred_label, m, "train_ligands", output_dir, fig_dir)
            for m, d in oof.groupby("model")
        }).T.add_prefix("pooled_")
        metric_cols = ["accuracy", "balanced_accuracy", "precision", "recall", "f1",
                       "f1_weighted", "inner_best_score"]
        summary = per_fold.groupby("model")[metric_cols].agg(["mean", "std"])
        summary.columns = [f"{m}_{s}" for m, s in summary.columns]
        summary = summary.join(pooled).reset_index()

        # Tie-break the majority vote by mean inner f1_macro, then mean outer f1.
        self.selected_params_ = {}
        for model_name, rows in per_fold.groupby("model", sort=False):
            ranking = (
                rows.groupby("best_params")
                .agg(
                    count=("inner_best_score", "size"),
                    mean_inner=("inner_best_score", "mean"),
                    mean_f1=("f1", "mean"),
                )
                .sort_values(["count", "mean_inner", "mean_f1"], ascending=False)
            )
            self.selected_params_[model_name] = {
                "chosen_params": json.loads(ranking.index[0]),
                "selection_rule": "majority_vote_then_mean_inner_f1_then_outer_f1",
                "outer_fold_params": [json.loads(p) for p in rows["best_params"]],
            }
        with open(output_dir / "selected_params.json", "w") as f:
            json.dump(self.selected_params_, f, indent=2, default=str)
        return summary

    def refit(self, genes: list[str]) -> None:
        """Refit every model with its selected params on the full panel restricted to genes."""
        for model_name, model in self.models.items():
            pipe = Pipeline([
                *preprocessing_pipeline().steps,
                ("select_genes", ColumnSelector(genes)),
                ("clf", clone(model)),
            ]).set_output(transform="pandas")
            pipe.set_params(**self.selected_params_[model_name]["chosen_params"])
            pipe.fit(self.X, self.y_enc, **(
                {"clf__sample_weight": compute_sample_weight("balanced", self.y_enc)}
                if model_name == "XGBoost" else {}
            ))
            self.trained_models[model_name] = pipe
