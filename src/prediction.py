"""Model prediction and evaluation on new data and multiple runs."""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Optional
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.metrics import classification_report, confusion_matrix

from .model_training import ModelTrainer, ModelFactory, make_score
from .config import SUBSET_CLASS_ORDERS, order_labels


class ModelPredictor:
    """Make predictions on new samples with trained models."""

    def __init__(self, trainer: ModelTrainer):
        """Initialize predictor with trained models.

        Args:
            trainer: Trained ModelTrainer instance with fitted models.
        """
        if not hasattr(trainer, "trained_models"):
            raise AttributeError("trainer must have trained_models attribute")
        if len(trainer.trained_models) == 0:
            raise ValueError("trainer must have at least one trained model")
        self.trainer = trainer
        self.predictions = {}
        self.probabilities = {}
        self.true_labels = None

    def predict_samples(
        self, X_test: pd.DataFrame, sample_names: np.ndarray = None, y_test=None
    ) -> Dict[str, pd.DataFrame]:
        """Predict labels for test samples using all trained models.

        Args:
            X_test: Test feature matrix.
            sample_names: Optional sample identifiers.
            y_test: Optional true labels for samples.

        Returns:
            Dictionary mapping model names to prediction DataFrames.
        """
        predictions_dict = {}

        if y_test is not None:
            self.true_labels = y_test

        for model_name in self.trainer.trained_models.keys():
            model = self.trainer.trained_models[model_name]
            y_pred_encoded = model.predict(X_test)
            y_pred = self.trainer.decode_predictions(y_pred_encoded)

            pred_df = pd.DataFrame(
                {
                    "sample": (
                        sample_names
                        if sample_names is not None
                        else np.arange(len(y_pred))
                    ),
                    "prediction": y_pred,
                }
            )

            predictions_dict[model_name] = pred_df
            self.predictions[model_name] = pred_df

            proba = model.predict_proba(X_test)
            proba_df = pd.DataFrame(
                proba,
                columns=self.trainer.label_encoder.classes_,
                index=(sample_names),
            )
            self.probabilities[model_name] = proba_df

        return predictions_dict

    def save_predictions(self, output_dir: Path):
        """Save predictions to CSV files.

        Args:
            output_dir: Directory to save predictions.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        for model_name, pred_df in self.predictions.items():
            output_path = output_dir / f"{model_name}_predictions.csv"
            pred_df.to_csv(output_path, index=False)
            print(f"Saved predictions to {output_path}")

        for model_name, proba_df in self.probabilities.items():
            output_path = output_dir / f"{model_name}_probabilities.csv"
            proba_df.to_csv(output_path)
            print(f"Saved probabilities to {output_path}")

    def evaluate(self, output_dir: Path, subset: str = "main_ligands") -> pd.DataFrame:
        """Score predictions against true labels, identically to nested-CV training.

        Reuses make_score (accuracy, balanced accuracy, macro precision/recall/f1,
        weighted f1) and writes a per-model classification report and confusion
        matrices, mirroring ModelTrainer.tune_nested.

        Args:
            output_dir: Directory to save score files.

        Returns:
            DataFrame with one row of metrics per model.
        """
        if self.true_labels is None:
            raise ValueError("true labels required; call predict_samples with y_test")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        y_true = np.asarray(self.true_labels)

        rows = []
        for model_name, pred_df in self.predictions.items():
            y_pred = pred_df["prediction"].to_numpy()
            rows.append({"model": model_name, **make_score(y_true, y_pred)})

            labels = order_labels(set(y_true) | set(y_pred), subset)
            report = classification_report(
                y_true, y_pred, labels=labels, zero_division=0, output_dict=True
            )
            pd.DataFrame(report).transpose().to_csv(
                output_dir / f"{model_name}_classification_report.csv"
            )

            cm = confusion_matrix(y_true, y_pred, labels=labels)
            cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True).clip(min=1)
            pd.DataFrame(cm, index=labels, columns=labels).to_csv(
                output_dir / f"{model_name}_confusion_matrix.csv"
            )
            pd.DataFrame(cm_norm, index=labels, columns=labels).to_csv(
                output_dir / f"{model_name}_confusion_matrix_normalized.csv"
            )
            self.trainer._save_confusion_matrix(
                cm_norm, labels, model_name, output_dir, subset=subset,
            )

        summary = pd.DataFrame(rows)
        summary.to_csv(output_dir / "test_scores_summary.csv", index=False)
        print(f"Saved scores to {output_dir / 'test_scores_summary.csv'}")
        return summary

    def create_probability_heatmaps(
        self, output_dir: Path, subset: str = "main_ligands", all_controls: bool = False
    ):
        """Create heatmap visualizations of prediction probabilities.

        Args:
            output_dir: Directory to save figures.
            subset: Key into SUBSET_CLASS_ORDERS ("main_ligands", "additional_ligands", "bacteria_ligands").
            all_controls: If True, plot every negative_control and LPS sample instead of one each.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        class_order = SUBSET_CLASS_ORDERS.get(subset, SUBSET_CLASS_ORDERS["main_ligands"])

        for model_name, proba_df in self.probabilities.items():
            proba_df_ordered = proba_df.copy()

            available_classes = [c for c in class_order if c in proba_df.columns]
            proba_df_ordered = proba_df_ordered[available_classes]

            if self.true_labels is not None:
                rng = np.random.default_rng(42)
                nc_idx = self.true_labels[self.true_labels == "negative_control"].index
                lps_idx = self.true_labels[self.true_labels == "LPS"].index
                if all_controls:
                    rand_nc = list(nc_idx)
                    rand_lps = list(lps_idx)
                else:
                    rand_nc = list(rng.choice(nc_idx, size=1, replace=False))
                    rand_lps = list(rng.choice(lps_idx, size=1, replace=False))
                remaining = [
                    i
                    for cls in class_order
                    if cls not in ("negative_control", "LPS")
                    for i in self.true_labels[self.true_labels == cls].index
                ]
                ordered_idx = rand_nc + rand_lps + remaining
                proba_df_ordered = proba_df_ordered.loc[ordered_idx]
                proba_df_ordered.index = self.true_labels.loc[ordered_idx].values

            plt.figure(figsize=(12, 8))
            sns.heatmap(
                proba_df_ordered, cmap="YlGnBu", cbar_kws={"label": "Predicted class probability"}
            )
            plt.title(f"{model_name} - Prediction Probabilities")
            plt.xlabel("Class")
            plt.ylabel("Sample")
            plt.tight_layout()

            output_path = output_dir / f"{model_name}_probabilities_heatmap.png"
            try:
                plt.savefig(output_path, dpi=300, bbox_inches="tight")
                print(f"Saved heatmap to {output_path}")
            finally:
                plt.close()

