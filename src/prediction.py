"""Model prediction and evaluation on new data and multiple runs."""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Optional
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.metrics import classification_report, confusion_matrix

from .model_training import ModelTrainer, ModelFactory, make_score
from .config import CLASS_ORDER
from .visualization import plot_probability_heatmap, subset_display, order_labels


class ModelPredictor:
    def __init__(self, trainer: ModelTrainer):
        if not hasattr(trainer, "trained_models"):
            raise AttributeError("trainer must have trained_models attribute")
        if len(trainer.trained_models) == 0:
            raise ValueError("trainer must have at least one trained model")
        self.trainer = trainer
        self.predictions = {}
        self.probabilities = {}
        self.y_test = None

    def predict_samples(
        self, X_test: pd.DataFrame, y_test: np.ndarray, sample_names: np.ndarray
    ) -> Dict[str, pd.DataFrame]:
        """Predict labels for test samples using all trained models.
        """
        predictions_dict = {}
        self.y_test = y_test

        for model_name in self.trainer.trained_models.keys():
            model = self.trainer.trained_models[model_name]
            y_pred_encoded = model.predict(X_test)
            y_pred = self.trainer.decode_predictions(y_pred_encoded)

            pred_df = pd.DataFrame(
                {
                    "sample": sample_names,
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
        """Score predictions against true labels
        """
        if self.y_test is None:
            raise ValueError("True labels required; call predict_samples")

        rows = []
        for model_name, pred_df in self.predictions.items():
            y_pred = pred_df["prediction"].to_numpy()
            rows.append({"model": model_name, **make_score(y_test, y_pred)})

            labels = order_labels(set(y_test) | set(y_pred), subset)
            report = classification_report(
                y_test, y_pred, labels=labels, zero_division=0, output_dict=True
            )
            pd.DataFrame(report).transpose().to_csv(
                output_dir / f"{model_name}_classification_report.csv"
            )

            cm = confusion_matrix(y_test, y_pred, labels=labels, normalize="true")
            pd.DataFrame(cm, index=labels, columns=labels).to_csv(
                output_dir / f"{model_name}_confusion_matrix.csv"
            )
            self.trainer._save_confusion_matrix(
                cm, labels, model_name, output_dir, subset=subset,
            )

        summary = pd.DataFrame(rows)
        summary.to_csv(output_dir / "test_scores_summary.csv", index=False)
        print(f"Saved scores to {output_dir / 'test_scores_summary.csv'}")
        return summary

    def create_probability_heatmaps(
        self, output_dir: Path, subset: str = "main_ligands", all_controls: bool = False
    ):
        """Create heatmap visualizations of prediction probabilities.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        class_order = CLASS_ORDER.get(subset, CLASS_ORDER["main_ligands"])

        for model_name, proba_df in self.probabilities.items():
            output_path = plot_probability_heatmap(
                proba_df, class_order,
                title=f"{model_name} Prediction Probabilities {subset_display(subset)}",
                true_labels=self.true_labels, all_controls=all_controls,
                output_dir=output_dir,
                filename=f"{model_name}_probabilities_heatmap.png",
            )
            print(f"Saved heatmap to {output_path}")

