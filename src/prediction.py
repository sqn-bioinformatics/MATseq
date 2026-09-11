import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Optional
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.metrics import classification_report, confusion_matrix
from .model_training import ModelTrainer, ModelFactory, make_score
from .config import CLASS_ORDER

class ModelPredictor:
    def __init__(self, trainer: ModelTrainer):
        self.trainer = trainer
        self.predictions = {}
        self.probabilities = {}
        self.y_test = None

    def predict_samples(
        self, X_test: pd.DataFrame, y_test: pd.Series, sample_names: np.ndarray
    ) -> Dict[str, pd.DataFrame]:
        """Predict labels for test samples.
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

    def evaluate(self, output_dir: Path, subset: str = "train_ligands") -> pd.DataFrame:
        """Score predictions against true labels
        """
        rows = []
        for model_name, pred_df in self.predictions.items():
            y_pred = pred_df["prediction"].to_numpy()
            rows.append({"model": model_name, **make_score(self.y_test, y_pred)})

            present = set(self.y_test) | set(y_pred)
            order = CLASS_ORDER[subset]
            labels = [c for c in order if c in present] + [c for c in present if c not in order]
            report = classification_report(
                self.y_test, y_pred, labels=labels, zero_division=0, output_dict=True
            )
            pd.DataFrame(report).transpose().to_csv(
                output_dir / f"{model_name}_classification_report.csv"
            )

            cm = confusion_matrix(self.y_test, y_pred, labels=labels, normalize="true")
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

