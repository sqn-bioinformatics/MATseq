from pathlib import Path

import pandas as pd

from .config import CLASS_ORDER, SUBSET_DISPLAY_NAMES
from .model_training import ModelTrainer, evaluate
from .visualization import plot_probability_heatmap


def predict_samples(trainer: ModelTrainer, X: pd.DataFrame, y: pd.Series, subset: str,
                    output_dir: Path, all_controls: bool) -> pd.DataFrame:
    """Predict X with every trained model, save predictions/probabilities/heatmaps, score against y."""
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for model_name, model in trainer.trained_models.items():
        y_pred = trainer.label_encoder.inverse_transform(model.predict(X))
        pd.DataFrame({"sample": X.index, "prediction": y_pred}).to_csv(
            output_dir / f"{model_name}_predictions.csv", index=False
        )
        proba = pd.DataFrame(
            model.predict_proba(X), columns=trainer.label_encoder.classes_, index=X.index
        )
        proba.to_csv(output_dir / f"{model_name}_probabilities.csv")
        plot_probability_heatmap(
            proba, CLASS_ORDER[subset],
            title=f"{model_name} Prediction Probabilities {SUBSET_DISPLAY_NAMES[subset]}",
            true_labels=y, all_controls=all_controls,
            output_path=output_dir,
            output_filename=f"{model_name}_probabilities_heatmap.png",
        )
        rows.append({"model": model_name, **evaluate(y, y_pred, model_name, subset,
                                                     output_dir, output_dir)})
    summary = pd.DataFrame(rows)
    summary.to_csv(output_dir / "test_scores_summary.csv", index=False)
    return summary
