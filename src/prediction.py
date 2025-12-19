"""Model prediction and evaluation on new data and multiple runs."""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Tuple, Optional
import matplotlib.pyplot as plt
import seaborn as sns

from .model_training import ModelFactory, ModelTrainer
from .evaluation import make_score
from .feature_engineering import create_feature_pipeline
from .utils import save_csv


class ModelPredictor:
    """Make predictions on new samples with trained models."""

    def __init__(self, trainer: ModelTrainer):
        """Initialize predictor with trained models.

        Args:
            trainer: Trained ModelTrainer instance with fitted models.
        """
        self.trainer = trainer
        self.predictions = {}
        self.probabilities = {}

    def predict_samples(
        self, X_test: pd.DataFrame, sample_names: np.ndarray = None
    ) -> Dict[str, pd.DataFrame]:
        """Predict labels for test samples using all trained models.

        Args:
            X_test: Test feature matrix.
            sample_names: Optional sample identifiers.

        Returns:
            Dictionary mapping model names to prediction DataFrames.
        """
        predictions_dict = {}

        for model_name in self.trainer.trained_models.keys():
            y_pred_encoded = self.trainer.predict(X_test, model_name)
            y_pred = self.trainer.decode_predictions(y_pred_encoded)

            # Create prediction DataFrame
            pred_df = pd.DataFrame(
                {"sample": sample_names if sample_names is not None else np.arange(len(y_pred)),
                 "prediction": y_pred}
            )

            predictions_dict[model_name] = pred_df
            self.predictions[model_name] = pred_df

            # Get probabilities if available
            if hasattr(self.trainer.trained_models[model_name], "predict_proba"):
                proba = self.trainer.predict_proba(X_test, model_name)
                proba_df = pd.DataFrame(
                    proba,
                    columns=self.trainer.label_encoder.classes_,
                    index=sample_names if sample_names is not None else np.arange(len(y_pred)),
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

    def create_probability_heatmaps(self, output_dir: Path):
        """Create heatmap visualizations of prediction probabilities.

        Args:
            output_dir: Directory to save figures.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        for model_name, proba_df in self.probabilities.items():
            plt.figure(figsize=(12, 8))
            sns.heatmap(proba_df, cmap="YlGnBu", cbar_kws={"label": "Probability"})
            plt.title(f"{model_name} - Prediction Probabilities")
            plt.xlabel("Class")
            plt.ylabel("Sample")
            plt.tight_layout()

            output_path = output_dir / f"{model_name}_probabilities_heatmap.png"
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
            plt.close()
            print(f"Saved heatmap to {output_path}")


class ModelComparator:
    """Compare model performance across multiple runs with different seeds."""

    def __init__(
        self,
        X_train: pd.DataFrame,
        y_train: np.ndarray,
        X_test: pd.DataFrame,
        y_test: np.ndarray,
        feature_pipeline=None,
        n_runs: int = 5,
        random_states: Optional[list] = None,
    ):
        """Initialize comparator.

        Args:
            X_train: Training feature matrix.
            y_train: Training labels.
            X_test: Test feature matrix.
            y_test: Test labels.
            feature_pipeline: Feature selection pipeline. If None, uses raw features.
            n_runs: Number of independent runs.
            random_states: List of random states for runs. If None, uses 0 to n_runs-1.
        """
        self.X_train = X_train
        self.y_train = y_train
        self.X_test = X_test
        self.y_test = y_test
        self.feature_pipeline = feature_pipeline
        self.n_runs = n_runs
        self.random_states = random_states or list(range(n_runs))

        self.results = []
        self.summary = None

    def run_comparisons(self) -> pd.DataFrame:
        """Run model training and evaluation multiple times.

        Returns:
            DataFrame with results from all runs.
        """
        print(f"Running model comparisons with {self.n_runs} random seeds...")

        results_list = []

        for run_idx, random_state in enumerate(self.random_states):
            print(f"  Run {run_idx + 1}/{self.n_runs} (seed={random_state})...")

            # Prepare data
            X_train = self.X_train.copy()
            X_test = self.X_test.copy()

            if self.feature_pipeline is not None:
                X_train = self.feature_pipeline.fit_transform(X_train, self.y_train)
                X_test = self.feature_pipeline.transform(X_test)

            # Train models
            models = ModelFactory.create_models(random_state=random_state)
            trainer = ModelTrainer(
                X_train,
                self.y_train,
                models=models,
                apply_smote=True,
                random_state=random_state,
            )
            trainer.train_all_models()

            # Evaluate on test set
            for model_name in trainer.trained_models.keys():
                y_pred = trainer.predict(X_test, model_name)
                y_pred_decoded = trainer.decode_predictions(y_pred)

                scores = make_score(self.y_test, y_pred_decoded)

                for score_name, score_value in scores.items():
                    results_list.append({
                        "run": run_idx,
                        "seed": random_state,
                        "model": model_name,
                        "metric": score_name,
                        "value": score_value,
                    })

        self.results = pd.DataFrame(results_list)
        return self.results

    def summarize_results(self) -> pd.DataFrame:
        """Create summary table with mean ± std for each model and metric.

        Returns:
            DataFrame with mean and std for each model/metric combination.
        """
        if self.results.empty:
            raise ValueError("Run comparisons first")

        summary_data = []

        for model_name in self.results["model"].unique():
            for metric in self.results["metric"].unique():
                subset = self.results[
                    (self.results["model"] == model_name)
                    & (self.results["metric"] == metric)
                ]

                if len(subset) > 0:
                    mean_val = subset["value"].mean()
                    std_val = subset["value"].std()

                    summary_data.append({
                        "Model": model_name,
                        "Metric": metric,
                        "Mean": mean_val,
                        "Std": std_val,
                        "Score": f"{mean_val:.2f} ± {std_val:.2f}",
                    })

        self.summary = pd.DataFrame(summary_data)
        return self.summary

    def save_results(self, output_dir: Path):
        """Save comparison results to CSV files.

        Args:
            output_dir: Directory for output files.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Save raw results
        results_path = output_dir / "model_comparison_raw.csv"
        self.results.to_csv(results_path, index=False)
        print(f"Saved raw results to {results_path}")

        # Save summary
        if self.summary is not None:
            summary_path = output_dir / "model_comparison_summary.csv"
            self.summary.to_csv(summary_path, index=False)
            print(f"Saved summary to {summary_path}")

    def create_comparison_table(self) -> pd.DataFrame:
        """Create pivot table matching the example.docx format.

        Returns:
            DataFrame in the format: Preprocessing | Model | Metric | Score
        """
        if self.summary is None:
            self.summarize_results()

        # Create the table in the specified format
        table_data = []
        for _, row in self.summary.iterrows():
            table_data.append({
                "Preprocessing": "",  # Will be filled by caller
                "Model": row["Model"],
                "Score Name": row["Metric"],
                "Score": row["Score"],
            })

        return pd.DataFrame(table_data)

    def create_comparison_figure(
        self, output_dir: Path, output_filename: str = "model_comparison.png"
    ) -> Path:
        """Create visualization of model performance comparisons.

        Args:
            output_dir: Directory to save figure.
            output_filename: Name for output figure.

        Returns:
            Path to saved figure.
        """
        if self.results.empty:
            raise ValueError("Run comparisons first")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Create figure with subplots for each metric
        metrics = sorted(self.results["metric"].unique())
        n_metrics = len(metrics)

        fig, axes = plt.subplots(
            n_metrics, 1, figsize=(12, 4 * n_metrics), sharex=True
        )
        if n_metrics == 1:
            axes = [axes]

        models = sorted(self.results["model"].unique())
        x_pos = np.arange(len(models))

        for ax_idx, metric in enumerate(metrics):
            ax = axes[ax_idx]

            # Get data for this metric
            metric_data = self.results[self.results["metric"] == metric]

            # Calculate mean and std for each model
            means = []
            stds = []
            for model in models:
                model_data = metric_data[metric_data["model"] == model]
                means.append(model_data["value"].mean())
                stds.append(model_data["value"].std())

            # Create bar plot
            bars = ax.bar(x_pos, means, yerr=stds, capsize=5, alpha=0.7, color="steelblue")

            ax.set_ylabel(metric.capitalize(), fontsize=11, fontweight="bold")
            ax.set_ylim(0, 1.05)
            ax.grid(True, alpha=0.3, axis="y")

            if ax_idx == len(metrics) - 1:
                ax.set_xticks(x_pos)
                ax.set_xticklabels(models, rotation=45, ha="right")
                ax.set_xlabel("Model", fontsize=11, fontweight="bold")
            else:
                ax.set_xticks([])

        plt.suptitle("Model Performance Comparison (5 runs, mean ± std)", fontsize=13, fontweight="bold")
        plt.tight_layout()

        output_path = output_dir / output_filename
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Saved comparison figure to {output_path}")

        return output_path
