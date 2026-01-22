"""Model prediction and evaluation on new data and multiple runs."""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Optional
import matplotlib.pyplot as plt
import seaborn as sns

from .model_training import ModelTrainer, ModelFactory, make_score
from .config import SUBSET_CLASS_ORDERS


class ModelPredictor:
    """Make predictions on new samples with trained models."""

    def __init__(self, trainer: ModelTrainer):
        """Initialize predictor with trained models.

        Args:
            trainer: Trained ModelTrainer instance with fitted models.
        """
        assert hasattr(trainer, 'trained_models'), "trainer must have trained_models attribute"
        assert len(trainer.trained_models) > 0, "trainer must have at least one trained model"
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

            # Create prediction DataFrame
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

            # Get probabilities if available
            if hasattr(model, "predict_proba"):
                proba = model.predict_proba(X_test)
                proba_df = pd.DataFrame(
                    proba,
                    columns=self.trainer.label_encoder.classes_,
                    index=(
                        sample_names
                        if sample_names is not None
                        else np.arange(len(y_pred))
                    ),
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

    def create_probability_heatmaps(self, output_dir: Path, class_order: list = None):
        """Create heatmap visualizations of prediction probabilities.

        Args:
            output_dir: Directory to save figures.
            class_order: Optional class order from config.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if class_order is None:
            class_order = SUBSET_CLASS_ORDERS.get("training")

        for model_name, proba_df in self.probabilities.items():
            proba_df_ordered = proba_df.copy()
            if class_order is not None:
                available_classes = [c for c in class_order if c in proba_df.columns]
                proba_df_ordered = proba_df_ordered[available_classes]

            if self.true_labels is not None:
                proba_df_ordered.index = self.true_labels.values

            plt.figure(figsize=(12, 8))
            sns.heatmap(
                proba_df_ordered, cmap="YlGnBu", cbar_kws={"label": "Probability"}
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

            if self.true_labels is not None:
                mask = self.true_labels.isin(["LPS", "negative_control"])
                if mask.any():
                    lps_idx = self.true_labels[self.true_labels == "LPS"].index
                    nc_idx = self.true_labels[self.true_labels == "negative_control"].index
                    if len(lps_idx) > 0 and len(nc_idx) > 0:
                        sampled_lps = np.random.choice(lps_idx, size=1)
                        subset_idx = np.concatenate([sampled_lps, nc_idx.values])

                        proba_subset = proba_df.loc[subset_idx].copy()
                        if class_order is not None:
                            available_classes = [c for c in class_order if c in proba_subset.columns]
                            proba_subset = proba_subset[available_classes]

                        proba_subset.index = self.true_labels.loc[subset_idx].values

                        plt.figure(figsize=(12, 8))
                        sns.heatmap(
                            proba_subset, cmap="YlGnBu", cbar_kws={"label": "Probability"}
                        )
                        plt.title(f"{model_name}")
                        plt.xlabel("Class")
                        plt.ylabel("Sample")
                        plt.tight_layout()

                        output_path_subset = output_dir / f"{model_name}_probabilities_heatmap_subset.png"
                        try:
                            plt.savefig(output_path_subset, dpi=300, bbox_inches="tight")
                            print(f"Saved subset heatmap to {output_path_subset}")
                        finally:
                            plt.close()


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
        assert len(X_train) > 0, "X_train cannot be empty"
        assert len(y_train) > 0, "y_train cannot be empty"
        assert len(X_test) > 0, "X_test cannot be empty"
        assert len(y_test) > 0, "y_test cannot be empty"
        assert len(X_train) == len(y_train), "X_train and y_train must have same length"
        assert len(X_test) == len(y_test), "X_test and y_test must have same length"
        assert n_runs > 0, f"n_runs must be positive, got {n_runs}"

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
            X_test_array = X_test.values if isinstance(X_test, pd.DataFrame) else X_test
            for model_name in trainer.trained_models.keys():
                y_pred = trainer.predict(X_test_array, model_name)
                y_pred_decoded = trainer.decode_predictions(y_pred)

                scores = make_score(self.y_test, y_pred_decoded)

                for score_name, score_value in scores.items():
                    results_list.append(
                        {
                            "run": run_idx,
                            "seed": random_state,
                            "model": model_name,
                            "metric": score_name,
                            "value": score_value,
                        }
                    )

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

                    summary_data.append(
                        {
                            "Model": model_name,
                            "Metric": metric,
                            "Mean": mean_val,
                            "Std": std_val,
                            "Score": f"{mean_val:.2f} ± {std_val:.2f}",
                        }
                    )

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
            table_data.append(
                {
                    "Preprocessing": "",  # Will be filled by caller
                    "Model": row["Model"],
                    "Score Name": row["Metric"],
                    "Score": row["Score"],
                }
            )

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

        fig, axes = plt.subplots(n_metrics, 1, figsize=(12, 4 * n_metrics), sharex=True)
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
            bars = ax.bar(
                x_pos, means, yerr=stds, capsize=5, alpha=0.7, color="steelblue"
            )

            ax.set_ylabel(metric.capitalize(), fontsize=11, fontweight="bold")
            ax.set_ylim(0, 1.05)
            ax.grid(True, alpha=0.3, axis="y")

            if ax_idx == len(metrics) - 1:
                ax.set_xticks(x_pos)
                ax.set_xticklabels(models, rotation=45, ha="right")
                ax.set_xlabel("Model", fontsize=11, fontweight="bold")
            else:
                ax.set_xticks([])

        plt.suptitle(
            "Model Performance Comparison (5 runs, mean ± std)",
            fontsize=13,
            fontweight="bold",
        )
        plt.tight_layout()

        output_path = output_dir / output_filename
        try:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
            print(f"Saved comparison figure to {output_path}")
        finally:
            plt.close()

        return output_path
