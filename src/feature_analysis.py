"""Feature selection analysis and Venn diagram generation."""

import pickle
import pandas as pd
import numpy as np
from functools import partial
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Union
from collections import Counter
from matplotlib_venn import venn2, venn3
from joblib import Parallel, delayed
from sklearn.model_selection import GridSearchCV
from sklearn.feature_selection import SelectKBest, mutual_info_classif
from sklearn.ensemble import ExtraTreesClassifier

from .feature_engineering import (
    create_feature_pipeline,
    create_preprocessing_pipeline,
    create_selection_pipeline,
)
from .config import FEATURE_SELECTION_CONFIG


def _run_one_selection(X_pre, y, random_state):
    sel = create_selection_pipeline(
        **FEATURE_SELECTION_CONFIG, random_state=int(random_state)
    ).set_output(transform="pandas")
    sel.fit(X_pre, y)
    return sel.get_feature_names_out()


def _elbow_index(scores) -> int:
    y = np.asarray(scores, dtype=float)
    x = np.arange(len(y), dtype=float)
    xn = (x - x.min()) / (np.ptp(x) or 1)
    yn = (y - y.min()) / (np.ptp(y) or 1)
    num = np.abs(
        (yn[-1] - yn[0]) * xn
        - (xn[-1] - xn[0]) * yn
        + xn[-1] * yn[0]
        - yn[-1] * xn[0]
    )
    return int(np.argmax(num)) + 1


def feature_count_analysis(
    X: pd.DataFrame,
    y: Union[np.ndarray, pd.Series],
    k_best: int,
    candidate_counts: List[int],
    n_runs: int = 200,
    n_estimators: int = 250,
    max_depth: int = 5,
    random_state: int = 42,
) -> Dict:
    """Model-independent analysis of how many features to select.

    Combines an intrinsic-signal elbow (sorted mutual-information /
    ExtraTrees-importance curves) with a selection-stability curve (mean pairwise
    Jaccard of reseeded top-N selections). Both are derived from one reseeded
    ExtraTrees pass on the SelectKBest(k_best) matrix; no classifier performance
    is used. Returns tidy DataFrames and the suggested counts; writes nothing.
    """
    pre = create_preprocessing_pipeline(
        drop_low_count=True, scale=False
    ).set_output(transform="pandas")
    X_pre = pre.fit_transform(X, y)

    mi = np.mean(
        [
            mutual_info_classif(X_pre, y, random_state=random_state + i)
            for i in range(3)
        ],
        axis=0,
    )
    mi_series = pd.Series(mi, index=X_pre.columns).sort_values(ascending=False)
    mi_elbow = _elbow_index(mi_series.to_numpy())

    k = min(k_best, X_pre.shape[1])
    X_k = (
        SelectKBest(partial(mutual_info_classif, random_state=random_state), k=k)
        .set_output(transform="pandas")
        .fit_transform(X_pre, y)
    )
    k_features = list(X_k.columns)

    rankings = []
    importance_sum = np.zeros(X_k.shape[1])
    for i in range(n_runs):
        et = ExtraTreesClassifier(
            n_estimators=n_estimators,
            max_depth=max_depth,
            random_state=random_state + i,
        )
        et.fit(X_k, y)
        imp = et.feature_importances_
        importance_sum += imp
        order = np.argsort(imp)[::-1]
        rankings.append([k_features[j] for j in order])

    mean_imp = pd.Series(importance_sum / n_runs, index=k_features).sort_values(
        ascending=False
    )
    importance_elbow = _elbow_index(mean_imp.to_numpy())

    counts = [c for c in candidate_counts if c <= len(k_features)]
    jaccard, core_50, core_75, core_90 = [], [], [], []
    for n in counts:
        topn = [set(r[:n]) for r in rankings]
        pair_jac = [
            len(a & b) / len(a | b)
            for i, a in enumerate(topn)
            for b in topn[i + 1 :]
            if a | b
        ]
        jaccard.append(float(np.mean(pair_jac)) if pair_jac else 0.0)
        freq = Counter()
        for s in topn:
            freq.update(s)
        core_50.append(sum(1 for v in freq.values() if v / n_runs >= 0.50))
        core_75.append(sum(1 for v in freq.values() if v / n_runs >= 0.75))
        core_90.append(sum(1 for v in freq.values() if v / n_runs >= 0.90))

    stable_count = counts[_elbow_index(jaccard) - 1] if counts else None

    ref = min(importance_elbow, len(k_features))
    ref_freq = Counter()
    for r in rankings:
        ref_freq.update(r[:ref])
    core_df = (
        pd.Series({g: c / n_runs for g, c in ref_freq.items()})
        .sort_values(ascending=False)
        .rename_axis("gene")
        .reset_index(name="selection_frequency")
    )

    scores_df = pd.DataFrame(
        {"rank": np.arange(1, len(mi_series) + 1), "mi_sorted": mi_series.to_numpy()}
    ).merge(
        pd.DataFrame(
            {
                "rank": np.arange(1, len(mean_imp) + 1),
                "importance_sorted": mean_imp.to_numpy(),
            }
        ),
        on="rank",
        how="left",
    )

    stability_df = pd.DataFrame(
        {
            "n_features": counts,
            "mean_jaccard": jaccard,
            "core_50pct": core_50,
            "core_75pct": core_75,
            "core_90pct": core_90,
        }
    )

    return {
        "mi_elbow": mi_elbow,
        "importance_elbow": importance_elbow,
        "stable_count": stable_count,
        "scores": scores_df,
        "stability": stability_df,
        "core": core_df,
    }


class FeatureSelectionAnalyzer:
    """Analyze feature selection variability across multiple runs."""

    def __init__(self, output_dir: Path = None):
        """Initialize analyzer.

        Args:
            output_dir: Directory for output files. Defaults to results/feature_analysis.
        """
        if output_dir is None:
            output_dir = Path.cwd() / "results" / "feature_analysis"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.gene_frequency = None
        self.feature_sets = []

    def run_multiple_selections(
        self,
        X: pd.DataFrame,
        y: Union[np.ndarray, pd.Series],
        n_runs: int = 1000,
    ) -> Dict[int, set]:
        """Run feature selection pipeline multiple times with different seeds.

        Args:
            X: Feature matrix (samples x features).
            y: Target labels.
            n_runs: Number of runs to perform.

        Returns:
            Dictionary mapping run number to selected gene set.
        """
        if len(X) == 0:
            raise ValueError("X cannot be empty")
        if len(y) == 0:
            raise ValueError("y cannot be empty")
        if n_runs <= 0:
            raise ValueError(f"n_runs must be positive, got {n_runs}")

        rng = np.random.default_rng(125)
        random_states = rng.integers(0, 500000, size=n_runs)

        pre = create_preprocessing_pipeline(
            drop_low_count=True, scale=False
        ).set_output(transform="pandas")
        X_pre = pre.fit_transform(X, y)

        print(f"Running feature selection {n_runs} times...")
        results = Parallel(n_jobs=-1, verbose=10, return_as="generator")(
            delayed(_run_one_selection)(X_pre, y, rs) for rs in random_states
        )

        self.feature_sets = []
        gene_counts = Counter()
        for selected_genes in results:
            self.feature_sets.append(selected_genes)
            gene_counts.update(selected_genes)

        self.gene_frequency = gene_counts
        print(f"Completed {n_runs} feature selection runs")
        print(f"Unique genes across all runs: {len(gene_counts)}")

        return dict(enumerate(self.feature_sets))

    def create_gene_frequency_table(
        self, de_genes: set = None, output_filename: str = "gene_frequency_table.csv"
    ) -> pd.DataFrame:
        """Create table of gene frequencies with DE status.

        Args:
            de_genes: Set of differentially expressed genes.
            output_filename: Name for output CSV file.

        Returns:
            DataFrame with columns: Gene, Count, Differentially Expressed.
        """
        if self.gene_frequency is None:
            raise ValueError("Run feature selection first")

        data = []
        for gene, count in self.gene_frequency.most_common():
            is_de = gene in de_genes if de_genes else False
            data.append(
                {
                    "Gene": gene,
                    "Count": count,
                    "Differentially Expressed": is_de,
                }
            )

        df = pd.DataFrame(data)

        output_path = self.output_dir / output_filename
        df.to_csv(output_path, index=False)
        print(f"Saved gene frequency table to {output_path}")

        return df

    def save_feature_sets(self, output_filename: str = "list_of_gene_name_set_1000"):
        """Save feature sets to pickle file.

        Args:
            output_filename: Name for output pickle file (without extension).
        """
        if not self.feature_sets:
            raise ValueError("Run feature selection first")

        output_path = self.output_dir / output_filename
        with open(output_path, "wb") as f:
            pickle.dump(self.feature_sets, f)

        print(f"Saved {len(self.feature_sets)} feature sets to {output_path}")

    def load_feature_sets(self, input_path: Path):
        """Load previously saved feature sets.

        Args:
            input_path: Path to pickle file with feature sets.
        """
        with open(input_path, "rb") as f:
            self.feature_sets = pickle.load(f)

        gene_counts = Counter()
        for gene_set in self.feature_sets:
            gene_counts.update(gene_set)

        self.gene_frequency = gene_counts
        print(f"Loaded {len(self.feature_sets)} feature sets")



class PipelineParamTuner:
    """Grid search over k_best and max_features pipeline parameters."""

    def __init__(self, output_dir: Path = None):
        if output_dir is None:
            output_dir = Path.cwd() / "results" / "feature_analysis"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.cv_results_ = None

    def run(
        self,
        X: pd.DataFrame,
        y: Union[np.ndarray, pd.Series],
        k_best_values: List[int],
        max_features_values: List[int],
        scoring: str = "balanced_accuracy",
        cv: int = 5,
    ) -> pd.DataFrame:
        param_grid = {
            "select_k_best__k": k_best_values,
            "select_forest__max_features": max_features_values,
        }
        pipe = create_feature_pipeline(**FEATURE_SELECTION_CONFIG)
        gs = GridSearchCV(pipe, param_grid, scoring=scoring, cv=cv, n_jobs=-1, refit=False)
        gs.fit(X, y)

        results = pd.DataFrame(gs.cv_results_)
        self.cv_results_ = results
        return results

    def plot(self, output_filename: str = "param_tuning_heatmap.png") -> Path:
        if self.cv_results_ is None:
            raise ValueError("Run grid search first")

        pivot = self.cv_results_.pivot(
            index="param_select_k_best__k",
            columns="param_select_forest__max_features",
            values="mean_test_score",
        )

        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(pivot.values, aspect="auto")
        plt.colorbar(im, ax=ax, label="mean test score")

        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels(pivot.columns)
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels(pivot.index)
        ax.set_xlabel("max_features")
        ax.set_ylabel("k_best")

        for i in range(len(pivot.index)):
            for j in range(len(pivot.columns)):
                ax.text(j, i, f"{pivot.values[i, j]:.3f}", ha="center", va="center", fontsize=8)

        output_path = self.output_dir / output_filename
        try:
            fig.savefig(output_path, dpi=150, bbox_inches="tight")
        finally:
            plt.close(fig)

        return output_path


class VennDiagramGenerator:
    """Generate Venn diagrams for gene set comparisons."""

    def __init__(self, output_dir: Path = None):
        """Initialize generator.

        Args:
            output_dir: Directory for output figures. Defaults to results/figures/venn.
        """
        if output_dir is None:
            output_dir = Path.cwd() / "results" / "figures" / "venn"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def plot_venn(
        self,
        de_genes: set,
        fs_genes: set,
        output_filename: str,
    ) -> Path:
        """Create 2-way Venn diagram comparing DE and feature-selected genes.

        Args:
            de_genes: Set of differentially expressed genes.
            fs_genes: Set of feature-selected genes.
            output_filename: Name for output file.

        Returns:
            Path to saved figure.
        """
        plt.figure(figsize=(8, 8))

        venn2(
            [de_genes, fs_genes],
            set_labels=("Differentially Expressed Genes", "Feature Selection Genes"),
        )
        output_path = self.output_dir / output_filename
        try:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
            print(f"Saved Venn diagram to {output_path}")
        finally:
            plt.close()

        return output_path

    def plot_venn3(
        self,
        sets: List[set],
        set_labels: tuple,
        output_filename: str,
        title: str = None,
    ) -> Path:
        """Create a 3-way Venn diagram comparing three gene sets.

        Args:
            sets: List of three gene sets.
            set_labels: Labels for the three sets.
            output_filename: Name for output file.
            title: Optional figure title.

        Returns:
            Path to saved figure.
        """
        plt.figure(figsize=(8, 8))

        venn3(sets, set_labels=set_labels)
        if title:
            plt.title(title)
        output_path = self.output_dir / output_filename
        try:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
            print(f"Saved Venn diagram to {output_path}")
        finally:
            plt.close()

        return output_path
