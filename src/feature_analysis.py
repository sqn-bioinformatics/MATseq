"""Feature selection analysis and Venn diagram generation."""

import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, Union
from collections import Counter
from matplotlib_venn import venn2

from .feature_engineering import create_feature_pipeline
from .config import FEATURE_SELECTION_CONFIG


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

        self.feature_sets = []
        gene_counts = Counter()

        print(f"Running feature selection {n_runs} times...")
        for i, random_state in enumerate(random_states):
            if (i + 1) % 100 == 0:
                print(f"  Completed {i + 1}/{n_runs} runs...")

            pipe = create_feature_pipeline(
                **FEATURE_SELECTION_CONFIG, random_state=random_state
            ).set_output(transform="pandas")
            pipe.fit_transform(X, y)
            selected_genes = pipe[:-1].get_feature_names_out()
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
