"""Feature selection analysis and Venn diagram generation."""

import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Tuple
from collections import Counter
from matplotlib_venn import venn2, venn3

from .feature_engineering import create_feature_pipeline


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

        self.figures_dir = Path.cwd() / "results" / "figures" / "feature_analysis"
        self.figures_dir.mkdir(parents=True, exist_ok=True)

        self.gene_frequency = None
        self.feature_sets = []

    def run_multiple_selections(
        self,
        X: pd.DataFrame,
        y: np.ndarray,
        n_runs: int = 1000,
        random_states: List[int] = None,
    ) -> Dict[int, set]:
        """Run feature selection pipeline multiple times with different seeds.

        Args:
            X: Feature matrix (samples x features).
            y: Target labels.
            n_runs: Number of runs to perform.
            random_states: List of random states. If None, uses range(n_runs).

        Returns:
            Dictionary mapping run number to selected gene set.
        """
        if random_states is None:
            random_states = list(range(n_runs))

        self.feature_sets = []
        gene_counts = Counter()

        print(f"Running feature selection {n_runs} times...")
        for i, random_state in enumerate(random_states):
            if (i + 1) % 100 == 0:
                print(f"  Completed {i + 1}/{n_runs} runs...")

            # Create pipeline with specific random state
            pipe = create_feature_pipeline(random_state=random_state)

            # Get selected features
            X_selected = pipe.fit_transform(X, y)

            # Store feature names
            selected_genes = set(X_selected.columns if isinstance(X_selected, pd.DataFrame)
                                else X.columns[pipe[-1].get_support()])

            self.feature_sets.append(selected_genes)

            # Count gene frequencies
            for gene in selected_genes:
                gene_counts[gene] += 1

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
            data.append({
                "Gene": gene,
                "Count": count,
                "Differentially Expressed": is_de,
            })

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

        # Recalculate gene frequency
        gene_counts = Counter()
        for gene_set in self.feature_sets:
            for gene in gene_set:
                gene_counts[gene] += 1

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

    def plot_de_vs_feature_selection(
        self,
        de_genes: set,
        feature_selected_genes: set,
        title: str = "DE genes vs Feature Selected genes",
        output_filename: str = "venn_de_vs_fs.png",
    ) -> Path:
        """Create 2-way Venn diagram comparing DE and feature-selected genes.

        Args:
            de_genes: Set of differentially expressed genes.
            feature_selected_genes: Set of feature-selected genes.
            title: Title for the diagram.
            output_filename: Name for output file.

        Returns:
            Path to saved figure.
        """
        plt.figure(figsize=(8, 8))

        venn2(
            [de_genes, feature_selected_genes],
            set_labels=("DE genes", "Feature Selected"),
        )

        plt.title(title, fontsize=14, fontweight="bold")
        plt.tight_layout()

        output_path = self.output_dir / output_filename
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"Saved Venn diagram to {output_path}")
        plt.close()

        return output_path

    def plot_three_way_comparison(
        self,
        set1: set,
        set2: set,
        set3: set,
        labels: Tuple[str, str, str] = ("Set 1", "Set 2", "Set 3"),
        title: str = "Three-way Venn diagram",
        output_filename: str = "venn_3way.png",
    ) -> Path:
        """Create 3-way Venn diagram.

        Args:
            set1: First gene set.
            set2: Second gene set.
            set3: Third gene set.
            labels: Labels for the three sets.
            title: Title for the diagram.
            output_filename: Name for output file.

        Returns:
            Path to saved figure.
        """
        plt.figure(figsize=(10, 10))

        venn3([set1, set2, set3], set_labels=labels)

        plt.title(title, fontsize=14, fontweight="bold")
        plt.tight_layout()

        output_path = self.output_dir / output_filename
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"Saved 3-way Venn diagram to {output_path}")
        plt.close()

        return output_path


class DownstreamGOAnalysis:
    """Consolidate GO enrichment results across multiple DE analyses."""

    def __init__(self, output_dir: Path = None):
        """Initialize GO analysis.

        Args:
            output_dir: Directory for output files. Defaults to results/differential_gene_expression.
        """
        if output_dir is None:
            output_dir = Path.cwd() / "results" / "differential_gene_expression"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def merge_go_tables(
        self,
        go_files: List[Path],
        output_filename: str = "Supplementary_Table_3_genes_GO.csv",
    ) -> pd.DataFrame:
        """Merge GO enrichment results from multiple ligand analyses.

        Args:
            go_files: List of paths to GO result CSV files.
            output_filename: Name for merged output file.

        Returns:
            Merged DataFrame with all GO terms.
        """
        if not go_files:
            raise ValueError("No GO files provided")

        merged_df = None

        for go_file in go_files:
            ligand_name = go_file.stem.split("_")[0]
            df = pd.read_csv(go_file, index_col=0)

            # Add ligand identifier
            df.insert(0, "Ligand", ligand_name)

            if merged_df is None:
                merged_df = df
            else:
                merged_df = pd.concat([merged_df, df], axis=0)

        # Sort by FDR
        merged_df = merged_df.sort_values("fdr", ascending=True)

        output_path = self.output_dir / output_filename
        merged_df.to_csv(output_path)
        print(f"Saved merged GO table to {output_path}")

        return merged_df

    def create_fs_de_go_table(
        self,
        de_genes: set,
        fs_genes: set,
        go_df: pd.DataFrame = None,
        output_filename: str = "Supplementary_Table_4_fsde_GO.csv",
    ) -> pd.DataFrame:
        """Create table of GO terms for genes in both DE and feature selection.

        Args:
            de_genes: Set of differentially expressed genes.
            fs_genes: Set of feature-selected genes.
            go_df: GO enrichment results DataFrame. If None, creates stub.
            output_filename: Name for output file.

        Returns:
            DataFrame with GO terms for DE+FS genes.
        """
        # Get intersection of DE and FS genes
        fsde_genes = de_genes & fs_genes

        if go_df is None:
            # Create simple table if no GO data available
            data = {
                "Gene": list(fsde_genes),
                "In DE": [True] * len(fsde_genes),
                "In FS": [True] * len(fsde_genes),
            }
            fsde_table = pd.DataFrame(data)
        else:
            # Filter GO results to FSDE genes
            fsde_table = go_df[go_df["gene_symbols"].apply(
                lambda x: any(gene in fsde_genes for gene in x.split(", "))
            )]

        output_path = self.output_dir / output_filename
        fsde_table.to_csv(output_path, index=False)
        print(f"Saved FS+DE GO table to {output_path}")

        return fsde_table
