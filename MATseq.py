#!/usr/bin/env python
"""MAT-seq pipeline orchestration script.

Main entry point for running the complete MAT-seq analysis pipeline:
1. Data preprocessing and normalization
2. Feature selection and model training
3. DESeq2 differential expression analysis
4. Downstream GO enrichment and gene set analysis
5. Venn diagram generation and feature stability analysis

This script uses disk-based caching to skip completed steps on re-runs.
"""

import sys
from pathlib import Path
import argparse
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

# Add parent directory to path for src imports
sys.path.insert(0, str(Path(__file__).parent))

from src import (
    # Configuration
    CUSTOM_PALETTE_6,
    SUBSET_PALETTES,
    SUBSET_CLASS_ORDERS,
    TRAINING_LIGANDS,
    ADDITIONAL_TLR_LIGANDS,
    DESEQ2_CONFIG,
    FEATURE_SELECTION_CONFIG,
    MODEL_TRAINING_CONFIG,
    # Preprocessing
    merge_counts,
    extract_subset,
    normalize_rpm,
    # Feature engineering
    create_feature_pipeline,
    # Model training
    ModelFactory,
    ModelTrainer,
    # Evaluation
    make_score,
    get_confusion_matrix,
    make_probability_matrix,
    # Visualization
    plot_pca_for_pandas,
    plot_gene_expression_by_class,
    plot_tlr_hek_blue,
    # DESeq2 analysis
    AnalysisPipeline,
    # Feature analysis
    FeatureSelectionAnalyzer,
    VennDiagramGenerator,
    DownstreamGOAnalysis,
    # Prediction and comparison
    ModelPredictor,
    ModelComparator,
)

from src.cache import PipelineCache


class MATseqPipeline:
    """Complete MAT-seq analysis pipeline with caching."""

    def __init__(self, cache_dir: Path = None, force_recompute: bool = False):
        """Initialize pipeline with caching.

        Args:
            cache_dir: Directory for cached intermediate files.
            force_recompute: Skip cache and recompute all steps.
        """
        self.cache = PipelineCache(cache_dir)
        self.force_recompute = force_recompute
        self.results_dir = Path.cwd() / "results"
        self.results_dir.mkdir(exist_ok=True, parents=True)

    def run_preprocessing(self) -> pd.DataFrame:
        """Load and preprocess count data.

        Returns:
            Merged count matrix.
        """
        def _merge():
            print("\n--- STEP 1: DATA PREPROCESSING ---")
            print("Merging featurecounts files...")
            df = merge_counts()
            print(f"Merged data shape: {df.shape}")
            return df

        return self.cache.cached_call(
            _merge,
            name="merged_counts",
            force_recompute=self.force_recompute,
        )

    def run_deseq2_analysis(self, raw_counts: pd.DataFrame) -> tuple:
        """Run DESeq2 differential expression analysis.

        Args:
            raw_counts: Raw count matrix (samples as rows, genes as columns).

        Returns:
            Tuple of (pipeline, de_genes_set).
        """
        print("\n--- STEP 2: DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS ---")

        pipeline = AnalysisPipeline(
            raw_counts=raw_counts,
            padj_threshold=DESEQ2_CONFIG.get("padj_threshold", 0.05),
            log2fc_threshold=DESEQ2_CONFIG.get("log2fc_threshold", 2),
            n_cpus=DESEQ2_CONFIG.get("n_cpus", 42),
        )

        # Run analysis for training ligands
        print(f"\nAnalyzing training ligands: {TRAINING_LIGANDS}")
        pipeline.run_analysis(TRAINING_LIGANDS, negative_control="IMDM")

        # Save DE genes
        pipeline.save_de_genes()

        return pipeline, pipeline.de_genes

    def run_feature_selection_analysis(
        self, X: pd.DataFrame, y: np.ndarray, n_runs: int = 1000
    ) -> tuple:
        """Run feature selection multiple times for stability analysis.

        Args:
            X: Feature matrix.
            y: Target labels.
            n_runs: Number of selection runs.

        Returns:
            Tuple of (analyzer, fs_genes_set).
        """
        print("\n--- STEP 3: FEATURE SELECTION STABILITY ANALYSIS ---")
        print(f"Running feature selection {n_runs} times...")

        analyzer = FeatureSelectionAnalyzer()

        feature_sets = analyzer.run_multiple_selections(
            X=X,
            y=y,
            n_runs=n_runs,
        )

        analyzer.save_feature_sets()
        analyzer.create_gene_frequency_table()

        # Get union of all selected features
        all_fs_genes = set()
        for gene_set in analyzer.feature_sets:
            all_fs_genes.update(gene_set)

        return analyzer, all_fs_genes

    def run_venn_analysis(
        self, de_genes: set, fs_genes: set, de_genes_path: Path
    ) -> VennDiagramGenerator:
        """Generate Venn diagrams for gene set comparisons.

        Args:
            de_genes: Set of differentially expressed genes.
            fs_genes: Set of feature-selected genes.
            de_genes_path: Path to DE genes file.

        Returns:
            VennDiagramGenerator instance.
        """
        print("\n--- STEP 4: VENN DIAGRAM GENERATION ---")

        venn_gen = VennDiagramGenerator()

        # DE vs Feature Selection
        venn_gen.plot_de_vs_feature_selection(
            de_genes,
            fs_genes,
            title="DE genes vs Feature Selected genes",
            output_filename="venn_de_vs_fs.png",
        )

        print(f"DE genes: {len(de_genes)}")
        print(f"Feature-selected genes: {len(fs_genes)}")
        print(f"Overlap: {len(de_genes & fs_genes)}")

        return venn_gen

    def process_subset(
        self, df: pd.DataFrame, subset_name: str
    ) -> tuple[pd.DataFrame, pd.Series]:
        """Process a data subset.

        Args:
            df: Full count matrix.
            subset_name: Name of subset.

        Returns:
            Tuple of (features_df, labels_series).
        """
        print(f"\nProcessing {subset_name} subset...")

        subset_data, labels = extract_subset(df, name=subset_name)
        X = subset_data.drop(columns=["label"], errors="ignore")

        # Normalize
        X_rpm = normalize_rpm(X.copy())

        # Scale
        X_scaled = StandardScaler().fit_transform(X_rpm)
        X_scaled = pd.DataFrame(X_scaled, index=X.index, columns=X.columns)

        # Plot before feature selection
        palette = SUBSET_PALETTES.get(subset_name, CUSTOM_PALETTE_6)
        hue_order = SUBSET_CLASS_ORDERS.get(subset_name)

        print(f"Creating PCA plots for {subset_name}...")
        plot_pca_for_pandas(
            name=f"{subset_name}_normalized",
            df=X_scaled,
            labels=labels,
            with_sample_names=False,
            output_filename=f"{subset_name}_pca_before_feature_selection.png",
            palette=palette,
            hue_order=hue_order,
        )

        plot_pca_for_pandas(
            name=f"{subset_name}_normalized_labeled",
            df=X_scaled,
            labels=labels,
            with_sample_names=True,
            output_filename=f"{subset_name}_pca_normalized_labeled.png",
            palette=palette,
            hue_order=hue_order,
        )

        return X_scaled, labels

    def run_feature_selection(
        self, X: pd.DataFrame, y: pd.Series, subset_name: str
    ) -> tuple[np.ndarray, object]:
        """Run feature selection pipeline.

        Args:
            X: Feature matrix.
            y: Labels.
            subset_name: Name of subset.

        Returns:
            Tuple of (selected_features, pipeline).
        """
        print(f"Running feature selection for {subset_name}...")

        def _select():
            pipe = create_feature_pipeline(**FEATURE_SELECTION_CONFIG)
            X_selected = pipe.fit_transform(X, y)
            return X_selected, pipe

        cache_key = f"feature_selection_{subset_name}"
        cached = self.cache.get(cache_key, params={"subset": subset_name})

        if cached is not None and not self.force_recompute:
            X_selected, pipe = cached
        else:
            X_selected, pipe = _select()
            self.cache.set(
                cache_key,
                (X_selected, pipe),
                params={"subset": subset_name},
                description=f"Feature-selected data for {subset_name}",
            )

        palette = SUBSET_PALETTES.get(subset_name, CUSTOM_PALETTE_6)
        hue_order = SUBSET_CLASS_ORDERS.get(subset_name)

        plot_pca_for_pandas(
            name=f"{subset_name}_selected",
            df=X_selected,
            labels=y,
            with_sample_names=False,
            output_filename=f"{subset_name}_pca_selected.png",
            palette=palette,
            hue_order=hue_order,
        )

        plot_pca_for_pandas(
            name=f"{subset_name}_selected_labeled",
            df=X_selected,
            labels=y,
            with_sample_names=True,
            output_filename=f"{subset_name}_pca_selected_labeled.png",
            palette=palette,
            hue_order=hue_order,
        )

        return X_selected, pipe

    def train_models(
        self, X: np.ndarray, y: pd.Series, subset_name: str
    ) -> ModelTrainer:
        """Train classification models.

        Args:
            X: Feature matrix (feature-selected).
            y: Labels.
            subset_name: Name of subset.

        Returns:
            Trained ModelTrainer instance.
        """
        print(f"Training models for {subset_name}...")

        def _train():
            models = ModelFactory.create_models(**MODEL_TRAINING_CONFIG)
            trainer = ModelTrainer(X, y, models=models, **MODEL_TRAINING_CONFIG)
            trainer.train_all_models()
            return trainer

        cache_key = f"trained_models_{subset_name}"
        cached = self.cache.get(cache_key, params={"subset": subset_name})

        if cached is not None and not self.force_recompute:
            trainer = cached
        else:
            trainer = _train()
            self.cache.set(
                cache_key,
                trainer,
                params={"subset": subset_name},
                description=f"Trained models for {subset_name}",
            )

        return trainer

    def run_pipeline(self):
        """Execute complete pipeline."""
        print("=" * 80)
        print("MAT-seq Analysis Pipeline")
        print("=" * 80)

        # Step 1: Preprocessing
        df = self.run_preprocessing()

        # Step 2: DESeq2 Analysis
        deseq2_pipeline, de_genes = self.run_deseq2_analysis(df)

        # Step 3: Feature Selection Analysis (on training subset)
        print("\n--- PREPARING DATA FOR FEATURE SELECTION ANALYSIS ---")
        X_train, y_train = self.process_subset(df, "training")
        fs_analyzer, fs_genes = self.run_feature_selection_analysis(
            X_train, y_train, n_runs=1000
        )

        # Step 4: Venn Diagram Generation
        de_genes_path = (
            Path.cwd() / "results" / "differential_gene_expression" / "de_genes.txt"
        )
        venn_gen = self.run_venn_analysis(de_genes, fs_genes, de_genes_path)

        # Step 5: Model Training (on training subset)
        print("\n--- STEP 5: MODEL TRAINING ---")
        X_fs, pipe = self.run_feature_selection("training", X_train, y_train)
        trainer = self.train_models(X_fs, y_train, "training")

        print("\n" + "=" * 80)
        print("PIPELINE COMPLETED SUCCESSFULLY")
        print("=" * 80)
        print(f"Results saved to: {self.results_dir.absolute()}")
        print(f"Cache saved to: {self.cache.cache_dir.absolute()}")

        return {
            "deseq2": deseq2_pipeline,
            "feature_selection": fs_analyzer,
            "venn": venn_gen,
            "models": trainer,
        }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="MAT-seq analysis pipeline with caching and full analysis"
    )
    parser.add_argument(
        "--force-recompute",
        action="store_true",
        help="Skip cache and recompute all steps",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Directory for cached files (default: results/cache)",
    )

    args = parser.parse_args()

    pipeline = MATseqPipeline(
        cache_dir=args.cache_dir, force_recompute=args.force_recompute
    )

    results = pipeline.run_pipeline()

    return results


if __name__ == "__main__":
    main()
