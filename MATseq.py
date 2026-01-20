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
    CUSTOM_PALETTE_9,
    SUBSET_PALETTES,
    SUBSET_CLASS_ORDERS,
    TRAINING_LIGANDS,
    ADDITIONAL_LIGANDS,
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

    def run_snakemake_preprocessing(
        self,
        snakemake_dir: Path = None,
        fastq_dir: Path = None,
        genome_dir: Path = None,
    ) -> bool:
        """Run Snakemake preprocessing pipeline to generate featurecounts.

        Args:
            snakemake_dir: Directory containing Snakemake pipeline (default: pipeline/).
            fastq_dir: Override FASTQ input directory.
            genome_dir: Override genome reference directory.

        Returns:
            True if successful, False otherwise.
        """
        if snakemake_dir is None:
            snakemake_dir = Path.cwd() / "pipeline"

        run_script = snakemake_dir / "run.sh"
        if not run_script.exists():
            print(f"Error: Snakemake pipeline not found at {snakemake_dir}")
            return False

        print("\n--- RUNNING SNAKEMAKE PREPROCESSING ---")
        print(f"Pipeline directory: {snakemake_dir}")

        try:
            import subprocess

            result = subprocess.run(
                ["bash", str(run_script)],
                cwd=str(snakemake_dir),
                capture_output=False,
            )

            if result.returncode == 0:
                print("✓ Snakemake preprocessing completed successfully")
                return True
            else:
                print(f"✗ Snakemake preprocessing failed with code {result.returncode}")
                return False
        except Exception as e:
            print(f"✗ Error running Snakemake: {e}")
            return False

    def run_preprocessing(self) -> pd.DataFrame:
        """Load and preprocess count data.

        Returns:
            Merged count matrix.
        """

        def _merge():
            print("Merging featurecounts files...")
            df = merge_counts()
            print(f"Merged data shape: {df.shape}")
            return df

        return self.cache.cached_call(
            _merge,
            name="merged_counts",
            force_recompute=self.force_recompute,
        )

    def run_venn_feature_selection_multiple_times(
        self, X: pd.DataFrame, y: pd.Series, n_runs: int = 1000
    ) -> tuple:
        """Run feature selection multiple times to see which genes are more prevalent.

        Args:
            X: Feature matrix.
            y: Target labels.
            n_runs: Number of selection runs.

        Returns:
            Tuple of (analyzer, fs_genes_set).
        """

        def _venn_selection():
            print("\n--- STEP 3: FEATURE SELECTION STABILITY ANALYSIS ---")
            print(f"Running feature selection {n_runs} times...")

            analyzer = FeatureSelectionAnalyzer()

            analyzer.run_multiple_selections(
                X=X,
                y=y,
                n_runs=n_runs,
            )

            analyzer.create_gene_frequency_table()

            # Get union of all selected features
            all_fs_genes = set()
            for gene_set in analyzer.feature_sets:
                all_fs_genes.update(gene_set)

            return analyzer, all_fs_genes

        cache_key = f"venn_feature_selection_{n_runs}"
        cached = self.cache.get(cache_key)

        if cached is not None and not self.force_recompute:
            analyzer, all_fs_genes = cached
        else:
            analyzer, all_fs_genes = _venn_selection()
            self.cache.set(cache_key, (analyzer, all_fs_genes))

        return analyzer, all_fs_genes

    def run_venn_analysis(self, de_genes: set, fs_genes: set) -> VennDiagramGenerator:
        """Generate Venn diagrams for gene set comparisons.

        Args:
            de_genes: Set of differentially expressed genes.
            fs_genes: Set of feature-selected genes.

        Returns:
            VennDiagramGenerator instance.
        """
        print("\n--- STEP 4: VENN DIAGRAM GENERATION ---")

        venn_gen = VennDiagramGenerator()
        venn_gen.plot_venn(
            de_genes,
            fs_genes,
            output_filename="venn_de_vs_fs.png",
        )

        return venn_gen

    def extract_and_run_pca_before_pipeline(
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

        subset_data = extract_subset(df, name=subset_name)
        y = subset_data["label"]
        X = subset_data.drop(columns=["label"])

        # Following processing steps are only to plot PCA
        X_rpm = normalize_rpm(X.copy())
        ss = StandardScaler()
        ss.fit(X_rpm)
        ss.set_output(transform="pandas")
        X_scaled = ss.transform(X_rpm)

        palette = SUBSET_PALETTES.get(subset_name, CUSTOM_PALETTE_9)
        hue_order = SUBSET_CLASS_ORDERS.get(subset_name)

        print(f"Creating PCA plots for {subset_name}...")
        common_kwargs = dict(
            df=X_scaled,
            labels=y,
            palette=palette,
            hue_order=hue_order,
        )

        for with_names, suffix in [
            (False, ""),
            (True, "_labeled"),
        ]:
            plot_pca_for_pandas(
                name=f"{subset_name}_selected{suffix}",
                with_sample_names=with_names,
                output_filename=f"{subset_name}_pca_selected{suffix}.png",
                **common_kwargs,
            )

        return X, y

    def run_feature_selection(
        self, X: pd.DataFrame, y: pd.Series, subset_name: str
    ) -> tuple[object, np.ndarray]:
        """Run feature selection pipeline.

        Args:
            X: Feature matrix.
            y: Labels.
            subset_name: Name of subset.

        Returns:
            Tuple of (selected_features, pipeline).
        """

        def _feature_select():
            print(f"Running feature selection for {subset_name}...")
            pipe = create_feature_pipeline(**FEATURE_SELECTION_CONFIG)
            pipe.fit(X, y)
            pipe.set_output(transform="pandas")
            X_selected = pipe.transform(X)
            return pipe, X_selected

        cache_key = f"feature_selection_{subset_name}"
        cached = self.cache.get(cache_key)

        if cached is not None and not self.force_recompute:
            pipe, X_selected = cached
        else:
            pipe, X_selected = _feature_select()
            self.cache.set(cache_key, (pipe, X_selected))

        palette = SUBSET_PALETTES.get(subset_name, CUSTOM_PALETTE_9)
        hue_order = SUBSET_CLASS_ORDERS.get(subset_name)

        common_kwargs = dict(
            df=X_selected,
            labels=y,
            palette=palette,
            hue_order=hue_order,
        )

        for with_names, suffix in [
            (False, ""),
            (True, "_labeled"),
        ]:
            plot_pca_for_pandas(
                name=f"{subset_name}_selected{suffix}",
                with_sample_names=with_names,
                output_filename=f"{subset_name}_pca_selected{suffix}.png",
                **common_kwargs,
            )

        return pipe, X_selected

    def train_models(
        self, X: pd.DataFrame, y: pd.Series, subset_name: str
    ) -> ModelTrainer:  #! Implement
        """Train classification models.

        Args:
            X: Feature matrix (feature-selected).
            y: Labels.
            subset_name: Name of subset.

        Returns:
            Trained ModelTrainer instance.
        """

        def _train():
            models = ModelFactory.create_models(**MODEL_TRAINING_CONFIG)
            trainer = ModelTrainer(X, y, models=models, **MODEL_TRAINING_CONFIG)
            trainer.train_all_models()
            return trainer

        cache_key = f"trained_models_{subset_name}"
        cached = self.cache.get(cache_key)

        if cached is not None and not self.force_recompute:
            trainer = cached
        else:
            trainer = _train()
            self.cache.set(
                cache_key,
                trainer,
            )

        return trainer

    def run_pipeline(self):
        print("=" * 80)
        print("MAT-seq Analysis Pipeline")
        print("=" * 80)

        # Step 1: Preprocessing
        print("\n--- STEP 1: DATA PREPROCESSING ---")
        df = self.run_preprocessing()

        # Step 2: DESeq2 Analysis (on training subset)
        print(
            "\n--- STEP 2: DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS ON TRAINING SUBSET ---"
        )
        X_train, y_train = self.extract_and_run_pca_before_pipeline(df, "training")

        deseq2_pipeline_training = AnalysisPipeline(
            raw_counts=X_train,
            sample_labels=y_train,
            padj_threshold=DESEQ2_CONFIG.get("padj_threshold", 0.05),
            log2fc_threshold=DESEQ2_CONFIG.get("log2fc_threshold", 2),
            n_cpus=DESEQ2_CONFIG.get("n_cpus", 42),
        )
        deseq2_pipeline_training.run_analysis(
            TRAINING_LIGANDS, negative_control="negative_control"
        )

        # Step 3: Venn Diagram of Feature Selection Analysis (on training subset)
        print(
            "\n--- STEP 3: GENERATING VENN DIAGRAM OF POST-FEATURE SELECTION AND DIFFERENTIALLY EXPRESSED GENES IN TRAINING SUBSET ---"
        )

        de_genes = deseq2_pipeline_training.get_de_genes()
        fs_analyzer, fs_genes = self.run_venn_feature_selection_multiple_times(
            X_train, y_train, n_runs=1000
        )
        venn_gen = VennDiagramGenerator()
        venn_gen.plot_venn(
            de_genes,
            fs_genes,
            output_filename="venn_de_vs_fs.png",
        )

        # Step 4: Model Training (on training subset)
        print("\n--- STEP 4: MODEL TRAINING ---")
        pipe, X_fs = self.run_feature_selection(
            X_train,
            y_train,
            "training",
        )
        trainer = self.train_models(X_fs, y_train, "training")

        # Step 4 DESeq2 Analysis (on other ligands and bacterial subsets)
        print(
            "\n--- STEP 5: DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS ON REMAINING LIGANDS ---"
        )
        X_other_ligands, y_other_ligands = self.extract_and_run_pca_before_pipeline(
            df, "other_ligands"
        )
        deseq2_pipeline_other_ligand = AnalysisPipeline(
            raw_counts=X_other_ligands,
            sample_labels=y_other_ligands,
            padj_threshold=DESEQ2_CONFIG.get("padj_threshold", 0.05),
            log2fc_threshold=DESEQ2_CONFIG.get("log2fc_threshold", 2),
            n_cpus=DESEQ2_CONFIG.get("n_cpus", 42),
        )
        deseq2_pipeline_other_ligand.run_analysis(
            ADDITIONAL_LIGANDS, negative_control="negative_control"
        )

        X_bacterial, y_bacterial = self.extract_and_run_pca_before_pipeline(
            df, "bacterial"
        )
        deseq2_pipeline_bacterial = AnalysisPipeline(
            raw_counts=X_bacterial,
            sample_labels=y_bacterial,
            padj_threshold=DESEQ2_CONFIG.get("padj_threshold", 0.05),
            log2fc_threshold=DESEQ2_CONFIG.get("log2fc_threshold", 2),
            n_cpus=DESEQ2_CONFIG.get("n_cpus", 42),
        )
        deseq2_pipeline_bacterial.run_analysis(
            ["HK E.coli", "HK S.aureus"],
            negative_control="negative_control",
        )

        print("\n--- STEP 7: CLASS PREDICTION ON ADDITIONAL AND BACTERIAL LIGANDS ---")

        predictor = ModelPredictor(trainer)
        predictions_dir = self.results_dir / "predictions"

        X_other_fs = pipe.transform(X_other_ligands)
        predictor.predict_samples(
            X_other_fs, sample_names=X_other_ligands.index.to_numpy()
        )
        predictor.save_predictions(predictions_dir / "other_ligands")
        predictor.create_probability_heatmaps(predictions_dir / "other_ligands")

        X_bacterial_fs = pipe.transform(X_bacterial)
        predictor.predict_samples(
            X_bacterial_fs, sample_names=X_bacterial.index.to_numpy()
        )
        predictor.save_predictions(predictions_dir / "bacterial")
        predictor.create_probability_heatmaps(predictions_dir / "bacterial")

        print("\n" + "=" * 80)
        print("PIPELINE COMPLETED SUCCESSFULLY")
        print("=" * 80)
        print(f"Results saved to: {self.results_dir.absolute()}")
        print(f"Cache saved to: {self.cache.cache_dir.absolute()}")

        return {
            "deseq2_training": deseq2_pipeline_training,
            "feature_selection": pipe,
            "feature_selection_for_venn": fs_analyzer,
            "venn": venn_gen,
            "models": trainer,
            "predictor": predictor,
            "deseq2_other_ligand": deseq2_pipeline_other_ligand,
            "deseq2_pipeline_bacterial": deseq2_pipeline_bacterial,
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
    parser.add_argument(
        "--run-snakemake",
        action="store_true",
        help="Run Snakemake preprocessing pipeline first to generate featurecounts",
    )
    parser.add_argument(
        "--snakemake-dir",
        type=Path,
        default=None,
        help="Directory containing Snakemake pipeline (default: pipeline/)",
    )
    parser.add_argument(
        "--fastq-dir",
        type=Path,
        default=None,
        help="Override FASTQ input directory for Snakemake",
    )
    parser.add_argument(
        "--genome-dir",
        type=Path,
        default=None,
        help="Override genome reference directory for Snakemake",
    )

    args = parser.parse_args()

    pipeline = MATseqPipeline(
        cache_dir=args.cache_dir, force_recompute=args.force_recompute
    )

    # Run Snakemake preprocessing if requested
    if args.run_snakemake:
        print("=" * 80)
        print("MAT-seq Analysis Pipeline - WITH SNAKEMAKE PREPROCESSING")
        print("=" * 80)
        success = pipeline.run_snakemake_preprocessing(
            snakemake_dir=args.snakemake_dir,
            fastq_dir=args.fastq_dir,
            genome_dir=args.genome_dir,
        )
        if not success:
            print("Snakemake preprocessing failed. Aborting pipeline.")
            return None

    results = pipeline.run_pipeline()

    return results


if __name__ == "__main__":
    main()
