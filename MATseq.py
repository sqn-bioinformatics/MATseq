#!/usr/bin/env python
"""MAT-seq pipeline orchestration script."""

import sys
from pathlib import Path
import argparse
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).parent))

from src import (
    # Configuration
    CUSTOM_PALETTE_9,
    SUBSET_PALETTES,
    SUBSET_CLASS_ORDERS,
    TRAINING_LIGANDS,
    ADDITIONAL_LIGANDS,
    BACTERIAL_LIGANDS,
    TRAINING_LIGANDS_WO_FLAPA,
    DESEQ2_CONFIG,
    FEATURE_SELECTION_CONFIG,
    MODEL_FACTORY_CONFIG,
    MODEL_TRAINING_CONFIG,
    # Preprocessing
    merge_counts,
    extract_subset,
    normalize_rpm,
    load_tlr_data,
    # Feature engineering
    create_feature_pipeline,
    # Model training
    ModelFactory,
    ModelTrainer,
    # Evaluation
    make_score,
    get_confusion_matrix,
    # Visualization
    plot_pca_for_pandas,
    plot_gene_expression_by_class,
    plot_tlr_hek_blue,
    make_probability_matrix,
    # DESeq2 analysis
    AnalysisPipeline,
    # Feature analysis
    FeatureSelectionAnalyzer,
    VennDiagramGenerator,
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
        fastq_dir: Path = None,
        genome_dir: Path = None,
        dry_run: bool = False,
    ) -> bool:
        """Run Snakemake preprocessing pipeline to generate featurecounts.

        Args:
            fastq_dir: Override FASTQ input directory (default: from config.json).
            genome_dir: Override genome reference directory (required).
            dry_run: If True, only validate Snakemake without running.

        Returns:
            True if successful, False otherwise.
        """
        import subprocess
        import json

        snakemake_dir = Path.cwd() / "pipeline"

        config_path = Path.cwd() / "config.json"
        with open(config_path) as f:
            config = json.load(f)

        snakemake_config = config.get("snakemake", {})

        if fastq_dir is None:
            fastq_dir = snakemake_config.get("sample_dir", "")
            if fastq_dir:
                fastq_dir = Path(fastq_dir).expanduser()

        if genome_dir is None:
            genome_dir = snakemake_config.get("genome_dir", "")
            if genome_dir:
                genome_dir = Path(genome_dir).expanduser()

        if not genome_dir or not Path(genome_dir).exists():
            print("Error: genome_dir is required and must exist")
            print(f"Current value: '{genome_dir}'")
            print("Set via --genome-dir or in config.json snakemake.genome_dir")
            return False

        if not fastq_dir or not Path(fastq_dir).exists():
            print("Error: sample_dir (FASTQ directory) is required and must exist")
            print(f"Current value: '{fastq_dir}'")
            print("Set via --fastq-dir or in config.json snakemake.sample_dir")
            return False

        # Write validated paths back to config
        config["snakemake"]["sample_dir"] = str(fastq_dir)
        config["snakemake"]["genome_dir"] = str(genome_dir)
        with open(config_path, "w") as f:
            json.dump(config, f, indent=2)

        snakefile = snakemake_dir / "0_MATseq.smk"

        if not snakefile.exists():
            print(f"Error: Snakemake pipeline not found at {snakefile}")
            return False

        print("\n--- RUNNING SNAKEMAKE PREPROCESSING ---")
        print(f"Pipeline directory: {snakemake_dir}")
        print(f"FASTQ directory: {fastq_dir}")
        print(f"Genome directory: {genome_dir}")

        work_dir = Path(
            snakemake_config.get("work_dir", "~/MATseq/pipeline/results")
        ).expanduser()
        threads = snakemake_config.get("threads", 42)

        cmd = [
            "poetry",
            "run",
            "snakemake",
            "--use-conda",
            "--cores",
            str(threads),
            "--snakefile",
            str(snakefile),
            "--directory",
            str(work_dir),
            "--config",
            f"SampleDir={fastq_dir}",
            f"GenomeDir={genome_dir}",
            "--latency-wait",
            "60",
        ]

        if dry_run:
            cmd.append("--dry-run")
            print("Dry run mode: validating pipeline...")

        try:
            result = subprocess.run(
                cmd,
                cwd=str(Path.cwd()),
                capture_output=False,
            )

            if result.returncode == 0:
                print("Snakemake preprocessing completed successfully")
                return True
            else:
                print(f"Snakemake preprocessing failed with code {result.returncode}")
                return False
        except FileNotFoundError:
            print(
                "Error: poetry or snakemake not found. Install with: poetry add snakemake"
            )
            return False
        except Exception as e:
            print(f"Error running Snakemake: {e}")
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
        self, X: pd.DataFrame, y: pd.Series, de_genes: set, n_runs: int = 1000
    ) -> tuple:
        """Run feature selection multiple times to see which genes are more prevalent.

        Args:
            X: Feature matrix.
            y: Target labels.
            n_runs: Number of selection runs.
            de_genes: Set of differentially expressed genes for frequency table.

        Returns:
            Tuple of (analyzer, fs_genes_set).
        """

        def _venn_selection():
            print(f"Running feature selection {n_runs} times...")

            analyzer = FeatureSelectionAnalyzer()

            analyzer.run_multiple_selections(
                X=X,
                y=y,
                n_runs=n_runs,
                pipeline_config=FEATURE_SELECTION_CONFIG,
            )

            analyzer.create_gene_frequency_table(de_genes=de_genes)

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
        X = subset_data.drop(columns="label")
        y = subset_data["label"]

        # Concat and standardization here is only for PCA plotting
        print(f"Creating PCA plots for {subset_name}...")
        if subset_name != "training":
            training_data = extract_subset(df, name="training")
            X = pd.concat([X, training_data.drop(columns="label")])
            y = pd.concat([y, training_data["label"]])

        X_rpm = normalize_rpm(X.copy())
        ss = StandardScaler()
        ss.fit(X_rpm)
        ss.set_output(transform="pandas")
        X_scaled = ss.transform(X_rpm)

        palette = SUBSET_PALETTES.get(subset_name, CUSTOM_PALETTE_9)
        hue_order = SUBSET_CLASS_ORDERS.get(subset_name)

        common_kwargs = dict(
            X=X_scaled,
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
                output_filename=f"{subset_name}_pca{suffix}.png",
                **common_kwargs,
            )

        return subset_data.drop(columns="label"), subset_data["label"]

    def run_feature_selection(
        self, X: pd.DataFrame, y: pd.Series, subset_name: str
    ) -> tuple[object, pd.DataFrame]:
        """Run feature selection pipeline.

        Args:
            X: Feature matrix.
            y: Labels.
            subset_name: Name of subset.

        Returns:
            Tuple of (pipeline, selected_features_df).
        """

        def _feature_select():
            print(f"Running feature selection for {subset_name}...")
            pipe = create_feature_pipeline(**FEATURE_SELECTION_CONFIG).set_output(
                transform="pandas"
            )
            X_selected = pipe.fit_transform(X, y)
            feature_names = pipe[:-1].get_feature_names_out()
            X_fs = pd.DataFrame(X_selected, columns=feature_names, index=X.index)

            return pipe, X_fs

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
            X=X_selected,
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
            models = ModelFactory.create_models(**MODEL_FACTORY_CONFIG)
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
            cache=self.cache,
            force_recompute=self.force_recompute,
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
            X_train, y_train, n_runs=1000, de_genes=de_genes
        )
        venn_gen = VennDiagramGenerator()
        venn_gen.plot_venn(
            de_genes,
            fs_genes,
            output_filename="venn_de_vs_fs.png",
        )

        # Step 4: Model Training
        print("\n--- STEP 4: MODEL TRAINING ---")
        pipe, X_fs = self.run_feature_selection(
            X_train,
            y_train,
            "training",
        )

        trainer = self.train_models(X_fs, y_train, "training")
        models_dir = self.results_dir / "models"
        trainer.save_models(models_dir)

        print("\n--- EVALUATING MODELS ---")
        eval_dir = self.results_dir / "model_evaluation"

        # Evaluation 1: Full feature set
        print("Evaluation on full feature set (X_train)...")
        trainer.evaluate(X_train, y_train, eval_dir=eval_dir, cv=5, eval_name="full")

        # Evaluation 2: Feature-selected genes
        print("Evaluation on feature-selected genes (X_fs)...")
        trainer.evaluate(X_fs, y_train, eval_dir=eval_dir, cv=5, eval_name="fs")

        # Evaluation 3: X_train subsetted to genes in X_fs ∪ de_genes
        print("X_fs:", X_fs)
        fs_genes = set(X_fs.columns)
        print(f"Feature Selections genes: {fs_genes}")
        print(f"Differentially expressed genes: {de_genes}")
        union_genes = list(fs_genes | de_genes)
        print("Union genes:", union_genes)

        X_fs_de = X_train[union_genes]
        print(f"Evaluation on FS ∪ DE genes ({len(union_genes)} genes)...")
        trainer.evaluate(X_fs_de, y_train, eval_dir=eval_dir, cv=5, eval_name="fs_de")

        print("Model evaluation complete. Results saved to results/model_evaluation")

        # Step 5: DESeq2 Analysis (on other ligands and bacterial subsets)
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
            cache=self.cache,
            force_recompute=self.force_recompute,
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
            cache=self.cache,
            force_recompute=self.force_recompute,
        )
        deseq2_pipeline_bacterial.run_analysis(
            BACTERIAL_LIGANDS,
            negative_control="negative_control",
        )

        # Step 6: Predict classes of other ligands
        print("\n--- STEP 6: CLASS PREDICTION ON ADDITIONAL AND BACTERIAL LIGANDS ---")

        # Added LPS for plotting
        lps_mask_train = y_train == "LPS"
        X_train_lps = X_train[lps_mask_train]
        y_train_lps = y_train[lps_mask_train]
        X_other_with_lps = pd.concat([X_other_ligands, X_train_lps])
        y_other_with_lps = pd.concat([y_other_ligands, y_train_lps])
        X_bacterial_with_lps = pd.concat([X_bacterial, X_train_lps])
        y_bacterial_with_lps = pd.concat([y_bacterial, y_train_lps])

        predictor = ModelPredictor(trainer)
        predictions_dir = self.results_dir / "predictions"

        X_other_fs = pipe.transform(X_other_with_lps)
        predictor.predict_samples(
            X_other_fs,
            sample_names=X_other_with_lps.index.to_numpy(),
            y_test=y_other_with_lps,
        )
        predictor.save_predictions(predictions_dir / "other_ligands")
        predictor.create_probability_heatmaps(predictions_dir / "other_ligands")

        X_bacterial_fs = pipe.transform(X_bacterial_with_lps)
        predictor.predict_samples(
            X_bacterial_fs,
            sample_names=X_bacterial_with_lps.index.to_numpy(),
            y_test=y_bacterial_with_lps,
        )
        predictor.save_predictions(predictions_dir / "bacterial")
        predictor.create_probability_heatmaps(predictions_dir / "bacterial")

        # Step 7: Model Training on training_wo_flapa and Prediction on Additional and Bacterial Ligands
        print(
            "\n--- STEP 7: MODEL TRAINING ON TRAINING_WO_FLAPA AND PREDICTION ON ADDITIONAL AND BACTERIAL LIGANDS ---"
        )
        X_train_wo_flapa, y_train_wo_flapa = self.extract_and_run_pca_before_pipeline(
            df, "training_wo_flapa"
        )
        pipe_wo_flapa, X_fs_wo_flapa = self.run_feature_selection(
            X_train_wo_flapa,
            y_train_wo_flapa,
            "training_wo_flapa",
        )
        trainer_wo_flapa = self.train_models(
            X_fs_wo_flapa, y_train_wo_flapa, "training_wo_flapa"
        )

        predictor_wo_flapa = ModelPredictor(trainer_wo_flapa)
        predictions_dir_wo_flapa = (
            self.results_dir / "predictions" / "training_wo_flapa"
        )

        X_other_fs_wo_flapa = pipe_wo_flapa.transform(X_other_with_lps)
        predictor_wo_flapa.predict_samples(
            X_other_fs_wo_flapa,
            sample_names=X_other_with_lps.index.to_numpy(),
            y_test=y_other_with_lps,
        )
        predictor_wo_flapa.save_predictions(predictions_dir_wo_flapa / "other_ligands")
        predictor_wo_flapa.create_probability_heatmaps(
            predictions_dir_wo_flapa / "other_ligands"
        )

        X_bacterial_fs_wo_flapa = pipe_wo_flapa.transform(X_bacterial_with_lps)
        predictor_wo_flapa.predict_samples(
            X_bacterial_fs_wo_flapa,
            sample_names=X_bacterial_with_lps.index.to_numpy(),
            y_test=y_bacterial_with_lps,
        )
        predictor_wo_flapa.save_predictions(predictions_dir_wo_flapa / "bacterial")
        predictor_wo_flapa.create_probability_heatmaps(
            predictions_dir_wo_flapa / "bacterial"
        )

        # Step 8: TLR HEK visualization
        print("\n--- STEP 8: TLR HEK BLUE VISUALIZATION ---")
        tlr2_df, tlr4_df, flapa_data = load_tlr_data()

        plot_tlr_hek_blue(
            tlr2_df,
            tlr4_df,
            flapa_data,
            output_filename="tlr_hek_blue.png",
        )

        print("\n" + "=" * 80)
        print("PIPELINE COMPLETED SUCCESSFULLY")
        print("=" * 80)
        print(f"Results saved to: {self.results_dir.absolute()}")
        print(f"Cache saved to: {self.cache.cache_dir.absolute()}")

        return {
            "feature_selection": pipe,
            "feature_analysis": fs_analyzer,
            "models": trainer,
            "predictor": predictor,
            "models_wo_flapa": trainer_wo_flapa,
            "predictor_wo_flapa": predictor_wo_flapa,
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
    parser.add_argument(
        "--snakemake-dry-run",
        action="store_true",
        help="Validate Snakemake pipeline without running (dry run)",
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
            fastq_dir=args.fastq_dir,
            genome_dir=args.genome_dir,
            dry_run=args.snakemake_dry_run,
        )
        if not success:
            print("Snakemake preprocessing failed. Aborting pipeline.")
            return None

    results = pipeline.run_pipeline()

    return results


if __name__ == "__main__":
    main()
