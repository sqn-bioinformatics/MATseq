#!/usr/bin/env python
"""MAT-seq pipeline orchestration script."""

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).parent))

from src import (
    CUSTOM_PALETTE_9,
    SUBSET_PALETTES,
    SUBSET_CLASS_ORDERS,
    DESEQ2_CONFIG,
    FEATURE_SELECTION_CONFIG,
    MODEL_FACTORY_CONFIG,
    MODEL_TRAINING_CONFIG,
    prepare_counts,
    extract_subset,
    normalize_rpm,
    load_tlr_data,
    create_feature_pipeline,
    ModelFactory,
    ModelTrainer,
    plot_pca_pandas,
    plot_tlr_hek_blue,
    DESeq2,
    FeatureSelectionAnalyzer,
    PipelineParamTuner,
    VennDiagramGenerator,
    create_fs_de_go_table,
    ModelPredictor,
)
from src.cache import PipelineCache


class MATseqPipeline:
    """MAT-seq analysis pipeline with disk caching."""

    def __init__(self, cache_dir: Path = None, force_recompute: bool = False):
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
        """Run Snakemake preprocessing pipeline to generate featurecounts."""
        from src.config import get_sample_dir, get_genome_dir, get_work_dir, get_config

        snakemake_dir = Path.cwd() / "pipeline"
        fastq_dir = fastq_dir or get_sample_dir()
        genome_dir = genome_dir or get_genome_dir()

        if not genome_dir or not Path(genome_dir).exists():
            print(
                f"Error: genome_dir is required and must exist. Current: '{genome_dir}'"
            )
            return False
        if not fastq_dir or not Path(fastq_dir).exists():
            print(
                f"Error: sample_dir is required and must exist. Current: '{fastq_dir}'"
            )
            return False

        snakefile = snakemake_dir / "0_MATseq.smk"
        if not snakefile.exists():
            print(f"Error: Snakemake pipeline not found at {snakefile}")
            return False

        work_dir = get_work_dir()
        threads = get_config("snakemake.threads")

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
        ]
        if dry_run:
            cmd.append("--dry-run")

        print(
            f"\n--- RUNNING SNAKEMAKE PREPROCESSING ---\nFASTQ: {fastq_dir}\nGenome: {genome_dir}"
        )
        result = subprocess.run(cmd, cwd=str(Path.cwd()))
        if result.returncode == 0:
            print("Snakemake preprocessing completed successfully")
            return True
        print(f"Snakemake preprocessing failed with code {result.returncode}")
        return False

    def _venn_feature_selection(
        self, X: pd.DataFrame, y: pd.Series, de_genes: set, n_runs: int = 1000
    ):
        def _run():
            analyzer = FeatureSelectionAnalyzer()
            analyzer.run_multiple_selections(X=X, y=y, n_runs=n_runs)
            analyzer.create_gene_frequency_table(de_genes=de_genes)
            return analyzer, set().union(*analyzer.feature_sets)

        return self.cache.cached_call(
            _run,
            name=f"venn_feature_selection_{n_runs}",
            force_recompute=self.force_recompute,
        )

    def _pca_plot(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        subset_name: str,
        file_suffix: str = "",
        tag: str | None = None,
        prescaled: bool = False,
    ) -> None:
        if prescaled:
            X_scaled = X
        else:
            X_rpm = normalize_rpm(X)
            scaler = StandardScaler().set_output(transform="pandas")
            X_scaled = scaler.fit_transform(X_rpm)

        palette = SUBSET_PALETTES.get(subset_name, CUSTOM_PALETTE_9)
        hue_order = SUBSET_CLASS_ORDERS.get(subset_name)
        label_base = tag or subset_name

        for with_names, label_suffix in [(False, ""), (True, "_labeled")]:
            plot_pca_pandas(
                X=X_scaled,
                labels=y,
                palette=palette,
                hue_order=hue_order,
                name=f"{label_base}{file_suffix}{label_suffix}",
                with_sample_names=with_names,
                output_filename=f"{label_base}_pca{file_suffix}{label_suffix}.png",
            )

    def _fs_and_train(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        subset_name: str,
        tag: str | None = None,
    ) -> tuple[object, pd.DataFrame, ModelTrainer]:
        cache_key = tag or subset_name

        def _fit_fs():
            pipe = create_feature_pipeline(
                **FEATURE_SELECTION_CONFIG,
                random_state=MODEL_TRAINING_CONFIG["random_state"],
            ).set_output(transform="pandas")
            X_selected = pipe.fit_transform(X, y)
            feature_names = pipe[:-1].get_feature_names_out()
            return pipe, pd.DataFrame(X_selected, columns=feature_names, index=X.index)

        pipe, X_fs = self.cache.cached_call(
            _fit_fs,
            name=f"feature_selection_{cache_key}",
            force_recompute=self.force_recompute,
        )

        self._pca_plot(
            X_fs, y, subset_name, file_suffix="_selected", tag=tag, prescaled=True
        )

        def _fit_models():
            models = ModelFactory.create_models(**MODEL_FACTORY_CONFIG)
            trainer = ModelTrainer(X_fs, y, models=models, **MODEL_TRAINING_CONFIG)
            trainer.train_all_models()
            return trainer

        trainer = self.cache.cached_call(
            _fit_models,
            name=f"trained_models_{cache_key}",
            force_recompute=self.force_recompute,
        )
        return pipe, X_fs, trainer

    def _predict(
        self,
        predictor: ModelPredictor,
        pipe,
        subset_key: str,
        X: pd.DataFrame,
        y: pd.Series,
        output_dir: Path,
    ) -> None:
        X_fs = pipe.transform(X)
        predictor.predict_samples(X_fs, sample_names=X.index.to_numpy(), y_test=y)
        predictor.save_predictions(output_dir / subset_key)
        predictor.create_probability_heatmaps(
            output_dir / subset_key, subset=subset_key
        )

    def run_pipeline(self):
        print("=" * 80)
        print("MAT-seq Analysis Pipeline")
        print("=" * 80)

        print("\n--- STEP 1: DATA PREPROCESSING ---")
        features, labels = self.cache.cached_call(
            prepare_counts, name="prepare_counts", force_recompute=self.force_recompute
        )
        print(f"Count df shape: {features.shape}")

        print("\n--- STEP 2: DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS ---")
        subset_names = ["main_ligands", "additional_ligands", "bacteria_ligands"]
        subset_xy: dict[str, tuple[pd.DataFrame, pd.Series]] = {}
        deseq2_pipes: dict[str, DESeq2] = {}

        for subset in subset_names:
            print(f"\nProcessing {subset} subset...")
            X_sub, y_sub = extract_subset(features, labels, subset)
            subset_xy[subset] = (X_sub, y_sub)

            if subset == "main_ligands":
                X_pca, y_pca = X_sub, y_sub
            else:
                X_main, y_main = subset_xy["main_ligands"]
                X_pca = pd.concat([X_sub, X_main])
                y_pca = pd.concat([y_sub, y_main])
            self._pca_plot(X_pca, y_pca, subset)

            deseq2 = DESeq2(
                raw_counts=X_sub,
                sample_labels=y_sub,
                padj_threshold=DESEQ2_CONFIG.get("padj_threshold", 0.05),
                log2fc_threshold=DESEQ2_CONFIG.get("log2fc_threshold", 2.0),
                baseMean_threshold=DESEQ2_CONFIG.get("baseMean_threshold", 10.0),
                n_cpus=DESEQ2_CONFIG.get("n_cpus", 42),
                cache=self.cache,
                force_recompute=self.force_recompute,
            )
            deseq2.run_analysis(
                SUBSET_CLASS_ORDERS[subset], negative_control="negative_control"
            )
            deseq2_pipes[subset] = deseq2

        X_train, y_train = subset_xy["main_ligands"]
        X_other, y_other = subset_xy["additional_ligands"]
        X_bact, y_bact = subset_xy["bacteria_ligands"]
        deseq2_main = deseq2_pipes["main_ligands"]

        print("\n--- STEP 2b: PIPELINE PARAMETER TUNING ---")
        tuner = PipelineParamTuner()
        tuner.run(
            X_train,
            y_train,
            k_best_values=[100, 500, 1000, 2000],
            max_features_values=[50, 100, 250, 500],
        )
        tuner.plot()

        print("\n--- STEP 3: VENN OF FS vs DE GENES (training subset) ---")
        de_genes = deseq2_main.get_de_genes()
        fs_analyzer, fs_genes = self._venn_feature_selection(
            X_train, y_train, de_genes=de_genes, n_runs=1000
        )
        VennDiagramGenerator().plot_venn(
            de_genes, fs_genes, output_filename="venn_de_vs_fs.png"
        )
        goeaobj, geneid_symbol_mapper = deseq2_main.get_go_objects()
        create_fs_de_go_table(
            de_genes=de_genes,
            fs_genes=fs_genes,
            goeaobj=goeaobj,
            geneid_symbol_mapper=geneid_symbol_mapper,
        )

        print("\n--- STEP 4: MODEL TRAINING (main_ligands) ---")
        pipe, _, trainer = self._fs_and_train(X_train, y_train, "main_ligands")
        trainer.save_models(self.results_dir / "models")

        print("\n--- EVALUATING MODELS (per-fold FS via ModelTrainer.evaluate) ---")
        trainer.evaluate(
            X_train,
            y_train,
            eval_dir=self.results_dir / "model_evaluation",
            cv=10,
            eval_name="main_ligands",
        )

        print("\n--- STEP 5: CLASS PREDICTION ON ADDITIONAL AND BACTERIA LIGANDS ---")
        # Append training LPS samples as positive control for the held-out ligand sets.
        predictor = ModelPredictor(trainer)
        pred_dir = self.results_dir / "predictions"
        self._predict(predictor, pipe, "additional_ligands", X_other, y_other, pred_dir)
        self._predict(predictor, pipe, "bacteria_ligands", X_bact, y_bact, pred_dir)

        print("\n--- STEP 6: TLR HEK BLUE VISUALIZATION ---")
        tlr2_df, tlr4_df, flapa_data = load_tlr_data(
            data_dir=Path(__file__).parent / "data" / "supplementary_data"
        )
        plot_tlr_hek_blue(
            tlr2_df, tlr4_df, flapa_data, output_filename="tlr_hek_blue.png"
        )

        print("\n--- STEP 7: MODEL TRAINING (main_ligands without Fla-PA) ---")
        flapa_mask = y_train != "Fla-PA"
        X_wo, y_wo = X_train[flapa_mask], y_train[flapa_mask]
        pipe_wo, _, trainer_wo = self._fs_and_train(
            X_wo, y_wo, "main_ligands", tag="main_ligands_no_flapa"
        )

        predictor_wo = ModelPredictor(trainer_wo)
        pred_dir_wo = self.results_dir / "predictions" / "main_ligands_no_flapa"
        self._predict(
            predictor_wo,
            pipe_wo,
            "additional_ligands",
            X_other,
            y_other,
            pred_dir_wo,
        )
        self._predict(
            predictor_wo,
            pipe_wo,
            "bacteria_ligands",
            X_bact,
            y_bact,
            pred_dir_wo,
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
            "models_wo_flapa": trainer_wo,
            "predictor_wo_flapa": predictor_wo,
        }


def main():
    parser = argparse.ArgumentParser(
        description="MAT-seq analysis pipeline with caching and full analysis"
    )
    parser.add_argument(
        "--force-recompute",
        action="store_true",
        help="Skip cache and recompute all steps",
    )
    parser.add_argument(
        "--cache-dir", type=Path, default=None, help="Directory for cached files"
    )
    parser.add_argument(
        "--fastq-dir", type=Path, default=None, help="Override FASTQ input directory"
    )
    parser.add_argument(
        "--genome-dir",
        type=Path,
        default=None,
        help="Override genome reference directory",
    )
    parser.add_argument(
        "--snakemake",
        choices=["dry-run", "run"],
        default=None,
        help="Run Snakemake preprocessing: 'dry-run' to validate, 'run' to execute",
    )

    args = parser.parse_args()
    pipeline = MATseqPipeline(
        cache_dir=args.cache_dir, force_recompute=args.force_recompute
    )

    if args.snakemake:
        ok = pipeline.run_snakemake_preprocessing(
            fastq_dir=args.fastq_dir,
            genome_dir=args.genome_dir,
            dry_run=args.snakemake == "dry-run",
        )
        if not ok or args.snakemake == "dry-run":
            return None

    return pipeline.run_pipeline()


if __name__ == "__main__":
    main()
