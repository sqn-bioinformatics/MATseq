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
    MODEL_FACTORY_CONFIG,
    MODEL_TRAINING_CONFIG,
    HYPERPARAMETER_GRIDS,
    prepare_counts,
    extract_subset,
    normalize_rpm,
    load_tlr_data,
    ModelFactory,
    ModelTrainer,
    plot_pca_pandas,
    plot_tlr_hek_blue,
    DESeq2,
    FeatureSelectionAnalyzer,
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
            "--rerun-incomplete",
        ]
        if dry_run:
            cmd.append("--dry-run")

        print(
            f"\n--- RUNNING SNAKEMAKE PREPROCESSING ---\nFASTQ: {fastq_dir}\nGenome: {genome_dir}"
        )
        subprocess.run(cmd + ["--unlock"], cwd=str(Path.cwd()))
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
    ) -> None:
        X_rpm = normalize_rpm(X)
        scaler = StandardScaler().set_output(transform="pandas")
        X_scaled = scaler.fit_transform(X_rpm)

        palette = SUBSET_PALETTES.get(subset_name, CUSTOM_PALETTE_9)
        hue_order = SUBSET_CLASS_ORDERS.get(subset_name)

        for with_names, label_suffix in [(False, ""), (True, "_labeled")]:
            plot_pca_pandas(
                X=X_scaled,
                labels=y,
                palette=palette,
                hue_order=hue_order,
                name=f"{subset_name}{label_suffix}",
                with_sample_names=with_names,
                output_filename=f"{subset_name}_pca{label_suffix}.png",
            )

    def _tune_subset(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        eval_name: str,
        output_subdir: str = "hyperparameter_tuning",
    ) -> ModelTrainer:
        def _run():
            models = ModelFactory.create_models(**MODEL_FACTORY_CONFIG)
            trainer = ModelTrainer(X, y, models=models, **MODEL_TRAINING_CONFIG)
            trainer.tune_nested(
                X,
                y,
                param_grids=HYPERPARAMETER_GRIDS,
                output_dir=self.results_dir / output_subdir,
                outer_cv=5,
                inner_cv=3,
                eval_name=eval_name,
            )
            return trainer

        return self.cache.cached_call(
            _run,
            name=f"tuned_models_{eval_name}",
            force_recompute=self.force_recompute,
        )

    def _predict(
        self,
        predictor: ModelPredictor,
        subset_key: str,
        X: pd.DataFrame,
        y: pd.Series,
        output_dir: Path,
    ) -> None:
        predictor.predict_samples(X, sample_names=X.index.to_numpy(), y_test=y)
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

        print("\n--- STEP 4: NESTED CV TUNING + DEPLOYMENT REFIT (main_ligands) ---")
        trainer = self._tune_subset(X_train, y_train, eval_name="main_ligands")
        trainer.save_models(self.results_dir / "models")

        print("\n--- STEP 5: CLASS PREDICTION ON ADDITIONAL AND BACTERIA LIGANDS ---")
        predictor = ModelPredictor(trainer)
        pred_dir = self.results_dir / "predictions"
        self._predict(predictor, "additional_ligands", X_other, y_other, pred_dir)
        self._predict(predictor, "bacteria_ligands", X_bact, y_bact, pred_dir)

        print("\n--- STEP 6: TLR HEK BLUE VISUALIZATION ---")
        tlr2_df, tlr4_df, flapa_data = load_tlr_data(
            data_dir=Path(__file__).parent / "data" / "supplementary_data"
        )
        plot_tlr_hek_blue(
            tlr2_df, tlr4_df, flapa_data, output_filename="tlr_hek_blue.png"
        )

        print(
            "\n--- STEP 7: NESTED CV TUNING + DEPLOYMENT REFIT (main_ligands without Fla-PA) ---"
        )
        flapa_mask = y_train != "Fla-PA"
        X_wo, y_wo = X_train[flapa_mask], y_train[flapa_mask]
        trainer_wo = self._tune_subset(
            X_wo,
            y_wo,
            eval_name="main_ligands_no_flapa",
            output_subdir="hyperparameter_tuning_no_flapa",
        )
        trainer_wo.save_models(self.results_dir / "models" / "no_flapa")

        predictor_wo = ModelPredictor(trainer_wo)
        pred_dir_wo = self.results_dir / "predictions" / "main_ligands_no_flapa"
        self._predict(predictor_wo, "additional_ligands", X_other, y_other, pred_dir_wo)
        self._predict(predictor_wo, "bacteria_ligands", X_bact, y_bact, pred_dir_wo)

        print("\n" + "=" * 80)
        print("PIPELINE COMPLETED SUCCESSFULLY")
        print("=" * 80)
        print(f"Results saved to: {self.results_dir.absolute()}")
        print(f"Cache saved to: {self.cache.cache_dir.absolute()}")

        return {
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
