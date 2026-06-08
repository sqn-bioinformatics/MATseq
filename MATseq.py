#!/usr/bin/env python
"""MAT-seq pipeline orchestration script."""

import argparse
import os
import subprocess
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from src import (
    CUSTOM_PALETTE_9,
    SUBSET_PALETTES,
    SUBSET_CLASS_ORDERS,
    DESEQ2_CONFIG,
    FEATURE_SELECTION_CONFIG,
    MODEL_FACTORY_CONFIG,
    MODEL_TRAINING_CONFIG,
    HYPERPARAMETER_GRIDS,
    LIGAND_ALIASES,
    prepare_counts,
    write_count_summary,
    extract_subset,
    create_feature_pipeline,
    create_preprocessing_pipeline,
    load_tlr_data,
    ModelFactory,
    ModelTrainer,
    plot_pca_pandas,
    plot_feature_count_analysis,
    plot_tlr_hek_blue,
    DESeq2,
    FeatureSelectionAnalyzer,
    VennDiagramGenerator,
    feature_count_analysis,
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
        work_dir: Path = None,
        dry_run: bool = False,
    ) -> bool:
        """Run Snakemake preprocessing pipeline to generate featurecounts."""
        from src.config import get_sample_dir, get_genome_dir, get_work_dir, get_config

        snakemake_dir = Path.cwd() / "pipeline"
        fastq_dir = fastq_dir or get_sample_dir()
        genome_dir = genome_dir or get_genome_dir()
        work_dir = work_dir or get_work_dir()

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

        threads = get_config("snakemake.threads")

        cmd = [
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
            f"WorkDir={work_dir}",
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
            key_inputs={
                "feature_selection": FEATURE_SELECTION_CONFIG,
                "n_runs": n_runs,
                "n_samples": int(len(X)),
                "labels": sorted(set(y.tolist())),
            },
        )

    def _pca_plot(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        subset_name: str,
        fs=None,
        suffix: str = "",
    ) -> None:
        if fs is None:
            pre = create_preprocessing_pipeline(
                drop_low_count=False, scale=True
            ).set_output(transform="pandas")
            X_scaled = pre.fit_transform(X)
        else:
            X_scaled = fs.transform(X)

        palette = SUBSET_PALETTES.get(subset_name, CUSTOM_PALETTE_9)
        hue_order = SUBSET_CLASS_ORDERS.get(subset_name)

        for with_names, label_suffix in [(False, ""), (True, "_labeled")]:
            plot_pca_pandas(
                X=X_scaled,
                labels=y,
                palette=palette,
                hue_order=hue_order,
                name=f"{subset_name}{suffix}{label_suffix}",
                with_sample_names=with_names,
                output_filename=f"{subset_name}{suffix}_pca{label_suffix}.png",
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
            key_inputs={
                "eval_name": eval_name,
                "hyperparameter_grids": HYPERPARAMETER_GRIDS,
                "feature_selection": FEATURE_SELECTION_CONFIG,
                "model_factory": MODEL_FACTORY_CONFIG,
                "model_training": MODEL_TRAINING_CONFIG,
                "labels": sorted(set(y.tolist())),
                "n_samples": int(len(X)),
            },
        )

    def _predict(
        self,
        predictor: ModelPredictor,
        subset_key: str,
        X: pd.DataFrame,
        y: pd.Series,
        output_dir: Path,
        subset: str = None,
        all_controls: bool = False,
    ) -> None:
        out_dir = output_dir / subset_key
        predictor.predict_samples(X, sample_names=X.index.to_numpy(), y_test=y)
        predictor.save_predictions(out_dir)
        predictor.evaluate(out_dir)
        predictor.create_probability_heatmaps(
            out_dir, subset=subset or subset_key, all_controls=all_controls
        )

    def tune_feature_count(self):
        print("\n--- FEATURE-COUNT TUNING (main_ligands) ---")
        features, labels = self.cache.cached_call(
            prepare_counts,
            name="prepare_counts",
            force_recompute=self.force_recompute,
            key_inputs={"ligand_aliases": LIGAND_ALIASES},
        )
        X, y = extract_subset(features, labels, "main_ligands")

        result = feature_count_analysis(
            X,
            y,
            k_best=FEATURE_SELECTION_CONFIG["k_best"],
            candidate_counts=[50, 100, 150, 200, 250, 300, 400, 500],
            n_estimators=FEATURE_SELECTION_CONFIG["n_estimators"],
            max_depth=FEATURE_SELECTION_CONFIG["max_depth"],
        )

        out_dir = self.results_dir / "feature_selection"
        out_dir.mkdir(parents=True, exist_ok=True)
        result["scores"].to_csv(out_dir / "feature_count_scores.csv", index=False)
        result["stability"].to_csv(
            out_dir / "feature_count_stability.csv", index=False
        )
        result["core"].to_csv(out_dir / "feature_count_core_genes.csv", index=False)

        fig_dir = self.results_dir / "figures" / "feature_selection"
        fig_path = plot_feature_count_analysis(result, output_path=fig_dir)

        print(f"MI elbow: {result['mi_elbow']}")
        print(f"ExtraTrees-importance elbow: {result['importance_elbow']}")
        print(f"Selection-stability count: {result['stable_count']}")
        print(f"Wrote scores/stability CSVs to {out_dir}")
        print(f"Wrote figure to {fig_path}")
        return result

    def run_pipeline(self):
        print("=" * 80)
        print("MAT-seq Analysis Pipeline")
        print("=" * 80)

        print("\n--- STEP 1: DATA PREPROCESSING ---")
        features, labels = self.cache.cached_call(
            prepare_counts,
            name="prepare_counts",
            force_recompute=self.force_recompute,
            key_inputs={"ligand_aliases": LIGAND_ALIASES},
        )
        write_count_summary(features, labels, self.results_dir / "counts")
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
                fs_pca = create_feature_pipeline(
                    **FEATURE_SELECTION_CONFIG
                ).set_output(transform="pandas")
                fs_pca.fit(X_sub, y_sub)
            else:
                X_main, y_main = subset_xy["main_ligands"]
                X_pca = pd.concat([X_sub, X_main])
                y_pca = pd.concat([y_sub, y_main])
            self._pca_plot(X_pca, y_pca, subset)
            self._pca_plot(X_pca, y_pca, subset, fs=fs_pca, suffix="_selected")

            deseq2 = DESeq2(
                raw_counts=X_sub,
                sample_labels=y_sub,
                padj_threshold=DESEQ2_CONFIG.get("padj_threshold", 0.05),
                log2fc_threshold=DESEQ2_CONFIG.get("log2fc_threshold", 2.0),
                n_cpus=DESEQ2_CONFIG.get("n_cpus", 42),
                cache=self.cache,
                force_recompute=self.force_recompute,
                name=subset,
            )
            deseq2.run_analysis(
                SUBSET_CLASS_ORDERS[subset], negative_control="negative_control"
            )
            deseq2_pipes[subset] = deseq2

        X_train, y_train = subset_xy["main_ligands"]
        X_other, y_other = subset_xy["additional_ligands"]
        X_bact, y_bact = subset_xy["bacteria_ligands"]
        deseq2_main = deseq2_pipes["main_ligands"]

        print("\n--- STEP 2b: VENN OF DE GENES ACROSS SUBSETS ---")
        de_sets = {s: deseq2_pipes[s].get_de_genes() for s in subset_names}
        for s, genes in de_sets.items():
            print(f"  {s}: {len(genes)} DE genes")
        shared = set.intersection(*de_sets.values())
        print(f"  shared across all three subsets: {len(shared)}")
        de_dir = self.results_dir / "differential_gene_expression"
        de_dir.mkdir(parents=True, exist_ok=True)
        for s, genes in de_sets.items():
            pd.Series(sorted(genes), name="gene").to_csv(
                de_dir / f"de_genes_{s}.csv", index=False
            )
        pd.Series(sorted(shared), name="gene").to_csv(
            de_dir / "de_genes_shared_all_subsets.csv", index=False
        )
        VennDiagramGenerator().plot_venn3(
            [de_sets[s] for s in subset_names],
            set_labels=tuple(subset_names),
            output_filename="venn_de_across_subsets.png",
            title="DE genes across ligand subsets",
        )

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

        from src.config import get_test_work_dir, get_test_name

        print(
            "\n--- STEP 4: NESTED CV TUNING + DEPLOYMENT REFIT (main_ligands, with and without Fla-PA) ---"
        )
        trainer = self._tune_subset(X_train, y_train, eval_name="main_ligands")
        trainer.save_models(self.results_dir / "models")

        flapa_mask = y_train != "Fla-PA"
        X_wo, y_wo = X_train[flapa_mask], y_train[flapa_mask]
        trainer_wo = self._tune_subset(
            X_wo,
            y_wo,
            eval_name="main_ligands_no_flapa",
            output_subdir="hyperparameter_tuning_no_flapa",
        )
        trainer_wo.save_models(self.results_dir / "models" / "no_flapa")

        print("\n--- STEP 5: MODEL VALIDATION ON EXTERNAL TEST BATCH ---")
        test_fc = get_test_work_dir() / "featurecounts"
        if not test_fc.is_dir() or not any(test_fc.glob("*.txt")):
            print(f"  Skipping validation: no featurecounts at {test_fc}")
        else:
            test_counts, test_labels = prepare_counts(featurecounts_dir=test_fc)
            Xv, yv = extract_subset(test_counts, test_labels, "main_ligands")
            deseq2_ext = DESeq2(
                raw_counts=Xv,
                sample_labels=yv,
                padj_threshold=DESEQ2_CONFIG.get("padj_threshold", 0.05),
                log2fc_threshold=DESEQ2_CONFIG.get("log2fc_threshold", 2.0),
                n_cpus=DESEQ2_CONFIG.get("n_cpus", 42),
                cache=self.cache,
                force_recompute=self.force_recompute,
                name="main_external_test",
            )
            deseq2_ext.run_analysis(
                SUBSET_CLASS_ORDERS["main_ligands"],
                negative_control="negative_control",
            )
            
            val_dir = self.results_dir / "validation" / get_test_name()
            self._predict(
                ModelPredictor(trainer),
                "main_ligands",
                Xv,
                yv,
                val_dir,
                all_controls=True,
            )
            wo = yv != "Fla-PA"
            self._predict(
                ModelPredictor(trainer_wo),
                "no_flapa",
                Xv[wo],
                yv[wo],
                val_dir,
                subset="main_ligands",
                all_controls=True,
            )

            self._pca_plot(Xv, yv, "main_ligands", suffix="_external_test")
            self._pca_plot(
                Xv, yv, "main_ligands", fs=fs_pca, suffix="_external_test_selected"
            )



        print("\n--- STEP 6: CLASS PREDICTION ON ADDITIONAL AND BACTERIA LIGANDS ---")
        predictor = ModelPredictor(trainer)
        pred_dir = self.results_dir / "predictions"
        self._predict(predictor, "additional_ligands", X_other, y_other, pred_dir)
        self._predict(predictor, "bacteria_ligands", X_bact, y_bact, pred_dir)

        print("\n--- STEP 7: TLR HEK BLUE VISUALIZATION ---")
        tlr2_df, tlr4_df, flapa_data = load_tlr_data(
            data_dir=Path(__file__).parent / "data" / "supplementary_data"
        )
        plot_tlr_hek_blue(
            tlr2_df, tlr4_df, flapa_data, output_filename="tlr_hek_blue.png"
        )

        print(
            "\n--- STEP 8: CLASS PREDICTION ON ADDITIONAL AND BACTERIA LIGANDS (without Fla-PA) ---"
        )
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
    parser.add_argument(
        "--tune-features",
        action="store_true",
        help="Run model-independent feature-count analysis and exit",
    )

    args = parser.parse_args()
    pipeline = MATseqPipeline(
        cache_dir=args.cache_dir, force_recompute=args.force_recompute
    )

    if args.tune_features:
        return pipeline.tune_feature_count()

    if args.snakemake:
        from src.config import get_test_sample_dir, get_test_work_dir

        dry_run = args.snakemake == "dry-run"
        batches = [
            (args.fastq_dir, None),
            (get_test_sample_dir(), get_test_work_dir()),
        ]
        for fastq_dir, work_dir in batches:
            ok = pipeline.run_snakemake_preprocessing(
                fastq_dir=fastq_dir,
                work_dir=work_dir,
                genome_dir=args.genome_dir,
                dry_run=dry_run,
            )
            if not ok:
                return None
        if dry_run:
            return None

    return pipeline.run_pipeline()


if __name__ == "__main__":
    main()
