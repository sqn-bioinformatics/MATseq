#!/usr/bin/env python
"""MAT-seq pipeline orchestration script."""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.pipeline import Pipeline as SkPipeline
from sklearn.utils.class_weight import compute_sample_weight

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
    ColumnSelector,
    FeatureSelector,
    load_tlr_data,
    ModelFactory,
    ModelTrainer,
    plot_pca_pandas,
    plot_feature_count_analysis,
    plot_venn,
    plot_tlr_hek_blue,
    DESeq2,
    create_fs_de_go_table,
    ModelPredictor,
)
from src.config import FOREST_SELECTION_GRID, update_config, primary_geneset_name
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

    def _select_features(self, X: pd.DataFrame, y: pd.Series) -> list:
        """Fit the tuned feature pipeline once and return its selected gene list."""

        def _run():
            fs = FeatureSelector.feature_pipeline(
                **FEATURE_SELECTION_CONFIG
            ).set_output(transform="pandas")
            fs.fit(X, y)
            return list(fs.get_feature_names_out())

        return self.cache.cached_call(
            _run,
            name="selected_features",
            force_recompute=self.force_recompute,
            key_inputs={
                "feature_selection": FEATURE_SELECTION_CONFIG,
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
            pre = FeatureSelector.preprocessing_pipeline().set_output(
                transform="pandas"
            )
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
                equal_aspect=True,
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

    @staticmethod
    def _best_two_models(summary: pd.DataFrame, rank_col: str = "pooled_f1") -> list:
        col = rank_col if rank_col in summary.columns else "f1_mean"
        return summary.sort_values(col, ascending=False)["model"].head(2).tolist()

    def _deploy_geneset_trainer(
        self,
        X_train: pd.DataFrame,
        y_train: pd.Series,
        gene_list: list,
        model_names: list,
        tuned_params: dict,
        label_encoder,
    ) -> ModelTrainer:
        """Refit the chosen models on a frozen gene list with their tuned params."""
        base = ModelFactory.create_models(**MODEL_FACTORY_CONFIG)
        y_array = y_train.values if isinstance(y_train, pd.Series) else y_train
        y_enc = label_encoder.transform(y_array)

        trained = {}
        for name in model_names:
            pipe = SkPipeline(
                [
                    *FeatureSelector.preprocessing_pipeline().steps,
                    ("select_genes", ColumnSelector(gene_list)),
                    ("clf", clone(base[name])),
                ]
            ).set_output(transform="pandas")
            params = tuned_params.get(name, {}).get("chosen_params", {})
            if params:
                pipe.set_params(**params)
            fit_kwargs = (
                {"clf__sample_weight": compute_sample_weight("balanced", y_enc)}
                if name == "XGBoost"
                else {}
            )
            pipe.fit(X_train, y_enc, **fit_kwargs)
            trained[name] = pipe

        return ModelTrainer.from_trained_models(trained, label_encoder)

    def _endpoint_predictions(
        self, trainers: dict, base_dir: Path, X_other, y_other, X_bact, y_bact
    ) -> None:
        for gs in (primary_geneset_name(), "de_overlap"):
            predictor = ModelPredictor(trainers[gs])
            self._predict(predictor, "additional_ligands", X_other, y_other, base_dir / gs)
            self._predict(predictor, "bacteria_ligands", X_bact, y_bact, base_dir / gs)

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
        predictor.evaluate(out_dir, subset=subset or subset_key)
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

        result = FeatureSelector.count_analysis(X, y)

        out_dir = self.results_dir / "feature_selection"
        out_dir.mkdir(parents=True, exist_ok=True)
        result["scores"].to_csv(out_dir / "feature_count_scores.csv", index=False)

        fig_dir = self.results_dir / "figures" / "feature_selection"
        fig_path = plot_feature_count_analysis(result, output_path=fig_dir)

        print(f"Per-seed MI elbows: {result['per_run_elbows']}")
        print(f"Mean MI elbow (k_best): {result['mi_elbow']}")
        print(f"Wrote scores CSV to {out_dir}")
        print(f"Wrote figure to {fig_path}")
        return result

    def tune_forest_features(self):
        print("\n--- FOREST FEATURE-COUNT TUNING BY KMEANS SEPARATION (main_ligands) ---")
        features, labels = self.cache.cached_call(
            prepare_counts,
            name="prepare_counts",
            force_recompute=self.force_recompute,
            key_inputs={"ligand_aliases": LIGAND_ALIASES},
        )
        X, y = extract_subset(features, labels, "main_ligands")

        result = FeatureSelector.forest_kmeans_separation(
            X,
            y,
            k_best=FEATURE_SELECTION_CONFIG["k_best"],
            grid=FOREST_SELECTION_GRID,
        )

        out_dir = self.results_dir / "feature_selection"
        out_dir.mkdir(parents=True, exist_ok=True)
        result["scan"].to_csv(out_dir / "forest_kmeans_separation.csv", index=False)

        best = result["best"]
        print(result["scan"].to_string(index=False))
        print(
            f"Best by ARI: n_estimators={best['n_estimators']}, "
            f"max_depth={best['max_depth']}, n_selected={best['n_selected']} "
            f"(ARI={best['ari']:.3f})"
        )
        print(f"Wrote scan to {out_dir / 'forest_kmeans_separation.csv'}")
        return result

    def replot(self):
        """Regenerate figures from cache + saved CSVs (NO model recompute).

        Rebuilds the equal-aspect PCA scatters and re-renders every saved
        confusion matrix in the Nature 'Blues' style, then assembles three
        review collages (Training / Test / Additional ligands), each a 2x3
        grid of one PCA panel + five model confusion matrices.
        """
        from functools import partial
        import numpy as np
        from src.config import (
            get_test_name,
            get_test_work_dir,
            primary_geneset_name,
            subset_display,
            confusion_title,
        )
        from src.visualization import (
            draw_pca,
            draw_confusion_matrix,
            assemble_panel_collage,
        )
        from sklearn.decomposition import PCA

        print("=" * 80)
        print("MAT-seq REPLOT (figures only; no model recompute)")
        print("=" * 80)

        models = ["LogisticRegression", "SGDClassifier", "LinearSVC",
                  "RandomForest", "XGBoost"]
        gs = primary_geneset_name()
        test_name = get_test_name()
        fig_root = self.results_dir / "figures"
        collage_dir = fig_root / "collages"
        collage_dir.mkdir(parents=True, exist_ok=True)

        # --- Data from cache (cheap) + frozen selected-gene PCA pipeline -----
        features, labels = self.cache.cached_call(
            prepare_counts, name="prepare_counts",
            force_recompute=False, key_inputs={"ligand_aliases": LIGAND_ALIASES},
        )
        X_main, y_main = extract_subset(features, labels, "main_ligands")
        fs_pca = FeatureSelector.feature_pipeline(
            **FEATURE_SELECTION_CONFIG
        ).set_output(transform="pandas")
        fs_pca.fit(X_main, y_main)

        subset_xy = {}
        for subset in ["main_ligands", "additional_ligands", "bacteria_ligands"]:
            X_sub, y_sub = extract_subset(features, labels, subset)
            subset_xy[subset] = (X_sub, y_sub)
            if subset == "main_ligands":
                X_pca, y_pca = X_sub, y_sub
            else:
                X_pca = pd.concat([X_sub, X_main])
                y_pca = pd.concat([y_sub, y_main])
            self._pca_plot(X_pca, y_pca, subset)
            self._pca_plot(X_pca, y_pca, subset, fs=fs_pca, suffix="_selected")

        def _cm(ax, csv_path, model, subset, show_cbar=False):
            df = pd.read_csv(csv_path, index_col=0)
            draw_confusion_matrix(
                ax, df.values, list(df.index),
                title=confusion_title(model, subset), show_cbar=show_cbar,
                show_tick_labels=False,
            )

        def _pca_panel(ax, X, y, subset_key, title):
            Xt = fs_pca.transform(X)
            coords = PCA(n_components=2).fit_transform(Xt)
            draw_pca(
                ax, coords, np.asarray(y),
                palette=SUBSET_PALETTES.get(subset_key, CUSTOM_PALETTE_9),
                hue_order=SUBSET_CLASS_ORDERS.get(subset_key),
                equal_aspect=True, title=title,
            )

        def build_collage(pca_draw, cm_dir, subset_label_key, out_name, title,
                          file_prefix=""):
            panels = [pca_draw]
            for m in models:
                csv = (Path(cm_dir)
                       / f"{file_prefix}{m}_confusion_matrix_normalized.csv")
                panels.append(
                    partial(_cm, csv_path=csv, model=m, subset=subset_label_key)
                )
            out = assemble_panel_collage(
                panels, collage_dir / out_name, ncols=3, title=title,
            )
            print(f"  Collage saved: {out}")
            return out

        # 1. Training Ligands (internal nested-CV OOF confusion matrices).
        X_tr, y_tr = subset_xy["main_ligands"]
        build_collage(
            partial(_pca_panel, X=X_tr, y=y_tr, subset_key="main_ligands",
                    title="PCA Training Ligands"),
            fig_root / "model_evaluation",
            "main_ligands",
            "training_ligands_collage.png",
            "Training Ligands",
            file_prefix="main_ligands_",
        )

        # 2. Test Ligands (external batch 7086)
        test_fc = get_test_work_dir() / "featurecounts"
        test_cm_dir = self.results_dir / "validation" / test_name / gs / "main_ligands"
        if test_fc.is_dir() and any(test_fc.glob("*.txt")):
            test_counts, test_labels = prepare_counts(featurecounts_dir=test_fc)
            Xv, yv = extract_subset(test_counts, test_labels, "main_ligands")
            test_pca = partial(_pca_panel, X=Xv, y=yv, subset_key="main_ligands",
                               title="PCA Test Ligands")
        else:
            print("  (no external featurecounts; Test PCA panel skipped)")
            test_pca = lambda ax: ax.axis("off")
        build_collage(
            test_pca, test_cm_dir, "external_test",
            "test_ligands_collage.png", "Test Ligands",
        )

        # 3. Additional Ligands (prediction on unseen ligands)
        X_add, y_add = subset_xy["additional_ligands"]
        add_cm_dir = self.results_dir / "predictions" / gs / "additional_ligands"
        build_collage(
            partial(_pca_panel, X=pd.concat([X_add, X_main]),
                    y=pd.concat([y_add, y_main]), subset_key="additional_ligands",
                    title="PCA Additional Ligands"),
            add_cm_dir, "additional_ligands",
            "additional_ligands_collage.png", "Additional Ligands",
        )

        # 4. Bacterial Ligands (heat-killed E. coli / S. aureus)
        X_bact, y_bact = subset_xy["bacteria_ligands"]
        bact_cm_dir = self.results_dir / "predictions" / gs / "bacteria_ligands"
        build_collage(
            partial(_pca_panel, X=pd.concat([X_bact, X_main]),
                    y=pd.concat([y_bact, y_main]), subset_key="bacteria_ligands",
                    title="PCA Bacterial Ligands"),
            bact_cm_dir, "bacteria_ligands",
            "bacterial_ligands_collage.png", "Bacterial Ligands",
        )

        print("\nReplot complete. Collages in:", collage_dir)
        return collage_dir

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

        print("\n--- STEP 1b: FEATURE-NUMBER EVALUATION (MI elbow + KMeans forest scan) ---")
        count_result = self.tune_feature_count()
        mi_elbow = int(count_result["mi_elbow"])
        k_best_prev = FEATURE_SELECTION_CONFIG["k_best"]
        k_best = int(round(mi_elbow))
        FEATURE_SELECTION_CONFIG["k_best"] = k_best
        print(f"  MI elbow={mi_elbow}; using k_best={k_best} (was {k_best_prev}).")

        fs_dir = self.results_dir / "feature_selection"
        fs_dir.mkdir(parents=True, exist_ok=True)
        with open(fs_dir / "mi_elbow.json", "w") as f:
            json.dump(
                {
                    "mi_elbow": mi_elbow,
                    "per_run_elbows": count_result["per_run_elbows"],
                    "k_best": k_best,
                },
                f,
                indent=2,
            )

        forest_result = self.tune_forest_features()
        best = forest_result["best"]
        update_config(
            {
                "k_best": k_best,
                "max_features": best["n_selected"],
                "n_estimators": best["n_estimators"],
                "max_depth": best["max_depth"],
            }
        )
        print(
            f"  Persisted to config.json: k_best={k_best}, "
            f"max_features={best['n_selected']}, n_estimators={best['n_estimators']}, "
            f"max_depth={best['max_depth']}."
        )

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
                fs_pca = FeatureSelector.feature_pipeline(
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
        plot_venn(
            [de_sets[s] for s in subset_names],
            set_labels=tuple(subset_names),
            output_filename="venn_de_across_subsets.png",
            title="DE genes across ligand subsets",
        )

        print("\n--- STEP 3: VENN OF FS vs DE GENES (training subset) ---")
        de_genes = deseq2_main.get_de_genes()
        fs_genes = set(self._select_features(X_train, y_train))
        plot_venn(
            [de_genes, fs_genes],
            set_labels=("Differentially Expressed Genes", "Feature Selection Genes"),
            output_filename="venn_de_vs_fs.png",
        )
        goeaobj, geneid_symbol_mapper = deseq2_main.get_go_objects()
        create_fs_de_go_table(
            de_genes=de_genes,
            fs_genes=fs_genes,
            goeaobj=goeaobj,
            geneid_symbol_mapper=geneid_symbol_mapper,
        )

        print("\n--- STEP 3b: FINAL GENE SETS (selected genes and DE overlap) ---")
        selected_genes = sorted(fs_genes)
        overlap_genes = sorted(fs_genes & set(de_genes))
        union_genes = sorted(fs_genes | set(de_genes))
        primary_gs = primary_geneset_name()
        genesets = {
            primary_gs: selected_genes,
            "de_overlap": overlap_genes,
            "union_stable_de": union_genes,
        }

        fs_dir = self.results_dir / "feature_selection"
        fs_dir.mkdir(parents=True, exist_ok=True)
        pd.Series(selected_genes, name="gene").to_csv(
            fs_dir / "selected_genes.csv", index=False
        )
        pd.Series(overlap_genes, name="gene").to_csv(
            fs_dir / "selected_genes_de_overlap.csv", index=False
        )
        pd.DataFrame(
            {"gene": selected_genes, "in_de": [g in de_genes for g in selected_genes]}
        ).to_csv(fs_dir / "selected_vs_de_overlap_table.csv", index=False)
        print(
            f"  selected: {len(selected_genes)}; DE overlap: {len(overlap_genes)}"
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

        all_models = trainer.nested_cv_summary_["model"].tolist()
        best_two_cv = self._best_two_models(trainer.nested_cv_summary_)
        tuned_params = trainer.selected_params_
        tuned_params_wo = trainer_wo.selected_params_
        print(f"  Best two models by nested CV (at ceiling, tie): {best_two_cv}")

        tables_dir = self.results_dir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)
        trainer.nested_cv_summary_.to_csv(
            tables_dir / "supp_nested_cv_main.csv", index=False
        )
        trainer_wo.nested_cv_summary_.to_csv(
            tables_dir / "supp_nested_cv_no_flapa.csv", index=False
        )

        print(
            "\n--- STEP 4b: DEPLOY ALL MODELS ON FINAL GENE SETS ---"
        )
        geneset_trainers = {
            gs: self._deploy_geneset_trainer(
                X_train, y_train, genes, all_models, tuned_params, trainer.label_encoder
            )
            for gs, genes in genesets.items()
        }
        geneset_trainers_wo = {
            gs: self._deploy_geneset_trainer(
                X_wo, y_wo, genes, all_models, tuned_params_wo, trainer_wo.label_encoder
            )
            for gs, genes in genesets.items()
        }
        for gs, gtr in geneset_trainers.items():
            gtr.save_models(self.results_dir / "models" / gs)
            geneset_trainers_wo[gs].save_models(
                self.results_dir / "models" / "no_flapa" / gs
            )

        print("\n--- STEP 5: MODEL VALIDATION ON EXTERNAL TEST BATCH ---")
        best_two = best_two_cv
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
            table2_rows = []
            for gs, gtr in geneset_trainers.items():
                self._predict(
                    ModelPredictor(gtr), "main_ligands", Xv, yv, val_dir / gs,
                    all_controls=True,
                )
                summ = pd.read_csv(
                    val_dir / gs / "main_ligands" / "test_scores_summary.csv"
                )
                summ.insert(0, "gene_set", gs)
                table2_rows.append(summ)
            external_perf_df = pd.concat(table2_rows, ignore_index=True)
            external_perf_csv = tables_dir / "external_validation_performance.csv"
            external_perf_df.to_csv(external_perf_csv, index=False)
            print(f"  External-validation performance table saved: {external_perf_csv}")
            table2_df = external_perf_df

            primary_gs = primary_geneset_name()
            primary = table2_df[table2_df["gene_set"] == primary_gs]
            best_two = primary.sort_values("f1", ascending=False)["model"].head(2).tolist()
            print(f"  Best two models by EXTERNAL F1 ({primary_gs}): {best_two}")

            self._pca_plot(Xv, yv, "main_ligands", suffix="_external_test")
            self._pca_plot(
                Xv, yv, "main_ligands", fs=fs_pca, suffix="_external_test_selected"
            )

        print(
            "\n--- STEP 6: ENDPOINT PREDICTION ON ADDITIONAL AND BACTERIA LIGANDS ---"
        )
        self._endpoint_predictions(
            geneset_trainers, self.results_dir / "predictions",
            X_other, y_other, X_bact, y_bact,
        )

        print("\n--- STEP 7: TLR HEK BLUE VISUALIZATION ---")
        tlr2_df, tlr4_df, flapa_data = load_tlr_data(
            data_dir=Path(__file__).parent / "data" / "supplementary_data"
        )
        plot_tlr_hek_blue(
            tlr2_df, tlr4_df, flapa_data, output_filename="tlr_hek_blue.png"
        )

        print(
            "\n--- STEP 8: ENDPOINT PREDICTION ON ADDITIONAL AND BACTERIA LIGANDS (without Fla-PA) ---"
        )
        self._endpoint_predictions(
            geneset_trainers_wo, self.results_dir / "predictions" / "no_flapa",
            X_other, y_other, X_bact, y_bact,
        )

        print("\n" + "=" * 80)
        print("PIPELINE COMPLETED SUCCESSFULLY")
        print("=" * 80)
        print(f"Results saved to: {self.results_dir.absolute()}")
        print(f"Cache saved to: {self.cache.cache_dir.absolute()}")

        return {
            "selected_genes": selected_genes,
            "models": trainer,
            "models_wo_flapa": trainer_wo,
            "best_two": best_two,
            "geneset_trainers": geneset_trainers,
            "geneset_trainers_wo": geneset_trainers_wo,
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
    parser.add_argument(
        "--tune-forest",
        action="store_true",
        help="Grid-search the forest selector by KMeans cluster separation (ARI) and exit",
    )
    parser.add_argument(
        "--replot",
        action="store_true",
        help="Regenerate figures/collages from cache + saved CSVs (no model recompute) and exit",
    )

    args = parser.parse_args()
    pipeline = MATseqPipeline(
        cache_dir=args.cache_dir, force_recompute=args.force_recompute
    )

    if args.replot:
        return pipeline.replot()

    if args.tune_features:
        return pipeline.tune_feature_count()

    if args.tune_forest:
        return pipeline.tune_forest_features()

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
