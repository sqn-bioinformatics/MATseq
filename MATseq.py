#!/usr/bin/env python
"""MAT-seq pipeline orchestration script."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
from sklearn.base import clone
from sklearn.pipeline import Pipeline as SkPipeline
from sklearn.utils.class_weight import compute_sample_weight
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))

from src import (
    CUSTOM_PALETTE_9,
    CLASS_ORDER,
    SUBSET_PALETTES,
    ColumnSelector,
    DESEQ2_CONFIG,
    FEATURE_SELECTION_CONFIG,
    HYPERPARAMETER_GRIDS,
    MODEL_FACTORY_CONFIG,
    MODEL_TRAINING_CONFIG,
    DESeq2,
    mutual_information,
    feature_pipeline,
    forest_kmeans,
    preprocessing_pipeline,
    selected_with_importance,
    ModelFactory,
    ModelPredictor,
    ModelTrainer,
    create_fs_de_go_table,
    extract_subset,
    load_tlr_data,
    plot_mutual_information,
    plot_forest_ari_sweep,
    plot_pca,
    plot_probability_heatmap,
    plot_tlr_hek_blue,
    plot_venn,
    prepare_counts,
    subset_display,
)
from src.config import (
    FOREST_SELECTION_GRID,
    get_config,
    get_genome_dir,
    get_sample_dir,
    get_work_dir,
    primary_geneset_name,
    update_config,
)
from src.compose_figures import compose_figures
from src.make_tables import assemble_supplementary_tables, format_table2

RESULTS_DIR = Path.cwd() / "results"

def run_snakemake_preprocessing(
    fastq_dir: Path | None = None,
    genome_dir: Path | None = None,
    work_dir: Path | None = None,
    dry_run: bool = False,
) -> bool:
    """Run Snakemake preprocessing to generate featureCounts outputs."""
    fastq_dir = fastq_dir or get_sample_dir()
    genome_dir = genome_dir or get_genome_dir()
    work_dir = work_dir or get_work_dir()
    snakefile = Path.cwd() / "pipeline" / "0_MATseq.smk"

    if not Path(genome_dir).exists():
        print(f"Error: genome_dir is required and must exist. Current: '{genome_dir}'")
        return False
    if not Path(fastq_dir).exists():
        print(f"Error: sample_dir is required and must exist. Current: '{fastq_dir}'")
        return False
    if not snakefile.exists():
        print(f"Error: Snakemake pipeline not found at {snakefile}")
        return False

    cmd = [
        "poetry",
        "run",
        "snakemake",
        "--use-conda",
        "--cores",
        str(get_config("snakemake.threads")),
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

    print(f"FASTQ: {fastq_dir}")
    print(f"Genome: {genome_dir}")
    result = subprocess.run(cmd, cwd=str(Path.cwd()))
    if result.returncode == 0:
        print("Snakemake preprocessing completed successfully")
        return True
    print(f"Snakemake preprocessing failed with code {result.returncode}")
    return False

def run_pipeline(
    snakemake: str | None = None,
    fastq_dir: Path | None = None,
    genome_dir: Path | None = None,
) -> None:
    RESULTS_DIR.mkdir(exist_ok=True, parents=True)

    print("=" * 80)
    print("MAT-seq Analysis Pipeline")
    print("=" * 80)

    print("\n--- STEP 0: SNAKEMAKE ---")
    if snakemake is None:
        print("Skipping Snakemake preprocessing; no --snakemake flag provided.")
    else:
        dry_run = snakemake == "dry-run"
        ok = run_snakemake_preprocessing(
            fastq_dir=fastq_dir,
            genome_dir=genome_dir,
            dry_run=dry_run,
        )
        if not ok:
            return None
        if dry_run:
            return None

    print("\n--- STEP 1: COUNT TABLE GENERATION ---")
    counts_dir = RESULTS_DIR / "counts"
    counts_dir.mkdir(parents=True, exist_ok=True)

    count_df = prepare_counts()
    count_df.to_csv(counts_dir / "MATseq_count_summary.csv")
    batch_specs = {"main_dataset": "7128", "external_test": "7086"}

    for split, batch in batch_specs.items():
        df = count_df.loc[count_df.index.str.contains(f"_{batch}_", regex=False)]
        features = df.drop(columns="label")
        labels = df["label"]
        if split == "external_test":
            X_test, y_test = extract_subset(features, labels, "main_ligands")
        else:
            X_train, y_train = extract_subset(features, labels, "main_ligands")
            X_other, y_other = extract_subset(features, labels, "additional_ligands")
            X_bact, y_bact = extract_subset(features, labels, "bacterial_ligands")

    print("\n--- STEP 2: DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS ---")
    de_dir = RESULTS_DIR / "differential_gene_expression"
    de_dir.mkdir(parents=True, exist_ok=True)

    subset_xy: dict[str, tuple[pd.DataFrame, pd.Series]] = {
        "train_ligands": (X_train, y_train),
        "test_ligands": (X_test, y_test),
        "additional_ligands": (X_other, y_other),
        "bacterial_ligands": (X_bact, y_bact),
    }
    pre = preprocessing_pipeline().set_output(transform="pandas")
    fs_pipe = feature_pipeline(**FEATURE_SELECTION_CONFIG).set_output(
        transform="pandas"
    )
    fs_pipe.fit(X_train, y_train)

    deseq2_tables: dict[str, DESeq2] = {}
    for subset, (X_sub, y_sub) in subset_xy.items():
        print(f"\nProcessing {subset} subset...")
        if subset in ["train_ligands", "test_ligands"]:
            X_pca, y_pca = X_sub, y_sub
        else:
            X_pca = pd.concat([X_sub, X_train])
            y_pca = pd.concat([y_sub, y_train])

        X_pca_pre = pre.fit_transform(X_pca)
        X_pca_selected = fs_pipe.transform(X_pca)
        palette = SUBSET_PALETTES.get(subset, CUSTOM_PALETTE_9)
        hue_order = CLASS_ORDER.get(subset)
        for with_names, label_suffix in [(False, ""), (True, "_labeled")]:
            plot_pca(
                X=X_pca_pre,
                labels=y_pca,
                palette=palette,
                hue_order=hue_order,
                name=f"{subset}{label_suffix}",
                with_sample_names=with_names,
                output_filename=f"{subset}_pca{label_suffix}.png",
                equal_aspect=True,
            )
            plot_pca(
                X=X_pca_selected,
                labels=y_pca,
                palette=palette,
                hue_order=hue_order,
                name=f"{subset}_feature_selected{label_suffix}",
                with_sample_names=with_names,
                output_filename=(
                    f"{subset}_feature_selected{label_suffix}.png"
                ),
                equal_aspect=True,
            )

        deseq2 = DESeq2(
            raw_counts=X_sub,
            sample_labels=y_sub,
            padj_threshold=DESEQ2_CONFIG.get("padj_threshold", 0.05),
            log2fc_threshold=DESEQ2_CONFIG.get("log2fc_threshold", 2.0),
            n_cpus=DESEQ2_CONFIG.get("n_cpus", 42),
            name=subset,
        )
        deseq2.run_analysis(
            CLASS_ORDER[subset], negative_control="negative_control"
        )
        deseq2_tables[subset] = deseq2

        pd.Series(sorted(deseq2.get_de_genes()), name="gene").to_csv(
            de_dir / f"de_genes_{subset}.csv", index=False
        )

    deseq2_train = deseq2_tables["train_ligands"]

    print("\n--- STEP 3: FINAL GENE NUMBER DETERMINATION ---")
    out_dir = RESULTS_DIR / "feature_selection"
    fig_dir = RESULTS_DIR / "figures" / "feature_selection"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    mi_result = mutual_information(fs_pipe, X_train, y_train)
    mi_elbow = mi_result["mi_elbow"]
    print(f"  MI elbow (mean curve): {mi_elbow} genes")

    kmeans_result = forest_kmeans(
        fs_pipe,
        X_train,
        y_train,
        grid=FOREST_SELECTION_GRID,
        k_best=mi_elbow,
    )

    mi_result["scores"].to_csv(out_dir / "mutual_information.csv", index=False)
    kmeans_result["scan"].to_csv(out_dir / "forest_kmeans.csv", index=False)
    
    best = kmeans_result["best"]
    n_genes = kmeans_result["plateau"]

    plot_mutual_information(mi_result, fig_dir)
    plot_forest_ari_sweep(kmeans_result["scan"], fig_dir, selected=n_genes)

    print(
        f"  Best by ARI: n_estimators={best['n_estimators']}, "
        f"max_depth={best['max_depth']} (ARI={best['ari']:.3f})"
    )
    print(f"  Plateau gene count (max_features): {n_genes}")

    update_config(
        {
            "k_best": mi_elbow,
            "n_estimators": best["n_estimators"],
            "max_depth": best["max_depth"],
            "max_features": n_genes,
        }
    )

    print("\n--- STEP 4: FS vs DE VENN AND GO ---")
    fig_dir = RESULTS_DIR / "figures" / "venn"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = RESULTS_DIR / "fs_de_genesets"
    tables_dir.mkdir(parents=True, exist_ok=True)
    
    de_genes = deseq2_train.get_de_genes()
    fs = feature_pipeline(**FEATURE_SELECTION_CONFIG).set_output(
        transform="pandas"
    )
    fs.fit(X_train, y_train)
    fs_ranked = selected_with_importance(fs)
    fs_genes = set(fs_ranked["gene"])
    plot_venn(
        [de_genes, fs_genes],
        set_labels=("Differentially Expressed Genes", "Feature Selection Genes"),
        output_path=fig_dir,
        output_filename="venn_de_vs_fs.png",
        title="Differential expression vs. feature selection",
    )
    goeaobj, geneid_symbol_mapper = deseq2_train.get_go_objects()
    create_fs_de_go_table(
        de_genes=de_genes,
        fs_genes=fs_genes,
        goeaobj=goeaobj,
        geneid_symbol_mapper=geneid_symbol_mapper,
    )

    selected_genes = sorted(fs_genes)
    overlap_genes = sorted(fs_genes & set(de_genes))
    union_genes = sorted(fs_genes | set(de_genes))
    
    fs_ranked.to_csv(tables_dir / "fs_genes_ranked.csv", index=False)
    pd.Series(selected_genes, name="gene").to_csv(
        tables_dir / "fs_gene_names.csv", index=False
    )

    overlap_ranked = fs_ranked.copy()
    overlap_ranked["in_de"] = overlap_ranked["gene"].isin(de_genes)
    overlap_ranked[["gene", "in_de", "rank", "importance"]].to_csv(
        tables_dir / "selected_vs_de_overlap_table.csv", index=False
    )
 
    print("\n--- STEP 5: BEST PARMETER DETERMINATION WITH NESTED CV ---")
    model_dir = RESULTS_DIR / "models"
    model_dir.mkdir(parents=True, exist_ok=True)
    model_dir_wo = RESULTS_DIR / "models" / "no_flapa"
    model_dir_wo.mkdir(parents=True, exist_ok=True)
    tables_dir = RESULTS_DIR / "nested_cv"
    tables_dir.mkdir(parents=True, exist_ok=True)
    hp_dir = RESULTS_DIR / "hyperparameter_tuning"
    hp_dir.mkdir(parents=True, exist_ok=True)
    hp_dir_wo = RESULTS_DIR / "hyperparameter_tuning_no_flapa"
    hp_dir_wo.mkdir(parents=True, exist_ok=True)

    # Include the Fla-Pa class
    models = ModelFactory.create_models(**MODEL_FACTORY_CONFIG)
    trainer = ModelTrainer(X_train, y_train, models=models, **MODEL_TRAINING_CONFIG)
    trainer.tune_nested(
        X_train,
        y_train,
        param_grids=HYPERPARAMETER_GRIDS,
        output_dir=hp_dir,
        outer_cv=5,
        inner_cv=3,
    )
    trainer.save_models(model_dir)

    # Exclude the Fla-Pa class to check if training improves
    flapa_mask = y_train != "Fla-PA"
    X_wo, y_wo = X_train[flapa_mask], y_train[flapa_mask]
    models_wo = ModelFactory.create_models(**MODEL_FACTORY_CONFIG)
    trainer_wo = ModelTrainer(X_wo, y_wo, models=models_wo, **MODEL_TRAINING_CONFIG)
    trainer_wo.tune_nested(
        X_wo,
        y_wo,
        param_grids=HYPERPARAMETER_GRIDS,
        output_dir=hp_dir_wo,
        outer_cv=5,
        inner_cv=3,
    )
    trainer_wo.save_models(model_dir_wo)

    all_models = trainer.nested_cv_summary_["model"].tolist()
    tuned_params = trainer.selected_params_
    tuned_params_wo = trainer_wo.selected_params_

    trainer.nested_cv_summary_.to_csv(tables_dir / "supp_nested_cv_main.csv", index=False)
    trainer_wo.nested_cv_summary_.to_csv(
        tables_dir / "supp_nested_cv_no_flapa.csv", index=False
    )
    
    print("\n--- STEP 6: REFIT THE MODELS WITH BEST PARAMETERS ON GENESETS  ---")
    geneset_model_dir = RESULTS_DIR / "models"
    geneset_model_dir.mkdir(parents=True, exist_ok=True)
    geneset_model_dir_wo = RESULTS_DIR / "models" / "no_flapa"
    geneset_model_dir_wo.mkdir(parents=True, exist_ok=True)

    primary_gs = primary_geneset_name()
    genesets = {
        primary_gs: selected_genes,
        "de_overlap": overlap_genes,
        "union_stable_de": union_genes,
    }
    fs_wo = feature_pipeline(**FEATURE_SELECTION_CONFIG).set_output(
        transform="pandas"
    )
    fs_wo.fit(X_wo, y_wo)
    fs_wo_ranked = selected_with_importance(fs_wo)
    fs_genes_wo = set(fs_wo_ranked["gene"])
    de_genes_wo = {
        gene
        for ligand_name, result in deseq2_train.results.items()
        if ligand_name != "Fla-PA"
        for gene in result["significant"].index
    }
    genesets_wo = {
        primary_gs: sorted(fs_genes_wo),
        "de_overlap": sorted(fs_genes_wo & de_genes_wo),
        "union_stable_de": sorted(fs_genes_wo | de_genes_wo),
    }

    y_enc = trainer.label_encoder.transform(
        y_train.values if isinstance(y_train, pd.Series) else y_train
    )
    y_enc_wo = trainer_wo.label_encoder.transform(
        y_wo.values if isinstance(y_wo, pd.Series) else y_wo
    )

    geneset_trainers = {}
    geneset_trainers_wo = {}
    for gene_set, genes in genesets.items():
        base = ModelFactory.create_models(**MODEL_FACTORY_CONFIG)
        trained = {}
        for name in all_models:
            pipe = SkPipeline(
                [
                    *preprocessing_pipeline().steps,
                    ("select_genes", ColumnSelector(genes)),
                    ("clf", clone(base[name])),
                ]
            ).set_output(transform="pandas")
            params = tuned_params.get(name, {}).get("chosen_params", {})
            pipe.set_params(**params)
            fit_kwargs = (
                {"clf__sample_weight": compute_sample_weight("balanced", y_enc)}
                if name == "XGBoost"
                else {}
            )
            pipe.fit(X_train, y_enc, **fit_kwargs)
            trained[name] = pipe

        geneset_trainers[gene_set] = ModelTrainer.from_trained_models(
            trained, trainer.label_encoder
        )

        base_wo = ModelFactory.create_models(**MODEL_FACTORY_CONFIG)
        trained_wo = {}
        genes_wo = genesets_wo[gene_set]
        for name in all_models:
            pipe_wo = SkPipeline(
                [
                    *preprocessing_pipeline().steps,
                    ("select_genes", ColumnSelector(genes_wo)),
                    ("clf", clone(base_wo[name])),
                ]
            ).set_output(transform="pandas")
            params_wo = tuned_params_wo.get(name, {}).get("chosen_params", {})
            pipe_wo.set_params(**params_wo)
            fit_kwargs_wo = (
                {"clf__sample_weight": compute_sample_weight("balanced", y_enc_wo)}
                if name == "XGBoost"
                else {}
            )
            pipe_wo.fit(X_wo, y_enc_wo, **fit_kwargs_wo)
            trained_wo[name] = pipe_wo

        geneset_trainers_wo[gene_set] = ModelTrainer.from_trained_models(
            trained_wo, trainer_wo.label_encoder
        )

    for gene_set, gene_trainer in geneset_trainers.items():
        gene_trainer.save_models(geneset_model_dir / gene_set)
        geneset_trainers_wo[gene_set].save_models(geneset_model_dir_wo / gene_set)

    print("\n--- STEP 7: VALIDATION ON TEST SET ---")
    val_dir = RESULTS_DIR / "validation" / "test_set"
    val_dir_wo = RESULTS_DIR / "validation" / "test_set_no_flapa"
    validation_tables_dir = RESULTS_DIR / "validation"
    validation_tables_dir.mkdir(parents=True, exist_ok=True)

    validation_rows = []
    for gene_set, gene_trainer in geneset_trainers.items():
        out_dir = val_dir / gene_set / "test_ligands"
        out_dir.mkdir(parents=True, exist_ok=True)
        predictor = ModelPredictor(gene_trainer)
        predictor.predict_samples(
            X_test, sample_names=X_test.index.to_numpy(), y_test=y_test
        )
        for model_name, pred_df in predictor.predictions.items():
            pred_df.to_csv(out_dir / f"{model_name}_predictions.csv", index=False)
        for model_name, proba_df in predictor.probabilities.items():
            proba_df.to_csv(out_dir / f"{model_name}_probabilities.csv")
        summary = predictor.evaluate(out_dir, subset="test_ligands")
        class_order = CLASS_ORDER.get("test_ligands", CLASS_ORDER["test_ligands"])

        for model_name, proba_df in predictor.probabilities.items():
            plot_probability_heatmap(
                proba_df, class_order,
                title=f"{model_name} Prediction Probabilities {subset_display('test_ligands')}",
                true_labels=predictor.y_test, all_controls=True,
                output_dir=out_dir,
                filename=f"{model_name}_probabilities_heatmap.png",
            )
        summary.insert(0, "gene_set", gene_set)
        validation_rows.append(summary)

    external_perf_df = pd.concat(validation_rows, ignore_index=True)
    external_perf_csv = validation_tables_dir / "external_validation_performance.csv"
    external_perf_df.to_csv(external_perf_csv, index=False)

    flapa_mask_ext = y_test != "Fla-PA"
    X_test_wo, y_test_wo = X_test[flapa_mask_ext], y_test[flapa_mask_ext]
    validation_wo_rows = []
    
    for gene_set, gene_trainer in geneset_trainers_wo.items():
        out_dir_wo = val_dir_wo / gene_set / "test_ligands"
        out_dir_wo.mkdir(parents=True, exist_ok=True)
        predictor_wo = ModelPredictor(gene_trainer)
        predictor_wo.predict_samples(
            X_test_wo, sample_names=X_test_wo.index.to_numpy(), y_test=y_test_wo
        )
        for model_name, pred_df in predictor_wo.predictions.items():
            pred_df.to_csv(out_dir_wo / f"{model_name}_predictions.csv", index=False)
        for model_name, proba_df in predictor_wo.probabilities.items():
            proba_df.to_csv(out_dir_wo / f"{model_name}_probabilities.csv")
        summary_wo = predictor_wo.evaluate(out_dir_wo, subset="test_ligands")
        class_order_wo = CLASS_ORDER.get("test_ligands", CLASS_ORDER["test_ligands"])
        for model_name, proba_df in predictor_wo.probabilities.items():
            plot_probability_heatmap(
                proba_df, class_order_wo,
                title=f"{model_name} Prediction Probabilities {subset_display('test_ligands')}",
                true_labels=predictor_wo.y_test, all_controls=True,
                output_dir=out_dir_wo,
                filename=f"{model_name}_probabilities_heatmap.png",
            )
        summary_wo.insert(0, "gene_set", gene_set)
        validation_wo_rows.append(summary_wo)
    
    external_perf_wo_df = pd.concat(validation_wo_rows, ignore_index=True)
    external_perf_wo_csv = (
        validation_tables_dir / "external_validation_no_flapa_performance.csv"
    )
    external_perf_wo_df.to_csv(external_perf_wo_csv, index=False)


    print("\n--- STEP 8: PREDICTIONS ON OTHER LIGANDS---")
    predictions_dir = RESULTS_DIR / "predictions"
    predictions_dir.mkdir(parents=True, exist_ok=True)
    predictions_dir_wo = RESULTS_DIR / "predictions" / "no_flapa"
    predictions_dir_wo.mkdir(parents=True, exist_ok=True)
    endpoint_subsets = {
        "additional_ligands": (X_other, y_other),
        "bacterial_ligands": (X_bact, y_bact),
    }

    for trainers, base_dir in (
        (geneset_trainers, predictions_dir),
        (geneset_trainers_wo, predictions_dir_wo),
    ):
        for gene_set in (primary_gs, "de_overlap"):
            gene_trainer = trainers[gene_set]
            for subset_key, (X_end, y_end) in endpoint_subsets.items():
                out_dir = base_dir / gene_set / subset_key
                out_dir.mkdir(parents=True, exist_ok=True)
                predictor = ModelPredictor(gene_trainer)
                predictor.predict_samples(
                    X_end, sample_names=X_end.index.to_numpy(), y_test=y_end
                )
                for model_name, pred_df in predictor.predictions.items():
                    pred_df.to_csv(
                        out_dir / f"{model_name}_predictions.csv", index=False
                    )
                for model_name, proba_df in predictor.probabilities.items():
                    proba_df.to_csv(out_dir / f"{model_name}_probabilities.csv")
                predictor.evaluate(out_dir, subset=subset_key)
                class_order_end = CLASS_ORDER.get(subset_key, CLASS_ORDER["train_ligands"])
                for model_name, proba_df in predictor.probabilities.items():
                    plot_probability_heatmap(
                        proba_df, class_order_end,
                        title=f"{model_name} Prediction Probabilities {subset_display(subset_key)}",
                        true_labels=predictor.y_test, all_controls=False,
                        output_dir=out_dir,
                        filename=f"{model_name}_probabilities_heatmap.png",
                    )

    print("\n--- STEP 9: TLR VISUALIZATION ---")
    tlr2_df, tlr4_df, flapa_data = load_tlr_data(
        data_dir=Path(__file__).parent / "data" / "supplementary_data"
    )
    plot_tlr_hek_blue(tlr2_df, tlr4_df, flapa_data, output_filename="tlr_hek_blue.png")

    print("\n--- STEP 10: ASSEMBLE COMPOSITE TABLES AND FIGURE COLLAGES ---")
    manuscript_tables_dir = RESULTS_DIR / "tables"
    manuscript_tables_dir.mkdir(parents=True, exist_ok=True)
    composite_figures_dir = RESULTS_DIR / "tables"
    composite_figures_dir.mkdir(parents=True, exist_ok=True)
    format_table2(
        RESULTS_DIR / "nested_cv" / "supp_nested_cv_main.csv",
        output_dir=manuscript_tables_dir,
    )
    assemble_supplementary_tables(RESULTS_DIR, output_dir=manuscript_tables_dir)
    compose_figures()

    print("\n" + "=" * 80)
    print("PIPELINE COMPLETED SUCCESSFULLY")
    print("=" * 80)
    print(f"Results saved to: {RESULTS_DIR.absolute()}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="MAT-seq analysis pipeline with full publication analysis"
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

    return run_pipeline(
        snakemake=args.snakemake,
        fastq_dir=args.fastq_dir,
        genome_dir=args.genome_dir,
    )


if __name__ == "__main__":
    main()
