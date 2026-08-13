#!/usr/bin/env python
"""MAT-seq pipeline orchestration script."""

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.model_selection import StratifiedKFold
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
    assemble_supplementary_tables,
    format_table2,
    make_score,
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
 
    print("\n--- STEP 5: GRID SEARCH FOR BEST PARAMETERS (ALL GENES) ---")
    model_dir = RESULTS_DIR / "models"
    model_dir.mkdir(parents=True, exist_ok=True)
    hp_dir = RESULTS_DIR / "hyperparameter_tuning"
    hp_dir.mkdir(parents=True, exist_ok=True)

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
    tuned_params = trainer.selected_params_

    print("\n--- STEP 6: NESTED CV ACROSS GENE-SET CONDITIONS ---")
    tables_dir = RESULTS_DIR / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    primary_gs = primary_geneset_name()
    rng = np.random.default_rng(MODEL_TRAINING_CONFIG.get("random_state", 42))
    random_genes = sorted(
        rng.choice(
            np.asarray(X_train.columns),
            size=min(len(selected_genes), X_train.shape[1]),
            replace=False,
        ).tolist()
    )
    gene_set_conditions = {
        "all_genes": None,
        "feature_selection": selected_genes,
        "fs_plus_de": union_genes,
        "random_selected": random_genes,
    }

    y_enc = trainer.label_encoder.transform(
        y_train.values if isinstance(y_train, pd.Series) else y_train
    )
    outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=trainer.random_state)
    folds = [
        (X_train.iloc[tr], X_train.iloc[te], y_enc[tr], y_enc[te])
        for tr, te in outer.split(np.arange(len(X_train)), y_enc)
    ]
    metric_keys = ["accuracy", "balanced_accuracy", "precision", "recall", "f1", "f1_weighted"]
    table2_rows = []
    for condition, genes in gene_set_conditions.items():
        print(f"\n===== CONDITION: {condition} =====", flush=True)
        for model_name, model in trainer.models.items():
            params = tuned_params.get(model_name, {}).get("chosen_params", {})
            per_fold = {k: [] for k in metric_keys}
            fit_seconds = []
            n_genes = []
            pooled_true, pooled_pred = [], []
            for X_tr, X_te, y_tr, y_te in folds:
                n_genes.append(X_train.shape[1] if genes is None else len(genes))
                steps = list(preprocessing_pipeline().steps)
                if genes is not None:
                    steps.append(("select_genes", ColumnSelector(genes)))
                steps.append(("clf", clone(model)))
                pipe = SkPipeline(steps).set_output(transform="pandas")
                pipe.set_params(**params)
                fit_kwargs = (
                    {"clf__sample_weight": compute_sample_weight("balanced", y_tr)}
                    if model_name == "XGBoost"
                    else {}
                )
                start = time.perf_counter()
                pipe.fit(X_tr, y_tr, **fit_kwargs)
                fit_seconds.append(time.perf_counter() - start)
                y_pred = pipe.predict(X_te)
                scores = make_score(y_te, y_pred)
                for k in metric_keys:
                    per_fold[k].append(scores[k])
                pooled_true.extend(y_te.tolist())
                pooled_pred.extend(y_pred.tolist())
            row = {
                "condition": condition, "model": model_name,
                "n_genes_mean": float(np.mean(n_genes)),
                "fit_seconds_mean": float(np.mean(fit_seconds)),
                "fit_seconds_std": float(np.std(fit_seconds, ddof=1)),
            }
            for k in metric_keys:
                row[f"{k}_mean"] = float(np.mean(per_fold[k]))
                row[f"{k}_std"] = float(np.std(per_fold[k], ddof=1))
            row["pooled_f1"] = make_score(pooled_true, pooled_pred)["f1"]
            table2_rows.append(row)
            print(f"  [{model_name}] f1={row['f1_mean']:.3f}±{row['f1_std']:.3f} "
                  f"acc={row['accuracy_mean']:.3f}±{row['accuracy_std']:.3f}", flush=True)
            pd.DataFrame(table2_rows).to_csv(
                tables_dir / "table2_feature_set_benchmark.csv", index=False
            )

    print("\n--- STEP 7: REFIT DEPLOYMENT MODELS (FEATURE-SELECTED GENES) ---")
    for model_name, model in trainer.models.items():
        steps = [
            *preprocessing_pipeline().steps,
            ("select_genes", ColumnSelector(selected_genes)),
            ("clf", clone(model)),
        ]
        pipe = SkPipeline(steps).set_output(transform="pandas")
        pipe.set_params(**tuned_params.get(model_name, {}).get("chosen_params", {}))
        fit_kwargs = (
            {"clf__sample_weight": compute_sample_weight("balanced", y_enc)}
            if model_name == "XGBoost"
            else {}
        )
        pipe.fit(X_train, y_enc, **fit_kwargs)
        trainer.trained_models[model_name] = pipe
    trainer.save_models(model_dir / primary_gs)

    model_dir_wo = RESULTS_DIR / "models_noflapa"
    model_dir_wo.mkdir(parents=True, exist_ok=True)
    flapa_mask = y_train != "Fla-PA"
    X_wo, y_wo = X_train[flapa_mask], y_train[flapa_mask]
    trainer_wo = ModelTrainer(
        X_wo, y_wo, models=ModelFactory.create_models(**MODEL_FACTORY_CONFIG),
        **MODEL_TRAINING_CONFIG,
    )
    y_enc_wo = trainer_wo.label_encoder.fit_transform(
        y_wo.values if isinstance(y_wo, pd.Series) else y_wo
    )
    for model_name, model in trainer_wo.models.items():
        steps = [
            *preprocessing_pipeline().steps,
            ("select_genes", ColumnSelector(selected_genes)),
            ("clf", clone(model)),
        ]
        pipe_wo = SkPipeline(steps).set_output(transform="pandas")
        pipe_wo.set_params(**tuned_params.get(model_name, {}).get("chosen_params", {}))
        fit_kwargs_wo = (
            {"clf__sample_weight": compute_sample_weight("balanced", y_enc_wo)}
            if model_name == "XGBoost"
            else {}
        )
        pipe_wo.fit(X_wo, y_enc_wo, **fit_kwargs_wo)
        trainer_wo.trained_models[model_name] = pipe_wo
    trainer_wo.save_models(model_dir_wo / primary_gs)

    print("\n--- STEP 8: VALIDATION ON TEST SET ---")
    val_dir = RESULTS_DIR / "validation" / "test_set"
    val_dir_wo = RESULTS_DIR / "validation" / "test_set_no_flapa"
    validation_tables_dir = RESULTS_DIR / "validation"
    validation_tables_dir.mkdir(parents=True, exist_ok=True)

    out_dir = val_dir / primary_gs / "test_ligands"
    out_dir.mkdir(parents=True, exist_ok=True)
    predictor = ModelPredictor(trainer)
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
    summary.insert(0, "gene_set", primary_gs)
    external_perf_csv = validation_tables_dir / "external_validation_performance.csv"
    summary.to_csv(external_perf_csv, index=False)

    flapa_mask_ext = y_test != "Fla-PA"
    X_test_wo, y_test_wo = X_test[flapa_mask_ext], y_test[flapa_mask_ext]
    out_dir_wo = val_dir_wo / primary_gs / "test_ligands"
    out_dir_wo.mkdir(parents=True, exist_ok=True)
    predictor_wo = ModelPredictor(trainer_wo)
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
    summary_wo.insert(0, "gene_set", primary_gs)
    external_perf_wo_csv = (
        validation_tables_dir / "external_validation_no_flapa_performance.csv"
    )
    summary_wo.to_csv(external_perf_wo_csv, index=False)

    print("\n--- STEP 9: PREDICTIONS ON OTHER LIGANDS---")
    predictions_dir = RESULTS_DIR / "predictions"
    predictions_dir.mkdir(parents=True, exist_ok=True)
    predictions_dir_wo = RESULTS_DIR / "predictions" / "no_flapa"
    predictions_dir_wo.mkdir(parents=True, exist_ok=True)
    endpoint_subsets = {
        "additional_ligands": (X_other, y_other),
        "bacterial_ligands": (X_bact, y_bact),
    }

    for gene_trainer, base_dir in (
        (trainer, predictions_dir),
        (trainer_wo, predictions_dir_wo),
    ):
        for subset_key, (X_end, y_end) in endpoint_subsets.items():
            out_dir = base_dir / primary_gs / subset_key
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

    print("\n--- STEP 10: TLR VISUALIZATION ---")
    tlr2_df, tlr4_df, flapa_data = load_tlr_data(
        data_dir=Path(__file__).parent / "data" / "supplementary_data"
    )
    plot_tlr_hek_blue(tlr2_df, tlr4_df, flapa_data, output_filename="tlr_hek_blue.png")

    print("\n--- STEP 11: NO-FLA-PA NESTED CV ACROSS GENE-SET CONDITIONS ---")
    nested_cv_dir = RESULTS_DIR / "nested_cv"
    nested_cv_dir.mkdir(parents=True, exist_ok=True)

    outer_wo = StratifiedKFold(n_splits=5, shuffle=True, random_state=trainer_wo.random_state)
    folds_wo = [
        (X_wo.iloc[tr], X_wo.iloc[te], y_enc_wo[tr], y_enc_wo[te])
        for tr, te in outer_wo.split(np.arange(len(X_wo)), y_enc_wo)
    ]
    no_flapa_rows = []
    for condition, genes in gene_set_conditions.items():
        print(f"\n===== CONDITION: {condition} =====", flush=True)
        for model_name, model in trainer_wo.models.items():
            params = tuned_params.get(model_name, {}).get("chosen_params", {})
            per_fold = {k: [] for k in metric_keys}
            fit_seconds = []
            n_genes = []
            pooled_true, pooled_pred = [], []
            for X_tr, X_te, y_tr, y_te in folds_wo:
                n_genes.append(X_wo.shape[1] if genes is None else len(genes))
                steps = list(preprocessing_pipeline().steps)
                if genes is not None:
                    steps.append(("select_genes", ColumnSelector(genes)))
                steps.append(("clf", clone(model)))
                pipe = SkPipeline(steps).set_output(transform="pandas")
                pipe.set_params(**params)
                fit_kwargs = (
                    {"clf__sample_weight": compute_sample_weight("balanced", y_tr)}
                    if model_name == "XGBoost"
                    else {}
                )
                start = time.perf_counter()
                pipe.fit(X_tr, y_tr, **fit_kwargs)
                fit_seconds.append(time.perf_counter() - start)
                y_pred = pipe.predict(X_te)
                scores = make_score(y_te, y_pred)
                for k in metric_keys:
                    per_fold[k].append(scores[k])
                pooled_true.extend(y_te.tolist())
                pooled_pred.extend(y_pred.tolist())
            row = {
                "condition": condition, "model": model_name,
                "n_genes_mean": float(np.mean(n_genes)),
                "fit_seconds_mean": float(np.mean(fit_seconds)),
                "fit_seconds_std": float(np.std(fit_seconds, ddof=1)),
            }
            for k in metric_keys:
                row[f"{k}_mean"] = float(np.mean(per_fold[k]))
                row[f"{k}_std"] = float(np.std(per_fold[k], ddof=1))
            row["pooled_f1"] = make_score(pooled_true, pooled_pred)["f1"]
            no_flapa_rows.append(row)
            print(f"  [{model_name}] f1={row['f1_mean']:.3f}±{row['f1_std']:.3f} "
                  f"acc={row['accuracy_mean']:.3f}±{row['accuracy_std']:.3f}", flush=True)
            pd.DataFrame(no_flapa_rows).to_csv(
                nested_cv_dir / "supp_nested_cv_no_flapa.csv", index=False
            )

    print("\n--- STEP 12: ASSEMBLE COMPOSITE TABLES AND FIGURE COLLAGES ---")
    format_table2(tables_dir / "table2_feature_set_benchmark.csv", output_dir=tables_dir)
    assemble_supplementary_tables(RESULTS_DIR, tables_dir)

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
