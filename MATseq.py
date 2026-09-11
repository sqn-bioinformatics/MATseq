#!/usr/bin/env python
"""MAT-seq pipeline orchestration script."""

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from src import (
    CUSTOM_PALETTE_9,
    CLASS_ORDER,
    SUBSET_PALETTES,
    DESEQ2_CONFIG,
    FEATURE_SELECTION_CONFIG,
    HYPERPARAMETER_GRIDS,
    MODEL_TRAINING_CONFIG,
    DESeq2,
    mutual_information,
    feature_pipeline,
    forest_kmeans,
    preprocessing_pipeline,
    selected_with_importance,
    ModelTrainer,
    create_fs_de_go_table,
    extract_subset,
    load_tlr_data,
    plot_mutual_information,
    plot_forest_ari_sweep,
    plot_pca,
    plot_tlr_hek_blue,
    plot_venn,
    predict_samples,
    prepare_counts,
)
from src.config import (
    get_config,
    get_genome_dir,
    get_sample_dir,
    get_work_dir,
    primary_geneset_name,
)
from src.compose_figures import compose_figures
from src.make_tables import assemble_supplementary_tables, format_table2

RESULTS_DIR = Path(__file__).parent / "results"

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

    subset_xy: dict[str, tuple[pd.DataFrame, pd.Series]] = {
        "train_ligands": (X_train, y_train),
        "test_ligands": (X_test, y_test),
        "additional_ligands": (X_other, y_other),
        "bacterial_ligands": (X_bact, y_bact),
    }

    print("\n--- STEP 2: DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS ---")
    de_dir = RESULTS_DIR / "differential_gene_expression"
    de_dir.mkdir(parents=True, exist_ok=True)
    go_data_dir = Path(__file__).parent / "data" / "go_terms_support"

    for subset, (X_sub, y_sub) in subset_xy.items():
        deseq2 = DESeq2(
            raw_counts=X_sub,
            sample_labels=y_sub,
            output_dir=de_dir / subset,
            figures_dir=RESULTS_DIR / "figures" / "deseq2" / subset,
            go_terms_dir=RESULTS_DIR / "go_terms" / subset,
            go_fig_dir=RESULTS_DIR / "figures" / "go" / subset,
            go_data_dir=go_data_dir,
            **DESEQ2_CONFIG,
            name=subset,
        )
        deseq2.run_analysis(
            CLASS_ORDER[subset], class_to_compare_to="negative_control"
        )
        if subset == "train_ligands":
            deseq2_train = deseq2
            de_genes = deseq2_train.get_de_genes()
        pd.Series(sorted(deseq2.get_de_genes()), name="gene").to_csv(
            de_dir / f"de_genes_{subset}.csv", index=False
        )

    print("\n--- STEP 3: FEATURE ENGINEERING GENE NUMBER DETERMINATION ---")
    out_dir = RESULTS_DIR / "feature_selection"
    fig_dir = RESULTS_DIR / "figures" / "feature_selection"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    pre_pipe = preprocessing_pipeline().set_output(transform="pandas")
    X_train_pre = pre_pipe.fit_transform(X_train, y_train)

    mi_result = mutual_information(X_train_pre, y_train)
    mi_elbow = mi_result["mi_elbow"]
    print(f"  MI elbow (mean curve): {mi_elbow} genes")

    ari_scan = forest_kmeans(
        X_train_pre,
        y_train,
        k_best=mi_elbow,
        n_estimators=FEATURE_SELECTION_CONFIG["n_estimators"],
        max_depth=FEATURE_SELECTION_CONFIG["max_depth"],
    )
    mi_result["scores"].to_csv(out_dir / "mutual_information.csv", index=False)
    ari_scan.to_csv(out_dir / "forest_kmeans.csv", index=False)

    plot_mutual_information(mi_result, fig_dir)
    plot_forest_ari_sweep(ari_scan, fig_dir)

    print("\n --- STEP 4: PLOT PCA GRAPHS ---")
    pca_dir = RESULTS_DIR / "figures" / "pca"
    pca_dir.mkdir(parents=True, exist_ok=True)

    fs_pipe = feature_pipeline(**FEATURE_SELECTION_CONFIG).set_output(
        transform="pandas"
    )
    fs_pipe.fit(X_train, y_train)

    for subset, (X_sub, y_sub) in subset_xy.items():
        print(f"\nProcessing {subset} subset...")
        if subset in ["train_ligands", "test_ligands"]:
            X_pca, y_pca = X_sub, y_sub
        else:
            # Showing where the other ligands clusters land compared to train
            X_pca = pd.concat([X_sub, X_train])
            y_pca = pd.concat([y_sub, y_train])

        X_pca_pre = pre_pipe.fit_transform(X_pca) # Refitting on each subset

        # The fs_pipe is fit once to train to keep parameters constant
        X_pca_selected = fs_pipe.transform(X_pca) 
        palette = SUBSET_PALETTES.get(subset, CUSTOM_PALETTE_9)
        hue_order = CLASS_ORDER.get(subset)
        for with_names, label_suffix in [(False, ""), (True, "_labeled")]:
            plot_pca(
                X=X_pca_pre,
                labels=y_pca,
                palette=palette,
                hue_order=hue_order,
                with_sample_names=with_names,
                output_path=pca_dir,
                output_filename=f"pca_{subset}{label_suffix}.png",
            )
            plot_pca(
                X=X_pca_selected,
                labels=y_pca,
                palette=palette,
                hue_order=hue_order,
                with_sample_names=with_names,
                output_path=pca_dir,
                output_filename=(
                    f"pca_{subset}_fs{label_suffix}.png"
                ),
            )

    print("\n--- STEP 5: FS vs DE VENN AND GO ---")
    fig_dir = RESULTS_DIR / "figures" / "venn"
    tables_dir = RESULTS_DIR / "fs_de_genesets"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    fs_ranked = selected_with_importance(fs_pipe)
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
        output_dir=RESULTS_DIR / "go_terms",
        fig_dir=RESULTS_DIR / "figures" / "go",
    )

    fs_ranked.to_csv(tables_dir / "fs_genes_ranked.csv", index=False)
    pd.Series(sorted(fs_genes), name="gene").to_csv(
        tables_dir / "fs_gene_names.csv", index=False
    )

    overlap_ranked = fs_ranked.copy()
    overlap_ranked["in_de"] = overlap_ranked["gene"].isin(de_genes)
    overlap_ranked[["gene", "in_de", "rank", "importance"]].to_csv(
        tables_dir / "selected_vs_de_overlap_table.csv", index=False
    )
 
    print("\n--- STEP 6: NESTED CV, GENESET REFIT, VALIDATION AND PREDICTIONS ---")
    primary_gs = primary_geneset_name()  # e.g. 'selected_130'
    endpoint_subsets = {
        "additional_ligands": (X_other, y_other),
        "bacterial_ligands": (X_bact, y_bact),
    }
    mask, mask_test = y_train != "Fla-PA", y_test != "Fla-PA"
    fs_wo = feature_pipeline(**FEATURE_SELECTION_CONFIG).set_output(transform="pandas")
    fs_wo.fit(X_train[mask], y_train[mask])
    de_genes_wo = {
        gene
        for ligand_name, result in deseq2_train.results.items()
        if ligand_name != "Fla-PA"
        for gene in result["significant"].index
    }
    nested_dir = RESULTS_DIR / "nested_cv"
    val_root = RESULTS_DIR / "validation"
    nested_dir.mkdir(parents=True, exist_ok=True)
    panels = {
        "main": {
            "X": X_train, "y": y_train, "X_test": X_test, "y_test": y_test,
            "fs_genes": fs_genes, "de_genes": de_genes,
            "hp_dir": RESULTS_DIR / "hyperparameter_tuning",
            "fig_dir": RESULTS_DIR / "figures" / "model_evaluation",
            "model_dir": RESULTS_DIR / "models",
            "nested_csv": nested_dir / "supp_nested_cv_main.csv",
            "val_dir": val_root / "test_set",
            "perf_csv": val_root / "external_validation_performance.csv",
            "pred_dir": RESULTS_DIR / "predictions",
        },
        "no_flapa": {
            "X": X_train[mask], "y": y_train[mask],
            "X_test": X_test[mask_test], "y_test": y_test[mask_test],
            "fs_genes": set(selected_with_importance(fs_wo)["gene"]),
            "de_genes": de_genes_wo,
            "hp_dir": RESULTS_DIR / "hyperparameter_tuning_no_flapa",
            "fig_dir": RESULTS_DIR / "figures" / "model_evaluation" / "no_flapa",
            "model_dir": RESULTS_DIR / "models" / "no_flapa",
            "nested_csv": nested_dir / "supp_nested_cv_no_flapa.csv",
            "val_dir": val_root / "test_set_no_flapa",
            "perf_csv": val_root / "external_validation_no_flapa_performance.csv",
            "pred_dir": RESULTS_DIR / "predictions" / "no_flapa",
        },
    }

    for panel_name, panel in panels.items():
        print(f"\nPanel: {panel_name}")
        trainer = ModelTrainer(panel["X"], panel["y"], **MODEL_TRAINING_CONFIG)
        trainer.tune_nested(
            HYPERPARAMETER_GRIDS, panel["hp_dir"], panel["fig_dir"], outer_cv=5, inner_cv=3
        ).to_csv(panel["nested_csv"], index=False)

        fs, de = panel["fs_genes"], panel["de_genes"]
        genesets = {
            primary_gs: sorted(fs),
            "de_overlap": sorted(fs & de),
            "union_stable_de": sorted(fs | de),
        }
        validation = {}
        for gene_set, genes in genesets.items():
            trainer.refit(genes)
            trainer.save_models(panel["model_dir"] / gene_set)
            validation[gene_set] = predict_samples(
                trainer, panel["X_test"], panel["y_test"], "test_ligands",
                panel["val_dir"] / gene_set / "test_ligands", all_controls=True,
            )
            if gene_set == "union_stable_de":
                continue
            for subset, (X_end, y_end) in endpoint_subsets.items():
                predict_samples(
                    trainer, X_end, y_end, subset,
                    panel["pred_dir"] / gene_set / subset, all_controls=False,
                )
        pd.concat(validation, names=["gene_set"]).reset_index(level=0).to_csv(
            panel["perf_csv"], index=False
        )

    print("\n--- STEP 7: TLR VISUALIZATION ---")
    tlr2_df, tlr4_df, flapa_data = load_tlr_data(
        data_dir=Path(__file__).parent / "data" / "supplementary_data"
    )
    plot_tlr_hek_blue(
        tlr2_df, tlr4_df, flapa_data,
        output_path=RESULTS_DIR / "figures" / "supplementary",
        output_filename="tlr_hek_blue.png",
    )

    print("\n--- STEP 8: ASSEMBLE COMPOSITE TABLES AND FIGURE COLLAGES ---")
    manuscript_tables_dir = RESULTS_DIR / "tables"
    manuscript_tables_dir.mkdir(parents=True, exist_ok=True)
    composite_figures_dir = Path(__file__).parent / "paper" / "paper_updated" / "figures"
    composite_figures_dir.mkdir(parents=True, exist_ok=True)
    format_table2(
        RESULTS_DIR / "nested_cv" / "supp_nested_cv_main.csv",
        output_dir=manuscript_tables_dir,
    )
    assemble_supplementary_tables(RESULTS_DIR, output_dir=manuscript_tables_dir)
    compose_figures(RESULTS_DIR, composite_figures_dir)

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
