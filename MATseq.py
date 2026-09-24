#!/usr/bin/env python
"""MAT-seq pipeline orchestration script."""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from src import (
    CLASS_ORDER,
    SUBSET_PALETTES,
    SUBSET_DISPLAY_NAMES,
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
    initialize_go,
    load_tlr_data,
    plot_mutual_information,
    plot_forest_ari_sweep,
    plot_pca,
    plot_tlr_hek_blue,
    plot_venn,
    predict_samples,
    prepare_counts,
    assemble_supplementary_tables,
    compose_figures,
    format_table2,
)
from src.config import (
    ADDITIONAL_LIGANDS,
    BACTERIAL_LIGANDS,
    MAIN_LIGANDS,
    get_config,
    get_genome_dir,
    get_sample_dir,
    get_work_dir,
)

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
        if not ok or dry_run:
            return

    print("\n--- STEP 1: COUNT TABLE GENERATION ---")
    counts_dir = RESULTS_DIR / "counts"
    counts_dir.mkdir(parents=True, exist_ok=True)

    count_df = prepare_counts()
    count_df.to_csv(counts_dir / "MATseq_count_summary.csv")

    main = count_df.loc[count_df.index.str.contains("_7128_", regex=False)]
    external = count_df.loc[count_df.index.str.contains("_7086_", regex=False)]
    main_features, main_labels = main.drop(columns="label"), main["label"]
    X_train, y_train = extract_subset(main_features, main_labels, "main_ligands")
    X_other, y_other = extract_subset(main_features, main_labels, "additional_ligands")
    X_bact, y_bact = extract_subset(main_features, main_labels, "bacterial_ligands")
    X_test, y_test = extract_subset(
        external.drop(columns="label"), external["label"], "main_ligands"
    )

    subset_xy: dict[str, tuple[pd.DataFrame, pd.Series]] = {
        "train_ligands": (X_train, y_train),
        "test_ligands": (X_test, y_test),
        "additional_ligands": (X_other, y_other),
        "bacterial_ligands": (X_bact, y_bact),
    }

    print("\n--- STEP 2: DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS ---")
    de_dir = RESULTS_DIR / "differential_gene_expression"
    deseq2_fig_dir = RESULTS_DIR / "figures" / "deseq2"
    go_dir = RESULTS_DIR / "go_terms"
    go_fig_dir = RESULTS_DIR / "figures" / "go"
    goeaobj, geneid_symbol_mapper = initialize_go(
        Path(__file__).parent / "data" / "go_terms_support"
    )

    deseq2_runs: dict[str, DESeq2] = {}
    for subset, (X_sub, y_sub) in subset_xy.items():
        deseq2 = DESeq2(
            raw_counts=X_sub,
            sample_labels=y_sub,
            output_dir=de_dir / subset,
            figures_dir=deseq2_fig_dir / subset,
            go_terms_dir=go_dir / subset,
            go_fig_dir=go_fig_dir / subset,
            goeaobj=goeaobj,
            geneid_symbol_mapper=geneid_symbol_mapper,
            **DESEQ2_CONFIG,
            name=subset,
        )
        deseq2.run_analysis(
            CLASS_ORDER[subset], class_to_compare_to="negative_control"
        )
        deseq2_runs[subset] = deseq2
        pd.Series(sorted(deseq2.get_de_genes()), name="gene").to_csv(
            de_dir / f"de_genes_{subset}.csv", index=False
        )

    deseq2_train = deseq2_runs["train_ligands"]
    de_genes = deseq2_train.get_de_genes()

    print("\n--- STEP 3: FEATURE ENGINEERING GENE NUMBER DETERMINATION ---")
    out_dir = RESULTS_DIR / "feature_selection"
    fs_fig_dir = RESULTS_DIR / "figures" / "feature_selection"
    out_dir.mkdir(parents=True, exist_ok=True)
    fs_fig_dir.mkdir(parents=True, exist_ok=True)

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

    plot_mutual_information(mi_result, fs_fig_dir)
    plot_forest_ari_sweep(ari_scan, fs_fig_dir)

    print("\n --- STEP 4: PLOT PCA GRAPHS ---")
    pca_dir = RESULTS_DIR / "figures" / "pca"

    fs_pipe = feature_pipeline(**FEATURE_SELECTION_CONFIG, k_best=mi_elbow).set_output(
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
        palette = SUBSET_PALETTES[subset]
        hue_order = CLASS_ORDER[subset]
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
    venn_dir = RESULTS_DIR / "figures" / "venn"
    tables_dir = RESULTS_DIR / "fs_de_genesets"
    venn_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    fs_ranked = selected_with_importance(fs_pipe)
    fs_genes = set(fs_ranked["gene"])
    plot_venn(
        [de_genes, fs_genes],
        set_labels=("Differentially Expressed Genes", "Feature Selection Genes"),
        output_path=venn_dir,
        output_filename="venn_de_vs_fs.png",
        title="Differential expression vs. feature selection",
    )
    create_fs_de_go_table(
        de_genes=de_genes,
        fs_genes=fs_genes,
        goeaobj=goeaobj,
        geneid_symbol_mapper=geneid_symbol_mapper,
        output_dir=go_dir,
        fig_dir=go_fig_dir,
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

    print("\n--- STEP 6: NESTED CV, REFIT, VALIDATION AND PREDICTIONS ---")
    endpoint_subsets = {
        "additional_ligands": (X_other, y_other),
        "bacterial_ligands": (X_bact, y_bact),
    }
    mask, mask_test = y_train != "Fla-PA", y_test != "Fla-PA"
    fs_wo = feature_pipeline(**FEATURE_SELECTION_CONFIG, k_best=mi_elbow).set_output(
        transform="pandas"
    )
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
        },
        "no_flapa": {
            "X": X_train[mask], "y": y_train[mask],
            "X_test": X_test[mask_test], "y_test": y_test[mask_test],
            "fs_genes": set(selected_with_importance(fs_wo)["gene"]),
            "de_genes": de_genes_wo,
        },
    }

    for panel_name, panel in panels.items():
        print(f"\nPanel: {panel_name}")
        panel.update(
            hp_dir=RESULTS_DIR / "hyperparameter_tuning" / panel_name,
            fig_dir=RESULTS_DIR / "figures" / "model_evaluation" / panel_name,
            model_dir=RESULTS_DIR / "models" / panel_name,
            nested_csv=nested_dir / f"supp_nested_cv_{panel_name}.csv",
            val_dir=val_root / f"test_set_{panel_name}" / "test_ligands",
            val_fig_dir=RESULTS_DIR / "figures" / "validation" / panel_name / "test_ligands",
            perf_csv=val_root / f"external_validation_{panel_name}_performance.csv",
            pred_dir=RESULTS_DIR / "predictions" / panel_name,
            pred_fig_dir=RESULTS_DIR / "figures" / "predictions" / panel_name,
        )
        cache_dir = panel["hp_dir"] / "pipeline_cache"
        trainer = ModelTrainer(panel["X"], panel["y"], **MODEL_TRAINING_CONFIG)
        shutil.rmtree(cache_dir, ignore_errors=True)
        trainer.tune_nested(
            HYPERPARAMETER_GRIDS, panel["hp_dir"], panel["fig_dir"], cache_dir,
            k_best=mi_elbow, fs_genes=panel["fs_genes"], de_genes=panel["de_genes"],
            outer_cv=5, inner_cv=3,
        ).to_csv(panel["nested_csv"], index=False)
        shutil.rmtree(cache_dir, ignore_errors=True)

        trainer.refit(sorted(panel["fs_genes"]))
        trainer.save_models(panel["model_dir"])
        predict_samples(
            trainer, panel["X_test"], panel["y_test"], "test_ligands",
            panel["val_dir"], panel["val_fig_dir"],
        ).to_csv(panel["perf_csv"], index=False)
        for subset, (X_end, y_end) in endpoint_subsets.items():
            predict_samples(
                trainer, X_end, y_end, subset,
                panel["pred_dir"] / subset, panel["pred_fig_dir"] / subset,
            )

    print("\n--- STEP 7: TLR VISUALIZATION ---")
    supp_data_dir = Path(__file__).parent / "data" / "supplementary_data"
    tlr2_df, tlr4_df, flapa_data = load_tlr_data(data_dir=supp_data_dir)
    plot_tlr_hek_blue(
        tlr2_df, tlr4_df, flapa_data,
        output_path=RESULTS_DIR / "figures" / "supplementary",
        output_filename="tlr_hek_blue.png",
    )

    print("\n--- STEP 8: ASSEMBLE COMPOSITE TABLES AND FIGURE COLLAGES ---")
    manuscript_tables_dir = RESULTS_DIR / "paper"/ "tables"
    composite_figures_dir = RESULTS_DIR/ "paper" / "figures"
    format_table2(
        panels["main"]["nested_csv"],
        manuscript_tables_dir / "table2_formatted.csv",
        manuscript_tables_dir / "Table_2.xlsx",
    )
    format_table2(
        panels["no_flapa"]["nested_csv"],
        manuscript_tables_dir / "Supplementary_Table_8.csv",
    )
    assemble_supplementary_tables(
        de_dir, go_dir, tables_dir, out_dir, supp_data_dir,
        output_dir=manuscript_tables_dir,
    )
    compose_figures(
        [deseq2_fig_dir / "train_ligands" / "LPS_volcano.png",
         deseq2_fig_dir / "train_ligands" / "LPS_histogram.png"],
        composite_figures_dir / "Figure2.png",
        ncols=2,
        margins=(0.45, 0.4, 0.3, 0.3),
    )
    figure3_panels = [
        fs_fig_dir / "mutual_information.png",
        fs_fig_dir / "forest_ari_sweep.png",
        venn_dir / "venn_de_vs_fs.png",
        pca_dir / "pca_train_ligands.png",
        pca_dir / "pca_train_ligands_fs.png",
        go_fig_dir / "de_intersect_fs_go.png",
    ]
    compose_figures(figure3_panels, composite_figures_dir / "Figure3.png")

    main_panel = panels["main"]
    for filename, subset, cm_dir in [
        ("Figure4.png", "train_ligands", main_panel["fig_dir"] / "feature_selection"),
        ("Figure5.png", "test_ligands", main_panel["val_fig_dir"]),
        ("Figure6.png", "additional_ligands",
         main_panel["pred_fig_dir"] / "additional_ligands"),
        ("Figure7.png", "bacterial_ligands",
         main_panel["pred_fig_dir"] / "bacterial_ligands"),
    ]:
        compose_figures(
            [pca_dir / f"pca_{subset}_fs.png"]
            + [cm_dir / f"Confusion_Matrix_{model}.png" for model in HYPERPARAMETER_GRIDS],
            composite_figures_dir / filename,
            ncols=3,
            max_size=(9.84, 6.69),
            title=SUBSET_DISPLAY_NAMES[subset],
        )

    supp1_ligands = [
        (subset, ligand)
        for subset, ligands in [
            ("train_ligands", MAIN_LIGANDS),
            ("additional_ligands", ADDITIONAL_LIGANDS),
            ("bacterial_ligands", BACTERIAL_LIGANDS),
        ]
        for ligand in CLASS_ORDER[subset]
        if ligand in ligands and ligand not in ("negative_control", "LPS")
    ]
    for page, start in enumerate(range(0, len(supp1_ligands), 2)):
        compose_figures(
            [deseq2_fig_dir / subset / f"{ligand}_{kind}.png"
             for subset, ligand in supp1_ligands[start:start + 2]
             for kind in ("volcano", "histogram")],
            composite_figures_dir / f"Supplementary_Figure1p{page + 1}.png",
            first_letter=start * 2,
            margins=(0.45, 0.4, 0.3, 0.3),
        )

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

    run_pipeline(
        snakemake=args.snakemake,
        fastq_dir=args.fastq_dir,
        genome_dir=args.genome_dir,
    )


if __name__ == "__main__":
    main()
