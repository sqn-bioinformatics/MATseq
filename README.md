<h1 align="center">
  <img src="assets/logo_MATseq.png" width="180"><br>
  MATseq
</h1>

<p align="center" style="margin-top: 6px;">
  <i>Machine-learning classification of pyrogen-induced monocyte transcriptomic signatures</i>
</p>

## Abstract

The monocyte activation test is an in vitro pyrogenicity assessment method that can utilise human peripheral blood mononuclear cells to detect pyrogens in injectable drugs, providing a binary outcome that indicates the presence or absence of a pyrogen. The added ability to distinguish between different types of pyrogens would greatly expand the applicability of the test, for example, by allowing to pinpoint the source of a contaminating pyrogen in pharmaceutical products. Pyrogens activate a unique set of pattern recognition receptors (PRRs), which contribute to inflammation, yielding distinct transcriptomic activation signatures. In this paper, we capture the unique expression signatures of activated monocytes through bulk RNA sequencing and introduce a data preprocessing pipeline that allows the training of a machine-learning model to classify pyrogenic contaminants. Using a dataset of 108 samples stimulated with five classes of PRR agonists, we could differentiate between these classes with more than 97% F1 on test data. We further demonstrate the model's capacity to generalise on the previously unseen data using different ligands for the same PRRs as well as heat-killed *Escherichia coli* and *Staphylococcus aureus*.

MATseq runs the full analysis end to end: Snakemake read preprocessing, DESeq2 differential expression with GO enrichment, data-driven selection of the classifier gene panel, nested cross-validated model tuning, external-batch validation, and assembly of the manuscript tables and composite figures.

## Installation

Requires Python >=3.10,<3.13 and [Poetry](https://python-poetry.org/). Snakemake rules pull their own tools through Conda (`pipeline/environment.yml`).

```bash
git clone https://github.com/sqn-bioinformatics/MATseq.git MATseq
cd MATseq
poetry install
```

## Usage

```bash
poetry run python MATseq.py --snakemake dry-run   # validate the Snakemake DAG, then exit
poetry run python MATseq.py --snakemake run       # FASTQ -> featureCounts, then exit
poetry run python MATseq.py                       # analysis on existing featureCounts
```

| Option | Description |
| --- | --- |
| `--snakemake {dry-run,run}` | Run read preprocessing; omit to start from existing featureCounts output |
| `--fastq-dir PATH` | Override `snakemake.sample_dir` |
| `--genome-dir PATH` | Override `snakemake.genome_dir` |

All output locations are defined in `MATseq.py` and passed into `src/`; only the inputs are configurable.

## Input data

| Input | Source |
| --- | --- |
| Raw reads (`fastq.gz`) | [GEO GSE313994](https://www.ncbi.nlm.nih.gov/geo/) → `data/raw/` |
| Reference genome (GRCh38_GCA_000001405.15) | [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/genome/) → `snakemake.genome_dir` |
| `go-basic.obo` | https://current.geneontology.org/ontology/go-basic.obo → `data/go_terms_support/` |
| `gene2go.gz` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz → `data/go_terms_support/` |
| HEK-Blue TLR reporter readouts | `data/supplementary_data/` |

The GO support files are downloaded automatically if absent. Sample labels are parsed from the third underscore-separated token of the FASTQ name; batch `7128` is the main dataset and batch `7086` the external test batch.

## Configuration

`config.json`:

| Key | Description |
| --- | --- |
| `snakemake.sample_dir`, `.genome_dir`, `.work_dir`, `.threads` | Read-preprocessing inputs and resources |
| `paths.featurecounts_dir` | featureCounts tables consumed by the analysis |
| `deseq2.padj_threshold`, `.log2fc_threshold`, `.n_cpus` | Significance calls (default padj < 0.05, \|log2FC\| > 2) |
| `feature_selection.n_estimators`, `.max_depth`, `.max_features` | ExtraTrees selector; `k_best` is not set here but derived from the MI elbow |
| `model_training.random_state` | Seed for every split, selector and classifier |
| `hyperparameter_grids` | Per-model grids searched by the inner CV |
| `ligands.*`, `ligand_aliases` | Subset membership and sample-name normalisation |
| `class_order_for_plotting`, `class_display_names`, `subset_display_names`, `colors`, `subset_palettes` | Figure labelling and palettes |

## Pipeline

**Step 0 — Read preprocessing (Snakemake, `pipeline/`).** FastQC → fastp trimming → STAR alignment → SAMtools merge/sort/index → UMI-tools deduplication → featureCounts gene-level quantification.

**Step 1 — Count table.** featureCounts tables are merged into a samples × genes matrix, samples below 1M reads are dropped, labels are derived from sample names, and the matrix is split into the `train_ligands`, `test_ligands` (batch 7086), `additional_ligands` and `bacterial_ligands` subsets.

**Step 2 — Differential expression.** DESeq2 for each ligand versus `negative_control` within every subset, with volcano plots, clustered heatmaps of the top 50 genes, and per-ligand GO enrichment merged per subset.

**Step 3 — Panel size.** The mutual-information curve is averaged over three seeds and its elbow sets `k_best`. An ExtraTrees importance ranking within those genes is scanned with k-means, scoring ARI against the true ligand labels at every panel size.

**Step 4 — PCA.** Per subset, before and after feature selection, labelled and unlabelled. The selector is fit once on the training subset so all subsets share one gene panel.

**Step 5 — Feature selection versus differential expression.** Venn of the selected panel against the DE gene set, GO enrichment of DE ∩ FS and FS \ DE, and the importance-ranked gene tables.

**Step 6 — Modelling.** For the `main` and `no_flapa` panels: nested stratified CV (5 outer × 3 inner) tunes LinearSVC, SGDClassifier, LogisticRegression, RandomForest and XGBoost on macro F1. Feature selection sits inside the CV pipeline and is refit on each outer training fold only, so no test-fold information reaches the panel and the reported metrics are leakage-free. Class imbalance is handled with balanced sample weights. Hyperparameters are chosen by majority vote across outer folds (ties broken by mean inner then outer F1) and refit on the full panel for three gene sets — `selected_<max_features>`, `de_overlap` (FS ∩ DE) and `union_stable_de` (FS ∪ DE). Each refit is validated on the external batch and applied to the additional and bacterial ligands.

**Step 7 — TLR reporter figure.** HEK-Blue TLR2 (Pam3) and TLR4 (LPS) dose-response with the Fla-PA reference.

## Output

```
results/
├── counts/MATseq_count_summary.csv
├── differential_gene_expression/
│   ├── de_genes_{subset}.csv
│   └── {subset}/{ligand}_deseq2_results.csv
├── go_terms/
│   ├── {de_intersect_fs,fs_only}_go_terms.csv
│   └── {subset}/{{ligand}_go_terms.csv, GO_merged_results.csv}
├── feature_selection/{mutual_information.csv, forest_kmeans.csv}
├── fs_de_genesets/
│   ├── fs_genes_ranked.csv                  # selected genes by ExtraTrees importance
│   ├── fs_gene_names.csv
│   └── selected_vs_de_overlap_table.csv
├── nested_cv/supp_nested_cv_{main,no_flapa}.csv
├── hyperparameter_tuning/{panel}/
│   ├── nested_cv_per_fold.csv
│   ├── oof_predictions.csv                  # pooled out-of-fold predictions
│   ├── selected_params.json
│   ├── {model}_{classification_report,confusion_matrix}.csv
│   └── inner_cv_results/{model}_fold_{n}.csv
├── models/{panel}/{gene_set}/{label_encoder.pkl, {model}.pkl}
├── validation/
│   ├── external_validation_{panel}_performance.csv
│   └── test_set_{panel}/{gene_set}/test_ligands/
├── predictions/{panel}/{gene_set}/{additional,bacterial}_ligands/
│   ├── {model}_{predictions,probabilities}.csv
│   ├── {model}_probabilities_heatmap.png
│   └── test_scores_summary.csv
├── tables/{table2_formatted.csv, Table_2.xlsx, Supplementary_Table_{1..12}.csv}
└── figures/
    ├── deseq2/{subset}/{ligand}_{volcano,histogram}.png
    ├── go/{subset}/{ligand}_go.png
    ├── feature_selection/{mutual_information,forest_ari_sweep}.png
    ├── pca/pca_{subset}[_fs][_labeled].png
    ├── venn/venn_de_vs_fs.png
    ├── model_evaluation/{panel}/Confusion_Matrix_{model}.png
    └── supplementary/tlr_hek_blue.png

```

`{panel}` is `main` or `no_flapa`; `{gene_set}` is `selected_<max_features>`, `de_overlap` or `union_stable_de`.

## Repository layout

```
MATseq/
├── MATseq.py                  # orchestration; defines every output path
├── config.json
├── pipeline/                  # 0_MATseq.smk .. 6_count_reads_featurecounts.smk, environment.yml
├── src/
│   ├── config.py              # config.json loader and path expansion
│   ├── preprocessing.py       # featureCounts loading, filtering, labelling, RPM
│   ├── feature_engineering.py # MI/ExtraTrees selection, MI elbow, ARI sweep
│   ├── pydeseq2.py            # DESeq2 wrapper
│   ├── go_term_analysis.py    # GO enrichment
│   ├── model_training.py      # nested CV tuning and per-gene-set refit
│   ├── prediction.py          # prediction, probabilities, scoring
│   ├── visualization.py       # PCA, volcano, heatmap, GO, Venn, confusion matrices
│   ├── tlr_analysis.py        # HEK-Blue reporter figure
│   ├── make_tables.py         # manuscript tables from raw CSVs
│   └── compose_figures.py     # multi-panel manuscript figures
├── data/
│   ├── raw/                   # FASTQ (required)
│   ├── go_terms_support/      # go-basic.obo, gene2go.gz (required)
│   └── supplementary_data/    # HEK-Blue reporter readouts
└── results/
```

## Citation

If you use this work, please cite:

> **Identifying pyrogenic contaminants using transcriptomic profiling of monocyte activation test with machine learning**
> Tess AV Afanasyeva, Bruno FM de Albuquerque, Paulien Doodeman, Miranda C Dieker-Meijer, Marijke Molenaar-de Backer, Teunis JP van Dam, Anja ten Brinke
> bioRxiv 2025.08.13.670109
> https://doi.org/10.1101/2025.08.13.670109
