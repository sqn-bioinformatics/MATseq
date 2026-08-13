<h1 align="center">
  <img src="assets/logo_MATseq.png" width="180"><br>
  MATseq
</h1>

<p align="center" style="margin-top: 6px;">
  <i>Machine-learning classification of pyrogen-induced monocyte transcriptomic signatures</i>
</p>

## Abstract

The monocyte activation test is an in vitro pyrogenicity assessment method that can utilise human peripheral blood mononuclear cells to detect pyrogens in injectable drugs, providing a binary outcome that indicates the presence or absence of a pyrogen. The added ability to distinguish between different types of pyrogens would greatly expand the applicability of the test, for example, by allowing to pinpoint the source of a contaminating pyrogen in pharmaceutical products. Pyrogens activate a unique set of pattern recognition receptors (PRRs), which contribute to inflammation, yielding distinct transcriptomic activation signatures. In this paper, we capture the unique expression signatures of activated monocytes through bulk RNA sequencing and introduce a data preprocessing pipeline that allows the training of a machine-learning model to classify pyrogenic contaminants. Using a dataset of 108 samples stimulated with five classes of PRR agonists, we could differentiate between these classes with more than 97% F1 on test data. We further demonstrate the model's capacity to generalise on the previously unseen data using different ligands for the same PRRs as well as heat-killed Escherichia coli and Staphylococcus aureus.

MATseq integrates differential expression analysis (DESeq2), feature selection, and machine learning models to classify ligand-induced transcriptomic signatures in monocyte samples. The pipeline processes raw RNA-seq counts through quality control, normalization, feature engineering, and model training to predict ligand classes with high accuracy.


## Installation

### Prerequisites

- Python >=3.10, <3.13
- [Poetry](https://python-poetry.org/) for dependency management


### Setup

```bash
# Clone repository
git clone https://github.com/sqn-bioinformatics/MATseq.git MATseq
cd MATseq

# Install dependencies via Poetry
poetry install

# Validate Snakemake pipeline without running
poetry run python MATseq.py --snakemake dry-run

## Run Snakemake preprocessing
poetry run python MATseq.py --snakemake run

## Run Analysis
poetry run python MATseq.py
```

Optional arguments:
- `--snakemake {dry-run,run}`: Validate (`dry-run`) or execute (`run`) Snakemake preprocessing
- `--genome-dir PATH`: Genome reference directory for Snakemake
- `--fastq-dir PATH`: Raw data (zipped FASTQ) directory for Snakemake
- `--force-recompute`: Skip cache and recompute all steps
- `--cache-dir PATH`: Directory for cached files (default: results/cache)

## Pipeline Steps

### Snakemake Preprocessing

Snakemake pipeline processes raw FASTQ into MATseq_count_summary.csv. It must be run before the remainder of the pipeline. Snakemake will request a path to the raw data in fastq.gz format. Download the raw data from [NCBI GEO database](https://www.ncbi.nlm.nih.gov/geo/), GEO accession number GSE313994. Provide it with --genome-dir or in the config.json. It will also request a path to a reference genome; we used GRCh38_GCA_000001405.15. Download a reference genome from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/). Provide it with --fastq-dir or in the config.json.

Run with `--snakemake run` to execute the RNA-seq preprocessing:
1. **Quality Control** - FastQC on raw reads
2. **Trimming** - Adapter removal with fastp
3. **Alignment** - STAR alignment to reference genome
4. **Merge/Sort/Index** - SAMtools BAM processing
5. **Deduplication** - UMI-tools deduplication
6. **Quantification** - featureCounts gene-level counting

### Configuration

The configuration is managed through `config.json`:

**Snakemake preprocessing:**
- `snakemake.sample_dir`: Raw FASTQ files directory
- `snakemake.work_dir`: Snakemake working directory
- `snakemake.genome_dir`: Reference genome directory
- `snakemake.threads`: Number of threads

**Paths:**
- `paths.data_dir`: Input data directory
- `paths.go_terms_support_dir`: GO enrichment support files directory
- `paths.results_dir`: Output results directory
- `paths.featurecounts_dir`: featureCounts output location

**GO Enrichment Files:**
The following files must be in `data/go_terms_support/`:
- `go-basic.obo` - Download from: https://current.geneontology.org/ontology/go-basic.obo
- `gene2go.gz` - Download from: https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz

MATseq will attempt to download these files automatically if missing (may be slow).

**DESeq2 parameters:**
- `deseq2.padj_threshold`: Adjusted p-value threshold (default: 0.05)
- `deseq2.log2fc_threshold`: Log2 fold-change threshold (default: 2)
- `deseq2.n_cpus`: Number of CPUs for parallel processing

**Feature selection:**
- `feature_selection.k_best`: Number of top genes kept by mutual-information selection
- `feature_selection.n_estimators`: Trees in the ExtraTrees selector
- `feature_selection.max_depth`: Maximum depth of the ExtraTrees selector
- `feature_selection.max_features`: Maximum genes kept by the ExtraTrees selector

**Model training:**
- `model_training.random_state`: Reproducibility seed for training and CV splits
- `model_factory.calibrate`: Probability-calibrate the LinearSVC (CalibratedClassifierCV)
- `hyperparameter_grids`: Per-model search grids used by the nested cross-validation

### Main Pipeline

The pipeline runs as eight steps in `MATseq.py` (`run_pipeline`). Intermediate
results (counts, DESeq2 contrasts, feature-selection runs, tuned models) are
cached under `results/cache`; rerun with `--force-recompute` to ignore the cache.

1. **Data Preprocessing** (`prepare_counts`)
   - Merge featureCounts outputs into a samples × genes matrix
   - Filter samples by read-count threshold (>1M reads)
   - Derive ligand labels from sample names (with `ligand_aliases`)
   - Write `counts/MATseq_count_summary.csv`

2. **DESeq2 Differential Expression** (per subset: main, additional, bacteria)
   - PCA plot per subset (RPM → log1p → standardised), labelled and unlabelled
   - Each ligand vs negative control; significant genes (padj < 0.05, |log2FC| > 2)
   - Volcano plots and clustered heatmaps per ligand
   - GO enrichment per ligand, merged into `GO_merged_results.csv`

3. **Feature-selection vs DE Venn** (training subset)
   - Run the selection pipeline 1000× with different seeds (parallelised across cores)
   - Compare the union of selected genes against the DESeq2 DE genes (Venn diagram)
   - Gene-frequency table; GO enrichment on DE ∩ FS and FS \ DE gene sets

4. **Nested CV tuning + deployment refit** (main panel, with and without Fla-PA)
   - Nested stratified CV (5 outer, 3 inner) tuning LinearSVC, SGDClassifier, LogisticRegression, RandomForest, XGBoost
   - Feature selection is embedded in the cross-validation pipeline and re-fit on each outer training fold only (no test-fold information reaches selection), so reported metrics are leakage-free
   - Inner `GridSearchCV` tunes classifier hyperparameters on macro F1
   - Class imbalance handled with balanced class/sample weights (no SMOTE)
   - Deployment hyperparameters chosen by majority vote across outer folds, then refit on the full panel
   - Writes per-fold and pooled out-of-fold metrics, confusion matrices, and `selected_params.json`
   - A second model set is trained on the main panel excluding Fla-PA

4c. **Table 2 feature-set benchmark** (`ModelTrainer.benchmark_feature_set_conditions`)
   - Re-runs the leakage-free nested CV (5 outer, 3 inner) for every model across four feature-set conditions: all genes, feature selection, feature selection + DE, and a randomly selected gene set of matching size
   - The gene set for each condition is derived on the outer-training fold only (feature selection and DESeq2 re-run per fold), frozen, then scored once on the outer-test fold
   - Writes the raw long-form `results/tables/table2_feature_set_benchmark.csv`; `src/make_tables.py` (`format_table2`) reshapes it into the publication table (`table2_formatted.csv` and `table2.tex`)

5. **Model validation on external test batch**
   - Apply the deployed models to an independent sequencing batch (`test.work_dir/featurecounts`)
   - Predictions, probabilities and probability heatmaps, with and without Fla-PA
   - PCA before and after feature selection (reusing the main-panel palette and class order)
   - DESeq2 of each ligand vs negative control on the external batch
   - Skipped automatically if no external featureCounts are present

6. **Prediction on additional and bacterial ligands**
   - Apply the deployed main-panel models to held-out ligands (LTA, MPLA, Pam2) and heat-killed bacteria
   - Generate predictions, probabilities and probability heatmaps

7. **TLR Reporter Visualization** (`src/tlr_analysis.py`)
   - HEK-Blue TLR2 (Pam3) and TLR4 (LPS) dose-response with Fla-PA reference

8. **Prediction without Fla-PA**
   - Apply the no-Fla-PA models to the additional and bacterial ligands


## Output Structure

```
results/
├── cache/                                  # Cached intermediates + manifest.json
├── counts/
│   └── MATseq_count_summary.csv
├── differential_gene_expression/
│   └── {ligand}_deseq2_results.csv
├── go_terms/
│   ├── {ligand}_go_terms.csv               # Per-ligand GO enrichment
│   ├── GO_merged_results.csv               # All per-ligand terms merged
│   ├── de_intersect_fs_go_terms.csv        # GO for DE ∩ FS gene set
│   └── fs_only_go_terms.csv                # GO for FS \ DE gene set
├── tables/
│   ├── supp_nested_cv_main.csv             # Nested CV summary (main panel)
│   ├── supp_nested_cv_no_flapa.csv         # Nested CV summary (no Fla-PA)
│   ├── table2_feature_set_benchmark.csv    # Raw feature-set benchmark (long form)
│   ├── table2_formatted.csv                # Publication Table 2 (mean ± SD grid)
│   └── table2.tex                          # Publication Table 2 (booktabs LaTeX)
├── feature_analysis/
│   └── gene_frequency_table.csv            # Gene selection frequency across 1000 runs
├── models/                                 # Deployed main-panel models
│   ├── label_encoder.pkl
│   ├── {model}.pkl
│   └── no_flapa/                           # Models trained without Fla-PA
│       ├── label_encoder.pkl
│       └── {model}.pkl
├── hyperparameter_tuning/                  # Nested CV outputs (main panel)
│   ├── nested_cv_per_fold.csv
│   ├── nested_cv_summary.csv
│   ├── {prefix}oof_predictions.csv         # Pooled out-of-fold predictions
│   ├── {prefix}{model}_classification_report.csv
│   ├── selected_params.json                # Deployment hyperparameters
│   └── inner_cv_results/{model}_fold_{n}.csv
├── hyperparameter_tuning_no_flapa/         # Same, for the no-Fla-PA models
├── figures/
│   ├── deseq2/
│   │   ├── {ligand}_volcano.png
│   │   └── {ligand}_histogram.png
│   ├── go/
│   │   ├── {ligand}_go.png
│   │   ├── de_intersect_fs_go.png
│   │   └── fs_only_go.png
│   ├── pca/
│   │   ├── {subset}_pca.png
│   │   └── {subset}_pca_labeled.png
│   ├── venn/
│   │   └── venn_de_vs_fs.png
│   ├── model_evaluation/
│   │   ├── {prefix}{model}_confusion_matrix.csv
│   │   ├── {prefix}{model}_confusion_matrix_normalized.csv
│   │   └── Confusion_Matrix_{prefix}{model}.png
│   └── supplementary/
│       └── tlr_hek_blue.png
├── validation/{test_name}/                 # External test batch (e.g. 7086)
│   ├── main_ligands/
│   │   ├── {model}_predictions.csv
│   │   ├── {model}_probabilities.csv
│   │   └── {model}_probabilities_heatmap.png
│   └── no_flapa/
└── predictions/
    ├── additional_ligands/
    │   ├── {model}_predictions.csv
    │   ├── {model}_probabilities.csv
    │   └── {model}_probabilities_heatmap.png
    ├── bacteria_ligands/
    └── main_ligands_no_flapa/              # No-Fla-PA model predictions
        ├── additional_ligands/
        └── bacteria_ligands/
```

## Project Structure

```
MATseq/
├── MATseq.py           
├── config.json         
├── pyproject.toml      
├── pipeline/           
│   ├── 0_MATseq.smk                       
│   ├── 1_control_quality_fastqc.smk       
│   ├── 2_trim_fastp.smk                   
│   ├── 3_align_star.smk                   
│   ├── 4_merge_sort_index_samtools.smk    
│   ├── 5_deduplicate_umitools.smk         
│   ├── 6_count_reads_featurecounts.smk    
│   └── environment.yml                    
├── src/
│   ├── cache.py               
│   ├── config.py              
│   ├── preprocessing.py       
│   ├── feature_engineering.py 
│   ├── model_training.py      
│   ├── make_tables.py          # Formats raw pipeline CSVs into publication tables
│   ├── prediction.py          
│   ├── pydeseq2.py           
│   ├── visualization.py       
│   └── go_term_analysis.py    
├── scripts/            # Standalone helper/exploratory scripts (PCA comparisons, figure composition)
├── data/               
│   ├── go_terms_support/      
│   │   ├── go-basic.obo                                                        # (required)
│   │   ├── gene2go.gz                                                          # (required)
│   │   └──  gene_result_ncbi_human_proteincoding.txt 
│   ├── reference_genome/
│   │   └── GRCh38_GCA_000001405.15/
│   │   │   ├── GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation     # (required)
│   │   │   └── GCA_000001405.15_GRCh38_no_alt_analysis_set                     # (required)
│   └── raw/                                                                    # FASTQ (required)
└── results/            
```

## Citation

If you use this work, please cite:

> **Identifying pyrogenic contaminants using transcriptomic profiling of monocyte activation test with machine learning**
> Tess AV Afanasyeva, Bruno FM de Albuquerque, Paulien Doodeman, Miranda C Dieker-Meijer, Marijke Molenaar-de Backer, Teunis JP van Dam, Anja ten Brinke
> bioRxiv 2025.08.13.670109
> https://doi.org/10.1101/2025.08.13.670109
