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

## Usage

### Running the Pipeline

Execute the main pipeline:
```bash
poetry run python MATseq.py
```

Optional arguments:

- `--run-snakemake`: Run Snakemake preprocessing before pipeline
- `--snakemake-dir PATH`: Directory containing Snakemake pipeline
- `--fastq-dir PATH`: Raw data (zipped FASTQ) directory for Snakemake
- `--genome-dir PATH`: Genome reference directory for Snakemake
- `--snakemake-dry-run`: Validate Snakemake pipeline without running
- `--force-recompute`: Skip cache and recompute all steps
- `--cache-dir PATH`: Directory for cached files (default: results/cache)


### Configuration

The configuration is managed through `config.json`:

**Snakemake preprocessing:**
- `snakemake.sample_dir`: Raw FASTQ files directory
- `snakemake.work_dir`: Snakemake working directory
- `snakemake.genome_dir`: Reference genome directory
- `snakemake.threads`: Number of threads

**Paths:**
- `paths.data_dir`: Input data directory
- `paths.deseq2_dir`: DESeq2 reference files (gene2go, go-basic.obo should download automatically)
- `paths.results_dir`: Output results directory
- `paths.featurecounts_dir`: featureCounts output location

**DESeq2 parameters:**
- `deseq2.padj_threshold`: Adjusted p-value threshold (default: 0.05)
- `deseq2.log2fc_threshold`: Log2 fold-change threshold (default: 2)
- `deseq2.n_cpus`: Number of CPUs for parallel processing

**Feature selection:**
- `feature_selection.k_best`: Number of top features to select
- `feature_selection.n_estimators`: Random forest estimators
- `feature_selection.random_state`: Reproducibility seed

**Model training:**
- `model_training.apply_smote`: Enable SMOTE oversampling (used in this paper)
- `model_training.smote_k_neighbors`: SMOTE k-neighbors parameter

## Pipeline Steps

### Snakemake Preprocessing
Snakemake will request a path to the raw data in fastq.gz format. Download the raw data from NCBI, GEO accession number GSE313994. Provide it with --genome-dir or in the config.json.
Snakemake will request a path to a reference genome; we used GRCh38_GCA_000001405.15. Provide it with --fastq-dir or in the config.json.

Run with `--run-snakemake` to execute the RNA-seq preprocessing:
1. **Quality Control** - FastQC on raw reads
2. **Trimming** - Adapter removal with fastp
3. **Alignment** - STAR alignment to reference genome
4. **Merge/Sort/Index** - SAMtools BAM processing
5. **Deduplication** - UMI-tools deduplication
6. **Quantification** - featureCounts gene-level counting

### Main Pipeline

1. **Data Preprocessing** (`merge_counts`)
   - Merge featureCounts outputs
   - Filter samples by read count threshold (>1M reads)
   - Extract training/test/bacterial subsets with proper labels

2. **DESeq2 Differential Expression** (`deseq2_training`)
   - Differential expression analysis on training data
   - Identify significant genes (padj < 0.05, |log2FC| > 2)
   - Generate volcano plots, heatmaps, PCA, GO enrichment

3. **Feature Selection Analysis** (`venn_diagram`)
   - Run feature selection 1000x with different random seeds
   - Compare feature-selected genes vs DESeq2 differentially expressed genes
   - Generate Venn diagrams and gene frequency tables

4. **Model Training** (`model_training`)
   - Train classification models (LinearSVC, SGDClassifier, RandomForest, XGBoost)
   - Apply SMOTE for class balancing
   - Evaluate with cross-validation
   - Generate performance metrics and confusion matrices

5. **Additional Ligand and Bacterial Analysis** (`deseq2_other_ligands`, `deseq2_bacterial`)
   - DESeq2 analysis on additional ligands (LTA, MPLA, Pam2)
   - DESeq2 analysis on bacterial samples (HK E.coli, HK S.aureus)
   - Generate visualization outputs

6. **Prediction on New Data** (`predict_other_ligands`, `predict_bacterial`)
   - Apply trained models to additional ligands and bacterial samples
   - Generate probability heatmaps and predictions

7. **TLR Reporter Visualization** (`tlr_hek_blue`)
   - Supplementary Figure 2

8. **Extended Training** (`training_wo_flapa`)
   - Retrain models on dataset without Fla-Pa samples
   - Supplementary Figure 3


## Project Structure

```
MATseq/
├── MATseq.py           
├── config.json         
├── pyproject.toml      # Poetry dependencies
├── pipeline/           # Snakemake preprocessing
│   ├── 0_MATseq.smk                       # Main workflow
│   ├── 1_control_quality_fastqc.smk       # FastQC quality control
│   ├── 2_trim_fastp.smk                   # Fastp adapter trimming
│   ├── 3_align_star.smk                   # STAR alignment
│   ├── 4_merge_sort_index_samtools.smk    # SAMtools processing
│   ├── 5_deduplicate_umitools.smk         # UMI-tools deduplication
│   ├── 6_count_reads_featurecounts.smk    # featureCounts quantification
│   └── environment.yml                    # Conda environment spec
├── src/
│   ├── cache.py               # Caching system
│   ├── config.py              # Configuration loading
│   ├── preprocessing.py       # Data loading and filtering
│   ├── feature_engineering.py # Feature selection pipeline
│   ├── feature_analysis.py    # Feature selection analysis
│   ├── model_training.py      # ML model training
│   ├── prediction.py          # Model prediction
│   ├── pydeseq2.py           # DESeq2 wrapper
│   └── visualization.py       # Plotting functions
├── data/               # Input data
└── results/            # Output directory (generated)
```


## Output Structure

```
results/
├── cache/                          # Cached pipeline step outputs
├── counts/                         # Processed count matrices
│   └── MATseq_count_summary.csv
├── subsets/                        # Training/test data subsets
│   ├── training_data_with_labels.csv
│   ├── other_ligands_data_with_labels.csv
│   └── bacterial_data_with_labels.csv
├── differential_gene_expression/   # DESeq2 results
│   ├── {ligand}_deseq2_results.csv
│   └── merged_results.csv
├── figures/
│   ├── deseq2/                     # DESeq2 visualizations
│   │   ├── {ligand}_volcano.png
│   │   ├── {ligand}_histogram.png
│   │   └── {ligand}_go.png
│   ├── pca/                        # PCA plots
│   ├── venn/                       # Venn diagrams
│   └── supplementary/              # Step 8: TLR reporter and gene expression
│       └── tlr_hek_blue.png
├── go_terms/                       # GO enrichment results
│   └── {ligand}_go_terms.csv
├── models/                         # Trained models
├── model_evaluation/               # Model performance metrics
│   ├── confusion_matrices/
│   └── evaluation_metrics.csv
└── predictions/                    
    ├── other_ligands/             
    │   ├── {model}_predictions.csv
    │   ├── {model}_probabilities.csv
    │   └── {model}_probabilities_heatmap.png
    ├── bacterial/                 
    │   ├── {model}_predictions.csv
    │   ├── {model}_probabilities.csv
    │   └── {model}_probabilities_heatmap.png
    └── training_wo_flapa/          # Step 7
        ├── other_ligands/
        └── bacterial/
```


## Citation

If you use this work, please cite:

> **Identifying pyrogenic contaminants using transcriptomic profiling of monocyte activation test with machine learning**
> Tess AV Afanasyeva, Bruno FM de Albuquerque, Paulien Doodeman, Miranda C Dieker-Meijer, Marijke Molenaar-de Backer, Teunis JP van Dam, Anja ten Brinke
> bioRxiv 2025.08.13.670109
> https://doi.org/10.1101/2025.08.13.670109
