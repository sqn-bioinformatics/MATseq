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

1. **Data Preprocessing** (`prepare_counts`)
   - Merge featureCounts outputs
   - Filter samples by read count threshold (>1M reads)
   - Extract main, additional, and bacteria ligand subsets with labels
   - Generate PCA plots for each subset

2. **DESeq2 Differential Expression** (`deseq2_training`)
   - Differential expression analysis on training data
   - Identify significant genes (padj < 0.05, |log2FC| > 2)
   - Generate volcano plots and heatmaps per ligand
   - Run GO enrichment per ligand; merge into single GO table

3. **Feature Selection Analysis** (`venn_diagram`)
   - Run feature selection 1000x with different random seeds
   - Compare feature-selected genes vs DESeq2 differentially expressed genes
   - Generate Venn diagrams and gene frequency tables
   - Run GO enrichment on DE ∩ FS and FS \ DE gene sets

4. **Model Training** (`model_training`)
   - Nested cross-validation (5 outer, 3 inner folds) tuning LinearSVC, LogisticRegression, RandomForest, and XGBoost
   - Balance classes with class/sample weights (no SMOTE)
   - Select deployment hyperparameters by majority vote and refit on the full panel
   - Write per-fold and pooled out-of-fold metrics, confusion matrices, and `selected_params.json`

5. **Additional Ligand and Bacterial Analysis** (`deseq2_other_ligands`, `deseq2_bacterial`)
   - DESeq2 analysis on additional ligands (LTA, MPLA, Pam2)
   - DESeq2 analysis on bacterial samples (HK E.coli, HK S.aureus)
   - Generate visualization outputs

6. **Prediction on New Data** (`predict_other_ligands`, `predict_bacterial`)
   - Apply trained models to additional ligands and bacterial samples
   - Generate probability heatmaps and predictions

7. **TLR Reporter Visualization** (`tlr_hek_blue`)
   - Loading and plotting handled in `src/tlr_analysis.py`
   (`load_tlr_data`, `plot_tlr_hek_blue`, `plot_tlr_panel`)

8. **Extended Training** (`training_wo_flapa`)
   - Retrain models on dataset without Fla-Pa samples


## Output Structure

```
results/
├── cache/                             
├── counts/
│   └── MATseq_count_summary.csv
├── subsets/
│   ├── training_data_with_labels.csv
│   ├── other_ligands_data_with_labels.csv
│   ├── bacterial_data_with_labels.csv
│   └── training_wo_flapa_data_with_labels.csv
├── differential_gene_expression/
│   └── {ligand}_deseq2_results.csv
├── go_terms/
│   ├── {ligand}_go_terms.csv           # Per-ligand GO enrichment
│   ├── GO_merged_results.csv           # All per-ligand terms merged
│   ├── de_intersect_fs_go_terms.csv    # GO for DE ∩ FS gene set
│   └── fs_only_go_terms.csv            # GO for FS \ DE gene set
├── feature_analysis/
│   ├── gene_frequency_table.csv        # Gene selection frequency across 1000 runs
│   └── list_of_gene_name_set_1000      # Pickle: all 1000 feature sets
├── models/
│   ├── label_encoder.pkl
│   └── {model_name}.pkl
├── model_evaluation/
│   ├── full_model_evaluation.csv           # Metrics on all genes
│   ├── full_model_evaluation_per_fold.csv
│   ├── fs_model_evaluation.csv             # Metrics on feature-selected genes
│   ├── fs_model_evaluation_per_fold.csv
│   ├── fs_de_model_evaluation.csv          # Metrics on FS ∪ DE genes
│   └── fs_de_model_evaluation_per_fold.csv
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
│   │   ├── {subset}_pca_labeled.png
│   │   ├── {subset}_pca_selected.png
│   │   └── {subset}_pca_selected_labeled.png
│   ├── venn/
│   │   └── venn_de_vs_fs.png
│   ├── model_evaluation/
│   │   └── Confusion_Matrix_{prefix}{model}_fold_{n}.png
│   └── supplementary/
│       └── tlr_hek_blue.png
└── predictions/
    ├── other_ligands/
    │   ├── {model}_predictions.csv
    │   ├── {model}_probabilities.csv
    │   └── {model}_probabilities_heatmap.png
    ├── bacterial/
    │   ├── {model}_predictions.csv
    │   ├── {model}_probabilities.csv
    │   └── {model}_probabilities_heatmap.png
    └── training_wo_flapa/
        ├── other_ligands/
        └── bacterial/
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
│   ├── feature_analysis.py    
│   ├── model_training.py      
│   ├── prediction.py          
│   ├── pydeseq2.py           
│   ├── visualization.py       
│   └── go_term_analysis.py    
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
