<h1 align="center">
  <img src="assets/logo_MATseq.png" width="180"><br>
  MATseq
</h1>

<p align="center" style="margin-top: 6px;">
  <i>Machine-learning classification of pyrogen-induced monocyte transcriptomic signatures</i>
</p>

## Abstract

The monocyte activation test is an in vitro pyrogenicity assessment method that can utilise human peripheral blood mononuclear cells to detect pyrogens in injectable drugs, providing a binary outcome that indicates the presence or absence of a pyrogen. The added ability to distinguish between different types of pyrogens would greatly expand the applicability of the test, for example, by allowing to pinpoint the source of a contaminating pyrogen in pharmaceutical products. Pyrogens activate a unique set of pattern recognition receptors (PRRs), which contribute to inflammation, yielding distinct transcriptomic activation signatures. In this paper, we capture the unique expression signatures of activated monocytes through bulk RNA sequencing and introduce a data preprocessing pipeline that allows the training of a machine-learning model to classify pyrogenic contaminants. We further demonstrate the model's capacity to generalise on the previously unseen data using different ligands for the same PRRs as well as heat-killed *Escherichia coli* and *Staphylococcus aureus*.

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


## Input data

| Input | Source |
| --- | --- |
| Raw reads (`fastq.gz`) | [GEO GSE313994](https://www.ncbi.nlm.nih.gov/geo/) → `data/raw/` |
| Reference genome (GRCh38_GCA_000001405.15) | [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/genome/) → `snakemake.genome_dir` |
| `go-basic.obo` | https://current.geneontology.org/ontology/go-basic.obo → `data/go_terms_support/` |
| `gene2go.gz` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz → `data/go_terms_support/` |
| HEK-Blue TLR reporter readouts | `data/supplementary_data/` |

The GO support files are downloaded automatically if absent. Sample labels are parsed from the third underscore-separated token of the FASTQ name; batch `7128` is the main dataset and batch `7086` the external test batch.

## Citation

If you use this work, please cite:

> **Identifying pyrogenic contaminants using transcriptomic profiling of monocyte activation test with machine learning**
> Tess AV Afanasyeva, Bruno FM de Albuquerque, Paulien Doodeman, Miranda C Dieker-Meijer, Marijke Molenaar-de Backer, Teunis JP van Dam, Anja ten Brinke
> bioRxiv 2025.08.13.670109
> https://doi.org/10.1101/2025.08.13.670109
