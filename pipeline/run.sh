#!/bin/bash

set -e

CONFIG="../config.json"

WORKDIR=$(python3 -c "import json; print(json.load(open('$CONFIG'))['snakemake']['work_dir'])" | sed "s|~|$HOME|")
THREADS=$(python3 -c "import json; print(json.load(open('$CONFIG'))['snakemake']['threads'])")
SAMPLEDIR=$(python3 -c "import json; print(json.load(open('$CONFIG'))['snakemake']['sample_dir'])" | sed "s|~|$HOME|")
GENOMEDIR=$(python3 -c "import json; print(json.load(open('$CONFIG'))['snakemake']['genome_dir'])" | sed "s|~|$HOME|")

if [ -z "$GENOMEDIR" ] || [ ! -d "$GENOMEDIR" ]; then
    echo "Error: genome_dir must be set in config.json and point to a valid directory"
    echo "Current value: '$GENOMEDIR'"
    echo "Please set snakemake.genome_dir in config.json to your reference genome path"
    exit 1
fi

if [ -z "$SAMPLEDIR" ] || [ ! -d "$SAMPLEDIR" ]; then
    echo "Error: sample_dir must be set in config.json and point to a valid directory"
    echo "Current value: '$SAMPLEDIR'"
    echo "Please set snakemake.sample_dir in config.json to your FASTQ files path"
    exit 1
fi

snakemake --directory ${WORKDIR} \
          --cores ${THREADS} \
          --snakefile workflow/rules/0_MATseq.smk \
          --config SampleDir=${SAMPLEDIR} GenomeDir=${GENOMEDIR} \
          --latency-wait 60 \
          --rerun-incomplete
