#!/bin/bash

function collect_yaml_value () {
  echo $(grep ${1} ${CONFIG} | cut -d ' ' -f2 | tr -d "\'\"")
}

CONFIG=config.yml
WORKDIR=$(collect_yaml_value "^WorkDir:")
THREADS=$(collect_yaml_value "^Threads:")

snakemake --use-conda \
	        --directory ${WORKDIR} \
          --cores ${THREADS} \
          --configfile ${CONFIG} \
          --snakefile library/scripts/0_MATseq.smk \
          --delete-all-output
 