# !/bin/bash

function get_yaml_value () {
  echo $(grep ${1} ${CONFIG} | cut -d ' ' -f2 | tr -d "\'\"")
}

CONFIG=config.yml
WORKDIR=$(get_yaml_value "^WorkDir:")
THREADS=$(get_yaml_value "^Threads:")

snakemake --use-conda \
	        --directory ${WORKDIR} \
          --cores ${THREADS} \
          --configfile ${CONFIG} \
          --snakefile workflow/rules/0_MATseq.smk \
          --latency-wait 60 \
          --rerun-incomplete \
          
 

 
        
          
  
   

          # --rerun-incomplete
          # --list \
          # --summary \
          # --debug-dag \
          # --dry-run \
          # --keep-going \
          # --batch RULE=BATCH/BATCHES
          # -F \
          # --lint  
          # --dag | dot -Tsvg > dag.svg  \
          





 