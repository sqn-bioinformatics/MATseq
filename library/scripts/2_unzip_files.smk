from glob import glob
import os
import re

files = glob(config["SampleDir"]+'/*.fastq.gz')

file_names = sorted(list(set([re.split('_R1|_R2', os.path.basename(file))[0] for file in files])))
 
# start here maybe not needed step
# https://github.com/sqn-bioinformatics/snakemake-workflows/blob/main/rnaseq_pipeline/rules/qc.smk


# rule all:
#   input:
# 	  expand("temp/raw_fastq/{sample}", sample=file_names)
        
rule gunzip:
  input:
    "experiment/raw_data/{sample}.gz"
  output:
    "temp/raw_fastq/{sample}"  
  shell:
    "gunzip -c {input} > {output}"