# from glob import glob
# import os
# import re

# files = glob(config["SampleDir"]+"/*.fastq.gz")

# file_names = sorted(list(set([re.split("_R1|_R2", os.path.basename(file))[0] for file in files])))
 
# start here maybe not needed step
# https://github.com/sqn-bioinformatics/snakemake-workflows/blob/main/rnaseq_pipeline/rules/qc.smk


# rule all:
#   input:
# 	  expand("temp/raw_fastq/{sample}", sample=file_names)
        
# rule gunzip:
#   input:
#     "experiment/raw_data/{sample}.gz"
#   output:
#     "temp/raw_fastq/{sample}"  
#   shell:
#     "gunzip -c {input} > {output}"


# rule all:
# 	input:
# 		expand("experiment/results/fastp/{sample}.html", sample=file_names)

# rule fastqc:
#     input:
#         read=config["SampleDir"]+"/{sample}.fastq.gz"
#     output:
#         html="experiment/results/fastqc/{sample}.html",
#         zip="experiment/results/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#     shell:
#         "mkdir -p fastqc/tmp && \
#         --fastqc \
#         --threads {threads} \
#         --outdir fastqc \
#         --dir fastqc/tmp \
#         {input.read}"

#     params: "--quiet"
#     conda:
#         "environment.yml"
#     threads: 1
#     resources:
#     mem_mb=2000
#     log:
#         "logs/fastp/{sample}.log"
#     wrapper:
#         "v1.7.1/bio/fastqc"




#### Sanquin Bioinformatics ####
# description: RNA-seq pipeline
# author: Margo Schuller
# date: 2021-03-10


rule fastqc:
# FastQC - A high throughput sequence QC analysis tool
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    input:
        read = config["SampleDir"]+"/{sample}_{paired}.fastq.gz"
    output:
        html = "fastqc/{sample}_{paired}_fastqc.html",
        zip = "fastqc/{sample}_{paired}_fastqc.zip"
    # conda:
    #     "environment.yml"
    threads: 1
    resources:
        mem_mb = 2000
    message:
        "Rule {rule} checking the quality of the .fastq file."
    shell:
        "mkdir -p fastqc/tmp && \
         fastqc \
         --threads {threads} \
         --outdir fastqc \
         --dir fastqc/tmp \
         {input.read}"


# use rule fastqc as fastqc_fastp with:
#     input:
#         read="fastp/reads/{sample}_{paired}_trimmed_paired.fastq.gz"
#     output:
#         html="fastqc/{sample}_{paired}_trimmed_paired_fastqc.html",
#         zip="fastqc/{sample}_{paired}_trimmed_paired_fastqc.zip"


# rule qualimap:
# # qualimap - RNA-seq QC reports quality control metrics and bias estimations
# # which are specific for whole transcriptome sequencing
# # http://qualimap.conesalab.org/
#     input:
#         sam="star/{sample}_Aligned.sortedByCoord.out.bam",
#         GTFfile=glob(config["GenomeDir"]+"/*.gtf")
#     output:
#         counts="qualimap/{sample}.counts"
#     conda:
#         "environment.yml"
#     threads: 1
#     resources:
#         mem_mb=8000
#     shell:
#         "unset DISPLAY && \
# 	 qualimap rnaseq \
#          -bam {input.sam} \
#          -gtf {input.GTFfile} \
#          -oc  {output.counts} \
# 	 -pe \
#          -outdir qualimap/{wildcards.sample}/ \
# 	 --java-mem-size={resources.mem_mb}M"