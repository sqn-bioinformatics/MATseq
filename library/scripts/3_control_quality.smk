import os

rule all:
	input:
		expand("experiment/results/fastp/{sample}.html", sample=file_names)

rule fastqc:
    input:
        read=config['SampleDir']+'/{sample}.fastq.gz'
    output:
        html="experiment/results/fastqc/{sample}.html",
        zip="experiment/results/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    shell:
        'mkdir -p fastqc/tmp && \
        --fastqc \
        --threads {threads} \
        --outdir fastqc \
        --dir fastqc/tmp \
        {input.read}'

    params: "--quiet"
    conda:
        'environment.yml'
    threads: 1
    resources:
    mem_mb=2000
    log:
        "logs/fastp/{sample}.log"
    wrapper:
        "v1.7.1/bio/fastqc"



   input:
        read=config['SampleDir']+'/{sample}_{paired}.fastq.gz'
    output:
        html='fastqc/{sample}_{paired}_fastqc.html',
        zip='fastqc/{sample}_{paired}_fastqc.zip'
    conda:
        'environment.yml'
    threads: 1
    resources:
        mem_mb=2000
    shell:
        'mkdir -p fastqc/tmp && \
         fastqc \
         --threads {threads} \
         --outdir fastqc \
         --dir fastqc/tmp \
         {input.read}'

use rule fastqc as fastqc_fastp with:
    input:
        read='fastp/reads/{sample}_{paired}_trimmed_paired.fastq.gz'
    output:
        html='fastqc/{sample}_{paired}_trimmed_paired_fastqc.html',
        zip='fastqc/{sample}_{paired}_trimmed_paired_fastqc.zip'