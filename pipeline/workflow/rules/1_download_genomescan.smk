#### Sanquin Bioinformatics ####
# name: MATseq pipeline
# description: Monocyte activation test transcriptomes analysis pipeline
# author: Tess Afanasyeva
# date: 2024-01-01


def required_files(wildcards):
    FILES = glob(config["SampleDir"]+"/*.fastq.gz")
    for file in FILES:
        if wildcards.sample in file and wildcards.paired in file:
            return file


rule copy_and_rename:
	input: required_files
	output: "fastqs/{sample}_{paired}.fastq.gz"
	shell:"cp -v {input} {output}"


rule fastqc:
# FastQC - A high throughput sequence QC analysis tool
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    input:
        read = "fastqs/{sample}_{paired}.fastq.gz"
    output:
        html = "fastqc/{sample}_{paired}.html",
        zip = "fastqc/{sample}_{paired}.zip"
    threads: 18
    resources:
        mem_mb= 32000
    shell:
        "fastqc \
        --threads {threads} \
        --outdir fastqc \
        {input.read}"