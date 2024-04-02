#### Sanquin Bioinformatics ####
# name: MATseq pipeline
# description: Monocyte activation test transcriptomes analysis pipeline
# author: Tess Afanasyeva
# date: 2024-01-01


def required_files(wildcards):
    FILES = glob(config["SampleDir"]+"/*.fastq.gz") # HGHVMDSX3_104972-001-010_R for whatever reason is doubled, will cause program fail if run
    for file in FILES:
        if wildcards.sample in file and wildcards.paired in file:
            print(f'Found {file}')
            return file
        else:
            print(f"No: {wildcards.sample}, {wildcards.paired }")

rule copy_and_rename:
    input: required_files
    output: config["WorkDir"]+"/fastqs/{sample}_{paired}.fastq.gz"
    benchmark: "benchmarks/copy_and_rename_{sample}_{paired}.txt"
    threads: 16
    shell: 'cp -v "{input}" {output}'


rule fastqc:
# FastQC - A high throughput sequence QC analysis tool
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    input:
        read = rules.copy_and_rename.output
    output:
        html = "fastqc/{sample}_{paired}_fastqc.html",
        zip = "fastqc/{sample}_{paired}_fastqc.zip"
    benchmark: "benchmarks/fastqc_{sample}_{paired}_.txt"
    threads: 18
    resources: mem_mb= 2000
    shell:
        'fastqc \
        --threads {threads} \
        --outdir fastqc \
        {input.read}'
