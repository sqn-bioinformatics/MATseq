#### Sanquin Bioinformatics ####
# name: MATseq pipeline
# description: Monocyte activation test transcriptomes analysis pipeline
# author: Tess Afanasyeva
# date: 2024-01-01


rule fastp:
# fastp - A tool designed to provide fast all-in-one preprocessing for FastQ files.
# https://github.com/OpenGene/fastp#quality-filter
    input:
        R1 =  os.path.join(config["WorkDir"], "fastqs/{sample}_R1.fastq.gz"),
        R2 = os.path.join(config["WorkDir"], "fastqs/{sample}_R2.fastq.gz")
    output:
        R1_tp = "fastp/reads/{sample}_R1_trimmed.fastq",
        R2_tp = "fastp/reads/{sample}_R2_trimmed.fastq",
        json = "fastp/report/{sample}.json",
        html = "fastp/report/{sample}.html",
    benchmark: "benchmarks/fastp_{sample}.txt"
    priority: 50
    params:
        length_required = 36,
        quality_score = 15,
    threads: 16
    resources:
        mem_mb= 5000
    log: "logs/fastp/{sample}.log"
    shell:
        "fastp \
        --in1 {input.R1} \
        --in2 {input.R2} \
        --out1 {output.R1_tp} \
        --out2 {output.R2_tp} \
        --html {output.html} \
        --json {output.json} \
        --thread {threads} \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --trim_poly_x \
        --length_required {params.length_required} \
        --low_complexity_filter \
        --cut_tail \
        --cut_tail_window_size 4 \
        --cut_tail_mean_quality {params.quality_score} \
        --dont_eval_duplication "

