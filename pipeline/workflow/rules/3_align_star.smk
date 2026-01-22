#### Sanquin Bioinformatics ####
# name: MATseq pipeline
# description: Monocyte activation test transcriptomes analysis pipeline
# author: Tess Afanasyeva
# date: 2024-01-01

rule star_index:
# STAR - Spliced Transcripts Alignment to a Reference (c) Alexander Dobin, 2009-2020
# https://github.com/alexdobin/STAR
    input:
        fasta=glob(config["GenomeDir"]+"/*.fna"),
        annotation=glob(config["GenomeDir"]+"/*.gtf")
    output:
        genome_dir=directory(config["GenomeDir"]+"/index")
    log: "sm_logs/star_index/star_index.log"
    message: "Rule {rule} is indexing reference genome."
    benchmark: "sm_benchmarks/star_index.txt"
    priority: 100
    params: overhang = 100
    threads: 8
    resources:
        mem_mb=32000
    shell:
        'STAR \
        --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir {output.genome_dir} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.annotation} \
        --sjdbOverhang {params.overhang}'
    
checkpoint star_alignment:
    input:
        R1="sm_fastp/reads/{sample}_R1_trimmed.fastq",
        R2="sm_fastp/reads/{sample}_R2_trimmed.fastq",
        annotation=glob(config["GenomeDir"]+"/*.gtf")
    output:
        bam = "sm_star/{sample}_Aligned.out.bam",
    benchmark: "sm_benchmarks/star_alignment_{sample}.txt"
    params: 
        overhang = 100,
        genome_dir=directory(config["GenomeDir"]+"/index")
    threads: 8
    resources:
        mem_mb=32000
    shell:
        "STAR \
        --runThreadN {threads} \
        --outFileNamePrefix sm_star/{wildcards.sample}_ \
        --genomeDir {params.genome_dir} \
        --readFilesIn {input.R1} {input.R2} \
        --sjdbGTFfile {input.annotation} \
        --sjdbOverhang {params.overhang} \
        --outSAMtype BAM Unsorted"

checkpoint star_all:
    input: expand("sm_star/{sample}_Aligned.out.bam", sample=SAMPLES)
    output: "sm_star_done.txt"
    shell: "touch sm_star_done.txt" 
