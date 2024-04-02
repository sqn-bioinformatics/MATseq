#### Sanquin Bioinformatics ####
# name: MATseq pipeline
# description: Monocyte activation test transcriptomes analysis pipeline
# author: Tess Afanasyeva
# date: 2024-01-01


from glob import glob

def required_files(wildcards):
    STAR_FILES = glob(config["WorkDir"]+"/star/*.bam")
    outfiles= []
    for file in STAR_FILES:
        if wildcards.sample_id not in file:
            continue
        else:
            outfiles.append(file)
    return outfiles


rule samtools_merge: 
    input:
        bams = required_files,
        summary="star_done.txt"
    output: 
        combined_bam = "samtools_merged/{sample_id}.merged.bam",
    benchmark: "benchmarks/samtools_merge_{sample_id}.txt"
    threads: 8
    shell: "samtools merge -@ {threads} -o {output.combined_bam} {input.bams}"


rule samtools_sort:
# samtools - sort
# http://www.htslib.org/doc/samtools-sort.html
    input:
        aligned_bam = "samtools_merged/{sample_id}.merged.bam",
    output:
        sorted_bam = "samtools/{sample_id}.sorted.bam"
    message: "Rule {rule} sorts {sample_id}."
    benchmark: "benchmarks/samtools_sort_{sample_id}.txt"
    threads: 8
    shell:
        " samtools sort \
        -@ {threads} \
        -T samtools/ \
        -O bam \
        -o {output.sorted_bam} \
        --no-PG \
        {input.aligned_bam} "


rule samtools_index:
# samtools - index
# http://www.htslib.org/doc/samtools-index.html
    input:
        bam = "samtools/{sample_id}.sorted.bam"
    output:
        "samtools/{sample_id}.sorted.bam.bai"
    message:
        "Rule {rule} indexes {sample_id}."
    benchmark: "benchmarks/samtools_index_{sample_id}.txt"
    threads: 8
    shell:
        "samtools index {input.bam}"
