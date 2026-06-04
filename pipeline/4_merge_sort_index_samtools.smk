#### Sanquin Bioinformatics ####
# name: MATseq pipeline
# description: Monocyte activation test transcriptomes analysis pipeline
# author: Tess Afanasyeva


from glob import glob

def get_bam_files(wildcards):
    STAR_FILES = glob(str(Path(WORKDIR) / "sm_star/*_Aligned.out.bam"))
    outfiles= []
    for file in STAR_FILES:
        sample = Path(file).name.split("_Aligned")[0]
        if sample.split("_", 1)[1] != wildcards.sample_id:
            continue
        else:
            outfiles.append(file)
    return outfiles


rule samtools_merge:
    input:
        bams = get_bam_files,
        summary="sm_star_done.txt"
    output:
        combined_bam = "sm_samtools_merged/{sample_id}.merged.bam",
    benchmark: "sm_benchmarks/samtools_merge_{sample_id}.txt"
    threads: 8
    conda: "environment.yml"
    shell: "samtools merge -@ {threads} -o {output.combined_bam} {input.bams}"


rule samtools_sort:
# samtools - sort
# http://www.htslib.org/doc/samtools-sort.html
    input:
        aligned_bam = "sm_samtools_merged/{sample_id}.merged.bam",
    output:
        sorted_bam = "sm_samtools/{sample_id}.sorted.bam"
    message: "Rule {rule} sorts {sample_id}."
    benchmark: "sm_benchmarks/samtools_sort_{sample_id}.txt"
    threads: 8
    conda: "environment.yml"
    shell:
        " samtools sort \
        -@ {threads} \
        -T sm_samtools/ \
        -O bam \
        -o {output.sorted_bam} \
        --no-PG \
        {input.aligned_bam} "


rule samtools_index:
# samtools - index
# http://www.htslib.org/doc/samtools-index.html
    input:
        bam = "sm_samtools/{sample_id}.sorted.bam"
    output:
        "sm_samtools/{sample_id}.sorted.bam.bai"
    message:
        "Rule {rule} indexes {sample_id}."
    benchmark: "sm_benchmarks/samtools_index_{sample_id}.txt"
    threads: 8
    conda: "environment.yml"
    shell:
        "samtools index {input.bam}"
