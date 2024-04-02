#### Sanquin Bioinformatics ####
# name: MATseq pipeline
# description: Monocyte activation test transcriptomes analysis pipeline
# author: Tess Afanasyeva
# date: 2024-01-01


rule featurecounts:
# https://subread.sourceforge.net/SubreadUsersGuide.pdf
    input:
        bam="umitools/{sample_id}.deduped.bam",
        annotation = glob(config["GenomeDir"]+"/*.gtf"),
    output:
        counts = "featurecounts/{sample_id}.txt",
        statistics = "featurecounts/{sample_id}.txt.summary"
    log: "logs/featurecounts/{sample_id}.log"
    message: "Rule {rule} counts gene reads in {sample_id}."
    benchmark: "benchmarks/featurecounts_{sample_id}.txt"
    threads: 42
    shell: "featureCounts -T {threads} -p --countReadPairs -C -B -t exon -g gene_id -a {input.annotation} -o {output.counts} {input.bam}"

# # Trim columns to leave those containing Geneid and sample_id name 
# rule trim:
#     input:
#         rules.featurecounts.output.counts
#     output:
#         "featurecounts/{sample_id}_clean.txt"
#     log: "logs/featurecounts_trim/{sample_id}.log"
#     benchmark: "benchmarks/featurecounts_trim/{sample_id}.txt"
#     run: 
#         shell("cut -f1,7 {input} > {output}")
#         shell("cut -f1 {input} > genes.txt")
#         shell("genes.txt > output.txt")

# # Merge genecounts into one dataset
# rule merge:
#     input:
#         rules.trim.output
#     output:
#         temp="{gene_id}.txt",
#     run:
#         shell("cp {input} {output.temp}"),
#         shell("paste {input} > genereads_summary/output.txt")








