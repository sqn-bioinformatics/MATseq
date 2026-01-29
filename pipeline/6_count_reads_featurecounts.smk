#### Sanquin Bioinformatics ####
# name: MATseq pipeline
# description: Monocyte activation test transcriptomes analysis pipeline
# author: Tess Afanasyeva


rule featurecounts:
# https://subread.sourceforge.net/SubreadUsersGuide.pdf
    input:
        bam="sm_umitools/{sample_id}.deduped.bam",
        annotation = glob(config["GenomeDir"]+"/*.gtf"),
    output:
        counts = "featurecounts/{sample_id}.txt",
        statistics = "featurecounts/{sample_id}.txt.summary"
    log: "sm_logs/featurecounts/{sample_id}.log"
    message: "Rule {rule} counts gene reads in {sample_id}."
    benchmark: "sm_benchmarks/featurecounts_{sample_id}.txt"
    threads: 42
    conda: "environment.yml"
    shell: "featureCounts -T {threads} -p --countReadPairs -C -B -t exon -g gene_id -a {input.annotation} -o {output.counts} {input.bam}"







