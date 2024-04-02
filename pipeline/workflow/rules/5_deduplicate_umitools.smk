#### Sanquin Bioinformatics ####
# name: MATseq pipeline
# description: Monocyte activation test transcriptomes analysis pipeline
# author: Tess Afanasyeva
# date: 2024-01-01


rule umitools_deduplication:
    input:
        bam = rules.samtools_sort.output,
        bai = rules.samtools_index.output
    output:
        deduplicated = "umitools/{sample_id}.deduped.bam"
    message: "Rule {rule} removes duplicates in {sample_id}."
    benchmark: "benchmarks/umitools_{sample_id}.txt"
    shell: "umi_tools dedup \
    --paired \
    --method=directional \
    --output-stats= umitools/{wildcards.sample_id} \
    --umi-separator=':' \
    --chimeric-pairs=discard \
    --unpaired-reads=discard \
    -I {input.bam} \
    -S {output.deduplicated} "
