#### Sanquin Bioinformatics ####
# name: MATseq pipeline
# description: Monocyte activation test transcriptomes analysis pipeline
# author: Tess Afanasyeva


rule umitools_deduplication:
    input:
        bam = rules.samtools_sort.output,
        bai = rules.samtools_index.output
    output:
        deduplicated = "sm_umitools/{sample_id}.deduped.bam"
    message: "Rule {rule} removes duplicates in {sample_id}."
    benchmark: "sm_benchmarks/umitools_{sample_id}.txt"
    conda: "environment.yml"
    shell: "umi_tools dedup \
    --paired \
    --method=directional \
    --output-stats= sm_umitools/{wildcards.sample_id} \
    --umi-separator=':' \
    --chimeric-pairs=discard \
    --unpaired-reads=discard \
    -I {input.bam} \
    -S {output.deduplicated} "
