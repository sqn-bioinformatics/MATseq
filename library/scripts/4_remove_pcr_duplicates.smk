#### Sanquin Bioinformatics ####
# description: RNA-seq pipeline
# author: Bruno Albuquerque
# updated: Tess Afanasyeva
# date: 2023-11-30

# from removes_pcr_duplicates import main

rule remove_pcr_duplicates:
# In house PCR duplicate removal script
    input:
        read = config["WorkDir"] + "/fastp/reads/{sample}_{paired}_trimmed_paired.fastq.gz"
    output:
        deduped = "unique_fastq/{sample}_{paired}.fastq.gz"
    script:
        "removes_pcr_duplicates.py"


