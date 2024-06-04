#### Sanquin Bioinformatics ####
# name: MATseq pipeline
# description: Monocyte activation test transcriptomes analysis pipeline
# author: Tess Afanasyeva
# date: 2024-01-01

import os
import random
import csv
import datetime
from glob import glob

def compile_benchmarks(benchmark: str, stats: str):
    """Aggregate all the benchmark files into one and put it in stats"""
    if not os.path.exists(benchmark):
        os.makedirs(benchmark)
    benchmarks = os.listdir(benchmark)
    if not benchmarks:
        print("No benchmark files found")
        return None
    headers = ["rule"]
    with open(os.path.join(benchmark, benchmarks[0]), "r") as f:
        reader = csv.reader(f, delimiter="\t")
        headers += next(reader)

    if not os.path.exists(stats):
        os.makedirs(stats)
    stats_file = os.path.join(
        stats,
        f"{str(int(datetime.datetime.now().timestamp() * 1000))}_benchmarks.tsv",
    )
    with open(stats_file, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)
        for fp in benchmarks:	
            with open(os.path.join(benchmark, fp), "r") as g:
                reader = csv.reader(g, delimiter="\t")
                next(reader)  # Headers line
                writer.writerow([fp[:-4]] + next(reader))

# Getting sample names and _R1 and _R2 strings
SAMPLES, _, _, PAIRED = map(set, glob_wildcards(os.path.join(config['SampleDir'], "{samples}_{tag}_{other}_{paired}.fastq.gz")))
SAMPLE_RUN_DICTIONARY = {}
for sample in SAMPLES:
    sample_id = sample.split("_")[1]
    run = sample.split("_")[0]
    SAMPLE_RUN_DICTIONARY.setdefault(sample_id, []).append(run)

SAMPLE_IDS = SAMPLE_RUN_DICTIONARY.keys()

"""
Target rules

"""

onsuccess:
    compile_benchmarks(benchmark=os.path.join(config['WorkDir'], "benchmarks"), stats=os.path.join(config['WorkDir'], "stats"))
    print("Workflow finished successfully.")

ruleorder:
    star_alignment > samtools_merge

		
rule all:
	input: 
		expand([config["WorkDir"]+"/fastqs/{sample}_{paired}.fastq.gz",
				"fastqc/{sample}_{paired}_fastqc.html",
				"fastp/report/{sample}.html",
				"star/{sample}_Aligned.out.bam",
				"samtools_merged/{sample_id}.merged.bam",
				"samtools/{sample_id}.sorted.bam",
				"samtools/{sample_id}.sorted.bam.bai",
				"umitools/{sample_id}.deduped.bam",
				"featurecounts/{sample_id}.txt",
				], sample=SAMPLES, sample_id=SAMPLE_IDS, paired=PAIRED),


"""
Load rules 

"""
# # include: "1_download_files.smk"
include: "1_control_quality_fastqc.smk"
include: "2_trim_fastp.smk"
include: "3_align_star.smk"
include: "4_merge_sort_index_samtools.smk"
include: "5_deduplicate_umitools.smk"
include: "6_count_reads_featurecounts.smk"
