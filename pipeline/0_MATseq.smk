#### Sanquin Bioinformatics ####
# name: MATseq pipeline
# description: Monocyte activation test transcriptomes analysis pipeline
# author: Tess Afanasyeva

import sys
import random
import csv
import datetime
from pathlib import Path
from glob import glob

sys.path.append(str(Path(__file__).parent.parent.parent.parent))
from src.config import get_work_dir, get_sample_dir, get_genome_dir

WORKDIR = get_work_dir()
SAMPLEDIR = get_sample_dir()
GENOMEDIR = get_genome_dir()

if not Path(SAMPLEDIR).exists():
    raise FileNotFoundError(
        f"sample_dir not found: {SAMPLEDIR}\n"
        "Provide via --fastq-dir or set snakemake.sample_dir in config.json"
    )

if not Path(GENOMEDIR).exists():
    raise FileNotFoundError(
        f"genome_dir not found: {GENOMEDIR}\n"
        "Provide via --genome-dir or set snakemake.genome_dir in config.json"
    )

print("WORKDIR:", WORKDIR)
print("SAMPLEDIR:", SAMPLEDIR)
print("GENOMEDIR:", GENOMEDIR)

def compile_benchmarks(benchmark: str, stats: str):
    """Aggregate all the benchmark files into one and put it in stats"""
    benchmark = Path(benchmark)
    if not benchmark.exists():
        benchmark.mkdir(parents=True, exist_ok=True)
    benchmarks = list(benchmark.iterdir())
    if not benchmarks:
        print("No benchmark files found")
        return None
    headers = ["rule"]
    with open(benchmarks[0], "r") as f:
        reader = csv.reader(f, delimiter="\t")
        headers += next(reader)

    stats = Path(stats)
    if not stats.exists():
        stats.mkdir(parents=True, exist_ok=True)
    stats_file = stats / f"{str(int(datetime.datetime.now().timestamp() * 1000))}_benchmarks.tsv"
    with open(stats_file, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)
        for fp in benchmarks:
            with open(fp, "r") as g:
                reader = csv.reader(g, delimiter="\t")
                next(reader)  # Headers line
                writer.writerow([fp.stem] + next(reader))
# Getting sample names and _R1 and _R2 strings
SAMPLES, _, _, PAIRED = map(set, glob_wildcards(str(Path(SAMPLEDIR) / "{samples}_{tag}_{other}_{paired}.fastq.gz")))
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
    compile_benchmarks(benchmark=str(Path(WORKDIR) / "sm_benchmarks"), stats=str(Path(WORKDIR) / "sm_stats"))
    print("Workflow finished successfully.")

ruleorder:
    star_alignment > samtools_merge

		
rule all:
	input:
		expand([WORKDIR/"sm_fastqs/{sample}_{paired}.fastq.gz",
				"sm_fastqc/{sample}_{paired}_fastqc.html",
				"sm_fastp/report/{sample}.html",
				"sm_star/{sample}_Aligned.out.bam",
				"sm_samtools_merged/{sample_id}.merged.bam",
				"sm_samtools/{sample_id}.sorted.bam",
				"sm_samtools/{sample_id}.sorted.bam.bai",
				"sm_umitools/{sample_id}.deduped.bam",
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
