import os
from parses_experiment_info import get_experiment_info
import datetime


# Sets path to ~/MATseq
path = ""

log_filename=datetime.datetime.now().strftime("%Y%m%d_%H%M%S")+".log"

exp_info, sample_dict = get_experiment_info()

project_number = exp_info["SQN_project_number"]

file_names=[]
for file_name in(os.listdir("temp/trimmed_fastq/")):
	if file_name.endswith("_R1.fastq"):
		file_names.append(file_name.replace("_R1.fastq",""))


# This tells snakemake what all the final files are to avoid wild card error

rule all:  # In future make the rule that removes everything form temp   
	input:
		expand(["temp/sam/{sample}.sam", 
		"temp/bed/{sample}.bed",
		"temp/read_counts/{sample}_genecount.csv",
		"logs/bowtie2/{sample}.log"], sample=SAMPLES),
		f"experiment/results/read_counts_analysis/{project_number}_raw_reads.csv", 
		f"experiment/results/read_counts_analysis/{project_number}_raw_reads.xlsx",
		f"experiment/results/stats/{project_number}_seqstats.xlsx",
		f"experiment/results/MLPC_prediction/{project_number}_MLPC_prediction.xlsx",
		f"experiment/results/MLPC_prediction/{project_number}_MLPC_prediction.csv",

rule Bowtie2_map:
	input:
		right="temp/trimmed_fastq/{sample}_R1.fastq",
		left="temp/trimmed_fastq/{sample}_R2.fastq"
	output:
		"temp/sam/{sample}.sam"
	log:
		"logs/bowtie2/{sample}.log"
	shell:
		"bowtie2 -x library/genomes/GRCh38_cds -1 {input.right} -2 {input.left} -p 16 --fr  --no-discordant --no-mixed -3 50 -5 10 -S {output} --no-unal 2>{log}"

rule Seq_stats: # This rule assumes that the demultiplexing was made prior to this pipeline
	input:
		deduped_stats="experiment/results/stats/deduping_stats.csv",
		seq_logs = expand("logs/bowtie2/{sample}.log", sample=SAMPLES),
	output:
		"experiment/results/stats/"+SQN_project_number+"_seqstats.xlsx"
	script:
		"parses_sequencing_stats.py"

rule sam_to_bed:
	input:
		"temp/sam/{sample}.sam"
	output:
		"temp/bed/{sample}.bed"
	script:
		"sam_to_bed.py"

rule count_genes_from_bed:
	input:
		"temp/bed/{sample}.bed"
	output:
		"temp/read_counts/{sample}_genecount.csv"
	script:
		"count_genes_from_bed.py"
		
		
rule merge_all_gene_counts:
    input:
      filenames=expand("temp/read_counts/{sample}_genecount.csv", sample=file_names),       
    output:
        "experiment/results/read_counts_analysis/"+SQN_project_number+"_raw_reads.xlsx",
        "experiment/results/read_counts_analysis/"+SQN_project_number+"_raw_reads.csv",
    script:
        "merge_all_gene_counts.py"
        

rule MLPC_classification:
    input:
      "experiment/results/read_counts_analysis/"+SQN_project_number+"_raw_reads.csv",
    output:
        "experiment/results/MLPC_prediction/"+SQN_project_number+"_MLPC_prediction.xlsx",
        "experiment/results/MLPC_prediction/"+SQN_project_number+"_MLPC_prediction.csv",
    script:
        "MLPC_classification.py"