import os
import re
from glob import glob
# from parses_experiment_info import get_experiment_info

# Path to files
FILES = glob(config["SampleDir"]+"/*.fastq.gz")

# Sample names
SAMPLES = sorted(list(set([ re.split("_R1|_R2", os.path.basename(x))[0] for x in FILES])))



# exp_info, _ = get_experiment_info()
# project_number = exp_info["SQN_project_number"]



# #### Target rules ####
rule all:
	input:
		expand(["fastqc/{sample}_{paired}_fastqc.html", "fastp/report/{sample}.html", "unique_fastq/{sample}_{paired}.fastq.gz"], sample=SAMPLES, paired=config["Paired"]),

# rule all:
# 	input:
# 		expand(["fastqc/{sample}_{paired}_fastqc.html", "fastp/report/{sample}.html", "unique_fastq/{sample}_{paired}.fastq.gz"], sample=SAMPLES, paired=config["Paired"]),

"""
Load rules 

"""
# include: "1_download_files.smk"
include: "2_control_quality.smk"
include: "3_trim_polyg_tails.smk"
include: "4_remove_pcr_duplicates.smk"



# hic sunt dracones
# 6_map_gene_counts.smk is being deconstructed into 4 files to have separete functions


# include: "5_trim_polyg_tails.smk"
# include: "6_aligment.smk"
# include: "7_gene_counts.smk"
# include: "8_classification.smk"

# include: "6_map_gene_counts.smk"


        # expand(["temp/sam/{sample}.sam", 
		# "temp/bed/{sample}.bed",
		# "temp/read_counts/{sample}_genecount.csv",
		# "logs/bowtie2/{sample}.log"], sample=SAMPLES
