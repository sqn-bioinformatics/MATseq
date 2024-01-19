# import os

# class InvalidFileFormat(Exception):
#     pass

# # Sets path to ~/MATseq
# path = ""

# file_names = []
# for file_name in os.listdir(os.path.join(path, "temp/unique_fastq")):
#     if not (file_name.endswith("_R1.fastq") or (file_name.endswith("_R2.fastq"))):
#         raise InvalidFileFormat(
#             f"Invalid filename extension for {file_name}, requires 2 files of format "<file_name>_R1.fastq", "<file_name>_R2.fastq""
#         )
#     elif file_name.endswith("_R1.fastq"):
#         file_names.append(file_name.replace("_R1.fastq", ""))


# rule all:
# 	input:
# 		expand("temp/trimmed_fastq/{sample}_R1.fastq", sample=file_names),
#         expand("temp/trimmed_fastq/{sample}_R2.fastq", sample=file_names),

# rule fastp_pe:
#     input:
#         sample=["temp/unique_fastq/{sample}_R1.fastq", "temp/unique_fastq/{sample}_R2.fastq"]
#     output:
#         trimmed=["temp/trimmed_fastq/{sample}_R1.fastq", "temp/trimmed_fastq/{sample}_R2.fastq"],
#         html="/home/t.afanasyeva/MATseq/experiment/results/fastp/{sample}.html",
#         json="/home/t.afanasyeva/MATseq/experiment/results/fastp/{sample}.json"
#     log:
#         "logs/fastp/{sample}.log"
#     params:
#         "--trim_poly_g --trim_poly_x"
#     threads: 16
#     wrapper:
#         "v2.6.0/bio/fastp"


#### Sanquin Bioinformatics ####
# description: RNA-seq pipeline
# author: Margo Schuller
# date: 2021-03-10


rule fastp:
# fastp - A tool designed to provide fast all-in-one preprocessing for FastQ files.
# https://github.com/OpenGene/fastp#quality-filter
    input:
        R1 = config["SampleDir"]+"/{sample}_R1.fastq.gz",
        R2 = config["SampleDir"]+"/{sample}_R2.fastq.gz"
    output:
        R1_tp = "fastp/reads/{sample}_R1_trimmed_paired.fastq.gz",
        R2_tp = "fastp/reads/{sample}_R2_trimmed_paired.fastq.gz",
        json = "fastp/report/{sample}.json",
        html = "fastp/report/{sample}.html",
    params:
        length_required = config["LengthRequired"],
        quality_score = config["QualityScore"]
    # conda:
    #     "environment.yml"
    threads: 2
    resources:
        mem_mb = 4000
    shell:
        "fastp \
         --in1 {input.R1} \
         --in2 {input.R2} \
         --out1 {output.R1_tp} \
         --out2 {output.R2_tp} \
         --html {output.html} \
         --thread {threads} \
         --dont_overwrite \
         --detect_adapter_for_pe \
         --trim_poly_g \
         --trim_poly_x \
         --length_required {params.length_required} \
     	 --low_complexity_filter \
         --disable_quality_filtering \
         --cut_right \
         --cut_right_window_size 4 \
         --cut_right_mean_quality {params.quality_score}"