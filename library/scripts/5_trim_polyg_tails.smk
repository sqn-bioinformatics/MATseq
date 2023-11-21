import os

class InvalidFileFormat(Exception):
    pass

# Sets path to ~/MATseq
path = ""

file_names = []
for file_name in os.listdir(os.path.join(path, "temp/unique_fastq")):
    if not (file_name.endswith("_R1.fastq") or (file_name.endswith("_R2.fastq"))):
        raise InvalidFileFormat(
            f'Invalid filename extension for {file_name}, requires 2 files of format "<file_name>_R1.fastq", "<file_name>_R2.fastq"'
        )
    elif file_name.endswith("_R1.fastq"):
        file_names.append(file_name.replace("_R1.fastq", ""))


rule all:
	input:
		expand("temp/trimmed_fastq/{sample}_R1.fastq", sample=file_names),
        expand("temp/trimmed_fastq/{sample}_R2.fastq", sample=file_names),

rule fastp_pe:
    input:
        sample=["temp/unique_fastq/{sample}_R1.fastq", "temp/unique_fastq/{sample}_R2.fastq"]
    output:
        trimmed=["temp/trimmed_fastq/{sample}_R1.fastq", "temp/trimmed_fastq/{sample}_R2.fastq"],
        html="/home/t.afanasyeva/MATseq/experiment/results/fastp/{sample}.html",
        json="/home/t.afanasyeva/MATseq/experiment/results/fastp/{sample}.json"
    log:
        "logs/fastp/{sample}.log"
    params:
        "--trim_poly_g --trim_poly_x"
    threads: 16
    wrapper:
        "v2.6.0/bio/fastp"
