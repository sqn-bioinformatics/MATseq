import os

class InvalidFileFormat(Exception):
    pass

# Sets path to ~/MATseq
path = ""

file_names = []
for file_name in os.listdir(os.path.join(path, "temp/raw_fastq")):
    if not (file_name.endswith("_R1.fastq") or (file_name.endswith("_R2.fastq"))):
        raise InvalidFileFormat(
            f'Invalid filename extension for {file_name}, requires 2 files of format "<file_name>_R1.fastq", "<file_name>_R2.fastq"'
        )
    elif file_name.endswith("_R1.fastq"):
        file_names.append(file_name.replace("_R1.fastq", ""))

rule all:
    input:
        expand("temp/unique_fastq/{sample}.fastq", sample=file_names)

rule PCR_duplicate_remove:
    input:
        "temp/raw_fastq/{sample}.fastq"
    output:
        "temp/unique_fastq/{sample}.fastq",
        "experiment/results/stats/deduping_stats.csv"
    script:
        "removes_PCR_duplicates.py"