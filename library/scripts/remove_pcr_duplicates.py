import os
import time
from mylogger import get_logger
import gzip
from Bio import SeqIO


# Gets logger instance
logger = get_logger(__name__)


def remove_duplicated_runs(file_path, unique_runs_per_sample):
    reads_dict = {}

    for file_name in unique_runs_per_sample:
        path = os.path.join(file_path, "temp/raw_fastq", file_name + "_R1.fastq")

        with gzip.open(path, mode="rt") as f:
            for record in SeqIO.parse(f, "fastq"):
                seq = record.seq
                umi = record.id.split("_")[-1]
                reads_dict[record.id] = umi + seq

        print(len(reads_dict))
        swap_dict = {v: k for k, v in reads_dict.items()}

        logger.info(
            "number of reads before deduping for "
            # + file_name[:-23] get sample id
            + " is "
            + str(len(reads_dict))
        )
        logger.info(
            "number of reads after deduping for "
            # + file_name[:-23]
            + " is "
            + str(len(swap_dict))
        )

    return swap_dict.keys()


def save_deduplicated(file_path, save_path, unique_runs_per_sample, unique_keys):
    for file_name in unique_runs_per_sample:
        path = os.path.join(file_path, "temp/raw_fastq", file_name + "_R1.fastq")

        with gzip.open(path, mode="rt") as f:
            for record in SeqIO.parse(f, "fastq"):
                if record.id in unique_keys:
                    SeqIO.write(record, save_path, "fastq")


def main():
    files_path = os.path.dirname(os.path.realpath(snakemake.input[0]))
    file_names = os.path.listdir(files_path)

    # Parses out the sample and run ids
    file_names_split = [filename.split("_") for filename in file_names]
    runs = {name[0] for name in file_names_split}
    sample_names = {name[1] for name in file_names_split}

    # Retrives the number of unique runs per sample
    for name in file_names_split:
        sample_id = name[1]

        logger.info("starting to dedup sample " + sample_id)

        unique_runs_per_sample = {
            "_".join(name[:-1]) for name in file_names_split if name[1] == sample_id
        }

        unique_runs_per_sample_length = len(unique_runs_per_sample)

        logger.info(
            f"found {unique_runs_per_sample_length} runs for sample " + sample_id
        )

        unique_keys = remove_duplicated_runs(
            files_path, unique_runs_per_sample, sample_id
        )
        save_path = "/home/t.afanasyeva/MATseq/temp/unique_fastq"
        save_deduplicated(files_path, save_path, unique_runs_per_sample, sample_id)


if __name__ == "__main__":
    main()
    # import cProfile

    # cProfile.run("main()")
