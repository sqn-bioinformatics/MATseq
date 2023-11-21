import os
import time
import csv
from mylogger import get_logger
from itertools import chain


# Sets path to ~/MATseq
path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Gets logger instance
logger = get_logger(__name__)


class FileNotFoundForRun(FileNotFoundError):
    pass


class UMIMismatchError(ValueError):
    pass


# Load reads
def load_reads(filepath):
    with open(filepath) as file:
        return file.readlines()


# Save reads
def save_reads(filepath, reads):
    with open(filepath, "w") as file:
        file.writelines(reads)


# Find unique reads
def uniquereads_finder(lst):
    seen = set()
    seen_add = seen.add
    return [x for x in lst if not (x in seen or seen_add(x))]


def remove_duplicates_runs_loader(unique_runs_per_sample, has_many_runs):
    R1_reads, R2_reads = [], []
    deduped_R1_list, deduped_R2_list = [], []
    initial_read_number_total = 0

    for file_name in unique_runs_per_sample:
        if file_name is None:
            logger.error(f"file not found for {file_name}")
            raise FileNotFoundForRun(f"File not found for {file_name}")

        else:
            R1_path = os.path.join(path, "temp/raw_fastq", file_name + "_R1.fastq")
            R2_path = os.path.join(path, "temp/raw_fastq", file_name + "_R2.fastq")

            R1_reads.extend(load_reads(R1_path))
            logger.info(f"{file_name}_R1.fastq loaded")
            R2_reads.extend(load_reads(R2_path))
            logger.info(f"{file_name}_R2.fastq loaded")

            logger.info("removing duplicates " + file_name)

            deduped_R1, deduped_R2, initial_read_number = remove_duplicates(
                R1_reads, R2_reads
            )

            logger.info(
                "number of reads before deduping for "
                + file_name[:-23]
                + " is "
                + str(initial_read_number)
            )
            logger.info(
                "number of reads after deduping for "
                + file_name[:-23]
                + " is "
                + str(len(deduped_R1) // 4)
            )

            deduped_R1_list.append(deduped_R1)
            deduped_R2_list.append(deduped_R2)
            initial_read_number_total += initial_read_number

    if has_many_runs == False:
        return (
            deduped_R1_list,
            deduped_R2_list,
            initial_read_number_total,
            len(deduped_R1_list) // 4,
        )

    else:
        logger.info("removing duplicates from combined runs")
        deduped_R1_list, deduped_R2_list = list(
            chain.from_iterable(deduped_R1_list)
        ), list(chain.from_iterable(deduped_R2_list))

        deduped_R1_all_runs, deduped_R2_all_runs, _ = remove_duplicates(
            deduped_R1_list, deduped_R2_list
        )

        return (
            deduped_R1_all_runs,
            deduped_R2_all_runs,
            initial_read_number_total,
            len(deduped_R1_all_runs) // 4,
        )


def remove_duplicates(R1_reads, R2_reads):
    start_time = time.time()

    initial_read_number = len(R1_reads) // 4

    umi_R1_R2_concatenated = [
        R1_reads[i * 4].split()[0].split(":")[-1]
        + R1_reads[i * 4 + 1][25:50]
        + R2_reads[i * 4 + 1][25:50]
        if R1_reads[i * 4].split()[0].split(":")[-1]
        == R2_reads[i * 4].split()[0].split(":")[-1]
        else UMIMismatchError(f"UMI mismatch at index {i * 4}")
        for i in range(initial_read_number)
    ]

    logger.info("pairs created")

    unique_reads = uniquereads_finder(umi_R1_R2_concatenated)

    final_read_number = len(unique_reads)

    coordinate_list = []
    count = 0

    # Traverses through the indexes of original reads to look up reads equal to the unique one until the end of the unique reads list
    for i in range(initial_read_number):
        if count < final_read_number:
            if umi_R1_R2_concatenated[i] == unique_reads[count]:
                coordinate_list.append(i * 4)
                coordinate_list.append(i * 4 + 1)
                coordinate_list.append(i * 4 + 2)
                coordinate_list.append(i * 4 + 3)
                count += 1

    logger.info("reads deduped")
    logger.info(str(initial_read_number) + " total reads")
    logger.info(str(final_read_number) + " unique reads")

    # Determines the position of each unique read in the initial fastq files
    dedupedR1 = [R1_reads[i] for i in coordinate_list]
    dedupedR2 = [R2_reads[i] for i in coordinate_list]

    logger.info(f"run deduped in {time.time() - start_time:.2f} seconds")

    return dedupedR1, dedupedR2, initial_read_number


def main():
    # # Checks/creates saving folder
    # if "unique_fastq" not in os.listdir(os.path.join(path, "temp")):
    #     os.mkdir(os.path.join(path, "temp/unique_fastq"))

    file_names = os.listdir(os.path.join(path, "temp/raw_fastq"))

    # file_names = snakemake.input
    # print(file_names)

    # Parses out the sample and run ids
    runs_samples = [filename.split("_") for filename in file_names]

    # Generates list to store each sample statistics dictionary
    sample_stats_list, done_samples_list = [], []

    # # Retrives the number of unique runs per sample
    for run_sample in runs_samples[:1]:
        sample_id = run_sample[1]
        logger.info("starting to dedup sample " + sample_id)

        unique_runs_per_sample = {
            "_".join(run_sample[:-1])
            for run_sample in runs_samples
            if run_sample[1] == sample_id
        }
        unique_runs_per_sample_length = len(unique_runs_per_sample)
        logger.info(
            f"found {unique_runs_per_sample_length} runs for sample " + sample_id
        )

        (
            dedupedR1,
            dedupedR2,
            ninitial,
            nfinal,
        ) = remove_duplicates_runs_loader(
            unique_runs_per_sample,
            lambda has_many_runs: False if unique_runs_per_sample_length == 1 else True,
        )

        save_reads(
            os.path.join(path, "temp/unique_fastq", sample_id + "_R1.fastq"),
            dedupedR1,
        )
        save_reads(
            os.path.join(path, "temp/unique_fastq", sample_id + "_R2.fastq"),
            dedupedR2,
        )

        logger.info(sample_id + " deduped and saved!")

        # Create a dictionary for the current sample's statistics
        sample_stats = {
            "sample": sample_id,
            "initial": ninitial,
            "unique": nfinal,
        }

        # Append the sample's statistics dictionary to the list
        sample_stats_list.append(sample_stats)
        done_samples_list.append(sample_id)

    # # Generates the statistics directory and file path
    # if "stats" not in os.listdir(os.path.join(path, "experiment/results/")):
    #     os.mkdir(os.path.join(path, "experiment/results/stats/"))

    csv_file_path = os.path.join(path, "experiment/results/stats/deduping_stats.csv")
    sample_stats = dict(zip(done_samples_list, sample_stats_list))

    # Writes the list of sample statistics dictionaries to the CSV file
    with open(csv_file_path, "w", newline="") as f:
        fieldnames = ["sample", "initial", "unique"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for sample_name, stats in sample_stats.items():
            writer.writerow(
                {
                    "sample": sample_name,
                    "initial": stats["initial"],
                    "unique": stats["unique"],
                }
            )

    logger.info(f"Statistics saved to experiment/results/stats/")


if __name__ == "__main__":
    import cProfile

    cProfile.run("main()")
