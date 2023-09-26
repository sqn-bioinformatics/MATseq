import os
import time
import datetime
from functools import wraps
import re
import csv

# Sets path to ~/MATseq
path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class FileNotFoundForRun(FileNotFoundError):
    pass


class UMIMismatchError(ValueError):
    pass


# Creates logs
def log_entry(text, show, log_filename):
    if "logs" not in os.listdir("./"):
        os.mkdir("logs")

    text = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S: ") + text

    with open("logs/" + log_filename, "a") as file_object:
        file_object.write(text + "\n")
    if show:
        print(text)


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


# Decorator function for loading reads if there are multiple runs
def load_reads_if_multiple_runs(func):
    @wraps(func)
    def wrapper(sample, unique_runs_per_sample, has_many_runs):
        if has_many_runs:
            # If there's only one run, uses the existing ddR1, ddR2
            R1_reads, R2_reads = ddR1, ddR2
        else:
            # Locates runs for sample
            R1_reads, R2_reads = [], []
            for run in unique_runs_per_sample:
                pattern = rf"{run}_{sample}_[A-Z-]{{17}}_\D\d{{3}}_R\d[.]fastq"
                file_name = next(
                    (
                        name[:-9]
                        for name in file_names
                        if run in name
                        and sample in name
                        and bool(re.match(pattern, name))
                    ),
                    None,
                )

                if file_name is None:
                    raise FileNotFoundForRun(
                        f"File not found for run '{run}' and sample '{sample}'."
                    )
                R1_path = os.path.join(path, "temp/raw_fastq", file_name + "_R1.fastq")
                R2_path = os.path.join(path, "temp/raw_fastq", file_name + "_R2.fastq")

                # Load reads
                R1_reads.extend(load_reads(R1_path))
                log_entry(f"{file_name}_R1.fastq loaded", True, log_filename)
                R2_reads.extend(load_reads(R2_path))
                log_entry(f"{file_name}_R2.fastq loaded", True, log_filename)

            return func(R1_reads, R2_reads, sample)

    return wrapper


@load_reads_if_multiple_runs
def remove_duplicates(R1_reads, R2_reads, sample):
    start_time = time.time()
    log_entry("removing duplicates " + sample, True, log_filename)

    ninitial = len(R1_reads) // 4

    #
    generator = (
        R1_reads[i * 4].split()[0].split(":")[-1]
        if R1_reads[i * 4].split()[0].split(":")[-1]
        == R2_reads[i * 4].split()[0].split(":")[-1]
        else False
        for i in range(ninitial)
    )

    # Creates reads and UMI pairs, throws an error if there
    mergedreads = []
    for i, match in enumerate(generator):
        if not match:
            raise UMIMismatchError(f"UMI mismatch at index {i}.")
        else:
            mergedreads.append(
                match + R1_reads[i * 4 + 1][25:50] + R2_reads[i * 4 + 1][25:50]
            )

    log_entry("pairs created", True, log_filename)

    # Removes duplicates
    uniquereads = uniquereads_finder(mergedreads)
    nfinal = len(uniquereads)
    coordinate_list = []
    for i, mergedread in enumerate(mergedreads):
        if i < nfinal and mergedread == uniquereads[i]:
            coordinate_list.extend([i * 4, i * 4 + 1, i * 4 + 2, i * 4 + 3])

    log_entry("reads deduped", True, log_filename)
    log_entry(str(ninitial) + " total reads", True, log_filename)
    log_entry(str(nfinal) + " unique reads", True, log_filename)

    # Determines the position of each unique read in the initial fastq files
    dedupedR1 = [R1_reads[i] for i in coordinate_list]
    dedupedR2 = [R2_reads[i] for i in coordinate_list]

    log_entry(
        f"run deduped in {time.time() - start_time:.2f} seconds", True, log_filename
    )

    return dedupedR1, dedupedR2, ninitial, nfinal


def main():
    global log_filename
    log_filename = datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + ".log"

    # Checks/creates saving folder
    if "unique_fastq" not in os.listdir(os.path.join(path, "temp")):
        os.mkdir(os.path.join(path, "/temp/unique_fastq"))

    global file_names
    file_names = os.listdir(path + "/temp/raw_fastq")

    # Parses out the sample and run ids
    samples, runs = [], []

    for filename in file_names:
        run_ref, sample_name = filename.split("_")[:2]
        samples.append(sample_name)
        runs.append(run_ref)

    runs_samples_dict = dict(zip(runs, samples))

    # Generates list to store each sample statistics dictionary
    sample_stats_list = []

    # Retrives the number of unique runs per sample
    for sample_id in runs_samples_dict.values():
        log_entry("starting to dedup sample " + sample_id, True, log_filename)
        unique_runs_per_sample = {
            i for i in runs_samples_dict if runs_samples_dict[i] == sample_id
        }
        log_entry(
            f"found {len(unique_runs_per_sample)} runs for sample " + sample_id,
            True,
            log_filename,
        )

        global ddR1
        global ddR2
        ddR1 = []
        ddR2 = []

        if len(unique_runs_per_sample) > 1 and ddR1:
            (
                dedupedR1_temp,
                dedupedR2_temp,
                ninitial,
                nfinal,
            ) = remove_duplicates(sample_id, unique_runs_per_sample, has_many_runs=True)
        else:
            (
                dedupedR1_temp,
                dedupedR2_temp,
                ninitial,
                nfinal,
            ) = remove_duplicates(
                sample_id, unique_runs_per_sample, has_many_runs=False
            )

        ddR1.extend(dedupedR1_temp)
        ddR2.extend(dedupedR2_temp)

        log_entry("saving deduped fastq files", True, log_filename)

        save_reads(
            os.path.join(path, "temp/unique_fastq", sample_id + "_R1.fastq"), ddR1
        )
        save_reads(
            os.path.join(path, "temp/unique_fastq", sample_id + "_R2.fastq"), ddR2
        )

        log_entry("deduped fastq files saved", True, log_filename)
        log_entry(sample_id + " deduped!", True, log_filename)

        # Create a dictionary for the current sample's statistics
        sample_stats = {
            "sample": sample_id,
            "initial": ninitial,
            "unique": nfinal,
        }

        # Append the sample's statistics dictionary to the list
        sample_stats_list.append(sample_stats)

    # Generates the statistics directory and file path
    if "stats" not in os.listdir(os.path.join(path, "experiment/results")):
        os.mkdir(os.path.join(path, "experiment/results/stats/"))
    csv_file_path = f"experiment/results/stats/deduping_stats.csv"

    # Writes the list of sample statistics dictionaries to the CSV file
    with open(csv_file_path, "a", newline="") as f:
        fieldnames = ["sample", "initial", "unique"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)

        # Writes the CSV header if the file is empty
        if os.stat(os.path.join(path, "experiment/results/stats")).st_size == 0:
            writer.writeheader()

        for sample_stats in sample_stats_list:
            writer.writerow(sample_stats)

    print(f"Statistics saved to {os.path.join(path, 'experiment/results/stats')}")


if __name__ == "__main__":
    main()
