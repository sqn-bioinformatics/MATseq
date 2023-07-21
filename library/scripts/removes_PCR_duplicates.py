import os
import time
import datetime
import pandas as pd
from log_entry import log_entry
from functools import wraps
import re

path = ""


class FileNotFoundForRun(FileNotFoundError):
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


def change_header(fastq, run):
    for i in range(len(fastq)):
        fastq[i] = fastq[i].replace("\n", "_" + run + "\n")
    return fastq


def load_reads_if_multiple_runs(func):
    @wraps(func)
    def wrapper(sample, list_of_runs, is_merged):
        if is_merged:
            # If there's only one run, use the existing ddR1, ddR2, ddR3
            R1_reads, R2_reads, R3_reads = ddR1, ddR2, ddR3
        else:
            # Locates runs for sample
            R1_reads, R2_reads, R3_reads = [], [], []
            for run in list_of_runs:
                filenamelist = os.listdir(path + "temp/raw_fastq")
                pattern = rf"{run}_{sample}_[A-Z-]{{17}}_\D\d{{3}}_R\d[.]fastq"
                filename = next(
                    (
                        file_name[:-9]
                        for file_name in filenamelist
                        if run in file_name
                        and sample in file_name
                        and bool(re.match(pattern, file_name))
                    ),
                    None,
                )

                if filename is None:
                    raise FileNotFoundForRun(
                        f"File not found for run '{run}' and sample '{sample}'."
                    )
                R1_path = os.path.join(path, "temp/raw_fastq", filename + "_R1.fastq")
                R2_path = os.path.join(path, "temp/raw_fastq", filename + "_R2.fastq")
                R3_path = os.path.join(path, "temp/raw_fastq", filename + "_R3.fastq")

                # Load reads
                R1_reads.extend(load_reads(R1_path))
                log_entry(R1_path + " loaded", True, log_filename)
                R2_reads.extend(load_reads(R2_path))
                log_entry(R2_path + " loaded", True, log_filename)
                R3_reads.extend(load_reads(R3_path))
                log_entry(R3_path + " loaded", True, log_filename)

        return func(R1_reads, R2_reads, R3_reads, sample, run)

    return wrapper


@load_reads_if_multiple_runs
def remove_duplicates(R1_reads, R2_reads, R3_reads, sample, run):
    start_time = time.time()
    log_entry("remove_duplicates is called ", True, log_filename)
    log_entry("removing duplicates " + sample, True, log_filename)

    ninitial = len(R1_reads) // 4

    # Creates reads and UMI pairs
    mergedreads = [
        R1_reads[i * 4 + 1][25:50] + R3_reads[i * 4 + 1][25:50] + R2_reads[i * 4 + 1]
        for i in range(ninitial)
    ]
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
    log_entry("coordinate list created", True, log_filename)

    del mergedreads
    del uniquereads

    # Determines the position of each unique read in the initial fastq files
    dedupedR1 = [R1_reads[i] for i in coordinate_list]
    dedupedR2 = [R2_reads[i] for i in coordinate_list]
    dedupedR3 = [R3_reads[i] for i in coordinate_list]

    del R1_reads
    del R2_reads
    del R3_reads

    dedupedR1 = change_header(dedupedR1, run)
    dedupedR2 = change_header(dedupedR2, run)
    dedupedR3 = change_header(dedupedR3, run)

    log_entry(
        f"run deduped in {time.time() - start_time:.2f} seconds", True, log_filename
    )

    return dedupedR1, dedupedR2, dedupedR3, ninitial, nfinal


log_filename = datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + ".log"

# Loads all fastq filenames into a pandas dataframe
filenamelist = os.listdir(path + "temp/raw_fastq")
columns = ["run_ref", "sample_name"]
filename_df = pd.DataFrame(
    [dict(zip(columns, filename.split("_")[:2])) for filename in filenamelist]
)

# Creates a list with sample names
samplelist = filename_df["sample_name"].unique().tolist()

# Checks/creates saving folder
if "unique_fastq" not in os.listdir(os.path.join(path, "temp")):
    os.mkdir(os.path.join(path, "temp/unique_fastq"))


stats_columns = ["sample", "initial", "unique"]
stats_df = pd.DataFrame(columns=stats_columns)

for sample in samplelist:
    log_entry("starting to dedup sample " + sample, True, log_filename)
    # Creates a list with the different runs.
    # each run is composed of 3 files (R1 - forward read, R2 - UMI and R3 - reverse read)
    list_of_runs = filename_df[filename_df["sample_name"] == sample]["run_ref"].unique()
    log_entry(
        f"{len(list_of_runs)}" + " runs found for sample " + sample, True, log_filename
    )

    if len(list_of_runs) > 1 and ddR1:
        (
            dedupedR1_temp,
            dedupedR2_temp,
            dedupedR3_temp,
            ninitial,
            nfinal,
        ) = remove_duplicates(sample, list_of_runs, is_merged=True)
    else:
        (
            dedupedR1_temp,
            dedupedR2_temp,
            dedupedR3_temp,
            ninitial,
            nfinal,
        ) = remove_duplicates(sample, list_of_runs, is_merged=False)

    # Add deduplicated reads to current working fastq
    ddR1 = []
    ddR2 = []
    ddR3 = []

    ddR1.extend(dedupedR1_temp)
    ddR2.extend(dedupedR2_temp)
    ddR3.extend(dedupedR3_temp)

    del dedupedR1_temp, dedupedR2_temp, dedupedR3_temp

    # Adds stats to stats data frame
    stats_df = pd.concat(
        [
            stats_df,
            pd.DataFrame(
                {"sample": sample, "initial": ninitial, "unique": nfinal}, index=[0]
            ),
        ],
        ignore_index=True,
    )


log_entry("saving deduped fastq files", True, log_filename)

save_reads(os.path.join(path, "temp/unique_fastq", sample + "_R1.fastq"), ddR1)
save_reads(os.path.join(path, "temp/unique_fastq", sample + "_R2.fastq"), ddR2)
save_reads(os.path.join(path, "temp/unique_fastq", sample + "_R3.fastq"), ddR3)

del ddR1, ddR2, ddR3

log_entry("deduped fastq files saved", True, log_filename)
log_entry(sample + " deduped!", True, log_filename)

results_dir = "experiment/results"
stats_dir = os.path.join(results_dir, "stats")

# Create directories if they don't exist
os.makedirs(results_dir, exist_ok=True)
os.makedirs(stats_dir, exist_ok=True)

stats_df.to_csv(os.path.join(stats_dir, "deduping_stats.csv"))
