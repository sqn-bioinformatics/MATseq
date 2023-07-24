import os
import time
import datetime
import pandas as pd
import gc
from log_entry import log_entry
from collections import OrderedDict

# from scratch.BA064.library.scripts.log_entry import log_entry #because Rstudio has path wrong
# path='scratch/BA064/' #because Rstudio has path wrong

path = ""  # because Rstudio has path wrong
log_filename = datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + ".log"


def removedups_with_loading(R1_path, R2_path, R3_path, sample, run):
    start_time = time.time()
    log_entry("starting removedups " + sample + " " + run, True, log_filename)

    # loads forward fastq
    R1_reads = open(R1_path).readlines()
    log_entry(R1_path + " loaded", True, log_filename)
    ninitial = int(len(R1_reads) / 4)

    R2_reads = open(R2_path).readlines()
    log_entry(R2_path + " loaded", True, log_filename)

    R3_reads = open(R3_path).readlines()
    log_entry(R3_path + " loaded", True, log_filename)

    # Creates reads and UMI pairs
    mergedreads = []
    for i in range(ninitial):
        mergedreads.append(
            R1_reads[i * 4 + 1][25:50]
            + R3_reads[i * 4 + 1][25:50]
            + R2_reads[i * 4 + 1]
        )
    log_entry("pairs created", True, log_filename)

    # Removes duplicates
    uniquereads = list(OrderedDict.fromkeys(mergedreads))

    nfinal = len(uniquereads)

    log_entry("reads dedduped", True, log_filename)
    log_entry(str(ninitial) + " total reads", True, log_filename)
    log_entry(str(nfinal) + " unique reads", True, log_filename)

    # Determines the position of each unique read in the initial fastq files
    rcount = 0
    coordinate_list = []
    for i in range(int(ninitial)):
        if rcount < len(uniquereads):
            if mergedreads[i] == uniquereads[rcount]:
                coordinate_list.append(i * 4)
                coordinate_list.append(i * 4 + 1)
                coordinate_list.append(i * 4 + 2)
                coordinate_list.append(i * 4 + 3)
                rcount = rcount + 1
    log_entry("coordenate list created", True, log_filename)
    del mergedreads
    del uniquereads

    dedupedR1 = []
    dedupedR2 = []
    dedupedR3 = []

    for i in coordinate_list:
        dedupedR1.append(R1_reads[i])
        dedupedR2.append(R2_reads[i])
        dedupedR3.append(R3_reads[i])

    del R1_reads
    del R2_reads
    del R3_reads

    # Changes header
    dedupedR1 = change_header(dedupedR1, run)
    dedupedR2 = change_header(dedupedR2, run)
    dedupedR3 = change_header(dedupedR3, run)

    log_entry(
        "run deduped in " + ("--- %s seconds ---" % (time.time() - start_time)),
        True,
        log_filename,
    )

    return (dedupedR1, dedupedR2, dedupedR3, ninitial, nfinal)


# Adds run info to each read header
def change_header(fastq, run):
    for i in range(int(len(fastq) / 4)):
        fastq[i * 4] = fastq[i * 4].replace("\n", "_" + run + "\n")
    return fastq


def removedups_NO_loading(ddR1, ddR2, ddR3, sample, run):
    start_time = time.time()
    log_entry("starting removedups " + sample + " " + run, True, log_filename)

    ninitial = int(len(ddR1) / 4)

    # Creates reads and UMI pairs
    mergedreads = []
    for i in range(ninitial):
        mergedreads.append(
            ddR1[i * 4 + 1][25:50] + ddR3[i * 4 + 1][25:50] + ddR2[i * 4 + 1]
        )
    log_entry("pairs created", True, log_filename)

    # Removes duplicates
    uniquereads = list(OrderedDict.fromkeys(mergedreads))

    nfinal = len(uniquereads)

    log_entry("reads dedduped", True, log_filename)
    log_entry(str(ninitial) + " total reads", True, log_filename)
    log_entry(str(nfinal) + " unique reads", True, log_filename)

    # Determines the position of each unique read in the initial fastq files
    rcount = 0
    coordinate_list = []

    for i in range(int(ninitial)):
        if rcount < len(uniquereads):
            if mergedreads[i] == uniquereads[rcount]:
                coordinate_list.append(i * 4)
                coordinate_list.append(i * 4 + 1)
                coordinate_list.append(i * 4 + 2)
                coordinate_list.append(i * 4 + 3)
                rcount = rcount + 1
    log_entry("coordenate list created", True, log_filename)

    del mergedreads
    del uniquereads

    dedupedR1 = []
    dedupedR2 = []
    dedupedR3 = []

    for i in coordinate_list:
        dedupedR1.append(ddR1[i])
        dedupedR2.append(ddR2[i])
        dedupedR3.append(ddR3[i])

    del ddR1
    del ddR2
    del ddR3

    log_entry(
        f"run deduped in {time.time() - start_time:.2f} seconds", True, log_filename
    )

    return (dedupedR1, dedupedR2, dedupedR3, ninitial, nfinal)


# Loads all fastq filenames into a pandas dataframe
filenamelist = os.listdir(path + "temp/raw_fastq")
columns = ["run_ref", "sample_name", "barcode", "L?", "seq_num"]
filename_DF = []

for filename in filenamelist:
    filename_DF.append(dict(zip(columns, filename.split("_"))))
filename_DF = pd.DataFrame(filename_DF)
filename_DF.drop("seq_num", axis=1, inplace=True)


# Creates a list with sample names
samplelist = list(dict.fromkeys(filename_DF.sample_name.values.tolist()))

"""Checks/creates saving folder

Important to save files, otherwise will raise error. 
"""
if "unique_fastq" not in os.listdir("temp"):
    os.mkdir(path + "temp/unique_fastq")

stats = []
stats_columns = ["sample", "run", "initial", "unique"]

for sample in samplelist:
    log_entry("starting to dedup sample " + sample, True, log_filename)

    ddR1 = []
    ddR2 = []
    ddR3 = []

    # Creates a dataframe with the multiple filenames that correspond to a specific sample
    tempDF = filename_DF[filename_DF["sample_name"] == sample]

    """Creates a list with the different runs.
    Each run is composed of 3 files (R1 - forward read, R2 - UMI and R3 - reverse read)
    """
    runlist = list(set(tempDF["run_ref"].values.tolist()))
    log_entry(
        str(len(runlist)) + " runs found for sample " + sample, True, log_filename
    )

    # Iterates every run to remove duplicates
    for run in runlist:
        runDF = tempDF[tempDF["run_ref"] == run]
        filename = ("_").join(runDF.iloc[0].values.tolist())

        # Determines file paths
        R1_path = path + "temp/raw_fastq/" + filename + "_R1.fastq"
        R2_path = path + "temp/raw_fastq/" + filename + "_R2.fastq"
        R3_path = path + "temp/raw_fastq/" + filename + "_R3.fastq"

        (ddR1_temp, ddR2_temp, ddR3_temp, ninitial, nfinal) = removedups_with_loading(
            R1_path, R2_path, R3_path, sample, run
        )
        gc.collect()

        # Adds dedduped reads to current working fastq
        ddR1.extend(ddR1_temp)
        ddR2.extend(ddR2_temp)
        ddR3.extend(ddR3_temp)
        log_entry("run added to final fastq files", True, log_filename)

        del ddR1_temp
        del ddR2_temp
        del ddR3_temp

        # Adds stats to statsDF
        stats.append(dict(zip(stats_columns, [sample, run, ninitial, nfinal])))

    # Removes duplicate from all merged runs
    if len(runlist) > 1:
        (ddR1, ddR2, ddR3, ninitial, nfinal) = removedups_NO_loading(
            ddR1, ddR2, ddR3, sample, "merged"
        )
        stats.append(dict(zip(stats_columns, [sample, "merged", ninitial, nfinal])))
        gc.collect()
    else:
        stats.append(dict(zip(stats_columns, [sample, "merged", ninitial, nfinal])))
        gc.collect()

    log_entry("saving deduped fastq files", True, log_filename)

    open("temp/unique_fastq/" + sample + "_R1.fastq", "w").writelines(ddR1)
    del ddR1

    open("temp/unique_fastq/" + sample + "_R2.fastq", "w").writelines(ddR2)
    del ddR2

    open("temp/unique_fastq/" + sample + "_R3.fastq", "w").writelines(ddR3)
    del ddR3

    log_entry("deduped fastq files saved", True, log_filename)
    gc.collect()

    log_entry(sample + " deduped!", True, log_filename)


stats_DF = pd.DataFrame(stats)

if "results" not in os.listdir("experiment"):
    os.mkdir("experiment/results")
if "stats" not in os.listdir("experiment/results"):
    os.mkdir("experiment/results/stats")

stats_DF.to_csv("experiment/results/stats/" + "deduping_stats.csv")
