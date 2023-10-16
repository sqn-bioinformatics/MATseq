import pandas as pd
from mylogger import get_logger


def extract_sample_name(log_path):
    return log_path.replace("logs/bowtie2/", "").replace(".log", "")


def calculate_mapped_reads(log_lines, logger):
    numbers = [
        int(line.split("(")[0].replace(" ", ""))
        for line in log_lines
        if "aligned concordantly exactly 1 time" in line
        or "aligned concordantly >1 times" in line
    ]

    mapped_reads = sum(numbers)

    return mapped_reads


def process_sample(sample, deduped_stats, logger):
    try:
        sample_stats = deduped_stats.loc[sample]
        unique_reads = sample_stats["unique"]
        initial_reads = sample_stats["initial"]

        log_path = f"logs/bowtie2/{sample}.log"
        with open(log_path, "r") as log_file:
            log_lines = log_file.readlines()
        mapped_reads = calculate_mapped_reads(log_lines, logger=logger)

        sample_statistics = {
            "sample": sample,
            "initial_reads": initial_reads,
            "unique_reads": unique_reads,
            "mapped_reads": mapped_reads,
        }

        return sample_statistics

    except KeyError:
        logger.error(f"file not found for {sample}")
        return None


def main():
    logger = get_logger(__name__)
    bowtie_log_paths = snakemake.input.seq_logs
    sample_list = [extract_sample_name(path) for path in bowtie_log_paths]

    deduped_stats = pd.read_csv(snakemake.input.deduped_stats, index_col="sample")

    stats = []
    for sample in sample_list:
        sample_statistics = process_sample(sample, deduped_stats, logger)
        if sample_statistics is None:
            logger.error(f"failed to generate statistics for {sample}")
        else:
            stats.append(sample_statistics)

    stats_df = pd.DataFrame(stats)
    stats_df.set_index("sample", inplace=True)

    stats_df.to_excel(snakemake.output[0])


if __name__ == "__main__":
    main()
