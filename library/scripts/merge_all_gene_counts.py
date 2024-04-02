import pandas as pd
from parses_experiment_info import get_experiment_info
from mylogger import get_logger

# Gets logger instance
logger = get_logger(__name__)

_, sample_dict = get_experiment_info()
file_names = snakemake.input

logger.info("loaded bed files")
# loads all diferent gene counts files and creates a DF
reads_DF = pd.DataFrame()

for file_name in file_names:
    in_DF = pd.read_csv(file_name).set_index("gene")
    sample_name = file_name.split("/")[2]
    sample_name = sample_dict[sample_name.replace("_genecount.csv", "")]
    reads_DF[sample_name] = in_DF["raw_reads"]

reads_DF.columns.name = "sample"

# Creates rawreadsDF with ordered columns
sorted_sample_names = sorted(reads_DF.columns, key=str.lower)
output_DF = pd.DataFrame()

for sample_name in sorted_sample_names:
    output_DF[sample_name] = reads_DF[sample_name]

# Saves rawreadsDF as excel and csv
output_DF.to_excel(snakemake.output[0])
output_DF.to_csv(snakemake.output[1])
