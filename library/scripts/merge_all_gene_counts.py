import pandas as pd
from load_exp_info import load_exp_info

log_filename = snakemake.log[0]

exp_info, sample_dict = load_exp_info(log_filename)

filenamelist = snakemake.input

# loads all diferent gene counts files and creates a DF
reads_DF = pd.DataFrame()

for filename in filenamelist:
    in_DF = pd.read_csv(filename).set_index("gene")
    sample_name = filename.split("/")[2]
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
