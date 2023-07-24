import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns


"""
Translation from R to python of DESeq2 normalisation. R code here:
https://scienceparkstudygroup.github.io/research-data-management-lesson/median_of_ratios_manual_normalization/index.html

"""


def normalization(data):
    import numpy as np
    import pandas as pd

    # Take the log
    log_data = np.log(data)

    # Calculate the pseudo-reference sample for each gene
    log_data["pseudo_reference"] = log_data.mean(axis=1)

    # Filter out genes with -Inf as their average
    filtered_log_data = log_data[log_data["pseudo_reference"] != float("-inf")]

    # Subtract the gene pseudo-references from log counts
    ratio_data = filtered_log_data.iloc[:, :-1].sub(
        filtered_log_data["pseudo_reference"], axis=0
    )

    # Find the median of the ratios for each sample
    sample_medians = ratio_data.median(axis=0)

    # Convert medians to scaling factors
    scaling_factors = np.exp(sample_medians)

    # Divide the original counts by the scaling factors
    manually_normalized = data.div(scaling_factors)

    # Calculate the mean count for each gene (row)
    mean_counts = manually_normalized.mean(axis=1)

    # Filter out rows where mean count is greater than 10
    filtered_normalized = manually_normalized[mean_counts >= 10]

    return filtered_normalized


current_script_path = os.path.realpath(__file__)
path = os.path.abspath(os.path.join(current_script_path, os.pardir))
filenamelist = os.listdir(path + "/input")

# Read the CSV file into a DataFrame
for file in filenamelist:
    file_path = os.path.join(path, "input", file)
    df_counts = pd.read_csv(file_path, index_col=0)
    normalized_counts = normalization(df_counts)
    print(normalized_counts.head())
