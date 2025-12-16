"""Data preprocessing functions for RNA-seq data."""

from pathlib import Path
import numpy as np
import pandas as pd


def process_featurecounts_files(
    featurecounts_dir: str, output_file: str = None
) -> pd.DataFrame:
    """Process featurecounts files by trimming and merging them.

    This function:
    - Trims each featurecounts file to keep only Geneid (column 1) and counts (column 7)
    - Merges all files into a single DataFrame with genes as rows and samples as columns

    Args:
        featurecounts_dir: Path to directory containing featurecounts .txt files
        output_file: Optional path to save the merged output file

    Returns:
        DataFrame: Merged gene counts with Geneid as index and samples as columns
    """
    featurecounts_path = Path(featurecounts_dir)

    # Find all .txt files (excluding .summary files)
    txt_files = [
        f for f in featurecounts_path.glob("*.txt") if not f.name.endswith(".summary")
    ]

    if not txt_files:
        raise ValueError(f"No .txt files found in {featurecounts_dir}")

    merged_df = None

    for file_path in sorted(txt_files):
        # Read the file, skip comment line (line 0), use line 1 as header
        df = pd.read_csv(file_path, sep="\t", comment="#", skiprows=1)

        sample_name = file_path.stem

        # Keep only Geneid and the count column (last column, which is column 7 in original)
        trimmed_df = df[["Geneid", df.columns[-1]]].copy()
        trimmed_df.columns = ["Geneid", sample_name]

        if merged_df is None:
            merged_df = trimmed_df
        else:
            merged_df = merged_df.merge(trimmed_df, on="Geneid", how="outer")

    merged_df.set_index("Geneid", inplace=True)

    if output_file:
        merged_df.to_csv(output_file, sep="\t")
    else:
        # Save to parent directory (results/) with default name
        default_output = featurecounts_path.parent / "MATseq_count_summary.csv"
        merged_df.to_csv(default_output)

    return merged_df


def clean_counts(
    data: str | pd.DataFrame,
    mapping_path: str = None,
    min_reads: int = 1000000,
) -> pd.DataFrame:
    """Load RNA-seq counts and perform initial cleaning.

    Args:
        data: Either a file path (str) to counts file (tab-separated) or a pandas DataFrame.
        mapping_path: Path to the sample names mapping CSV. Optional if data is already a DataFrame.
        min_reads: Minimum total read count threshold for filtering samples.

    Returns:
        DataFrame: Cleaned counts data with samples as rows and genes as columns.
    """
    # Load data if file path is provided
    if isinstance(data, str):
        data = pd.read_csv(data, sep="\t", header=1, index_col="Geneid")
    elif isinstance(data, pd.DataFrame):
        # Make a copy to avoid modifying the original
        data = data.copy()
    else:
        raise TypeError("data must be either a file path (str) or a pandas DataFrame")

    # Transpose and handle duplicates
    data = data.T
    data = data[~data.index.duplicated(keep="last")]
    data = data[data.index.notnull()]
    data.sort_index(inplace=True)

    # Filter samples with insufficient read counts
    mask = data.sum(axis=1) > min_reads
    data = data[mask]

    data.index.name = "samples"
    data.reset_index(inplace=True)
    data.set_index("samples", inplace=True)

    return data


def prepare_training(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Extract training subset and labels from main dataset.

    Args:
        df: Main DataFrame with all samples.

    Returns:
        tuple: (training_data, labels) where labels is a DataFrame with 'label' column.
    """

    training_classes = [
        "_IMDM_",
        "_LPS_",
        "_Fla-PA_",
        "_PGN_",
        "_Pam3_",
        "_R848_",
    ]

    training_subset = [
        index
        for index, sample_id in enumerate(df.index)
        for c in training_classes
        if c in sample_id
    ]

    training_data = df.iloc[training_subset]

    labels = pd.DataFrame(index=training_data.index)
    labels["label"] = [i.split("_")[2] for i in training_data.index]

    return training_data, labels


def normalize_rpm(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize gene counts to reads per million (RPM).

    Args:
        df: DataFrame with gene counts (samples as rows, genes as columns).

    Returns:
        DataFrame: Normalized counts in RPM.
    """
    df = df.apply(
        lambda x: (x / (np.sum(x) if np.sum(x) != 0 else 1)) * 1000000, axis=1
    )
    return df
