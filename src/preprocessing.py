"""Data preprocessing functions for RNA-seq data."""

from pathlib import Path
import pandas as pd
import numpy as np


def process_featurecounts_files(
    featurecounts_dir: str, output_path: str | None = None
) -> pd.DataFrame:
    """Process featurecounts files by trimming and merging them.

    This function:
    - Trims each featurecounts file to keep only Geneid (column 1) and counts (column 7)
    - Merges all files into a single DataFrame with genes as rows and samples as columns

    Args:
        featurecounts_dir: Path to directory containing featurecounts .txt files
        output_path: Optional path to save the merged output file

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

    dfs = []

    for f in txt_files:
        df = pd.read_csv(f, sep="\t", comment="#", skiprows=1, index_col=0)
        last_col = df.columns[-1]
        dfs.append(df.loc[:, [last_col]].rename(columns={last_col: f.stem}))

    counts_df = pd.concat(dfs, axis=1)
    filtered_counts_df = filtered_counts(counts_df)
    filtered_counts_df.index.rename("samples")

    # Set default output path
    if output_path is None:
        output_path = featurecounts_path.parent / "results" / "counts"
    output_path.mkdir(parents=True, exist_ok=True)
    filtered_counts_df.to_csv(output_path / "MATseq_count_summary.csv")

    return filtered_counts_df


def filtered_counts(
    df: pd.DataFrame,
    min_reads: int = 1000000,
) -> pd.DataFrame:
    """Load RNA-seq counts and perform initial filtering.

    Args:
        df: pandas DataFrame with raw counts.
        min_reads: Minimum total read count threshold for filtering samples.

    Returns:
        DataFrame: Filtered counts data with samples as rows and genes as columns.
    """

    # Transpose and handle duplicates
    df = df.rename_axis(index="samples", columns=None).T
    df = df[~df.index.duplicated(keep="last")]
    df = df[df.index.notnull()]

    # Filter samples with insufficient read counts
    mask = df.sum(axis=1) > min_reads
    df = df[mask]

    return df


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
    row_sums = df.sum(axis=1).replace(0, 1)  # Avoid division by zero
    return df.div(row_sums, axis=0) * 1000000


def load_tlr_data(data_dir: Path = None) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Load TLR2 (Pam3) and TLR4 (LPS) data from supplementary tables.

    Args:
        data_dir: Path to the supplementary_data directory.
                  Defaults to data/supplementary_data relative to project root.

    Returns:
        Tuple of (tlr2_df, tlr4_df, fla_pa_data) where fla_pa_data contains
        Fla-PA measurements for each TLR.
    """
    if data_dir is None:
        data_dir = Path(__file__).parent.parent / "data" / "supplementary_data"

    # Load TLR4/LPS data (Supplementary Table 5)
    tlr4_raw = pd.read_csv(data_dir / "Supplementary_Table_5.csv")
    tlr4_lps_mask = tlr4_raw["OD630nm_LPS_Replicate1"].notna()
    tlr4_lps = tlr4_raw[tlr4_lps_mask]
    tlr4_df = pd.DataFrame(
        {
            "Concentration_EU_mL": tlr4_lps["Concentration_(EU_mL)"],
            "Average": tlr4_lps[
                ["OD630nm_LPS_Replicate1", "OD630nm_LPS_Replicate2"]
            ].mean(axis=1),
        }
    )

    # Extract Fla-PA data for TLR4 (row with Fla-PA values)
    tlr4_fla_mask = tlr4_raw["OD630nm_Fla-PA_Replicate1"].notna()
    tlr4_fla = tlr4_raw[tlr4_fla_mask].iloc[0]

    # Load TLR2/Pam3Csk4 data (Supplementary Table 6)
    tlr2_raw = pd.read_csv(data_dir / "Supplementary_Table_6.csv")
    tlr2_pam_mask = tlr2_raw["OD630nm_Pam3_Replicate1"].notna()
    tlr2_pam = tlr2_raw[tlr2_pam_mask]
    tlr2_df = pd.DataFrame(
        {
            "Concentration_ng_mL": tlr2_pam["Concentration_(ng_mL)"],
            "Average": tlr2_pam[
                ["OD630nm_Pam3_Replicate1", "OD630nm_Pam3_Replicate2"]
            ].mean(axis=1),
        }
    )

    # Extract Fla-PA data for TLR2 (row with Fla-PA values)
    tlr2_fla_mask = tlr2_raw["OD630nm_Fla-PA_Replicate1"].notna()
    tlr2_fla = tlr2_raw[tlr2_fla_mask].iloc[0]

    fla_pa_data = {
        "tlr4": {
            "concentration": tlr4_fla["Concentration_(EU_mL)"],
            "average": np.mean(
                [
                    tlr4_fla["OD630nm_Fla-PA_Replicate1"],
                    tlr4_fla["OD630nm_Fla-PA_Replicate2"],
                ]
            ),
        },
        "tlr2": {
            "concentration": tlr2_fla["Concentration_(ng_mL)"],
            "average": np.mean(
                [
                    tlr2_fla["OD630nm_Fla-PA_Replicate1"],
                    tlr2_fla["OD630nm_Fla-PA_Replicate2"],
                ]
            ),
        },
    }

    return tlr2_df, tlr4_df, fla_pa_data
