"""Data preprocessing functions for RNA-seq data."""

from pathlib import Path
import pandas as pd
import numpy as np
from .config import (
    TRAINING_LIGANDS,
    ADDITIONAL_LIGANDS,
    BACTERIAL_LIGANDS,
    TRAINING_LIGANDS_WO_FLAPA,
    BACTERIAL_LIGANDS_ORIGINAL_NAMES,
    get_featurecounts_dir,
)


def merge_counts(
    featurecounts_dir: str | None = None, output_path: str | None = None
) -> pd.DataFrame:
    """Process featurecounts files by trimming and merging them.

    This function:
    - Trims each featurecounts file to keep only Geneid (column 1) and counts (column 7)
    - Merges all files into a single DataFrame with genes as rows and samples as columns

    Args:
        featurecounts_dir: Path to directory containing featurecounts.txt files
        output_path: Optional path to save the merged output file

    Returns:
        DataFrame: Merged gene counts with Geneid as index and samples as columns
    """

    if featurecounts_dir is None:
        featurecounts_path = get_featurecounts_dir()
    else:
        featurecounts_path = Path(featurecounts_dir).expanduser()

    if not featurecounts_path.exists():
        raise FileNotFoundError(f"Directory not found: {featurecounts_path}")
    if not featurecounts_path.is_dir():
        raise NotADirectoryError(f"Not a directory: {featurecounts_path}")

    txt_files = [
        f for f in featurecounts_path.glob("*.txt") if not f.name.endswith(".summary")
    ]
    if not txt_files:
        raise ValueError(f"No .txt files found in {featurecounts_path}")

    dfs = []

    for f in txt_files:
        df = pd.read_csv(f, sep="\t", comment="#", skiprows=1, index_col=0)
        last_col = df.columns[-1]
        dfs.append(df.loc[:, [last_col]].rename(columns={last_col: f.stem}))

    counts_df = pd.concat(dfs, axis=1)
    counts_df = counts_df.rename_axis(index="samples", columns=None)
    filtered_counts_df = filter_counts(counts_df)

    if output_path is None:
        output_path = Path(__file__).parent.parent / "results" / "counts"
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    csv_path = output_path / "MATseq_count_summary.csv"
    filtered_counts_df.to_csv(csv_path)
    print(f"CSV saved: {csv_path}")

    return filtered_counts_df


def filter_counts(
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
    if not isinstance(df, pd.DataFrame):
        raise TypeError(f"df must be DataFrame, got {type(df)}")
    if len(df) == 0:
        raise ValueError("df cannot be empty")
    if min_reads <= 0:
        raise ValueError(f"min_reads must be positive, got {min_reads}")

    df = df.T
    df = df[~df.index.duplicated(keep="last")]
    df = df[df.index.notnull()]

    mask = df.sum(axis=1) > min_reads
    df = df[mask]
    if len(df) == 0:
        raise ValueError(f"No samples remain after filtering with min_reads={min_reads}")

    return df


def extract_subset(
    df: pd.DataFrame, name: str = "training", output_path: str | None = None
) -> pd.DataFrame:
    """Extract subset and labels from main dataset based on subset name.

    Args:
        df: Main DataFrame with all samples.
        name: Subset name - "training", "other_ligands", "bacterial", or "training_wo_flapa".
        output_path: Optional path to save the output file.

    Returns:
        DataFrame with samples as rows, genes as columns, and 'label' column.
    """

    # Select classes based on name
    if name == "training":
        selected_classes = [f"_{ligand}_" for ligand in TRAINING_LIGANDS + ["IMDM"]]
    elif name == "other_ligands":
        selected_classes = [f"_{ligand}_" for ligand in ADDITIONAL_LIGANDS + ["IMDM"]]
    elif name == "bacterial":
        selected_classes = [
            f"_{ligand}_" for ligand in BACTERIAL_LIGANDS_ORIGINAL_NAMES + ["IMDM"]
        ]
    elif name == "training_wo_flapa":
        selected_classes = [
            f"_{ligand}_" for ligand in TRAINING_LIGANDS_WO_FLAPA + ["IMDM"]
        ]
    else:
        raise ValueError(
            f"Invalid name '{name}'. Must be 'training', 'other_ligands', 'bacterial', or 'training_wo_flapa'."
        )

    subset_indices = [
        index
        for index, sample_id in enumerate(df.index)
        for c in selected_classes
        if c in sample_id
    ]
    if not subset_indices:
        raise ValueError(f"No samples found for subset '{name}'")

    subset_data = df.iloc[subset_indices].copy()

    labels = []
    for sample_id in subset_data.index:
        parts = sample_id.split("_")
        if len(parts) < 3:
            raise ValueError(
                f"Sample ID '{sample_id}' does not match expected format "
                "'<run>_<batch>_<ligand>_...'"
            )
        labels.append(parts[2])
    subset_data["label"] = labels

    # Replace label names
    label_mapping = {
        "HKSA": "HK S.aureus",
        "HKEB": "HK E.coli",
        "IMDM": "negative_control",
    }
    subset_data["label"] = subset_data["label"].replace(label_mapping)

    if output_path is None:
        output_path = Path(__file__).parent.parent / "results" / "subsets"
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    csv_path = output_path / f"{name}_data_with_labels.csv"
    subset_data.to_csv(csv_path)
    print(f"CSV saved: {csv_path}")

    return subset_data


def normalize_rpm(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize gene counts to reads per million (RPM).

    Args:
        df: DataFrame with gene counts (samples as rows, genes as columns).

    Returns:
        DataFrame: Normalized counts in RPM.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError(f"df must be DataFrame, got {type(df)}")
    if len(df) == 0:
        raise ValueError("df cannot be empty")

    row_sums = df.sum(axis=1).replace(0, 1)
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

    data_dir = Path(data_dir)
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")
    for fname in ("Supplementary_Table_5.csv", "Supplementary_Table_6.csv"):
        if not (data_dir / fname).exists():
            raise FileNotFoundError(f"Missing {fname} in {data_dir}")

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

    flapa_data = {
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

    return tlr2_df, tlr4_df, flapa_data
