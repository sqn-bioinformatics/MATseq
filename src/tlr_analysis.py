"""TLR reporter assay data loading and analysis."""

from pathlib import Path
import numpy as np
import pandas as pd


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
    tlr4_lps = tlr4_raw[tlr4_raw["OD630nm_LPS_Replicate1"].notna()]
    tlr4_df = pd.DataFrame(
        {
            "Concentration_EU_mL": tlr4_lps["Concentration_(EU_mL)"],
            "Average": tlr4_lps[
                ["OD630nm_LPS_Replicate1", "OD630nm_LPS_Replicate2"]
            ].mean(axis=1),
        }
    )

    tlr4_fla = tlr4_raw[tlr4_raw["OD630nm_Fla-PA_Replicate1"].notna()].iloc[0]

    tlr2_raw = pd.read_csv(data_dir / "Supplementary_Table_6.csv")
    tlr2_pam = tlr2_raw[tlr2_raw["OD630nm_Pam3_Replicate1"].notna()]
    tlr2_df = pd.DataFrame(
        {
            "Concentration_ng_mL": tlr2_pam["Concentration_(ng_mL)"],
            "Average": tlr2_pam[
                ["OD630nm_Pam3_Replicate1", "OD630nm_Pam3_Replicate2"]
            ].mean(axis=1),
        }
    )

    tlr2_fla = tlr2_raw[tlr2_raw["OD630nm_Fla-PA_Replicate1"].notna()].iloc[0]

    flapa_data = {
        "tlr4": {
            "concentration": tlr4_fla["Concentration_(EU_mL)"],
            "average": np.mean(
                [tlr4_fla["OD630nm_Fla-PA_Replicate1"], tlr4_fla["OD630nm_Fla-PA_Replicate2"]]
            ),
        },
        "tlr2": {
            "concentration": tlr2_fla["Concentration_(ng_mL)"],
            "average": np.mean(
                [tlr2_fla["OD630nm_Fla-PA_Replicate1"], tlr2_fla["OD630nm_Fla-PA_Replicate2"]]
            ),
        },
    }

    return tlr2_df, tlr4_df, flapa_data
