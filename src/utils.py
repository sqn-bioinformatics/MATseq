"""Utility functions for data analysis and visualization."""

from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# Custom color palette
CUSTOM_PALETTE_6 = [
    "#1f77b4",  # Muted Blue
    "#ff7f0e",  # Soft Orange
    "#2ca02c",  # Green
    "#d62728",  # Red
    "#9467bd",  # Purple
    "#17becf",  # Teal/Cyan
]


def get_output_path() -> Path:
    """Get the output directory path based on current date.

    Returns:
        Path: Directory path for saving outputs (YYMMDD_output format).
    """
    p = Path().cwd().parent
    date = datetime.today().strftime("%Y%m%d")[2:]
    output_path = p / f"results/{date}"
    output_path.mkdir(parents=True, exist_ok=True)
    return output_path


def save_fig(
    fig_name: str,
    tight_layout: bool = True,
    fig_extension: str = "png",
    resolution: int = 300,
    output_path: Path = None,
) -> None:
    """Save figure to file.

    Args:
        fig_name: Name of the figure file (without extension).
        tight_layout: Whether to apply tight layout before saving.
        fig_extension: File extension/format for the figure.
        resolution: DPI resolution for the saved figure.
        output_path: Optional custom output path. If None, uses get_output_path().
    """
    if output_path is None:
        output_path = get_output_path()

    path = output_path / f"{fig_name}.{fig_extension}"
    if tight_layout:
        plt.tight_layout()
    plt.savefig(path, format=fig_extension, dpi=resolution)


def save_csv(table: pd.DataFrame, table_name: str, output_path: Path = None) -> None:
    """Save DataFrame to CSV file.

    Args:
        table: DataFrame to save.
        table_name: Name for the CSV file (without extension).
        output_path: Optional custom output path. If None, uses get_output_path().
    """
    if output_path is None:
        output_path = get_output_path()

    path = output_path / f"{table_name}.csv"
    table.to_csv(path)
