"""Featurecounts loading, sample filtering, labelling, subsetting, RPM normalisation."""

from pathlib import Path

import numpy as np
import pandas as pd

from .config import (
    get_featurecounts_dir,
    MAIN_LIGANDS,
    ADDITIONAL_LIGANDS,
    BACTERIA_LIGANDS,
)

_BACTERIA_RAW_TO_DISPLAY = {"E.coli": "HK E.coli", "S.aureus": "HK S.aureus"}

_SUBSET_LIGANDS = {
    "main_ligands": MAIN_LIGANDS,
    "additional_ligands": ADDITIONAL_LIGANDS,
    "bacteria_ligands": BACTERIA_LIGANDS,
}


def load_featurecounts(featurecounts_dir: str | Path | None = None) -> pd.DataFrame:
    """Concatenate per-sample featureCounts .txt files into a (samples x genes) matrix."""
    data_dir = (
        Path(featurecounts_dir).expanduser() if featurecounts_dir else get_featurecounts_dir()
    )
    if not data_dir.is_dir():
        raise NotADirectoryError(data_dir)

    files = sorted(
        f for f in data_dir.glob("*.txt") if not f.name.endswith(".summary")
    )
    if not files:
        raise ValueError(f"No .txt files found in {data_dir}")

    per_sample = []
    for f in files:
        df = pd.read_csv(f, sep="\t", comment="#", skiprows=1, index_col=0)
        last_col = df.columns[-1]
        per_sample.append(df.loc[:, [last_col]].rename(columns={last_col: f.stem}))

    counts = pd.concat(per_sample, axis=1).T
    counts.index.name = "samples"
    counts = counts[~counts.index.duplicated(keep="last")]
    return counts


def filter_low_read_samples(
    counts: pd.DataFrame, min_reads: float = 1_000_000.0
) -> pd.DataFrame:
    """Drop samples whose total read count is below min_reads."""
    return counts[counts.sum(axis=1) > min_reads]


def label_samples(counts: pd.DataFrame) -> pd.Series:
    """Derive labels from sample names; the third underscore-separated token is the ligand."""
    labels = [
        _BACTERIA_RAW_TO_DISPLAY.get(name.split("_")[2], name.split("_")[2])
        for name in counts.index
    ]
    return pd.Series(labels, index=counts.index, name="label")


def prepare_counts(
    featurecounts_dir: str | Path | None = None,
    output_path: str | Path | None = None,
    min_reads: float = 1_000_000.0,
) -> tuple[pd.DataFrame, pd.Series]:
    """Load featureCounts, filter low-read samples, derive labels, persist summary."""
    counts = load_featurecounts(featurecounts_dir)
    counts = filter_low_read_samples(counts, min_reads=min_reads)
    labels = label_samples(counts)

    output_dir = (
        Path(output_path) if output_path else Path(__file__).parent.parent / "results" / "counts"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    counts.assign(label=labels).to_csv(output_dir / "MATseq_count_summary.csv")

    return counts, labels


def extract_subset(
    features: pd.DataFrame, labels: pd.Series, name: str
) -> tuple[pd.DataFrame, pd.Series]:
    """Filter (features, labels) to the samples whose label belongs to the named subset."""
    if name not in _SUBSET_LIGANDS:
        raise KeyError(
            f"Unknown subset '{name}'. Valid: {sorted(_SUBSET_LIGANDS)}"
        )
    ligands = _SUBSET_LIGANDS[name]
    mask = labels.isin(ligands)
    return features.loc[mask], labels.loc[mask]


def normalize_rpm(x: np.ndarray | pd.DataFrame) -> np.ndarray | pd.DataFrame:
    """Normalise per-sample counts to reads-per-million."""
    is_df = isinstance(x, pd.DataFrame)
    values = x.to_numpy(dtype=float) if is_df else np.asarray(x, dtype=float)

    if values.ndim != 2:
        raise ValueError(f"x must be 2D, got shape {values.shape}")
    if values.shape[0] == 0:
        raise ValueError("x cannot be empty")

    totals = values.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    rpm = (values / totals) * 1_000_000.0

    if is_df:
        return pd.DataFrame(rpm, index=x.index, columns=x.columns)
    return rpm
