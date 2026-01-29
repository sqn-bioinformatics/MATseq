"""Configuration loader for MAT-seq pipeline."""

import json
from pathlib import Path
from functools import lru_cache
from typing import Any, Dict


def _find_config_path() -> Path:
    """Find config.json in project root or current working directory."""
    candidates = [
        Path(__file__).parent.parent / "config.json",
        Path.cwd() / "config.json",
    ]
    for path in candidates:
        if path.exists():
            return path
    return candidates[0]


@lru_cache
def _load_config() -> dict[str, Any]:
    config_path = _find_config_path()
    try:
        with config_path.open() as f:
            return json.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"Config file not found: {config_path}")
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in config file: {e}")


def get_config(key: str) -> Any:
    """Get a configuration value using dot notation.

    Example:
        get_config("paths.data_dir")
    """
    config = _load_config()
    for part in key.split("."):
        try:
            config = config[part]
        except (TypeError, KeyError):
            raise KeyError(f"Config key '{key}' not found")
    return config


def _get_project_root() -> Path:
    """Get project root directory."""
    src_dir = Path(__file__).parent
    if src_dir.name == "src":
        return src_dir.parent
    return Path.cwd()


def expand_path(path_str: str) -> Path:
    """Expand path relative to project root or home directory.

    Args:
        path_str: Path string (may contain ~ or relative path).

    Returns:
        Expanded absolute Path object.
    """
    path = Path(path_str)
    if path_str.startswith("~"):
        return path.expanduser()
    if path_str.startswith("~/MATseq"):
        relative = path_str.replace("~/MATseq/", "")
        return _get_project_root() / relative
    if not path.is_absolute():
        return _get_project_root() / path
    return path


def get_work_dir() -> Path:
    """Get Snakemake working directory."""
    return expand_path(get_config("snakemake.work_dir"))


def get_sample_dir() -> Path:
    """Get directory containing input samples."""
    return expand_path(get_config("snakemake.sample_dir"))


def get_genome_dir() -> Path:
    """Get genome reference directory."""
    return expand_path(get_config("snakemake.genome_dir"))


def get_data_dir() -> Path:
    """Get data directory path."""
    return expand_path(get_config("paths.data_dir"))


def get_go_terms_support_dir() -> Path:
    """Get GO terms support data directory path."""
    return expand_path(get_config("paths.go_terms_support_dir"))


def get_results_dir() -> Path:
    """Get results directory path."""
    return expand_path(get_config("paths.results_dir"))


def get_cache_dir() -> Path:
    """Get cache directory path."""
    return expand_path(get_config("paths.cache_dir"))


def get_figures_dir() -> Path:
    """Get figures directory path."""
    return expand_path(get_config("paths.figures_dir"))


def get_featurecounts_dir() -> Path:
    """Get featurecounts directory path."""
    return expand_path(get_config("paths.featurecounts_dir"))


# Load and export configuration as module-level constants for backward compatibility
_config = _load_config()

CUSTOM_PALETTE_6 = _config["colors"]["palette_6"]
CUSTOM_PALETTE_8 = _config["colors"]["palette_8"]
CUSTOM_PALETTE_9 = _config["colors"]["palette_9"]

CLASS_ORDER_TRAINING = _config["class_orders"]["training"]
CLASS_ORDER_OTHER_LIGANDS = _config["class_orders"]["other_ligands"]
CLASS_ORDER_BACTERIAL = _config["class_orders"]["bacterial"]

SUBSET_PALETTES = {
    "training": CUSTOM_PALETTE_6,
    "other_ligands": CUSTOM_PALETTE_9,
    "bacterial": CUSTOM_PALETTE_8,
}

SUBSET_CLASS_ORDERS = {
    "training": CLASS_ORDER_TRAINING,
    "other_ligands": CLASS_ORDER_OTHER_LIGANDS,
    "bacterial": CLASS_ORDER_BACTERIAL,
}

DESEQ2_CONFIG = _config["deseq2"]
FEATURE_SELECTION_CONFIG = _config["feature_selection"]
MODEL_FACTORY_CONFIG = _config["model_factory"]
MODEL_TRAINING_CONFIG = _config["model_training"]
CLASS_LABELS = _config["class_labels"]

TRAINING_LIGANDS = _config["ligands"]["training"]
ADDITIONAL_LIGANDS = _config["ligands"]["additional"]
BACTERIAL_LIGANDS_ORIGINAL_NAMES = _config["ligands"]["bacterial"]
BACTERIAL_LIGANDS = _config["ligands"]["bacterial_renamed"]
TRAINING_LIGANDS_WO_FLAPA = _config["ligands"]["training_wo_flapa"]
