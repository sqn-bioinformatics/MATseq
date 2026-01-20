"""Configuration loader for MAT-seq pipeline."""

import json
from pathlib import Path
from functools import lru_cache
from typing import Any, Dict

_CONFIG_FILE = Path(__file__).parent / "config.json"
_config_cache: Dict[str, Any] = None


@lru_cache
def _load_config() -> dict[str, Any]:
    config_path = Path(__file__).parent.parent / "config.json"
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


def expand_path(path_str: str) -> Path:
    """Expand home directory and convert to Path.

    Args:
        path_str: Path string (may contain ~).

    Returns:
        Expanded Path object.
    """
    return Path(path_str).expanduser()


def get_data_dir() -> Path:
    """Get data directory path."""
    return expand_path(get_config("paths.data_dir"))


def get_deseq2_dir() -> Path:
    """Get DESeq2 data directory path."""
    return expand_path(get_config("paths.deseq2_dir"))


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
