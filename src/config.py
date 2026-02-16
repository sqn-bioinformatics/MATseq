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
    if path_str.startswith("~/MATseq"):
        relative = path_str.replace("~/MATseq/", "")
        return _get_project_root() / relative
    if path_str.startswith("~"):
        return path.expanduser()
    if not path.is_absolute():
        return _get_project_root() / path
    return path


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

CLASS_ORDER_TRAINING_WO_FLAPA = _config["class_orders"]["training_wo_flapa"]

SUBSET_PALETTES = {
    "training": CUSTOM_PALETTE_6,
    "other_ligands": CUSTOM_PALETTE_9,
    "bacterial": CUSTOM_PALETTE_8,
    "training_wo_flapa": CUSTOM_PALETTE_6[:len(CLASS_ORDER_TRAINING_WO_FLAPA)],
}

SUBSET_CLASS_ORDERS = {
    "training": CLASS_ORDER_TRAINING,
    "other_ligands": CLASS_ORDER_OTHER_LIGANDS,
    "bacterial": CLASS_ORDER_BACTERIAL,
    "training_wo_flapa": CLASS_ORDER_TRAINING_WO_FLAPA,
}

DESEQ2_CONFIG = _config["deseq2"]
FEATURE_SELECTION_CONFIG = _config["feature_selection"]
MODEL_FACTORY_CONFIG = _config["model_factory"]
MODEL_TRAINING_CONFIG = _config["model_training"]
TRAINING_LIGANDS = _config["ligands"]["training"]
ADDITIONAL_LIGANDS = _config["ligands"]["additional"]
BACTERIAL_LIGANDS_ORIGINAL_NAMES = _config["ligands"]["bacterial"]
BACTERIAL_LIGANDS = _config["ligands"]["bacterial_renamed"]
TRAINING_LIGANDS_WO_FLAPA = _config["ligands"]["training_wo_flapa"]
