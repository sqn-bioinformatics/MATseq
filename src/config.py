"""Configuration loader for MAT-seq pipeline."""

import json
from pathlib import Path
from functools import lru_cache
from typing import Any


def _find_config_path() -> Path:
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
    """Get a configuration value using dot notation, e.g. 'paths.data_dir'."""
    config = _load_config()
    for part in key.split("."):
        try:
            config = config[part]
        except (TypeError, KeyError):
            raise KeyError(f"Config key '{key}' not found")
    return config


def _project_root() -> Path:
    src_dir = Path(__file__).parent
    return src_dir.parent if src_dir.name == "src" else Path.cwd()


def expand_path(path_str: str) -> Path:
    """Expand path relative to project root or home directory."""
    path = Path(path_str)
    if path_str.startswith("~/MATseq"):
        return _project_root() / path_str.replace("~/MATseq/", "")
    if path_str.startswith("~"):
        return path.expanduser()
    if not path.is_absolute():
        return _project_root() / path
    return path


def get_featurecounts_dir() -> Path:
    return expand_path(get_config("paths.featurecounts_dir"))


def get_work_dir() -> Path:
    return expand_path(get_config("snakemake.work_dir"))


def get_sample_dir() -> Path:
    return expand_path(get_config("snakemake.sample_dir"))


def get_genome_dir() -> Path:
    return expand_path(get_config("snakemake.genome_dir"))


def get_test_sample_dir() -> Path:
    return expand_path(get_config("test.sample_dir"))


def get_test_work_dir() -> Path:
    return expand_path(get_config("test.work_dir"))


def get_test_name() -> str:
    return get_config("test.name")


_config = _load_config()

CUSTOM_PALETTE_6 = _config["colors"]["palette_6"]
CUSTOM_PALETTE_8 = _config["colors"]["palette_8"]
CUSTOM_PALETTE_9 = _config["colors"]["palette_9"]

_PALETTES = {
    "palette_6": CUSTOM_PALETTE_6,
    "palette_8": CUSTOM_PALETTE_8,
    "palette_9": CUSTOM_PALETTE_9,
}

LIGAND_ALIASES = _config.get("ligand_aliases", {})

MAIN_LIGANDS = _config["ligands"]["main_ligands"]
ADDITIONAL_LIGANDS = _config["ligands"]["additional_ligands"]
BACTERIA_LIGANDS = _config["ligands"]["bacteria_ligands"]

CLASS_ORDER_MAIN_LIGANDS = _config["class_order_for_plotting"]["main_ligands"]
CLASS_ORDER_ADDITIONAL_LIGANDS = _config["class_order_for_plotting"]["additional_ligands"]
CLASS_ORDER_BACTERIA_LIGANDS = _config["class_order_for_plotting"]["bacteria_ligands"]

SUBSET_PALETTES = {
    subset: _PALETTES[palette_name]
    for subset, palette_name in _config["subset_palettes"].items()
}

SUBSET_CLASS_ORDERS = {
    "main_ligands": CLASS_ORDER_MAIN_LIGANDS,
    "additional_ligands": CLASS_ORDER_ADDITIONAL_LIGANDS,
    "bacteria_ligands": CLASS_ORDER_BACTERIA_LIGANDS,
}


def order_labels(present, subset="main_ligands"):
    order = SUBSET_CLASS_ORDERS.get(subset, SUBSET_CLASS_ORDERS["main_ligands"])
    present = list(present)
    ordered = [c for c in order if c in present]
    return ordered + [c for c in present if c not in ordered]
CLASS_DISPLAY_NAMES = {
    "negative_control": "Negative Control",
}
SUBSET_DISPLAY_NAMES = {
    "main_ligands": "Training Ligands",
    "main_ligands_no_flapa": "Training Ligands (w/o Fla-PA)",
    "additional_ligands": "Additional Ligands",
    "bacteria_ligands": "Bacterial Ligands",
    "external_test": "Test Ligands",
    "no_flapa": "w/o Fla-PA",
}


def display_label(label):
    """Human-readable form of a single class label (no underscores)."""
    if label in CLASS_DISPLAY_NAMES:
        return CLASS_DISPLAY_NAMES[label]
    return str(label).replace("_", " ")


def display_labels(labels):
    """Human-readable forms of an iterable of class labels."""
    return [display_label(x) for x in labels]


def subset_display(subset):
    """Human-readable subset/panel name (e.g. main_ligands -> 'Training Ligands')."""
    if subset in SUBSET_DISPLAY_NAMES:
        return SUBSET_DISPLAY_NAMES[subset]
    return str(subset).replace("_", " ").title()


def confusion_title(model, subset):
    """Figure title in the '<Model> Confusion Matrix <Subset>' format."""
    return f"{model} Confusion Matrix {subset_display(subset)}"


DESEQ2_CONFIG = _config["deseq2"]
FEATURE_SELECTION_CONFIG = _config["feature_selection"]
FOREST_SELECTION_GRID = _config["forest_selection_grid"]
MODEL_FACTORY_CONFIG = _config["model_factory"]
MODEL_TRAINING_CONFIG = _config["model_training"]
HYPERPARAMETER_GRIDS = _config["hyperparameter_grids"]


def update_config(updates: dict) -> None:
    """Merge updates into config.json's feature_selection block and write it back."""
    config_path = _find_config_path()
    with config_path.open() as f:
        config = json.load(f)
    config.setdefault("feature_selection", {}).update(updates)
    with config_path.open("w") as f:
        json.dump(config, f, indent=2)
    FEATURE_SELECTION_CONFIG.update(updates)


def primary_geneset_name() -> str:
    """Geneset key for the full tuned selection, e.g. 'selected_130' (= n_selected)."""
    return f"selected_{FEATURE_SELECTION_CONFIG['max_features']}"

if FEATURE_SELECTION_CONFIG["k_best"] < FEATURE_SELECTION_CONFIG["max_features"]:
    raise ValueError(
        f"feature_selection infeasible: k_best ({FEATURE_SELECTION_CONFIG['k_best']}) "
        f"< max_features ({FEATURE_SELECTION_CONFIG['max_features']}); SelectFromModel "
        f"requires k_best >= max_features."
    )
