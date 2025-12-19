"""Configuration constants for MAT-seq pipeline."""

# Color palettes for different subsets
CUSTOM_PALETTE_6 = [
    "#1f77b4",  # negative_control - Muted Blue
    "#ff7f0e",  # LPS - Soft Orange
    "#2ca02c",  # Pam3 - Green
    "#d62728",  # R848 - Red
    "#9467bd",  # PGN - Purple
    "#17becf",  # Fla-PA - Teal/Cyan
]

CUSTOM_PALETTE_8 = [
    "#1f77b4",  # negative_control - Muted Blue
    "#ff7f0e",  # LPS - Soft Orange
    "#2ca02c",  # Pam3 - Green
    "#d62728",  # R848 - Red
    "#9467bd",  # PGN - Purple
    "#17becf",  # Fla-PA - Teal/Cyan
    "#e377c2",  # HK E.coli - Pink
    "#bcbd22",  # HK S.aureus - Yellow-Green
]

CUSTOM_PALETTE_9 = [
    "#1f77b4",  # negative_control - Muted Blue
    "#ff7f0e",  # LPS - Soft Orange
    "#2ca02c",  # Pam3 - Green
    "#d62728",  # R848 - Red
    "#9467bd",  # PGN - Purple
    "#17becf",  # Fla-PA - Teal/Cyan
    "#8c564b",  # LTA - Brown
    "#00fa9a",  # MPLA - Medium Spring Green
    "#e377c2",  # Pam2 - Pink
]

# Class ordering for different data subsets
CLASS_ORDER_TRAINING = [
    "negative_control",
    "LPS",
    "Pam3",
    "R848",
    "PGN",
    "Fla-PA",
]

CLASS_ORDER_OTHER_LIGANDS = [
    "negative_control",
    "LPS",
    "Pam3",
    "R848",
    "PGN",
    "Fla-PA",
    "LTA",
    "MPLA",
    "Pam2",
]

CLASS_ORDER_BACTERIAL = [
    "negative_control",
    "LPS",
    "Pam3",
    "R848",
    "PGN",
    "Fla-PA",
    "HK E.coli",
    "HK S.aureus",
]

# Mapping of subset names to color palettes
SUBSET_PALETTES = {
    "training": CUSTOM_PALETTE_6,
    "other_ligands": CUSTOM_PALETTE_9,
    "bacterial": CUSTOM_PALETTE_8,
}

# Mapping of subset names to class orders
SUBSET_CLASS_ORDERS = {
    "training": CLASS_ORDER_TRAINING,
    "other_ligands": CLASS_ORDER_OTHER_LIGANDS,
    "bacterial": CLASS_ORDER_BACTERIAL,
}

# DESeq2 Analysis parameters
DESEQ2_CONFIG = {
    "padj_threshold": 0.05,
    "log2fc_threshold": 2,
    "baseMean_threshold": 10,
    "n_cpus": 42,
}

# Feature selection parameters
FEATURE_SELECTION_CONFIG = {
    "k_best": 1000,
    "n_estimators": 250,
    "max_depth": 5,
    "max_features": 250,
    "feature_threshold": 0.001,
    "random_state": 42,
}

# Model training parameters
MODEL_TRAINING_CONFIG = {
    "random_state": 42,
    "apply_smote": True,
    "smote_sampling_strategy": "not majority",
    "smote_k_neighbors": 1,
}

# Class labels mapping
CLASS_LABELS = {
    "IMDM": "negative_control",
    "nc": "negative_control",
}

# DESeq2 ligand lists for different analyses
TRAINING_LIGANDS = ["Fla-PA", "LPS", "PGN", "R848", "Pam3", "PGN"]

ADDITIONAL_TLR_LIGANDS = ["LTA", "MPLA", "Pam2", "HKEB", "HKSA"]

# Negative control identifier
NEGATIVE_CONTROL = "nc"
