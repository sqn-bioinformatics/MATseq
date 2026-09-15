"""MAT-seq analysis package."""

__version__ = "0.1.0"

from .config import (
    SUBSET_PALETTES,
    CLASS_ORDER,
    SUBSET_DISPLAY_NAMES,
    DESEQ2_CONFIG,
    FEATURE_SELECTION_CONFIG,
    MODEL_TRAINING_CONFIG,
    HYPERPARAMETER_GRIDS,
)
from .preprocessing import prepare_counts, extract_subset
from .feature_engineering import (
    feature_pipeline,
    preprocessing_pipeline,
    selected_with_importance,
    mutual_information,
    forest_kmeans,
)
from .model_training import ModelTrainer
from .pydeseq2 import DESeq2
from .visualization import (
    plot_mutual_information,
    plot_forest_ari_sweep,
    plot_venn,
    plot_pca,
)
from .go_term_analysis import create_fs_de_go_table, initialize_go
from .tlr_analysis import load_tlr_data, plot_tlr_hek_blue
from .prediction import predict_samples
from .make_tables import format_table2, assemble_supplementary_tables
