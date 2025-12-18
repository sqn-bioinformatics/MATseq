"""MAT-seq analysis package."""

__version__ = "0.1.0"

# Import commonly used functions and classes
from .preprocessing import (
    process_featurecounts_files,
    normalize_rpm,
    prepare_training,
)
from .feature_engineering import LibraryLengthNormalizer, create_feature_pipeline
from .evaluation import (
    multiclass_roc_auc_score,
    make_score,
    get_confusion_matrix,
    make_probability_matrix,
)
from .utils import get_output_path, save_fig, save_csv, CUSTOM_PALETTE_6

# Optional import for visualization (requires goatools)
try:
    from .visualization import (
        Plotter,
        plot_gene_expression_by_class,
        load_tlr_data,
        plot_tlr_hek_blue,
    )
except ImportError:
    Plotter = None
    plot_gene_expression_by_class = None
    load_tlr_data = None
    plot_tlr_hek_blue = None

__all__ = [
    # Preprocessing
    "prepare_training",
    "normalize_rpm",
    "process_featurecounts_files",
    # Feature engineering
    "LibraryLengthNormalizer",
    "create_feature_pipeline",
    # Evaluation
    "multiclass_roc_auc_score",
    "make_score",
    "get_confusion_matrix",
    "make_probability_matrix",
    # Utils
    "get_output_path",
    "save_fig",
    "save_csv",
    "CUSTOM_PALETTE_6",
    # Visualization
    "Plotter",
    "plot_gene_expression_by_class",
    "load_tlr_data",
    "plot_tlr_hek_blue",
]
