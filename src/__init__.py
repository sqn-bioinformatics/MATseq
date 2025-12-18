"""MAT-seq analysis package."""

__version__ = "0.1.0"

# Import commonly used functions and classes
from .preprocessing import (
    merge_counts,
    filter_counts,
    extract_subset,
    normalize_rpm,
    load_tlr_data,
)
# Optional import for DESeq2 (requires pydeseq2)
try:
    from .deseq2 import DataProcessor
except ImportError:
    DataProcessor = None

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
        plot_tlr_hek_blue,
        plot_pca_for_pandas,
    )
except ImportError:
    Plotter = None
    plot_gene_expression_by_class = None
    plot_tlr_hek_blue = None
    plot_pca_for_pandas = None

__all__ = [
    # Preprocessing
    "merge_counts",
    "filter_counts",
    "extract_subset",
    "normalize_rpm",
    "load_tlr_data",
    # DESeq2
    "DataProcessor",
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
    "plot_tlr_hek_blue",
    "plot_pca_for_pandas",
]
