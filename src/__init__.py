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

from .feature_engineering import LibraryLengthNormalizer, create_feature_pipeline

from .evaluation import (
    multiclass_roc_auc_score,
    make_score,
    get_confusion_matrix,
    make_probability_matrix,
)

from .config import (
    CUSTOM_PALETTE_6,
    CUSTOM_PALETTE_8,
    CUSTOM_PALETTE_9,
    CLASS_ORDER_TRAINING,
    CLASS_ORDER_OTHER_LIGANDS,
    CLASS_ORDER_BACTERIAL,
    SUBSET_PALETTES,
    SUBSET_CLASS_ORDERS,
    DESEQ2_CONFIG,
    FEATURE_SELECTION_CONFIG,
    MODEL_TRAINING_CONFIG,
    CLASS_LABELS,
    TRAINING_LIGANDS,
    ADDITIONAL_TLR_LIGANDS,
    NEGATIVE_CONTROL,
)

from .utils import get_output_path, save_fig, save_csv

from .model_training import ModelFactory, ModelTrainer

# Optional import for DESeq2 (requires pydeseq2)
try:
    from .pydeseq2 import DataProcessor as DESeq2DataProcessor, Plotter
except ImportError:
    DESeq2DataProcessor = None
    Plotter = None

# Optional import for visualization (requires goatools)
try:
    from .visualization import (
        plot_gene_expression_by_class,
        plot_tlr_hek_blue,
        plot_pca_for_pandas,
    )
except ImportError:
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
    # Feature engineering
    "LibraryLengthNormalizer",
    "create_feature_pipeline",
    # Model training
    "ModelFactory",
    "ModelTrainer",
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
    # DESeq2 analysis
    "DESeq2DataProcessor",
    "Plotter",
    # Visualization
    "plot_gene_expression_by_class",
    "plot_tlr_hek_blue",
    "plot_pca_for_pandas",
]
