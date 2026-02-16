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

from .model_training import (
    multiclass_roc_auc_score,
    make_score,
    get_confusion_matrix,
    ModelFactory,
    ModelTrainer,
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
    MODEL_FACTORY_CONFIG,
    MODEL_TRAINING_CONFIG,
    TRAINING_LIGANDS,
    TRAINING_LIGANDS_WO_FLAPA,
    ADDITIONAL_LIGANDS,
    BACTERIAL_LIGANDS,
    BACTERIAL_LIGANDS_ORIGINAL_NAMES,
)

from .pydeseq2 import (
    DataProcessor as DESeq2DataProcessor,
    AnalysisPipeline,
)

from .visualization import (
    plot_gene_expression_by_class,
    plot_tlr_hek_blue,
    plot_pca_pandas,
    plot_volcano,
    plot_heatmap,
    plot_pca_deseq2,
    plot_go,
)

from .go_term_analysis import (
    initialize_go,
    generate_go_table,
    run_go_analysis,
    merge_go_tables,
)

# Feature analysis
from .feature_analysis import (
    FeatureSelectionAnalyzer,
    VennDiagramGenerator,
)

# Prediction and model comparison
from .prediction import ModelPredictor, ModelComparator

__all__ = [
    # Configuration
    "CUSTOM_PALETTE_6",
    "CUSTOM_PALETTE_8",
    "CUSTOM_PALETTE_9",
    "CLASS_ORDER_TRAINING",
    "CLASS_ORDER_OTHER_LIGANDS",
    "CLASS_ORDER_BACTERIAL",
    "SUBSET_PALETTES",
    "SUBSET_CLASS_ORDERS",
    "DESEQ2_CONFIG",
    "FEATURE_SELECTION_CONFIG",
    "MODEL_FACTORY_CONFIG",
    "MODEL_TRAINING_CONFIG",
    "TRAINING_LIGANDS",
    "ADDITIONAL_LIGANDS",
    "BACTERIAL_LIGANDS",
    "BACTERIAL_LIGANDS_ORIGINAL_NAMES",
    "TRAINING_LIGANDS_WO_FLAPA",
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
    # DESeq2 analysis
    "DESeq2DataProcessor",
    "AnalysisPipeline",
    # Visualization
    "plot_gene_expression_by_class",
    "plot_tlr_hek_blue",
    "plot_pca_pandas",
    "plot_volcano",
    "plot_heatmap",
    "plot_pca_deseq2",
    "plot_go",
    # GO term analysis
    "initialize_go",
    "generate_go_table",
    "run_go_analysis",
    "merge_go_tables",
    # Feature analysis
    "FeatureSelectionAnalyzer",
    "VennDiagramGenerator",
    # Prediction and model comparison
    "ModelPredictor",
    "ModelComparator",
]
