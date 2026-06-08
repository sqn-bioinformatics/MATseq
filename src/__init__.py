"""MAT-seq analysis package."""

__version__ = "0.1.0"

# Import commonly used functions and classes
from .preprocessing import (
    prepare_counts,
    write_count_summary,
    load_featurecounts,
    filter_low_read_samples,
    label_samples,
    extract_subset,
    normalize_rpm,
)

from .tlr_analysis import load_tlr_data, plot_tlr_hek_blue

from .feature_engineering import (
    LibraryLengthNormalizer,
    create_feature_pipeline,
    create_preprocessing_pipeline,
    create_selection_pipeline,
)

from .model_training import (
    make_score,
    ModelFactory,
    ModelTrainer,
)

from .config import (
    CUSTOM_PALETTE_6,
    CUSTOM_PALETTE_8,
    CUSTOM_PALETTE_9,
    CLASS_ORDER_MAIN_LIGANDS,
    CLASS_ORDER_ADDITIONAL_LIGANDS,
    CLASS_ORDER_BACTERIA_LIGANDS,
    SUBSET_PALETTES,
    SUBSET_CLASS_ORDERS,
    DESEQ2_CONFIG,
    FEATURE_SELECTION_CONFIG,
    MODEL_FACTORY_CONFIG,
    MODEL_TRAINING_CONFIG,
    HYPERPARAMETER_GRIDS,
    LIGAND_ALIASES,
    MAIN_LIGANDS,
    ADDITIONAL_LIGANDS,
    BACTERIA_LIGANDS,
)

from .pydeseq2 import (
    DataProcessor as DESeq2DataProcessor,
    DESeq2,
)

from .visualization import (
    plot_gene_expression_by_class,
    plot_feature_count_analysis,
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
    create_fs_de_go_table,
)

# Feature analysis
from .feature_analysis import (
    FeatureSelectionAnalyzer,
    PipelineParamTuner,
    VennDiagramGenerator,
    feature_count_analysis,
)

# Prediction and model comparison
from .prediction import ModelPredictor

__all__ = [
    # Configuration
    "CUSTOM_PALETTE_6",
    "CUSTOM_PALETTE_8",
    "CUSTOM_PALETTE_9",
    "CLASS_ORDER_MAIN_LIGANDS",
    "CLASS_ORDER_ADDITIONAL_LIGANDS",
    "CLASS_ORDER_BACTERIA_LIGANDS",
    "SUBSET_PALETTES",
    "SUBSET_CLASS_ORDERS",
    "DESEQ2_CONFIG",
    "FEATURE_SELECTION_CONFIG",
    "MODEL_FACTORY_CONFIG",
    "MODEL_TRAINING_CONFIG",
    "HYPERPARAMETER_GRIDS",
    "LIGAND_ALIASES",
    "MAIN_LIGANDS",
    "ADDITIONAL_LIGANDS",
    "BACTERIA_LIGANDS",
    # Preprocessing
    "prepare_counts",
    "write_count_summary",
    "load_featurecounts",
    "filter_low_read_samples",
    "label_samples",
    "extract_subset",
    "normalize_rpm",
    "load_tlr_data",
    # Feature engineering
    "LibraryLengthNormalizer",
    "create_feature_pipeline",
    "create_preprocessing_pipeline",
    "create_selection_pipeline",
    # Model training
    "ModelFactory",
    "ModelTrainer",
    # Evaluation
    "make_score",
    # DESeq2 analysis
    "DESeq2DataProcessor",
    "DESeq2",
    # Visualization
    "plot_gene_expression_by_class",
    "plot_feature_count_analysis",
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
    "create_fs_de_go_table",
    # Feature analysis
    "FeatureSelectionAnalyzer",
    "PipelineParamTuner",
    "VennDiagramGenerator",
    "feature_count_analysis",
    # Prediction and model comparison
    "ModelPredictor",
]
