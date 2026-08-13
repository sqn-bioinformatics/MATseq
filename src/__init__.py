"""MAT-seq analysis package."""

__version__ = "0.1.0"

# Import commonly used functions and classes
from .preprocessing import (
    prepare_counts,
    load_featurecounts,
    filter_low_read_samples,
    label_samples,
    extract_subset,
    normalize_rpm,
)

from .tlr_analysis import load_tlr_data, plot_tlr_hek_blue

from .feature_engineering import (
    LibraryLengthNormalizer,
    ColumnSelector,
    feature_pipeline,
    preprocessing_pipeline,
    selected_with_importance,
    write_selected_gene_tables,
    mutual_information,
    forest_kmeans,
)

from .model_training import (
    make_score,
    ModelFactory,
    ModelTrainer,
)

from .config import (
    CUSTOM_PALETTE_9,
    SUBSET_PALETTES,
    CLASS_ORDER,
    DESEQ2_CONFIG,
    FEATURE_SELECTION_CONFIG,
    MODEL_FACTORY_CONFIG,
    MODEL_TRAINING_CONFIG,
    HYPERPARAMETER_GRIDS,
    LIGAND_ALIASES,
    MAIN_LIGANDS,
    ADDITIONAL_LIGANDS,
    BACTERIAL_LIGANDS,
)

from .pydeseq2 import (
    DataProcessor as DESeq2DataProcessor,
    DESeq2,
)

from .visualization import (
    plot_mutual_information,
    plot_forest_ari_sweep,
    plot_venn,
    plot_pca,
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

# Prediction and model comparison
from .prediction import ModelPredictor

__all__ = [
    # Configuration
    "CUSTOM_PALETTE_9",
    "SUBSET_PALETTES",
    "CLASS_ORDER",
    "DESEQ2_CONFIG",
    "FEATURE_SELECTION_CONFIG",
    "MODEL_FACTORY_CONFIG",
    "MODEL_TRAINING_CONFIG",
    "HYPERPARAMETER_GRIDS",
    "LIGAND_ALIASES",
    "MAIN_LIGANDS",
    "ADDITIONAL_LIGANDS",
    "BACTERIAL_LIGANDS",
    # Preprocessing
    "prepare_counts",
    "load_featurecounts",
    "filter_low_read_samples",
    "label_samples",
    "extract_subset",
    "normalize_rpm",
    "load_tlr_data",
    # Feature engineering
    "LibraryLengthNormalizer",
    "ColumnSelector",
    "feature_pipeline",
    "preprocessing_pipeline",
    "selected_with_importance",
    "write_selected_gene_tables",
    "mutual_information",
    "forest_kmeans",
    # Model training
    "ModelFactory",
    "ModelTrainer",
    # Evaluation
    "make_score",
    # DESeq2 analysis
    "DESeq2DataProcessor",
    "DESeq2",
    # Visualization
    "plot_mutual_information",
    "plot_forest_ari_sweep",
    "plot_venn",
    "plot_tlr_hek_blue",
    "plot_pca",
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
    # Prediction and model comparison
    "ModelPredictor",
]
