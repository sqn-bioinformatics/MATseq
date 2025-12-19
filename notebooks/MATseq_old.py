# Import from src modules
from sklearn.preprocessing import StandardScaler
from src import (
    # Preprocessing
    merge_counts,
    extract_subset,
    normalize_rpm,
    # Feature engineering
    create_feature_pipeline,
    # Visualization
    plot_pca_for_pandas,
    CUSTOM_PALETTE_6,
)


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

# Custom palettes matching class order
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

# Map subset names to palettes and class orders
SUBSET_PALETTES = {
    "training": CUSTOM_PALETTE_6,
    "other_ligands": CUSTOM_PALETTE_9,
    "bacterial": CUSTOM_PALETTE_8,
}

SUBSET_CLASS_ORDERS = {
    "training": CLASS_ORDER_TRAINING,
    "other_ligands": CLASS_ORDER_OTHER_LIGANDS,
    "bacterial": CLASS_ORDER_BACTERIAL,
}


def extract_subset(df, subset_name: str):

    subset_data, labels = extract_subset(df, name=subset_name)
    X = subset_data.drop(columns=["label"])
    X_rpm = normalize_rpm(X)
    X_scaled = StandardScaler().fit_transform(X_rpm)

    # Get appropriate palette and class order for this subset
    palette = SUBSET_PALETTES.get(subset_name, CUSTOM_PALETTE_6)
    hue_order = SUBSET_CLASS_ORDERS.get(subset_name)

    plot_pca_for_pandas(
        name=f"{subset_name}_normalized",
        df=X_scaled,
        labels=labels,
        with_sample_names=False,
        output_filename=f"{subset_name}_pca_before_feature_selection.png",
        palette=palette,
        hue_order=hue_order,
    )
    plot_pca_for_pandas(
        name=f"{subset_name}_pca_before_feature_selection_labeled",
        df=X_scaled,
        labels=labels,
        with_sample_names=True,
        output_filename=f"{subset_name}_pca_normalized_labeled.png",
        palette=palette,
        hue_order=hue_order,
    )

    return subset_data, labels


def run_feature_selection(X, y, labels, subset_name: str):
    """Run feature selection pipeline and create PCA plots on selected features."""
    print(f"\nRunning feature selection for {subset_name}...")

    # Create and fit the feature selection pipeline
    pipe = create_feature_pipeline()
    X_selected = pipe.fit_transform(X, y)
    palette = SUBSET_PALETTES.get(subset_name, CUSTOM_PALETTE_6)
    hue_order = SUBSET_CLASS_ORDERS.get(subset_name)

    print(f"Creating PCA plots for {subset_name} (feature-selected)...")
    plot_pca_for_pandas(
        name=f"{subset_name}_selected",
        df=X_selected,
        labels=labels,
        with_sample_names=False,
        output_filename=f"{subset_name}_pca_selected.png",
        palette=palette,
        hue_order=hue_order,
    )
    plot_pca_for_pandas(
        name=f"{subset_name}_selected_labeled",
        df=X_selected,
        labels=labels,
        with_sample_names=True,
        output_filename=f"{subset_name}_pca_selected_labeled.png",
        palette=palette,
        hue_order=hue_order,
    )

    return X_selected, pipe


def main():
    # Step 1: Merge featurecounts files into a single DataFrame
    print("Step 1: Merging featurecounts files...")
    df = merge_counts()
    print(f"Merged data shape: {df.shape}")

    subsets = ["training", "other_ligands", "bacterial"]
    results = {}

    for subset_name in subsets:
        X, labels = process_subset(df, subset_name)
        X_selected, pipe = run_feature_selection(X, y, labels, subset_name)
        results[subset_name] = {
            "X": X,
            "labels": labels,
        }

    return results


if __name__ == "__main__":
    main()
