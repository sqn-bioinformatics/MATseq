"""Exploratory: external validation-set PCA using the deployed ML selector.

The deployed selector = FeatureSelector.feature_pipeline(**FEATURE_SELECTION_CONFIG)
(preprocess -> SelectKBest(MI, k=1000) -> SelectFromModel(ExtraTrees, max_features=200)),
fit on the training main_ligands subset, then .transform applied to the external 7086
batch -- identical to MATseq.py STEP 5. Reuses src pipeline/plot; main pipeline untouched.
"""

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from src import (
    prepare_counts,
    extract_subset,
    FeatureSelector,
    plot_pca_pandas,
    SUBSET_PALETTES,
    SUBSET_CLASS_ORDERS,
    FEATURE_SELECTION_CONFIG,
)
from src.config import get_test_work_dir
from de_vs_random_pca import RESULTS, load_training_counts

OUT_SUBDIR = "explore_ml"


def main():
    (RESULTS / "figures" / "pca" / OUT_SUBDIR).mkdir(parents=True, exist_ok=True)

    features, labels = load_training_counts()
    Xm, ym = extract_subset(features, labels, "main_ligands")

    fs = FeatureSelector.feature_pipeline(**FEATURE_SELECTION_CONFIG).set_output(
        transform="pandas"
    )
    fs.fit(Xm, ym)

    test_fc = get_test_work_dir() / "featurecounts"
    tf, tl = prepare_counts(featurecounts_dir=test_fc)
    Xe, ye = extract_subset(tf, tl, "main_ligands")
    Xe = Xe.reindex(columns=Xm.columns, fill_value=0)

    palette = SUBSET_PALETTES["main_ligands"]
    hue_order = SUBSET_CLASS_ORDERS["main_ligands"]

    def plot(name, X):
        for with_names, label_suffix in [(False, ""), (True, "_labeled")]:
            plot_pca_pandas(
                name=f"{name}{label_suffix}",
                X=X,
                labels=ye,
                palette=palette,
                hue_order=hue_order,
                with_sample_names=with_names,
                output_filename=f"{OUT_SUBDIR}/{name}_pca{label_suffix}.png",
            )

    X_before = FeatureSelector.preprocessing_pipeline().set_output(
        transform="pandas"
    ).fit_transform(Xe)
    print(f"Before selection: {X_before.shape[1]} genes on {X_before.shape[0]} external samples")
    plot("external_test_before", X_before)

    X_sel = fs.transform(Xe)
    print(f"Deployed ML selector: {X_sel.shape[1]} genes on {X_sel.shape[0]} external samples")
    plot("external_test_ml", X_sel)

    print(f"Figures: {RESULTS / 'figures' / 'pca' / OUT_SUBDIR}")


if __name__ == "__main__":
    main()
