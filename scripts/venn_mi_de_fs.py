"""Exploratory: 3-way Venn of MI-elbow, main DE, and the tuned forest selection.

The feature-selection circle is the deployed selection: the full feature_pipeline
(SelectKBest -> ExtraTrees SelectFromModel, using the tuned FEATURE_SELECTION_CONFIG) fit
once on the main labels, taking its selected genes. MI elbow = top genes by mutual
information with the main labels. DE (main) = union of significant genes over the per-ligand
main_ligands DESeq2 contrasts. Reuses src and the existing exploratory loaders; main pipeline
untouched.
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from src import FeatureSelector, plot_venn
from src.config import FEATURE_SELECTION_CONFIG
from de_vs_random_pca import load_main_de_genes, build_conditions
from mi_vs_random_pca import mi_ranked_genes, MI_TOP


def main():
    de_genes = load_main_de_genes()
    conditions, Xm = build_conditions()
    ym = conditions["main_ligands"][2]

    mi_genes = set(mi_ranked_genes(Xm, ym).index[:MI_TOP])

    fs = FeatureSelector.feature_pipeline(**FEATURE_SELECTION_CONFIG).set_output(
        transform="pandas"
    )
    fs.fit(Xm, ym)
    fs_genes = set(fs.get_feature_names_out())

    print(f"DE (main) = {len(de_genes)}")
    print(f"MI elbow (top {MI_TOP}) = {len(mi_genes)}")
    print(f"FS tuned selection = {len(fs_genes)}")

    plot_venn(
        [de_genes, mi_genes, fs_genes],
        set_labels=("DE (main)", f"MI elbow (top {MI_TOP})", "FS tuned selection"),
        output_filename="venn_mi_de_fs.png",
        title="Main DE genes vs MI elbow vs tuned forest selection",
    )


if __name__ == "__main__":
    main()
