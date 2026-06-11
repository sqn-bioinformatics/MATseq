"""Exploratory: 3-way Venn of MI-elbow, main DE, and the forest selection (max_features=100).

Mirrors results/figures/venn/venn_de_vs_fs_methods.png but replaces the stability-knee circle
with the deployed feature selection: the full feature_pipeline (SelectKBest k=4063 -> ExtraTrees
SelectFromModel max_features=100, from FEATURE_SELECTION_CONFIG) reseeded 1000 times, taking the
union of selected genes. MI elbow = top 4063 genes by mutual information with the main labels.
DE (main) = union of significant genes over the per-ligand main_ligands DESeq2 contrasts.
Reuses src and the existing exploratory loaders; main pipeline untouched.
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from src import FeatureSelector, VennDiagramGenerator
from de_vs_random_pca import load_main_de_genes, build_conditions
from mi_vs_random_pca import mi_ranked_genes, MI_TOP

N_RUNS = 1000


def main():
    de_genes = load_main_de_genes()
    conditions, Xm = build_conditions()
    ym = conditions["main_ligands"][2]

    mi_genes = set(mi_ranked_genes(Xm, ym).index[:MI_TOP])

    selector = FeatureSelector()
    selector.run_multiple_selections(X=Xm, y=ym, n_runs=N_RUNS)
    fs_genes = set().union(*selector.feature_sets)

    print(f"DE (main) = {len(de_genes)}")
    print(f"MI elbow (top {MI_TOP}) = {len(mi_genes)}")
    print(f"FS forest 100, union of {N_RUNS} runs = {len(fs_genes)}")

    VennDiagramGenerator().plot_venn3(
        [de_genes, mi_genes, fs_genes],
        set_labels=("DE (main)", f"MI elbow (top {MI_TOP})", f"FS forest 100 ({N_RUNS} runs)"),
        output_filename="venn_mi_de_fs.png",
        title="Main DE genes vs MI elbow vs forest selection (max_features=100)",
    )


if __name__ == "__main__":
    main()
