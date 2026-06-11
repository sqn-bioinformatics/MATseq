"""Exploratory: PCA on main-subset DE genes vs equally-many random genes.

Selection = union of significant genes across the per-ligand main_ligands DESeq2
contrasts (padj < threshold & |log2FC| > threshold), i.e. the same rule as
DESeq2.de_genes in src/pydeseq2.py. For every condition (main, additional+main,
bacteria+main, external test) it draws a "selected" PCA on those DE genes and,
as a control, on 5 random N-gene sets sampled from QC-passing genes.

Reuses src loaders/pipeline/plotting; reads the DE CSVs you pointed at. Outputs
to results/figures/pca/explore_de_random/. Nothing here touches the main pipeline.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import (
    prepare_counts,
    extract_subset,
    FeatureSelector,
    plot_pca_pandas,
    SUBSET_PALETTES,
    SUBSET_CLASS_ORDERS,
    CUSTOM_PALETTE_9,
    DESEQ2_CONFIG,
)
from src.config import get_test_work_dir

RESULTS = Path.cwd() / "results"
DE_DIR = RESULTS / "differential_gene_expression" / "main_ligands"
COUNT_SUMMARY = RESULTS / "counts" / "MATseq_count_summary.csv"
OUT_SUBDIR = "explore_de_random"
RANDOM_SEEDS = [42, 43, 44, 45, 46]


def load_main_de_genes() -> set:
    """Union of significant genes over the per-ligand main_ligands contrast CSVs."""
    padj = DESEQ2_CONFIG.get("padj_threshold", 0.05)
    lfc = DESEQ2_CONFIG.get("log2fc_threshold", 2.0)
    files = sorted(DE_DIR.glob("*_deseq2_results.csv"))
    if not files:
        raise FileNotFoundError(f"No *_deseq2_results.csv in {DE_DIR}")
    genes: set = set()
    for f in files:
        df = pd.read_csv(f, index_col=0)
        sig = df[(df["padj"] < padj) & (df["log2FoldChange"].abs() > lfc)]
        print(f"  {f.name}: {len(sig)} significant")
        genes.update(sig.index)
    return genes


def make_pca(cond_key: str, palette_key: str, X, y, genes, tag: str,
             out_subdir: str = OUT_SUBDIR, equal_aspect: bool = False) -> int:
    """Preprocess X, keep `genes`, write labelled+unlabelled selected PCA. Returns N used."""
    pre = FeatureSelector.preprocessing_pipeline().set_output(transform="pandas")
    Xs = pre.fit_transform(X)
    cols = [g for g in genes if g in Xs.columns]
    Xs = Xs[cols]

    palette = SUBSET_PALETTES.get(palette_key, CUSTOM_PALETTE_9)
    hue_order = SUBSET_CLASS_ORDERS.get(palette_key)
    for with_names, label_suffix in [(False, ""), (True, "_labeled")]:
        plot_pca_pandas(
            name=f"{cond_key}_{tag}{label_suffix}",
            X=Xs,
            labels=y,
            palette=palette,
            hue_order=hue_order,
            with_sample_names=with_names,
            output_filename=f"{out_subdir}/{cond_key}_{tag}_pca{label_suffix}.png",
            equal_aspect=equal_aspect,
        )
    return len(cols)


def load_training_counts():
    """Reuse the pipeline's samples x genes (+label) summary CSV instead of re-reading featureCounts."""
    df = pd.read_csv(COUNT_SUMMARY, index_col=0)
    return df.drop(columns=["label"]), df["label"]


def build_conditions():
    """(condition_key -> (palette_key, X, y)) mirroring the pipeline's PCA grouping."""
    features, labels = load_training_counts()
    Xm, ym = extract_subset(features, labels, "main_ligands")

    conditions = {"main_ligands": ("main_ligands", Xm, ym)}
    for sub in ["additional_ligands", "bacteria_ligands"]:
        Xs, ys = extract_subset(features, labels, sub)
        conditions[sub] = (sub, pd.concat([Xs, Xm]), pd.concat([ys, ym]))

    test_fc = get_test_work_dir() / "featurecounts"
    if test_fc.is_dir() and any(test_fc.glob("*.txt")):
        tf, tl = prepare_counts(featurecounts_dir=test_fc)
        Xe, ye = extract_subset(tf, tl, "main_ligands")
        conditions["external_test"] = ("main_ligands", Xe, ye)
    else:
        print(f"  external test skipped: no featurecounts at {test_fc}")

    return conditions, Xm


def main():
    (RESULTS / "figures" / "pca" / OUT_SUBDIR).mkdir(parents=True, exist_ok=True)

    print("Loading main-subset DE genes...")
    de_genes = load_main_de_genes()

    conditions, Xm = build_conditions()
    de_present = sorted(de_genes & set(Xm.columns))
    n = len(de_present)
    print(f"\n=== N = {n} DE genes (main subset, present in count matrix) ===\n")

    qc_pool = set(Xm.columns[Xm.sum(axis=0) >= 10])
    if "external_test" in conditions:
        Xe = conditions["external_test"][1]
        qc_pool &= set(Xe.columns[Xe.sum(axis=0) >= 10])
    qc_genes = sorted(qc_pool)
    print(f"Random pool (QC-passing in both batches): {len(qc_genes)} genes")
    random_sets = {
        seed: list(np.random.default_rng(seed).choice(qc_genes, size=n, replace=False))
        for seed in RANDOM_SEEDS
    }

    summary = []
    for cond_key, (palette_key, X, y) in conditions.items():
        used = make_pca(cond_key, palette_key, X, y, de_present, "de")
        print(f"{cond_key}: DE genes used = {used}")
        summary.append({"condition": cond_key, "gene_set": "de", "n_genes": used})
        for seed, genes in random_sets.items():
            used_r = make_pca(cond_key, palette_key, X, y, genes, f"random_seed{seed}")
            summary.append(
                {"condition": cond_key, "gene_set": f"random_seed{seed}", "n_genes": used_r}
            )

    out = RESULTS / "figures" / "pca" / OUT_SUBDIR / "gene_counts.csv"
    pd.DataFrame(summary).to_csv(out, index=False)
    print(f"\nN (DE genes) = {n}")
    print(f"Figures: {RESULTS / 'figures' / 'pca' / OUT_SUBDIR}")
    print(f"Gene-count summary: {out}")


if __name__ == "__main__":
    main()
