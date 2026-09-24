import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import textwrap
from itertools import cycle
from pathlib import Path
from matplotlib.patches import Patch
from matplotlib_venn import venn2
from sklearn.decomposition import PCA
from adjustText import adjust_text

from .config import CLASS_DISPLAY_NAMES


def plot_confusion_matrix(cm, class_names, output_path: Path, output_filename: str,
                          title: str | None = None) -> Path:
    """Render a confusion matrix normalized over the true classes (rows)."""
    fig, ax = plt.subplots(figsize=(6.5, 6))
    cm = np.asarray(cm, dtype=float)
    labels = [CLASS_DISPLAY_NAMES.get(c, c) for c in class_names]
    n = cm.shape[0]
    annot_fs = 15 if n <= 6 else (13 if n == 7 else 12)
    tick_fs = 18 if n <= 7 else 16
    im = ax.imshow(cm, cmap="Blues", vmin=0.0, vmax=1.0, aspect="auto")
    ax.set_box_aspect(1)

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=tick_fs)
    ax.set_yticklabels(labels, fontsize=tick_fs)
    ax.set_xlabel("Predicted class", fontsize=20)
    ax.set_ylabel("True class", fontsize=20)
    if title:
        ax.set_title(title, fontsize=20, pad=8)

    for i in range(n):
        for j in range(n):
            v = cm[i, j]
            # Drop the decimal for a full 100% so it never overruns the cell.
            txt = "100%" if v >= 0.9995 else f"{v * 100:.1f}%"
            ax.text(j, i, txt, ha="center", va="center",
                    color="white" if v > 0.5 else "#222222", fontsize=annot_fs)

    ax.tick_params(which="both", length=0)
    ax.spines[:].set_visible(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(length=0, labelsize=16)

    output_path.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    save_path = output_path / output_filename
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return save_path


def plot_probability_heatmap(proba_df, class_order, true_labels, output_path: Path,
                             output_filename: str, title: str | None = None,
                             seed: int = 42) -> Path:
    """Render a per-sample prediction-probability heatmap, one sample per control class."""
    rng = np.random.default_rng(seed)
    control_idx = []
    for cls in ("negative_control", "LPS"):
        idx = true_labels.index[true_labels == cls]
        control_idx += list(rng.choice(idx, size=1, replace=False))
    remaining = [
        i
        for cls in class_order
        if cls not in ("negative_control", "LPS")
        for i in true_labels[true_labels == cls].index
    ]
    ordered_idx = control_idx + remaining
    available_classes = [c for c in class_order if c in proba_df.columns]

    fig, ax = plt.subplots(figsize=(12, 8))
    mat = proba_df.loc[ordered_idx, available_classes].to_numpy(dtype=float)
    col_labels = [CLASS_DISPLAY_NAMES.get(c, c) for c in available_classes]
    row_labels = [CLASS_DISPLAY_NAMES.get(c, c) for c in true_labels.loc[ordered_idx]]
    nrow, ncol = mat.shape

    im = ax.imshow(mat, cmap="YlGnBu", vmin=0.0, vmax=1.0, aspect="auto")
    ax.set_box_aspect(1)

    tick_fs = 9 if max(nrow, ncol) <= 8 else 7
    ax.set_xticks(range(ncol))
    ax.set_yticks(range(nrow))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=tick_fs)
    ax.set_yticklabels(row_labels, fontsize=tick_fs)

    ax.set_xlabel("Reference class", fontsize=10)
    ax.set_ylabel("Sample", fontsize=10)
    if title:
        ax.set_title(title, fontsize=11, pad=8)

    ax.set_xticks(np.arange(-0.5, ncol, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, nrow, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="both", length=0)
    ax.spines[:].set_visible(False)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Predicted class probability")
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(length=0)

    output_path.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    save_path = output_path / output_filename
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return save_path


def plot_venn(
    sets: list,
    set_labels: tuple,
    output_path: Path,
    output_filename: str = "venn.png",
    title: str | None = None,
) -> Path:
    """Plot a 2-set Venn diagram and save it."""
    with plt.rc_context({"font.size": 22}):
        fig = plt.figure(figsize=(8, 8))
        venn2([set(s) for s in sets],
              set_labels=[textwrap.fill(lb, 16) for lb in set_labels])
    if title:
        plt.title(title, fontsize=22)

    save_path = output_path / output_filename
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Figure saved: {save_path}")
    plt.close(fig)

    return save_path


def plot_mutual_information(
    result: dict,
    output_path: Path,
    output_filename: str = "mutual_information.png",
) -> Path:
    """Plot the sorted mutual-information curve with its elbow."""
    mi_elbow = result["mi_elbow"]
    scores = result["scores"]

    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax1.plot(scores["rank"], scores["mi_sorted"])
    ax1.axvline(
        mi_elbow, color="r", ls="--", label=f"mean elbow = {mi_elbow}"
    )
    ax1.set_xlabel("Gene rank", fontsize=19)
    ax1.set_ylabel("Mutual information", fontsize=19)
    ax1.tick_params(labelsize=17)
    ax1.legend(frameon=False, loc="upper right", fontsize=17)
    ax1.spines[["top", "right"]].set_visible(False)

    save_path = output_path / output_filename
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Figure saved: {save_path}")
    plt.close(fig)

    return save_path


def plot_forest_ari_sweep(
    scan: pd.DataFrame,
    output_path: Path,
    output_filename: str = "forest_ari_sweep.png",
    title: str = "k-means separation per gene rank",
) -> Path:
    """Mean and seed range of k-means/ligand ARI over the ExtraTrees gene rank."""
    ari_cols = [c for c in scan.columns if c.startswith("ari_seed_")]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.fill_between(
        scan["n_selected"], scan[ari_cols].min(axis=1), scan[ari_cols].max(axis=1),
        color="#1f77b4", alpha=0.2, linewidth=0,
        label=f"range over {len(ari_cols)} seeds",
    )
    ax.plot(
        scan["n_selected"], scan["ari_mean"], color="#1f77b4", linewidth=1.5,
        label="mean",
    )
    ax.axhline(
        scan["ari_mean"].max(), color="grey", ls=":", linewidth=1,
        label="max mean ARI",
    )
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_title(title, fontsize=20, pad=8)
    ax.set_xlabel("Gene rank (ExtraTrees importance)", fontsize=18)
    ax.set_ylabel("Adjusted Rand Index\n(k-means vs. ligand class)", fontsize=18)
    ax.tick_params(labelsize=16)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="best", fontsize=15)

    save_path = output_path / output_filename
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Figure saved: {save_path}")
    plt.close(fig)

    return save_path


def plot_pca(
    X: pd.DataFrame,
    labels: pd.Series,
    output_filename: str,
    output_path: Path,
    with_sample_names: bool = False,
    palette: str | list | None = None,
    hue_order: list | None = None,
) -> Path:
    """Create PCA visualization for pandas DataFrame data."""
    X_reduced = PCA(n_components=2).fit_transform(X)

    with plt.rc_context({"figure.facecolor": "white", "axes.facecolor": "white"}):
        fig, ax = plt.subplots(figsize=(20, 20) if with_sample_names else (6, 6))
        sns.scatterplot(
            x=X_reduced[:, 0],
            y=X_reduced[:, 1],
            hue=labels,
            hue_order=hue_order,
            s=200 if with_sample_names else 80,
            alpha=0.6,
            palette=palette,
            ax=ax,
        )
        ax.set_xlabel("PC1", fontsize=18)
        ax.set_ylabel("PC2", fontsize=18)
        ax.tick_params(axis="both", labelsize=16)

        if with_sample_names:
            texts = [
                ax.text(x, y, name, ha="left", va="bottom", alpha=0.8, fontsize=12)
                for (x, y), name in zip(X_reduced, X.index)
            ]
            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle="->", color="black"))

        handles, labels_txt = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles,
                      [CLASS_DISPLAY_NAMES.get(lb, lb) for lb in labels_txt],
                      loc="upper left", bbox_to_anchor=(1.02, 1.0),
                      borderaxespad=0, ncol=1, fontsize=15,
                      frameon=False, handletextpad=0.4)

        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_linewidth(1.5)
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_box_aspect(1)
        plt.tight_layout()

        output_path.mkdir(parents=True, exist_ok=True)
        save_path = output_path / output_filename
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Figure saved to: {save_path.absolute()}")
        plt.close(fig)

    return save_path.absolute()


def plot_volcano(
    res: pd.DataFrame,
    analysis_name: str,
    output_path: Path,
    log2foldchange: float = 2.0,
) -> Path:
    """Create volcano plot showing differentially expressed genes."""
    grapher = res.assign(
        padj_log=-np.log10(res["padj"].replace(0, 1e-300)),
        color="no_expression_change",
    )

    grapher.loc[grapher["log2FoldChange"] > log2foldchange, "color"] = "overexpressed"
    grapher.loc[grapher["log2FoldChange"] < -log2foldchange, "color"] = "underexpressed"

    grapher_subset = grapher[grapher["color"].isin(["overexpressed", "underexpressed"])]

    sorted_grapher_padj_log = grapher_subset.sort_values(by="padj_log", ascending=False)
    sorted_grapher_log2foldchange = grapher_subset.sort_values(
        by="log2FoldChange", ascending=True
    )

    annotation_subset = pd.concat(
        [
            sorted_grapher_padj_log.head(12),
            sorted_grapher_log2foldchange.head(4),
            sorted_grapher_log2foldchange.tail(4),
        ]
    ).drop_duplicates()

    fig = plt.figure(figsize=(5.5, 5.55))
    ax = sns.scatterplot(
        data=grapher,
        x="log2FoldChange",
        y="padj_log",
        hue="color",
        hue_order=["no_expression_change", "overexpressed", "underexpressed"],
        palette=["grey", "orange", "purple"],
        alpha=0.7,
    )

    ax.axhline(1.3, color="black", linestyle="--", linewidth=1)
    ax.axvline(log2foldchange, color="black", linestyle="--", linewidth=1)
    ax.axvline(-log2foldchange, color="black", linestyle="--", linewidth=1)

    texts = [
        plt.text(x=row.log2FoldChange, y=row.padj_log, s=row.name,
                 weight="bold", size=7)
        for _, row in annotation_subset.iterrows()
    ]

    plt.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False,
               prop={"size": 9})
    plt.xticks(size=10)
    plt.yticks(size=10)
    plt.xlabel("$log_{2}$ fold change", fontsize=11)
    plt.ylabel("-$log_{10}$ FDR", fontsize=11)
    plt.ylim(-2, grapher["padj_log"].max() * 1.15)
    plt.title(
        f"{analysis_name} Differentially Expressed Genes",
        fontsize=12,
        pad=20,
    )

    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(1.5)

    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="k"))

    fig.patch.set_facecolor("white")

    output_path.mkdir(parents=True, exist_ok=True)
    save_path = output_path / f"{analysis_name}_volcano.png"
    fig.savefig(save_path, dpi=300, bbox_inches="tight", pad_inches=0.25)
    print(f"Figure saved: {save_path}")
    plt.close(fig)

    return save_path


def plot_heatmap(
    dds,
    sigs: pd.DataFrame,
    analysis_name: str,
    output_path: Path,
) -> Path:
    """Create hierarchical clustering heatmap of the top 50 significant genes."""
    dds_sigs = dds[:, sigs.sort_values("padj").index[:50]]
    grapher = pd.DataFrame(
        np.log1p(dds_sigs.layers["normed_counts"]).T,
        index=dds_sigs.var_names,
        columns=dds_sigs.obs.condition,
    )

    other_colors = cycle(["green", "tab:blue", "tab:orange", "tab:purple", "tab:brown"])
    lut = {}
    for cond in dict.fromkeys(dds_sigs.obs.condition):
        c = str(cond).lower()
        lut[cond] = "magenta" if "negative" in c or c == "control" else next(other_colors)
    col_colors = list(dds_sigs.obs.condition.map(lut))
    g = sns.clustermap(
        figsize=(5.5, 5.5),
        data=grapher,
        cmap="RdYlBu_r",
        z_score=0,
        dendrogram_ratio=(0.1, 0.1),
        cbar_pos=(0.93, 0.2, 0.03, 0.45),
        cbar_kws=dict(
            location="left",
            orientation="vertical",
            pad=2,
            label="Row z-score (log1p normed counts)",
        ),
        col_colors=col_colors,
    )

    vmin, vmax = g.ax_heatmap.collections[0].get_clim()
    g.ax_cbar.set_yticks([vmin, 0, vmax])
    g.ax_cbar.set_yticklabels([f"{vmin:.1f}", "0", f"{vmax:.1f}"], fontsize=5)
    g.ax_cbar.yaxis.label.set_size(5)

    handles = [Patch(facecolor=lut[name], label=name) for name in lut]

    g.ax_col_dendrogram.legend(
        handles=handles,
        bbox_to_anchor=(1.0, 1.0),
        loc="upper left",
        fontsize=7,
        fancybox=True,
        frameon=True,
        facecolor="white",
        edgecolor="black",
    )

    g.ax_heatmap.set_xticklabels([])
    g.ax_heatmap.tick_params(bottom=False)
    g.ax_heatmap.set(xlabel=None)

    reordered_genes = grapher.index[g.dendrogram_row.reordered_ind]
    g.ax_heatmap.set_yticks(np.arange(len(reordered_genes)) + 0.5)
    g.ax_heatmap.set_yticklabels(reordered_genes, fontsize=5)

    g.ax_col_dendrogram.set_title(
        f"{analysis_name} Differentially Expressed Genes",
        fontsize=12,
        pad=20,
    )

    g.figure.subplots_adjust(hspace=0.01, right=0.82)
    g.figure.patch.set_facecolor("white")
    output_path.mkdir(parents=True, exist_ok=True)

    save_path = output_path / f"{analysis_name}_histogram.png"
    g.figure.savefig(save_path, dpi=300, bbox_inches="tight", pad_inches=0.25)
    print(f"Figure saved: {save_path}")
    plt.close(g.figure)

    return save_path


def plot_go(
    go_df: pd.DataFrame,
    output_path: Path,
    output_filename: str,
    condition: str,
    title: str = "Enriched GO Terms",
) -> Path:
    """Create horizontal bar plot of GO enrichment terms."""
    go_terms = go_df.head(15).sort_values("ratio_in_study", ascending=False)

    if len(go_terms) == 0:
        raise ValueError("No GO terms remaining after filtering")

    norm = mpl.colors.LogNorm(vmin=go_terms.fdr.min(), vmax=go_terms.fdr.max())
    color_mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.bwr_r)

    fig = plt.figure(figsize=(8, 10))

    ax = sns.barplot(
        data=go_terms,
        x=go_terms["n_genes"] / go_terms["n_study"],
        y="term",
        palette=list(color_mapper.to_rgba(go_terms.fdr.values)),
    )

    ax.set_yticklabels([textwrap.fill(term, 30) for term in go_terms["term"]])
    ax.set_xlabel("Gene Ratio (n_genes in term / n_study genes)", fontsize=18)
    ax.set_ylabel("")
    ax.tick_params(labelsize=17)
    ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
    ax.set_title(f"{condition} {title}", fontsize=20)

    cbar = fig.colorbar(
        color_mapper,
        ax=ax,
        orientation="vertical",
        pad=0.02,
        fraction=0.03,
        aspect=60,
    )
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(which="both", labelsize=16)
    cbar.set_label("FDR (adjusted p)", fontsize=18)

    ax.spines[["top", "right"]].set_visible(False)

    fig.patch.set_facecolor("white")
    plt.tight_layout()

    output_path.mkdir(parents=True, exist_ok=True)
    save_path = output_path / output_filename
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Figure saved: {save_path}")
    plt.close(fig)

    return save_path
