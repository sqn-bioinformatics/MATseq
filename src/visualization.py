"""Genomic visualization functions."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import textwrap
from pathlib import Path
from typing import Union
from matplotlib.patches import Patch, FancyArrowPatch
from matplotlib_venn import venn2, venn3
from sklearn.decomposition import PCA
from adjustText import adjust_text
import scanpy as sc

from .preprocessing import normalize_rpm
from .config import (
    CUSTOM_PALETTE_6,
    display_labels,
    subset_display,
    confusion_title,
)


def _savefig_with_arrow_fallback(save_path, **kwargs):
    """Save the current figure, retrying without adjust_text connector arrows.

    adjustText 0.8 with matplotlib >= 3.x can emit a degenerate FancyArrowPatch
    (text left at its anchor, so posA ~= posB) whose path cannot be clipped,
    raising StopIteration inside the draw call at savefig time. The connectors
    are cosmetic, so on that failure we drop them and re-save.
    """
    try:
        plt.savefig(save_path, **kwargs)
    except StopIteration:
        for ax in plt.gcf().axes:
            for patch in [p for p in ax.patches if isinstance(p, FancyArrowPatch)]:
                patch.remove()
            for txt in list(ax.texts):
                if getattr(txt, "arrow_patch", None) is None:
                    continue
                ax.text(
                    *txt.get_position(),
                    txt.get_text(),
                    ha=txt.get_ha(),
                    va=txt.get_va(),
                    size=txt.get_size(),
                    weight=txt.get_weight(),
                    color=txt.get_color(),
                    alpha=txt.get_alpha(),
                )
                txt.remove()
        plt.savefig(save_path, **kwargs)


def draw_confusion_matrix(ax, cm_norm, class_names, title=None, show_cbar=True,
                          show_tick_labels=True):
    """Render a row-normalized confusion matrix on ``ax`` in Nature 'Blues' style.

    Parameters
    ----------
    ax : matplotlib Axes
    cm_norm : 2D array-like, row-normalized (values in [0, 1])
    class_names : sequence of RAW class labels (display-cleaned internally)
    title : optional str (used verbatim; caller builds via confusion_title)
    show_cbar : whether to draw the colour bar (off for collage panels)
    show_tick_labels : whether to print the per-class ligand names on the axes.
        Off for collage panels, where the panel title already identifies the
        matrix and dropping the (repeated) labels tightens inter-panel spacing.
    """
    cm = np.asarray(cm_norm, dtype=float)
    labels = display_labels(class_names)
    n = cm.shape[0]
    ncol = cm.shape[1]

    # Fontsizes shrink as the matrix widens so "100.0%" always fits a cell.
    ncell = max(n, ncol)
    annot_fs = 9 if ncell <= 6 else (7 if ncell == 7 else 6)
    tick_fs = 9 if ncell <= 7 else 8

    # aspect="auto" + a fixed square box makes every matrix render at the SAME
    # physical size regardless of class count (5x5 / 6x6 / 7x7). With equal
    # aspect the larger matrices drew bigger and tight_layout left uneven
    # margins between collage panels; a fixed box keeps the panels uniform and
    # lets them pack tightly.
    im = ax.imshow(cm, cmap="Blues", vmin=0.0, vmax=1.0, aspect="auto")
    ax.set_box_aspect(1)

    ax.set_xticks(range(ncol))
    ax.set_yticks(range(n))
    if show_tick_labels:
        ax.set_xticklabels(labels[:ncol], rotation=45, ha="right",
                           fontsize=tick_fs)
        ax.set_yticklabels(labels[:n], fontsize=tick_fs)
    else:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    ax.set_xlabel("Predicted class", fontsize=10)
    ax.set_ylabel("True class", fontsize=10)
    if title:
        ax.set_title(title, fontsize=11, pad=8)

    for i in range(n):
        for j in range(ncol):
            v = cm[i, j]
            # Drop the decimal for a full 100% so it never overruns the cell.
            txt = "100%" if v >= 0.9995 else f"{v * 100:.1f}%"
            ax.text(j, i, txt, ha="center", va="center",
                    color="white" if v > 0.5 else "#222222", fontsize=annot_fs)

    # White minor gridlines between cells; hide spines and tick marks.
    ax.set_xticks(np.arange(-0.5, cm.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    if show_cbar:
        cbar = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.outline.set_visible(False)
        cbar.ax.tick_params(length=0)
    return im


def save_confusion_matrix(cm_norm, class_names, model, subset, output_dir,
                          filename=None):
    """Save a single Nature-style confusion matrix figure and return its path."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.5, 6))
    draw_confusion_matrix(ax, cm_norm, class_names,
                          title=confusion_title(model, subset))
    fig.tight_layout()
    if filename is None:
        filename = f"{subset}_{model}_confusion_matrix.png"
    save_path = output_dir / filename
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return save_path


def assemble_panel_collage(panels, output_path, ncols=3, panel_size=5.0,
                           title=None):
    """Compose an N-panel collage from a list of draw callables.

    Each entry in ``panels`` is a callable ``fn(ax)`` that renders one panel on
    the Axes it is given. Panels are laid out row-major into ceil(N/ncols) rows.
    """
    n = len(panels)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(panel_size * ncols, panel_size * nrows)
    )
    axes = np.atleast_1d(axes).ravel()
    for ax, draw in zip(axes, panels):
        draw(ax)
    for ax in axes[n:]:
        ax.axis("off")
    if title:
        fig.suptitle(title, fontsize=15, y=1.0)
    # Panels are fixed square boxes. The gutter is widened (wspace) so the PCA
    # panel's legend, anchored just outside its top-right corner, clears the
    # neighbouring confusion matrix without truncating any class name.
    fig.subplots_adjust(wspace=0.45, hspace=0.18,
                        left=0.04, right=0.98, top=0.93, bottom=0.05)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_venn(
    sets: list,
    set_labels: tuple,
    output_filename: str,
    output_dir: Path = None,
    title: str = None,
) -> Path:
    """Plot a 2- or 3-set Venn diagram and save it.

    Args:
        sets: List of 2 or 3 sets to compare.
        set_labels: Labels for the sets, same length as `sets`.
        output_filename: Name for the output figure.
        output_dir: Directory for the figure. Defaults to results/figures/venn.
        title: Optional figure title.
    """
    if len(sets) == 2:
        venn = venn2
    elif len(sets) == 3:
        venn = venn3
    else:
        raise ValueError(f"plot_venn supports 2 or 3 sets, got {len(sets)}")

    if output_dir is None:
        output_dir = Path.cwd() / "results" / "figures" / "venn"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(8, 8))
    venn([set(s) for s in sets], set_labels=set_labels)
    if title:
        plt.title(title)
    output_path = output_dir / output_filename
    try:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"Saved Venn diagram to {output_path}")
    finally:
        plt.close()

    return output_path


def plot_feature_count_analysis(
    result: dict,
    output_path: Path = None,
    output_filename: str = None,
) -> Path:
    """Plot the sorted mutual-information curve with per-seed and mean elbows."""
    mi_elbow = result["mi_elbow"]
    per_run_elbows = result["per_run_elbows"]

    fig, ax1 = plt.subplots(figsize=(8, 5))
    scores = result["scores"]
    ax1.plot(scores["rank"], scores["mi_sorted"], label="mutual information")
    for e in per_run_elbows:
        ax1.axvline(e, color="C1", ls=":", alpha=0.6)
    ax1.axvline(
        mi_elbow, color="C0", ls="--", label=f"MI elbow (mean {mi_elbow})"
    )
    ax1.set_xlabel("gene rank")
    ax1.set_ylabel("sorted mutual information")
    ax1.set_title(f"MI elbow (per-seed: {per_run_elbows})")
    ax1.legend()

    if output_path is None:
        output_path = (
            Path(__file__).parent.parent / "results" / "figures" / "feature_selection"
        )
    output_path.mkdir(parents=True, exist_ok=True)

    if output_filename is None:
        output_filename = "feature_count_analysis.png"

    save_path = output_path / output_filename
    try:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved: {save_path}")
    finally:
        plt.close(fig)

    return save_path


def plot_gene_expression_by_class(
    data: pd.DataFrame,
    gene: str = "IL6",
    output_filename: str = None,
    palette: list = CUSTOM_PALETTE_6,
    n_cols: int = 3,
) -> Path:
    """Create multipanel bar plots showing gene expression across ligand classes.

    Args:
        data: DataFrame with gene counts (samples as rows, genes as columns).
        gene: Gene name to plot.
        output_filename: Output file name.
        palette: Color palette.
        n_cols: Number of columns in layout.

    Returns:
        Path to saved figure.
    """
    normalized_df = normalize_rpm(data.copy())

    if gene not in normalized_df.columns:
        raise ValueError(f"Gene '{gene}' not found in data columns")

    gene_df = normalized_df[gene].to_frame()
    gene_df.columns = [gene]

    ligand_classes = sorted(set("_" + gene_df.index.str.split("_").str[2] + "_"))
    n_classes = len(ligand_classes)
    n_rows = (n_classes + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 6 * n_rows))
    axes = axes.flatten() if n_classes > 1 else [axes]

    for idx, ligand_class in enumerate(ligand_classes):
        ax = axes[idx]
        gene_of_ligand_class = gene_df[gene_df.index.str.contains(ligand_class)]
        gene_sorted = gene_of_ligand_class.sort_values(by=gene, ascending=False)

        ax.bar(
            range(len(gene_sorted)),
            gene_sorted[gene],
            color=palette[idx % len(palette)],
            alpha=0.8,
            edgecolor="black",
            linewidth=1.2,
        )

        ax.set_xticks(range(len(gene_sorted)))
        ax.set_xticklabels(gene_sorted.index, rotation=45, ha="right", fontsize=9)
        ax.set_ylabel(f"Normalized {gene} counts (RPM)", fontsize=11, fontweight="bold")
        ax.set_xlabel("Samples", fontsize=11, fontweight="bold")
        ax.set_title(
            f'{gene} expression in {ligand_class.strip("_")} samples (n={len(gene_sorted)})',
            fontsize=12,
            fontweight="bold",
            pad=10,
        )
        ax.yaxis.grid(True, alpha=0.3, linestyle="--", linewidth=0.8)
        ax.set_axisbelow(True)

        mean_val = gene_sorted[gene].mean()
        ax.axhline(
            mean_val,
            color="red",
            linestyle="--",
            linewidth=2,
            alpha=0.7,
            label=f"Mean: {mean_val:.1f}",
        )
        ax.legend(loc="upper right", fontsize=9)

        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
        for spine in ["left", "bottom"]:
            ax.spines[spine].set_linewidth(1.5)

    for idx in range(n_classes, len(axes)):
        fig.delaxes(axes[idx])

    plt.tight_layout()

    project_root = Path(__file__).parent.parent
    output_path = project_root / "results" / "figures" / "supplementary"
    output_path.mkdir(parents=True, exist_ok=True)

    if output_filename is None:
        output_filename = f"{gene.lower()}_expression.png"

    save_path = output_path / output_filename
    try:
        plt.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"Figure saved to: {save_path.absolute()}")
    finally:
        plt.close(fig)

    return save_path.absolute()


def draw_pca(ax, X_reduced, label_values, palette=CUSTOM_PALETTE_6,
             hue_order=None, with_sample_names=False, sample_names=None,
             equal_aspect=True, title=None, marker_size=80):
    """Render a 2-component PCA scatter on ``ax`` (equal aspect by default).

    Legend entries are display-cleaned (no underscores). ``X_reduced`` is the
    Nx2 PCA-projected array; ``label_values`` the per-sample class labels.
    """
    sns.scatterplot(
        x=X_reduced[:, 0],
        y=X_reduced[:, 1],
        hue=label_values,
        hue_order=hue_order,
        s=200 if with_sample_names else marker_size,
        alpha=0.6,
        palette=palette,
        ax=ax,
    )
    ax.set_xlabel("PC1", fontsize=12)
    ax.set_ylabel("PC2", fontsize=12)
    ax.tick_params(axis="both", labelsize=11)
    if title:
        ax.set_title(title, fontsize=11, pad=8)

    if with_sample_names and sample_names is not None:
        texts = [
            ax.text(X_reduced[i, 0], X_reduced[i, 1], sample_names[i],
                    ha="left", va="bottom", alpha=0.8, fontsize=12)
            for i in range(len(X_reduced))
        ]
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="->", color="black"))

    # Clean legend labels (strip underscores / apply display names).
    handles, labels_txt = ax.get_legend_handles_labels()
    if handles:
        # Place the legend OUTSIDE the axes, anchored to the top-right corner
        # of the (fixed-square) PCA box. The collage reserves gutter space to
        # the right of the PCA panel (see assemble_panel_collage wspace) so the
        # legend clears the neighbouring confusion matrix and no class name is
        # truncated, while never overlapping the scatter points.
        ax.legend(handles, display_labels(labels_txt),
                  loc="upper left", bbox_to_anchor=(1.02, 1.0),
                  borderaxespad=0, ncol=1, fontsize=9,
                  frameon=False, handletextpad=0.4)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_linewidth(1.5)

    if equal_aspect:
        # Equal data aspect inside a fixed square box: keeps PC1/PC2 on the
        # same scale while matching the fixed-square confusion-matrix panels so
        # all collage panels align and pack tightly.
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_box_aspect(1)
    return ax


def plot_pca_pandas(
    name: str,
    X: pd.DataFrame,
    labels: Union[pd.DataFrame, np.ndarray],
    with_sample_names: bool = False,
    output_filename: str = None,
    palette: list = CUSTOM_PALETTE_6,
    hue_order: list = None,
    equal_aspect: bool = False,
) -> Path:
    """Create PCA visualization for pandas DataFrame data.

    Args:
        name: Name for output file.
        X: DataFrame with features (samples as rows).
        labels: DataFrame with 'label' column or array with labels.
        with_sample_names: If True, add sample name labels.
        output_filename: Output file name.
        palette: Color palette.
        hue_order: Order of classes for legend.
        equal_aspect: If True, equal x/y scaling (no axis stretching).

    Returns:
        Path to saved figure.
    """
    label_values = (
        labels["label"].to_numpy() if isinstance(labels, pd.DataFrame) else labels
    )

    if hasattr(X, "index"):
        sample_names = X.index.to_numpy()
    elif isinstance(labels, pd.DataFrame):
        sample_names = labels.index.to_numpy()
    else:
        sample_names = np.arange(len(X))

    X_reduced = PCA(n_components=2).fit_transform(X)

    figsize = (20, 15) if with_sample_names else (6, 6)

    with plt.rc_context(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    ):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        draw_pca(
            ax,
            X_reduced,
            label_values,
            palette=palette,
            hue_order=hue_order,
            with_sample_names=with_sample_names,
            sample_names=sample_names,
            equal_aspect=equal_aspect,
        )
        plt.tight_layout()

        project_root = Path(__file__).parent.parent
        output_path = project_root / "results" / "figures" / "pca"
        output_path.mkdir(parents=True, exist_ok=True)

        if output_filename is None:
            output_filename = (
                f"{name}_pca_labeled.png" if with_sample_names else f"{name}_pca.png"
            )

        save_path = output_path / output_filename
        try:
            _savefig_with_arrow_fallback(save_path, dpi=300, bbox_inches="tight")
            print(f"Figure saved to: {save_path.absolute()}")
        finally:
            plt.close()

    return save_path.absolute()


def plot_volcano(
    res: pd.DataFrame,
    analysis_name: str,
    log2foldchange: float = 2.0,
    output_path: Path = None,
    output_filename: str = None,
) -> Path:
    """Create volcano plot showing differentially expressed genes.

    Args:
        res: DataFrame with DESeq2 results.
        analysis_name: Name for title and output file.
        log2foldchange: Threshold for log2 fold change.
        output_path: Directory to save plot.
        output_filename: Output file name.

    Returns:
        Path to saved figure.
    """
    grapher = res.assign(
        padj_log=res["padj"].apply(
            lambda x: -np.log10(x) if x != 0 else -np.log10(x + 1e-300)
        ),
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
            sorted_grapher_padj_log.head(20),
            sorted_grapher_log2foldchange.head(10),
            sorted_grapher_log2foldchange.tail(10),
        ]
    ).drop_duplicates()

    fig = plt.figure(figsize=(8, 10))
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
    ax.axvline(2, color="black", linestyle="--", linewidth=1)
    ax.axvline(-2, color="black", linestyle="--", linewidth=1)

    texts = []
    for i, row in annotation_subset.iterrows():
        texts.append(
            plt.text(
                x=row.log2FoldChange,
                y=row.padj_log,
                s=row.name,
                weight="bold",
                size=8,
            )
        )

    adjust_text(texts, arrowprops=dict(arrowstyle="-", color="k"))
    plt.legend(bbox_to_anchor=(1.4, 1), prop={"size": 10, "weight": "bold"})
    plt.xticks(size=10, weight="bold")
    plt.yticks(size=10, weight="bold")
    plt.xlabel("$log_{2}$ fold change", fontsize=12)
    plt.ylabel("-$log_{10}$ FDR", fontsize=12)
    plt.ylim(-2, grapher["padj_log"].max() + 5)
    plt.title(
        f"{analysis_name} Differentially Expressed Genes",
        fontsize=12,
        fontweight="bold",
    )

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_linewidth(1.5)

    fig.patch.set_facecolor("white")

    if output_path is None:
        output_path = Path(__file__).parent.parent / "results" / "figures" / "deseq2"
    output_path.mkdir(parents=True, exist_ok=True)

    if output_filename is None:
        output_filename = f"{analysis_name}_volcano.png"

    save_path = output_path / output_filename
    try:
        _savefig_with_arrow_fallback(save_path, dpi=300, bbox_inches="tight")
        print(f"Figure saved: {save_path}")
    finally:
        plt.close()

    return save_path


def plot_heatmap(
    dds,
    sigs: pd.DataFrame,
    analysis_name: str,
    num_top_sig: Union[int, str] = 50,
    output_path: Path = None,
    output_filename: str = None,
) -> Path:
    """Create hierarchical clustering heatmap of significant genes.

    Args:
        dds: DESeq2 AnnData object.
        sigs: DataFrame with significant genes.
        analysis_name: Name for title and output file.
        num_top_sig: Number of top genes to plot, or "all".
        output_path: Directory to save plot.
        output_filename: Output file name.

    Returns:
        Path to saved figure.
    """
    if num_top_sig != "all":
        sigs = sigs.sort_values("padj")[:num_top_sig]

    dds_sigs = dds[:, sigs.index].copy()
    dds_sigs.layers["log1p"] = np.log1p(dds_sigs.layers["normed_counts"])

    grapher = pd.DataFrame(
        dds_sigs.layers["log1p"].T,
        index=dds_sigs.var_names,
        columns=dds_sigs.obs.condition,
    )

    def _is_negative_control(cond):
        c = str(cond).lower().replace("-", "_")
        return "negative" in c or c == "control"

    other_palette = ["green", "tab:blue", "tab:orange", "tab:purple", "tab:brown"]
    lut, oi = {}, 0
    for cond in dict.fromkeys(dds_sigs.obs.condition):
        if _is_negative_control(cond):
            lut[cond] = "magenta"
        else:
            lut[cond] = other_palette[oi % len(other_palette)]
            oi += 1
    col_colors = list(dds_sigs.obs.condition.map(lut))
    g = sns.clustermap(
        figsize=(8, 10),
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
    g.ax_cbar.set_yticklabels([f"{vmin:.1f}", "0", f"{vmax:.1f}"])

    handles = [Patch(facecolor=lut[name], label=name) for name in lut]

    g.ax_col_dendrogram.legend(
        handles=handles,
        bbox_to_anchor=(1.0, 1.0),
        loc="upper left",
        fontsize=10,
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
    g.ax_heatmap.set_yticklabels(reordered_genes, fontsize=8)

    g.ax_col_dendrogram.set_title(
        f"{analysis_name} Differentially Expressed Genes",
        fontsize=12,
        fontweight="bold",
        pad=2,
    )

    g.figure.subplots_adjust(hspace=0.01, right=0.82)
    g.figure.patch.set_facecolor("white")
    if output_path is None:
        output_path = Path(__file__).parent.parent / "results" / "figures" / "deseq2"
    output_path.mkdir(parents=True, exist_ok=True)

    if output_filename is None:
        output_filename = f"{analysis_name}_histogram.png"

    save_path = output_path / output_filename
    try:
        g.figure.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Figure saved: {save_path}")
    finally:
        plt.close()

    return save_path


def plot_pca_deseq2(
    dds,
    analysis_name: str,
    with_text: bool = False,
    output_path: Path = None,
    output_filename: str = None,
) -> Path:
    """Create PCA visualization from DESeq2 results.

    Args:
        dds: DESeq2 AnnData object.
        analysis_name: Name for title and output file.
        with_text: If True, add sample name labels.
        output_path: Directory to save plot.
        output_filename: Output file name.

    Returns:
        Path to saved figure.
    """
    dds_copy = dds.copy()
    sc.tl.pca(dds_copy, n_comps=2)
    pca = dds_copy.obsm["X_pca"]
    sample_names = list(sc.get.obs_df(dds_copy).index)

    fig = sc.pl.pca(
        dds_copy,
        color="condition",
        size=300,
        show=False,
        title=analysis_name,
        return_fig=True,
    )
    fig.set_facecolor("white")

    if with_text:
        ax = fig.axes[0]
        for i in range(len(pca)):
            ax.text(
                pca[i][0],
                pca[i][1],
                sample_names[i],
                ha="left",
                va="bottom",
            )

    if output_path is None:
        output_path = Path(__file__).parent.parent / "results" / "figures" / "deseq2"
    output_path.mkdir(parents=True, exist_ok=True)

    if output_filename is None:
        output_filename = f"{analysis_name}_pca.png"

    save_path = output_path / output_filename
    try:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Figure saved: {save_path}")
    finally:
        plt.close()

    return save_path


def plot_go(
    go_df: pd.DataFrame,
    title: str = "Top 20 Significant GO Terms",
    output_path: Path = None,
    output_filename: str = None,
) -> Path:
    """Create horizontal bar plot of GO enrichment terms.

    Args:
        go_df: DataFrame with GO enrichment results (from generate_go_table or CSV).
        title: Plot title.
        output_path: Directory to save plot.
        output_filename: Output file name.

    Returns:
        Path to saved figure.
    """
    if go_df.empty:
        raise ValueError("No GO enrichment terms to plot")

    go_terms = go_df.head(20).copy()
    go_terms = go_terms.sort_values("ratio_in_study", ascending=False)

    if len(go_terms) == 0:
        raise ValueError("No GO terms remaining after filtering")

    norm = mpl.colors.Normalize(vmin=go_terms.fdr.min(), vmax=go_terms.fdr.max())
    color_mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.bwr_r)

    fig = plt.figure(figsize=(8, 10))

    ax = sns.barplot(
        data=go_terms,
        x=go_terms["n_genes"] / go_terms["n_study"],
        y="term",
        palette=list(color_mapper.to_rgba(go_terms.fdr.values)),
    )

    ax.set_yticklabels([textwrap.fill(term, 40) for term in go_terms["term"]])
    ax.set_xlabel("Gene Ratio (n_genes in term / n_study genes)", fontsize=12)
    ax.set_ylabel("")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.grid(True, alpha=0.3, linestyle="-", linewidth=0.5, axis="x")

    cbar = fig.colorbar(
        color_mapper,
        ax=ax,
        orientation="vertical",
        pad=0.01,
        format=mpl.ticker.LogFormatterSciNotation(),
    )
    cbar.ax.set_position([0.8, 0.5, 0.2, 0.3])
    cbar.ax.set_title("padj", loc="left", pad=4.0)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_linewidth(1.5)

    fig.patch.set_facecolor("white")
    plt.tight_layout()

    if output_path is None:
        output_path = Path(__file__).parent.parent / "results" / "figures" / "go"
    output_path.mkdir(parents=True, exist_ok=True)

    if output_filename is None:
        output_filename = "go_enrichment.png"

    save_path = output_path / output_filename
    try:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Figure saved: {save_path}")
    finally:
        plt.close()

    return save_path
