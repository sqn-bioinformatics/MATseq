"""Genomic visualization functions."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import textwrap
from pathlib import Path
from typing import Union
from matplotlib.patches import Patch
from sklearn.decomposition import PCA
from adjustText import adjust_text
import scanpy as sc

from .preprocessing import normalize_rpm
from .config import CUSTOM_PALETTE_6


def plot_tlr_panel(
    ax_main,
    ax_bar,
    df: pd.DataFrame,
    conc_col: str,
    fla_pa_val: float = None,
    xlabel: str = "Concentration",
    title: str = "TLR",
    label: str = "Ligand",
    color: str = "#1f77b4",
    xlim: tuple = (0.01, 100),
):
    """Plot a single TLR dose-response panel with optional Fla-PA bar.

    Args:
        ax_main: Main axis for dose-response curve.
        ax_bar: Axis for Fla-PA bar.
        df: DataFrame with concentration and Average columns.
        conc_col: Name of concentration column.
        fla_pa_val: Fla-PA average value (None to skip bar).
        xlabel: X-axis label.
        title: Plot title.
        label: Legend label for curve.
        color: Curve color.
        xlim: X-axis limits as (min, max).
    """
    assert conc_col in df.columns, f"Column '{conc_col}' not found in DataFrame"
    assert "Average" in df.columns, "Column 'Average' not found in DataFrame"

    df_plot = df[df[conc_col] > 0].copy()

    x = df_plot[conc_col].values
    y = df_plot["Average"].values

    sort_idx = np.argsort(x)
    x_sorted = x[sort_idx]
    y_sorted = y[sort_idx]

    ax_main.plot(
        x_sorted, y_sorted, "o-", linewidth=2, markersize=8, color=color, label=label
    )

    ax_main.set_xscale("log")
    ax_main.set_xlabel(xlabel, fontsize=12)
    ax_main.set_ylabel("OD (630 nm)", fontsize=12)
    ax_main.set_title(title, fontsize=14)
    ax_main.grid(True, alpha=0.3, linestyle="-", linewidth=0.5)
    ax_main.legend(fontsize=10)
    ax_main.set_xlim(xlim)

    if fla_pa_val is not None:
        ax_bar.bar(
            ["Fla-PA"],
            [fla_pa_val],
            color="lightblue",
            alpha=0.4,
            width=0.5,
            edgecolor="black",
            linewidth=1.2,
        )
        ax_bar.set_ylim(ax_main.get_ylim())
        ax_bar.set_ylabel("")
        ax_bar.tick_params(left=False, labelleft=False)
        ax_bar.grid(True, alpha=0.3, linestyle="-", linewidth=0.5, axis="y")
        ax_main.axhline(y=fla_pa_val, color="black", linestyle=":", linewidth=1.5)
    else:
        ax_bar.axis("off")


def plot_tlr_hek_blue(
    tlr2_df: pd.DataFrame,
    tlr4_df: pd.DataFrame,
    fla_pa_data: dict = None,
    output_path: Path = None,
    output_filename: str = "TLR_HEK_Blue.png",
) -> Path:
    """Create TLR2/TLR4 dose-response plots with Fla-PA bar.

    Args:
        tlr2_df: DataFrame with TLR2 data (Concentration_ng_mL, Average).
        tlr4_df: DataFrame with TLR4 data (Concentration_EU_mL, Average).
        fla_pa_data: Dictionary with Fla-PA measurements for each TLR.
        output_path: Directory to save plot.
        output_filename: Output file name.

    Returns:
        Path to saved figure.
    """
    fig, axes = plt.subplots(
        2, 2, figsize=(10, 10), gridspec_kw={"width_ratios": [4, 1]}
    )
    ax1, ax1_bar = axes[0]
    ax2, ax2_bar = axes[1]

    fla_pa_tlr4 = fla_pa_data["tlr4"]["average"] if fla_pa_data else None
    fla_pa_tlr2 = fla_pa_data["tlr2"]["average"] if fla_pa_data else None

    plot_tlr_panel(
        ax_main=ax1,
        ax_bar=ax1_bar,
        df=tlr4_df,
        conc_col="Concentration_EU_mL",
        fla_pa_val=fla_pa_tlr4,
        xlabel="Concentration (EU/ml)",
        title="HEK-Blue™ Reporter Line TLR4 LPS",
        label="LPS",
    )

    plot_tlr_panel(
        ax_main=ax2,
        ax_bar=ax2_bar,
        df=tlr2_df,
        conc_col="Concentration_ng_mL",
        fla_pa_val=fla_pa_tlr2,
        xlabel="Concentration (ng/mL)",
        title="HEK-Blue™ Reporter Line TLR2 Pam3",
        label="Pam3",
    )

    plt.tight_layout()

    if output_path is None:
        output_path = (
            Path(__file__).parent.parent / "results" / "figures" / "supplementary"
        )
    output_path.mkdir(parents=True, exist_ok=True)

    save_path = output_path / output_filename
    try:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Plot saved as '{save_path}'")
    finally:
        plt.close()

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


def plot_pca_pandas(
    name: str,
    X: pd.DataFrame,
    labels: Union[pd.DataFrame, np.ndarray],
    with_sample_names: bool = False,
    output_filename: str = None,
    palette: list = CUSTOM_PALETTE_6,
    hue_order: list = None,
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
    marker_size = 200 if with_sample_names else 80

    with plt.rc_context(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    ):
        fig = plt.figure(figsize=figsize)
        ax = sns.scatterplot(
            x=X_reduced[:, 0],
            y=X_reduced[:, 1],
            hue=label_values,
            hue_order=hue_order,
            s=marker_size,
            alpha=0.6,
            palette=palette,
        )

        ax.set_xlabel("1st Eigenvector", fontsize=14)
        ax.set_ylabel("2nd Eigenvector", fontsize=14)
        ax.tick_params(axis="both", labelsize=14)

        if with_sample_names:
            texts = [
                ax.text(
                    X_reduced[:, 0][i],
                    X_reduced[:, 1][i],
                    sample_names[i],
                    ha="left",
                    va="bottom",
                    alpha=0.8,
                    fontsize=12,
                )
                for i in range(len(X_reduced))
            ]
            adjust_text(texts, arrowprops=dict(arrowstyle="->", color="black"))
            ax.legend(bbox_to_anchor=(1, 1.0), ncol=1, fontsize=12)
        else:
            ax.get_legend().remove()

        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
        for spine in ["left", "bottom"]:
            ax.spines[spine].set_linewidth(1.5)

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
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
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
        size="baseMean",
        sizes=(20, 50),
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
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
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

    dds_sigs = dds[:, sigs.index]
    dds_sigs.layers["log1p"] = np.log1p(dds_sigs.layers["normed_counts"])

    grapher = pd.DataFrame(
        dds_sigs.layers["log1p"].T,
        index=dds_sigs.var_names,
        columns=dds_sigs.obs.condition,
    )

    lut = dict(zip(set(dds_sigs.obs.condition), "rgb"))
    col_colors = list(dds_sigs.obs.condition.map(lut))

    g = sns.clustermap(
        figsize=(8, 10),
        data=grapher,
        cmap="RdYlBu_r",
        z_score=None,
        dendrogram_ratio=(0.1, 0.1),
        cbar_pos=(0.93, 0.2, 0.03, 0.45),
        cbar_kws=dict(
            location="left",
            orientation="vertical",
            pad=2,
        ),
        col_colors=col_colors,
    )

    handles = [Patch(facecolor=lut[name]) for name in lut]
    plt.legend(
        handles,
        lut,
        bbox_transform=plt.gcf().transFigure,
        bbox_to_anchor=(1, 1),
        loc="upper right",
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
    plt.subplots_adjust(hspace=0.01)
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
    ax.set_xlabel("Gene Ratio", fontsize=12)
    ax.set_ylabel("")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.grid(True, alpha=0.3, linestyle="-", linewidth=0.5, axis="x")

    cbar = fig.colorbar(
        color_mapper, ax=ax, orientation="vertical", pad=0.01, format=mpl.ticker.LogFormatterSciNotation()
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
