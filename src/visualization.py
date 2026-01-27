"""Genomic visualization and GO enrichment analysis."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import textwrap
from pathlib import Path
from datetime import datetime
from typing import Union
from matplotlib.patches import Patch

from .preprocessing import normalize_rpm, load_tlr_data
from .cache import load_geneid2nt, load_go_dag, load_gene2go
from .config import get_deseq2_dir, CUSTOM_PALETTE_6

sc = None
adjust_text = None
GOEnrichmentStudyNS = None


def get_output_path() -> Path:
    """Get the output directory path based on current date."""
    p = Path().cwd().parent
    date = datetime.today().strftime("%Y%m%d")[2:]
    output_path = p / f"results/{date}"
    output_path.mkdir(parents=True, exist_ok=True)
    return output_path


def save_csv(table: pd.DataFrame, table_name: str, output_path: Path = None) -> None:
    """Save DataFrame to CSV file."""
    if output_path is None:
        output_path = get_output_path()

    path = output_path / f"{table_name}.csv"
    table.to_csv(path)
    print(f"CSV saved: {path}")


def load_scanpy():
    """Lazily load scanpy for PCA plotting."""
    global sc
    if sc is None:
        import scanpy as _sc

        sc = _sc
    return sc


def load_adjust_text():
    """Lazily load adjustText for volcano plot labels."""
    global adjust_text
    if adjust_text is None:
        from adjustText import adjust_text as _adjust_text

        adjust_text = _adjust_text
    return adjust_text


def load_goatools():
    """Lazily load goatools dependencies."""
    global GOEnrichmentStudyNS
    if GOEnrichmentStudyNS is None:
        from goatools.goea.go_enrichment_ns import (
            GOEnrichmentStudyNS as _GOEnrichmentStudyNS,
        )

        GOEnrichmentStudyNS = _GOEnrichmentStudyNS
    return GOEnrichmentStudyNS


AnnData = None


def _plot_tlr_panel(
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
        ax_main: Main axis for the dose-response curve.
        ax_bar: Axis for the Fla-PA bar.
        df: DataFrame with concentration and Average columns.
        conc_col: Name of the concentration column.
        fla_pa_val: Fla-PA average value (None to skip bar).
        xlabel: Label for x-axis.
        title: Plot title.
        label: Legend label for the curve.
        color: Color for the curve.
        xlim: X-axis limits as (min, max).
    """
    assert conc_col in df.columns, f"Column '{conc_col}' not found in DataFrame"
    assert "Average" in df.columns, "Column 'Average' not found in DataFrame"

    # Filter out zero concentrations for log scale
    df_plot = df[df[conc_col] > 0].copy()

    x = df_plot[conc_col].values
    y = df_plot["Average"].values

    # Sort by x for proper line connection
    sort_idx = np.argsort(x)
    x_sorted = x[sort_idx]
    y_sorted = y[sort_idx]

    # Plot dose-response curve
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

    # Fla-PA bar
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

        # Draw dotted line from y-axis across to the bar
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
    """Create separate TLR2/TLR4 dose-response plots with Fla-PA bar.

    Args:
        tlr2_df: DataFrame with TLR2 data (Concentration_ng_mL, Average).
        tlr4_df: DataFrame with TLR4 data (Concentration_EU_mL, Average).
        fla_pa_data: Dictionary with Fla-PA measurements for each TLR.
        output_path: Directory to save the plot. Defaults to results/figures/supplementary.
        output_filename: Name of the output file.

    Returns:
        Path to the saved figure.
    """
    # Create the plot with 2 rows, each with line plot and bar subplot
    fig, axes = plt.subplots(
        2, 2, figsize=(10, 10), gridspec_kw={"width_ratios": [4, 1]}
    )
    ax1, ax1_bar = axes[0]
    ax2, ax2_bar = axes[1]

    # Get Fla-PA values if available
    fla_pa_tlr4 = fla_pa_data["tlr4"]["average"] if fla_pa_data else None
    fla_pa_tlr2 = fla_pa_data["tlr2"]["average"] if fla_pa_data else None

    # TLR4/LPS plot (top)
    _plot_tlr_panel(
        ax_main=ax1,
        ax_bar=ax1_bar,
        df=tlr4_df,
        conc_col="Concentration_EU_mL",
        fla_pa_val=fla_pa_tlr4,
        xlabel="Concentration (EU/ml)",
        title="HEK-Blue™ Reporter Line TLR4 LPS",
        label="LPS",
    )

    # TLR2/Pam3 plot (bottom)
    _plot_tlr_panel(
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

    # Set default output path
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
        gene: Name of the gene to plot (default: "IL6").
        output_filename: Name of output file (default: "{gene}_expression.png").
        palette: Color palette to use (default: CUSTOM_PALETTE_6).
        n_cols: Number of columns in the multipanel layout (default: 3).

    Returns:
        Path: Absolute path to the saved figure.
    """
    # Normalize data
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

        # Sort by gene expression values
        gene_sorted = gene_of_ligand_class.sort_values(by=gene, ascending=False)

        bars = ax.bar(
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

    # Set default filename if not provided
    if output_filename is None:
        output_filename = f"{gene.lower()}_expression.png"

    save_path = output_path / output_filename
    try:
        plt.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"Figure saved to: {save_path.absolute()}")
    finally:
        plt.close(fig)

    return save_path.absolute()


def plot_pca_for_pandas(
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
        X: DataFrame with gene counts (samples as rows, genes as columns).
        labels: DataFrame with 'label' column or numpy array with labels for coloring.
        with_sample_names: If True, creates larger plot with sample name labels.
                          If False, creates smaller plot without labels.
        output_filename: Name of output file. Defaults to 'pca_plot.png'.
        palette: Color palette to use.
        hue_order: Order of classes for legend and color assignment.

    Returns:
        Path to the saved figure.
    """
    from sklearn.decomposition import PCA

    # Extract labels as array
    label_values = (
        labels["label"].to_numpy() if isinstance(labels, pd.DataFrame) else labels
    )

    # Get sample names - handle both DataFrame and numpy array
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
            adjust_text = load_adjust_text()
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


class Plotter:
    geneid_symbol_mapper_human = None
    goeaobj = None
    total_go_terms = None
    is_initialized = False

    def __init__(
        self,
        dds,
        res: pd.DataFrame,
        sigs: pd.DataFrame,
        analysis_name: str,
        figures_dir: Path = None,
    ):
        """Initialize Plotter with DESeq2 results and analysis name.

        Args:
            dds: AnnData object from DESeq2 analysis.
            res: DataFrame with all DESeq2 results.
            sigs: DataFrame containing significant gene signatures.
            analysis_name: Name for this analysis (used in plot titles and filenames).
            figures_dir: Directory for saving figures.
        """
        self.dds = dds
        self.res = res
        self.sigs = sigs
        self.analysis_name = analysis_name
        self.figures_dir = figures_dir or (
            Path.cwd() / "results" / "figures" / "deseq2"
        )
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        self.go_terms_dir = Path.cwd() / "results" / "go_terms"
        self.go_terms_dir.mkdir(parents=True, exist_ok=True)

        if not Plotter.is_initialized:
            load_goatools()
            Plotter._initialise_go()
            Plotter.is_initialized = True

    def make_figure(self, plot_type: str):
        """Create and save a figure based on the specified plot type.

        Args:
            plot_type: Type of plot to create. Options: 'volcano', 'histogram', 'pca', 'go'.

        Raises:
            ValueError: If plot_type is not one of the supported types.
        """
        if plot_type not in ["volcano", "histogram", "pca", "go"]:
            raise ValueError(
                "plot() takes 'volcano', 'histogram', 'pca' or 'go' as a condition."
            )

        fig_extension = "png"
        fname = self.figures_dir / f"{self.analysis_name}_{plot_type}.{fig_extension}"
        fontsize = 12

        if plot_type == "volcano":
            fig = self.plot_volcano()
            plt.title(
                f"{self.analysis_name} Differentially Expressed Genes",
                fontsize=fontsize,
                fontweight="bold",
            )

        if plot_type == "histogram":
            fig = self.plot_historgram()
            fig.ax_col_dendrogram.set_title(
                f"{self.analysis_name} Differentially Expressed Genes",
                fontsize=fontsize,
                fontweight="bold",
                pad=2,
            )
        if plot_type == "pca":
            fig = self.plot_pca()
            fig.axes[0].set_title(
                f"{self.analysis_name}", fontsize=fontsize, fontweight="bold"
            )

        if plot_type == "go":
            fig = self.plot_go()
            plt.title(
                f"{self.analysis_name} Top 20 Significant GO Terms",
                fontsize=fontsize,
                fontweight="bold",
            )

            go_df = self.generate_go_table()
            save_csv(
                go_df, f"{self.analysis_name}_go_terms", output_path=self.go_terms_dir
            )

        try:
            fig.figure.savefig(
                fname,
                format=fig_extension,
                dpi=300,
                bbox_inches="tight",
            )
            print(f"Figure saved: {fname}")
        finally:
            plt.close()

    @classmethod
    def _initialise_go(cls, path_to_supporting_files: Path = None):
        """Initialize GO terms and goatools classes.

        Args:
            path_to_supporting_files: Path to directory containing gene2go and go-basic.obo files.
                                     Defaults to ~/MATseq/data/deseq2 if not provided.
        """
        print("Initializing GO terms...")

        if path_to_supporting_files is None:
            path_to_supporting_files = get_deseq2_dir()

        GOEnrichmentStudyNS = load_goatools()

        # Load gene to GO associations (cached)
        genes = load_gene2go(
            path_to_supporting_files / "gene2go", taxids=[9606], namespaces={"BP"}
        )
        ns2assoc = genes.get_ns2assc()

        # Load GO ontologies (cached)
        obodag = load_go_dag(path_to_supporting_files / "go-basic.obo")

        # Load gene ID mapper (cached)
        GeneID2nt_hs = load_geneid2nt(
            path_to_supporting_files / "gene_result_ncbi_human_proteincoding.txt"
        )

        cls.goeaobj = GOEnrichmentStudyNS(
            GeneID2nt_hs.keys(),  # List of human protein-coding genes
            ns2assoc,  # Geneid/GO associations
            obodag,  # Ontologies
            propagate_counts=False,
            alpha=0.05,  # Default significance cut-off
            methods=["fdr_bh"],  # Default correction method for multiple testing
        )

        cls.geneid_symbol_mapper_human = {
            GeneID2nt_hs[key].Symbol: GeneID2nt_hs[key].GeneID for key in GeneID2nt_hs
        }

    def generate_go_table(self) -> pd.DataFrame:
        """Generate GO enrichment table for significant genes.

        Returns:
            DataFrame with GO enrichment results including terms, p-values, and gene lists.
        """
        sigs_list = [str(gene) for gene in self.sigs.index]
        sigs_ids = [
            int(Plotter.geneid_symbol_mapper_human[gene])
            for gene in sigs_list
            if gene in Plotter.geneid_symbol_mapper_human
        ]
        print(
            f"Mapped {len(sigs_ids)/len(sigs_list)*100:.2f}% of "
            "significantly differentially expressed gene symbols to gene IDs."
        )

        goea_results = Plotter.goeaobj.run_study(sigs_ids, prt=None)
        goea_results_sig = [r for r in goea_results if r.p_fdr_bh < 0.05]

        inverted_mapping_dictionary = {
            v: k for k, v in Plotter.geneid_symbol_mapper_human.items()
        }

        go_df = pd.DataFrame(
            list(
                map(
                    lambda x: [
                        x.GO,
                        x.goterm.name,
                        x.goterm.namespace,
                        x.p_uncorrected,
                        x.p_fdr_bh,
                        x.ratio_in_study[0],
                        x.ratio_in_study[1],
                        x.ratio_in_study[0] / x.ratio_in_study[1],
                        list(
                            map(lambda y: inverted_mapping_dictionary[y], x.study_items)
                        ),
                    ],
                    goea_results_sig,
                )
            ),
            columns=[
                "GO",
                "term",
                "class",
                "raw_pvalue",
                "fdr",
                "n_genes",
                "n_study",
                "ratio_in_study",
                "gene_symbols",
            ],
        )
        go_df["gene_symbols"] = go_df["gene_symbols"].apply(lambda x: ", ".join(x))
        go_df = go_df.sort_values("fdr", ascending=True)
        return go_df

    def plot_go(self):
        """Create horizontal bar plot of top 20 GO enrichment terms.

        Returns:
            Matplotlib figure object.
        """
        if self.sigs is None or len(self.sigs) == 0:
            raise ValueError("No significant genes available for GO enrichment")

        go_terms = self.generate_go_table()

        if go_terms.empty:
            raise ValueError("No GO enrichment terms to plot")

        go_terms = go_terms[:20]
        go_terms = go_terms.sort_values("ratio_in_study", ascending=False)

        if len(go_terms) == 0:
            raise ValueError("No GO terms remaining after filtering")

        norm = mpl.colors.Normalize(vmin=go_terms.fdr.min(), vmax=go_terms.fdr.max())
        color_mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.bwr_r)
        rc = {"axes.edgecolor": "black"}
        sns.set_theme(rc=rc)

        g = plt.figure(figsize=(8, 10))

        ax = sns.barplot(
            data=go_terms,
            x=go_terms["n_genes"] / go_terms["n_study"],
            y="term",
            palette=list(color_mapper.to_rgba(go_terms.fdr.values)),
        )

        ax.set_yticklabels([textwrap.fill(term, 40) for term in go_terms["term"]])
        ax.set_xlabel("Gene Ratio")
        cbar = g.colorbar(
            color_mapper, ax=ax, orientation="vertical", pad=0.01, format="{x:.4f}"
        )
        cbar.ax.set_position([0.8, 0.5, 0.2, 0.3])
        cbar.ax.set_title("padj", loc="left", pad=4.0)
        g.patch.set_facecolor("white")

        return g

    def plot_volcano(self, log2foldchange: float = 2):
        """Create volcano plot showing differentially expressed genes.

        Args:
            log2foldchange: Threshold for log2 fold change to classify genes.

        Returns:
            Matplotlib figure object.
        """
        adjust_text = load_adjust_text()

        # Add a small number to zeros to avoid inf in padj_log
        grapher = self.res.assign(
            padj_log=self.res["padj"].apply(
                lambda x: -np.log10(x) if x != 0 else -np.log10(x + 1e-300)
            ),
            color="no_expression_change",
        )

        grapher.loc[grapher["log2FoldChange"] > log2foldchange, "color"] = (
            "overexpressed"
        )
        grapher.loc[grapher["log2FoldChange"] < -log2foldchange, "color"] = (
            "underexpressed"
        )

        # Subset to only overexpressed and underexpressed
        grapher_subset = grapher[
            grapher["color"].isin(["overexpressed", "underexpressed"])
        ]

        print(f"Number of DE genes: {len(grapher_subset)}")

        # Sort by padj_log (descending) and log2FoldChange (ascending)
        sorted_grapher_padj_log = grapher_subset.sort_values(
            by="padj_log", ascending=False
        )
        sorted_grapher_log2foldchange = grapher_subset.sort_values(
            by="log2FoldChange", ascending=True
        )

        # Select genes for annotation (top by padj and extremes by fold change)
        annotation_subset = pd.concat(
            [
                sorted_grapher_padj_log.head(20),
                sorted_grapher_log2foldchange.head(10),
                sorted_grapher_log2foldchange.tail(10),
            ]
        ).drop_duplicates()

        g = plt.figure(figsize=(8, 10))
        rc = {
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.titlepad": 20,
            "font.size": 10,
            "font.family": "sans-serif",
            "legend.frameon": "False",
            "legend.loc": "upper right",
            "lines.linestyle": "--",
            "lines.linewidth": 1,
            "lines.color": "k",
            "axes.facecolor": "white",
        }
        sns.set_theme(rc=rc)

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

        # Draw threshold lines
        ax.axhline(1.3, zorder=1)
        ax.axvline(2, zorder=1)
        ax.axvline(-2, zorder=1)

        # Add gene name labels
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
        plt.xlabel("$log_{2}$ fold change")
        plt.ylabel("-$log_{10}$ FDR")
        plt.ylim(-2, grapher["padj_log"].max() + 5)

        g.figure.patch.set_facecolor("white")
        return g

    def plot_historgram(self, num_top_sig: Union[int, str] = 50):
        """Create hierarchical clustering heatmap of significant genes.

        Args:
            num_top_sig: Number of top significant genes to plot, or "all".

        Returns:
            Seaborn ClusterGrid object.
        """
        if num_top_sig != "all":
            sigs = self.sigs.sort_values("padj")[:num_top_sig]
        else:
            sigs = self.sigs

        dds_sigs = self.dds[:, sigs.index]
        dds_sigs.layers["log1p"] = np.log1p(dds_sigs.layers["normed_counts"])

        sns.set_theme(rc={"ytick.labelsize": 8})
        grapher = pd.DataFrame(
            dds_sigs.layers["log1p"].T,
            index=dds_sigs.var_names,
            columns=dds_sigs.obs.condition,
        )

        lut = dict(zip(set(dds_sigs.obs.condition), "rgb"))
        col_colors = list(dds_sigs.obs.condition.map(lut))

        plt.figure(figsize=(8, 10), layout="tight")
        ax = sns.clustermap(
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
        ax.ax_heatmap.set_xticklabels([])
        ax.ax_heatmap.set(xlabel=None)
        plt.subplots_adjust(hspace=0.01)
        ax.figure.patch.set_facecolor("white")

        return ax

    def plot_pca(self, w_text: bool = False):
        """Create PCA visualization of samples.

        Args:
            w_text: Whether to add sample name labels to points.

        Returns:
            Matplotlib figure object.
        """
        sc = load_scanpy()

        dds = self.dds.copy()
        sc.tl.pca(dds, n_comps=2)
        pca = dds.obsm["X_pca"]
        sample_names = list(sc.get.obs_df(dds).index)
        ax = sc.pl.pca(
            dds,
            color="condition",
            size=300,
            show=False,
            title=" ",
            return_fig=True,
        )
        ax.set_facecolor("white")

        if w_text:
            for i in range(len(pca)):
                ax.text(
                    pca[i][0],
                    pca[i][1],
                    sample_names[i],
                    ha="left",
                    va="bottom",
                )

        return ax


def make_probability_matrix(
    model_name: str, test_pred: pd.DataFrame, output_path: Path = None
) -> sns.matrix.ClusterGrid:
    """Create and save heatmap of model prediction probabilities.

    Args:
        model_name: Name of the model for the plot title.
        test_pred: DataFrame with prediction probabilities (samples x classes).
        output_path: Directory to save the figure. Defaults to results/figures.

    Returns:
        Seaborn heatmap object.
    """
    if output_path is None:
        project_root = Path(__file__).parent.parent
        output_path = project_root / "results" / "figures"
        output_path.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(12, 8))

    g = sns.heatmap(test_pred, cmap="YlGnBu")
    g.set_title(f"{model_name} predictions", fontsize=15)
    g.set_xlabel("Reference", fontsize=15)
    g.set_ylabel("Sample", fontsize=15, labelpad=10)
    plt.tight_layout()

    save_path = output_path / f"{model_name}_predictions.png"
    try:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Figure saved to: {save_path.absolute()}")
    finally:
        plt.close()

    return g
