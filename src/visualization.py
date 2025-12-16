"""Genomic visualization and GO enrichment analysis."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import textwrap
from pathlib import Path

from .utils import save_csv, get_output_path, CUSTOM_PALETTE_6
from .preprocessing import normalize_rpm

# Lazy imports for goatools (only loaded when needed)
GODag = None
Gene2GoReader = None
GOEnrichmentStudyNS = None


def _load_goatools():
    """Lazily load goatools dependencies."""
    global GODag, Gene2GoReader, GOEnrichmentStudyNS
    if GODag is None:
        from goatools.obo_parser import GODag as _GODag
        from goatools.anno.genetogo_reader import Gene2GoReader as _Gene2GoReader
        from goatools.goea.go_enrichment_ns import (
            GOEnrichmentStudyNS as _GOEnrichmentStudyNS,
        )

        GODag = _GODag
        Gene2GoReader = _Gene2GoReader
        GOEnrichmentStudyNS = _GOEnrichmentStudyNS


def load_tlr_data(data_dir: Path = None) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Load TLR2 (Pam3) and TLR4 (LPS) data from supplementary tables.

    Args:
        data_dir: Path to the supplementary_data directory.
                  Defaults to results/supplementary_data relative to project root.

    Returns:
        Tuple of (tlr2_df, tlr4_df, fla_pa_data) where fla_pa_data contains
        Fla-PA measurements for each TLR.
    """
    if data_dir is None:
        data_dir = Path(__file__).parent.parent / "results" / "supplementary_data"

    # Load TLR4/LPS data (Supplementary Table 5)
    tlr4_raw = pd.read_csv(data_dir / "Supplementary_Table_5.csv")

    # Filter rows with LPS data (non-NaN in LPS columns)
    tlr4_lps_mask = tlr4_raw["OD630nm_LPS_Replicate1"].notna()
    tlr4_lps = tlr4_raw[tlr4_lps_mask]
    tlr4_df = pd.DataFrame(
        {
            "Concentration_EU_mL": tlr4_lps["Concentration_(EU_mL)"],
            "Average": tlr4_lps[
                ["OD630nm_LPS_Replicate1", "OD630nm_LPS_Replicate2"]
            ].mean(axis=1),
        }
    )

    # Extract Fla-PA data for TLR4 (row with Fla-PA values)
    tlr4_fla_mask = tlr4_raw["OD630nm_Fla-PA_Replicate1"].notna()
    tlr4_fla = tlr4_raw[tlr4_fla_mask].iloc[0]

    # Load TLR2/Pam3Csk4 data (Supplementary Table 6)
    tlr2_raw = pd.read_csv(data_dir / "Supplementary_Table_6.csv")

    # Filter rows with Pam3 data (non-NaN in Pam3 columns)
    tlr2_pam_mask = tlr2_raw["OD630nm_Pam3_Replicate1"].notna()
    tlr2_pam = tlr2_raw[tlr2_pam_mask]
    tlr2_df = pd.DataFrame(
        {
            "Concentration_ng_mL": tlr2_pam["Concentration_(ng_mL)"],
            "Average": tlr2_pam[
                ["OD630nm_Pam3_Replicate1", "OD630nm_Pam3_Replicate2"]
            ].mean(axis=1),
        }
    )

    # Extract Fla-PA data for TLR2 (row with Fla-PA values)
    tlr2_fla_mask = tlr2_raw["OD630nm_Fla-PA_Replicate1"].notna()
    tlr2_fla = tlr2_raw[tlr2_fla_mask].iloc[0]

    # Build Fla-PA data dictionary
    fla_pa_data = {
        "tlr4": {
            "concentration": tlr4_fla["Concentration_(EU_mL)"],
            "average": np.mean(
                [
                    tlr4_fla["OD630nm_Fla-PA_Replicate1"],
                    tlr4_fla["OD630nm_Fla-PA_Replicate2"],
                ]
            ),
        },
        "tlr2": {
            "concentration": tlr2_fla["Concentration_(ng_mL)"],
            "average": np.mean(
                [
                    tlr2_fla["OD630nm_Fla-PA_Replicate1"],
                    tlr2_fla["OD630nm_Fla-PA_Replicate2"],
                ]
            ),
        },
    }

    return tlr2_df, tlr4_df, fla_pa_data


def plot_tlr_hek_blue(
    tlr2_df: pd.DataFrame,
    tlr4_df: pd.DataFrame,
    fla_pa_data: dict = None,
    output_path: Path = None,
    output_filename: str = "TLR_HEK_Blue.png",
) -> Path:
    """Create separate TLR2/TLR4 dose-response plots.

    Args:
        tlr2_df: DataFrame with TLR2 data (Concentration_ng_mL, Average, StdDev).
        tlr4_df: DataFrame with TLR4 data (Concentration_EU_mL, Average, StdDev).
        fla_pa_data: Dictionary with Fla-PA measurements for each TLR.
        output_path: Directory to save the plot. Defaults to results/figures/supplementary.
        output_filename: Name of the output file.

    Returns:
        Path to the saved figure.
    """
    # Remove the 0 concentration rows for log scale plotting
    tlr2_plot = tlr2_df[tlr2_df["Concentration_ng_mL"] > 0].copy()
    tlr4_plot = tlr4_df[tlr4_df["Concentration_EU_mL"] > 0].copy()

    # Create the plot with 2 subplots (ax1=top, ax2=bottom)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))

    # TLR4/LPS plot (top)
    ax1.errorbar(
        tlr4_plot["Concentration_EU_mL"],
        tlr4_plot["Average"],
        fmt="o-",
        capsize=5,
        capthick=2,
        linewidth=2,
        markersize=8,
        color="#1f77b4",
        label="LPS",
    )

    # Add Fla-PA point to TLR4 plot
    if fla_pa_data is not None:
        ax1.errorbar(
            fla_pa_data["tlr4"]["concentration"],
            fla_pa_data["tlr4"]["average"],
            fmt="o",
            capsize=5,
            capthick=2,
            markersize=10,
            color="magenta",
            label="Fla-PA",
        )

    ax1.set_xscale("log")
    ax1.set_xlabel("Concentration (EU/ml)", fontsize=12)
    ax1.set_ylabel("OD (630 nm)", fontsize=12)
    ax1.set_title("HEK-Blue™ Reporter Line TLR4 LPS", fontsize=14)
    ax1.grid(True, alpha=0.3, linestyle="-", linewidth=0.5)
    ax1.legend(fontsize=10)
    ax1.set_xlim(0.01, 100)

    # TLR2/Pam3 plot (bottom)
    ax2.errorbar(
        tlr2_plot["Concentration_ng_mL"],
        tlr2_plot["Average"],
        fmt="o-",
        capsize=5,
        capthick=2,
        linewidth=2,
        markersize=8,
        color="#1f77b4",
        label="Pam3",
    )

    # Add Fla-PA point to TLR2 plot
    if fla_pa_data is not None:
        ax2.errorbar(
            fla_pa_data["tlr2"]["concentration"],
            fla_pa_data["tlr2"]["average"],
            fmt="o",
            capsize=5,
            capthick=2,
            markersize=10,
            color="magenta",
            label="Fla-PA",
        )

    ax2.set_xscale("log")
    ax2.set_xlabel("Concentration (ng/mL)", fontsize=12)
    ax2.set_ylabel("OD (630 nm)", fontsize=12)
    ax2.set_title("HEK-Blue™ Reporter Line TLR2 Pam3", fontsize=14)
    ax2.grid(True, alpha=0.3, linestyle="-", linewidth=0.5)
    ax2.legend(fontsize=10)
    ax2.set_xlim(0.01, 100)

    plt.tight_layout()

    # Set default output path
    if output_path is None:
        output_path = (
            Path(__file__).parent.parent / "results" / "figures" / "supplementary"
        )
    output_path.mkdir(parents=True, exist_ok=True)

    save_path = output_path / output_filename
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Plot saved as '{save_path}'")
    plt.show()

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

    # Set up figure
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
    plt.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Figure saved to: {save_path.absolute()}")
    plt.show()
    # plt.close(fig)

    return save_path.absolute()


class Plotter:
    """Comprehensive visualization and analysis class for genomic data.

    This class provides methods for creating various plots and performing GO enrichment analysis.
    It maintains class-level GO term data that is initialized once and shared across instances.

    Attributes:
        geneid_symbol_mapper_human: Mapping from gene symbols to gene IDs.
        goeaobj: GO enrichment analysis object.
        is_initialized: Flag indicating if GO terms have been initialized.
    """

    geneid_symbol_mapper_human = None
    goeaobj = None
    total_go_terms = None
    is_initialized = False

    def __init__(
        self,
        sigs: pd.DataFrame,
        analysis_name: str,
    ):
        """Initialize Plotter with gene signatures and analysis name.

        Args:
            sigs: DataFrame containing gene signatures.
            analysis_name: Name for this analysis (used in plot titles and filenames).
        """
        self.sigs = sigs
        self.analysis_name = analysis_name

        if not Plotter.is_initialized:
            _load_goatools()  # Load goatools dependencies
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

        output_path = get_output_path()
        fig_extension = "png"
        fname = output_path / f"{self.analysis_name}_{plot_type}.{fig_extension}"
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
            save_csv(go_df, f"{self.analysis_name}_go_terms")

        fig.figure.savefig(
            fname,
            format=fig_extension,
            dpi=300,
            bbox_inches="tight",
        )

    @classmethod
    def _initialise_go(cls, path_to_supporting_files: Path = None):
        """Initialize GO terms and goatools classes.

        Args:
            path_to_supporting_files: Path to directory containing gene2go and go-basic.obo files.
                                     Defaults to notebooks/DeSeq2/deseq2 if not provided.
        """
        print("Initializing GO terms...")

        if path_to_supporting_files is None:
            path_to_supporting_files = Path().cwd() / "notebooks/DeSeq2/deseq2"

        # Load gene to GO associations
        genes = Gene2GoReader(
            path_to_supporting_files / "gene2go", taxids=[9606], namespaces={"BP"}
        )
        ns2assoc = genes.get_ns2assc()

        # Load GO ontologies
        obodag = GODag(path_to_supporting_files / "go-basic.obo")

        # Import gene ID mapper (requires DeSeq2 module)
        # This assumes the DeSeq2 module is in the path
        try:
            from genes_ncbi_homo_sapiens_proteincoding import GENEID2NT as GeneID2nt_hs
        except ImportError:
            import sys

            sys.path.insert(1, str(path_to_supporting_files.parent))
            from genes_ncbi_homo_sapiens_proteincoding import GENEID2NT as GeneID2nt_hs

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
        sigs_ids = [
            Plotter.geneid_symbol_mapper_human[gene]
            for gene in self.sigs
            if gene in Plotter.geneid_symbol_mapper_human
        ]
        print(
            f"Mapped {len(sigs_ids)/len(self.sigs)*100:.2f}% of "
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
        go_terms = self.generate_go_table()

        go_terms = go_terms[:20]
        go_terms = go_terms.sort_values("ratio_in_study", ascending=False)

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
        cbar = g.colorbar(
            color_mapper, ax=ax, orientation="vertical", pad=0.01, format="{x:.4f}"
        )
        cbar.ax.set_position([0.8, 0.5, 0.2, 0.3])
        cbar.ax.set_title("padj", loc="left", pad=4.0)

        return g

    def plot_volcano(self):
        """Create volcano plot (placeholder - implementation not shown in source).

        Returns:
            Matplotlib figure object.
        """
        raise NotImplementedError("Volcano plot method needs to be implemented")

    def plot_historgram(self):
        """Create hierarchical clustering heatmap (placeholder - implementation not shown in source).

        Returns:
            Matplotlib figure object.
        """
        raise NotImplementedError("Histogram plot method needs to be implemented")

    def plot_pca(self):
        """Create PCA visualization (placeholder - implementation not shown in source).

        Returns:
            Matplotlib figure object.
        """
        raise NotImplementedError("PCA plot method needs to be implemented")
