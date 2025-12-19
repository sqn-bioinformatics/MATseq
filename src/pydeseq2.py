"""PyDESeq2 analysis and visualization for differential expression studies."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import textwrap
from pathlib import Path
from typing import Union, Tuple
from matplotlib.patches import Patch
from anndata import AnnData
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

from .utils import save_csv, get_output_path, CUSTOM_PALETTE_6

# Lazy imports for optional dependencies
sc = None
adjust_text = None
GODag = None
Gene2GoReader = None
GOEnrichmentStudyNS = None
GeneID2nt_hs = None


def _load_scanpy():
    """Lazily load scanpy for PCA plotting."""
    global sc
    if sc is None:
        import scanpy as _sc
        sc = _sc


def _load_adjust_text():
    """Lazily load adjustText for volcano plot labels."""
    global adjust_text
    if adjust_text is None:
        from adjustText import adjust_text as _adjust_text
        adjust_text = _adjust_text


def _load_goatools():
    """Lazily load goatools dependencies."""
    global GODag, Gene2GoReader, GOEnrichmentStudyNS, GeneID2nt_hs
    if GODag is None:
        from goatools.obo_parser import GODag as _GODag
        from goatools.anno.genetogo_reader import Gene2GoReader as _Gene2GoReader
        from goatools.goea.go_enrichment_ns import (
            GOEnrichmentStudyNS as _GOEnrichmentStudyNS,
        )
        from .genes_ncbi_homo_sapiens_proteincoding import GENEID2NT as _GeneID2nt_hs

        GODag = _GODag
        Gene2GoReader = _Gene2GoReader
        GOEnrichmentStudyNS = _GOEnrichmentStudyNS
        GeneID2nt_hs = _GeneID2nt_hs


class DataProcessor:
    """Process raw count data and perform DESeq2 analysis."""

    def __init__(
        self,
        raw_counts: pd.DataFrame,
        classes: list[str],
        n_cpus: int,
        batches: list[str] = None,
    ):
        """Initialize data processor.

        Args:
            raw_counts: DataFrame with gene counts (multiindexed by sample/label).
            classes: List of class labels to compare (must be 2 for basic comparison).
            n_cpus: Number of CPUs for parallel processing.
            batches: Optional list of batch identifiers for batch effect analysis.
        """
        self.raw_counts = raw_counts
        self.classes = classes
        self.batches = batches
        self.n_cpus = n_cpus

    @staticmethod
    def make_pairs_with_negative_control(
        class_list: list[str], negative_control: str = "IMDM"
    ) -> list[list[str]]:
        """Create pairs of classes with negative control for comparison.

        Args:
            class_list: List of class names to pair with negative control.
            negative_control: Name of the negative control class.

        Returns:
            List of [class, negative_control] pairs.
        """
        return [[my_class, negative_control] for my_class in class_list]

    def prepare_metadata(self) -> pd.DataFrame:
        """Prepare metadata for DESeq2 analysis.

        Returns:
            DataFrame with metadata for comparison.

        Raises:
            ValueError: If invalid class or batch configuration.
        """
        if not self.batches and len(self.classes) != 2:
            raise ValueError(
                "prepare_metadata() takes a list with two classes (str) as a condition."
            )
        if self.batches and len(self.classes) != 2:
            raise ValueError(
                "prepare_metadata() takes a list with one class (str) when batches provided."
            )

        samples_to_compare = list(
            self.raw_counts[
                self.raw_counts.index.get_level_values("label").isin(self.classes)
            ].index
        )

        number_of_unique_classes = len(
            (
                self.raw_counts.loc[samples_to_compare]
                .index.get_level_values("label")
                .unique()
            )
        )
        if number_of_unique_classes <= 1:
            raise ValueError(
                f"Provided counts for {number_of_unique_classes} class(es). Please provide counts for at least two classes."
            )

        metadata = pd.DataFrame(samples_to_compare, columns=["index", "condition"])
        metadata.set_index("index", inplace=True)

        if self.batches:
            metadata["batches"] = metadata.index.map(
                lambda sample: (
                    sample.split("_")[1]
                    if sample.split("_")[1] in self.batches
                    else None
                )
            )
            metadata = metadata.dropna(subset=["batches"])
            if len(metadata["batches"].unique()) <= 1:
                raise ValueError(
                    f"Provided counts for {len(metadata['batches'].unique())} batch(es). Please provide counts for at least two batches."
                )

        return metadata

    def make_dds(self) -> AnnData:
        """Create and run DESeq2 analysis.

        Returns:
            AnnData object with DESeq2 results.
        """
        metadata = self.prepare_metadata()
        counts = self.raw_counts.loc[metadata.index]

        # Filter genes with low total counts
        counts = counts[counts.columns[counts.sum(axis=0) >= 10]]

        # Convert multi-indexed DataFrame to single-indexed
        counts.reset_index(level="label", drop=True, inplace=True)

        dds = DeseqDataSet(
            counts=counts,
            metadata=metadata,
            design_factors=metadata.columns,
            refit_cooks=True,
            inference=DefaultInference(self.n_cpus),
        )
        dds.deseq2()
        return dds

    def make_statistics(
        self,
        padj_value: float = 0.05,
        log2foldchange_value: int = 2,
    ) -> Tuple[AnnData, pd.DataFrame, pd.DataFrame]:
        """Calculate differential expression statistics.

        Args:
            padj_value: Adjusted p-value threshold for significance.
            log2foldchange_value: Log2 fold-change threshold.

        Returns:
            Tuple of (dds, results_df, significant_genes_df).
        """
        dds = self.make_dds()

        if not self.batches:
            conditions_to_contrast = self.classes
            design_factor = "condition"
        else:
            conditions_to_contrast = self.batches
            design_factor = "group"

        stat_res = DeseqStats(
            dds,
            contrast=(design_factor, *conditions_to_contrast),
            inference=DefaultInference(self.n_cpus),
            quiet=True,
        )

        stat_res.summary()

        res = stat_res.results_df
        res = res[res.baseMean >= 10]

        sigs = res[
            (res.padj < padj_value) & (abs(res.log2FoldChange) > log2foldchange_value)
        ]

        if sigs.empty:
            print("Warning: No significant genes found with current thresholds.")

        return dds, res, sigs


class Plotter:
    """Visualization and GO enrichment analysis for DESeq2 results.

    Maintains class-level GO term data for efficiency across multiple analyses.
    """

    geneid_symbol_mapper_human = None
    goeaobj = None
    is_initialized = False

    def __init__(
        self,
        dds: AnnData,
        res: pd.DataFrame,
        sigs: pd.DataFrame,
        analysis_name: str,
    ):
        """Initialize plotter with DESeq2 results.

        Args:
            dds: AnnData object from DESeq2 analysis.
            res: DataFrame with all DESeq2 results.
            sigs: DataFrame with significant genes.
            analysis_name: Name for this analysis (used in filenames).
        """
        self.dds = dds
        self.res = res
        self.sigs = sigs
        self.analysis_name = analysis_name

        if not Plotter.is_initialized:
            self._initialise_go()
            Plotter.is_initialized = True

    @classmethod
    def _initialise_go(cls):
        """Initialize GO terms database and gene mappings.

        Requires go-basic.obo and gene2go files in working directory or data path.
        """
        print("Initializing GO terms...")
        _load_goatools()

        # Determine data path
        path_to_supporting_files = Path.cwd() / "deseq2"
        if not path_to_supporting_files.exists():
            path_to_supporting_files = Path(__file__).parent.parent / "data" / "deseq2"

        # Check if required files exist
        gene2go_file = path_to_supporting_files / "gene2go"
        obo_file = path_to_supporting_files / "go-basic.obo"

        if not gene2go_file.exists() or not obo_file.exists():
            print(f"Warning: GO data files not found at {path_to_supporting_files}")
            print("GO enrichment analysis will be skipped.")
            return

        try:
            genes = Gene2GoReader(
                gene2go_file, taxids=[9606], namespaces={"BP"}
            )
            ns2assoc = genes.get_ns2assc()
            obodag = GODag(obo_file)

            cls.goeaobj = GOEnrichmentStudyNS(
                GeneID2nt_hs.keys(),
                ns2assoc,
                obodag,
                propagate_counts=False,
                alpha=0.05,
                methods=["fdr_bh"],
            )

            cls.geneid_symbol_mapper_human = {
                GeneID2nt_hs[key].Symbol: GeneID2nt_hs[key].GeneID
                for key in GeneID2nt_hs
            }
            print("GO terms initialized successfully.")
        except Exception as e:
            print(f"Warning: Could not initialize GO terms - {e}")
            print("GO enrichment analysis will be skipped.")

    def make_figure(self, plot_type: str):
        """Generate and save a figure.

        Args:
            plot_type: Type of plot ('volcano', 'histogram', 'pca', or 'go').

        Raises:
            ValueError: If plot_type is invalid or no significant genes found.
        """
        if plot_type not in ["volcano", "histogram", "pca", "go"]:
            raise ValueError(
                "plot_type must be 'volcano', 'histogram', 'pca', or 'go'."
            )

        if plot_type in ["volcano", "histogram", "go"] and self.sigs.empty:
            print(f"Skipping {plot_type} plot: No significant genes found.")
            return

        fontsize = 12
        fig_extension = "png"
        output_path = get_output_path()
        fname = output_path / f"{self.analysis_name}_{plot_type}.{fig_extension}"

        if plot_type == "volcano":
            fig = self.plot_volcano()
            plt.title(
                f"{self.analysis_name} Differentially Expressed Genes",
                fontsize=fontsize,
                fontweight="bold",
            )

        elif plot_type == "histogram":
            fig = self.plot_histogram()
            fig.ax_col_dendrogram.set_title(
                f"{self.analysis_name} Differentially Expressed Genes",
                fontsize=fontsize,
                fontweight="bold",
                pad=2,
            )

        elif plot_type == "pca":
            fig = self.plot_pca()
            fig.axes[0].set_title(
                f"{self.analysis_name}", fontsize=fontsize, fontweight="bold"
            )

        elif plot_type == "go":
            fig = self.plot_go()
            plt.title(
                f"{self.analysis_name} Top 20 Significant GO Terms",
                fontsize=fontsize,
                fontweight="bold",
            )
            go_df = self.generate_go_table()
            save_csv(go_df, f"{self.analysis_name}_go_terms")

        fig.figure.savefig(fname, format=fig_extension, dpi=300, bbox_inches="tight")
        print(f"Figure saved: {fname}")

    def plot_volcano(self, log2foldchange: float = 2):
        """Create volcano plot of differentially expressed genes.

        Args:
            log2foldchange: Log2 fold-change threshold for highlighting.

        Returns:
            Figure object.
        """
        _load_adjust_text()

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

        grapher_subset = grapher[
            grapher["color"].isin(["overexpressed", "underexpressed"])
        ]

        print(f"Number of DE genes: {len(grapher_subset)}")

        sorted_padj = grapher_subset.sort_values("padj_log", ascending=False)
        sorted_lfc = grapher_subset.sort_values("log2FoldChange", ascending=True)

        annotation_subset = pd.concat(
            [
                sorted_padj.head(20),
                sorted_lfc.head(10),
                sorted_lfc.tail(10),
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

        ax.axhline(1.3, zorder=1)
        ax.axvline(2, zorder=1)
        ax.axvline(-2, zorder=1)

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

        return g

    def plot_histogram(self, num_top_sig: Union[int, str] = 50):
        """Create heatmap of top differentially expressed genes.

        Args:
            num_top_sig: Number of top significant genes to plot ('all' for all).

        Returns:
            Clustermap figure object.
        """
        _load_scanpy()

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

        return ax

    def plot_pca(self, w_text: bool = False):
        """Create PCA plot of samples.

        Args:
            w_text: If True, add sample labels to plot.

        Returns:
            Scanpy figure object.
        """
        _load_scanpy()

        dds = self.dds.copy()
        sc.tl.pca(dds, n_comps=2)
        pca = dds.obsm["X_pca"]
        sample_names = list(sc.get.obs_df(dds).index)

        rc = {
            "axes.facecolor": "white",
            "axes.edgecolor": "black",
        }
        sns.set_theme(rc=rc)
        plt.figure(figsize=(8, 10))
        ax = sc.pl.pca(
            dds, color="condition", size=300, show=False, title=" ", return_fig=True
        )

        if w_text:
            texts = [
                plt.text(pca[i][0], pca[i][1], sample_names[i], ha="left", va="bottom")
                for i in range(len(pca))
            ]

        return ax

    def generate_go_table(self) -> pd.DataFrame:
        """Generate table of significant GO terms.

        Returns:
            DataFrame with GO term information and gene associations.

        Raises:
            RuntimeError: If GO enrichment was not initialized.
        """
        if not Plotter.goeaobj:
            raise RuntimeError("GO enrichment not initialized. Check data files.")

        sigs_ids = [
            Plotter.geneid_symbol_mapper_human[gene]
            for gene in self.sigs.index
            if gene in Plotter.geneid_symbol_mapper_human
        ]

        print(
            f"Mapped {len(sigs_ids)/len(self.sigs.index)*100:.2f}% of gene symbols."
        )

        goea_results = Plotter.goeaobj.run_study(sigs_ids, prt=None)
        goea_results_sig = [r for r in goea_results if r.p_fdr_bh < 0.05]

        inverted_mapping = {
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
                        list(map(lambda y: inverted_mapping[y], x.study_items)),
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

    def plot_go(self) -> plt.Figure:
        """Create barplot of top GO terms by enrichment.

        Returns:
            Figure object.
        """
        go_terms = self.generate_go_table()
        go_terms = go_terms[:20]
        go_terms = go_terms.sort_values("ratio_in_study", ascending=False)

        norm = mpl.colors.Normalize(vmin=go_terms.fdr.min(), vmax=go_terms.fdr.max())
        color_mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.bwr_r)

        rc = {"axes.edgecolor": "black"}
        sns.set_theme(rc=rc)

        g = plt.figure(figsize=(10, 12))

        ax = sns.barplot(
            data=go_terms,
            x=go_terms["n_genes"] / go_terms["n_study"],
            y="term",
            palette=list(color_mapper.to_rgba(go_terms.fdr.values)),
        )

        ax.set_yticklabels(
            [textwrap.fill(term, 40) for term in go_terms["term"]],
            fontsize=10,
        )

        cbar = g.colorbar(
            color_mapper, ax=ax, orientation="vertical", pad=0.01, format="{x:.2f}"
        )
        cbar.ax.set_position([0.8, 0.5, 0.2, 0.3])
        cbar.ax.set_title("padj", loc="left", pad=4.0)

        return g


class AnalysisPipeline:
    """Execute DESeq2 analysis workflow with caching and figure generation."""

    def __init__(
        self,
        raw_counts: pd.DataFrame,
        sample_labels: pd.Series = None,
        output_dir: Path = None,
        padj_threshold: float = 0.05,
        log2fc_threshold: int = 2,
        n_cpus: int = 42,
    ):
        """Initialize analysis pipeline.

        Args:
            raw_counts: DataFrame with gene counts (samples as rows, genes as columns).
            sample_labels: Series with sample labels. If None, extracts from index by splitting on "_".
            output_dir: Directory for output files. Defaults to results/differential_gene_expression.
            padj_threshold: Adjusted p-value threshold for significance.
            log2fc_threshold: Log2 fold-change threshold.
            n_cpus: Number of CPUs for parallel processing.
        """
        # Prepare raw counts with multiindex if needed
        if sample_labels is not None:
            self.raw_counts = raw_counts.copy()
            self.raw_counts.index = pd.MultiIndex.from_arrays(
                [self.raw_counts.index, sample_labels], names=["samples", "label"]
            )
        elif isinstance(raw_counts.index, pd.MultiIndex):
            self.raw_counts = raw_counts
        else:
            # Extract labels from sample names
            labels = [i.split("_")[2] for i in raw_counts.index]
            self.raw_counts = raw_counts.copy()
            self.raw_counts.index = pd.MultiIndex.from_arrays(
                [self.raw_counts.index, labels], names=["samples", "label"]
            )

        self.padj_threshold = padj_threshold
        self.log2fc_threshold = log2fc_threshold
        self.n_cpus = n_cpus

        if output_dir is None:
            output_dir = Path.cwd() / "results" / "differential_gene_expression"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.figures_dir = Path.cwd() / "results" / "figures" / "deseq2"
        self.figures_dir.mkdir(parents=True, exist_ok=True)

        self.results = {}
        self.de_genes = set()

    def run_analysis(
        self, class_list: list[str], negative_control: str = "nc"
    ) -> dict:
        """Run DESeq2 analysis for all class pairs with negative control.

        Args:
            class_list: List of ligand classes to analyze.
            negative_control: Name of negative control class.

        Returns:
            Dictionary with results for each ligand.
        """
        # Create pairs with negative control
        pairs = DataProcessor.make_pairs_with_negative_control(
            class_list, negative_control
        )

        results_list = []

        for class_pair in pairs:
            ligand_name = class_pair[0]
            print(f"Running analysis for {ligand_name}...")

            processor = DataProcessor(
                raw_counts=self.raw_counts,
                classes=class_pair,
                n_cpus=self.n_cpus,
            )

            dds, res, sigs = processor.make_statistics(
                padj_value=self.padj_threshold,
                log2foldchange_value=self.log2fc_threshold,
            )

            # Store results
            self.results[ligand_name] = {
                "dds": dds,
                "results": res,
                "significant": sigs,
            }

            # Save results CSV
            res_output = self.output_dir / f"{ligand_name}_results.csv"
            res.to_csv(res_output)
            print(f"Saved results to {res_output}")

            # Track DE genes
            self.de_genes.update(sigs.index)

            # Add results for merging
            res_copy = res.copy()
            res_copy.columns = [f"{col}_{ligand_name}" for col in res_copy.columns]
            results_list.append(res_copy)

            # Generate figures
            if not sigs.empty:
                self._generate_figures(ligand_name, dds, res, sigs)

        # Merge all results
        if results_list:
            merged = results_list[0]
            for res_df in results_list[1:]:
                merged = pd.merge(
                    merged,
                    res_df,
                    left_index=True,
                    right_index=True,
                    how="outer",
                )

            merged_output = self.output_dir / "merged_results.csv"
            merged.to_csv(merged_output)
            print(f"Saved merged results to {merged_output}")

        return self.results

    def _generate_figures(
        self, ligand_name: str, dds: AnnData, res: pd.DataFrame, sigs: pd.DataFrame
    ):
        """Generate visualization figures for analysis.

        Args:
            ligand_name: Name of ligand for filenames.
            dds: DESeq2 AnnData object.
            res: Full results DataFrame.
            sigs: Significant genes DataFrame.
        """
        plotter = Plotter(dds, res, sigs, ligand_name)

        # Generate individual plots
        plotter.make_figure("volcano")
        plotter.make_figure("histogram")
        plotter.make_figure("pca")

        # Try GO enrichment if available
        try:
            plotter.make_figure("go")
        except Exception as e:
            print(f"Warning: GO enrichment failed for {ligand_name}: {e}")

    def save_de_genes(self, filename: str = "de_genes.txt"):
        """Save list of all differentially expressed genes.

        Args:
            filename: Name of output file.
        """
        output_path = self.output_dir / filename
        with open(output_path, "w") as f:
            for gene in sorted(self.de_genes):
                f.write(f"{gene}\n")
        print(f"Saved {len(self.de_genes)} DE genes to {output_path}")

    def create_supplementary_figure_1a(
        self, comparison_pairs: list[Tuple[str, str]] = None
    ):
        """Create 1x2 supplementary figure with volcano and histogram.

        Args:
            comparison_pairs: List of (ligand, negative_control) pairs to plot.
                            If None, uses all computed results.
        """
        if comparison_pairs is None:
            comparison_pairs = [(k, "nc") for k in self.results.keys()]

        for ligand_name, nc_name in comparison_pairs:
            if ligand_name not in self.results:
                print(f"Warning: No results for {ligand_name}")
                continue

            dds = self.results[ligand_name]["dds"]
            res = self.results[ligand_name]["results"]
            sigs = self.results[ligand_name]["significant"]

            if sigs.empty:
                print(f"Skipping {ligand_name}: No significant genes")
                continue

            plotter = Plotter(dds, res, sigs, ligand_name)

            # Create 1x2 layout: volcano + histogram
            fig, axes = plt.subplots(1, 2, figsize=(16, 6))

            # Volcano plot (left)
            _plot_volcano_on_ax(axes[0], plotter, res)

            # Histogram/heatmap (right)
            # Use simplified version since we need subplot integration
            num_top_sig = 50
            if num_top_sig != "all":
                sigs_plot = sigs.sort_values("padj")[:num_top_sig]
            else:
                sigs_plot = sigs

            dds_sigs = dds[:, sigs_plot.index]
            dds_sigs.layers["log1p"] = np.log1p(dds_sigs.layers["normed_counts"])

            grapher = pd.DataFrame(
                dds_sigs.layers["log1p"].T,
                index=dds_sigs.var_names,
                columns=dds_sigs.obs.condition,
            )

            sns.heatmap(grapher, cmap="RdYlBu_r", ax=axes[1], cbar_kws=dict(label="log1p counts"))
            axes[1].set_title(f"{ligand_name} Top {len(sigs_plot)} DE Genes")
            axes[1].set_xlabel("Sample Condition")
            axes[1].set_ylabel("Gene")

            plt.suptitle(f"{ligand_name} vs {nc_name}", fontsize=14, fontweight="bold")
            plt.tight_layout()

            output_path = self.figures_dir / f"Supplementary_Figure_1A_{ligand_name}.png"
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
            print(f"Saved supplementary figure to {output_path}")
            plt.close(fig)


def _plot_volcano_on_ax(ax, plotter: Plotter, res: pd.DataFrame, log2foldchange: float = 2):
    """Plot volcano plot on existing matplotlib axis.

    Args:
        ax: Matplotlib axis to plot on.
        plotter: Plotter instance with data.
        res: Results DataFrame.
        log2foldchange: Log2 fold-change threshold.
    """
    _load_adjust_text()

    grapher = res.assign(
        padj_log=res["padj"].apply(
            lambda x: -np.log10(x) if x != 0 else -np.log10(x + 1e-300)
        ),
        color="no_expression_change",
    )

    grapher.loc[grapher["log2FoldChange"] > log2foldchange, "color"] = "overexpressed"
    grapher.loc[grapher["log2FoldChange"] < -log2foldchange, "color"] = "underexpressed"

    grapher_subset = grapher[grapher["color"].isin(["overexpressed", "underexpressed"])]

    sorted_padj = grapher_subset.sort_values("padj_log", ascending=False)
    sorted_lfc = grapher_subset.sort_values("log2FoldChange", ascending=True)

    annotation_subset = pd.concat(
        [sorted_padj.head(20), sorted_lfc.head(10), sorted_lfc.tail(10)]
    ).drop_duplicates()

    ax.scatter(
        grapher[grapher["color"] == "no_expression_change"]["log2FoldChange"],
        grapher[grapher["color"] == "no_expression_change"]["padj_log"],
        c="grey",
        alpha=0.5,
        s=20,
        label="No change",
    )
    ax.scatter(
        grapher[grapher["color"] == "overexpressed"]["log2FoldChange"],
        grapher[grapher["color"] == "overexpressed"]["padj_log"],
        c="orange",
        alpha=0.7,
        s=30,
        label="Overexpressed",
    )
    ax.scatter(
        grapher[grapher["color"] == "underexpressed"]["log2FoldChange"],
        grapher[grapher["color"] == "underexpressed"]["padj_log"],
        c="purple",
        alpha=0.7,
        s=30,
        label="Underexpressed",
    )

    ax.axhline(1.3, color="k", linestyle="--", linewidth=0.5)
    ax.axvline(2, color="k", linestyle="--", linewidth=0.5)
    ax.axvline(-2, color="k", linestyle="--", linewidth=0.5)

    texts = []
    for i, row in annotation_subset.iterrows():
        texts.append(
            ax.text(
                x=row.log2FoldChange,
                y=row.padj_log,
                s=row.name,
                fontsize=7,
                weight="bold",
            )
        )

    if texts:
        adjust_text(texts, arrowprops=dict(arrowstyle="-", color="k"), ax=ax)

    ax.set_xlabel("log2 Fold Change", fontsize=11)
    ax.set_ylabel("-log10 FDR", fontsize=11)
    ax.set_ylim(-2, grapher["padj_log"].max() + 5)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
