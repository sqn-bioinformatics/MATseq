"""PyDESeq2 analysis and visualization for differential expression studies."""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple, List, Optional
from anndata import AnnData
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

from .visualization import plot_volcano, plot_heatmap, plot_go
from .go_term_analysis import run_go_analysis
from .cache import PipelineCache


class DataProcessor:
    """Process raw count data and perform DESeq2 analysis."""

    def __init__(
        self,
        raw_counts: pd.DataFrame,
        sample_labels: pd.Series,
        classes: list[str],
        n_cpus: int,
        cache: Optional[PipelineCache] = None,
        force_recompute: bool = False,
    ):
        """Initialize data processor.

        Args:
            raw_counts: DataFrame with gene counts (samples as rows, genes as columns).
            sample_labels: Series with sample class labels.
            classes: List of class labels to compare (must be 2 for comparison).
            n_cpus: Number of CPUs for parallel processing.
            cache: PipelineCache instance for disk-based caching.
            force_recompute: Skip cache and recompute.
        """
        assert isinstance(
            raw_counts, pd.DataFrame
        ), f"raw_counts must be DataFrame, got {type(raw_counts)}"
        assert isinstance(
            sample_labels, pd.Series
        ), f"sample_labels must be Series, got {type(sample_labels)}"
        assert isinstance(classes, list), f"classes must be list, got {type(classes)}"
        assert isinstance(n_cpus, int), f"n_cpus must be int, got {type(n_cpus)}"

        if len(raw_counts) != len(sample_labels):
            raise ValueError(
                f"raw_counts ({len(raw_counts)}) and sample_labels ({len(sample_labels)}) must have same length"
            )
        self.raw_counts = raw_counts
        self.sample_labels = sample_labels
        self.classes = classes
        self.n_cpus = n_cpus
        self.cache = cache
        self.force_recompute = force_recompute

    def prepare_metadata(self) -> pd.DataFrame:
        """Prepare metadata for DESeq2 analysis.

        Returns:
            DataFrame with metadata (condition column) indexed by sample names.

        Raises:
            ValueError: If invalid class configuration.
        """
        if len(self.classes) != 2:
            raise ValueError(
                f"Expected exactly 2 classes for comparison, got {len(self.classes)}: {self.classes}"
            )

        # Filter samples and labels for the classes being compared
        mask = self.sample_labels.isin(self.classes)
        sample_names = self.raw_counts.index[mask]
        labels = self.sample_labels[mask]

        unique_labels = labels.unique()
        if len(unique_labels) < 2:
            raise ValueError(
                f"Must have at least 2 unique classes in data. Got {len(unique_labels)}: {unique_labels.tolist()}"
            )

        # Verify both classes exist in the filtered data
        for cls in self.classes:
            if cls not in unique_labels:
                raise ValueError(
                    f"Class '{cls}' not found in sample labels. Available: {unique_labels.tolist()}"
                )

        # Create metadata DataFrame with condition as categorical
        # This is required by pydeseq2 to properly set reference levels
        condition_cat = pd.Categorical(
            labels.values, categories=self.classes, ordered=False
        )

        metadata = pd.DataFrame({"condition": condition_cat}, index=sample_names)

        return metadata

    def make_dds(self) -> AnnData:
        """Create and run DESeq2 analysis.

        Returns:
            AnnData object with DESeq2 results.
        """
        metadata = self.prepare_metadata()
        counts = self.raw_counts.loc[metadata.index].copy()

        # Filter genes with low total counts
        counts = counts[counts.columns[counts.sum(axis=0) >= 10]]

        print(f"DESeq2 analysis: {counts.shape[0]} samples, {counts.shape[1]} genes")

        dds = DeseqDataSet(
            counts=counts,
            metadata=metadata,
            design_factors=["condition"],
            refit_cooks=True,
            inference=DefaultInference(self.n_cpus),
            quiet=True,
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
        class_key = "_".join(self.classes)
        cache_name = f"deseq2_{class_key}"
        cache_params = {
            "padj": padj_value,
            "log2fc": log2foldchange_value,
        }

        if self.cache is not None and not self.force_recompute:
            cached = self.cache.get(cache_name, cache_params)
            if cached is not None:
                print(f"Loading DESeq2 results for {class_key} from cache...")
                return cached

        dds = self.make_dds()

        tested_level = self.classes[0].replace("_", "-")
        ref_level = self.classes[1].replace("_", "-")
        contrast = ["condition", tested_level, ref_level]

        print(f"DESeq2 contrast: {tested_level} vs {ref_level}")

        stat_res = DeseqStats(
            dds,
            contrast=contrast,
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
            print(
                f"Warning: No significant genes found (padj < {padj_value}, |log2FC| > {log2foldchange_value})"
            )

        result = (dds, res, sigs)

        if self.cache is not None:
            self.cache.set(cache_name, result, cache_params, f"DESeq2 {class_key}")

        return result


class AnalysisPipeline:
    """Execute DESeq2 analysis workflow with caching and figure generation."""

    def __init__(
        self,
        raw_counts: pd.DataFrame,
        sample_labels: pd.Series,
        output_dir: Path = None,
        padj_threshold: float = 0.05,
        log2fc_threshold: int = 2,
        n_cpus: int = 42,
        cache: Optional[PipelineCache] = None,
        force_recompute: bool = False,
    ):
        """Initialize analysis pipeline.

        Args:
            raw_counts: DataFrame with gene counts (samples as rows, genes as columns).
            sample_labels: Series with sample class labels.
            output_dir: Directory for output files. Defaults to results/differential_gene_expression.
            padj_threshold: Adjusted p-value threshold for significance.
            log2fc_threshold: Log2 fold-change threshold.
            n_cpus: Number of CPUs for parallel processing.
            cache: PipelineCache instance for disk-based caching.
            force_recompute: Skip cache and recompute.
        """

        if len(raw_counts) != len(sample_labels):
            raise ValueError(
                f"raw_counts ({len(raw_counts)}) and sample_labels ({len(sample_labels)}) must have same length"
            )

        self.raw_counts = raw_counts.copy()
        self.sample_labels = pd.Series(sample_labels.values, index=raw_counts.index)

        self.padj_threshold = padj_threshold
        self.log2fc_threshold = log2fc_threshold
        self.n_cpus = n_cpus
        self.cache = cache
        self.force_recompute = force_recompute

        if output_dir is None:
            output_dir = Path.cwd() / "results" / "differential_gene_expression"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.figures_dir = Path.cwd() / "results" / "figures" / "deseq2"
        self.figures_dir.mkdir(parents=True, exist_ok=True)

        self.results = {}
        self.de_genes = set()

    def run_analysis(
        self, class_list: list[str], negative_control: str = "negative_control"
    ) -> dict:
        """Run DESeq2 analysis for all class pairs with negative control.

        Args:
            class_list: List of ligand classes to analyze.
            negative_control: Name of negative control class.

        Returns:
            Dictionary with results for each ligand.
        """
        filtered_list = [c for c in class_list if c != negative_control]
        pairs = [[my_class, negative_control] for my_class in filtered_list]
        results_list = []

        for class_pair in pairs:
            ligand_name = class_pair[0]
            print(f"Running analysis for {ligand_name}...")

            processor = DataProcessor(
                raw_counts=self.raw_counts,
                sample_labels=self.sample_labels,
                classes=class_pair,
                n_cpus=self.n_cpus,
                cache=self.cache,
                force_recompute=self.force_recompute,
            )

            dds, res, sigs = processor.make_statistics(
                padj_value=self.padj_threshold,
                log2foldchange_value=self.log2fc_threshold,
            )

            self.results[ligand_name] = {
                "dds": dds,
                "results": res,
                "significant": sigs,
            }

            res_output = self.output_dir / f"{ligand_name}_deseq2_results.csv"
            res.to_csv(res_output)
            print(f"Saved results to {res_output}")

            self.de_genes.update(sigs.index)

            res_copy = res.copy()
            res_copy.columns = [f"{col}_{ligand_name}" for col in res_copy.columns]
            results_list.append(res_copy)

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

        return processor, self.results

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
        plot_volcano(res, ligand_name, output_path=self.figures_dir)
        plot_heatmap(dds, sigs, ligand_name, output_path=self.figures_dir)

        try:
            go_terms_dir = Path.cwd() / "results" / "go_terms"
            go_terms_dir.mkdir(parents=True, exist_ok=True)
            go_df = run_go_analysis(sigs, ligand_name, output_dir=go_terms_dir)
            if not go_df.empty:
                go_fig_dir = Path.cwd() / "results" / "figures" / "go"
                plot_go(
                    go_df,
                    title=f"{ligand_name} Top 20 Significant GO Terms",
                    output_path=go_fig_dir,
                    output_filename=f"{ligand_name}_go.png",
                )
        except Exception as e:
            print(f"Warning: GO enrichment failed for {ligand_name}: {e}")

    def get_de_genes(self):
        """Provide list of all differentially expressed genes.

        Args:
            None.
        """
        return self.de_genes


def merge_go_tables(
    go_files: List[Path],
    output_dir: Path = None,
    output_filename: str = "Supplementary_Table_3_genes_GO.csv",
) -> pd.DataFrame:
    """Merge GO enrichment results from multiple ligand analyses.

    Args:
        go_files: List of paths to GO result CSV files.
        output_dir: Directory for output files. Defaults to results/differential_gene_expression.
        output_filename: Name for merged output file.

    Returns:
        Merged DataFrame with all GO terms.
    """
    if not go_files:
        raise ValueError("No GO files provided")

    if output_dir is None:
        output_dir = Path.cwd() / "results" / "differential_gene_expression"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    merged_df = None

    for go_file in go_files:
        ligand_name = go_file.stem.split("_")[0]
        df = pd.read_csv(go_file, index_col=0)

        df.insert(0, "Ligand", ligand_name)

        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.concat([merged_df, df], axis=0)

    merged_df = merged_df.sort_values("fdr", ascending=True)

    output_path = output_dir / output_filename
    merged_df.to_csv(output_path)
    print(f"Saved merged GO table to {output_path}")

    return merged_df
