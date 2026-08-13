"""PyDESeq2 analysis and visualization for differential expression studies."""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple
from anndata import AnnData
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

from .visualization import plot_volcano, plot_heatmap, plot_go
from .go_term_analysis import (
    initialize_go,
    generate_go_table,
    merge_go_tables,
)
class DataProcessor:
    """Process raw count data and perform DESeq2 analysis."""

    def __init__(
        self,
        raw_counts: pd.DataFrame,
        sample_labels: pd.Series,
        classes: list[str],
        n_cpus: int,
        name: str = None,
    ):
        """Initialize data processor.
        """
        if not isinstance(raw_counts, pd.DataFrame):
            raise TypeError(f"raw_counts must be DataFrame, got {type(raw_counts)}")
        if not isinstance(sample_labels, pd.Series):
            raise TypeError(f"sample_labels must be Series, got {type(sample_labels)}")
        if not isinstance(classes, list):
            raise TypeError(f"classes must be list, got {type(classes)}")
        if not isinstance(n_cpus, int):
            raise TypeError(f"n_cpus must be int, got {type(n_cpus)}")

        if len(raw_counts) != len(sample_labels):
            raise ValueError(
                f"raw_counts ({len(raw_counts)}) and sample_labels ({len(sample_labels)}) must have same length"
            )
        self.raw_counts = raw_counts
        self.sample_labels = sample_labels
        self.classes = classes
        self.n_cpus = n_cpus
        self.name = name

    def prepare_metadata(self) -> pd.DataFrame:
        """Prepare metadata for DESeq2 analysis.
        """
        if len(self.classes) != 2:
            raise ValueError(
                f"Expected exactly 2 classes for comparison, got {len(self.classes)}: {self.classes}"
            )

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
        """
        metadata = self.prepare_metadata()
        counts = self.raw_counts.loc[metadata.index].copy()

        # Keep genes with >= 10 counts in at least the smallest group's sample count
        min_group = int(metadata["condition"].value_counts().min())
        counts = counts.loc[:, (counts >= 10).sum(axis=0) >= min_group]

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
        padj_threshold: float = 0.05,
        log2fc_threshold: float = 2.0,
    ) -> Tuple[AnnData, pd.DataFrame, pd.DataFrame]:
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
        sigs = res[
            (res.padj < padj_threshold) & (abs(res.log2FoldChange) > log2fc_threshold)
        ]

        if sigs.empty:
            print(
                f"Warning: No significant genes found (padj < {padj_threshold}, |log2FC| > {log2fc_threshold})"
            )

        return dds, res, sigs


class DESeq2:
    def __init__(
        self,
        raw_counts: pd.DataFrame,
        sample_labels: pd.Series,
        output_dir: Path = None,
        padj_threshold: float = 0.05,
        log2fc_threshold: float = 2.0,
        n_cpus: int = 42,
        name: str = None,
    ):
        """Initialize analysis pipeline.
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
        self.name = name

        results = Path.cwd() / "results"
        if output_dir is None:
            output_dir = results / "differential_gene_expression"
            if name:
                output_dir = output_dir / name
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.figures_dir = results / "figures" / "deseq2"
        self.go_terms_dir = results / "go_terms"
        self.go_fig_dir = results / "figures" / "go"
        if name:
            self.figures_dir = self.figures_dir / name
            self.go_terms_dir = self.go_terms_dir / name
            self.go_fig_dir = self.go_fig_dir / name
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        self.go_terms_dir.mkdir(parents=True, exist_ok=True)
        self.go_fig_dir.mkdir(parents=True, exist_ok=True)

        self.results = {}
        self.de_genes = set()
        self._goeaobj = None
        self._geneid_symbol_mapper = None

    def run_analysis(
        self, class_list: list[str], negative_control: str = "negative_control"
    ) -> dict:
        """Run DESeq2 analysis for all class pairs with negative control.
        """
        present = set(self.sample_labels.unique())
        filtered_list = [
            c for c in class_list if c != negative_control and c in present
        ]
        pairs = [[my_class, negative_control] for my_class in filtered_list]

        for class_pair in pairs:
            ligand_name = class_pair[0]
            print(f"Running analysis for {ligand_name}...")

            processor = DataProcessor(
                raw_counts=self.raw_counts,
                sample_labels=self.sample_labels,
                classes=class_pair,
                n_cpus=self.n_cpus,
                name=self.name,
            )

            dds, res, sigs = processor.make_statistics(
                padj_threshold=self.padj_threshold,
                log2fc_threshold=self.log2fc_threshold,
            )

            self.results[ligand_name] = {
                "dds": dds,
                "results": res,
                "significant": sigs,
            }

            res_output = self.output_dir / f"{ligand_name}_deseq2_results.csv"
            res.to_csv(res_output, index_label="gene")
            print(f"Saved results to {res_output}")

            self.de_genes.update(sigs.index)

            res_copy = res.copy()
            res_copy.columns = [f"{col}_{ligand_name}" for col in res_copy.columns]

            if not sigs.empty:
                self._generate_figures(ligand_name, dds, res, sigs)

        go_files = [
            self.go_terms_dir / f"{ligand}_go_terms.csv"
            for ligand in filtered_list
            if (self.go_terms_dir / f"{ligand}_go_terms.csv").exists()
        ]
        if go_files:
            merge_go_tables(go_files, output_dir=self.go_terms_dir)

        return self.results

    def _generate_figures(
        self, ligand_name: str, dds: AnnData, res: pd.DataFrame, sigs: pd.DataFrame
    ):
        """Generate visualization figures for analysis.
        """
        plot_volcano(res, ligand_name, output_path=self.figures_dir)
        plot_heatmap(dds, sigs, ligand_name, output_path=self.figures_dir)

        try:
            go_df = self._go_table(ligand_name, set(sigs.index))
            go_df.to_csv(self.go_terms_dir / f"{ligand_name}_go_terms.csv")
            if not go_df.empty:
                plot_go(
                    go_df,
                    title=f"{ligand_name} Top 20 Significant GO Terms",
                    output_path=self.go_fig_dir,
                    output_filename=f"{ligand_name}_go.png",
                )
        except Exception as e:
            print(f"Warning: GO enrichment failed for {ligand_name}: {e}")

    def get_go_objects(self) -> tuple:
        """Return (goeaobj, geneid_symbol_mapper), initializing if needed."""
        if self._goeaobj is None:
            self._goeaobj, self._geneid_symbol_mapper = initialize_go()
        return self._goeaobj, self._geneid_symbol_mapper

    def _go_table(self, ligand_name: str, genes: set) -> pd.DataFrame:
        goeaobj, mapper = self.get_go_objects()
        return generate_go_table(set(genes), goeaobj, mapper)

    def get_de_genes(self):
        """Return set of all differentially expressed genes."""
        return self.de_genes
