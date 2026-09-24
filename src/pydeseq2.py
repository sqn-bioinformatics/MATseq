"""PyDESeq2 analysis and visualization for differential expression studies."""

import pandas as pd
from pathlib import Path
from anndata import AnnData
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

from .visualization import plot_volcano, plot_heatmap, plot_go
from .go_term_analysis import merge_go_tables, run_go_analysis

class DataProcessor:
    """Process raw count data and perform DESeq2 analysis."""

    def __init__(
        self,
        raw_counts: pd.DataFrame,
        sample_labels: pd.Series,
        classes: list[str],
        n_cpus: int,
    ):
        self.raw_counts = raw_counts
        self.sample_labels = sample_labels
        self.classes = classes
        self.n_cpus = n_cpus

    def prepare_metadata(self) -> pd.DataFrame:
        """Restrict the samples to the two compared classes and build DESeq2 metadata."""
        mask = self.sample_labels.isin(self.classes)
        # Categorical with explicit categories is what sets the pydeseq2 reference level.
        condition = pd.Categorical(
            self.sample_labels[mask].values, categories=self.classes, ordered=False
        )
        return pd.DataFrame(
            {"condition": condition}, index=self.raw_counts.index[mask]
        )

    def make_dds(self) -> AnnData:
        """Create and run DESeq2 analysis."""
        metadata = self.prepare_metadata()
        counts = self.raw_counts.loc[metadata.index].copy()
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
    ) -> tuple[AnnData, pd.DataFrame, pd.DataFrame]:
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
        output_dir: Path,
        figures_dir: Path,
        go_terms_dir: Path,
        go_fig_dir: Path,
        goeaobj,
        geneid_symbol_mapper: dict,
        padj_threshold: float = 0.05,
        log2fc_threshold: float = 2.0,
        n_cpus: int = 42,
        name: str | None = None,
    ):
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

        self.output_dir = output_dir
        self.figures_dir = figures_dir
        self.go_terms_dir = go_terms_dir
        self.go_fig_dir = go_fig_dir
        self.goeaobj = goeaobj
        self.geneid_symbol_mapper = geneid_symbol_mapper
        for directory in (output_dir, figures_dir, go_terms_dir, go_fig_dir):
            directory.mkdir(parents=True, exist_ok=True)

        self.results = {}
        self.de_genes = set()

    def run_analysis(
        self, class_list: list[str], class_to_compare_to: str = "negative_control"
    ) -> dict:
        """Run DESeq2 analysis against a condition/class."""
        present = set(self.sample_labels.unique())
        filtered_list = [
            c for c in class_list if c != class_to_compare_to and c in present
        ]
        pairs = [[my_class, class_to_compare_to] for my_class in filtered_list]

        for class_pair in pairs:
            ligand_name = class_pair[0]
            print(f"Running analysis for {ligand_name}...")

            processor = DataProcessor(
                raw_counts=self.raw_counts,
                sample_labels=self.sample_labels,
                classes=class_pair,
                n_cpus=self.n_cpus,
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
        """Generate visualization figures for analysis."""
        plot_volcano(res, ligand_name, output_path=self.figures_dir)
        plot_heatmap(dds, sigs, ligand_name, output_path=self.figures_dir)

        try:
            go_df = run_go_analysis(
                set(sigs.index), ligand_name, self.go_terms_dir,
                self.goeaobj, self.geneid_symbol_mapper,
            )
            if not go_df.empty:
                plot_go(
                    go_df,
                    condition=ligand_name,
                    output_path=self.go_fig_dir,
                    output_filename=f"{ligand_name}_go.png",
                )
        except Exception as e:
            print(f"Warning: GO enrichment failed for {ligand_name}: {e}")

    def get_de_genes(self) -> set[str]:
        """Return set of all differentially expressed genes."""
        return self.de_genes
