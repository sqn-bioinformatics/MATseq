"""DESeq2 analysis utilities."""

import pandas as pd
from anndata import AnnData
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference


class DataProcessor:
    def __init__(
        self,
        raw_counts: pd.DataFrame,
        classes: list[str],
        n_cpus: int,
        batches: list[str] = None,
    ):
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
            negative_control: Name of the negative control class (default: "nc").

        Returns:
            List of [class, negative_control] pairs.
        """
        return [[my_class, negative_control] for my_class in class_list]

    def prepare_metadata(self) -> pd.DataFrame:
        """
        A method that checks if the list of batches was provided.
        If there are no batches, it checks the number of classes and subsets samples that belong to the classes,
        then it will prepare filtered counts and metadata. Otherwise, it will subset by class and batch and check
        if only one class was provided. Multiple batch comparisons are implemented here as only per one class!

        Returns tuple of counts and metadata

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
            )  # Retrieve batch number if the sample name as TA001_batch_class_replica
            metadata = metadata.dropna(subset=["batches"])
            if len(metadata["batches"].unique()) <= 1:
                raise ValueError(
                    f"Provided counts for {len(metadata['batches'].unique())} batch(es). Please provide counts for at least two batches."
                )

        return metadata

    def make_dds(self) -> AnnData:
        metadata = self.prepare_metadata()
        counts = self.raw_counts.loc[metadata.index]

        # Filtering genes that have a total sum of counts lower than 10
        counts = counts[counts.columns[counts.sum(axis=0) >= 10]]

        # Convert the multi-indexed DataFrame to a single-indexed one by resetting the index
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
    ) -> tuple[AnnData, pd.DataFrame, pd.DataFrame]:

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

        # LImit the p value to minimal pvalue detected

        sigs = res[
            (res.padj < padj_value) & (abs(res.log2FoldChange) > log2foldchange_value)
        ]

        try:
            if sigs.empty:
                raise ValueError("No significant gene found.")
        except ValueError as e:
            print(e)

        return dds, res, sigs
