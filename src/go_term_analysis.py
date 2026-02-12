"""GO term enrichment analysis using goatools."""

import gzip
import shutil
import pandas as pd
from pathlib import Path
from collections import namedtuple
from typing import Dict, Any, List


from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS


def _get_data_dir() -> Path:
    """Get data directory for GO files.

    Returns:
        Path to data/go_terms_support directory.
    """
    project_root = Path(__file__).parent.parent
    data_dir = project_root / "data" / "go_terms_support"
    if not data_dir.exists():
        data_dir = Path.cwd() / "data" / "go_terms_support"
    data_dir.mkdir(parents=True, exist_ok=True)
    return data_dir


def _parse_geneid2nt(gene_file: Path) -> Dict[int, Any]:
    """Parse NCBI gene result file to namedtuple dictionary.

    Args:
        gene_file: Path to gene_result_ncbi_human_proteincoding.txt.

    Returns:
        Dictionary mapping GeneID to namedtuple with gene info.
    """
    ntncbi = namedtuple(
        "ntncbi",
        "tax_id Org_name GeneID CurrentID Status Symbol Aliases description "
        "other_designations map_location chromosome genomic_nucleotide_accession_version "
        "start_position_on_the_genomic_accession end_position_on_the_genomic_accession "
        "orientation exon_count OMIM",
    )

    geneid2nt = {}
    with open(gene_file, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) >= 17:
                try:
                    gene_id = int(fields[2])
                    genomic_acc = fields[11].replace(".", "_")
                    entry = ntncbi(
                        tax_id=fields[0],
                        Org_name=fields[1],
                        GeneID=gene_id,
                        CurrentID=fields[3],
                        Status=fields[4],
                        Symbol=fields[5],
                        Aliases=fields[6],
                        description=fields[7],
                        other_designations=fields[8],
                        map_location=fields[9],
                        chromosome=fields[10],
                        genomic_nucleotide_accession_version=genomic_acc,
                        start_position_on_the_genomic_accession=fields[12],
                        end_position_on_the_genomic_accession=fields[13],
                        orientation=fields[14],
                        exon_count=fields[15],
                        OMIM=fields[16],
                    )
                    geneid2nt[gene_id] = entry
                except (ValueError, IndexError):
                    continue
    return geneid2nt


def initialize_go(data_dir: Path = None) -> tuple:
    """Initialize GO enrichment analysis objects.

    Loads GO ontology, gene associations, and creates enrichment study object.

    Args:
        data_dir: Directory containing GO data files.

    Returns:
        Tuple of (goeaobj, geneid_symbol_mapper) for GO enrichment analysis.
    """
    print("Initializing GO terms...")

    if data_dir is None:
        data_dir = _get_data_dir()

    gene2go_file = data_dir / "gene2go"
    gene2go_gz = data_dir / "gene2go.gz"

    if not gene2go_file.exists():
        if gene2go_gz.exists():
            print("gene2go.gz found, unzipping...")
        else:
            print("gene2go.gz not found, downloading from NCBI...")
            downloaded = download_ncbi_associations()
            shutil.move(downloaded, gene2go_gz)

        print("Extracting gene2go.gz...")
        with gzip.open(gene2go_gz, "rb") as f_in:
            with open(gene2go_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"gene2go extracted to {gene2go_file}")

    obo_file = data_dir / "go-basic.obo"
    if not obo_file.exists():
        print("go-basic.obo not found, downloading...")
        downloaded = download_go_basic_obo()
        shutil.move(downloaded, obo_file)
        print(f"go-basic.obo saved to {obo_file}")

    gene_file = data_dir / "gene_result_ncbi_human_proteincoding.txt"
    if not gene_file.exists():
        raise FileNotFoundError(
            f"Gene file not found: {gene_file}. "
            "Download from NCBI Gene with filter: "
            '"9606"[Taxonomy ID] AND alive[property] AND genetype protein coding[Properties]'
        )

    obodag = GODag(str(obo_file))
    genes = Gene2GoReader(str(gene2go_file), taxids=[9606])
    ns2assoc = genes.get_ns2assc()

    geneid2nt = _parse_geneid2nt(gene_file)

    goeaobj = GOEnrichmentStudyNS(
        geneid2nt.keys(),
        ns2assoc,
        obodag,
        propagate_counts=False,
        alpha=0.05,
        methods=["fdr_bh"],
    )

    geneid_symbol_mapper = {
        geneid2nt[key].Symbol: geneid2nt[key].GeneID for key in geneid2nt
    }

    print("GO initialization complete.")
    return goeaobj, geneid_symbol_mapper


def generate_go_table(
    sigs: pd.DataFrame, goeaobj, geneid_symbol_mapper: dict
) -> pd.DataFrame:
    """Generate GO enrichment table for significant genes.

    Args:
        sigs: DataFrame with significant genes (gene symbols as index).
        goeaobj: GO enrichment study object.
        geneid_symbol_mapper: Dictionary mapping gene symbols to gene IDs.

    Returns:
        DataFrame with GO enrichment results.
    """
    sigs_list = [str(gene) for gene in sigs.index]
    sigs_ids = [
        int(geneid_symbol_mapper[gene])
        for gene in sigs_list
        if gene in geneid_symbol_mapper
    ]
    print(
        f"Mapped {len(sigs_ids)/len(sigs_list)*100:.2f}% of "
        "significantly differentially expressed gene symbols to gene IDs."
    )

    goea_results = goeaobj.run_study(sigs_ids, prt=None)
    goea_results_sig = [r for r in goea_results if r.p_fdr_bh < 0.05]

    inverted_mapping = {v: k for k, v in geneid_symbol_mapper.items()}

    go_df = pd.DataFrame(
        [
            [
                r.GO,
                r.goterm.name,
                r.goterm.namespace,
                r.p_uncorrected,
                r.p_fdr_bh,
                r.ratio_in_study[0],
                r.ratio_in_study[1],
                r.ratio_in_study[0] / r.ratio_in_study[1] if r.ratio_in_study[1] else 0,
                list(map(lambda y: inverted_mapping.get(y, str(y)), r.study_items)),
            ]
            for r in goea_results_sig
        ],
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


def run_go_analysis(
    sigs: pd.DataFrame,
    analysis_name: str,
    output_dir: Path = None,
    data_dir: Path = None,
    goeaobj=None,
    geneid_symbol_mapper: dict = None,
) -> pd.DataFrame:
    """Run GO enrichment analysis and save results to CSV.

    Args:
        sigs: DataFrame with significant genes.
        analysis_name: Name for output files.
        output_dir: Directory for output CSV.
        data_dir: Directory containing GO data files.
        goeaobj: Pre-initialized GO enrichment object (avoids re-initialization).
        geneid_symbol_mapper: Pre-built gene symbol to ID mapping.

    Returns:
        DataFrame with GO enrichment results.
    """
    if goeaobj is None or geneid_symbol_mapper is None:
        goeaobj, geneid_symbol_mapper = initialize_go(data_dir)
    go_df = generate_go_table(sigs, goeaobj, geneid_symbol_mapper)

    if output_dir is None:
        output_dir = Path(__file__).parent.parent / "results" / "go_terms"
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_path = output_dir / f"{analysis_name}_go_terms.csv"
    go_df.to_csv(csv_path)
    print(f"GO table saved: {csv_path}")

    return go_df


def merge_go_tables(
    go_files: List[Path],
    output_dir: Path = None,
    output_filename: str = "GO_merged_results.csv",
) -> pd.DataFrame:
    """Merge GO enrichment results from multiple ligand analyses.

    Args:
        go_files: List of paths to GO result CSV files.
        output_dir: Directory for output files. Defaults to results/go_terms.
        output_filename: Name for merged output file.

    Returns:
        Merged DataFrame with all GO terms.
    """
    if not go_files:
        raise ValueError("No GO files provided")

    merged_df = None

    for go_file in go_files:
        ligand_name = go_file.stem.split("_")[0]
        df = pd.read_csv(go_file, index_col=0)

        df.insert(0, "Ligand", ligand_name)

        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.concat([merged_df, df], axis=0)

    merged_df = merged_df.sort_values("fdr", ascending=True).reset_index(drop=True)

    if output_dir is None:
        output_dir = Path(__file__).parent.parent / "results" / "go_terms"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / output_filename
    merged_df.to_csv(output_path)
    print(f"Saved merged GO table to {output_path}")

    return merged_df
