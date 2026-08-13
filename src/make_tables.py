"""Assemble and format the manuscript's tables from raw pipeline CSVs.

The modelling code (MATseq.py -> src/model_training.py, src/pydeseq2.py,
src/go_term_analysis.py) writes *raw* per-ligand / per-condition CSVs of
computed numbers. This module is the single place that reshapes those raw
outputs into the final, correctly-named manuscript tables. No numbers are
produced here: only selection, joining, column naming and formatting, so a
table can be re-assembled or re-styled without re-running the pipeline.

Two groups of functions:

1. Main-text Table 2 -- ``format_table2`` reads the nested-CV summary written by
   MATseq.py at results/nested_cv/supp_nested_cv_main.csv and turns it into the
   manuscript table.

2. Supplementary Tables S1-S12 -- ``assemble_supplementary_table_{1..12}`` and the
   convenience wrapper ``assemble_supplementary_tables``. These read the raw
   DESeq2 / GO / feature-selection CSVs and emit the named Supplementary_Table_N
   files. Assembling them here (rather than by hand) fixes three long-standing
   naming problems at source:
     * S1's LPS block was unlabelled -- every per-ligand block is now suffixed,
       including ``_LPS``, and the gene column is named ``gene``.
     * S2/S4 carried a stray ``Unnamed: 0`` index -- all writes use index=False.
     * the gene identifier column was unnamed on read-back -- set explicitly.

Raw inputs consumed (relative to the results directory):
     Table 2: nested_cv/supp_nested_cv_main.csv
     S1: differential_gene_expression/<subset>/<ligand>_deseq2_results.csv
     S2: go_terms/<subset>/<ligand>_go_terms.csv   (per-ligand GO)
     S3: fs_de_genesets/selected_vs_de_overlap_table.csv  (gene, in_de),
         ordered by feature-selection importance rank
     S4: go_terms/de_intersect_fs_go_terms.csv, go_terms/fs_only_go_terms.csv
     S5-S6: ../data/supplementary_data/Supplementary_Table_{5,6}.csv
     S7: feature_selection/mutual_information.csv
     S8: feature_selection/forest_kmeans.csv
     S9-S10: nested_cv/supp_nested_cv_{main,no_flapa}.csv
     S11-S12: validation/external_validation{,_no_flapa}_performance.csv
"""
from pathlib import Path

import pandas as pd

# ----------------------------------------------------------------------------
# Supplementary table configuration
# ----------------------------------------------------------------------------
# Column order within each per-ligand DESeq2 block (DESeq2 native order).
DESEQ2_STATS = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]

# Display order of the ligand blocks across the wide S1 table. Each entry is
# (subset_dir, file_key, column_suffix): file_key is the stem in
# "<file_key>_deseq2_results.csv" inside differential_gene_expression/<subset_dir>/;
# column_suffix is what every stat column in that block is suffixed with, so the
# LPS block is explicitly labelled (previously it was left unsuffixed). LPS is
# taken from the train_ligands run; bacterial aliases map to HKEB/HKSA.
S1_LIGAND_ORDER = [
    ("train_ligands", "Fla-PA", "Fla-PA"),
    ("train_ligands", "LPS", "LPS"),
    ("train_ligands", "PGN", "PGN"),
    ("train_ligands", "R848", "R848"),
    ("train_ligands", "Pam3", "Pam3"),
    ("additional_ligands", "LTA", "LTA"),
    ("additional_ligands", "MPLA", "MPLA"),
    ("additional_ligands", "Pam2", "Pam2"),
    ("bacterial_ligands", "HK E.coli", "HKEB"),
    ("bacterial_ligands", "HK S.aureus", "HKSA"),
]

# Per-ligand GO files merged into S2 (the five main-ligand classes), given as
# (subset_dir, file_key, ligand_label).
S2_LIGAND_ORDER = [
    ("train_ligands", "Fla-PA", "Fla-PA"),
    ("train_ligands", "LPS", "LPS"),
    ("train_ligands", "PGN", "PGN"),
    ("train_ligands", "R848", "R848"),
    ("train_ligands", "Pam3", "Pam3"),
]

# Human-readable labels for the S4 by-condition GO table, mapped to the raw
# GO result file stem each is built from.
S4_CONDITIONS = [
    ("Intersect Feature Selection and Differentially Expressed Genes",
     "de_intersect_fs_go_terms.csv"),
    ("Feature Selection Genes not Differentially Expressed",
     "fs_only_go_terms.csv"),
]

# GO tables share this column layout (the leading index is dropped on write).
GO_COLUMNS = ["GO", "term", "class", "raw_pvalue", "fdr",
              "n_genes", "n_study", "ratio_in_study", "gene_symbols"]

# Display order and human-readable labels for legacy condition-shaped inputs.
CONDITION_ORDER = ["all_genes", "feature_selection", "fs_plus_de", "random_selected"]
CONDITION_LABELS = {
    "all_genes": "All genes",
    "feature_selection": "Feature selection",
    "fs_plus_de": "Feature selection + DE",
    "random_selected": "Randomly selected",
}
MODEL_ORDER = ["LinearSVC", "SGDClassifier", "LogisticRegression",
               "RandomForest", "XGBoost"]
METRICS = [("accuracy", "Accuracy"), ("f1", "F1"),
           ("precision", "Precision"), ("recall", "Recall")]
TABLE2_COLUMNS = ["Condition", "Model", "Accuracy", "F1", "Precision", "Recall",
                  "Training time"]
SUPP_TABLE_7_COLUMNS = ["gene_rank", "MI_seed_1043", "MI_seed_3669",
                        "MI_seed_6360", "MI_mean"]
SUPP_TABLE_8_COLUMNS = ["n_selected_genes", "ARI_seed_1043", "ARI_seed_3669",
                        "ARI_seed_6360", "ARI_mean", "ARI_std"]
SUPP_TABLE_9_COLUMNS = ["model", "accuracy_mean", "accuracy_std",
                        "balanced_accuracy_mean", "balanced_accuracy_std",
                        "f1_mean", "f1_std", "f1_weighted_mean",
                        "f1_weighted_std", "pooled_f1"]
SUPP_TABLE_10_COLUMNS = ["model", "accuracy_mean", "accuracy_std",
                         "balanced_accuracy_mean", "balanced_accuracy_std",
                         "precision_mean", "precision_std", "recall_mean",
                         "recall_std", "f1_mean", "f1_std",
                         "f1_weighted_mean", "f1_weighted_std",
                         "inner_best_score_mean", "inner_best_score_std",
                         "pooled_accuracy", "pooled_balanced_accuracy",
                         "pooled_precision", "pooled_recall", "pooled_f1",
                         "pooled_f1_weighted"]
SUPP_TABLE_VALIDATION_COLUMNS = ["gene_set", "model", "accuracy",
                                 "balanced_accuracy", "precision", "recall",
                                 "f1", "f1_weighted"]


def _cell(mean, std):
    return f"{mean:.3f} ± {std:.3f}"


def _condition_label(cond, n_genes):
    base = CONDITION_LABELS.get(cond, cond)
    if n_genes is None:
        return base
    return f"{base} ({int(round(n_genes)):,} genes)"


def format_table2(
    raw_csv,
    output_dir=None,
    bold_scope="column",
):
    """Format the MATseq.py nested-CV summary into Table 2.

    Args:
        raw_csv: path to results/nested_cv/supp_nested_cv_main.csv.
        output_dir: where to write outputs (defaults to the CSV's directory).
        bold_scope: 'column' bolds the single best mean per metric across the
            table; for legacy condition-shaped input, 'group' bolds the best per
            metric within each condition.

    Writes:
        table2_formatted.csv  -- plain 'mean ± SD' cells, ** around bolded cells.
        Table_2.xlsx          -- Excel table with the same columns as the paper example.
        table2.tex            -- booktabs LaTeX, \textbf{} on bolded cells.
    Returns the formatted DataFrame that backs both outputs.
    """
    raw_csv = Path(raw_csv)
    df = pd.read_csv(raw_csv)
    if output_dir:
        output_dir = Path(output_dir)
    elif raw_csv.parent.name == "nested_cv":
        output_dir = raw_csv.parent.parent / "tables"
    else:
        output_dir = raw_csv.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    required = ["model"]
    for key, _ in METRICS:
        required.extend([f"{key}_mean", f"{key}_std"])
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Table 2 input is missing columns: {missing}")

    has_condition = "condition" in df.columns
    if not has_condition:
        df.insert(0, "condition", "feature_selection")
        has_condition = True
    model_categories = [m for m in MODEL_ORDER if m in set(df["model"])]
    model_categories.extend(sorted(set(df["model"]) - set(model_categories)))
    df["model"] = pd.Categorical(df["model"], model_categories, ordered=True)
    if has_condition:
        df["condition"] = pd.Categorical(
            df["condition"], CONDITION_ORDER, ordered=True
        )
        df = df.sort_values(["condition", "model"]).reset_index(drop=True)
    else:
        df = df.sort_values("model").reset_index(drop=True)

    best_mask = {}
    for key, _ in METRICS:
        col = f"{key}_mean"
        if bold_scope == "group" and has_condition:
            best_mask[key] = df.groupby("condition", observed=True)[col].transform(
                "max"
            )
        else:
            best_mask[key] = pd.Series(df[col].max(), index=df.index)

    rows = []
    for _, r in df.iterrows():
        row = {"Model": r["model"]}
        if has_condition:
            row = {
                "Condition": _condition_label(r["condition"], r.get("n_genes_mean")),
                "Model": r["model"],
            }
        for key, header in METRICS:
            txt = _cell(r[f"{key}_mean"], r[f"{key}_std"])
            is_best = abs(r[f"{key}_mean"] - best_mask[key].iloc[r.name]) < 1e-9
            row[header] = f"**{txt}**" if is_best else txt
            row[f"_{key}_best"] = bool(is_best)
        rows.append(row)
    fmt = pd.DataFrame(rows)

    fmt["Training time"] = ""
    time_cols = ["Training time", "training_time", "training_time_seconds"]
    for time_col in time_cols:
        if time_col in df.columns:
            fmt["Training time"] = df[time_col].fillna("").astype(str).values
            break

    display = fmt.copy()
    display["Condition"] = display["Condition"].mask(
        display["Condition"].duplicated(), ""
    )
    display[TABLE2_COLUMNS].to_csv(output_dir / "table2_formatted.csv", index=False)
    display[TABLE2_COLUMNS].to_excel(output_dir / "Table_2.xlsx", index=False)

    _write_latex(fmt, output_dir / "table2.tex")
    return fmt


def _write_latex(fmt, out_path):
    headers = [h for _, h in METRICS] + ["Training time"]
    has_condition = "Condition" in fmt.columns
    lines = [
        r"\begin{table}[t]",
        r"\centering",
        r"\caption{Classification performance from MATseq.py nested "
        r"cross-validation. Cells show mean $\pm$ standard deviation across "
        r"the five outer folds; the best value per metric is in bold.}",
        r"\label{tab:table2}",
        r"\begin{tabular}{" + ("ll" if has_condition else "l")
        + "c" * len(headers) + "}",
        r"\toprule",
        (("Condition & " if has_condition else "") + "Model & "
         + " & ".join(headers) + r" \\"),
        r"\midrule",
    ]
    prev_cond = None
    for _, r in fmt.iterrows():
        cells = []
        for _, header in METRICS:
            v = r[header]
            if v.startswith("**") and v.endswith("**"):
                v = r"\textbf{" + v.strip("*") + "}"
            cells.append(v.replace("±", r"$\pm$"))
        cells.append(str(r.get("Training time", "")))
        if has_condition:
            if prev_cond is not None and r["Condition"] != prev_cond:
                lines.append(r"\midrule")
            cond = r["Condition"] if r["Condition"] != prev_cond else ""
            lines.append(
                f"{cond} & {r['Model']} & " + " & ".join(cells) + r" \\"
            )
            prev_cond = r["Condition"]
        else:
            lines.append(f"{r['Model']} & " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    Path(out_path).write_text("\n".join(lines) + "\n")


# ============================================================================
# Supplementary tables S1-S12
# ============================================================================
def _read_go_table(path):
    """Read a raw GO CSV, dropping a leading unnamed index column if present."""
    df = pd.read_csv(path)
    unnamed = [c for c in df.columns if str(c).startswith("Unnamed:")]
    if unnamed:
        df = df.drop(columns=unnamed)
    return df


def assemble_supplementary_table_1(results_dir, output_dir=None,
                                   filename="Supplementary_Table_1.csv"):
    """Assemble the wide per-ligand DESeq2 table (S1).

    Reads differential_gene_expression/<subset>/<ligand>_deseq2_results.csv for
    each ligand in S1_LIGAND_ORDER, suffixes every statistic column with the
    ligand label (so the LPS block is no longer unlabelled), and outer-joins the
    blocks on the gene identifier. The gene column is named ``gene``.

    Returns the assembled DataFrame.
    """
    results_dir = Path(results_dir)
    de_dir = results_dir / "differential_gene_expression"
    merged = None
    for subset, file_key, suffix in S1_LIGAND_ORDER:
        fpath = de_dir / subset / f"{file_key}_deseq2_results.csv"
        if not fpath.exists():
            raise FileNotFoundError(f"S1: missing DESeq2 file {fpath}")
        block = pd.read_csv(fpath)
        # The raw file's first column is the gene id (header 'Geneid'); make it
        # a named 'gene' column and suffix the six statistic columns.
        block = block.rename(columns={block.columns[0]: "gene"})
        block = block.rename(columns={c: f"{c}_{suffix}" for c in DESEQ2_STATS
                                      if c in block.columns})
        merged = block if merged is None else merged.merge(block, on="gene", how="outer")

    out_dir = Path(output_dir) if output_dir else results_dir / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_dir / filename, index=False)
    return merged


def assemble_supplementary_table_2(results_dir, output_dir=None,
                                   filename="Supplementary_Table_2.csv"):
    """Assemble the per-ligand GO enrichment table (S2).

    Concatenates the five main-ligand GO tables (S2_LIGAND_ORDER) with a leading
    ``Ligand`` column, sorted by FDR. No stray index column is written.

    Returns the assembled DataFrame.
    """
    results_dir = Path(results_dir)
    go_dir = results_dir / "go_terms"
    frames = []
    for subset, file_key, label in S2_LIGAND_ORDER:
        fpath = go_dir / subset / f"{file_key}_go_terms.csv"
        if not fpath.exists():
            raise FileNotFoundError(f"S2: missing GO file {fpath}")
        df = _read_go_table(fpath)
        df.insert(0, "Ligand", label)
        frames.append(df)
    merged = pd.concat(frames, axis=0, ignore_index=True)
    merged = merged.sort_values("fdr", ascending=True).reset_index(drop=True)

    out_dir = Path(output_dir) if output_dir else results_dir / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_dir / filename, index=False)
    return merged


def assemble_supplementary_table_3(results_dir, output_dir=None,
                                   filename="Supplementary_Table_3.csv",
                                   include_rank=True):
    """Assemble the selected-genes table (S3).

    Reads fs_de_genesets/selected_vs_de_overlap_table.csv (columns
    ``gene``, ``in_de`` and, from the current pipeline, ``rank`` and
    ``importance``) and emits ``Gene`` and ``Differentially_Expressed``.
    When the source carries a ``rank`` column the rows are ordered by
    ExtraTrees feature importance (rank 1 = most important) and the optional
    ``Rank`` column reflects that true importance rank; otherwise Rank falls
    back to row order. Set include_rank=False to omit the Rank column.

    Returns the assembled DataFrame.
    """
    results_dir = Path(results_dir)
    fpath = results_dir / "fs_de_genesets" / "selected_vs_de_overlap_table.csv"
    if not fpath.exists():
        legacy_path = results_dir / "feature_selection" / "selected_vs_de_overlap_table.csv"
        fpath = legacy_path
    if not fpath.exists():
        raise FileNotFoundError(f"S3: missing feature-selection file {fpath}")
    src = pd.read_csv(fpath)
    if "rank" in src.columns:
        src = src.sort_values("rank").reset_index(drop=True)
    out = pd.DataFrame({
        "Gene": src["gene"],
        "Differentially_Expressed": src["in_de"].astype(bool),
    })
    if include_rank:
        if "rank" in src.columns:
            out.insert(0, "Rank", src["rank"].astype(int).values)
        else:
            out.insert(0, "Rank", range(1, len(out) + 1))

    out_dir = Path(output_dir) if output_dir else results_dir / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / filename, index=False)
    return out


def assemble_supplementary_table_4(results_dir, output_dir=None,
                                   filename="Supplementary_Table_4.csv"):
    """Assemble the by-condition GO enrichment table (S4).

    Concatenates the two gene-set GO tables (S4_CONDITIONS) with a leading
    ``Condition`` column. No stray index column is written.

    Returns the assembled DataFrame.
    """
    results_dir = Path(results_dir)
    go_dir = results_dir / "go_terms"
    frames = []
    for label, fname in S4_CONDITIONS:
        fpath = go_dir / fname
        if not fpath.exists():
            raise FileNotFoundError(f"S4: missing GO file {fpath}")
        df = _read_go_table(fpath)
        df.insert(0, "Condition", label)
        frames.append(df)
    merged = pd.concat(frames, axis=0, ignore_index=True)

    out_dir = Path(output_dir) if output_dir else results_dir / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_dir / filename, index=False)
    return merged


def assemble_supplementary_table_5(results_dir, output_dir=None,
                                   filename="Supplementary_Table_5.csv"):
    """Assemble the TLR4 reporter assay source table (S5)."""
    results_dir = Path(results_dir)
    fpath = results_dir.parent / "data" / "supplementary_data" / filename
    if not fpath.exists():
        fpath = Path("data") / "supplementary_data" / filename
    if not fpath.exists():
        raise FileNotFoundError(f"S5: missing TLR table {fpath}")
    out = pd.read_csv(fpath).rename(
        columns={"Concentration_(EU_mL)": "Concentration_LPS_(eq_mL)"}
    )
    out_dir = Path(output_dir) if output_dir else results_dir / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / filename, index=False)
    return out


def assemble_supplementary_table_6(results_dir, output_dir=None,
                                   filename="Supplementary_Table_6.csv"):
    """Assemble the TLR2 reporter assay source table (S6)."""
    results_dir = Path(results_dir)
    fpath = results_dir.parent / "data" / "supplementary_data" / filename
    if not fpath.exists():
        fpath = Path("data") / "supplementary_data" / filename
    if not fpath.exists():
        raise FileNotFoundError(f"S6: missing TLR table {fpath}")
    out = pd.read_csv(fpath).rename(
        columns={"Concentration_(ng_mL)": "Concentration_Pam3_(ng_mL)"}
    )
    out_dir = Path(output_dir) if output_dir else results_dir / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / filename, index=False)
    return out


def assemble_supplementary_table_7(results_dir, output_dir=None,
                                   filename="Supplementary_Table_7.csv"):
    """Assemble the mutual-information curve table (S7)."""
    results_dir = Path(results_dir)
    fpath = results_dir / "feature_selection" / "mutual_information.csv"
    if not fpath.exists():
        raise FileNotFoundError(f"S7: missing mutual-information file {fpath}")
    src = pd.read_csv(fpath)
    out = pd.DataFrame(index=src.index)
    out["gene_rank"] = src["rank"] if "rank" in src.columns else range(1, len(src) + 1)
    for col in SUPP_TABLE_7_COLUMNS[1:4]:
        raw_col = col.replace("MI_seed_", "mi_seed_")
        out[col] = src[raw_col] if raw_col in src.columns else pd.NA
    out["MI_mean"] = src["mi_sorted"] if "mi_sorted" in src.columns else src.get("MI_mean")
    out = out[SUPP_TABLE_7_COLUMNS]
    out_dir = Path(output_dir) if output_dir else results_dir / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / filename, index=False)
    return out


def assemble_supplementary_table_8(results_dir, output_dir=None,
                                   filename="Supplementary_Table_8.csv"):
    """Assemble the forest/KMeans ARI scan table (S8)."""
    results_dir = Path(results_dir)
    fpath = results_dir / "feature_selection" / "forest_kmeans.csv"
    if not fpath.exists():
        raise FileNotFoundError(f"S8: missing forest/KMeans file {fpath}")
    src = pd.read_csv(fpath)
    if {"n_estimators", "max_depth", "ari_mean"}.issubset(src.columns):
        best = src.loc[src["ari_mean"].idxmax()]
        src = src[
            (src["n_estimators"] == best["n_estimators"])
            & (src["max_depth"] == best["max_depth"])
        ]
    src = src.sort_values("n_selected").reset_index(drop=True)
    out = pd.DataFrame({"n_selected_genes": src["n_selected"]})
    for col in SUPP_TABLE_8_COLUMNS[1:4]:
        raw_col = col.replace("ARI_seed_", "ari_seed_")
        out[col] = src[raw_col] if raw_col in src.columns else pd.NA
    out["ARI_mean"] = src["ari_mean"]
    out["ARI_std"] = src["ari_std"]
    out = out[SUPP_TABLE_8_COLUMNS]
    out_dir = Path(output_dir) if output_dir else results_dir / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / filename, index=False)
    return out


def assemble_supplementary_table_9(results_dir, output_dir=None,
                                   filename="Supplementary_Table_9.csv"):
    """Assemble the main-panel nested-CV summary table (S9)."""
    results_dir = Path(results_dir)
    fpath = results_dir / "nested_cv" / "supp_nested_cv_main.csv"
    if not fpath.exists():
        fpath = results_dir / "tables" / "supp_nested_cv_main.csv"
    if not fpath.exists():
        raise FileNotFoundError(f"S9: missing nested-CV summary {fpath}")
    src = pd.read_csv(fpath)
    out = src.reindex(columns=SUPP_TABLE_9_COLUMNS)
    out_dir = Path(output_dir) if output_dir else results_dir / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / filename, index=False)
    return out


def assemble_supplementary_table_10(results_dir, output_dir=None,
                                    filename="Supplementary_Table_10.csv"):
    """Assemble the no-Fla-PA nested-CV summary table (S10)."""
    results_dir = Path(results_dir)
    fpath = results_dir / "nested_cv" / "supp_nested_cv_no_flapa.csv"
    if not fpath.exists():
        fpath = results_dir / "tables" / "supp_nested_cv_no_flapa.csv"
    if not fpath.exists():
        raise FileNotFoundError(f"S10: missing nested-CV summary {fpath}")
    src = pd.read_csv(fpath)
    out = src.reindex(columns=SUPP_TABLE_10_COLUMNS)
    out_dir = Path(output_dir) if output_dir else results_dir / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / filename, index=False)
    return out


def assemble_supplementary_table_11(results_dir, output_dir=None,
                                    filename="Supplementary_Table_11.csv"):
    """Assemble the external-validation performance table (S11)."""
    results_dir = Path(results_dir)
    fpath = results_dir / "validation" / "external_validation_performance.csv"
    if not fpath.exists():
        fpath = results_dir / "tables" / "external_validation_performance.csv"
    if not fpath.exists():
        raise FileNotFoundError(f"S11: missing validation summary {fpath}")
    src = pd.read_csv(fpath)
    out = src.reindex(columns=SUPP_TABLE_VALIDATION_COLUMNS)
    out_dir = Path(output_dir) if output_dir else results_dir / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / filename, index=False)
    return out


def assemble_supplementary_table_12(results_dir, output_dir=None,
                                    filename="Supplementary_Table_12.csv"):
    """Assemble the no-Fla-PA external-validation performance table (S12)."""
    results_dir = Path(results_dir)
    fpath = results_dir / "validation" / "external_validation_no_flapa_performance.csv"
    if not fpath.exists():
        fpath = results_dir / "tables" / "external_validation_no_flapa_performance.csv"
    if not fpath.exists():
        raise FileNotFoundError(f"S12: missing validation summary {fpath}")
    src = pd.read_csv(fpath)
    out = src.reindex(columns=SUPP_TABLE_VALIDATION_COLUMNS)
    out_dir = Path(output_dir) if output_dir else results_dir / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / filename, index=False)
    return out


def assemble_supplementary_tables(results_dir, output_dir=None):
    """Regenerate Supplementary Tables S1-S12 from raw pipeline outputs.

    Returns a dict mapping table name to the assembled DataFrame.
    """
    return {
        "S1": assemble_supplementary_table_1(results_dir, output_dir),
        "S2": assemble_supplementary_table_2(results_dir, output_dir),
        "S3": assemble_supplementary_table_3(results_dir, output_dir),
        "S4": assemble_supplementary_table_4(results_dir, output_dir),
        "S5": assemble_supplementary_table_5(results_dir, output_dir),
        "S6": assemble_supplementary_table_6(results_dir, output_dir),
        "S7": assemble_supplementary_table_7(results_dir, output_dir),
        "S8": assemble_supplementary_table_8(results_dir, output_dir),
        "S9": assemble_supplementary_table_9(results_dir, output_dir),
        "S10": assemble_supplementary_table_10(results_dir, output_dir),
        "S11": assemble_supplementary_table_11(results_dir, output_dir),
        "S12": assemble_supplementary_table_12(results_dir, output_dir),
    }


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "supp":
        results_dir = sys.argv[2] if len(sys.argv) > 2 else "results"
        tables = assemble_supplementary_tables(results_dir)
        for name, df in tables.items():
            print(f"{name}: {df.shape}  cols={list(df.columns)[:6]}")
    else:
        raw = sys.argv[1] if len(sys.argv) > 1 else "results/nested_cv/supp_nested_cv_main.csv"
        out = format_table2(raw)
        print(out[TABLE2_COLUMNS].to_string(index=False))
