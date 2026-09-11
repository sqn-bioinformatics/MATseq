"""Assemble the per-panel pipeline figures into the manuscript composite pages.

Run as MATseq.py step 10, once every panel figure (DESeq2, PCA, GO, venn,
prediction and validation heatmaps) has been written. No numbers are produced
here: panels are only cropped, ordered and lettered.

Model order and the external-batch ligand order come from the config
(``hyperparameter_grids``, ``class_order_for_plotting``); the per-ligand DESeq2
pages follow ``make_tables.S1_LIGAND_ORDER``, so Supplementary Figure 1 and
Supplementary Table S1 present their blocks in the same order. Most pages hold
at most four panels so they drop cleanly onto
one Word page; the model-comparison pages carry PCA plus all five models. Panel
white margins are auto-cropped to remove inter-panel whitespace, panel letters
run in sequence across the split pages of a figure, and missing panels are
flagged in place rather than breaking the run.
"""

from __future__ import annotations

from pathlib import Path
from string import ascii_uppercase

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import FancyArrow, FancyBboxPatch

from .config import (
    CLASS_ORDER,
    FEATURE_SELECTION_CONFIG,
    HYPERPARAMETER_GRIDS,
    primary_geneset_name,
)
from .feature_engineering import SEEDS
from .make_tables import S1_LIGAND_ORDER

MODELS = list(HYPERPARAMETER_GRIDS)
MODEL_TITLES = {
    "LinearSVC": "LinearSVC",
    "SGDClassifier": "SGD",
    "LogisticRegression": "Logistic Regression",
    "RandomForest": "Random Forest",
    "XGBoost": "XGBoost",
}


def _load_cropped(path: Path) -> np.ndarray:
    """Read an image and trim its uniform white border."""
    pad = 6
    img = mpimg.imread(str(path))
    a = img.astype(np.float32)
    if a.max() > 1.0:
        a /= 255.0
    rgb = a[..., :3] if a.ndim == 3 else np.dstack([a] * 3)
    nonwhite = (rgb < 0.97).any(axis=2)
    if a.ndim == 3 and a.shape[2] == 4:
        nonwhite &= a[..., 3] > 0.05
    if not nonwhite.any():
        return img
    rows = np.where(nonwhite.any(axis=1))[0]
    cols = np.where(nonwhite.any(axis=0))[0]
    r0, r1 = max(0, rows[0] - pad), min(img.shape[0], rows[-1] + 1 + pad)
    c0, c1 = max(0, cols[0] - pad), min(img.shape[1], cols[-1] + 1 + pad)
    return img[r0:r1, c0:c1]


def _aspect(path: Path) -> float:
    if not Path(path).is_file():
        return 1.0
    h, w = _load_cropped(path).shape[:2]
    return w / h


def _place_image(ax: Axes, path: Path, title: str | None = None) -> None:
    ax.axis("off")
    if not Path(path).is_file():
        ax.text(0.5, 0.5, f"missing:\n{Path(path).name}", ha="center", va="center",
                fontsize=7, color="0.6", transform=ax.transAxes)
        return
    ax.imshow(_load_cropped(path))
    if title:
        ax.set_title(title, fontsize=9)


def _letter(ax: Axes, letter: str | None) -> None:
    if letter:
        ax.text(-0.01, 1.01, letter, transform=ax.transAxes, fontsize=15,
                fontweight="bold", ha="right", va="bottom")


def _save(fig: Figure, out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=250, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"  Composite page saved: {out_path}")
    return out_path


def _grid_page(paths: list[Path], titles: list[str | None], letters: list[str],
               out_path: Path, ncols: int, height: float = 4.4) -> Path:
    """Lay panels out row-major; column widths track the widest cropped panel."""
    nrows = -(-len(paths) // ncols)
    wr = [max(_aspect(p) for p in paths[c::ncols]) for c in range(ncols)]
    fig = plt.figure(figsize=(height * sum(wr), height * nrows))
    gs = fig.add_gridspec(nrows, ncols, width_ratios=wr, wspace=0.04, hspace=0.12)
    for i, (p, t, l) in enumerate(zip(paths, titles, letters)):
        r, c = divmod(i, ncols)
        ax = fig.add_subplot(gs[r, c])
        _place_image(ax, p, t)
        _letter(ax, l)
    return _save(fig, out_path)


def _draw_pipeline_schematic(ax: Axes, results_dir: Path) -> None:
    """Feature-selection flow chart, with the gene counts of the current run."""
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    overlap = pd.read_csv(
        results_dir / "fs_de_genesets" / "selected_vs_de_overlap_table.csv"
    )
    steps = [
        "Library-size\nnormalisation\n+ log1p, z-score",
        f"Mutual information\nranking over\n{len(SEEDS)} seeds",
        f"MI elbow\n({FEATURE_SELECTION_CONFIG['k_best']:,} genes)",
        f"ExtraTrees ranking\n({FEATURE_SELECTION_CONFIG['n_estimators']:,} trees);"
        "\nk-means ARI per rank",
        f"{FEATURE_SELECTION_CONFIG['max_features']} most\nimportant genes",
        f"Overlap with DE\n({int(overlap['in_de'].sum())}-gene set)",
    ]
    n = len(steps)
    gap = 0.012
    w = (1.0 - gap * (n - 1)) / n
    y, h = 0.30, 0.40
    for i, text in enumerate(steps):
        x = i * (w + gap)
        ax.add_patch(FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.01",
            linewidth=1.1, edgecolor="#2c5f8a",
            facecolor="#e8f0f7" if i < n - 1 else "#dcecdc"))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=7.5)
        if i < n - 1:
            ax.add_patch(FancyArrow(
                x + w, y + h / 2, gap, 0, width=0.004, head_width=0.05,
                head_length=gap * 0.9, length_includes_head=True, color="#2c5f8a"))


def _model_grid_page(pca_png: Path, pca_title: str, pred_dir: Path,
                     letters: list[str], out_path: Path) -> Path:
    """PCA plus every model's probability heatmap on one page."""
    return _grid_page(
        [pca_png] + [pred_dir / f"{m}_probabilities_heatmap.png" for m in MODELS],
        [pca_title] + [MODEL_TITLES.get(m, m) for m in MODELS],
        letters, out_path, ncols=3,
    )


def _prediction_pages(base: str, results_dir: Path, pred_root: Path,
                      out_dir: Path) -> list[Path]:
    """One page per prediction subset; each shows PCA plus every model."""
    specs = [
        ("unseen ligands", "additional_ligands"),
        ("heat-killed bacteria", "bacterial_ligands"),
    ]
    n_panels = len(MODELS) + 1
    return [
        _model_grid_page(
            results_dir / "figures" / "pca" / f"{subset}_feature_selected.png", f"PCA ({label})",
            pred_root / subset,
            list(ascii_uppercase[i * n_panels:(i + 1) * n_panels]),
            out_dir / f"{base}_p{i + 1}.png",
        )
        for i, (label, subset) in enumerate(specs)
    ]


def _paginate_de(rows: list[tuple[str, str]], base: str, results_dir: Path,
                 out_dir: Path) -> list[Path]:
    """Two ligands (four panels) per page; rows are (ligand, DESeq2 subset)."""
    chunks = [rows[i:i + 2] for i in range(0, len(rows), 2)]
    de0 = results_dir / "figures" / "deseq2" / rows[0][1]
    wr = [_aspect(de0 / f"{rows[0][0]}_volcano.png"),
          _aspect(de0 / f"{rows[0][0]}_histogram.png")]
    letters = iter(ascii_uppercase)
    outs = []
    for ci, chunk in enumerate(chunks):
        fig = plt.figure(figsize=(12, 5.4 * len(chunk)))
        gs = fig.add_gridspec(len(chunk), 2, width_ratios=wr, hspace=0.16, wspace=0.05)
        for r, (ligand, subset) in enumerate(chunk):
            de = results_dir / "figures" / "deseq2" / subset
            ax0 = fig.add_subplot(gs[r, 0])
            ax1 = fig.add_subplot(gs[r, 1])
            _place_image(ax0, de / f"{ligand}_volcano.png")
            _place_image(ax1, de / f"{ligand}_histogram.png")
            _letter(ax0, next(letters))
            _letter(ax1, next(letters))
        fname = f"{base}.png" if len(chunks) == 1 else f"{base}_p{ci + 1}.png"
        outs.append(_save(fig, out_dir / fname))
    return outs


def compose_figure2(results_dir: Path, out_dir: Path) -> list[Path]:
    """LPS DESeq2 on the training batch: A) volcano, B) clustered heatmap."""
    de = results_dir / "figures" / "deseq2" / "train_ligands"
    return [_grid_page(
        [de / "LPS_volcano.png", de / "LPS_histogram.png"],
        [None, None], ["A", "B"], out_dir / "Figure2.png", ncols=2, height=5.2)]


def compose_figure3(results_dir: Path, out_dir: Path) -> list[Path]:
    """Feature selection, split into two <=4-panel pages (letters run A-E)."""
    fig = plt.figure(figsize=(11, 7))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.55, 1], hspace=0.1, wspace=0.05)
    ax_a = fig.add_subplot(gs[0, :])
    _draw_pipeline_schematic(ax_a, results_dir)
    _letter(ax_a, "A")
    ax_b1 = fig.add_subplot(gs[1, 0])
    ax_b2 = fig.add_subplot(gs[1, 1])
    _place_image(ax_b1, results_dir / "figures" / "pca" / "train_ligands_pca.png", "before selection")
    _place_image(ax_b2, results_dir / "figures" / "pca" / "train_ligands_feature_selected.png",
                 "after selection")
    _letter(ax_b1, "B")
    p1 = _save(fig, out_dir / "Figure3_p1.png")

    fig = plt.figure(figsize=(11, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.12, wspace=0.05)
    ax_c = fig.add_subplot(gs[0, 0])
    ax_d = fig.add_subplot(gs[0, 1])
    ax_e = fig.add_subplot(gs[1, :])
    _place_image(ax_c, results_dir / "figures" / "venn" / "venn_de_vs_fs.png")
    _place_image(ax_d, results_dir / "figures" / "go" / "de_intersect_fs_go.png")
    _place_image(ax_e, results_dir / "figures" / "go" / "fs_only_go.png")
    _letter(ax_c, "C")
    _letter(ax_d, "D")
    _letter(ax_e, "E")
    p2 = _save(fig, out_dir / "Figure3_p2.png")
    return [p1, p2]


def compose_figure4(results_dir: Path, out_dir: Path) -> list[Path]:
    """Exploratory predictions on unseen ligands and heat-killed bacteria."""
    return _prediction_pages(
        "Figure4", results_dir, results_dir / "predictions" / primary_geneset_name(),
        out_dir,
    )


def compose_figure4a(results_dir: Path, out_dir: Path) -> list[Path]:
    """Figure 4 composition applied to the external test batch (all models)."""
    val = (results_dir / "validation" / "test_set" / primary_geneset_name()
           / "test_ligands")
    return [_model_grid_page(
        results_dir / "figures" / "pca" / "test_ligands_feature_selected.png",
        "PCA (external test batch)", val,
        list(ascii_uppercase[:len(MODELS) + 1]),
        out_dir / "Figure4a_external_test.png")]


def compose_supp_figure1(results_dir: Path, out_dir: Path) -> list[Path]:
    """DESeq2 volcano and heatmap for every stimulus except LPS (Figure 2)."""
    rows = [
        (ligand, subset)
        for subset, ligand, _ in S1_LIGAND_ORDER
        if ligand != "LPS"
    ]
    return _paginate_de(rows, "SupplementaryFigure1", results_dir, out_dir)


def compose_supp_figure2(results_dir: Path, out_dir: Path) -> list[Path]:
    """HEK-Blue TLR2/TLR4 reporter dose-response."""
    fig = plt.figure(figsize=(8, 9))
    ax = fig.add_subplot(111)
    _place_image(ax, results_dir / "figures" / "supplementary" / "tlr_hek_blue.png")
    return [_save(fig, out_dir / "SupplementaryFigure2.png")]


def compose_supp_figure3(results_dir: Path, out_dir: Path) -> list[Path]:
    """Figure 4 composition for the models trained without Fla-PA."""
    return _prediction_pages(
        "SupplementaryFigure3",
        results_dir,
        results_dir / "predictions" / "no_flapa" / primary_geneset_name(),
        out_dir,
    )


def compose_supp_figure4_external_de(results_dir: Path, out_dir: Path) -> list[Path]:
    """DESeq2 volcano and heatmap for the external test batch."""
    rows = [
        (ligand, "test_ligands")
        for ligand in CLASS_ORDER["test_ligands"]
        if ligand != "negative_control"
    ]
    return _paginate_de(rows, "SupplementaryFigure4_external_DE", results_dir, out_dir)


def compose_supp_figure5(results_dir: Path, out_dir: Path) -> list[Path]:
    """Gene-number selection: A) MI elbow curve, B) forest/k-means ARI sweep."""
    fs = results_dir / "figures" / "feature_selection"
    return [_grid_page(
        [fs / "mutual_information.png", fs / "forest_ari_sweep.png"],
        [None, None], ["A", "B"], out_dir / "SupplementaryFigure5.png",
        ncols=2, height=4.6)]


def compose_figures(results_dir: Path, out_dir: Path) -> list[Path]:
    """Compose every manuscript figure page; returns the written PNG paths."""
    builders = [
        compose_figure2, compose_figure3, compose_figure4, compose_figure4a,
        compose_supp_figure1, compose_supp_figure2, compose_supp_figure3,
        compose_supp_figure4_external_de, compose_supp_figure5,
    ]
    paths = []
    for build in builders:
        paths.extend(build(results_dir, out_dir))
    return paths
