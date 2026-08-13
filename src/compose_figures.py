"""Assemble the per-panel pipeline figures into the manuscript composite pages.

Run as MATseq.py step 10, once every panel figure (DESeq2, PCA, GO, venn,
prediction and validation heatmaps) has been written. No numbers are produced
here: panels are only cropped, ordered and lettered.

Ligand and model order come from the config (``class_order_for_plotting``,
``ligands`` and ``hyperparameter_grids``), so a change there propagates to the
composed pages. Most pages hold at most four panels so they drop cleanly onto
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
    ADDITIONAL_LIGANDS,
    BACTERIAL_LIGANDS,
    CLASS_ORDER,
    FEATURE_SELECTION_CONFIG,
    HYPERPARAMETER_GRIDS,
    MAIN_LIGANDS,
    primary_geneset_name,
)
from .feature_engineering import SEEDS

PROJECT_ROOT = Path(__file__).parent.parent
RESULTS = PROJECT_ROOT / "results"
FIG = RESULTS / "figures"
OUT = PROJECT_ROOT / "paper" / "paper_updated" / "figures"

MODELS = list(HYPERPARAMETER_GRIDS)
MODEL_TITLES = {
    "LinearSVC": "LinearSVC",
    "SGDClassifier": "SGD",
    "LogisticRegression": "Logistic Regression",
    "RandomForest": "Random Forest",
    "XGBoost": "XGBoost",
}
DE_SUBSET_LIGANDS = {
    "train_ligands": MAIN_LIGANDS,
    "test_ligands": MAIN_LIGANDS,
    "additional_ligands": ADDITIONAL_LIGANDS,
    "bacterial_ligands": BACTERIAL_LIGANDS,
}


def _load_cropped(path: Path, pad: int = 6) -> np.ndarray:
    """Read an image and trim its uniform white border."""
    img = mpimg.imread(str(path))
    a = img.astype(float)
    if a.max() > 1.0:
        a = a / 255.0
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


def _row_page(paths: list[Path], titles: list[str | None], letters: list[str | None],
              out_path: Path, height: float = 4.6) -> Path:
    """One tight row; column widths track each cropped image so panels abut."""
    wr = [_aspect(p) for p in paths]
    fig = plt.figure(figsize=(height * sum(wr), height))
    gs = fig.add_gridspec(1, len(paths), width_ratios=wr, wspace=0.03)
    for i, (p, t, l) in enumerate(zip(paths, titles, letters)):
        ax = fig.add_subplot(gs[0, i])
        _place_image(ax, p, t)
        _letter(ax, l)
    return _save(fig, out_path)


def _draw_pipeline_schematic(ax: Axes) -> None:
    """Feature-selection flow chart, with the gene counts of the current run."""
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    overlap = pd.read_csv(
        RESULTS / "fs_de_genesets" / "selected_vs_de_overlap_table.csv"
    )
    steps = [
        "Library-size\nnormalisation\n+ log1p, z-score",
        f"Mutual information\nranking over\n{len(SEEDS)} seeds",
        f"MI elbow\n({FEATURE_SELECTION_CONFIG['k_best']:,} genes)",
        "ExtraTrees / k-means\nARI scan over\nforest settings",
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
                     letters: list[str], out_path: Path, ncols: int = 3,
                     height: float = 4.4) -> Path:
    """PCA plus every model's probability heatmap on one page (2 rows x ncols)."""
    paths = [pca_png] + [pred_dir / f"{m}_probabilities_heatmap.png" for m in MODELS]
    titles = [pca_title] + [MODEL_TITLES[m] for m in MODELS]
    nrows = -(-len(paths) // ncols)
    wr = []
    for c in range(ncols):
        col = [_aspect(paths[r * ncols + c]) for r in range(nrows)
               if r * ncols + c < len(paths)]
        wr.append(max(col) if col else 1.0)
    fig = plt.figure(figsize=(height * sum(wr), height * nrows))
    gs = fig.add_gridspec(nrows, ncols, width_ratios=wr, wspace=0.04, hspace=0.12)
    for i, (p, t, l) in enumerate(zip(paths, titles, letters)):
        r, c = divmod(i, ncols)
        ax = fig.add_subplot(gs[r, c])
        _place_image(ax, p, t)
        _letter(ax, l)
    return _save(fig, out_path)


def _prediction_pages(base: str, pred_root: Path, out_dir: Path) -> list[Path]:
    """One page per prediction subset; each shows PCA plus all five models."""
    specs = [
        ("unseen ligands", "additional_ligands", list("ABCDEF"), f"{base}_p1.png"),
        ("heat-killed bacteria", "bacterial_ligands", list("GHIJKL"), f"{base}_p2.png"),
    ]
    return [
        _model_grid_page(
            FIG / "pca" / f"{subset}_feature_selected.png", f"PCA ({label})",
            pred_root / subset, letters, out_dir / fname,
        )
        for label, subset, letters, fname in specs
    ]


def _paginate_de(rows: list[tuple[str, str]], base: str, out_dir: Path) -> list[Path]:
    """Two ligands (four panels) per page; rows are (ligand, DESeq2 subset)."""
    chunks = [rows[i:i + 2] for i in range(0, len(rows), 2)]
    de0 = FIG / "deseq2" / rows[0][1]
    wr = [_aspect(de0 / f"{rows[0][0]}_volcano.png"),
          _aspect(de0 / f"{rows[0][0]}_histogram.png")]
    outs = []
    for ci, chunk in enumerate(chunks):
        fig = plt.figure(figsize=(12, 5.4 * len(chunk)))
        gs = fig.add_gridspec(len(chunk), 2, width_ratios=wr, hspace=0.16, wspace=0.05)
        for r, (ligand, subset) in enumerate(chunk):
            de = FIG / "deseq2" / subset
            ax0 = fig.add_subplot(gs[r, 0])
            ax1 = fig.add_subplot(gs[r, 1])
            _place_image(ax0, de / f"{ligand}_volcano.png")
            _place_image(ax1, de / f"{ligand}_histogram.png")
            _letter(ax0, ascii_uppercase[4 * ci + 2 * r])
            _letter(ax1, ascii_uppercase[4 * ci + 2 * r + 1])
        fname = f"{base}.png" if len(chunks) == 1 else f"{base}_p{ci + 1}.png"
        outs.append(_save(fig, out_dir / fname))
    return outs


def compose_figure2(out_dir: Path = OUT) -> list[Path]:
    """LPS DESeq2 on the training batch: A) volcano, B) clustered heatmap."""
    de = FIG / "deseq2" / "train_ligands"
    return [_row_page(
        [de / "LPS_volcano.png", de / "LPS_histogram.png"],
        [None, None], ["A", "B"], out_dir / "Figure2.png", height=5.2)]


def compose_figure3(out_dir: Path = OUT) -> list[Path]:
    """Feature selection, split into two <=4-panel pages (letters run A-E)."""
    fig = plt.figure(figsize=(11, 7))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.55, 1], hspace=0.1, wspace=0.05)
    ax_a = fig.add_subplot(gs[0, :])
    _draw_pipeline_schematic(ax_a)
    _letter(ax_a, "A")
    ax_b1 = fig.add_subplot(gs[1, 0])
    ax_b2 = fig.add_subplot(gs[1, 1])
    _place_image(ax_b1, FIG / "pca" / "train_ligands_pca.png", "before selection")
    _place_image(ax_b2, FIG / "pca" / "train_ligands_feature_selected.png",
                 "after selection")
    _letter(ax_b1, "B")
    p1 = _save(fig, out_dir / "Figure3_p1.png")

    fig = plt.figure(figsize=(11, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.12, wspace=0.05)
    ax_c = fig.add_subplot(gs[0, 0])
    ax_d = fig.add_subplot(gs[0, 1])
    ax_e = fig.add_subplot(gs[1, :])
    _place_image(ax_c, FIG / "venn" / "venn_de_vs_fs.png")
    _place_image(ax_d, FIG / "go" / "de_intersect_fs_go.png")
    _place_image(ax_e, FIG / "go" / "fs_only_go.png")
    _letter(ax_c, "C")
    _letter(ax_d, "D")
    _letter(ax_e, "E")
    p2 = _save(fig, out_dir / "Figure3_p2.png")
    return [p1, p2]


def compose_figure4(out_dir: Path = OUT) -> list[Path]:
    """Exploratory predictions on unseen ligands and heat-killed bacteria."""
    return _prediction_pages(
        "Figure4", RESULTS / "predictions" / primary_geneset_name(), out_dir
    )


def compose_figure4a(out_dir: Path = OUT) -> list[Path]:
    """Figure 4 composition applied to the external test batch (all models)."""
    val = (RESULTS / "validation" / "test_set" / primary_geneset_name()
           / "test_ligands")
    return [_model_grid_page(
        FIG / "pca" / "test_ligands_feature_selected.png",
        "PCA (external test batch)", val, list("ABCDEF"),
        out_dir / "Figure4a_external_test.png")]


def compose_supp_figure1(out_dir: Path = OUT) -> list[Path]:
    """DESeq2 volcano and heatmap for every stimulus except LPS (Figure 2)."""
    rows = [
        (ligand, subset)
        for subset in ("train_ligands", "additional_ligands", "bacterial_ligands")
        for ligand in CLASS_ORDER[subset]
        if ligand in DE_SUBSET_LIGANDS[subset]
        and ligand not in ("negative_control", "LPS")
    ]
    return _paginate_de(rows, "SupplementaryFigure1", out_dir)


def compose_supp_figure2(out_dir: Path = OUT) -> list[Path]:
    """HEK-Blue TLR2/TLR4 reporter dose-response."""
    fig = plt.figure(figsize=(8, 9))
    ax = fig.add_subplot(111)
    _place_image(ax, FIG / "supplementary" / "tlr_hek_blue.png")
    return [_save(fig, out_dir / "SupplementaryFigure2.png")]


def compose_supp_figure3(out_dir: Path = OUT) -> list[Path]:
    """Figure 4 composition for the models trained without Fla-PA."""
    return _prediction_pages(
        "SupplementaryFigure3",
        RESULTS / "predictions" / "no_flapa" / primary_geneset_name(),
        out_dir,
    )


def compose_supp_figure4_external_de(out_dir: Path = OUT) -> list[Path]:
    """DESeq2 volcano and heatmap for the external test batch."""
    rows = [
        (ligand, "test_ligands")
        for ligand in CLASS_ORDER["test_ligands"]
        if ligand in DE_SUBSET_LIGANDS["test_ligands"] and ligand != "negative_control"
    ]
    return _paginate_de(rows, "SupplementaryFigure4_external_DE", out_dir)


def compose_supp_figure5(out_dir: Path = OUT) -> list[Path]:
    """Gene-number selection: A) MI elbow curve, B) forest/k-means ARI sweep."""
    fs = FIG / "feature_selection"
    return [_row_page(
        [fs / "mutual_information.png", fs / "forest_ari_sweep.png"],
        [None, None], ["A", "B"], out_dir / "SupplementaryFigure5.png")]


def compose_figures(out_dir: Path = OUT) -> list[Path]:
    """Compose every manuscript figure page; returns the written PNG paths."""
    builders = [
        compose_figure2, compose_figure3, compose_figure4, compose_figure4a,
        compose_supp_figure1, compose_supp_figure2, compose_supp_figure3,
        compose_supp_figure4_external_de, compose_supp_figure5,
    ]
    paths = []
    for build in builders:
        paths.extend(build(out_dir))
    return paths


if __name__ == "__main__":
    compose_figures()
