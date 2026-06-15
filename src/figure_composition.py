"""Assemble the per-panel pipeline outputs into the manuscript composite figures.

Each emitted page holds at most four panels so it drops cleanly onto one Word page.
Panel white margins are auto-cropped to remove inter-panel whitespace. Original panel
letters are preserved across split pages so manuscript text references stay valid.

Run standalone (`python src/figure_composition.py`) or call `compose_all` from the
pipeline. Missing panels are flagged in place rather than breaking the run.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrow, FancyBboxPatch

RESULTS = Path(__file__).resolve().parent.parent / "results"
FIG = RESULTS / "figures"
PRED = RESULTS / "predictions" / "selected_100"
OUT = Path(__file__).resolve().parent.parent / "paper" / "paper_updated" / "figures"


# --------------------------------------------------------------------------- #
# panel helpers
# --------------------------------------------------------------------------- #
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


def _place_image(ax, path: Path, title: str | None = None) -> None:
    ax.axis("off")
    if not Path(path).is_file():
        ax.text(0.5, 0.5, f"missing:\n{Path(path).name}", ha="center", va="center",
                fontsize=7, color="0.6", transform=ax.transAxes)
        return
    ax.imshow(_load_cropped(path))
    if title:
        ax.set_title(title, fontsize=9)


def _letter(ax, letter: str | None) -> None:
    if letter:
        ax.text(-0.01, 1.01, letter, transform=ax.transAxes, fontsize=15,
                fontweight="bold", ha="right", va="bottom")


def _img(path, title=None):
    return lambda ax: _place_image(ax, path, title)


def _save(fig, out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=250, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return out_path


def _row_page(paths, titles, letters, out_path: Path, height: float = 4.6) -> Path:
    """One tight row; column widths track each cropped image so panels abut."""
    wr = [_aspect(p) for p in paths]
    fig = plt.figure(figsize=(height * sum(wr), height))
    gs = fig.add_gridspec(1, len(paths), width_ratios=wr, wspace=0.03)
    for i, (p, t, l) in enumerate(zip(paths, titles, letters)):
        ax = fig.add_subplot(gs[0, i])
        _place_image(ax, p, t)
        _letter(ax, l)
    return _save(fig, out_path)


# --------------------------------------------------------------------------- #
# drawn panels
# --------------------------------------------------------------------------- #
def _draw_pipeline_schematic(ax) -> None:
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    steps = [
        "Expression filter\n+ library-size\nnormalisation",
        "Mutual information\nranking",
        "MI elbow\n(~2,900 genes)",
        "ExtraTrees over\n1000 reseeded runs",
        "100 most\nstable genes",
        "Overlap with DE\n(58-gene set)",
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


def _read_gene_set(path: Path) -> set[str]:
    if not Path(path).is_file():
        return set()
    import csv

    with open(path) as fh:
        rows = list(csv.reader(fh))
    return {r[0] for r in rows[1:] if r}


def _draw_stable_de_venn(ax) -> None:
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    stable = _read_gene_set(RESULTS / "feature_selection" / "selected_genes_100.csv")
    de = _read_gene_set(
        RESULTS / "differential_gene_expression" / "de_genes_main_ligands.csv")
    inter = len(stable & de)
    only_s, only_de = len(stable) - inter, len(de) - inter
    ax.add_patch(Circle((0.40, 0.5), 0.30, facecolor="#e8a0a0", alpha=0.55,
                        edgecolor="#b05050", linewidth=1.2))
    ax.add_patch(Circle((0.62, 0.5), 0.34, facecolor="#9ec79e", alpha=0.55,
                        edgecolor="#508050", linewidth=1.2))
    ax.text(0.27, 0.5, str(only_s), ha="center", va="center", fontsize=12)
    ax.text(0.50, 0.5, str(inter), ha="center", va="center", fontsize=12,
            fontweight="bold")
    ax.text(0.78, 0.5, str(only_de), ha="center", va="center", fontsize=12)
    ax.text(0.30, 0.12, "Stable-100\nsignature", ha="center", va="center", fontsize=8.5)
    ax.text(0.74, 0.10, "DE genes\n(main set)", ha="center", va="center", fontsize=8.5)


# --------------------------------------------------------------------------- #
# main-text figures
# --------------------------------------------------------------------------- #
def compose_figure2(out_dir: Path = OUT) -> list[Path]:
    """LPS DESeq2: A) volcano, B) clustered expression heatmap."""
    de = FIG / "deseq2" / "main_ligands"
    return [_row_page(
        [de / "LPS_volcano.png", de / "LPS_histogram.png"],
        [None, None], ["A", "B"], out_dir / "Figure2.png", height=5.2)]


def compose_figure3(out_dir: Path = OUT) -> list[Path]:
    """Feature-selection figure, split into two <=4-panel pages (letters preserved)."""
    # page 1: A) schematic, B) PCA before / after
    fig = plt.figure(figsize=(11, 7))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.55, 1], hspace=0.1, wspace=0.05)
    ax_a = fig.add_subplot(gs[0, :])
    _draw_pipeline_schematic(ax_a)
    _letter(ax_a, "A")
    ax_b1 = fig.add_subplot(gs[1, 0])
    ax_b2 = fig.add_subplot(gs[1, 1])
    _place_image(ax_b1, FIG / "pca" / "main_ligands_pca.png", "before selection")
    _place_image(ax_b2, FIG / "pca" / "main_ligands_selected_pca.png", "after selection")
    _letter(ax_b1, "B")
    p1 = _save(fig, out_dir / "Figure3_p1.png")

    # page 2: C) venn, D) GO of selected genes, E) GO of selected non-DE genes
    fig = plt.figure(figsize=(11, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.12, wspace=0.05)
    ax_c = fig.add_subplot(gs[0, 0])
    ax_d = fig.add_subplot(gs[0, 1])
    ax_e = fig.add_subplot(gs[1, :])
    _draw_stable_de_venn(ax_c)
    _place_image(ax_d, FIG / "go" / "de_intersect_fs_go.png")
    _place_image(ax_e, FIG / "go" / "fs_only_go.png")
    _letter(ax_c, "C")
    _letter(ax_d, "D")
    _letter(ax_e, "E")
    p2 = _save(fig, out_dir / "Figure3_p2.png")
    return [p1, p2]


def _prediction_pages(base: str, pred_root: Path, out_dir: Path):
    """Two pages (unseen ligands / heat-killed bacteria); each a tight PCA+SVC+LogReg row."""
    specs = [
        ("unseen ligands", FIG / "pca" / "additional_ligands_selected_pca.png",
         pred_root / "additional_ligands", ["A", "B", "C"], f"{base}_p1.png"),
        ("heat-killed bacteria", FIG / "pca" / "bacteria_ligands_selected_pca.png",
         pred_root / "bacteria_ligands", ["D", "E", "F"], f"{base}_p2.png"),
    ]
    outs = []
    for label, pca_png, pdir, letters, fname in specs:
        outs.append(_row_page(
            [pca_png, pdir / "LinearSVC_probabilities_heatmap.png",
             pdir / "LogisticRegression_probabilities_heatmap.png"],
            [f"PCA ({label})", "LinearSVC", "Logistic Regression"],
            letters, out_dir / fname))
    return outs


def compose_figure4(out_dir: Path = OUT) -> list[Path]:
    """Exploratory predictions on unseen ligands and heat-killed bacteria."""
    return _prediction_pages("Figure4", PRED, out_dir)


def _external_test_name() -> str:
    val = RESULTS / "validation"
    names = [p.name for p in val.iterdir() if p.is_dir()] if val.is_dir() else []
    return names[0] if names else "7086"


def compose_figure4a(out_dir: Path = OUT) -> list[Path]:
    """Figure 4a: same composition as Figure 4 but on the external test batch."""
    val = RESULTS / "validation" / _external_test_name() / "selected_100" / "main_ligands"
    return [_row_page(
        [FIG / "pca" / "main_ligands_external_test_selected_pca.png",
         val / "LinearSVC_probabilities_heatmap.png",
         val / "LogisticRegression_probabilities_heatmap.png"],
        ["PCA (external test batch)", "LinearSVC", "Logistic Regression"],
        ["A", "B", "C"], out_dir / "Figure4a_external_test.png")]


# --------------------------------------------------------------------------- #
# supplementary figures
# --------------------------------------------------------------------------- #
def _paginate_de(rows, base: Path, out_dir: Path) -> list[Path]:
    """Two ligands (4 panels) per page; rows = (ligand, subset, vol_letter, heat_letter)."""
    chunks = [rows[i:i + 2] for i in range(0, len(rows), 2)]
    outs = []
    wr = None
    for ci, chunk in enumerate(chunks):
        n = len(chunk)
        fig = plt.figure(figsize=(12, 5.4 * n))
        if wr is None:
            de0 = FIG / "deseq2" / chunk[0][1]
            wr = [_aspect(de0 / f"{chunk[0][0]}_volcano.png"),
                  _aspect(de0 / f"{chunk[0][0]}_histogram.png")]
        gs = fig.add_gridspec(n, 2, width_ratios=wr, hspace=0.16, wspace=0.05)
        for r, (ligand, subset, vl, hl) in enumerate(chunk):
            de = FIG / "deseq2" / subset
            ax0 = fig.add_subplot(gs[r, 0])
            ax1 = fig.add_subplot(gs[r, 1])
            _place_image(ax0, de / f"{ligand}_volcano.png")
            _place_image(ax1, de / f"{ligand}_histogram.png")
            _letter(ax0, vl)
            _letter(ax1, hl)
        name = base if len(chunks) == 1 else base.with_name(
            f"{base.stem}_p{ci + 1}{base.suffix}")
        outs.append(_save(fig, out_dir / name.name))
    return outs


def compose_supp_figure1(out_dir: Path = OUT) -> list[Path]:
    rows = [
        ("Pam3", "main_ligands", "A", "B"),
        ("R848", "main_ligands", "C", "D"),
        ("PGN", "main_ligands", "E", "F"),
        ("Fla-PA", "main_ligands", "G", "I"),
        ("LTA", "additional_ligands", "J", "K"),
        ("MPLA", "additional_ligands", "L", "M"),
        ("Pam2", "additional_ligands", "N", "O"),
        ("HK E.coli", "bacteria_ligands", "P", "Q"),
        ("HK S.aureus", "bacteria_ligands", "R", "S"),
    ]
    return _paginate_de(rows, out_dir / "SupplementaryFigure1.png", out_dir)


def compose_supp_figure2(out_dir: Path = OUT) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(8, 9))
    ax = fig.add_subplot(111)
    _place_image(ax, FIG / "supplementary" / "tlr_hek_blue.png")
    return [_save(fig, out_dir / "SupplementaryFigure2.png")]


def compose_supp_figure3(out_dir: Path = OUT) -> list[Path]:
    pred = RESULTS / "predictions" / "no_flapa" / "selected_100"
    return _prediction_pages("SupplementaryFigure3", pred, out_dir)


def compose_supp_figure4_external_de(out_dir: Path = OUT) -> list[Path]:
    rows = [
        ("LPS", "main_external_test", "A", "B"),
        ("Fla-PA", "main_external_test", "C", "D"),
        ("Pam3", "main_external_test", "E", "F"),
        ("R848", "main_external_test", "G", "H"),
    ]
    return _paginate_de(rows, out_dir / "SupplementaryFigure4_external_DE.png", out_dir)


def compose_all(out_dir: Path = OUT) -> list[Path]:
    builders = [
        compose_figure2, compose_figure3, compose_figure4, compose_figure4a,
        compose_supp_figure1, compose_supp_figure2, compose_supp_figure3,
        compose_supp_figure4_external_de,
    ]
    paths = []
    for b in builders:
        paths.extend(b(out_dir))
    for p in paths:
        print(f"  Composite page saved: {p}")
    return paths


if __name__ == "__main__":
    compose_all()
