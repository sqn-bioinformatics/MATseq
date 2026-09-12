"""Compose the manuscript's multi-panel figure from rendered panel PNGs."""
from pathlib import Path
from string import ascii_uppercase

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np


def compose_figures(
    panel_paths: list[Path],
    output_path: Path,
    ncols: int = 2,
    max_size: tuple[float, float] = (6.69, 8.66),
) -> Path:
    """Lay panels out row by row, labelled A, B, ..., within max_size inches.

    The default max_size is an A4 text block (170 mm wide) with ~37 mm of
    height left for the caption. Panels share one column width, are aligned
    to the top-left of their cell and are cropped only of their white border.
    """
    missing = [p for p in panel_paths if not p.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing panel images: {missing}")
    images = []
    for path in panel_paths:
        img = mpimg.imread(path)
        ink = (img[..., :3] < 0.98).any(axis=-1)
        rows = np.flatnonzero(ink.any(axis=1))
        cols = np.flatnonzero(ink.any(axis=0))
        images.append(img[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1])
    gap, left, top, right, bottom = 0.3, 0.3, 0.25, 0.05, 0.05
    nrows = -(-len(images) // ncols)
    width = (max_size[0] - left - right - (ncols - 1) * gap) / ncols
    heights = [
        max(width * im.shape[0] / im.shape[1]
            for im in images[r * ncols:(r + 1) * ncols])
        for r in range(nrows)
    ]
    scale = min(
        1.0, (max_size[1] - top - bottom - (nrows - 1) * gap) / sum(heights)
    )
    width, heights = width * scale, [h * scale for h in heights]
    fig_w = left + ncols * width + (ncols - 1) * gap + right
    fig_h = top + sum(heights) + (nrows - 1) * gap + bottom
    fig = plt.figure(figsize=(fig_w, fig_h))
    for i, img in enumerate(images):
        r, c = divmod(i, ncols)
        h = width * img.shape[0] / img.shape[1]
        x0 = left + c * (width + gap)
        y_top = fig_h - top - sum(heights[:r]) - r * gap
        ax = fig.add_axes(
            (x0 / fig_w, (y_top - h) / fig_h, width / fig_w, h / fig_h)
        )
        ax.imshow(img)
        ax.axis("off")
        ax.annotate(
            ascii_uppercase[i], (0, 1), xycoords="axes fraction",
            xytext=(-6, 3), textcoords="offset points", ha="right",
            va="bottom", fontsize=12, fontweight="bold",
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    print(f"Figure saved: {output_path}")
    return output_path
