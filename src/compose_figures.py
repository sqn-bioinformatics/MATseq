"""Compose the manuscript's multi-panel figure from rendered panel PNGs."""
import textwrap
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
    title: str | None = None,
    first_letter: int = 0,
    margins: tuple[float, float, float, float] = (0.3, 0.3, 0.05, 0.05),
    spans: list[int] | None = None,
    width_ratios: tuple[float, ...] | None = None,
) -> Path:
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
    gap, left, right, bottom = margins
    heading = textwrap.fill(title, 90) if title else ""
    top = 0.25 + (0.22 * (heading.count("\n") + 1) if title else 0.0)
    spans = spans or [1] * len(images)
    ratios = width_ratios or (1.0,) * ncols
    avail = max_size[0] - left - right - (ncols - 1) * gap
    widths = [avail * r / sum(ratios) for r in ratios]
    next_free = [0] * ncols
    cells = []
    for span in spans:
        c = next_free.index(min(next_free))
        cells.append((next_free[c], c, span))
        next_free[c] += span
    nrows = max(next_free)
    heights = [0.0] * nrows
    for im, (r, c, span) in zip(images, cells):
        if span == 1:
            heights[r] = max(heights[r], widths[c] * im.shape[0] / im.shape[1])
    for im, (r, c, span) in zip(images, cells):
        if span > 1:
            deficit = (
                widths[c] * im.shape[0] / im.shape[1]
                - sum(heights[r:r + span]) - (span - 1) * gap
            )
            heights[r + span - 1] += max(deficit, 0.0)
    scale = min(
        1.0, (max_size[1] - top - bottom - (nrows - 1) * gap) / sum(heights)
    )
    widths, heights = [w * scale for w in widths], [h * scale for h in heights]
    fig_w = left + sum(widths) + (ncols - 1) * gap + right
    fig_h = top + sum(heights) + (nrows - 1) * gap + bottom
    fig = plt.figure(figsize=(fig_w, fig_h))
    if title:
        fig.text(0.5, 1 - 0.1 / fig_h, heading, ha="center", va="top", fontsize=13)
    for i, (img, (r, c, span)) in enumerate(zip(images, cells)):
        box_h = sum(heights[r:r + span]) + (span - 1) * gap
        w = widths[c]
        h = w * img.shape[0] / img.shape[1]
        if h > box_h:
            h, w = box_h, box_h * img.shape[1] / img.shape[0]
        x0 = left + sum(widths[:c]) + c * gap
        y_top = fig_h - top - sum(heights[:r]) - r * gap
        ax = fig.add_axes(
            (x0 / fig_w, (y_top - h) / fig_h, w / fig_w, h / fig_h)
        )
        ax.imshow(img)
        ax.axis("off")
        ax.annotate(
            ascii_uppercase[first_letter + i], (0, 1), xycoords="axes fraction",
            xytext=(-8, 6), textcoords="offset points", ha="right",
            va="bottom", fontsize=12, fontweight="bold",
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    print(f"Figure saved: {output_path}")
    return output_path
