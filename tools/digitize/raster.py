"""Fallback for figures that are embedded images: trace a curve in a rendered raster.

The page (or the image inside it) is rendered to greyscale at a stated resolution. Pixels
darker than a threshold are ink. Starting from a seed point on the curve, the tracer walks
column by column to the right and then to the left; in each column it takes the run of ink
whose centre is nearest to where the curve was heading, and it bridges gaps (dashes) of up
to a stated width. It returns run centres and run heights, so the caller can state the
error: half the run height plus half a pixel, which is much larger than for vector paths.

Pure numpy; :func:`render_gray` imports pypdfium2 only when called.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = ["Trace", "ink_runs", "render_gray", "trace_curve"]


@dataclass(frozen=True)
class Trace:
    """A traced curve in pixel coordinates (column, row; row 0 at the top)."""

    cols: np.ndarray
    rows: np.ndarray  # run centres
    heights: np.ndarray  # run heights, in pixels


def ink_runs(column: np.ndarray) -> list[tuple[int, int]]:
    """``(start, stop)`` row ranges of consecutive ``True`` values in a boolean column."""
    padded = np.concatenate([[False], column.astype(bool), [False]])
    diff = np.diff(padded.astype(np.int8))
    starts = np.flatnonzero(diff == 1)
    stops = np.flatnonzero(diff == -1)
    return list(zip(starts.tolist(), stops.tolist()))


def _walk(ink, col, row, slope, step, max_jump, max_gap, max_height):
    cols, rows, heights = [], [], []
    width = ink.shape[1]
    gap = 0
    c = col
    while True:
        c += step
        if not 0 <= c < width:
            break
        predicted = row + slope * (gap + 1) * step
        runs = [(start, stop) for start, stop in ink_runs(ink[:, c]) if stop - start <= max_height]
        best = None
        for start, stop in runs:
            centre = (start + stop - 1) / 2
            distance = abs(centre - predicted)
            if distance <= max_jump * (gap + 1) and (best is None or distance < best[0]):
                best = (distance, centre, stop - start)
        if best is None:
            gap += 1
            if gap > max_gap:
                break
            continue
        _, centre, height = best
        slope = 0.7 * slope + 0.3 * (centre - row) / ((gap + 1) * step)
        row = centre
        gap = 0
        cols.append(c)
        rows.append(centre)
        heights.append(height)
    return cols, rows, heights


def trace_curve(
    ink: np.ndarray,
    seed: tuple[int, int],
    *,
    max_jump: float = 3.0,
    max_gap: int = 12,
    max_height: int = 40,
) -> Trace:
    """Follow the curve through ``seed = (col, row)`` across the boolean ``ink`` image.

    ``max_jump`` is how far (pixels) a run centre may be from the predicted row per column
    stepped; ``max_gap`` how many empty columns (dash gaps) may be bridged; runs taller than
    ``max_height`` (grid lines crossed lengthwise, text blocks) are ignored.
    """
    col, row = seed
    runs = ink_runs(ink[:, col])
    if not runs:
        raise ValueError(f"no ink in column {col} at the seed")
    start, stop = min(runs, key=lambda r: abs((r[0] + r[1] - 1) / 2 - row))
    centre = (start + stop - 1) / 2
    if abs(centre - row) > max_jump * 2:
        raise ValueError(f"the seed ({col}, {row}) is not on a curve (nearest ink row {centre})")
    right = _walk(ink, col, centre, 0.0, 1, max_jump, max_gap, max_height)
    left = _walk(ink, col, centre, 0.0, -1, max_jump, max_gap, max_height)
    cols = left[0][::-1] + [col] + right[0]
    rows = left[1][::-1] + [centre] + right[1]
    heights = left[2][::-1] + [stop - start] + right[2]
    return Trace(np.asarray(cols, float), np.asarray(rows, float), np.asarray(heights, float))


def mask_boxes(ink: np.ndarray, boxes_px: Sequence[tuple[int, int, int, int]]) -> np.ndarray:
    """Copy of ``ink`` with the given ``(col0, row0, col1, row1)`` boxes cleared."""
    out = ink.copy()
    for c0, r0, c1, r1 in boxes_px:
        out[max(r0, 0) : max(r1, 0), max(c0, 0) : max(c1, 0)] = False
    return out


def render_gray(pdf_path: Path | str, page_number: int, scale: float) -> np.ndarray:
    """Greyscale rendering (uint8, row 0 at the top) of a page at ``scale`` pixels per pt."""
    try:
        import pypdfium2 as pdfium
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise SystemExit("pypdfium2 is missing: run with `uv run --group digitize ...`") from exc
    document = pdfium.PdfDocument(str(pdf_path))
    try:
        image = document[page_number - 1].render(scale=scale, grayscale=True).to_pil()
        return np.asarray(image.convert("L"), dtype=np.uint8)
    finally:
        document.close()
