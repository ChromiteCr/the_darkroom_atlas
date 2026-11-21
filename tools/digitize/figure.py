"""Finding things in a chart: tick positions, tick labels, legend samples and curve paths,
and turning selected paths into ordered data points.

Works on :class:`~tools.digitize.page.Page` objects only, so it needs no PDF library.
"""

import itertools
import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from tools.digitize.calibrate import Axis
from tools.digitize.page import BBox, Page, VectorPath, Word, bbox_contains, dash_equal

__all__ = [
    "AssembleStats",
    "DigitizeError",
    "assemble",
    "find_phrase",
    "has_bar",
    "is_axis_aligned",
    "label_words",
    "legend_sample",
    "match_style",
    "nearest_to_annotation",
    "tick_candidates",
    "to_data",
]

_AXIS_ALIGNED_TOL = 0.05  # pt


class DigitizeError(Exception):
    """A recipe cannot be applied to the page as written."""


# --------------------------------------------------------------------------------------
# Ticks and labels
# --------------------------------------------------------------------------------------


def _edge(frame: BBox, side: str) -> float:
    return {"bottom": frame[1], "top": frame[3], "left": frame[0], "right": frame[2]}[side]


def tick_candidates(
    page: Page, frame: BBox, side: str, tolerance: float = 0.6
) -> list[tuple[float, str]]:
    """Positions of tick marks, grid lines and frame edges that meet one side of ``frame``.

    For ``bottom``/``top`` these are vertical straight segments touching that edge, returned
    as x coordinates; for ``left``/``right`` horizontal segments, returned as y. Each comes
    with a short description of the path it belongs to. Sorted, near-duplicates (within
    0.05 pt) merged.
    """
    edge = _edge(frame, side)
    horizontal_axis = side in ("bottom", "top")
    lo, hi = (frame[0], frame[2]) if horizontal_axis else (frame[1], frame[3])
    found: list[tuple[float, str]] = []
    for path in page.paths:
        if not path.stroke:
            continue
        for p, q in path.segments():
            if horizontal_axis:
                along, across = (p[0], q[0]), (p[1], q[1])
            else:
                along, across = (p[1], q[1]), (p[0], q[0])
            if abs(along[0] - along[1]) > _AXIS_ALIGNED_TOL:
                continue
            if abs(across[0] - across[1]) < 0.3:
                continue  # a dot, not a tick
            if not (min(across) - tolerance <= edge <= max(across) + tolerance):
                continue
            coord = (along[0] + along[1]) / 2
            if lo - tolerance <= coord <= hi + tolerance:
                kind = "frame edge" if path.is_rect else "tick/grid line"
                found.append((coord, f"{kind}, path #{path.index}"))
    found.sort()
    merged: list[tuple[float, str]] = []
    for coord, source in found:
        if merged and abs(coord - merged[-1][0]) <= 0.05:
            continue
        merged.append((coord, source))
    return merged


def label_words(page: Page, frame: BBox, side: str, band: float, margin: float) -> list[Word]:
    """Words in the band just outside one side of ``frame`` (where tick labels sit).

    The band is ``band`` points deep; along the side it extends ``margin`` points past the
    corners. Words whose centre is within 2 pt of the edge line are left out, which keeps
    the labels of the other axis out of the band.
    """
    x0, y0, x1, y1 = frame
    out = []
    for word in page.words:
        cx, cy = word.centre
        if side == "bottom":
            ok = y0 - band <= cy <= y0 - 2 and x0 - margin <= cx <= x1 + margin
        elif side == "top":
            ok = y1 + 2 <= cy <= y1 + band and x0 - margin <= cx <= x1 + margin
        elif side == "left":
            ok = x0 - band <= cx <= x0 - 2 and y0 - margin <= cy <= y1 + margin
        else:
            ok = x1 + 2 <= cx <= x1 + band and y0 - margin <= cy <= y1 + margin
        if ok:
            out.append(word)
    return out


def has_bar(word: Word, page: Page) -> bool:
    """True if a short horizontal stroke sits over the first glyph of ``word``.

    Kodak prints negative logarithms in bar notation: a bar over the characteristic, so
    that a label reading 3.0 with a bar over the 3 means -3.0. The bar is a path, not a
    glyph, so the text layer alone reads "3.0".
    """
    glyph = word.chars[0]
    width = glyph.x1 - glyph.x0
    height = glyph.y1 - glyph.y0
    for path in page.paths:
        if not path.stroke:
            continue
        # The bar may be its own path or one segment of a path that also draws ticks.
        for (ax, ay), (bx, by) in path.segments():
            if abs(ay - by) > _AXIS_ALIGNED_TOL:
                continue
            left, right = min(ax, bx), max(ax, bx)
            length = right - left
            if not 0.2 * width <= length <= 1.6 * width + 0.5:
                continue
            overlap = min(right, glyph.x1) - max(left, glyph.x0)
            if overlap < 0.5 * length:
                continue
            if glyph.y0 + 0.5 * height <= ay <= glyph.y1 + 0.5 * glyph.size:
                return True
    return False


def find_phrase(page: Page, phrase: str, region: BBox) -> BBox:
    """Box of the text run inside ``region`` that reads exactly ``phrase``.

    A text run is a stretch of words on one line with gaps of less than one em between
    them, so the legend entry "D-76" does not match inside "D-76 (1:1)". Words are compared
    with single spaces between them. Exactly one match must exist.
    """
    target = phrase.split()
    words = page.words_in(region)
    lines: list[list[Word]] = []
    for word in sorted(words, key=lambda w: (-w.centre[1], w.bbox[0])):
        for line in lines:
            if abs(line[0].centre[1] - word.centre[1]) <= 0.4 * max(line[0].size, word.size):
                line.append(word)
                break
        else:
            lines.append([word])
    runs: list[list[Word]] = []
    for line in lines:
        line.sort(key=lambda w: w.bbox[0])
        run = [line[0]]
        for word in line[1:]:
            if word.bbox[0] - run[-1].bbox[2] < max(word.size, run[-1].size):
                run.append(word)
            else:
                runs.append(run)
                run = [word]
        runs.append(run)
    matches = []
    for run in runs:
        if [w.text for w in run] == target:
            matches.append(
                (
                    min(w.bbox[0] for w in run),
                    min(w.bbox[1] for w in run),
                    max(w.bbox[2] for w in run),
                    max(w.bbox[3] for w in run),
                )
            )
    if len(matches) != 1:
        raise DigitizeError(
            f"legend text {phrase!r}: expected exactly one match inside the frame, "
            f"found {len(matches)}"
        )
    return matches[0]


def legend_sample(
    page: Page, label: BBox, size_hint: float = 8.0, max_gap: float = 12.0
) -> VectorPath:
    """The legend sample line drawn just left of a legend label.

    A sample is a stroked, straight, horizontal line 5-80 pt long whose right end lies at
    most ``max_gap`` pt (default 12) left of the label and whose height is within 0.6 em of
    the label's centre. The one nearest the label wins; two equally near samples are an
    error.
    """
    cy = (label[1] + label[3]) / 2
    em = max(label[3] - label[1], size_hint * 0.5)
    candidates = []
    for path in page.paths:
        line = path.straight_line()
        if line is None or not path.stroke:
            continue
        (ax, ay), (bx, by) = line
        if abs(ay - by) > _AXIS_ALIGNED_TOL:
            continue
        right = max(ax, bx)
        length = abs(bx - ax)
        if not 5.0 <= length <= 80.0:
            continue
        if not (label[0] - max_gap <= right <= label[0] + 0.5):
            continue
        if abs(ay - cy) > 0.6 * em:
            continue
        candidates.append((label[0] - right + abs(ay - cy), path))
    if not candidates:
        raise DigitizeError(f"no legend sample line found left of the label at {label}")
    candidates.sort(key=lambda item: item[0])
    if len(candidates) > 1 and math.isclose(candidates[0][0], candidates[1][0], abs_tol=0.05):
        raise DigitizeError(f"two legend sample lines are equally near the label at {label}")
    return candidates[0][1]


# --------------------------------------------------------------------------------------
# Curve paths
# --------------------------------------------------------------------------------------


def is_axis_aligned(path: VectorPath) -> bool:
    """True if every piece of ``path`` is a horizontal or vertical straight segment.

    Such paths are frames, ticks, grid lines and legend samples; data curves never are.
    """
    if any(op.kind in ("c", "v", "y") for sub in path.subpaths for op in sub):
        return False
    segments = path.segments()
    return all(
        abs(p[0] - q[0]) <= _AXIS_ALIGNED_TOL or abs(p[1] - q[1]) <= _AXIS_ALIGNED_TOL
        for p, q in segments
    )


def match_style(
    path: VectorPath,
    width: float,
    dash: Sequence[float],
    stroke_color: Sequence[float] | None = None,
    tolerance: float = 0.01,
) -> bool:
    """True if ``path`` is stroked with this width, dash array and (if given) colour."""
    if not path.stroke or not math.isclose(path.line_width, width, abs_tol=tolerance):
        return False
    if not dash_equal(path.dash, dash, tolerance):
        return False
    if stroke_color is not None:
        color = path.stroke_color
        if isinstance(color, str) or color is None or len(color) != len(stroke_color):
            return False
        if any(abs(a - b) > 0.01 for a, b in zip(color, stroke_color)):
            return False
    return True


def _vertical_gap(polylines: Sequence[np.ndarray], box: BBox) -> float:
    """Smallest vertical distance between the polylines and ``box``, over the box's width.

    Only the parts of the polylines whose x lies within the box's horizontal extent count;
    ``inf`` if no part does. 0 if a polyline passes through the box.
    """
    bx0, by0, bx1, by1 = box
    best = math.inf
    for poly in polylines:
        for (x0, y0), (x1, y1) in itertools.pairwise(poly):
            lo, hi = min(x0, x1), max(x0, x1)
            if hi < bx0 or lo > bx1:
                continue
            if x1 == x0:
                ends = [y0, y1]
            else:
                ends = []
                for x in (max(lo, bx0), min(hi, bx1)):
                    ends.append(y0 + (y1 - y0) * (x - x0) / (x1 - x0))
            low, high = min(ends), max(ends)
            if high >= by0 and low <= by1:
                return 0.0
            best = min(best, by0 - high if high < by0 else low - by1)
    return best


def nearest_to_annotation(
    paths: Sequence[VectorPath], box: BBox, tolerance: float, ratio: float = 2.0
) -> VectorPath:
    """The path that passes vertically nearest to an annotation printed beside a curve.

    Some figures label their curves with text placed next to each curve instead of a
    legend (e.g. "0.3 greater than D-min"). Over the horizontal extent of the text's box,
    the vertical gap between each path and the box is measured; the nearest path wins. The
    runner-up must be at least ``ratio`` times as far (plus 1 pt), or the choice is
    ambiguous and an error.
    """
    gaps = sorted(
        ((_vertical_gap(p.polylines(tolerance), box), p.index, p) for p in paths),
        key=lambda item: (item[0], item[1]),
    )
    if not gaps or math.isinf(gaps[0][0]):
        raise DigitizeError(f"no candidate path passes above or below the annotation at {box}")
    if len(gaps) > 1 and gaps[1][0] < ratio * gaps[0][0] + 1.0:
        raise DigitizeError(
            f"the annotation at {box} is about as near path #{gaps[0][1]} "
            f"({gaps[0][0]:.2f} pt) as path #{gaps[1][1]} ({gaps[1][0]:.2f} pt)"
        )
    return gaps[0][2]


def curve_candidates(
    page: Page, frame: BBox, exclude: Sequence[BBox], skip: set[int]
) -> list[VectorPath]:
    """Stroked paths inside ``frame`` that could be data curves.

    Leaves out rectangles, axis-aligned paths (ticks, grid lines, legend samples), paths
    shorter than 2 pt, paths listed in ``skip`` and paths inside an ``exclude`` box.
    """
    out = []
    for path in page.paths:
        if not path.stroke or path.index in skip or path.is_rect or is_axis_aligned(path):
            continue
        if not bbox_contains(frame, path.bbox, 1.0):
            continue
        if any(bbox_contains(box, path.bbox, 0.5) for box in exclude):
            continue
        x0, y0, x1, y1 = path.bbox
        if math.hypot(x1 - x0, y1 - y0) < 2.0:
            continue
        out.append(path)
    return out


# --------------------------------------------------------------------------------------
# Paths to data
# --------------------------------------------------------------------------------------


def to_data(polylines: Sequence[np.ndarray], x_axis: Axis, y_axis: Axis) -> list[np.ndarray]:
    """Convert page-coordinate polylines to data coordinates."""
    out = []
    for poly in polylines:
        xs = x_axis.to_data(poly[:, 0])
        ys = y_axis.to_data(poly[:, 1])
        out.append(np.column_stack([xs, ys]))
    return out


@dataclass(frozen=True)
class AssembleStats:
    """What :func:`assemble` did."""

    raw_points: int
    duplicates: int  # identical to the previous point after rounding
    same_x: int  # same rounded x as the previous point but a different y
    backtracking: int  # x smaller than the previous kept point
    reversed_pieces: int  # polylines drawn right to left, reversed


def assemble(
    polylines: Sequence[np.ndarray], decimals: tuple[int, int], policy: str
) -> tuple[np.ndarray, np.ndarray, AssembleStats]:
    """Join data-coordinate polylines into one curve with strictly increasing x.

    Each polyline is turned to run left to right, the pieces are put in order of their
    first x and joined, and values are rounded to ``decimals = (x, y)`` places. Then each
    point is kept only if its x exceeds the last kept x. Dropped points are counted by
    kind: exact duplicates and same-x points always drop; points whose x goes backwards
    drop only under ``policy="drop-backtracking"`` and raise :class:`DigitizeError` under
    ``policy="strict"``. No point is moved, averaged or interpolated.
    """
    if policy not in ("strict", "drop-backtracking"):
        raise ValueError(f"unknown monotonic policy {policy!r}")
    pieces = []
    reversed_pieces = 0
    for poly in polylines:
        if len(poly) == 0:
            continue
        if poly[-1, 0] < poly[0, 0]:
            poly = poly[::-1]
            reversed_pieces += 1
        pieces.append(poly)
    pieces.sort(key=lambda p: (float(p[0, 0]), float(p[0, 1])))
    if not pieces:
        raise DigitizeError("no points to assemble")
    points = np.concatenate(pieces)
    xs = np.round(points[:, 0], decimals[0])
    ys = np.round(points[:, 1], decimals[1])
    keep_x: list[float] = []
    keep_y: list[float] = []
    duplicates = same_x = backtracking = 0
    for x, y in zip(xs, ys):
        x, y = float(x), float(y)
        if not keep_x or x > keep_x[-1]:
            keep_x.append(x)
            keep_y.append(y)
        elif x == keep_x[-1] and y == keep_y[-1]:
            duplicates += 1
        elif x == keep_x[-1]:
            same_x += 1
        else:
            if policy == "strict":
                raise DigitizeError(
                    f"the curve goes backwards in x at ({x}, {y}) after ({keep_x[-1]}, "
                    f"{keep_y[-1]}); set monotonic: drop-backtracking in the recipe if "
                    "dropping such points is acceptable, and say why"
                )
            backtracking += 1
    if len(keep_x) < 2:
        raise DigitizeError("fewer than two points remain after ordering")
    stats = AssembleStats(len(points), duplicates, same_x, backtracking, reversed_pieces)
    return np.asarray(keep_x), np.asarray(keep_y), stats
