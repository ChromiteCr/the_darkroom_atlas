"""Path geometry: path operators, Bezier flattening and straight segments.

Pure numpy; no PDF library is needed. Coordinates are PDF page coordinates in points
(1/72 inch), origin at the bottom-left corner of the page, y pointing up.
"""

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

__all__ = [
    "Op",
    "Point",
    "flatten_cubic",
    "flatten_subpath",
    "polyline_deviation",
    "straight_segments",
]

Point = tuple[float, float]

_MAX_DEPTH = 24


@dataclass(frozen=True)
class Op:
    """One path construction operator with its points already in page coordinates.

    ``kind`` is ``m`` (move to), ``l`` (line to), ``c`` (cubic Bezier: two control points
    and the end point), ``v`` (cubic whose first control point is the current point: one
    control point and the end point), ``y`` (cubic whose second control point is the end
    point: one control point and the end point) or ``h`` (close the subpath, no points).
    """

    kind: str
    points: tuple[Point, ...]

    def __post_init__(self) -> None:
        expected = {"m": 1, "l": 1, "c": 3, "v": 2, "y": 2, "h": 0}
        if self.kind not in expected:
            raise ValueError(f"unknown path operator {self.kind!r}")
        if len(self.points) != expected[self.kind]:
            raise ValueError(
                f"operator {self.kind!r} takes {expected[self.kind]} point(s), "
                f"got {len(self.points)}"
            )


def _distance_to_chord(p: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    chord = b - a
    length = float(np.hypot(*chord))
    if length == 0.0:
        return float(np.hypot(*(p - a)))
    # |cross(chord, p - a)| / |chord|
    return abs(float(chord[0] * (p[1] - a[1]) - chord[1] * (p[0] - a[0]))) / length


def flatten_cubic(p0: Point, p1: Point, p2: Point, p3: Point, tolerance: float) -> list[Point]:
    """Points along the cubic Bezier ``p0 p1 p2 p3``, excluding ``p0``, ending with ``p3``.

    The curve is split at its midpoint (de Casteljau) until both control points of every
    piece lie within ``tolerance`` of that piece's chord. A cubic lies inside the convex
    hull of its control points, so the returned polyline deviates from the curve by at most
    ``tolerance`` (in the same units as the points). The subdivision is deterministic.
    """
    if not tolerance > 0.0:
        raise ValueError(f"tolerance must be positive, got {tolerance!r}")
    out: list[Point] = []

    def recurse(a, b, c, d, depth):
        flat = max(_distance_to_chord(b, a, d), _distance_to_chord(c, a, d)) <= tolerance
        if flat or depth >= _MAX_DEPTH:
            out.append((float(d[0]), float(d[1])))
            return
        ab, bc, cd = (a + b) / 2, (b + c) / 2, (c + d) / 2
        abc, bcd = (ab + bc) / 2, (bc + cd) / 2
        mid = (abc + bcd) / 2
        recurse(a, ab, abc, mid, depth + 1)
        recurse(mid, bcd, cd, d, depth + 1)

    recurse(*(np.asarray(p, dtype=np.float64) for p in (p0, p1, p2, p3)), 0)
    return out


def flatten_subpath(ops: Sequence[Op], tolerance: float) -> np.ndarray:
    """Polyline, shape ``(N, 2)``, through one subpath; Beziers flattened to ``tolerance``.

    The subpath must start with ``m``. A closing ``h`` adds the start point unless the
    subpath already ends there.
    """
    if not ops or ops[0].kind != "m":
        raise ValueError("a subpath must start with an 'm' operator")
    start = current = ops[0].points[0]
    points: list[Point] = [start]
    for op in ops[1:]:
        if op.kind == "m":
            raise ValueError("a subpath contains only one 'm' operator; split the path first")
        if op.kind == "l":
            current = op.points[0]
            points.append(current)
        elif op.kind == "c":
            points.extend(flatten_cubic(current, *op.points, tolerance))
            current = op.points[-1]
        elif op.kind == "v":
            points.extend(flatten_cubic(current, current, *op.points, tolerance))
            current = op.points[-1]
        elif op.kind == "y":
            control, end = op.points
            points.extend(flatten_cubic(current, control, end, end, tolerance))
            current = end
        elif op.kind == "h":
            if current != start:
                points.append(start)
            current = start
    return np.asarray(points, dtype=np.float64)


def straight_segments(ops: Sequence[Op]) -> list[tuple[Point, Point]]:
    """The straight pieces (``l`` operators and ``h`` closings) of one subpath."""
    segments: list[tuple[Point, Point]] = []
    start = current = None
    for op in ops:
        if op.kind == "m":
            start = current = op.points[0]
        elif op.kind == "l":
            segments.append((current, op.points[0]))
            current = op.points[0]
        elif op.kind in ("c", "v", "y"):
            current = op.points[-1]
        elif op.kind == "h":
            if current is not None and current != start:
                segments.append((current, start))
            current = start
    return segments


def polyline_deviation(polyline: np.ndarray, samples: np.ndarray) -> float:
    """Largest distance from any of ``samples`` (shape ``(M, 2)``) to ``polyline``."""
    a = polyline[:-1][None, :, :]
    b = polyline[1:][None, :, :]
    p = samples[:, None, :]
    ab = b - a
    denom = np.sum(ab * ab, axis=-1)
    t = np.where(denom > 0, np.sum((p - a) * ab, axis=-1) / np.where(denom > 0, denom, 1), 0)
    t = np.clip(t, 0.0, 1.0)
    nearest = a + t[..., None] * ab
    return float(np.max(np.min(np.hypot(*(p - nearest).transpose(2, 0, 1)), axis=1)))
