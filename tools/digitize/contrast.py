"""Contrast index read off a digitised characteristic curve, for checking CI-versus-time curves.

Run: ``uv run python -I -m tools.digitize.contrast CHAR.csv CI.csv MINUTES [CHAR.csv CI.csv
MINUTES ...]`` with reference CSV paths. For each triple it prints the contrast index read off the
characteristic curve, the value the CI curve gives at that development time (linear interpolation
between its points, nothing outside its span), and the difference.

The construction is Kodak's contrast index as a collinear chord: pick a centre on the fog line
(the density at the curve's left end, where the digitised curves are flat), mark where the curve
is 0.2 and 2.2 log units (distance in the log-exposure/density plane) from the centre, and slide
the centre along the fog line until the centre and both marks lie on one straight line. The
slope of that line is the contrast index. A simpler chord (from 0.1 above fog to the density 2.0
log exposure further on) is printed alongside; the two agree to a few hundredths.

This reads only the reference CSVs (no PDF) and needs numpy only. It changes nothing; the numbers it
prints are quoted in the ``ci_cross_check`` conditions of the F-4016 recipes.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

from tools.digitize.runner import read_staged

__all__ = ["chord_index", "contrast_index", "main"]

INNER = 0.2
OUTER = 2.2


def _first_hit(x: np.ndarray, y: np.ndarray, cx: float, cy: float, r: float):
    """First point, walking left to right along the polyline, at distance ``r`` from (cx, cy)
    where the distance passes from below ``r`` to at least ``r``; None if there is none."""
    for i in range(len(x) - 1):
        ax, ay = x[i] - cx, y[i] - cy
        dx, dy = x[i + 1] - x[i], y[i + 1] - y[i]
        da, db = math.hypot(ax, ay), math.hypot(ax + dx, ay + dy)
        if not (da < r <= db):
            continue
        # |a + t d| = r, solved for the root in [0, 1].
        qa = dx * dx + dy * dy
        qb = 2.0 * (ax * dx + ay * dy)
        qc = ax * ax + ay * ay - r * r
        root = math.sqrt(max(qb * qb - 4.0 * qa * qc, 0.0))
        t = (-qb + root) / (2.0 * qa)
        return x[i] + t * dx, y[i] + t * dy
    return None


def _misalignment(x, y, cx: float, fog: float):
    p1 = _first_hit(x, y, cx, fog, INNER)
    p2 = _first_hit(x, y, cx, fog, OUTER)
    if p1 is None or p2 is None:
        return None
    a1 = math.atan2(p1[1] - fog, p1[0] - cx)
    a2 = math.atan2(p2[1] - fog, p2[0] - cx)
    return a1 - a2, p1, p2


def contrast_index(x, y, fog: float | None = None, steps: int = 2000) -> float:
    """Collinear 0.2/2.2 contrast index of the curve (x log exposure, y density)."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    fog = float(y[0]) if fog is None else float(fog)
    centres = np.linspace(x[0] - OUTER, x[-1], steps)
    previous = None
    for c in centres:
        m = _misalignment(x, y, float(c), fog)
        if m is None:
            previous = None
            continue
        if previous is not None and previous[1] * m[0] <= 0.0:
            lo, hi, flo = previous[0], float(c), previous[1]
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                mm = _misalignment(x, y, mid, fog)
                if mm is None:
                    break
                if flo * mm[0] <= 0.0:
                    hi = mid
                else:
                    lo, flo = mid, mm[0]
            _, p1, p2 = _misalignment(x, y, 0.5 * (lo + hi), fog)
            return (p2[1] - p1[1]) / (p2[0] - p1[0])
        previous = (float(c), m[0])
    raise ValueError("no centre on the fog line makes the 0.2 and 2.2 marks collinear")


def chord_index(x, y, fog: float | None = None) -> float:
    """Slope from 0.1 above fog to the point 2.0 log exposure further on."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    fog = float(y[0]) if fog is None else float(fog)
    i = int(np.argmax(y >= fog + 0.1))
    if i == 0:
        raise ValueError("the curve never rises 0.1 above fog")
    xa = float(np.interp(fog + 0.1, y[i - 1 : i + 1], x[i - 1 : i + 1]))
    if xa + 2.0 > x[-1]:
        raise ValueError("the curve ends less than 2.0 log exposure above the 0.1 point")
    return (float(np.interp(xa + 2.0, x, y)) - (fog + 0.1)) / 2.0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -I -m tools.digitize.contrast")
    ap.add_argument("triples", nargs="+", help="CHAR.csv CI.csv MINUTES, repeated")
    args = ap.parse_args(argv)
    if len(args.triples) % 3:
        ap.error("give CHAR.csv CI.csv MINUTES triples")
    for k in range(0, len(args.triples), 3):
        char, ci, minutes = Path(args.triples[k]), Path(args.triples[k + 1]), args.triples[k + 2]
        x, y = read_staged(char)
        tx, tci = read_staged(ci)
        t = float(minutes)
        collinear = contrast_index(x, y)
        chord = chord_index(x, y)
        if tx[0] <= t <= tx[-1]:
            on_curve = float(np.interp(t, tx, tci))
            tail = f"{ci.name} at {t:g} min {on_curve:.3f}, difference {collinear - on_curve:+.3f}"
        else:
            tail = f"{ci.name} spans {tx[0]:g}-{tx[-1]:g} min, {t:g} min is outside: not compared"
        print(
            f"{char.name}: fog {y[0]:.4f}, CI collinear {collinear:.3f}, chord {chord:.3f}; {tail}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
