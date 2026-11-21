"""Axis calibration: map a page coordinate (points) to a data value and back.

An axis is calibrated from two or more reference points ``(coordinate, value)``, normally
labelled tick marks. On a ``linear`` axis the value itself is linear in the coordinate; on a
``log10`` axis its base-10 logarithm is. The "axis variable" ``t`` is the value (linear) or
``log10(value)`` (log10).

Two fits:

* ``lsq``: one straight line ``t = a + b * coordinate``, least squares over all reference
  points. The residuals measure how far the printed ticks are from a perfect axis.
* ``piecewise``: straight lines between neighbouring reference points (extended past the
  ends with the end segments), so the transform passes exactly through every reference
  point. Use it when the printed grid itself is uneven. The least-squares residual is still
  reported, because it measures that unevenness.

Pure numpy.
"""

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

__all__ = ["Axis", "Reference"]

Scale = Literal["linear", "log10"]
Fit = Literal["lsq", "piecewise"]


@dataclass(frozen=True)
class Reference:
    """One calibration point: a page coordinate (pt) and the data value printed there."""

    coord: float
    value: float
    source: str = ""


@dataclass(frozen=True)
class Axis:
    """A calibrated axis. Raises ``ValueError`` for unusable reference points."""

    references: tuple[Reference, ...]
    scale: Scale = "linear"
    fit: Fit = "lsq"
    _a: float = field(init=False, repr=False)
    _b: float = field(init=False, repr=False)
    _coords: np.ndarray = field(init=False, repr=False)
    _ts: np.ndarray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        refs = tuple(self.references)
        object.__setattr__(self, "references", refs)
        if self.scale not in ("linear", "log10"):
            raise ValueError(f"unknown axis scale {self.scale!r}")
        if self.fit not in ("lsq", "piecewise"):
            raise ValueError(f"unknown axis fit {self.fit!r}")
        if len(refs) < 2:
            raise ValueError("an axis needs at least two reference points")
        coords = np.array([r.coord for r in refs], dtype=np.float64)
        values = np.array([r.value for r in refs], dtype=np.float64)
        if not (np.all(np.isfinite(coords)) and np.all(np.isfinite(values))):
            raise ValueError("reference points must be finite")
        if self.scale == "log10" and np.any(values <= 0.0):
            raise ValueError("a log10 axis needs positive reference values")
        order = np.argsort(coords, kind="stable")
        coords, values = coords[order], values[order]
        if np.any(np.diff(coords) <= 0.0):
            raise ValueError("reference points must have distinct coordinates")
        ts = np.log10(values) if self.scale == "log10" else values
        if np.any(np.diff(ts) == 0.0):
            raise ValueError("reference points must have distinct values")
        steps = np.sign(np.diff(ts))
        if not (np.all(steps > 0) or np.all(steps < 0)):
            raise ValueError(
                "reference values must increase (or decrease) monotonically with the coordinate"
            )
        design = np.column_stack([np.ones_like(coords), coords])
        (a, b), *_ = np.linalg.lstsq(design, ts, rcond=None)
        object.__setattr__(self, "_a", float(a))
        object.__setattr__(self, "_b", float(b))
        object.__setattr__(self, "_coords", coords)
        object.__setattr__(self, "_ts", ts)

    # ---- transform ------------------------------------------------------------------

    def _t(self, coord) -> np.ndarray:
        c = np.asarray(coord, dtype=np.float64)
        if self.fit == "lsq":
            return self._a + self._b * c
        xs, ts = self._coords, self._ts
        inside = np.interp(c, xs, ts)
        lo = ts[0] + (c - xs[0]) * (ts[1] - ts[0]) / (xs[1] - xs[0])
        hi = ts[-1] + (c - xs[-1]) * (ts[-1] - ts[-2]) / (xs[-1] - xs[-2])
        return np.where(c < xs[0], lo, np.where(c > xs[-1], hi, inside))

    def to_data(self, coord) -> np.ndarray:
        """Data value at page coordinate ``coord`` (scalar or array)."""
        t = self._t(coord)
        return np.power(10.0, t) if self.scale == "log10" else t

    def to_page(self, value) -> np.ndarray:
        """Page coordinate of data value ``value`` (inverse of :meth:`to_data`)."""
        v = np.asarray(value, dtype=np.float64)
        t = np.log10(v) if self.scale == "log10" else v
        if self.fit == "lsq":
            return (t - self._a) / self._b
        xs, ts = self._coords, self._ts
        if ts[-1] < ts[0]:
            xs, ts = xs[::-1], ts[::-1]
        inside = np.interp(t, ts, xs)
        lo = xs[0] + (t - ts[0]) * (xs[1] - xs[0]) / (ts[1] - ts[0])
        hi = xs[-1] + (t - ts[-1]) * (xs[-1] - xs[-2]) / (ts[-1] - ts[-2])
        return np.where(t < ts[0], lo, np.where(t > ts[-1], hi, inside))

    # ---- quality --------------------------------------------------------------------

    @property
    def residuals(self) -> np.ndarray:
        """Least-squares residuals ``t_ref - t_fit`` in axis-variable units, by coordinate."""
        return self._ts - (self._a + self._b * self._coords)

    @property
    def residual_max(self) -> float:
        """Largest absolute least-squares residual, in axis-variable units.

        With exactly two reference points the fit is not over-determined and this is 0.
        """
        return float(np.max(np.abs(self.residuals)))

    @property
    def residual_rms(self) -> float:
        """Root-mean-square least-squares residual, in axis-variable units."""
        return float(np.sqrt(np.mean(self.residuals**2)))

    @property
    def max_slope(self) -> float:
        """Largest ``|dt / dcoordinate|`` of the transform, axis-variable units per point."""
        if self.fit == "lsq":
            return abs(self._b)
        return float(np.max(np.abs(np.diff(self._ts) / np.diff(self._coords))))

    def width_to_axis(self, width_pt: float) -> float:
        """A page distance ``width_pt`` expressed in axis-variable units (worst case)."""
        return abs(width_pt) * self.max_slope

    def describe_error(self, err_t: float) -> str:
        """Human-readable form of an error ``err_t`` in axis-variable units."""
        if self.scale == "log10":
            percent = (10.0**err_t - 1.0) * 100.0
            return f"{err_t:.4f} decades ({percent:.1f} %)"
        return f"{err_t:.4f}"


def fit_summary(axis: Axis) -> dict:
    """Plain-data description of a calibration, for reports."""
    return {
        "scale": axis.scale,
        "fit": axis.fit,
        "references": [
            {"coord": round(r.coord, 3), "value": r.value, "source": r.source}
            for r in axis.references
        ],
        "lsq_residual_max": round(axis.residual_max, 6),
        "lsq_residual_rms": round(axis.residual_rms, 6),
    }


def axis_from_pairs(
    pairs: Sequence[tuple[float, float]], scale: Scale = "linear", fit: Fit = "lsq"
) -> Axis:
    """Shorthand: an :class:`Axis` from ``(coord, value)`` pairs."""
    return Axis(tuple(Reference(c, v) for c, v in pairs), scale=scale, fit=fit)
