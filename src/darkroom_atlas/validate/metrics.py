"""Sensitometric metrics on characteristic curves (spec D1 sections 6.3-6.5, 6.10).

Curves are ``density(log_h)`` with ``log_h`` in log10(lux*s), strictly increasing; every
metric interpolates them with one monotone PCHIP (spec 6.3 step 1). A metric that cannot be
computed (curve too short or too flat, no bracket) returns NaN; :func:`iso_speed` also says
why in :attr:`SpeedResult.reason`.
"""

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from darkroom_atlas.validate.sensitometer import StripResult

__all__ = [
    "SpeedResult",
    "average_gradient",
    "chord_contrast_index",
    "contrast_index",
    "iso_bin",
    "iso_speed",
    "speed_point",
    "time_to_ci",
]

# ISO 6:1993 5.6 (spec 6.3): the speed point lies 0.10 above base plus fog, the gradient is
# taken over 1.30 log H from it, and development is chosen so that the density rise is 0.80.
ISO_SPEED_POINT_RISE = 0.10
ISO_LOG_H_RANGE = 1.30
ISO_DELTA_D = 0.80
# Kodak H-740 contrast index (spec 6.5): the marks lie 0.2 and 2.2 from the centre.
CI_INNER = 0.2
CI_OUTER = 2.2
# Spec 6.5: the centre is searched on [x_lo - 1, x_hi]; crossings are bracketed on a
# Delta x = 0.002 sampling of the PCHIP curve and refined with brentq.
CI_CENTRE_MARGIN = 1.0
CI_SAMPLE_DX = 0.002
# Step of the scan that brackets the collinear centre before brentq (implementation choice).
_CENTRE_SCAN_DX = 0.01
# Chord cross-check (spec 6.5, CAL 1.3.2): 0.10 above D_min, then 2.0 log H further on.
CHORD_RISE = 0.10
CHORD_LENGTH = 2.0
_XTOL = 1e-13

# ISO 6:1993 Table 1 (spec 6.3 step 6): bins 0.10 wide in log H_m, the ISO 400 bin is
# -2.75 ... -2.66 (two-decimal values); one bin per third of a stop.
_ISO_SERIES = (
    0.8, 1, 1.2, 1.6, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64, 80, 100, 125,
    160, 200, 250, 320, 400, 500, 640, 800, 1000, 1250, 1600, 2000, 2500, 3200, 4000, 5000,
    6400, 8000, 10000,
)  # fmt: skip
_ISO_400_INDEX = _ISO_SERIES.index(400)
_ISO_400_LOW_CENTI = -275  # lower end of the ISO 400 bin, in hundredths of log H_m


@dataclass(frozen=True)
class SpeedResult:
    """ISO 6 speed of a development series (spec 6.3).

    ``log_h_m`` (log10 lux*s) and ``delta_d`` are interpolated to the development time at
    which the density rise over 1.30 log H is 0.80, so ``delta_d`` is 0.80 when defined.
    ``speed = 0.80 / H_m``, ``speed_degrees = 1 + 10 log10(0.80 / H_m)``; ``iso_bin`` is the
    ISO 6 Table 1 arithmetic speed (``"400"``) or ``None``. All NaN with ``reason`` set when
    the speed cannot be determined.
    """

    log_h_m: float
    delta_d: float
    speed: float
    speed_degrees: float
    iso_bin: str | None
    reason: str = ""


def _curve(log_h: np.ndarray, density: np.ndarray) -> PchipInterpolator:
    x = np.asarray(log_h, dtype=np.float64)
    y = np.asarray(density, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1 or x.shape != y.shape:
        raise ValueError(
            f"log_h and density must be 1-D of the same length, got {x.shape} and {y.shape}"
        )
    if x.size < 2:
        raise ValueError(f"a curve needs at least 2 points, got {x.size}")
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        raise ValueError("log_h and density must be finite")
    if not np.all(np.diff(x) > 0.0):
        raise ValueError("log_h must be strictly increasing")
    return PchipInterpolator(x, y, extrapolate=False)


def _check_level(value: float, what: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{what} must be finite, got {value}")
    return value


def _first_rise(g: Callable[[np.ndarray], np.ndarray], xs: np.ndarray) -> float | None:
    """First x in [xs[0], xs[-1]] where g passes from below 0 to at least 0, else None.

    ``xs`` are increasing sample points; the bracket is refined with brentq. ``None`` also when
    g is already >= 0 at the first sample (the crossing is not inside the range).
    """
    values = g(xs)
    if not values[0] < 0.0:
        return None
    hits = np.flatnonzero((values[:-1] < 0.0) & (values[1:] >= 0.0))
    if hits.size == 0:
        return None
    i = int(hits[0])
    if values[i + 1] == 0.0:
        return float(xs[i + 1])
    return float(brentq(lambda t: float(g(np.array([t]))[0]), xs[i], xs[i + 1], xtol=_XTOL))


def _level_crossing(curve: PchipInterpolator, level: float) -> float | None:
    """First log H where the curve reaches ``level`` from below (PCHIP is monotone per
    interval, so the nodes bracket every crossing)."""
    return _first_rise(lambda x: curve(x) - level, curve.x)


def speed_point(log_h: np.ndarray, density: np.ndarray, base_fog: float) -> tuple[float, float]:
    """``(log H_m, delta_d)`` of one curve (spec 6.3 steps 3-4); NaN where undefined.

    ``log H_m`` is where the density first reaches ``base_fog + 0.10``; ``delta_d`` is
    ``D(log H_m + 1.30) - D(log H_m)``, NaN when the curve ends before ``log H_m + 1.30``.
    """
    curve = _curve(log_h, density)
    base_fog = _check_level(base_fog, "base_fog")
    target = base_fog + ISO_SPEED_POINT_RISE
    x_m = _level_crossing(curve, target)
    if x_m is None:
        return math.nan, math.nan
    x_far = x_m + ISO_LOG_H_RANGE
    if x_far > curve.x[-1]:
        return x_m, math.nan
    return x_m, float(curve(x_far)) - target


def average_gradient(log_h: np.ndarray, density: np.ndarray, base_fog: float) -> float:
    """ISO average gradient ``G_ISO = delta_d / 1.30`` (spec 6.4); NaN when undefined."""
    return speed_point(log_h, density, base_fog)[1] / ISO_LOG_H_RANGE


def iso_bin(log_h_m: float) -> str | None:
    """ISO 6 Table 1 arithmetic speed for ``log_h_m``, or ``None`` outside the table."""
    if not math.isfinite(log_h_m):
        return None
    centi = round(round(log_h_m, 2) * 100.0)
    steps = (centi - _ISO_400_LOW_CENTI) // 10  # positive: more exposure, slower
    index = _ISO_400_INDEX - steps
    if not 0 <= index < len(_ISO_SERIES):
        return None
    return f"{_ISO_SERIES[index]:g}"


def _nan_speed(reason: str) -> SpeedResult:
    return SpeedResult(math.nan, math.nan, math.nan, math.nan, None, reason)


def iso_speed(
    strip: StripResult,
    *,
    density: Literal["type2", "visual"] = "type2",
    base_fog: float | np.ndarray | None = None,
) -> SpeedResult:
    """ISO 6:1993 speed of a development series (spec 6.3).

    By default it uses Type 2 densities and the unexposed patch as base plus fog (6.3 steps
    1-2). ``base_fog`` overrides it with a scalar or one value per development time, e.g. the
    left-end density ``D_min*`` when comparing with Kodak curves. The series is interpolated
    linearly in time to ``delta_d = 0.80`` (step 5, ISO 6 Annex A).
    """
    if density not in ("type2", "visual"):
        raise ValueError(f"density must be 'type2' or 'visual', got {density!r}")
    curves = strip.type2 if density == "type2" else strip.visual
    n_t = strip.times_s.size
    if base_fog is None:
        fog = np.asarray(strip.base_fog_type2 if density == "type2" else strip.base_fog_visual)
    else:
        fog = np.asarray(base_fog, dtype=np.float64)
        if fog.ndim == 0:
            fog = np.full(n_t, float(fog))
    if fog.shape != (n_t,) or not np.all(np.isfinite(fog)):
        raise ValueError(f"base_fog must be finite with shape () or ({n_t},), got {fog.shape}")

    points = [speed_point(strip.log_h, curves[i], float(fog[i])) for i in range(n_t)]
    log_h_m = np.array([p[0] for p in points])
    delta_d = np.array([p[1] for p in points])
    times = np.asarray(strip.times_s)

    found: tuple[float, float] | None = None
    for i in range(n_t):
        if delta_d[i] == ISO_DELTA_D:
            found = (float(log_h_m[i]), float(delta_d[i]))
            break
        if i + 1 < n_t and np.all(np.isfinite(delta_d[i : i + 2])):
            a, b = delta_d[i] - ISO_DELTA_D, delta_d[i + 1] - ISO_DELTA_D
            if a * b < 0.0:
                w = a / (a - b)
                t_star = times[i] + w * (times[i + 1] - times[i])
                x_star = float(np.interp(t_star, times[i : i + 2], log_h_m[i : i + 2]))
                found = (x_star, ISO_DELTA_D)
                break
    if found is None:
        finite = delta_d[np.isfinite(delta_d)]
        if finite.size == 0:
            return _nan_speed(
                "no development time has a speed point and a point 1.30 log H beyond it"
            )
        return _nan_speed(
            f"the density rise over 1.30 log H ({finite.min():.3f} ... {finite.max():.3f}) "
            f"does not bracket {ISO_DELTA_D} between adjacent development times"
        )
    x_m, d_d = found
    h_m = 10.0**x_m
    return SpeedResult(
        log_h_m=x_m,
        delta_d=d_d,
        speed=ISO_DELTA_D / h_m,
        speed_degrees=1.0 + 10.0 * math.log10(ISO_DELTA_D / h_m),
        iso_bin=iso_bin(x_m),
    )


def _mark(curve: PchipInterpolator, x_c: float, d_min: float, r: float) -> float | None:
    """First log H > x_c where the curve's distance from (x_c, d_min) reaches r from below."""
    lo, hi = max(x_c, float(curve.x[0])), float(curve.x[-1])
    if not lo < hi:
        return None
    xs = np.union1d(np.append(np.arange(lo, hi, CI_SAMPLE_DX), hi), curve.x[curve.x > lo])
    return _first_rise(lambda x: np.hypot(x - x_c, curve(x) - d_min) - r, xs)


def _ci_marks(
    curve: PchipInterpolator, x_c: float, d_min: float
) -> tuple[float, float, float] | None:
    """(misalignment angle, x1, x2) for centre x_c, or None when a mark is missing."""
    x1 = _mark(curve, x_c, d_min, CI_INNER)
    x2 = _mark(curve, x_c, d_min, CI_OUTER)
    if x1 is None or x2 is None:
        return None
    a1 = math.atan2(float(curve(x1)) - d_min, x1 - x_c)
    a2 = math.atan2(float(curve(x2)) - d_min, x2 - x_c)
    return a1 - a2, x1, x2


def contrast_index(log_h: np.ndarray, density: np.ndarray, d_min: float) -> float:
    """Kodak H-740 contrast index by the collinear 0.2/2.2 construction (spec 6.5).

    The centre C = (x_C, d_min) slides along the fog line until C and the marks P1, P2 at
    distances 0.2 and 2.2 lie on one line; CI is the slope P1 P2. NaN when no centre works
    (curve too short or too flat); a non-positive slope is also NaN, since a curve lying on
    the fog line makes every centre trivially collinear.
    """
    curve = _curve(log_h, density)
    d_min = _check_level(d_min, "d_min")
    x_lo, x_hi = float(curve.x[0]), float(curve.x[-1])
    centres = np.append(np.arange(x_lo - CI_CENTRE_MARGIN, x_hi, _CENTRE_SCAN_DX), x_hi)
    previous: tuple[float, float] | None = None
    for c in centres:
        marks = _ci_marks(curve, float(c), d_min)
        if marks is None:
            previous = None
            continue
        if marks[0] == 0.0 and _slope(curve, marks[1], marks[2]) > 0.0:
            return _slope(curve, marks[1], marks[2])
        if previous is not None and previous[1] * marks[0] < 0.0:
            return _solve_centre(curve, d_min, previous[0], float(c))
        previous = (float(c), marks[0])
    return math.nan


def _misalignment(curve: PchipInterpolator, d_min: float, x_c: float) -> float:
    marks = _ci_marks(curve, x_c, d_min)
    if marks is None:
        raise _UndefinedMark
    return marks[0]


class _UndefinedMark(Exception):
    pass


def _solve_centre(curve: PchipInterpolator, d_min: float, lo: float, hi: float) -> float:
    try:
        x_c = brentq(lambda c: _misalignment(curve, d_min, c), lo, hi, xtol=_XTOL)
    except _UndefinedMark:
        return math.nan
    marks = _ci_marks(curve, float(x_c), d_min)
    if marks is None:
        return math.nan
    slope = _slope(curve, marks[1], marks[2])
    return slope if slope > 0.0 else math.nan


def _slope(curve: PchipInterpolator, x1: float, x2: float) -> float:
    return (float(curve(x2)) - float(curve(x1))) / (x2 - x1)


def chord_contrast_index(log_h: np.ndarray, density: np.ndarray, d_min: float) -> float:
    """Chord cross-check of the contrast index (spec 6.5, CAL 1.3.2).

    Slope from the point 0.10 above ``d_min`` to the point 2.0 log H further on; NaN when the
    curve never rises 0.10 or ends too soon. Not used for fitting or decisions.
    """
    curve = _curve(log_h, density)
    d_min = _check_level(d_min, "d_min")
    level = d_min + CHORD_RISE
    x_a = _level_crossing(curve, level)
    if x_a is None or x_a + CHORD_LENGTH > curve.x[-1]:
        return math.nan
    return (float(curve(x_a + CHORD_LENGTH)) - level) / CHORD_LENGTH


def time_to_ci(times_s: np.ndarray, ci: np.ndarray, target: float) -> float:
    """Development time (s) at which CI(t) first reaches ``target`` (spec 6.10).

    CI(t) is interpolated with monotone PCHIP through its finite values; NaN when the target
    is not reached inside the series. Extending the series (6.10) is the caller's job.
    """
    t = np.asarray(times_s, dtype=np.float64)
    c = np.asarray(ci, dtype=np.float64)
    if t.ndim != 1 or c.shape != t.shape:
        raise ValueError(f"times_s and ci must be 1-D of the same length, got {t.shape}, {c.shape}")
    if not (np.all(np.isfinite(t)) and np.all(np.diff(t) > 0.0)):
        raise ValueError("times_s must be finite and strictly increasing")
    target = _check_level(target, "target")
    keep = np.isfinite(c)
    t, c = t[keep], c[keep]
    if t.size == 0:
        return math.nan
    if c[0] == target:
        return float(t[0])
    if t.size < 2:
        return math.nan
    curve = PchipInterpolator(t, c, extrapolate=False)
    hit = _first_rise(lambda x: curve(x) - target, t)
    return math.nan if hit is None else hit
