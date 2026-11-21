"""Sensitometric metrics on characteristic curves (spec D1 sections 6.3-6.5, 6.9, 6.10).

Curves are ``density(log_h)`` with ``log_h`` in log10(lux*s), strictly increasing; every
metric interpolates them with one monotone PCHIP (spec 6.3 step 1).

Undefined values (spec 6.3 step 5, 6.5, 8.2). Malformed input (shapes, unsorted or
non-finite ``log_h``) raises ``ValueError``. A metric that cannot be computed from
well-formed input (curve too short or too flat, no bracket, non-finite densities from a
failed strip, a non-finite ``d_min``) returns NaN and never raises, so callers can count it
as a failed metric. The ``*_detail`` variants return the same value with the reason as a
:class:`MetricValue`; :func:`iso_speed` carries it in :attr:`SpeedResult.reason`.

Base plus fog (spec 6.3 step 2, 6.5): every comparison with a Kodak characteristic curve uses
the same operator on both sides, ``D_min* = D_ref(x_left)`` with ``x_left`` the reference
curve's leftmost log H (:func:`d_min_star`). Comparisons with Kodak CI(t) curves and time
tables, and prediction-only conditions, use the unexposed patch of the strip.
"""

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from darkroom_atlas._frozen import readonly_array
from darkroom_atlas.units import photon_energy_j
from darkroom_atlas.validate.sensitometer import StripResult

__all__ = [
    "MetricValue",
    "SpeedResult",
    "average_gradient",
    "average_gradient_detail",
    "chord_contrast_index",
    "chord_contrast_index_detail",
    "ci_at_times",
    "contrast_index",
    "contrast_index_detail",
    "d_min_star",
    "energy_exposure_erg_cm2",
    "iso_bin",
    "iso_speed",
    "speed_point",
    "time_to_ci",
    "time_to_ci_detail",
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
# Spec 6.9 step 3: E_b [erg/cm2] = N_b [photons/um2] * (h c / lambda_b) [J] * 1e7 * 1e8.
_ERG_PER_J = 1e7  # spec 6.9 step 3
_UM2_PER_CM2 = 1e8  # spec 6.9 step 3

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
class MetricValue:
    """A metric value and, when it is NaN, the reason (spec 6.5, 8.2)."""

    value: float
    reason: str = ""


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


def _nodes(log_h: np.ndarray, density: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Check a curve's shape and abscissa; raise ``ValueError`` on malformed input."""
    x = np.asarray(log_h, dtype=np.float64)
    y = np.asarray(density, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1 or x.shape != y.shape:
        raise ValueError(
            f"log_h and density must be 1-D of the same length, got {x.shape} and {y.shape}"
        )
    if x.size < 2:
        raise ValueError(f"a curve needs at least 2 points, got {x.size}")
    if not np.all(np.isfinite(x)):
        raise ValueError("log_h must be finite")
    if not np.all(np.diff(x) > 0.0):
        raise ValueError("log_h must be strictly increasing")
    return x, y


def _curve(log_h: np.ndarray, density: np.ndarray) -> tuple[PchipInterpolator | None, str]:
    """Monotone PCHIP of a curve, or ``(None, reason)`` when its densities are not finite."""
    x, y = _nodes(log_h, density)
    bad = int(np.count_nonzero(~np.isfinite(y)))
    if bad:
        return None, f"the density has {bad} non-finite value(s)"
    return PchipInterpolator(x, y, extrapolate=False), ""


def _level_reason(value: float, what: str) -> str:
    return "" if math.isfinite(value) else f"{what} is not finite ({value})"


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


def d_min_star(
    log_h_ref: np.ndarray,
    density_ref: np.ndarray,
    log_h_sim: np.ndarray,
    density_sim: np.ndarray,
) -> tuple[float, float]:
    """``(reference, simulated)`` D_min*, both read at ``x_left = log_h_ref[0]`` (spec 6.3).

    The reference side is ``density_ref[0]``; the simulated side is the PCHIP of the simulated
    curve at the same ``x_left``. The simulated value is NaN when ``x_left`` lies outside the
    simulated log H range or the simulated densities are not finite. The reference curve must
    be finite.
    """
    x_ref, y_ref = _nodes(log_h_ref, density_ref)
    if not np.all(np.isfinite(y_ref)):
        raise ValueError("density_ref must be finite")
    x_left = float(x_ref[0])
    curve, _ = _curve(log_h_sim, density_sim)
    sim = math.nan if curve is None else float(curve(x_left))
    return float(y_ref[0]), sim


def _speed_point_detail(
    log_h: np.ndarray, density: np.ndarray, base_fog: float
) -> tuple[float, float, str]:
    curve, reason = _curve(log_h, density)
    base_fog = float(base_fog)
    reason = reason or _level_reason(base_fog, "base_fog")
    if curve is None or reason:
        return math.nan, math.nan, reason
    target = base_fog + ISO_SPEED_POINT_RISE
    x_m = _level_crossing(curve, target)
    if x_m is None:
        return math.nan, math.nan, f"the curve does not rise through base_fog + 0.10 = {target:.3f}"
    x_far = x_m + ISO_LOG_H_RANGE
    if x_far > curve.x[-1]:
        return x_m, math.nan, "the curve ends before log H_m + 1.30"
    return x_m, float(curve(x_far)) - target, ""


def speed_point(log_h: np.ndarray, density: np.ndarray, base_fog: float) -> tuple[float, float]:
    """``(log H_m, delta_d)`` of one curve (spec 6.3 steps 3-4); NaN where undefined.

    ``log H_m`` is where the density first reaches ``base_fog + 0.10``; ``delta_d`` is
    ``D(log H_m + 1.30) - D(log H_m)``, NaN when the curve ends before ``log H_m + 1.30``.
    """
    x_m, delta_d, _ = _speed_point_detail(log_h, density, base_fog)
    return x_m, delta_d


def average_gradient_detail(log_h: np.ndarray, density: np.ndarray, base_fog: float) -> MetricValue:
    """:func:`average_gradient` with the reason when it is NaN."""
    _, delta_d, reason = _speed_point_detail(log_h, density, base_fog)
    return MetricValue(delta_d / ISO_LOG_H_RANGE, reason)


def average_gradient(log_h: np.ndarray, density: np.ndarray, base_fog: float) -> float:
    """ISO average gradient ``G_ISO = delta_d / 1.30`` (spec 6.4); NaN when undefined."""
    return average_gradient_detail(log_h, density, base_fog).value


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
    1-2). ``base_fog`` overrides it with a scalar or one value per development time; pass
    ``D_min*`` (:func:`d_min_star`) when comparing with Kodak curves, which gives log H_m*.
    The series is interpolated linearly in time to ``delta_d = 0.80`` (step 5, ISO 6 Annex
    A). A non-finite ``base_fog`` entry makes that development time undefined.
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
    if fog.shape != (n_t,):
        raise ValueError(f"base_fog must have shape () or ({n_t},), got {fog.shape}")

    points = [_speed_point_detail(strip.log_h, curves[i], float(fog[i])) for i in range(n_t)]
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
            causes = sorted({p[2] for p in points if p[2]})
            return _nan_speed(
                "no development time has a speed point and a point 1.30 log H beyond it"
                + (f" ({'; '.join(causes)})" if causes else "")
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


_NO_CENTRE = "no collinear 0.2/2.2 construction exists (curve too short or too flat)"


def contrast_index_detail(log_h: np.ndarray, density: np.ndarray, d_min: float) -> MetricValue:
    """:func:`contrast_index` with the reason when it is NaN."""
    curve, reason = _curve(log_h, density)
    d_min = float(d_min)
    reason = reason or _level_reason(d_min, "d_min")
    if curve is None or reason:
        return MetricValue(math.nan, reason)
    x_lo, x_hi = float(curve.x[0]), float(curve.x[-1])
    centres = np.append(np.arange(x_lo - CI_CENTRE_MARGIN, x_hi, _CENTRE_SCAN_DX), x_hi)
    previous: tuple[float, float] | None = None
    for c in centres:
        marks = _ci_marks(curve, float(c), d_min)
        if marks is None:
            previous = None
            continue
        if marks[0] == 0.0 and _slope(curve, marks[1], marks[2]) > 0.0:
            return MetricValue(_slope(curve, marks[1], marks[2]))
        if previous is not None and previous[1] * marks[0] < 0.0:
            return _solve_centre(curve, d_min, previous[0], float(c))
        previous = (float(c), marks[0])
    return MetricValue(math.nan, _NO_CENTRE)


def contrast_index(log_h: np.ndarray, density: np.ndarray, d_min: float) -> float:
    """Kodak H-740 contrast index by the collinear 0.2/2.2 construction (spec 6.5).

    The centre C = (x_C, d_min) slides along the fog line until C and the marks P1, P2 at
    distances 0.2 and 2.2 lie on one line; CI is the slope P1 P2. ``d_min`` is ``D_min*``
    (:func:`d_min_star`) when the result is compared with a CI read off a Kodak
    characteristic curve, and the unexposed patch otherwise (spec 6.5). NaN when no centre
    works (curve too short or too flat); a non-positive slope is also NaN, since a curve
    lying on the fog line makes every centre trivially collinear.
    """
    return contrast_index_detail(log_h, density, d_min).value


def _misalignment(curve: PchipInterpolator, d_min: float, x_c: float) -> float:
    marks = _ci_marks(curve, x_c, d_min)
    if marks is None:
        raise _UndefinedMark
    return marks[0]


class _UndefinedMark(Exception):
    pass


def _solve_centre(curve: PchipInterpolator, d_min: float, lo: float, hi: float) -> MetricValue:
    try:
        x_c = brentq(lambda c: _misalignment(curve, d_min, c), lo, hi, xtol=_XTOL)
    except _UndefinedMark:
        return MetricValue(math.nan, _NO_CENTRE)
    marks = _ci_marks(curve, float(x_c), d_min)
    if marks is None:
        return MetricValue(math.nan, _NO_CENTRE)
    slope = _slope(curve, marks[1], marks[2])
    if not slope > 0.0:
        return MetricValue(math.nan, "the collinear construction gives a non-positive slope")
    return MetricValue(slope)


def _slope(curve: PchipInterpolator, x1: float, x2: float) -> float:
    return (float(curve(x2)) - float(curve(x1))) / (x2 - x1)


def chord_contrast_index_detail(
    log_h: np.ndarray, density: np.ndarray, d_min: float
) -> MetricValue:
    """:func:`chord_contrast_index` with the reason when it is NaN."""
    curve, reason = _curve(log_h, density)
    d_min = float(d_min)
    reason = reason or _level_reason(d_min, "d_min")
    if curve is None or reason:
        return MetricValue(math.nan, reason)
    level = d_min + CHORD_RISE
    x_a = _level_crossing(curve, level)
    if x_a is None:
        return MetricValue(math.nan, f"the curve does not rise through d_min + 0.10 = {level:.3f}")
    if x_a + CHORD_LENGTH > curve.x[-1]:
        return MetricValue(math.nan, "the curve ends less than 2.0 log H after d_min + 0.10")
    return MetricValue((float(curve(x_a + CHORD_LENGTH)) - level) / CHORD_LENGTH)


def chord_contrast_index(log_h: np.ndarray, density: np.ndarray, d_min: float) -> float:
    """Chord cross-check of the contrast index (spec 6.5, CAL 1.3.2).

    Slope from the point 0.10 above ``d_min`` to the point 2.0 log H further on; NaN when the
    curve never rises 0.10 or ends too soon. Not used for fitting or decisions.
    """
    return chord_contrast_index_detail(log_h, density, d_min).value


def _series(times_s: np.ndarray, ci: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(times_s, dtype=np.float64)
    c = np.asarray(ci, dtype=np.float64)
    if t.ndim != 1 or c.shape != t.shape:
        raise ValueError(f"times_s and ci must be 1-D of the same length, got {t.shape}, {c.shape}")
    if t.size < 1:
        raise ValueError("times_s must not be empty")
    if not (np.all(np.isfinite(t)) and np.all(np.diff(t) > 0.0)):
        raise ValueError("times_s must be finite and strictly increasing")
    return t, c


def _finite_runs(c: np.ndarray) -> list[tuple[int, int]]:
    """Maximal runs ``[start, stop)`` of consecutive finite values of ``c``."""
    finite = np.concatenate(([False], np.isfinite(c), [False]))
    edges = np.flatnonzero(np.diff(finite.astype(np.int8)))
    return [(int(a), int(b)) for a, b in zip(edges[::2], edges[1::2], strict=True)]


def time_to_ci_detail(times_s: np.ndarray, ci: np.ndarray, target: float) -> MetricValue:
    """:func:`time_to_ci` with the reason when it is NaN."""
    t, c = _series(times_s, ci)
    target = float(target)
    if not math.isfinite(target):
        raise ValueError(f"target must be finite, got {target}")
    runs = _finite_runs(c)
    if not runs:
        return MetricValue(math.nan, "CI(t) has no finite value")
    below_before_gap = False
    for a, b in runs:
        if c[a] == target:
            return MetricValue(float(t[a]))
        if c[a] > target:
            if below_before_gap:
                return MetricValue(
                    math.nan, f"CI(t) reaches {target} inside a run of undefined CI values"
                )
            return MetricValue(math.nan, f"CI(t) is above {target} at the first defined time")
        if b - a >= 2:
            curve = PchipInterpolator(t[a:b], c[a:b], extrapolate=False)
            hit = _first_rise(lambda x, curve=curve: curve(x) - target, t[a:b])
            if hit is not None:
                return MetricValue(hit)
        below_before_gap = True
    return MetricValue(math.nan, f"CI(t) does not reach {target} within the series")


def time_to_ci(times_s: np.ndarray, ci: np.ndarray, target: float) -> float:
    """Development time (s) at which CI(t) first reaches ``target`` (spec 6.10).

    Report only, and for the "developed to CI 0.56" conditions of spec 6.6-6.9; it is not a
    fitting residual or a threshold (table cells use :func:`ci_at_times`). CI(t) is
    interpolated with monotone PCHIP over each run of consecutive finite values; NaN when the
    target is not reached inside the series or is crossed where CI is undefined. Extending the
    series (6.10) is the caller's job.
    """
    return time_to_ci_detail(times_s, ci, target).value


def ci_at_times(times_s: np.ndarray, ci: np.ndarray, at_s: np.ndarray) -> np.ndarray:
    """CI at the times ``at_s`` (s) by monotone PCHIP through CI(t) (spec 6.10, 7.2).

    Used to score table cells: ``CI_sim(t_cell) - CI_aim``. The result has the shape of
    ``at_s`` and is NaN where ``at_s`` lies outside ``[times_s[0], times_s[-1]]`` or in an
    interval touching an undefined (non-finite) CI value; undefined values are never bridged.
    Read-only float64.
    """
    t, c = _series(times_s, ci)
    at = np.asarray(at_s, dtype=np.float64)
    if not np.all(np.isfinite(at)):
        raise ValueError("at_s must be finite")
    out = np.full(at.shape, math.nan)
    for a, b in _finite_runs(c):
        inside = (at >= t[a]) & (at <= t[b - 1])
        if not np.any(inside):
            continue
        if b - a == 1:
            out[inside] = c[a]
        else:
            out[inside] = PchipInterpolator(t[a:b], c[a:b], extrapolate=False)(at[inside])
    return readonly_array(out, np.float64)


def energy_exposure_erg_cm2(photons_um2: float | np.ndarray, wavelength_nm: float | np.ndarray):
    """Energy exposure (erg/cm2) of ``photons_um2`` photons per um2 at ``wavelength_nm``.

    Spec 6.9 step 3: ``E_b = N_b * (h c / lambda_b) * 1e7 * 1e8`` with the photon energy in J
    from D0's :func:`~darkroom_atlas.units.photon_energy_j` (not the eV-valued E_ph,b of
    spec 2.1.2). At 550 nm, 1 erg/cm2 is about 2769 photons/um2 (LAT 3.3).
    """
    n = np.asarray(photons_um2, dtype=np.float64)
    energy = n * photon_energy_j(wavelength_nm) * _ERG_PER_J * _UM2_PER_CM2
    return float(energy) if np.ndim(energy) == 0 else np.asarray(energy, dtype=np.float64)
