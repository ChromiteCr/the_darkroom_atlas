"""Statistical grain mode over the whole frame (spec D1 2.7.1).

Per pixel and size class ``c``, from the developed fields (spec 2.0.3):

* ``mu_c = nu_c A_pix F_c`` developed grains on average (``A_pix = p^2``);
* ``x_c = X_c / F_c`` and ``<x^2>_c = Y_c / F_c``, the developed-fraction moments of the
  initiated grains, with ``X_c = silver_c / (1e12 m_Ag,c nu_c)``;
* ``phi = r_CP CP_ref silver_physical / sum_c CP_c silver_c`` enlarges the developed disks
  (``phi = 0`` and physical silver as a mean factor ``10^-(r_CP CP_ref silver_physical)``
  when there is no chemical silver);
* ``q_c = a_c (1 + phi) x_c / A_pix``.

Classes with ``p >= 4 r_c``, ``r_c = sqrt(a_c / pi)``, are *counted*: ``N_c ~ Poisson(mu_c)``
from :func:`~darkroom_atlas.grain.poisson.poisson_inverse` of one uniform per pixel (stream
``grain/emulsion/<c>``, ``c`` the 0-based class index). The others are *mean-only*: they
contribute the factor ``exp(-mu_c q_c)`` and no random draw. Where physical silver enlarges
a counted class to ``q_c >= 1`` (``(1 - q_c)^N`` would no longer average to
``exp(-mu_c q_c)``), that class is mean-only in that pixel; its uniform is still drawn, so the
streams do not change. The pixel transmittance is

``T_pix = prod_counted (1 - q_c)^N_c * prod_mean_only exp(-mu_c q_c) * (1 + eps)``,

truncated to [0, 1], with ``eps = sqrt(Var(eps)) z`` (``z = ndtri(u)``, stream
``grain/emulsion/overlap``) and ``Var(eps) = max(Var_bool - Var_count, 0) / (T^2 +
Var_count)``, ``Var_count = T^2 (exp(sum_counted mu_c q_c^2) - 1)``: the untruncated pixel
variance is the Boolean-model variance ``Var_bool`` of
:func:`~darkroom_atlas.grain.overlap.boolean_excess_variance` (intensity ``nu_c F_c``,
radius ``r~_c = sqrt(a_c (1 + phi) <x^2>_c^(1/2) / pi)``, all classes). Its mean is
``T = exp(-sum_c mu_c q_c) = 10^-D_Ag``, Nutting's transmittance (spec 2.5).

Random numbers are keyed by ``(state.seed, stream, absolute pixel)`` and every per-pixel
value is computed element-wise (sums over classes as explicit adds), so a pixel's sample
does not depend on the window it is computed in (spec 3.6).
"""

import itertools
import math
import warnings
from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
from scipy.special import ndtri

from darkroom_atlas.grain.overlap import boolean_excess_variance
from darkroom_atlas.grain.poisson import poisson_inverse
from darkroom_atlas.imaging.silver import (
    LN10,
    chemical_silver_density,
    imaging_constants,
    resolve_window,
)
from darkroom_atlas.rng import pixel_uniform
from darkroom_atlas.state import EmulsionState

__all__ = [
    "APERTURE_BIAS_COEFFICIENT",
    "APERTURE_BIAS_WARNING",
    "COUNTED_RADIUS_FACTOR",
    "OVERLAP_STREAM",
    "GrainBlockStats",
    "GrainStatistics",
    "aperture_bias",
    "class_stream",
    "counted_classes",
    "epsilon_variance",
    "exact_expansion",
    "grain_meta",
    "grain_statistics",
    "statistical_transmittance",
    "statistical_transmittance_with_stats",
]

COUNTED_RADIUS_FACTOR = 4.0
"""A class is counted when the pixel pitch is at least 4 r_c (spec 2.7.1)."""

APERTURE_BIAS_COEFFICIENT = 1.15
"""``b_ap = 1.15 r~_eff / D_ap``, the aperture-variance bias of the white-noise pixels
(spec 2.7.1; constant-radius numerical check)."""

APERTURE_BIAS_WARNING = 0.05
"""``b_ap`` above which a warning is raised (spec 2.7.1)."""

OVERLAP_STREAM = "grain/emulsion/overlap"
"""Random stream of the overlap term ``z`` (spec 2.7.1)."""

_STREAM_PREFIX = "grain/emulsion"

# Uniforms are clipped to [2^-53, 1 - 2^-53] before ndtri so z stays finite (as in poisson.py).
_U_EPS = 2.0**-53

# Pixels per boolean_excess_variance call: bounds its (pixels, 48) work arrays.
_VAR_CHUNK_PX = 8192

# Relative tolerance on F <= 1, X <= F and Y <= X for rounding in the developed fields.
_FRACTION_TOL = 1e-9

# Upper bound on the terms of an exact float expansion (a finite float64 sum needs < 45).
_EXPANSION_MAX_TERMS = 64


def class_stream(c: int) -> str:
    """Random stream of the Poisson count of class ``c`` (0-based): ``grain/emulsion/<c>``."""
    return f"{_STREAM_PREFIX}/{int(c)}"


def counted_classes(developed_areas_um2: np.ndarray, pixel_um: float) -> np.ndarray:
    """Which classes enter the Poisson count product: ``p >= 4 sqrt(a_c / pi)`` (spec 2.7.1).

    ``developed_areas_um2`` is ``a_c`` (``(K,)``, um^2), ``pixel_um`` the pixel pitch.
    Returns a ``(K,)`` bool array.
    """
    a = np.asarray(developed_areas_um2, dtype=np.float64)
    if a.ndim != 1:
        raise ValueError(f"developed_areas_um2 must have shape (K,), got {a.shape}")
    if not np.all(np.isfinite(a)) or np.any(a < 0.0):
        raise ValueError("developed_areas_um2 must be finite and non-negative")
    if isinstance(pixel_um, bool) or not np.isfinite(pixel_um) or not pixel_um > 0.0:
        raise ValueError(f"pixel_um must be positive and finite, got {pixel_um!r}")
    return float(pixel_um) >= COUNTED_RADIUS_FACTOR * np.sqrt(a / math.pi)


def _count_variance(
    t: np.ndarray, mu: np.ndarray, q: np.ndarray, counted: np.ndarray
) -> np.ndarray:
    """``Var_count = T^2 (exp(sum_{c counted} mu_c q_c^2) - 1)``, sums as element-wise adds.

    Where ``expm1`` overflows (or meets ``T^2 = 0``) the value is ``exp(2 ln T + s) - T^2``.
    """
    s = np.zeros(mu.shape[:-1], dtype=np.float64)
    for c in np.flatnonzero(counted):
        s = s + mu[..., c] * (q[..., c] * q[..., c])
    t2 = t * t
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        out = t2 * np.expm1(s)
        bad = ~np.isfinite(out)
        if np.any(bad):
            out = np.where(bad, np.exp(2.0 * np.log(t) + s) - t2, out)
    return out


def epsilon_variance(
    var_bool: np.ndarray,
    mean_t: np.ndarray,
    mu: np.ndarray,
    q: np.ndarray,
    counted: np.ndarray,
) -> np.ndarray:
    """``Var(eps) = max(Var_bool - Var_count, 0) / (T^2 + Var_count)`` (spec 2.7.1).

    ``Var_count = T^2 (exp(sum_{c counted} mu_c q_c^2) - 1)`` with ``q_c`` from the mean
    developed fraction. ``mu`` and ``q`` are ``(..., K)``, ``counted`` is ``(K,)`` bool,
    ``var_bool`` and ``mean_t`` broadcast to ``(...)``. Where ``T^2 + Var_count`` is 0 the
    result is 0.
    """
    mu = np.asarray(mu, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    counted = np.asarray(counted, dtype=bool)
    if mu.shape != q.shape or mu.ndim == 0:
        raise ValueError(f"mu and q must have the same shape (..., K), got {mu.shape}, {q.shape}")
    if counted.shape != mu.shape[-1:]:
        raise ValueError(f"counted must have shape ({mu.shape[-1]},), got {counted.shape}")
    t = np.asarray(mean_t, dtype=np.float64)
    v_bool = np.asarray(var_bool, dtype=np.float64)
    for name, value in (("var_bool", v_bool), ("mean_t", t), ("mu", mu), ("q", q)):
        if not np.all(np.isfinite(value)) or np.any(value < 0.0):
            raise ValueError(f"{name} must be finite and non-negative")
    t2 = t * t
    var_count = _count_variance(t, mu, q, counted)
    denominator = t2 + var_count
    excess = np.maximum(v_bool - var_count, 0.0)
    safe = np.where(denominator > 0.0, denominator, 1.0)
    return np.where(denominator > 0.0, excess / safe, 0.0)


@dataclass(frozen=True, eq=False)
class GrainStatistics:
    """Per-pixel grain statistics of a window (spec 2.7.1); arrays are read-only.

    ``mean_t`` is the silver transmittance ``T = 10^-D_Ag`` ``(h, w)``; ``mu`` and ``q`` are
    ``(h, w, K)``; ``count_q`` ``(h, w, K)`` is ``q_c`` where class ``c`` enters the count
    product in that pixel (counted and ``q_c < 1``) and 0 elsewhere; ``intensity_per_um2`` (``nu_c F_c``) and ``radius_um`` (``r~_c``) are the
    Boolean-model inputs, ``(h, w, K)``; ``mean_factor`` is the deterministic part
    ``prod_mean_only exp(-mu_c q_c)`` (demoted pixel-classes included) times the mean factor
    of physical silver without disks, ``(h, w)``; ``var_bool``, ``var_count`` and ``var_eps``
    are ``(h, w)``; ``counted`` is the frame-level ``(K,)`` bool.
    """

    mean_t: np.ndarray
    mu: np.ndarray
    q: np.ndarray
    count_q: np.ndarray
    intensity_per_um2: np.ndarray
    radius_um: np.ndarray
    mean_factor: np.ndarray
    var_bool: np.ndarray
    var_count: np.ndarray
    var_eps: np.ndarray
    counted: np.ndarray

    @property
    def var_pixel(self) -> np.ndarray:
        """Untruncated pixel variance ``Var_count + (T^2 + Var_count) Var(eps)`` (spec 6.6)."""
        return self.var_count + (self.mean_t * self.mean_t + self.var_count) * self.var_eps


@dataclass(frozen=True)
class GrainBlockStats:
    """What the statistical mode records for ``meta["grain"]`` of one window (spec 2.7.1).

    The sums ``sum_pix sum_c nu_c F_c pi r~_c^2`` and ``sum_pix sum_c nu_c F_c (pi r~_c^2)^2``
    are kept as exact float expansions (``disk_area_parts``, ``disk_area2_parts``: floats
    whose exact sum is the exact sum of the per-pixel terms), so merging windows in any order
    or partition gives the same bits (spec 3.6). ``demoted_pixels`` counts the pixels where a
    counted class had ``q_c >= 1`` and was taken as mean-only there.
    """

    mean_only_classes: tuple[int, ...]
    pixels: int
    truncated_pixels: int
    disk_area_parts: tuple[float, ...]
    disk_area2_parts: tuple[float, ...]
    demoted_pixels: int = 0

    def __post_init__(self) -> None:
        for name in ("disk_area_parts", "disk_area2_parts"):
            parts = tuple(float(v) for v in getattr(self, name))
            if not all(math.isfinite(v) for v in parts):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, parts)

    def merge(self, other: "GrainBlockStats") -> "GrainBlockStats":
        """Combine with the stats of another window of the same frame (order-independent)."""
        if other.mean_only_classes != self.mean_only_classes:
            raise ValueError("cannot merge grain stats with different mean-only classes")
        return GrainBlockStats(
            self.mean_only_classes,
            self.pixels + other.pixels,
            self.truncated_pixels + other.truncated_pixels,
            exact_expansion(self.disk_area_parts + other.disk_area_parts),
            exact_expansion(self.disk_area2_parts + other.disk_area2_parts),
            self.demoted_pixels + other.demoted_pixels,
        )

    @property
    def disk_area_sum(self) -> float:
        """``sum_pix sum_c nu_c F_c pi r~_c^2``, correctly rounded."""
        return math.fsum(self.disk_area_parts)

    @property
    def disk_area2_sum(self) -> float:
        """``sum_pix sum_c nu_c F_c (pi r~_c^2)^2``, correctly rounded."""
        return math.fsum(self.disk_area2_parts)

    @property
    def truncated_fraction(self) -> float:
        """Fraction of pixels whose ``T_pix`` was truncated to [0, 1]."""
        return self.truncated_pixels / self.pixels if self.pixels else 0.0

    @property
    def effective_radius_um(self) -> float:
        """``r~_eff = sqrt(sum nu F (pi r~^2)^2 / (pi sum nu F pi r~^2))``, 0 without disks."""
        if not self.disk_area_sum > 0.0:
            return 0.0
        return math.sqrt(self.disk_area2_sum / (math.pi * self.disk_area_sum))


def exact_expansion(values: Iterable[float]) -> tuple[float, ...]:
    """Floats ``r_1, r_2, ...`` whose exact sum is the exact sum of ``values``.

    ``r_1 = fsum(values)``, ``r_k = fsum(values - r_1 - ... - r_{k-1})``, until the remainder
    is 0. The result depends only on the exact sum, not on the order or grouping of the
    values, so expansions of disjoint windows can be concatenated and re-expanded.
    """
    terms = [float(v) for v in values]
    parts: list[float] = []
    for _ in range(_EXPANSION_MAX_TERMS):
        r = math.fsum(itertools.chain(terms, (-p for p in parts)))
        if r == 0.0:
            return tuple(parts)
        parts.append(r)
    raise RuntimeError("exact_expansion did not converge")  # pragma: no cover


def aperture_bias(effective_radius_um: float, aperture_um: float) -> float:
    """``b_ap = 1.15 r~_eff / D_ap`` for a circular aperture of diameter ``aperture_um``."""
    if isinstance(aperture_um, bool) or not np.isfinite(aperture_um) or not aperture_um > 0.0:
        raise ValueError(f"aperture_um must be positive and finite, got {aperture_um!r}")
    if not np.isfinite(effective_radius_um) or effective_radius_um < 0.0:
        raise ValueError(f"effective_radius_um must be finite and >= 0, got {effective_radius_um}")
    return APERTURE_BIAS_COEFFICIENT * float(effective_radius_um) / float(aperture_um)


def grain_meta(stats: GrainBlockStats, aperture_um: float | None = None) -> dict:
    """The ``meta["grain"]`` record of a statistical-mode readout (spec 2.7.1).

    Keys: ``mean_only_classes`` (0-based class indices), ``b_ap`` (``None`` without an
    aperture), ``effective_radius_um``, ``aperture_um``, ``truncated_fraction``,
    ``truncated_pixels``, ``pixels`` and ``demoted_pixels``. Warns when ``b_ap > 0.05``.
    """
    r_eff = stats.effective_radius_um
    b_ap = None if aperture_um is None else aperture_bias(r_eff, aperture_um)
    if b_ap is not None and b_ap > APERTURE_BIAS_WARNING:
        warnings.warn(
            f"statistical grain: aperture bias b_ap = {b_ap:.3f} > {APERTURE_BIAS_WARNING} "
            f"(r_eff = {r_eff:.3f} um, aperture {aperture_um} um); the aperture variance is "
            "biased high",
            stacklevel=2,
        )
    return {
        "mean_only_classes": list(stats.mean_only_classes),
        "b_ap": b_ap,
        "effective_radius_um": r_eff,
        "aperture_um": None if aperture_um is None else float(aperture_um),
        "truncated_fraction": stats.truncated_fraction,
        "truncated_pixels": stats.truncated_pixels,
        "pixels": stats.pixels,
        "demoted_pixels": stats.demoted_pixels,
    }


def _read_only(array: np.ndarray) -> np.ndarray:
    array = np.ascontiguousarray(array, dtype=np.float64)
    array.flags.writeable = False
    return array


def _fields(state: EmulsionState, window: tuple[slice, slice]):
    out = {}
    for name in ("silver", "initiated", "silver_m2", "silver_physical"):
        value = np.asarray(state.get("emulsion", name)[window], dtype=np.float64)
        if not np.all(np.isfinite(value)) or np.any(value < 0.0):
            raise ValueError(f"emulsion/{name} must be finite and non-negative")
        out[name] = value
    k_shape = out["silver"].shape
    for name in ("initiated", "silver_m2"):
        if out[name].shape != k_shape:
            raise ValueError(
                f"emulsion/{name} has shape {out[name].shape}; it must match emulsion/silver "
                f"{k_shape}"
            )
    if out["silver_physical"].shape != k_shape[:2]:
        raise ValueError("emulsion/silver_physical must have shape (H, W)")
    if np.any(out["initiated"] > 1.0 + _FRACTION_TOL):
        raise ValueError("emulsion/initiated must not exceed 1")
    return out


def _boolean_variance(
    lam: np.ndarray, radius: np.ndarray, mean_t: np.ndarray, pixel_area: float
) -> np.ndarray:
    """``boolean_excess_variance`` pixel by pixel, in chunks; one call for a uniform window."""
    h, w, k = lam.shape
    n = h * w
    lam2 = lam.reshape(n, k)
    rad2 = radius.reshape(n, k)
    t2 = mean_t.reshape(n)
    uniform = np.all(lam2 == lam2[0]) and np.all(rad2 == rad2[0]) and np.all(t2 == t2[0])
    if uniform:
        one = boolean_excess_variance(lam2[:1], rad2[:1], t2[:1], pixel_area)
        return np.full((h, w), one[0], dtype=np.float64)
    out = np.empty(n, dtype=np.float64)
    for start in range(0, n, _VAR_CHUNK_PX):
        stop = min(start + _VAR_CHUNK_PX, n)
        out[start:stop] = boolean_excess_variance(
            lam2[start:stop], rad2[start:stop], t2[start:stop], pixel_area
        )
    return out.reshape(h, w)


def grain_statistics(state: EmulsionState, rows: slice, cols: slice) -> GrainStatistics:
    """Per-pixel statistics of the statistical grain mode over a window (spec 2.7.1).

    Needs ``emulsion/silver``, ``initiated``, ``silver_m2`` and ``silver_physical``.
    Raises ``ValueError`` on negative or non-finite fields, on ``initiated > 1``, on developed
    silver above the initiated grains (``X_c > F_c``, including silver with ``F_c = 0``) and on
    ``Y_c > X_c``, each beyond a relative rounding tolerance of 1e-9.
    """
    window = resolve_window(state, rows, cols)
    const = imaging_constants(state)
    f = _fields(state, window)
    silver, initiated, m2, physical = (
        f["silver"],
        f["initiated"],
        f["silver_m2"],
        f["silver_physical"],
    )
    nu = np.asarray(const.classes.number_per_um2)
    a = np.asarray(const.developed_areas_um2)
    pitch = state.geometry.pixel_pitch_um
    pixel_area = pitch * pitch

    developed = silver / const.class_silver_g_m2  # X_c
    _check_not_above("emulsion/silver (X)", developed, "emulsion/initiated (F)", initiated)
    _check_not_above("emulsion/silver_m2 (Y)", m2, "emulsion/silver (X)", developed)
    initiated_pos = initiated > 0.0
    f_safe = np.where(initiated_pos, initiated, 1.0)
    x_mean = np.clip(np.where(initiated_pos, developed / f_safe, 0.0), 0.0, 1.0)
    x2_mean = np.clip(np.where(initiated_pos, m2 / f_safe, 0.0), 0.0, 1.0)

    chemical = chemical_silver_density(silver, const.covering_power)
    physical_density = const.physical_covering_ratio * const.covering_power_ref * physical
    has_disks = chemical > 0.0
    phi = np.where(has_disks, physical_density / np.where(has_disks, chemical, 1.0), 0.0)
    physical_mean = np.where(has_disks, 0.0, physical_density)
    scale = 1.0 + phi

    mu = nu * pixel_area * initiated
    q = a * scale[..., None] * x_mean / pixel_area
    lam = nu * initiated
    radius = np.sqrt(a * scale[..., None] * np.sqrt(x2_mean) / math.pi)
    counted = counted_classes(a, pitch)
    # Pointwise: a counted class with q_c >= 1 is mean-only in that pixel.
    in_count = counted & (q < 1.0)
    count_q = np.where(in_count, q, 0.0)

    total = np.zeros(chemical.shape, dtype=np.float64)
    mean_only = LN10 * physical_mean
    for c in range(mu.shape[-1]):
        mq = mu[..., c] * q[..., c]
        total = total + mq
        if counted[c]:
            mean_only = mean_only + np.where(in_count[..., c], 0.0, mq)
        else:
            mean_only = mean_only + mq
    mean_t = np.exp(-(total + LN10 * physical_mean))
    mean_factor = np.exp(-mean_only)

    var_bool = _boolean_variance(lam, radius, mean_t, pixel_area)
    var_count = _count_variance(mean_t, mu, count_q, counted)
    var_eps = epsilon_variance(var_bool, mean_t, mu, count_q, counted)
    counted_ro = np.array(counted, dtype=bool)
    counted_ro.flags.writeable = False
    return GrainStatistics(
        mean_t=_read_only(mean_t),
        mu=_read_only(mu),
        q=_read_only(q),
        count_q=_read_only(count_q),
        intensity_per_um2=_read_only(lam),
        radius_um=_read_only(radius),
        mean_factor=_read_only(mean_factor),
        var_bool=_read_only(var_bool),
        var_count=_read_only(var_count),
        var_eps=_read_only(var_eps),
        counted=counted_ro,
    )


def _check_not_above(name_lo: str, lo: np.ndarray, name_hi: str, hi: np.ndarray) -> None:
    """Raise unless ``lo <= hi`` within the relative tolerance ``_FRACTION_TOL``."""
    bad = lo > hi + _FRACTION_TOL * np.maximum(hi, lo)
    if np.any(bad):
        n = int(np.count_nonzero(bad))
        raise ValueError(
            f"{name_lo} exceeds {name_hi} in {n} pixel-class entries; the developed fields "
            "must satisfy Y <= X <= F (spec 2.0.3)"
        )


def _survival(q: np.ndarray, n: np.ndarray) -> np.ndarray:
    """``(1 - q)^n`` as ``exp(n log1p(-q))``, with ``1 - q`` clipped to [0, 1].

    Not ``np.power``: it squares a size-1 operand on a fast path, so a 1-pixel window would
    round differently from a larger one.
    """
    q = np.clip(q, 0.0, 1.0)
    opaque = q >= 1.0
    log_survive = np.log1p(-np.where(opaque, 0.0, q))
    out = np.exp(n * log_survive)
    return np.where(opaque, np.where(n > 0.0, 0.0, 1.0), out)


def statistical_transmittance_with_stats(
    state: EmulsionState, rows: slice, cols: slice
) -> tuple[np.ndarray, GrainBlockStats]:
    """Sampled silver transmittance ``T_pix`` of a window and its ``meta["grain"]`` sums.

    See the module docs; returns a read-only ``(h, w)`` array in [0, 1] and the
    :class:`GrainBlockStats` of the window (feed it to :func:`grain_meta`).
    """
    window = resolve_window(state, rows, cols)
    stats = grain_statistics(state, *window)
    geometry = state.geometry
    r0 = geometry.origin_px[0] + window[0].start
    c0 = geometry.origin_px[1] + window[1].start
    h = window[0].stop - window[0].start
    w = window[1].stop - window[1].start
    frame_width = geometry.frame_shape_hw[1]
    abs_rows, abs_cols = (r0, r0 + h), (c0, c0 + w)

    t = np.array(stats.mean_factor, dtype=np.float64)
    for c in np.flatnonzero(stats.counted):
        u = pixel_uniform(state.seed, class_stream(c), abs_rows, abs_cols, frame_width)
        n = poisson_inverse(u, stats.mu[..., c]).astype(np.float64)
        t = t * _survival(stats.count_q[..., c], n)
    u = pixel_uniform(state.seed, OVERLAP_STREAM, abs_rows, abs_cols, frame_width)
    z = ndtri(np.clip(u, _U_EPS, 1.0 - _U_EPS))
    t = t * (1.0 + np.sqrt(stats.var_eps) * z)
    truncated = int(np.count_nonzero((t < 0.0) | (t > 1.0)))
    t = np.clip(t, 0.0, 1.0)

    area = math.pi * stats.radius_um * stats.radius_um
    weighted = stats.intensity_per_um2 * area
    # Per-pixel class sums as explicit adds, then exact (order-free) sums over the pixels.
    disk = np.zeros((h, w), dtype=np.float64)
    disk2 = np.zeros((h, w), dtype=np.float64)
    for c in range(weighted.shape[-1]):
        disk = disk + weighted[..., c]
        disk2 = disk2 + weighted[..., c] * area[..., c]
    demoted = stats.counted & (stats.count_q != stats.q)
    block = GrainBlockStats(
        mean_only_classes=tuple(int(c) for c in np.flatnonzero(~stats.counted)),
        pixels=h * w,
        truncated_pixels=truncated,
        disk_area_parts=exact_expansion(disk.ravel().tolist()),
        disk_area2_parts=exact_expansion(disk2.ravel().tolist()),
        demoted_pixels=int(np.count_nonzero(np.any(demoted, axis=-1))),
    )
    return _read_only(t), block


def statistical_transmittance(state: EmulsionState, rows: slice, cols: slice) -> np.ndarray:
    """Sampled silver transmittance ``T_pix`` of a window, ``(h, w)``, keyed by ``state.seed``.

    The value of a pixel depends only on the state's fields there, the seed and the pixel's
    absolute position in the frame, so tiled and untiled results are identical bit for bit.
    """
    return statistical_transmittance_with_stats(state, rows, cols)[0]
