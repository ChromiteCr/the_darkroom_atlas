"""In-layer scatter point spread function (spec §2.1.6).

The PSF that belongs to Frieser's exponential line spread function:

    P(r) = (1 - w_P) delta(r) + w_P K0(r / a_P) / (2 pi a_P^2)        [um^-2]
    P^(nu) = 1 - w_P + w_P / (1 + (2 pi a_P nu)^2)                    (nu in cycles/um)

The mass of the K0 term inside radius R is ``1 - (R/a) K1(R/a)``.
"""

import math
from functools import lru_cache
from numbers import Integral, Real

import numpy as np
from scipy import optimize, special

__all__ = ["scatter_half_width", "scatter_kernel", "scatter_otf", "tail_radius_ratio"]


def _check(weight: float, length_um: float) -> tuple[float, float]:
    for name, value in (("weight", weight), ("length_um", length_um)):
        if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
            raise ValueError(f"{name} must be a finite number, got {value!r}")
    if not 0.0 <= weight <= 1.0:
        raise ValueError(f"weight must be in [0, 1], got {weight}")
    if length_um <= 0.0:
        raise ValueError(f"length_um must be positive, got {length_um}")
    return float(weight), float(length_um)


def scatter_otf(nu_per_um: np.ndarray, weight: float, length_um: float) -> np.ndarray:
    """Transfer function ``P^(nu)`` of the untruncated PSF (spec §2.1.6); ``nu`` in cycles/um."""
    weight, length_um = _check(weight, length_um)
    nu = np.asarray(nu_per_um, dtype=np.float64)
    return 1.0 - weight + weight / (1.0 + (2.0 * math.pi * length_um * nu) ** 2)


@lru_cache(maxsize=32)
def tail_radius_ratio(tail: float) -> float:
    """``x`` with ``x K1(x) = tail``: the truncation radius in units of ``a_P`` (spec §12.7)."""
    if not 0.0 < tail < 1.0:
        raise ValueError(f"tail must be in (0, 1), got {tail}")
    # x K1(x) falls monotonically from 1 at x = 0; it is below 1e-300 well before x = 700.
    return float(optimize.brentq(lambda x: x * special.k1(x) - tail, 1e-12, 700.0, xtol=1e-14))


def scatter_kernel(
    pitch_um: float, weight: float, length_um: float, *, tail: float = 1e-4, supersample: int = 16
) -> np.ndarray:
    """Pixel kernel of the scatter PSF: odd square, sums to 1 (spec §2.1.6).

    Each pixel holds the PSF mass over the pixel, from ``supersample`` x ``supersample``
    midpoint samples; the centre pixel uses the equal-area disc ``1 - (rho/a) K1(rho/a)``,
    ``rho = pitch/sqrt(pi)``. The K0 term is cut at ``R_P`` with ``(R_P/a) K1(R_P/a) = tail``
    and renormalised to ``w_P``; the delta part ``1 - w_P`` sits on the centre pixel.
    """
    weight, length_um = _check(weight, length_um)
    if isinstance(pitch_um, bool) or not isinstance(pitch_um, Real) or not pitch_um > 0.0:
        raise ValueError(f"pitch_um must be a positive number, got {pitch_um!r}")
    if not math.isfinite(pitch_um):
        raise ValueError(f"pitch_um must be finite, got {pitch_um!r}")
    if isinstance(supersample, bool) or not isinstance(supersample, Integral) or supersample < 1:
        raise ValueError(f"supersample must be a positive integer, got {supersample!r}")
    if isinstance(tail, bool) or not isinstance(tail, Real):
        raise TypeError(f"tail must be a number, got {tail!r}")
    pitch = float(pitch_um)
    radius = tail_radius_ratio(float(tail)) * length_um  # R_P, spec §2.1.6 / §12.7
    half = scatter_half_width(pitch, length_um, tail=float(tail))
    if weight == 0.0:
        kernel = np.zeros((2 * half + 1, 2 * half + 1))
        kernel[half, half] = 1.0
        kernel.flags.writeable = False
        return kernel
    ss = int(supersample)
    offsets = ((np.arange(ss) + 0.5) / ss - 0.5) * pitch
    cols = (np.arange(half + 1)[:, None] * pitch + offsets[None, :]).ravel()
    # One quadrant (pixel indices 0..half on both axes), one pixel row at a time, then mirror.
    quadrant = np.empty((half + 1, half + 1))
    for i in range(half + 1):
        r = np.hypot((i * pitch + offsets)[:, None], cols[None, :])
        inside = (r <= radius) & (r > 0.0)
        density = np.zeros_like(r)
        density[inside] = special.k0(r[inside] / length_um) / (2.0 * math.pi * length_um**2)
        quadrant[i] = density.reshape(ss, half + 1, ss).mean(axis=(0, 2)) * pitch**2
    rho = min(pitch / math.sqrt(math.pi), radius)
    quadrant[0, 0] = 1.0 - (rho / length_um) * special.k1(rho / length_um)
    scattered = _mirror(quadrant)
    scattered /= scattered.sum()
    kernel = weight * scattered
    kernel[half, half] += 1.0 - weight
    kernel.flags.writeable = False
    return kernel


def _mirror(quadrant: np.ndarray) -> np.ndarray:
    """Full odd square kernel from its quadrant ``[0..h, 0..h]`` (index 0 is the centre)."""
    quadrant = 0.5 * (quadrant + quadrant.T)  # exact x <-> y symmetry despite rounding
    rows = np.concatenate([quadrant[:0:-1], quadrant], axis=0)
    return np.concatenate([rows[:, :0:-1], rows], axis=1)


def scatter_half_width(pitch_um: float, length_um: float, *, tail: float = 1e-4) -> int:
    """Half width in pixels of :func:`scatter_kernel`: pixels that reach inside ``R_P``."""
    radius = tail_radius_ratio(float(tail)) * float(length_um)
    return math.floor(radius / float(pitch_um) + 0.5)
