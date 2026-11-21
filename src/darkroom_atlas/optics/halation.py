"""Halation: diffuse light reflected back to the emulsion by the base-air interface (§2.1.5).

Lambertian light entering the base at one point returns at ``r = 2 t_b tan(theta)``. Per unit
incident flux the returned irradiance is

    H(r) = R(theta) T_g(theta) cos^4(theta) / (4 pi t_b^2)            [um^-2]
    T_g(theta) = 10^(-2 (D_g + D_back) / cos(theta))
    f_H = integral_0^{pi/2} 2 sin(theta) cos(theta) R(theta) T_g(theta) d(theta)

with ``R`` the unpolarised Fresnel reflectance from the base (index ``n_b``) into air, 1 at
and beyond the critical angle ``theta_c = arcsin(1/n_b)``. Only the first return is kept.
"""

import math
from dataclasses import dataclass
from numbers import Real

import numpy as np
from scipy import integrate, optimize, special

from darkroom_atlas.optics.stack import BaseOptics

__all__ = [
    "HalationKernel",
    "halation_fraction",
    "halation_kernel",
    "halation_otf",
]

# Pixel supersampling (spec §2.1.5): 8 x 8 in general, 32 x 32 across r = r_c.
_SUPERSAMPLE = 8
_SUPERSAMPLE_CRITICAL = 32
_QUAD = {"epsabs": 1e-13, "epsrel": 1e-12, "limit": 500}


@dataclass(frozen=True, eq=False)
class HalationKernel:
    """Pixel halation kernel (spec §2.1.5).

    ``kernel`` is an odd square holding the returned fraction per pixel for unit flux into
    the base; ``fraction`` is the untruncated ``f_H``; ``truncated = fraction - kernel.sum()``;
    ``radius_um`` is the truncation radius ``R_H``; ``critical_radius_um`` is ``r_c``.
    """

    kernel: np.ndarray
    fraction: float
    truncated: float
    radius_um: float
    critical_radius_um: float


def _fresnel(theta: np.ndarray, n: float) -> np.ndarray:
    """Unpolarised Fresnel reflectance from index ``n`` into air; 1 beyond ``theta_c``."""
    theta = np.asarray(theta, dtype=np.float64)
    sin_t = n * np.sin(theta)
    cos_i = np.cos(theta)
    cos_t = np.sqrt(np.clip(1.0 - sin_t**2, 0.0, None))
    with np.errstate(invalid="ignore", divide="ignore"):
        r_s = ((n * cos_i - cos_t) / (n * cos_i + cos_t)) ** 2
        r_p = ((cos_i - n * cos_t) / (cos_i + n * cos_t)) ** 2
    return np.where(sin_t >= 1.0, 1.0, 0.5 * (r_s + r_p))


def _transmission(cos_theta: np.ndarray, base: BaseOptics) -> np.ndarray:
    """Two oblique passes through the grey base and backing, ``T_g`` (spec §2.1.5)."""
    density = base.grey_density + base.backing_density
    if density == 0.0:
        return np.ones_like(cos_theta)
    with np.errstate(divide="ignore", over="ignore"):
        exponent = -2.0 * density / np.asarray(cos_theta)
    return np.power(10.0, exponent)


def _angular(theta: float, base: BaseOptics) -> float:
    """Return density per unit angle, ``2 sin cos R T_g``."""
    c = math.cos(theta)
    return float(
        2.0 * math.sin(theta) * c * _fresnel(theta, base.refractive_index) * _transmission(c, base)
    )


def _fraction_parts(base: BaseOptics) -> tuple[float, float]:
    """Return ``(Fresnel part for theta < theta_c, total-reflection part)`` of ``f_H``."""
    theta_c = base.critical_angle_rad
    fresnel = integrate.quad(_angular, 0.0, theta_c, args=(base,), **_QUAD)[0]
    total = integrate.quad(_angular, theta_c, 0.5 * math.pi, args=(base,), **_QUAD)[0]
    return fresnel, total


def _check_base(base: BaseOptics) -> None:
    if not isinstance(base, BaseOptics):
        raise TypeError(f"base must be a BaseOptics, got {type(base).__name__}")


def halation_fraction(base: BaseOptics) -> float:
    """Return fraction ``f_H`` of the flux entering the base (spec §2.1.5)."""
    _check_base(base)
    return float(sum(_fraction_parts(base)))


def _cumulative(theta: float, base: BaseOptics) -> float:
    """Return fraction returned within angle ``theta``."""
    theta_c = base.critical_angle_rad
    value = integrate.quad(_angular, 0.0, min(theta, theta_c), args=(base,), **_QUAD)[0]
    if theta > theta_c:
        value += integrate.quad(_angular, theta_c, theta, args=(base,), **_QUAD)[0]
    return value


def _energy_radius(base: BaseOptics, share: float) -> float:
    """Radius (um) inside which ``share`` of ``f_H`` returns."""
    if not 0.0 < share < 1.0:
        raise ValueError(f"share must be in (0, 1), got {share}")
    target = share * halation_fraction(base)
    theta = optimize.brentq(
        lambda t: _cumulative(t, base) - target, 0.0, 0.5 * math.pi, xtol=1e-15, rtol=1e-14
    )
    return 2.0 * base.thickness_um * math.tan(theta)


def _irradiance(r: np.ndarray, base: BaseOptics) -> np.ndarray:
    """Return irradiance ``H(r)`` per unit flux (um^-2)."""
    t = base.thickness_um
    two_t = 2.0 * t
    cos_theta = two_t / np.hypot(r, two_t)
    theta = np.arctan2(r, two_t)
    return (
        _fresnel(theta, base.refractive_index)
        * _transmission(cos_theta, base)
        * cos_theta**4
        / (4.0 * math.pi * t**2)
    )


def _positive(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number, got {value!r}")
    if value <= 0.0:
        raise ValueError(f"{name} must be positive, got {value}")
    return float(value)


def halation_kernel(
    pitch_um: float, base: BaseOptics, *, energy_tol: float = 1e-3, max_radius_um: float = 3000.0
) -> HalationKernel:
    """Pixel kernel of ``H`` cut at ``R_H`` (spec §2.1.5).

    ``R_H`` is the radius holding ``(1 - energy_tol) f_H``, capped at ``max_radius_um``. Each
    pixel holds the integral of ``H`` over the pixel from 8 x 8 midpoint samples (32 x 32 for
    pixels crossing ``r_c``); samples beyond ``R_H`` are dropped. The kernel is not
    renormalised.
    """
    _check_base(base)
    pitch = _positive("pitch_um", pitch_um)
    max_radius = _positive("max_radius_um", max_radius_um)
    tol = _positive("energy_tol", energy_tol)
    if tol >= 1.0:
        raise ValueError(f"energy_tol must be below 1, got {energy_tol}")
    fraction = halation_fraction(base)
    radius = min(_energy_radius(base, 1.0 - tol), max_radius)
    r_c = base.critical_radius_um
    half = math.floor(radius / pitch + 0.5)
    quadrant = _quadrant(pitch, half, radius, base, _SUPERSAMPLE)
    # Pixels whose nearest and farthest points bracket r_c get the finer sampling.
    index = np.arange(half + 1)
    near_1d = np.maximum(index - 0.5, 0.0) * pitch
    far_1d = (index + 0.5) * pitch
    near = np.hypot(near_1d[:, None], near_1d[None, :])
    far = np.hypot(far_1d[:, None], far_1d[None, :])
    for i, j in zip(*np.nonzero((near < r_c) & (far > r_c)), strict=True):
        quadrant[i, j] = _pixel(pitch, int(i), int(j), radius, base, _SUPERSAMPLE_CRITICAL)
    kernel = _mirror(quadrant)
    kernel.flags.writeable = False
    return HalationKernel(
        kernel=kernel,
        fraction=fraction,
        truncated=float(fraction - kernel.sum()),
        radius_um=float(radius),
        critical_radius_um=float(r_c),
    )


def _offsets(pitch: float, ss: int) -> np.ndarray:
    return ((np.arange(ss) + 0.5) / ss - 0.5) * pitch


def _quadrant(pitch: float, half: int, radius: float, base: BaseOptics, ss: int) -> np.ndarray:
    """Pixel integrals for indices ``0..half`` on both axes, one pixel row at a time."""
    offsets = _offsets(pitch, ss)
    cols = (np.arange(half + 1)[:, None] * pitch + offsets[None, :]).ravel()
    out = np.empty((half + 1, half + 1))
    for i in range(half + 1):
        r = np.hypot((i * pitch + offsets)[:, None], cols[None, :])
        values = np.where(r <= radius, _irradiance(r, base), 0.0)
        out[i] = values.reshape(ss, half + 1, ss).mean(axis=(0, 2)) * pitch**2
    return out


def _pixel(pitch: float, i: int, j: int, radius: float, base: BaseOptics, ss: int) -> float:
    offsets = _offsets(pitch, ss)
    r = np.hypot((i * pitch + offsets)[:, None], (j * pitch + offsets)[None, :])
    values = np.where(r <= radius, _irradiance(r, base), 0.0)
    return float(values.mean() * pitch**2)


def _mirror(quadrant: np.ndarray) -> np.ndarray:
    """Full odd square kernel from its quadrant ``[0..h, 0..h]`` (index 0 is the centre)."""
    quadrant = 0.5 * (quadrant + quadrant.T)  # exact x <-> y symmetry despite rounding
    rows = np.concatenate([quadrant[:0:-1], quadrant], axis=0)
    return np.concatenate([rows[:, :0:-1], rows], axis=1)


def halation_otf(nu_per_um: np.ndarray, base: BaseOptics) -> np.ndarray:
    """Transfer function of the untruncated ``H`` (its Hankel transform); ``nu`` in cycles/um.

    ``H^(nu) = integral 2 sin cos R T_g J0(4 pi nu t_b tan(theta)) d(theta)``, so
    ``H^(0) = f_H``.

    For a clear base (``D_g + D_back == 0``) the total-reflection part uses the closed form
    ``integral_0^{pi/2} 2 sin cos J0(2 t_b k tan) = 2 t_b k K1(2 t_b k)``, ``k = 2 pi nu``,
    minus its ``theta < theta_c`` share. With an absorbing base the integral is done by
    quadrature; for very small densities (well below 0.1) the integrand oscillates without
    decaying out to large radii and ``quad`` may warn about accuracy.
    """
    _check_base(base)
    nu = np.asarray(nu_per_um, dtype=np.float64)
    if not np.all(np.isfinite(nu)):
        raise ValueError("nu_per_um must be finite")
    theta_c = base.critical_angle_rad
    scale = 4.0 * math.pi * base.thickness_um
    clear = base.grey_density + base.backing_density == 0.0

    def integrand(theta: float, frequency: float) -> float:
        return _angular(theta, base) * float(special.j0(scale * frequency * math.tan(theta)))

    def clear_tir(theta: float, frequency: float) -> float:
        return math.sin(2.0 * theta) * float(special.j0(scale * frequency * math.tan(theta)))

    out = np.empty(nu.shape)
    for idx, frequency in np.ndenumerate(nu):
        f = abs(float(frequency))
        value = integrate.quad(integrand, 0.0, theta_c, args=(f,), **_QUAD)[0]
        if clear:
            x = scale * f  # 2 t_b k
            full = 1.0 if x == 0.0 else float(x * special.k1(x))
            value += full - integrate.quad(clear_tir, 0.0, theta_c, args=(f,), **_QUAD)[0]
        else:
            value += integrate.quad(integrand, theta_c, 0.5 * math.pi, args=(f,), **_QUAD)[0]
        out[idx] = value
    return out
