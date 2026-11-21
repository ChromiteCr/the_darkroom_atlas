"""Boolean-model overlap of developed grains as disks (spec D1 section 2.7.1, OPT 5.4).

For a planar Poisson field of opaque disks the transmittance covariance at separation tau is
``T^2 (exp(gamma_mix(tau)) - 1)`` with ``gamma_mix = sum_c lambda_c gamma(tau; r_c)`` and
``gamma`` the intersection area of two disks. Integrated over the plane and divided by the
pixel area it gives the pixel variance ``Var_bool``.
"""

import math

import numpy as np
from numpy.polynomial.legendre import leggauss

__all__ = ["boolean_excess_variance", "disk_intersection_area", "kappa"]

# ln(10): net density D = lambda * pi r^2 / ln 10 for constant-radius disks (Nutting, spec 2.5).
_LN10 = math.log(10.0)


def _check_nodes(nodes: int) -> int:
    if isinstance(nodes, bool) or not isinstance(nodes, int | np.integer) or nodes < 1:
        raise ValueError(f"nodes must be a positive integer, got {nodes!r}")
    return int(nodes)


def disk_intersection_area(tau_um: np.ndarray, radius_um: np.ndarray) -> np.ndarray:
    """Intersection area (um^2) of two disks of radius ``radius_um`` whose centres are ``tau_um`` apart.

    ``gamma = 2 r^2 arccos(tau / 2r) - (tau / 2) sqrt(4 r^2 - tau^2)`` for ``tau < 2r``, else 0
    (OPT 5.4). Arguments broadcast. Raises ``ValueError`` on negative or non-finite input.
    """
    tau = np.asarray(tau_um, dtype=np.float64)
    r = np.asarray(radius_um, dtype=np.float64)
    if not np.all(np.isfinite(tau)) or np.any(tau < 0.0):
        raise ValueError("tau_um must be finite and non-negative")
    if not np.all(np.isfinite(r)) or np.any(r < 0.0):
        raise ValueError("radius_um must be finite and non-negative")
    tau, r = np.broadcast_arrays(tau, r)
    d = 2.0 * r
    inside = tau < d
    t = np.where(inside, tau, 0.0)
    safe_d = np.where(inside, d, 1.0)
    area = 2.0 * r * r * np.arccos(t / safe_d) - 0.5 * t * np.sqrt(np.maximum(d * d - t * t, 0.0))
    return np.where(inside, area, 0.0)


def boolean_excess_variance(
    intensity_per_um2: np.ndarray,
    radius_um: np.ndarray,
    mean_t: np.ndarray,
    pixel_area_um2: float,
    *,
    nodes: int = 48,
) -> np.ndarray:
    """Boolean-model pixel transmittance variance ``Var_bool`` (spec D1 section 2.7.1).

    ``Var_bool = (T^2 / A) * integral_0^{2 r_max} (exp(gamma_mix(tau)) - 1) 2 pi tau dtau``
    with ``gamma_mix(tau) = sum_c intensity_c * gamma(tau; r_c)``, evaluated with ``nodes``
    point Gauss-Legendre on ``[0, 2 r_max]`` where ``r_max`` is the largest radius of each
    pixel.

    ``intensity_per_um2`` (``nu_c F_c``, disks per um^2) and ``radius_um`` (``r~_c``) carry
    the size classes on their last axis and broadcast to ``(..., K)``; ``mean_t`` (``T``)
    broadcasts to ``(...)``. Returns an array of shape ``(...)``.
    """
    nodes = _check_nodes(nodes)
    lam = np.asarray(intensity_per_um2, dtype=np.float64)
    r = np.asarray(radius_um, dtype=np.float64)
    t_mean = np.asarray(mean_t, dtype=np.float64)
    if lam.ndim == 0 or r.ndim == 0:
        raise ValueError("intensity_per_um2 and radius_um need a trailing size-class axis")
    if not np.all(np.isfinite(lam)) or np.any(lam < 0.0):
        raise ValueError("intensity_per_um2 must be finite and non-negative")
    if not np.all(np.isfinite(r)) or np.any(r < 0.0):
        raise ValueError("radius_um must be finite and non-negative")
    if not np.all(np.isfinite(t_mean)) or np.any((t_mean < 0.0) | (t_mean > 1.0)):
        raise ValueError("mean_t must be finite and lie in [0, 1]")
    if isinstance(pixel_area_um2, bool) or not np.isfinite(pixel_area_um2) or pixel_area_um2 <= 0.0:
        raise ValueError(f"pixel_area_um2 must be positive and finite, got {pixel_area_um2!r}")
    try:
        lam, r = np.broadcast_arrays(lam, r)
        lead = np.broadcast_shapes(lam.shape[:-1], t_mean.shape)
    except ValueError as exc:
        raise ValueError(
            f"shapes do not broadcast: intensity {np.shape(intensity_per_um2)}, "
            f"radius {np.shape(radius_um)}, mean_t {t_mean.shape}"
        ) from exc
    lam = np.broadcast_to(lam, (*lead, lam.shape[-1]))
    r = np.broadcast_to(r, (*lead, r.shape[-1]))

    x, w = leggauss(nodes)
    upper = 2.0 * r.max(axis=-1)  # (...)
    half = 0.5 * upper[..., None]
    tau = half * (1.0 + x)  # (..., nodes)
    gamma_mix = np.zeros(tau.shape, dtype=np.float64)
    for c in range(lam.shape[-1]):
        gamma_mix += lam[..., c, None] * disk_intersection_area(tau, r[..., c, None])
    integral = np.sum(half * w * np.expm1(gamma_mix) * 2.0 * np.pi * tau, axis=-1)
    return np.broadcast_to(t_mean, lead) ** 2 / float(pixel_area_um2) * integral


def kappa(net_density: float, *, nodes: int = 48) -> float:
    """Boolean overlap factor ``Var_bool / Var_Siedentopf`` for constant-radius disks (OPT 5.4).

    With unit radius and ``lambda pi = D ln 10``, ``Var_Siedentopf`` is the linearised
    integral ``lambda (pi r^2)^2``; the ratio is independent of radius and pixel area.
    """
    nodes = _check_nodes(nodes)
    if isinstance(net_density, bool) or not np.isfinite(net_density) or net_density <= 0.0:
        raise ValueError(f"net_density must be positive and finite, got {net_density!r}")
    lam = float(net_density) * _LN10 / math.pi
    var_bool = boolean_excess_variance(np.array([lam]), np.array([1.0]), 1.0, 1.0, nodes=nodes)
    return float(var_bool) / (lam * math.pi**2)
