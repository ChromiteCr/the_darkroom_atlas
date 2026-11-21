"""Poisson inverse-CDF sampling from given uniforms (spec D1 section 2.7.1).

The sample depends only on ``(u, mean)`` element by element, so a tile of uniforms drawn
from :func:`darkroom_atlas.rng.pixel_uniform` gives the same counts whatever the tiling.
"""

import numpy as np
from scipy.special import ndtri

__all__ = ["poisson_inverse"]

# Uniforms are clipped to [2^-53, 1 - 2^-53] before ndtri so z stays finite (|z| < 8.3).
_U_EPS = 2.0**-53


def _upper_bound(mean: np.ndarray) -> np.ndarray:
    """Largest count returned for ``mean``: floor(mu + 10 sqrt(mu) + 10) (spec 2.7.1)."""
    return np.floor(mean + 10.0 * np.sqrt(mean) + 10.0)


def _exact(u: np.ndarray, mean: np.ndarray) -> np.ndarray:
    """Smallest k with CDF(k) >= u, by summing the pmf recursively up to the bound."""
    cap = _upper_bound(mean)
    k_max = int(cap.max()) if cap.size else 0
    pmf = np.exp(-mean)
    cdf = pmf.copy()
    k = np.zeros(mean.shape, dtype=np.int64)
    for j in range(k_max):
        k += u > cdf
        pmf = pmf * mean / (j + 1)
        cdf = cdf + pmf
    return np.minimum(k, cap.astype(np.int64))


def _cornish_fisher(u: np.ndarray, mean: np.ndarray) -> np.ndarray:
    """Cornish-Fisher quantile k = floor(mu + sqrt(mu) z + (z^2-1)/6 + (z^3-7z)/(72 sqrt(mu)) + 1/2).

    Spec 2.7.1; the result is clipped to [0, floor(mu + 10 sqrt(mu) + 10)].
    """
    z = ndtri(np.clip(u, _U_EPS, 1.0 - _U_EPS))
    s = np.sqrt(mean)
    k = np.floor(mean + s * z + (z * z - 1.0) / 6.0 + (z**3 - 7.0 * z) / (72.0 * s) + 0.5)
    return np.clip(k, 0.0, _upper_bound(mean)).astype(np.int64)


def poisson_inverse(u: np.ndarray, mean: np.ndarray, *, exact_below: float = 30.0) -> np.ndarray:
    """Poisson counts with means ``mean`` from uniforms ``u`` in ``[0, 1]``, as ``int64``.

    Where ``mean <= exact_below`` the exact inverse CDF is used (the smallest k with
    ``CDF(k) >= u``, capped at ``floor(mean + 10 sqrt(mean) + 10)``); above it the
    Cornish-Fisher expansion of spec D1 section 2.7.1. ``u`` and ``mean`` broadcast against
    each other; the result has the broadcast shape.

    Raises ``ValueError`` if ``u`` is outside ``[0, 1]`` or not finite, or if ``mean`` is
    negative or not finite.
    """
    u_arr = np.asarray(u, dtype=np.float64)
    mean_arr = np.asarray(mean, dtype=np.float64)
    try:
        u_arr, mean_arr = np.broadcast_arrays(u_arr, mean_arr)
    except ValueError as exc:
        raise ValueError(
            f"u of shape {u_arr.shape} and mean of shape {mean_arr.shape} do not broadcast"
        ) from exc
    if not np.all(np.isfinite(u_arr)) or np.any((u_arr < 0.0) | (u_arr > 1.0)):
        raise ValueError("u must be finite and lie in [0, 1]")
    if not np.all(np.isfinite(mean_arr)) or np.any(mean_arr < 0.0):
        raise ValueError("mean must be finite and non-negative")
    if isinstance(exact_below, bool) or not np.isfinite(exact_below) or exact_below < 0.0:
        raise ValueError(f"exact_below must be a non-negative finite number, got {exact_below!r}")

    out = np.zeros(u_arr.shape, dtype=np.int64)
    small = mean_arr <= exact_below
    if np.any(small):
        out[small] = _exact(u_arr[small], mean_arr[small])
    large = ~small
    if np.any(large):
        out[large] = _cornish_fisher(u_arr[large], mean_arr[large])
    return out
