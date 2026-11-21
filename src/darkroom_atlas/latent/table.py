"""Per-frame latent-image lookup table (spec D1 2.2.5).

For a fixed exposure time each size class gets ``p(n | m)`` on ``log10 m`` nodes spaced
evenly over ``[-6, 4]`` (201 nodes by default). Lookup interpolates each entry with a cubic
spline in ``log10 m`` (spec D1 3.1; linear interpolation misses the 1e-3 bound of the
interpolation test by about 2x on this grid), clamps to the first node below the range and
to the last node above it, clips the tiny spline undershoots at 0 and renormalises.
"""

import math
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline

from darkroom_atlas.emulsion.classes import SizeClasses
from darkroom_atlas.latent.chain import (
    ChainParams,
    _check_quadrature,
    _check_time,
    _distributions,
    efficiency,
)

__all__ = ["LatentTable", "latent_table"]


@dataclass(frozen=True, eq=False)
class LatentTable:
    """Cluster-size distributions per size class on a ``log10 m`` grid.

    ``log10_m`` (N_m,) is evenly spaced; ``probabilities`` (K, N_m, cap+1) is read-only,
    entry ``n = 1`` is 0 and each row sums to 1. ``size_class`` is the 0-based class index
    of the :class:`SizeClasses` the table was built from.
    """

    log10_m: np.ndarray  # (N_m,)
    probabilities: np.ndarray  # (K, N_m, cap+1), read-only
    exposure_time_s: float

    def __post_init__(self) -> None:
        grid = np.asarray(self.log10_m, dtype=np.float64)
        prob = np.asarray(self.probabilities, dtype=np.float64)
        if grid.ndim != 1 or grid.shape[0] < 2:
            raise ValueError("log10_m must be 1-D with at least 2 nodes")
        step = np.diff(grid)
        if not (np.all(step > 0.0) and np.allclose(step, step[0], rtol=1e-9, atol=0.0)):
            raise ValueError("log10_m must be evenly spaced and increasing")
        if prob.ndim != 3 or prob.shape[1] != grid.shape[0] or prob.shape[2] < 3:
            raise ValueError(
                f"probabilities must have shape (K, {grid.shape[0]}, cap+1), got {prob.shape}"
            )
        for name, array in (("log10_m", grid), ("probabilities", prob)):
            if array.flags.writeable:
                array = array.copy()
                array.flags.writeable = False
                object.__setattr__(self, name, array)
        # Piecewise cubic coefficients, (K, N_m - 1, 4, cap+1), highest power first.
        spline = CubicSpline(grid, prob, axis=1)
        coefficients = np.ascontiguousarray(np.moveaxis(spline.c, (0, 1), (2, 1)))
        coefficients.flags.writeable = False
        object.__setattr__(self, "_coefficients", coefficients)

    def lookup(self, mean_photons: np.ndarray, size_class: int) -> np.ndarray:
        """``p(n | m)`` for mean absorbed photons ``mean_photons`` (any shape) -> (..., cap+1).

        ``m = 0`` maps to the first node. Raises ``ValueError`` for negative or NaN ``m`` and
        for a class index outside ``0 .. K-1``.
        """
        k = self.probabilities.shape[0]
        if (
            isinstance(size_class, bool)
            or not isinstance(size_class, (int, np.integer))
            or not 0 <= size_class < k
        ):
            raise ValueError(f"size_class must be an integer from 0 to {k - 1}, got {size_class!r}")
        m = np.asarray(mean_photons, dtype=np.float64)
        if np.any(np.isnan(m)) or np.any(m < 0.0):
            raise ValueError("mean_photons must be >= 0 and not NaN")
        grid = self.log10_m
        n_m = grid.shape[0]
        with np.errstate(divide="ignore"):
            u = np.log10(m)
        position = np.clip((u - grid[0]) / (grid[-1] - grid[0]) * (n_m - 1), 0.0, n_m - 1.0)
        lower = np.minimum(np.floor(position).astype(np.intp), n_m - 2)
        x = ((position - lower) * ((grid[-1] - grid[0]) / (n_m - 1)))[..., None]
        c = self._coefficients[int(size_class)][lower]  # (..., 4, cap+1)
        p = ((c[..., 0, :] * x + c[..., 1, :]) * x + c[..., 2, :]) * x + c[..., 3, :]
        p = np.clip(p, 0.0, None)
        return p / p.sum(axis=-1, keepdims=True)


def latent_table(
    classes: SizeClasses,
    params: ChainParams,
    exposure_time_s: float,
    *,
    nodes: int = 201,
    log10_range: tuple[float, float] = (-6.0, 4.0),
    quadrature: int = 9,
) -> LatentTable:
    """Build the latent lookup table for ``classes`` and one exposure time (spec D1 2.2.5).

    One batched ``expm`` call covers all K x nodes x quadrature matrices. Caching per
    ``(material, t_exp, profile)`` is left to the caller.
    """
    t = _check_time(exposure_time_s)
    if isinstance(nodes, bool) or not isinstance(nodes, (int, np.integer)) or nodes < 2:
        raise ValueError(f"nodes must be an integer >= 2, got {nodes!r}")
    lo, hi = (float(v) for v in log10_range)
    if not (math.isfinite(lo) and math.isfinite(hi) and lo < hi):
        raise ValueError(f"log10_range must be finite with low < high, got {log10_range!r}")
    q = _check_quadrature(quadrature)

    log10_m = np.linspace(lo, hi, int(nodes))
    eta = np.asarray(efficiency(params, classes.ecd_um))  # (K,)
    electrons = eta[:, None] * 10.0 ** log10_m[None, :]  # (K, N_m)
    probabilities = np.ascontiguousarray(_distributions(electrons, t, params, q))
    probabilities.flags.writeable = False
    log10_m.flags.writeable = False
    return LatentTable(log10_m=log10_m, probabilities=probabilities, exposure_time_s=t)
