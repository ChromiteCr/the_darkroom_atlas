"""Explicit grain placement for the microscope mode (spec D1 2.7.2, step 1).

Grains of each size class form a Poisson process of surface density ``nu_c`` with fractional
depth uniform in ``[0, 1)``. The film plane is cut into square micro-cells of side
:data:`MICRO_CELL_UM`; the grains of a cell depend only on ``(seed, class, cell row, cell
column)``, so a cell holds the same grains whatever patch origin or size is requested.
"""

import math
from dataclasses import dataclass
from numbers import Integral

import numpy as np

from darkroom_atlas.emulsion.classes import SizeClasses
from darkroom_atlas.grain.poisson import poisson_inverse
from darkroom_atlas.rng import pixel_uniform

__all__ = ["MAX_PATCH_UM", "MICRO_CELL_UM", "GrainSet", "place_grains"]

MICRO_CELL_UM = 4.0
"""Side of a micro-cell, um: the 4 um chemistry cell of the microscope mode (spec D1 2.7.2)."""

MAX_PATCH_UM = 2000.0
"""Largest patch side, um: the microscope mode covers patches up to 2 mm (spec D1 2.7.2)."""

# Values drawn per chunk of cell rows, to bound memory on large patches.
_CHUNK_VALUES = 4_000_000


@dataclass(frozen=True, eq=False)
class GrainSet:
    """Explicit grains of a patch; all arrays have shape ``(G,)`` and are read-only.

    Positions are absolute film coordinates in um (``x`` along columns, ``y`` along rows).
    ``depth_index`` is the 0-based depth cell, ``size_class`` the 0-based size class.
    ``cluster_size`` is 0 (no developable cluster) or 2..cap, ``init_time_s`` is ``inf`` for
    a grain that never initiated, ``developed_fraction`` lies in ``[0, 1]``.
    """

    x_um: np.ndarray
    y_um: np.ndarray
    depth_index: np.ndarray
    size_class: np.ndarray
    cluster_size: np.ndarray
    init_time_s: np.ndarray
    developed_fraction: np.ndarray

    def __post_init__(self) -> None:
        floats = ("x_um", "y_um", "init_time_s", "developed_fraction")
        ints = ("depth_index", "size_class", "cluster_size")
        n = None
        for name in floats + ints:
            value = np.asarray(getattr(self, name))
            if value.ndim != 1:
                raise ValueError(f"GrainSet.{name} must be 1-D, got shape {value.shape}")
            if n is None:
                n = value.shape[0]
            elif value.shape[0] != n:
                raise ValueError(f"GrainSet.{name} has {value.shape[0]} entries, expected {n}")
            if name in ints:
                if value.size and not np.issubdtype(value.dtype, np.integer):
                    raise ValueError(f"GrainSet.{name} must be integer")
                value = np.array(value, dtype=np.int64)
                if np.any(value < 0):
                    raise ValueError(f"GrainSet.{name} must not be negative")
            else:
                value = np.array(value, dtype=np.float64)
            value.flags.writeable = False
            object.__setattr__(self, name, value)
        for name in ("x_um", "y_um"):
            if not np.all(np.isfinite(getattr(self, name))):
                raise ValueError(f"GrainSet.{name} must be finite")
        if np.any(np.isnan(self.init_time_s)) or np.any(self.init_time_s < 0.0):
            raise ValueError("GrainSet.init_time_s must be non-negative (inf if never initiated)")
        x = self.developed_fraction
        if not np.all(np.isfinite(x)) or np.any((x < 0.0) | (x > 1.0)):
            raise ValueError("GrainSet.developed_fraction must lie in [0, 1]")

    @property
    def count(self) -> int:
        """Number of grains ``G``."""
        return int(self.x_um.shape[0])


def _cap(mean: float) -> int:
    """Largest Poisson count drawn for ``mean``: floor(mu + 10 sqrt(mu) + 10) (spec D1 2.7.1)."""
    return math.floor(mean + 10.0 * math.sqrt(mean) + 10.0)


def _pair(name: str, value: object, *, positive: bool) -> tuple[float, float]:
    try:
        a, b = (float(v) for v in value)  # type: ignore[union-attr]
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a pair of numbers, got {value!r}") from None
    for v in (a, b):
        if not math.isfinite(v) or (v <= 0.0 if positive else v < 0.0):
            kind = "positive" if positive else "non-negative"
            raise ValueError(f"{name} must be {kind} and finite, got {value!r}")
    return a, b


def _cells(lo: float, extent: float) -> tuple[int, int]:
    """Half-open range of micro-cell indices that meet ``[lo, lo + extent)``."""
    first = math.floor(lo / MICRO_CELL_UM)
    stop = math.ceil((lo + extent) / MICRO_CELL_UM)
    return first, max(stop, first + 1)


def place_grains(
    classes: SizeClasses,
    patch_um: tuple[float, float],
    n_z: int,
    seed: int,
    origin_um: tuple[float, float] = (0.0, 0.0),
) -> GrainSet:
    """Place the grains of the patch ``[x0, x0 + w) x [y0, y0 + h)`` (spec D1 2.7.2, step 1).

    ``patch_um`` is ``(w, h)`` and ``origin_um`` is ``(x0, y0)``, both in um; the origin must
    be non-negative because micro-cells are indexed from the film origin. In each micro-cell
    the number of class ``c`` grains is Poisson with mean ``nu_c * MICRO_CELL_UM^2`` (stream
    ``micro/<c>``); each grain's ``x``, ``y`` and fractional depth are uniform in the cell
    (stream ``micro/position/<c>``), and ``depth_index = floor(zeta * n_z)``. Grains come in
    order of cell row, cell column, class and slot. Development fields start empty:
    ``cluster_size`` 0, ``init_time_s`` inf, ``developed_fraction`` 0.
    """
    w, h = _pair("patch_um", patch_um, positive=True)
    if w > MAX_PATCH_UM or h > MAX_PATCH_UM:
        raise ValueError(f"patch_um sides must not exceed {MAX_PATCH_UM} um, got {patch_um!r}")
    x0, y0 = _pair("origin_um", origin_um, positive=False)
    if isinstance(n_z, bool) or not isinstance(n_z, Integral) or n_z < 1:
        raise ValueError(f"n_z must be a positive integer, got {n_z!r}")
    if isinstance(seed, bool) or not isinstance(seed, Integral) or seed < 0:
        raise ValueError(f"seed must be a non-negative integer, got {seed!r}")
    n_z, seed = int(n_z), int(seed)
    density = np.asarray(classes.number_per_um2, dtype=np.float64)
    if density.ndim != 1 or not np.all(np.isfinite(density)) or np.any(density < 0.0):
        raise ValueError("classes.number_per_um2 must be a finite, non-negative 1-D array")

    c0, c1 = _cells(x0, w)
    r0, r1 = _cells(y0, h)
    n_cols = c1 - c0
    means = density * MICRO_CELL_UM * MICRO_CELL_UM
    caps = [_cap(m) for m in means]
    rows_per_chunk = max(1, _CHUNK_VALUES // (n_cols * (1 + 3 * max(caps, default=1))))

    parts: list[tuple[np.ndarray, ...]] = []
    for ra in range(r0, r1, rows_per_chunk):
        rb = min(ra + rows_per_chunk, r1)
        for c, (mean, cap) in enumerate(zip(means, caps, strict=True)):
            u = pixel_uniform(seed, f"micro/{c}", (ra, rb), (c0, c1), c1)
            counts = poisson_inverse(u, np.full(u.shape, mean))
            if not np.any(counts):
                continue
            pos = pixel_uniform(seed, f"micro/position/{c}", (ra, rb), (c0, c1), c1, (cap, 3))
            keep = np.arange(cap) < counts[..., None]  # (rows, cols, cap)
            ri, ci, _ = np.nonzero(keep)
            p = pos[keep]  # (n, 3)
            x = (c0 + ci + p[:, 0]) * MICRO_CELL_UM
            y = (ra + ri + p[:, 1]) * MICRO_CELL_UM
            depth = np.minimum(np.floor(p[:, 2] * n_z).astype(np.int64), n_z - 1)
            inside = (x >= x0) & (x < x0 + w) & (y >= y0) & (y < y0 + h)
            # Sort key: cell row, cell column, class, slot (np.nonzero is row-major).
            order_key = (ra + ri) * n_cols + ci
            parts.append(
                (order_key[inside], np.full(inside.sum(), c), x[inside], y[inside], depth[inside])
            )

    if parts:
        key, cls, x, y, depth = (np.concatenate(a) for a in zip(*parts, strict=True))
        order = np.lexsort((cls, key))  # stable: slots stay in order within a cell and class
        key, cls, x, y, depth = key[order], cls[order], x[order], y[order], depth[order]
    else:
        cls = depth = np.zeros(0, dtype=np.int64)
        x = y = np.zeros(0, dtype=np.float64)
    g = x.shape[0]
    return GrainSet(
        x_um=x,
        y_um=y,
        depth_index=depth,
        size_class=cls.astype(np.int64),
        cluster_size=np.zeros(g, dtype=np.int64),
        init_time_s=np.full(g, np.inf),
        developed_fraction=np.zeros(g, dtype=np.float64),
    )
