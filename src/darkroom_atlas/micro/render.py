"""Coverage raster of developed grains as opaque disks (spec D1 2.7.2, steps 4 and 5).

A developed grain is a disk of area ``areas_um2[size_class] * developed_fraction`` centred on
the grain. Each disk's share of a raster pixel is the exact area of the disk-square
intersection, so a lone disk's total coverage equals its area to rounding. Overlapping disks
cover once (Boolean model): a pixel's uncovered fraction is the product of ``1 - c_i`` over
the disks touching it, which keeps the mean uncovered fraction ``exp(-lambda a)`` exact for
a Poisson field at any pitch (Campbell's theorem).
"""

import math

import numpy as np

from darkroom_atlas.micro.grains import GrainSet

__all__ = ["render_coverage"]

# Corner evaluations per vectorised batch, to bound memory.
_BATCH_VALUES = 2_000_000
# Per-disk pixel fractions this close to 1 count as full cover; below _TINY they are dropped.
_FULL = 1.0 - 1e-12
_TINY = 1e-15


def _primitive(t: np.ndarray, r: np.ndarray) -> np.ndarray:
    """``integral_0^t sqrt(r^2 - s^2) ds`` for ``|t| <= r``."""
    t = np.clip(t, -r, r)
    return 0.5 * (t * np.sqrt(np.maximum(r * r - t * t, 0.0)) + r * r * np.arcsin(t / r))


def _half_plane(x: np.ndarray, r: np.ndarray) -> np.ndarray:
    """Area of the disk ``|X| <= r`` centred at 0 with ``X >= x``."""
    return 0.5 * math.pi * r * r - 2.0 * _primitive(x, r)


def _corner(x: np.ndarray, y: np.ndarray, r: np.ndarray) -> np.ndarray:
    """Area of the disk with ``X >= x, Y >= y`` for ``x, y >= 0``."""
    inside = x * x + y * y < r * r
    a = np.sqrt(np.maximum(r * r - y * y, 0.0))
    area = _primitive(a, r) - _primitive(x, r) - y * (a - x)
    return np.where(inside, np.maximum(area, 0.0), 0.0)


def _quadrant(x: np.ndarray, y: np.ndarray, r: np.ndarray) -> np.ndarray:
    """Area of the disk of radius ``r`` centred at 0 with ``X >= x, Y >= y`` (any signs)."""
    c = _corner(np.abs(x), np.abs(y), r)
    sx, sy = x < 0.0, y < 0.0
    full = math.pi * r * r
    return np.where(
        sx,
        np.where(sy, full - _half_plane(-x, r) - _half_plane(-y, r) + c, _half_plane(y, r) - c),
        np.where(sy, _half_plane(x, r) - c, c),
    )


def _disk_pixel_overlap(dx_edges: np.ndarray, dy_edges: np.ndarray, r: np.ndarray) -> np.ndarray:
    """Exact areas of disks ``(n,)`` within the pixels of their own grids.

    ``dx_edges`` ``(n, bx + 1)`` and ``dy_edges`` ``(n, by + 1)`` are pixel edges relative to
    each disk centre. Returns ``(n, by, bx)`` areas in um2.
    """
    rr = r[:, None, None]
    s = _quadrant(dx_edges[:, None, :], dy_edges[:, :, None], rr)  # (n, by + 1, bx + 1)
    return s[:, :-1, :-1] - s[:, :-1, 1:] - s[:, 1:, :-1] + s[:, 1:, 1:]


def _coverage(
    cx: np.ndarray,
    cy: np.ndarray,
    r: np.ndarray,
    origin: tuple[float, float],
    shape: tuple[int, int],
    pitch: float,
) -> np.ndarray:
    """Covered fraction ``(ny, nx)`` of disks ``(cx, cy, r)`` on the raster at ``origin``."""
    x0, y0 = origin
    ny, nx = shape
    keep = (
        (r > 0.0)
        & (cx + r > x0)
        & (cx - r < x0 + nx * pitch)
        & (cy + r > y0)
        & (cy - r < y0 + ny * pitch)
    )
    cx, cy, r = cx[keep], cy[keep], r[keep]
    i_lo = np.floor((cx - r - x0) / pitch).astype(np.int64)
    j_lo = np.floor((cy - r - y0) / pitch).astype(np.int64)
    span = np.maximum(
        np.floor((cx + r - x0) / pitch).astype(np.int64) - i_lo,
        np.floor((cy + r - y0) / pitch).astype(np.int64) - j_lo,
    )
    span += 1  # box side in pixels

    size = nx * ny
    full_count = np.zeros(size, dtype=np.int64)
    log_uncovered = np.zeros(size, dtype=np.float64)
    pending: list[tuple[np.ndarray, np.ndarray]] = []
    pending_n = 0

    def flush() -> None:
        nonlocal full_count, log_uncovered, pending_n
        if pending:
            flat = np.concatenate([p[0] for p in pending])
            frac = np.concatenate([p[1] for p in pending])
            full = frac >= _FULL
            full_count += np.bincount(flat[full], minlength=size)
            log_uncovered += np.bincount(
                flat[~full], weights=np.log1p(-frac[~full]), minlength=size
            )
            pending.clear()
            pending_n = 0

    for b in np.unique(span):
        members = np.flatnonzero(span == b)
        step = max(1, _BATCH_VALUES // int((b + 1) * (b + 1)))
        k = np.arange(b + 1)
        for start in range(0, members.size, step):
            m = members[start : start + step]
            dx = x0 + (i_lo[m, None] + k) * pitch - cx[m, None]
            dy = y0 + (j_lo[m, None] + k) * pitch - cy[m, None]
            frac = _disk_pixel_overlap(dx, dy, r[m]) / (pitch * pitch)
            cols = i_lo[m, None, None] + k[None, None, :-1]
            rows = j_lo[m, None, None] + k[None, :-1, None]
            valid = (cols >= 0) & (cols < nx) & (rows >= 0) & (rows < ny) & (frac > _TINY)
            pending.append(((rows * nx + cols)[valid], frac[valid]))
            pending_n += pending[-1][0].size
            # Amortise the raster-sized bincount over at least as many entries.
            if pending_n >= max(size, _BATCH_VALUES):
                flush()
    flush()
    covered = np.where(full_count > 0, 1.0, -np.expm1(log_uncovered))
    return covered.reshape(ny, nx)


def disk_radii(grains: GrainSet, areas_um2: np.ndarray) -> np.ndarray:
    """Disk radius of every grain, ``sqrt(areas_um2[size_class] * developed_fraction / pi)``."""
    areas = np.asarray(areas_um2, dtype=np.float64)
    if areas.ndim != 1 or not np.all(np.isfinite(areas)) or np.any(areas < 0.0):
        raise ValueError("areas_um2 must be a finite, non-negative 1-D array")
    if grains.count and int(grains.size_class.max()) >= areas.shape[0]:
        raise ValueError(
            f"areas_um2 has {areas.shape[0]} classes but a grain has size class "
            f"{int(grains.size_class.max())}"
        )
    return np.sqrt(areas[grains.size_class] * grains.developed_fraction / math.pi)


def render_coverage(
    grains: GrainSet,
    areas_um2: np.ndarray,
    window_um: tuple[float, float, float, float],
    pitch_um: float,
) -> np.ndarray:
    """Anti-aliased covered fraction per raster pixel (spec D1 2.7.2, step 5).

    ``window_um`` is ``(x0, y0, x1, y1)`` in film um; its sides must be whole multiples of
    ``pitch_um``. The result has shape ``(ny, nx)``, rows along ``y``; pixel ``[j, i]`` spans
    ``[x0 + i p, x0 + (i + 1) p) x [y0 + j p, y0 + (j + 1) p)``. Disks of area
    ``areas_um2[size_class] * developed_fraction`` (the caller folds in ``1 + phi``) cover
    once where they overlap. Disks centred outside the window still count where they reach in.
    """
    try:
        x0, y0, x1, y1 = (float(v) for v in window_um)
    except (TypeError, ValueError):
        raise ValueError(f"window_um must be (x0, y0, x1, y1), got {window_um!r}") from None
    if not all(math.isfinite(v) for v in (x0, y0, x1, y1)) or not (x1 > x0 and y1 > y0):
        raise ValueError(f"window_um must be finite with x1 > x0 and y1 > y0, got {window_um!r}")
    if isinstance(pitch_um, bool) or not math.isfinite(pitch_um) or pitch_um <= 0.0:
        raise ValueError(f"pitch_um must be positive and finite, got {pitch_um!r}")
    pitch = float(pitch_um)
    shape = []
    for side in (y1 - y0, x1 - x0):
        n = round(side / pitch)
        if n < 1 or abs(n * pitch - side) > 1e-6 * pitch:
            raise ValueError(f"window side {side} um is not a whole multiple of pitch {pitch}")
        shape.append(int(n))
    r = disk_radii(grains, areas_um2)
    return _coverage(grains.x_um, grains.y_um, r, (x0, y0), (shape[0], shape[1]), pitch)
