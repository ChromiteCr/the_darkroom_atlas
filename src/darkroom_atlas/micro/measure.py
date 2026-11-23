"""Aperture measurement on the microscope-mode disks (spec D1 2.7.2, step 5).

The transmittance through a circular aperture is its uncovered area fraction. By default the
covered area ``|aperture ∩ union of disks|`` is computed exactly with Green's theorem: the
boundary of the intersection is made of the disk arcs that lie inside the aperture and
outside every other disk, plus the aperture arcs that lie inside some disk. This is the
``p_r -> 0`` limit of the rendered raster. With an explicit ``pitch_um`` the coverage is
rendered instead (:mod:`darkroom_atlas.micro.render`) and each pixel is weighted by the exact
area of its intersection with the aperture circle.
"""

import math

import numpy as np
from scipy.spatial import cKDTree

from darkroom_atlas.micro.grains import GrainSet
from darkroom_atlas.micro.render import _coverage, _disk_pixel_overlap, disk_radii

__all__ = [
    "DEFAULT_APERTURE_UM",
    "MAX_RASTER_PITCH_UM",
    "aperture_transmittances",
    "default_pitch",
]

DEFAULT_APERTURE_UM = 48.0
"""Diameter of the measuring aperture, um (spec D1 2.7.2, step 5; ISO RMS granularity)."""

MAX_RASTER_PITCH_UM = 0.05
"""Upper bound of the raster pitch p_r = min(r_min / 4, 0.05 um) (spec D1 2.7.2, step 5)."""

_TWO_PI = 2.0 * math.pi
# Angle offset between interval groups when merging; larger than any angle in [0, 2 pi].
_GROUP_STRIDE = 8.0
# Apertures handled per batch of the exact path, and raster pixels shared by the raster path.
_APERTURES_PER_BATCH = 256
_DISKS_PER_BLOCK = 200_000
_MAX_SHARED_PIXELS = 8_000_000


def default_pitch(grains: GrainSet, areas_um2: np.ndarray) -> float | None:
    """``p_r = min(r_min / 4, 0.05 um)``, or ``None`` when no grain is developed.

    ``r_c = sqrt(areas_c * <x^2>_c^(1/2) / pi)`` with ``<x^2>_c`` the mean squared developed
    fraction of the developed class-``c`` grains (spec D1 2.7.1, 2.7.2); ``r_min`` is the
    smallest over classes that have developed grains.
    """
    areas = np.asarray(areas_um2, dtype=np.float64)
    disk_radii(grains, areas)  # validates areas against the grains
    x = grains.developed_fraction
    developed = (x > 0.0) & (areas[grains.size_class] > 0.0)
    if not np.any(developed):
        return None
    cls = grains.size_class[developed]
    k = areas.shape[0]
    n = np.bincount(cls, minlength=k)
    mean_x2 = np.bincount(cls, weights=x[developed] ** 2, minlength=k)[n > 0] / n[n > 0]
    r = np.sqrt(areas[n > 0] * np.sqrt(mean_x2) / math.pi)
    return float(min(r.min() / 4.0, MAX_RASTER_PITCH_UM))


# --- exact path ---------------------------------------------------------------------------


def _arcs(
    group: np.ndarray, theta: np.ndarray, cos_half: np.ndarray, full: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Split arcs ``theta +- arccos(cos_half)`` (or whole circles) into ``[lo, hi]`` in [0, 2 pi].

    Arcs with ``cos_half >= 1`` and not ``full`` are empty and dropped.
    """
    keep = full | (cos_half < 1.0)
    group, theta, cos_half, full = group[keep], theta[keep], cos_half[keep], full[keep]
    half = np.where(full, math.pi, np.arccos(np.clip(cos_half, -1.0, 1.0)))
    lo = np.mod(theta - half, _TWO_PI)
    hi = lo + 2.0 * half
    lo = np.where(full, 0.0, lo)
    hi = np.where(full, _TWO_PI, hi)
    wrap = hi > _TWO_PI
    return (
        np.concatenate([group, group[wrap]]),
        np.concatenate([lo, np.zeros(int(wrap.sum()))]),
        np.concatenate([np.minimum(hi, _TWO_PI), hi[wrap] - _TWO_PI]),
    )


def _sort_arcs(
    group: np.ndarray, lo: np.ndarray, hi: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Arcs ordered by group, then start angle."""
    order = np.argsort(group * _GROUP_STRIDE + lo)
    return group[order], lo[order], hi[order]


def _union_sums(
    group: np.ndarray, lo: np.ndarray, hi: np.ndarray, n_groups: int, *, presorted: bool = False
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per group, over the union of its arcs: total angle, sum of d(sin), sum of d(cos)."""
    if group.size == 0:
        return np.zeros(n_groups), np.zeros(n_groups), np.zeros(n_groups)
    if not presorted:
        group, lo, hi = _sort_arcs(group, lo, hi)
    offset = group * _GROUP_STRIDE
    reach = np.maximum.accumulate(hi + offset)  # offsets keep groups apart; only compared
    new = np.ones(group.size, dtype=bool)
    new[1:] = (group[1:] != group[:-1]) | (lo[1:] + offset[1:] > reach[:-1])
    starts = np.flatnonzero(new)
    g, a = group[starts], lo[starts]
    b = np.maximum.reduceat(hi, starts)
    return (
        np.bincount(g, weights=b - a, minlength=n_groups),
        np.bincount(g, weights=np.sin(b) - np.sin(a), minlength=n_groups),
        np.bincount(g, weights=np.cos(b) - np.cos(a), minlength=n_groups),
    )


def _coverage_arcs(
    xi: np.ndarray, yi: np.ndarray, ri: np.ndarray, xj: np.ndarray, yj: np.ndarray, rj: np.ndarray,
    j_first: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:  # fmt: skip
    """Arc of circle ``i`` inside disk ``j``: direction, cos of half-width, whole-circle flag.

    Coincident equal circles: only the one with ``j_first`` set is covered, so the union
    keeps one copy of the boundary.
    """
    dx, dy = xj - xi, yj - yi
    d = np.hypot(dx, dy)
    zero = d <= 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        cos_half = (ri * ri + d * d - rj * rj) / (2.0 * ri * d)
    cos_half = np.where(zero, 1.0, cos_half)
    full = np.where(zero, (ri < rj) | ((ri == rj) & j_first), cos_half <= -1.0)
    return np.arctan2(dy, dx), cos_half, full


def _exact(
    x: np.ndarray, y: np.ndarray, r: np.ndarray, centres: np.ndarray, radius: float
) -> np.ndarray:
    """Exact transmittances, apertures grouped in square tiles to bound memory."""
    reach = radius + float(r.max())
    lo_x, hi_x = float(x.min()), float(x.max())
    lo_y, hi_y = float(y.min()), float(y.max())
    density = x.size / max((hi_x - lo_x + 2.0 * reach) * (hi_y - lo_y + 2.0 * reach), 1e-300)
    tile = max(math.sqrt(_DISKS_PER_BLOCK / density) - 2.0 * reach, radius)
    key = np.floor((centres - centres.min(axis=0)) / tile).astype(np.int64)
    _, block = np.unique(key, axis=0, return_inverse=True)
    block = block.ravel()
    order = np.argsort(x, kind="stable")
    x, y, r = x[order], y[order], r[order]
    out = np.ones(centres.shape[0])
    for b in range(int(block.max()) + 1):
        sel = np.flatnonzero(block == b)
        cb = centres[sel]
        a, z = np.searchsorted(x, [cb[:, 0].min() - reach, cb[:, 0].max() + reach])
        ys = y[a:z]
        near = (ys > cb[:, 1].min() - reach) & (ys < cb[:, 1].max() + reach)
        out[sel] = _exact_block(x[a:z][near], ys[near], r[a:z][near], cb, radius)
    return out


def _exact_block(
    x: np.ndarray, y: np.ndarray, r: np.ndarray, centres: np.ndarray, radius: float
) -> np.ndarray:
    """Exact transmittances of apertures of ``radius`` at ``centres`` over disks ``(x, y, r)``."""
    m = centres.shape[0]
    out = np.ones(m)
    if x.size == 0:
        return out
    r_max = float(r.max())
    # Only disks that meet some aperture shape the measured boundary (their nearest centre).
    near, _ = cKDTree(centres).query(np.stack([x, y], axis=1), distance_upper_bound=radius + r_max)
    meets = near < radius + r
    x, y, r = x[meets], y[meets], r[meets]
    n = x.size
    if n == 0:
        return out

    # Arcs of each disk covered by the other disks (independent of the aperture).
    tree = cKDTree(np.stack([x, y], axis=1))
    pairs = tree.query_pairs(2.0 * r_max, output_type="ndarray")
    i, j = pairs[:, 0], pairs[:, 1]
    touch = np.hypot(x[j] - x[i], y[j] - y[i]) < r[i] + r[j]
    i, j = i[touch], j[touch]
    owner = np.concatenate([i, j])
    other = np.concatenate([j, i])
    theta, cos_half, full = _coverage_arcs(
        x[owner], y[owner], r[owner], x[other], y[other], r[other], other < owner
    )
    nb_group, nb_lo, nb_hi = _arcs(owner, theta, cos_half, full)
    nb_group, nb_lo, nb_hi = _sort_arcs(nb_group, nb_lo, nb_hi)
    nb_count = np.bincount(nb_group, minlength=n)
    nb_start = np.cumsum(nb_count) - nb_count
    cov_len, cov_sin, cov_cos = _union_sums(nb_group, nb_lo, nb_hi, n, presorted=True)
    # Uncovered arcs = whole circle minus the covered union.
    free_len, free_sin, free_cos = _TWO_PI - cov_len, -cov_sin, -cov_cos

    for q0 in range(0, m, _APERTURES_PER_BATCH):
        cq = centres[q0 : q0 + _APERTURES_PER_BATCH]
        lists = tree.query_ball_point(cq, radius + r_max)
        sizes = np.fromiter((len(v) for v in lists), dtype=np.int64, count=len(lists))
        if not sizes.any():
            continue
        disk = np.concatenate([np.asarray(v, dtype=np.int64) for v in lists])
        ap = np.repeat(np.arange(cq.shape[0]), sizes)
        ox, oy = cq[ap, 0], cq[ap, 1]
        rx, ry, rr = x[disk] - ox, y[disk] - oy, r[disk]
        d = np.hypot(rx, ry)
        hit = d < radius + rr
        disk, ap, ox, oy, rx, ry, rr, d = (v[hit] for v in (disk, ap, ox, oy, rx, ry, rr, d))
        inside = d + rr <= radius

        # Disks wholly inside: their free arcs, with the aperture centre as Green origin.
        k = disk[inside]
        val = (
            0.5
            * rr[inside]
            * (rr[inside] * free_len[k] + rx[inside] * free_sin[k] - ry[inside] * free_cos[k])
        )
        area = np.zeros(cq.shape[0])
        area += np.bincount(ap[inside], weights=val, minlength=cq.shape[0])

        # Straddling disks: free arcs that also lie inside the aperture.
        s = np.flatnonzero(~inside)
        if s.size:
            ks = disk[s]
            with np.errstate(divide="ignore", invalid="ignore"):
                cos_out = (radius * radius - d[s] ** 2 - rr[s] ** 2) / (2.0 * d[s] * rr[s])
            zero = d[s] <= 0.0
            cos_out = np.where(zero, 1.0, cos_out)
            full_out = np.where(zero, rr[s] > radius, cos_out <= -1.0)
            ex_group, ex_lo, ex_hi = _arcs(
                np.arange(s.size), np.arctan2(ry[s], rx[s]), cos_out, full_out
            )
            cnt = nb_count[ks]
            total = int(cnt.sum())
            rep = np.repeat(np.arange(s.size), cnt)
            idx = (
                np.repeat(nb_start[ks], cnt)
                + np.arange(total)
                - np.repeat(np.cumsum(cnt) - cnt, cnt)
            )
            group = np.concatenate([rep, ex_group])
            lo = np.concatenate([nb_lo[idx], ex_lo])
            hi = np.concatenate([nb_hi[idx], ex_hi])
            c_len, c_sin, c_cos = _union_sums(group, lo, hi, s.size)
            val = 0.5 * rr[s] * (
                rr[s] * (_TWO_PI - c_len) - rx[s] * c_sin + ry[s] * c_cos
            )  # fmt: skip
            area += np.bincount(ap[s], weights=val, minlength=cq.shape[0])

            # Aperture arcs inside some disk (Green origin at the centre: r^2 dtheta / 2).
            with np.errstate(divide="ignore", invalid="ignore"):
                cos_in = (radius * radius + d[s] ** 2 - rr[s] ** 2) / (2.0 * radius * d[s])
            cos_in = np.where(zero, 1.0, cos_in)
            full_in = np.where(zero, rr[s] > radius, cos_in <= -1.0)
            a_group, a_lo, a_hi = _arcs(ap[s], np.arctan2(ry[s], rx[s]), cos_in, full_in)
            a_len, _, _ = _union_sums(a_group, a_lo, a_hi, cq.shape[0])
            area += 0.5 * radius * radius * a_len

        covered = np.clip(area / (math.pi * radius * radius), 0.0, 1.0)
        out[q0 : q0 + cq.shape[0]] = 1.0 - covered
    return out


# --- raster path --------------------------------------------------------------------------


def _aperture_weights(
    edge_x: np.ndarray, edge_y: np.ndarray, radius: float, pitch: float
) -> np.ndarray:
    """Exact pixel fractions inside a circle; pixel edges are relative to its centre."""
    cx = 0.5 * (edge_x[:-1] + edge_x[1:])
    cy = 0.5 * (edge_y[:-1] + edge_y[1:])
    d = np.hypot(cx[None, :], cy[:, None])
    half_diag = pitch / math.sqrt(2.0)
    weights = (d + half_diag <= radius).astype(np.float64)
    ring = np.abs(d - radius) < half_diag
    jj, ii = np.nonzero(ring)
    if jj.size:
        dx = np.stack([edge_x[ii], edge_x[ii + 1]], axis=1)
        dy = np.stack([edge_y[jj], edge_y[jj + 1]], axis=1)
        area = _disk_pixel_overlap(dx, dy, np.full(jj.size, radius))[:, 0, 0]
        weights[jj, ii] = area / (pitch * pitch)
    return weights


def _rendered(
    x: np.ndarray, y: np.ndarray, r: np.ndarray, centres: np.ndarray, radius: float, pitch: float
) -> np.ndarray:
    """Transmittances from the coverage raster at ``pitch``."""
    m = centres.shape[0]
    gx0 = float(centres[:, 0].min()) - radius
    gy0 = float(centres[:, 1].min()) - radius
    i0 = np.floor((centres[:, 0] - radius - gx0) / pitch).astype(np.int64)
    j0 = np.floor((centres[:, 1] - radius - gy0) / pitch).astype(np.int64)
    i1 = np.ceil((centres[:, 0] + radius - gx0) / pitch).astype(np.int64)
    j1 = np.ceil((centres[:, 1] + radius - gy0) / pitch).astype(np.int64)

    order = np.argsort(y, kind="stable")
    gx, gy, gr = x[order], y[order], r[order]
    reach = float(gr.max()) if gr.size else 0.0

    # Apertures sorted by row share band rasters of at most _MAX_SHARED_PIXELS pixels.
    out = np.empty(m)
    band: list[int] = []
    ia = ib = ja = jb = 0
    for a in [*np.argsort(j0, kind="stable").tolist(), None]:
        if a is not None:
            box = (min(ia, i0[a]), max(ib, i1[a]), min(ja, j0[a]), max(jb, j1[a]))
            if not band:
                box = (i0[a], i1[a], j0[a], j1[a])
            if not band or (box[1] - box[0]) * (box[3] - box[2]) <= _MAX_SHARED_PIXELS:
                band.append(a)
                ia, ib, ja, jb = (int(v) for v in box)
                continue
        lo, hi = np.searchsorted(gy, [gy0 + ja * pitch - reach, gy0 + jb * pitch + reach])
        cov = _coverage(
            gx[lo:hi], gy[lo:hi], gr[lo:hi], (gx0 + ia * pitch, gy0 + ja * pitch),
            (jb - ja, ib - ia), pitch,
        )  # fmt: skip
        for b in band:
            cx, cy = centres[b]
            ex = gx0 + np.arange(i0[b], i1[b] + 1) * pitch - cx
            ey = gy0 + np.arange(j0[b], j1[b] + 1) * pitch - cy
            w = _aperture_weights(ex, ey, radius, pitch)
            window = cov[j0[b] - ja : j1[b] - ja, i0[b] - ia : i1[b] - ia]
            out[b] = float(np.sum(w * (1.0 - window)) / np.sum(w))
        if a is not None:
            band = [a]
            ia, ib, ja, jb = int(i0[a]), int(i1[a]), int(j0[a]), int(j1[a])
    return out


def aperture_transmittances(
    grains: GrainSet,
    areas_um2: np.ndarray,
    centres_um: np.ndarray,
    *,
    aperture_um: float = DEFAULT_APERTURE_UM,
    pitch_um: float | None = None,
) -> np.ndarray:
    """Transmittance through circular apertures of diameter ``aperture_um`` at ``centres_um``.

    ``centres_um`` is ``(M, 2)`` absolute ``(x, y)`` in um; returns ``(M,)``. The disks are
    those of :func:`darkroom_atlas.micro.render.render_coverage`; grains must be placed at
    least one disk radius beyond every aperture, otherwise the edge reads too clear.
    ``pitch_um=None`` computes the covered area exactly; a pitch renders the coverage at that
    pitch (:func:`default_pitch` gives the spec's ``p_r``). With no developed grain every
    transmittance is 1.
    """
    centres = np.asarray(centres_um, dtype=np.float64)
    if centres.ndim != 2 or centres.shape[1] != 2 or not np.all(np.isfinite(centres)):
        raise ValueError(f"centres_um must be a finite (M, 2) array, got shape {centres.shape}")
    if isinstance(aperture_um, bool) or not math.isfinite(aperture_um) or aperture_um <= 0.0:
        raise ValueError(f"aperture_um must be positive and finite, got {aperture_um!r}")
    if pitch_um is not None and (
        isinstance(pitch_um, bool) or not math.isfinite(pitch_um) or pitch_um <= 0.0
    ):
        raise ValueError(f"pitch_um must be positive and finite, got {pitch_um!r}")
    r = disk_radii(grains, areas_um2)
    m = centres.shape[0]
    drawn = r > 0.0
    if m == 0 or not np.any(drawn):
        return np.ones(m)
    x, y, r = grains.x_um[drawn], grains.y_um[drawn], r[drawn]
    radius = 0.5 * float(aperture_um)
    if pitch_um is None:
        return _exact(x, y, r, centres, radius)
    return _rendered(x, y, r, centres, radius, float(pitch_um))
