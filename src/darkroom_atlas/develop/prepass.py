"""Tank coupling of a whole frame: exposure histogram, prepass and its check (spec D1 3.7).

Tank depletion depends only on the frame's mean source term, which the exposure histogram
nearly fixes. So:

1. :func:`exposure_histogram`: key ``s = log10 m-bar`` with ``m-bar = (1/n_z) sum_z sum_c
   w_c (V_c e_V,z + S_c e_S,z)`` (``w_c`` the number fractions), 64 equal-width bins over
   the frame's range of ``s`` plus a zero-exposure bin (``m-bar < 1e-6``); integer counts and
   exposure sums, accumulated over canonical blocks in row-major order.
2. :func:`tank_prepass`: the 65 columns (bin means; the zero column has no exposure),
   weighted ``f_img count / total`` and the zero column ``+ (1 - f_img)``, exactly coupled
   to the tank. Gives the bulk trajectory and ``Jbar_zero``, the flux of unexposed film
   under it (one column run against the trajectory).
3. The image blocks develop against that trajectory and sum their inner-cell fluxes.
4. :func:`prepass_check`: re-integrate the tank from ``f_img Jbar_blk + (1 - f_img)
   Jbar_zero``, rerun the 65 columns under both trajectories and compare their developed
   silver ``M_b`` by ``delta`` (pass at ``delta <= 1e-3``).

Empty bins get the exposure of their bin centre (uniform in depth on the volume
descriptor); their weight is 0, so they never touch the tank. A frame whose non-zero pixels
share one key gets the bins ``[s - 0.5, s + 0.5]``. Strips use :func:`coupled_columns`, the
same exactly coupled column code (spec D1 3.7, last paragraph).
"""

import math
from dataclasses import dataclass

import numpy as np

from darkroom_atlas.chemistry.kinetics import Kinetics, positive_part
from darkroom_atlas.chemistry.solution import WorkingSolution
from darkroom_atlas.chemistry.speciation import solve_ph
from darkroom_atlas.develop import tank as tk
from darkroom_atlas.develop.ensemble import (
    BulkTrajectory,
    EnsembleResult,
    FilmColumnModel,
    Members,
    TankSpec,
    run_ensemble,
)
from darkroom_atlas.develop.schedule import Schedule
from darkroom_atlas.emulsion.classes import SizeClasses

__all__ = [
    "HISTOGRAM_BINS",
    "ZERO_EXPOSURE",
    "ExposureHistogram",
    "PrepassCheck",
    "bulk_ph",
    "coupled_columns",
    "developed_silver",
    "exposure_histogram",
    "prepass_check",
    "tank_prepass",
    "trajectory_array",
    "zero_column",
]

HISTOGRAM_BINS = 64
"""Equal-width bins of the exposure key (spec D1 3.7 step 1)."""

ZERO_EXPOSURE = 1e-6
"""Pixels with ``m-bar`` below this go to the zero-exposure bin (spec D1 3.7 step 1)."""


@dataclass(frozen=True, eq=False)
class ExposureHistogram:
    """Histogram of the exposure key ``log10 m-bar`` (spec D1 3.7, 10).

    ``edges`` ``(bins + 1,)``, ``counts`` ``(bins,)`` int64, ``mean_volume`` and
    ``mean_surface`` ``(bins, n_z)`` the bins' mean exposure descriptors (empty bins: the bin
    centre, see the module docs), ``zero_count`` the pixels of the zero-exposure bin.
    """

    edges: np.ndarray
    counts: np.ndarray
    mean_volume: np.ndarray
    mean_surface: np.ndarray
    zero_count: int

    def __post_init__(self) -> None:
        set_ = object.__setattr__
        set_(self, "edges", _readonly(np.asarray(self.edges, dtype=np.float64)))
        set_(self, "counts", _readonly(np.asarray(self.counts, dtype=np.int64)))
        set_(self, "mean_volume", _readonly(np.asarray(self.mean_volume, dtype=np.float64)))
        set_(self, "mean_surface", _readonly(np.asarray(self.mean_surface, dtype=np.float64)))
        set_(self, "zero_count", int(self.zero_count))
        b = self.counts.shape[0]
        if self.edges.shape != (b + 1,) or self.mean_volume.shape[0] != b:
            raise ValueError("edges, counts and mean exposures disagree on the number of bins")
        if self.mean_surface.shape != self.mean_volume.shape or self.mean_volume.ndim != 2:
            raise ValueError("mean_volume and mean_surface must share a shape (bins, n_z)")
        for name in ("edges", "mean_volume", "mean_surface"):
            if not np.all(np.isfinite(getattr(self, name))):
                raise ValueError(f"{name} must be finite")
        if np.any(self.counts < 0) or self.zero_count < 0:
            raise ValueError("counts must be >= 0")
        if int(self.counts.sum()) + self.zero_count < 1:
            raise ValueError("the histogram holds no pixel")

    @property
    def bins(self) -> int:
        """Number of bins without the zero bin."""
        return int(self.counts.shape[0])

    @property
    def total(self) -> int:
        """All pixels, zero bin included."""
        return int(self.counts.sum()) + self.zero_count


def _readonly(a: np.ndarray) -> np.ndarray:
    a = np.ascontiguousarray(a)
    a.flags.writeable = False
    return a


def _mean_hits(ev: np.ndarray, es: np.ndarray, classes: SizeClasses) -> np.ndarray:
    """``m-bar`` per pixel, ``(..., n_z) -> (...)``; element-wise sums only (spec D1 3.6)."""
    w = np.asarray(classes.number_fraction, dtype=np.float64)
    wv = float(np.sum(w * classes.volume_um3))
    ws = float(np.sum(w * classes.surface_um2))
    return np.sum(wv * ev + ws * es, axis=-1) / ev.shape[-1]


def exposure_histogram(
    exposure_volume: np.ndarray,
    exposure_surface: np.ndarray,
    classes: SizeClasses,
    *,
    bins: int = HISTOGRAM_BINS,
    block_px: int = 512,
) -> ExposureHistogram:
    """Histogram of ``log10 m-bar`` over a frame (spec D1 3.7 step 1).

    ``exposure_*`` ``(H, W, n_z)`` (memory maps are read block by block). Sums run over
    ``block_px`` blocks anchored at the origin in row-major order, so the result does not
    depend on how the fields were produced. ``ValueError`` on bad shapes, negative or
    non-finite exposures, ``bins < 1`` or ``block_px < 1``.
    """
    ev_all, es_all = exposure_volume, exposure_surface
    if not isinstance(classes, SizeClasses):
        raise TypeError(f"classes must be SizeClasses, got {type(classes).__name__}")
    if (
        getattr(ev_all, "ndim", 0) != 3
        or getattr(es_all, "shape", None) != ev_all.shape
        or ev_all.shape[0] < 1
        or ev_all.shape[1] < 1
        or ev_all.shape[2] < 1
    ):
        raise ValueError("exposure_volume and exposure_surface must share a shape (H, W, n_z)")
    for name, value in (("bins", bins), ("block_px", block_px)):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError(f"{name} must be an int >= 1, got {value!r}")
    h, w, n_z = ev_all.shape

    def blocks():
        for r0 in range(0, h, block_px):
            for c0 in range(0, w, block_px):
                sl = (slice(r0, min(r0 + block_px, h)), slice(c0, min(c0 + block_px, w)))
                ev = np.asarray(ev_all[sl], dtype=np.float64)
                es = np.asarray(es_all[sl], dtype=np.float64)
                if not (np.all(np.isfinite(ev)) and np.all(np.isfinite(es))):
                    raise ValueError("exposures must be finite")
                if np.any(ev < 0.0) or np.any(es < 0.0):
                    raise ValueError("exposures must be >= 0")
                yield ev, es, _mean_hits(ev, es, classes)

    lo, hi = math.inf, -math.inf
    for _, _, m in blocks():
        pos = m[m >= ZERO_EXPOSURE]
        if pos.size:
            s = np.log10(pos)
            lo, hi = min(lo, float(s.min())), max(hi, float(s.max()))
    if lo > hi:  # every pixel is unexposed
        lo, hi = 0.0, 1.0
    elif lo == hi:
        lo, hi = lo - 0.5, hi + 0.5
    edges = np.linspace(lo, hi, bins + 1)

    counts = np.zeros(bins, dtype=np.int64)
    sum_v = np.zeros((bins, n_z))
    sum_s = np.zeros((bins, n_z))
    zero = 0
    for ev, es, m in blocks():
        exposed = m >= ZERO_EXPOSURE
        zero += int(np.count_nonzero(~exposed))
        s = np.log10(m[exposed])
        k = np.clip(((s - lo) / (hi - lo) * bins).astype(np.int64), 0, bins - 1)
        counts += np.bincount(k, minlength=bins)
        ev_x, es_x = ev[exposed], es[exposed]
        for j in range(n_z):
            sum_v[:, j] += np.bincount(k, weights=ev_x[:, j], minlength=bins)
            sum_s[:, j] += np.bincount(k, weights=es_x[:, j], minlength=bins)

    filled = counts > 0
    mean_v = np.empty((bins, n_z))
    mean_s = np.zeros((bins, n_z))
    mean_v[filled] = sum_v[filled] / counts[filled, None]
    mean_s[filled] = sum_s[filled] / counts[filled, None]
    # Empty bins: the exposure whose key is the bin centre, all on the volume descriptor.
    w_c = np.asarray(classes.number_fraction, dtype=np.float64)
    wv = float(np.sum(w_c * classes.volume_um3))
    centre = 0.5 * (edges[:-1] + edges[1:])
    mean_v[~filled] = (10.0 ** centre[~filled] / wv)[:, None]
    return ExposureHistogram(edges, counts, mean_v, mean_s, zero)


def _columns(hist: ExposureHistogram) -> Members:
    """The 65 prepass columns: the bins, then the zero column (no exposure)."""
    n_z = hist.mean_volume.shape[1]
    ev = np.concatenate([hist.mean_volume, np.zeros((1, n_z))])
    es = np.concatenate([hist.mean_surface, np.zeros((1, n_z))])
    return Members(exposure_volume=ev, exposure_surface=es, layout="columns")


def _weights(hist: ExposureHistogram, image_area_fraction: float) -> np.ndarray:
    f = _fraction(image_area_fraction)
    total = float(hist.total)
    w = np.concatenate([hist.counts / total, [hist.zero_count / total]]) * f
    w[-1] += 1.0 - f
    return w


def _fraction(value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"image_area_fraction must be a number, got {value!r}")  # noqa: TRY004
    f = float(value)
    if not (0.0 < f <= 1.0):
        raise ValueError(f"image_area_fraction must be in (0, 1], got {value!r}")
    return f


def coupled_columns(
    film: FilmColumnModel,
    solution: WorkingSolution,
    kin: Kinetics,
    schedule: Schedule,
    members: Members,
    weights: np.ndarray,
    tank: TankSpec,
    presoaked: bool,
    *,
    method: str = "fixed",
    record_cells: int = 0,
) -> EnsembleResult:
    """Columns exactly coupled to ``tank`` (prepass and strips, spec D1 3.7)."""
    if not isinstance(tank, TankSpec):
        raise TypeError(f"tank must be a TankSpec, got {type(tank).__name__}")
    w = np.asarray(weights, dtype=np.float64)
    return run_ensemble(
        film,
        solution,
        kin,
        schedule,
        members,
        n_cells=int(w.size),
        cell_weights=w,
        tank=tank,
        presoaked=bool(presoaked),
        method=method,  # type: ignore[arg-type]
        record_cells=record_cells,
    )


def _given_columns(film, solution, kin, schedule, members, bulk, presoaked) -> EnsembleResult:
    n = members.exposure_volume.shape[0]
    return run_ensemble(
        film,
        solution,
        kin,
        schedule,
        members,
        n_cells=n,
        cell_weights=np.ones(n),
        bulk=bulk,
        presoaked=bool(presoaked),
        method="fixed",
    )


def zero_column(
    film: FilmColumnModel,
    solution: WorkingSolution,
    kin: Kinetics,
    schedule: Schedule,
    bulk: BulkTrajectory,
    presoaked: bool,
) -> EnsembleResult:
    """Unexposed film developed against ``bulk``: ``flux_mean`` is ``Jbar_zero``."""
    n_z = film.grid.n_z
    zero = np.zeros((1, n_z))
    members = Members(exposure_volume=zero, exposure_surface=zero, layout="columns")
    return _given_columns(film, solution, kin, schedule, members, bulk, presoaked)


def developed_silver(result: EnsembleResult) -> np.ndarray:
    """``M = sum_c silver_c + silver_physical`` per member, g/m2 (spec D1 3.8)."""
    return np.sum(result.silver, axis=1) + result.silver_physical


def trajectory_array(bulk: BulkTrajectory) -> np.ndarray:
    """``(S+1, K_c + 1)``: totals then TOTH, the layout of ``flux_mean``."""
    return np.concatenate(
        [np.asarray(bulk.totals, np.float64), np.asarray(bulk.toth, np.float64)[:, None]], axis=1
    )


def _trajectory(times: np.ndarray, values: np.ndarray) -> BulkTrajectory:
    return BulkTrajectory(
        times_s=_readonly(np.array(times, dtype=np.float64)),
        totals=_readonly(np.array(values[:, :-1])),
        toth=_readonly(np.array(values[:, -1])),
    )


def tank_prepass(
    hist: ExposureHistogram,
    film: FilmColumnModel,
    solution: WorkingSolution,
    kin: Kinetics,
    schedule: Schedule,
    tank: TankSpec,
    image_area_fraction: float,
    presoaked: bool,
) -> tuple[BulkTrajectory, np.ndarray]:
    """Bulk trajectory and ``Jbar_zero`` ``(S, K_c + 1)`` from the histogram columns (spec
    D1 3.7 step 2)."""
    if not isinstance(hist, ExposureHistogram):
        raise TypeError(f"hist must be an ExposureHistogram, got {type(hist).__name__}")
    if hist.mean_volume.shape[1] != film.grid.n_z:
        raise ValueError(
            f"histogram has n_z = {hist.mean_volume.shape[1]}, the film model {film.grid.n_z}"
        )
    weights = _weights(hist, image_area_fraction)
    run = coupled_columns(film, solution, kin, schedule, _columns(hist), weights, tank, presoaked)
    zero = zero_column(film, solution, kin, schedule, run.bulk, presoaked)
    return run.bulk, np.array(zero.flux_mean)


@dataclass(frozen=True)
class PrepassCheck:
    """Result of :func:`prepass_check` (spec D1 3.7 step 4, 10): ``delta`` the acceptance
    quantity, ``depletion_change`` and ``ph_change`` diagnostics, ``trajectory`` ``c'_bulk``
    (the next trajectory if the check fails)."""

    delta: float
    depletion_change: float
    ph_change: float
    trajectory: BulkTrajectory

    @property
    def passed(self) -> bool:
        """``delta <= epsilon_tank``."""
        return tk.prepass_passed(self.delta)


def bulk_ph(solution: WorkingSolution, values: np.ndarray) -> np.ndarray:
    """Bulk pH along a trajectory ``(S+1, K_c + 1)`` (totals then TOTH)."""
    v = np.asarray(values, dtype=np.float64)
    spec = solution.speciation
    idx = {cid: i for i, cid in enumerate(solution.component_ids)}
    cpos = positive_part(v[:, :-1])
    totals = np.stack([cpos[:, idx[f.component]] for f in spec.families])
    ph0 = np.full(v.shape[0], float(solution.ph))
    return solve_ph(spec, totals, v[:, -1], ph0=ph0)


def prepass_check(
    hist: ExposureHistogram,
    film: FilmColumnModel,
    solution: WorkingSolution,
    kin: Kinetics,
    schedule: Schedule,
    tank: TankSpec,
    image_area_fraction: float,
    presoaked: bool,
    current: BulkTrajectory,
    block_flux: np.ndarray,
    zero_flux: np.ndarray,
) -> PrepassCheck:
    """The output check of spec D1 3.7 step 4.

    ``block_flux`` ``(S, K_c + 1)`` is the area mean over all inner cells of the frame,
    ``zero_flux`` ``Jbar_zero`` under ``current``. The tank is re-integrated from their
    combination; ``delta`` compares the 65 columns' developed silver under the new and the
    current trajectory.
    """
    if not isinstance(current, BulkTrajectory):
        raise TypeError(f"current must be a BulkTrajectory, got {type(current).__name__}")
    f = _fraction(image_area_fraction)
    c = trajectory_array(current)
    times = np.asarray(current.times_s, dtype=np.float64)
    flux = tk.combine_flux(block_flux, zero_flux, f)
    coef = tk.tank_coefficient(film_area_um2=tank.film_area_um2, volume_l=tank.volume_l)
    c_new = tk.integrate_bulk(c[0], flux, times, coef)
    new = _trajectory(times, c_new)
    columns = _columns(hist)
    m_cur = developed_silver(
        _given_columns(film, solution, kin, schedule, columns, current, presoaked)
    )
    m_new = developed_silver(_given_columns(film, solution, kin, schedule, columns, new, presoaked))
    delta = tk.prepass_delta(np.maximum(m_new, 0.0), np.maximum(m_cur, 0.0))
    depletion = tk.depletion_change(c_new[:, :-1], c[:, :-1])
    ph = float(np.max(np.abs(bulk_ph(solution, c_new) - bulk_ph(solution, c))))
    return PrepassCheck(delta, depletion, ph, new)
