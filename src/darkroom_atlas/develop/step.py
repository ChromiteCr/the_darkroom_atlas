"""Presoak and development steps of the ``bw-mechanistic`` model, and the strip path (spec D1
2.0.3, 2.0.4, 2.4.9, 2.4.10, 3.6-3.8, 4.5).

``develop`` (``interaction_range_um=None``: the pipeline hands it the whole frame):

1. resolve the working solution, kinetics, film column model and the fixed step table of the
   run's profile (spec D1 3.1);
2. histogram the exposure and run the tank prepass (:mod:`darkroom_atlas.develop.prepass`);
3. develop canonical blocks of ``BLOCK_CELLS`` chemistry cells of 24 um, anchored at the
   frame origin, with a halo of ``ceil(halo_um / 24 um)`` cells (spec D1 3.6), each against
   the prescribed bulk trajectory, through :func:`~darkroom_atlas.blocks.run_blocks`; outside
   the frame is unexposed film. Workers write their owned pixels into the outputs and return
   only the block's flux sum;
4. check the prepass (spec D1 3.7) and recompute all blocks in place up to twice;
5. write ``emulsion/silver``, ``initiated``, ``silver_m2`` ``(H, W, K)``,
   ``silver_physical``, ``halide``, ``growth_rate`` ``(H, W)``, drop (and release) the two
   exposure fields and set ``meta["develop"]`` (spec D1 2.4.10).

Every block is a pure function of its padded input, the trajectory and the parameters, and
cross-block sums run in row-major block order, so the result does not depend on
``RunContext.tile_px`` or the number of workers.

Choices the spec leaves open: a develop step needs a temperature and a tank with
``film_area_cm2``; ``tank.rolls_before > 0`` (spec D1 5.6) is not supported yet; the
reference profile integrates adaptively (spec D1 3.1) on the calibration step table, the
other profiles use their fixed tables.
"""

from __future__ import annotations

import math
import shutil
import tempfile
import uuid
import warnings
from collections.abc import Sequence
from dataclasses import dataclass, replace
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from darkroom_atlas._frozen import freeze_mapping
from darkroom_atlas.blocks import (
    Block,
    block_grid,
    cell_count,
    overlapping_pixels,
    owned_pixels,
    padded_window,
    run_blocks,
    workers_for_budget,
)
from darkroom_atlas.chemistry.kinetics import Kinetics, resolve_kinetics
from darkroom_atlas.chemistry.solution import WorkingSolution, working_solution
from darkroom_atlas.data.models import BW_MECHANISTIC, StepKind
from darkroom_atlas.develop import lateral as lat
from darkroom_atlas.develop import prepass as pp
from darkroom_atlas.develop import tank as tk
from darkroom_atlas.develop.ensemble import (
    CHEM_CELL_UM,
    BulkTrajectory,
    EnsembleResult,
    FilmColumnModel,
    Lateral,
    Members,
    TankSpec,
    diffusion_factor,
    film_column_model,
    run_ensemble,
)
from darkroom_atlas.develop.schedule import Schedule, agitation_periods, step_schedule
from darkroom_atlas.fieldstore import FieldDescriptor, FieldStore, descriptor_of, open_field
from darkroom_atlas.registry import register_step
from darkroom_atlas.state import FieldSpec, FrameAssembler

if TYPE_CHECKING:
    from darkroom_atlas.data.models import Bath, Material, ProcessStep
    from darkroom_atlas.pipeline import RunContext, StepContext
    from darkroom_atlas.state import EmulsionState

__all__ = [
    "BLOCK_CELLS",
    "FORMAT_IMAGE_AREA_FRACTION",
    "MODEL",
    "PROFILE_SCHEDULES",
    "develop_bw",
    "develop_frame",
    "develop_strip",
    "halo_um",
    "presoak_bw",
]

MODEL = BW_MECHANISTIC
"""Material model the steps are registered for (spec D1 2.0.4)."""

BLOCK_CELLS = 256
"""Inner side of a canonical development block, chemistry cells (B_c, spec D1 3.6)."""

FORMAT_IMAGE_AREA_FRACTION = {"135": (36.0 * 24.0) / (38.0 * 35.0), "120": 0.80}
"""Default ``tank.image_area_fraction`` by film format (spec D1 3.7, assumed)."""

PROFILE_SCHEDULES: dict[str, tuple[float, float, float]] = {
    "image": (0.1, 1.6, 5.0),
    "calibration": (0.05, 1.4, 3.0),
    "reference": (0.05, 1.4, 3.0),
}
"""Step table ``(h0 s, ratio, h_max s)`` per profile (spec D1 3.1); the reference profile
subdivides the calibration table adaptively."""

_HALO_SIGMAS = 3.0  # 3 sigma of lateral liquid diffusion (spec D1 3.6)
_HALO_SCREEN = 4.6  # screening-length cap in units of delta_max (spec D1 3.6)
_PERIODIC = ("inversion", "lift")
_CONTINUOUS = ("rotary", "continuous", "tray")
_OUT_K = ("emulsion/silver", "emulsion/initiated", "emulsion/silver_m2")
_OUT_1 = ("emulsion/silver_physical", "emulsion/halide", "emulsion/growth_rate")
_EXPOSURE = ("emulsion/exposure_volume", "emulsion/exposure_surface")


# ---------------------------------------------------------------------------------------
# Conditions
# ---------------------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class _Conditions:
    film: FilmColumnModel
    solution: WorkingSolution
    kin: Kinetics
    schedule: Schedule
    tank: TankSpec
    image_area_fraction: float
    time_s: float
    temperature_c: float
    method: str


def _conditions(
    material: Material,
    film_format: str,
    exposure_time_s: float,
    profile: str,
    step: ProcessStep,
    bath: Bath | None,
    library,
    record_s: Sequence[float] = (),
) -> _Conditions:
    if step.kind is not StepKind.DEVELOP:
        raise ValueError(f"expected a develop step, got {step.kind.value!r}")
    if bath is None:
        raise ValueError("a develop step of the bw-mechanistic model needs a bath")
    time_s = step.time_s
    if time_s is None:
        raise ValueError("a develop step needs a time")
    if step.temperature is None:
        raise ValueError("a develop step of the bw-mechanistic model needs a temperature")
    temperature_c = float(step.temperature.to("C"))
    tank = step.tank
    if tank is None or tank.film_area_cm2 is None:
        raise ValueError(
            "a develop step of the bw-mechanistic model needs a tank with volume_ml and "
            "film_area_cm2 (spec D1 2.4.9, 3.7)"
        )
    if tank.rolls_before > 0:
        raise NotImplementedError("tank.rolls_before > 0 (spec D1 5.6) is not supported yet")
    if tank.image_area_fraction is not None:
        f_img = float(tank.image_area_fraction.value)
    elif film_format in FORMAT_IMAGE_AREA_FRACTION:
        f_img = FORMAT_IMAGE_AREA_FRACTION[film_format]
    else:
        raise ValueError(
            f"no default image_area_fraction for film format {film_format!r}; set "
            "tank.image_area_fraction"
        )
    if profile not in PROFILE_SCHEDULES:
        raise ValueError(f"unknown profile {profile!r}")
    species = library.species()
    kfile = library.kinetics()
    solution = working_solution(
        bath, species, kfile, dilution=step.dilution, temperature_c=temperature_c
    )
    kin = resolve_kinetics(kfile, solution, species)
    film = film_column_model(
        material, film_format, exposure_time_s, kin, solution, profile, species=species
    )
    h0, ratio, h_max = PROFILE_SCHEDULES[profile]
    periods = agitation_periods(step.agitation, time_s)
    schedule = step_schedule(
        periods, time_s, h0_s=h0, ratio=ratio, h_max_s=h_max, record_s=tuple(record_s)
    )
    return _Conditions(
        film=film,
        solution=solution,
        kin=kin,
        schedule=schedule,
        tank=TankSpec(volume_l=tank.volume_ml * 1e-3, film_area_um2=tank.film_area_cm2 * 1e8),
        image_area_fraction=f_img,
        time_s=float(time_s),
        temperature_c=temperature_c,
        method="adaptive" if profile == "reference" else "fixed",
    )


def halo_um(
    film: FilmColumnModel, solution: WorkingSolution, kin: Kinetics, step: ProcessStep
) -> float:
    """Halo of the canonical development blocks (spec D1 3.6), um.

    ``min(max_k 3 sqrt(2 D_k t_q,max / (1 - rho_k)), 4.6 delta_max)`` with the liquid
    ``D_w(T)`` (no ionic-strength factor) of every mobile component and of mobile TOTH at the
    bath's initial speciation; ``t_q,max`` the longest still interval (the development time
    under continuous agitation), ``rho_k`` from the shortest periodic agitation period.
    """
    time_s = step.time_s
    if time_s is None:
        raise ValueError("a develop step needs a time")
    t = float(solution.temperature_c)
    # D(T, I) -> D_w(T): drop the 1 / (1 + b_eta I) factor (spec D1 3.6, an upper bound).
    unscale = diffusion_factor(kin.transport, t, 0.0) / diffusion_factor(
        kin.transport, t, solution.ionic_strength
    )
    d_comp = np.asarray(film.diffusion_um2_s, dtype=np.float64) * unscale
    d_h, d_oh = (d * unscale for d in film.water_ion_diffusion_um2_s)
    # Lateral D of mobile TOTH (spec D1 3.5) in the initial bath.
    spec = solution.speciation
    idx = {cid: i for i, cid in enumerate(solution.component_ids)}
    fam_comp = [idx[f.component] for f in spec.families]
    mobile = [d_comp[c] > 0.0 for c in fam_comp]
    totals = np.array([[solution.totals[c]] for c in fam_comp], dtype=np.float64)
    beta = lat.carrier_capacities(spec, totals, np.array([solution.ph]), mobile)[:, 0]
    d_carrier = np.array(
        [d_comp[c] for c, m in zip(fam_comp, mobile, strict=True) if m] + [d_h, d_oh]
    )
    d_toth = float(np.sum(d_carrier * beta) / np.sum(beta))
    d_all = np.concatenate([d_comp[d_comp > 0.0], [d_toth]])

    periods = agitation_periods(step.agitation, time_s)
    method = None if step.agitation is None else step.agitation.method
    if method in _CONTINUOUS:
        t_q = float(time_s)
    else:
        edges = [0.0]
        for p in periods:
            edges.extend([p.start_s, p.stop_s])
        edges.append(float(time_s))
        t_q = max(edges[i + 1] - edges[i] for i in range(0, len(edges), 2))
    delta_agit, delta_max = film.grid.delta_agit_um, film.grid.delta_max_um
    if method in _PERIODIC and periods:
        t_a = min(p.stop_s - p.start_s for p in periods)
        rho = np.exp(-t_a * d_all * (math.pi / (2.0 * delta_agit)) ** 2)
    else:
        rho = np.zeros_like(d_all)
    spread = float(np.max(_HALO_SIGMAS * np.sqrt(2.0 * d_all * t_q / (1.0 - rho))))
    return min(spread, _HALO_SCREEN * delta_max)


# ---------------------------------------------------------------------------------------
# Canonical blocks
# ---------------------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class _BlockSetup:
    """Everything a block job needs besides its block (pickled into every task)."""

    film: FilmColumnModel
    solution: WorkingSolution
    kin: Kinetics
    schedule: Schedule
    bulk: BulkTrajectory
    presoaked: bool
    pitch_um: float
    frame_hw: tuple[int, int]
    halo_cells: int


@dataclass(frozen=True)
class _Summary:
    pixels: tuple[slice, slice]
    flux_sum: np.ndarray
    cells: int
    min_scaled: float


def _block_values(block: Block, setup: _BlockSetup, ev_src, es_src):
    """Develop one canonical block; returns (owned pixel region, values, summary)."""
    p = setup.pitch_um
    bp = block.padded_shape
    starts, lengths, origin = [], [], []
    for c0, n in zip(block.padded_origin, bp, strict=True):
        j0, _ = overlapping_pixels(slice(c0, c0 + n), p, CHEM_CELL_UM)
        starts.append(j0)
        lengths.append(math.ceil(n * CHEM_CELL_UM / p) + 1)  # fixed for every block
        origin.append(j0 * p - c0 * CHEM_CELL_UM)
    window = Block(
        index=block.index,
        inner=block.inner,
        padded_origin=(starts[0], starts[1]),
        padded_shape=(lengths[0], lengths[1]),
    )
    ev = padded_window(ev_src, window, 0.0)
    es = padded_window(es_src, window, 0.0)
    weights = np.zeros(bp)
    weights[block.inner_in_padded] = 1.0
    n_inner = int(np.count_nonzero(weights))
    res = run_ensemble(
        setup.film,
        setup.solution,
        setup.kin,
        setup.schedule,
        Members(ev, es, "raster", pixel_um=p, origin_um=(origin[0], origin[1])),
        n_cells=bp[0] * bp[1],
        cell_weights=weights.ravel(),
        lateral=Lateral(shape_cells=bp, spacing_um=CHEM_CELL_UM, boundary="neumann"),
        bulk=setup.bulk,
        presoaked=setup.presoaked,
        method="fixed",
    )
    # Pixels the engine returned (its own owned-pixel rule, window coordinates) and the
    # pixels this block owns in the frame (spec D1 3.6).
    local, own = [], []
    for axis in range(2):
        inner = block.inner_in_padded[axis]
        a = max(math.ceil((inner.start * CHEM_CELL_UM - origin[axis]) / p - 0.5), 0)
        b = min(math.ceil((inner.stop * CHEM_CELL_UM - origin[axis]) / p - 0.5), lengths[axis])
        frame = owned_pixels(block.inner[axis], p, setup.frame_hw[axis], CHEM_CELL_UM)
        lo, hi = frame.start - starts[axis] - a, frame.stop - starts[axis] - a
        if frame.stop > frame.start and (lo < 0 or hi > b - a):
            raise RuntimeError(f"block {block.index}: engine pixels do not cover the owned ones")
        local.append((b - a, slice(lo, max(lo, hi))))
        own.append(frame)
    (nr, rows), (nc, cols) = local
    k = res.silver.shape[1]
    values = {}
    for key, arr in zip(
        (*_OUT_K, *_OUT_1),
        (
            res.silver,
            res.initiated,
            res.silver_m2,
            res.silver_physical,
            res.halide,
            res.growth_rate,
        ),
        strict=True,
    ):
        shaped = arr.reshape((nr, nc, k) if arr.ndim == 2 else (nr, nc))
        values[key] = shaped[rows, cols]
    summary = _Summary(
        pixels=(own[0], own[1]),
        flux_sum=np.asarray(res.flux_mean) * float(n_inner),
        cells=n_inner,
        min_scaled=float(res.min_scaled_concentration),
    )
    return values, summary


def _worker_job(
    block: Block,
    *,
    setup: _BlockSetup,
    inputs: tuple[FieldDescriptor, FieldDescriptor],
    outputs: dict[str, FieldDescriptor],
) -> _Summary:
    """Worker side: develop one block, write its owned pixels, return its summary."""
    values, summary = _block_values(
        block, setup, open_field(inputs[0], "r"), open_field(inputs[1], "r")
    )
    if summary.pixels[0].stop > summary.pixels[0].start and summary.pixels[1].stop > (
        summary.pixels[1].start
    ):
        for key, v in values.items():
            out = open_field(outputs[key], "r+")
            out[summary.pixels] = v
            out.flush()
            del out
    return summary


def _empty(region: tuple[slice, slice]) -> bool:
    return region[0].stop <= region[0].start or region[1].stop <= region[1].start


class _Runner:
    """Runs all blocks of one frame against a trajectory, in- or out-of-process."""

    def __init__(self, ev, es, blocks, setup, workers, assembler, geometry, specs, store):
        self.ev, self.es, self.blocks, self.setup = ev, es, blocks, setup
        self.workers = min(workers, len(blocks))
        self.assembler = assembler
        self.geometry, self.specs, self.store = geometry, specs, store
        self.target = assembler
        self.transient: FieldStore | None = None
        self.scratch: Path | None = None
        self.owned_scratch = False
        self.inputs = None
        self.outputs = None
        if self.workers > 1:
            self._prepare_workers()

    def _prepare_workers(self) -> None:
        outputs = self.assembler.descriptors()
        mapped = set(outputs) == {s.key for s in self.specs}
        if self.store is not None and self.store.root is not None:
            self.store.root.mkdir(parents=True, exist_ok=True)
            self.scratch = self.store.root
        else:
            self.scratch = Path(tempfile.mkdtemp(prefix="atlas-develop-"))
            self.owned_scratch = True
        self.transient = FieldStore(self.scratch, ram_threshold_bytes=0)
        prefix = f"develop-{uuid.uuid4().hex[:8]}:"
        inputs = []
        for i, src in enumerate((self.ev, self.es)):
            try:
                inputs.append(descriptor_of(src))
            except ValueError:
                name = f"{prefix}in{i}"
                copy = self.transient.allocate(name, src.shape)
                copy[...] = src
                del copy
                self.transient.freeze(name)
                inputs.append(self.transient.descriptor(name))
        self.inputs = tuple(inputs)
        if not mapped:
            self.target = FrameAssembler(
                self.geometry, self.specs, self.transient, name_prefix=prefix
            )
            outputs = self.target.descriptors()
        self.outputs = outputs

    def run(self, bulk: BulkTrajectory, first: bool) -> list[_Summary]:
        setup = replace(self.setup, bulk=bulk)
        if not first:
            self.target.reset()
        if self.workers <= 1:
            out = []
            for block in self.blocks:
                values, summary = _block_values(block, setup, self.ev, self.es)
                if not _empty(summary.pixels):
                    for key, v in values.items():
                        self.target.write(key, summary.pixels, v)
                out.append(summary)
            return out
        job = partial(_worker_job, setup=setup, inputs=self.inputs, outputs=self.outputs)
        out = run_blocks(job, self.blocks, workers=self.workers)
        for summary in out:
            if not _empty(summary.pixels):
                self.target.mark(summary.pixels)
        return out

    def finish(self) -> dict[str, np.ndarray]:
        try:
            fields = self.target.finish()
            if self.target is not self.assembler:
                h, w = self.geometry.shape_hw
                for key, value in fields.items():
                    self.assembler.write(key, (slice(0, h), slice(0, w)), np.array(value))
                fields = self.assembler.finish()
            return fields
        finally:
            self.close()

    def close(self) -> None:
        if self.transient is not None:
            self.transient.close()
            self.transient = None
        if self.owned_scratch and self.scratch is not None:
            shutil.rmtree(self.scratch, ignore_errors=True)
            self.owned_scratch = False


# ---------------------------------------------------------------------------------------
# Steps
# ---------------------------------------------------------------------------------------


def _exposure_inputs(state: EmulsionState) -> tuple[np.ndarray, np.ndarray, dict]:
    missing = [k for k in _EXPOSURE if k not in state.fields]
    if missing:
        raise ValueError(
            f"the bw-mechanistic develop step needs the exposure fields; missing: "
            f"{', '.join(missing)} (expose with the bw-mechanistic model first)"
        )
    meta = state.meta.get("exposure")
    if meta is None or not all(k in meta for k in ("time_s", "film_format", "profile", "n_z")):
        raise ValueError(
            'the develop step needs meta["exposure"] with time_s, film_format, profile, n_z'
        )
    ev, es = (state.fields[k] for k in _EXPOSURE)
    if ev.ndim != 3 or es.shape != ev.shape or ev.shape[2] != int(meta["n_z"]):
        raise ValueError(
            f"exposure fields must share a shape (H, W, n_z = {meta['n_z']}), got "
            f"{ev.shape} and {es.shape}"
        )
    return ev, es, dict(meta)


def _presoaked(state: EmulsionState) -> bool:
    field = state.fields.get("frame/presoaked")
    if field is None:
        return False  # no presoak field: a dry film (spec D1 2.4.10)
    flat = np.asarray(field).reshape(field.shape[0], field.shape[1], -1)
    first = bool(flat[0, 0, 0])
    if not np.all(np.asarray(field) == first):
        raise ValueError("frame/presoaked must be uniform over the frame")
    return first


def _bulk_dict(solution: WorkingSolution, row: np.ndarray) -> dict[str, float]:
    out = {cid: float(v) for cid, v in zip(solution.component_ids, row[:-1], strict=True)}
    out["toth"] = float(row[-1])
    return out


def develop_frame(
    state: EmulsionState,
    step: ProcessStep,
    ctx: StepContext,
    *,
    block_cells: int = BLOCK_CELLS,
    halo_cells: int | None = None,
) -> EmulsionState:
    """The develop step with an explicit block size and halo (``None``: spec D1 3.6).

    :func:`develop_bw` is this with the spec values; tests use small blocks.
    """
    if state.material.model != MODEL:
        raise ValueError(f"material {state.material.id!r} is not a {MODEL!r} material")
    if state.geometry.origin_px != (0, 0) or state.geometry.frame_shape_hw != (
        state.geometry.shape_hw
    ):
        raise ValueError("the develop step needs the whole frame, not a window")
    for name, value, minimum in (("block_cells", block_cells, 1), ("halo_cells", halo_cells, 0)):
        if value is not None and (
            isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < minimum
        ):
            raise ValueError(f"{name} must be an int >= {minimum}, got {value!r}")
    ev, es, exp_meta = _exposure_inputs(state)
    profile = str(exp_meta["profile"])
    cond = _conditions(
        state.material,
        str(exp_meta["film_format"]),
        float(exp_meta["time_s"]),
        profile,
        step,
        ctx.bath,
        ctx.run.library,
    )
    film, solution, kin, schedule = cond.film, cond.solution, cond.kin, cond.schedule
    if film.grid.n_z != ev.shape[2]:
        raise ValueError(f"exposure has n_z = {ev.shape[2]}, profile {profile!r} {film.grid.n_z}")
    presoaked = _presoaked(state)
    halo = halo_um(film, solution, kin, step)
    h_c = math.ceil(halo / CHEM_CELL_UM) if halo_cells is None else int(halo_cells)

    hist = pp.exposure_histogram(ev, es, film.classes)
    bulk, j_zero = pp.tank_prepass(
        hist, film, solution, kin, schedule, cond.tank, cond.image_area_fraction, presoaked
    )

    geometry = state.geometry
    pitch = geometry.pixel_pitch_um
    h_px, w_px = geometry.shape_hw
    cells = (cell_count(h_px, pitch, CHEM_CELL_UM), cell_count(w_px, pitch, CHEM_CELL_UM))
    blocks = block_grid(cells, int(block_cells), h_c)
    k = film.classes.count
    specs = [FieldSpec(key, (k,)) for key in _OUT_K] + [FieldSpec(key) for key in _OUT_1]
    store = ctx.run.store
    assembler = FrameAssembler(geometry, specs, store, name_prefix=f"s{ctx.step_index}:")
    setup = _BlockSetup(
        film=film,
        solution=solution,
        kin=kin,
        schedule=schedule,
        bulk=bulk,
        presoaked=presoaked,
        pitch_um=float(pitch),
        frame_hw=(h_px, w_px),
        halo_cells=h_c,
    )
    workers = workers_for_budget(ctx.run.memory_budget_gb, ctx.run.workers)  # spec D1 4.6
    runner = _Runner(ev, es, blocks, setup, workers, assembler, geometry, specs, store)
    deltas: list[float] = []
    try:
        for n_pass in range(1 + tk.MAX_EXTRA_PASSES):
            summaries = runner.run(bulk, first=n_pass == 0)
            flux = np.zeros_like(summaries[0].flux_sum)
            for s in summaries:  # row-major block order (spec D1 3.6)
                flux = flux + s.flux_sum
            flux = flux / float(sum(s.cells for s in summaries))
            check = pp.prepass_check(
                hist,
                film,
                solution,
                kin,
                schedule,
                cond.tank,
                cond.image_area_fraction,
                presoaked,
                bulk,
                flux,
                j_zero,
            )
            deltas.append(check.delta)
            if check.passed or n_pass == tk.MAX_EXTRA_PASSES:
                break
            bulk = check.trajectory
            j_zero = np.array(
                pp.zero_column(film, solution, kin, schedule, bulk, presoaked).flux_mean
            )
        fields = runner.finish()
    finally:
        runner.close()
    if not check.passed:
        warnings.warn(
            f"tank prepass check failed after {len(deltas)} passes: delta = {deltas[-1]:.3g} "
            f"> {tk.TANK_TOLERANCE} (spec D1 3.7)",
            RuntimeWarning,
            stacklevel=2,
        )
    contraction = tk.contraction(deltas)
    if any(c > tk.CONTRACTION_WARN for c in contraction):
        warnings.warn(
            f"tank prepass contraction estimate {max(contraction):.3g} > "
            f"{tk.CONTRACTION_WARN} (spec D1 3.7)",
            RuntimeWarning,
            stacklevel=2,
        )

    zero = pp.zero_column(film, solution, kin, schedule, bulk, presoaked)
    traj = pp.trajectory_array(bulk)
    prepass_meta: dict[str, Any] = {
        "bins": hist.bins,
        "passes": len(deltas),
        "delta": tuple(deltas),
        "depletion_change": check.depletion_change,
        "ph_change": check.ph_change,
    }
    if contraction:
        prepass_meta["contraction"] = contraction
    agitation = None if step.agitation is None else step.agitation.model_dump()
    meta = {
        "bath": ctx.bath.id if ctx.bath is not None else None,
        "dilution": step.dilution,
        "temperature_c": cond.temperature_c,
        "time_s": cond.time_s,
        "agitation": agitation,
        "tank": {
            "volume_l": cond.tank.volume_l,
            "film_area_um2": cond.tank.film_area_um2,
            "image_area_fraction": cond.image_area_fraction,
        },
        "profile": profile,
        "chem_cell_um": CHEM_CELL_UM,
        "halo_um": float(halo),
        "block_cells": int(block_cells),
        "n_steps": schedule.n_steps,
        "min_scaled_concentration": min(s.min_scaled for s in summaries),
        "prepass": prepass_meta,
        "bulk_initial": _bulk_dict(solution, traj[0]),
        "bulk_final": _bulk_dict(solution, traj[-1]),
        "zero_exposure": {
            "silver": tuple(float(v) for v in zero.silver[0]),
            "initiated": tuple(float(v) for v in zero.initiated[0]),
            "silver_m2": tuple(float(v) for v in zero.silver_m2[0]),
            "silver_physical": float(zero.silver_physical[0]),
            "halide": float(zero.halide[0]),
        },
    }
    out = state.drop(*_EXPOSURE).with_fields(fields).with_meta(develop=freeze_mapping(meta))
    if store is not None:
        for field in (ev, es):
            try:
                store.release(store.name_of(field))
            except KeyError:
                pass
    return out


@register_step(MODEL, StepKind.DEVELOP, interaction_range_um=None)
def develop_bw(state: EmulsionState, step: ProcessStep, ctx: StepContext) -> EmulsionState:
    """Develop a ``bw-mechanistic`` frame (spec D1 2.4, 3.6-3.8)."""
    return develop_frame(state, step, ctx, block_cells=BLOCK_CELLS)


@register_step(MODEL, StepKind.PRESOAK, interaction_range_um=0.0)
def presoak_bw(state: EmulsionState, step: ProcessStep, ctx: StepContext) -> EmulsionState:
    """Presoak: the broadcast field ``frame/presoaked = True``, nothing else (spec D1 2.4.10)."""
    one = np.ones((1, 1), dtype=bool)
    one.flags.writeable = False
    return state.with_fields({"frame/presoaked": np.broadcast_to(one, state.geometry.shape_hw)})


def develop_strip(
    material: Material,
    film_format: str,
    exposures: tuple[np.ndarray, np.ndarray],
    exposure_time_s: float,
    step: ProcessStep,
    ctx: RunContext,
    *,
    record_s: Sequence[float] = (),
) -> EnsembleResult:
    """Develop a strip: one column per step, exactly coupled to the tank (spec D1 3.7).

    ``exposures`` ``(e_V, e_S)``, each ``(M, n_z)`` for ``ctx.profile``. The strip covers
    ``image_area_fraction`` of the tank's film area in equal steps, the rest is unexposed
    film (one more column, not returned). The physics code is the image path's.
    """
    try:
        ev, es = (np.asarray(a, dtype=np.float64) for a in exposures)
    except (TypeError, ValueError):
        raise ValueError("exposures must be a pair (e_V, e_S) of (M, n_z) arrays") from None
    if ev.ndim != 2 or ev.shape != es.shape or ev.shape[0] < 1:
        raise ValueError(f"exposures must share a shape (M, n_z), got {ev.shape} and {es.shape}")
    bath = None if step.bath is None else ctx.library.bath(step.bath)
    cond = _conditions(
        material, film_format, exposure_time_s, ctx.profile, step, bath, ctx.library, record_s
    )
    n_z = cond.film.grid.n_z
    if ev.shape[1] != n_z:
        raise ValueError(f"exposures have n_z = {ev.shape[1]}, profile {ctx.profile!r} {n_z}")
    m = ev.shape[0]
    f = cond.image_area_fraction
    weights = np.full(m, f / m)
    if f < 1.0:
        ev = np.concatenate([ev, np.zeros((1, n_z))])
        es = np.concatenate([es, np.zeros((1, n_z))])
        weights = np.concatenate([weights, [1.0 - f]])
    res = pp.coupled_columns(
        cond.film,
        cond.solution,
        cond.kin,
        cond.schedule,
        Members(ev, es, "columns"),
        weights,
        cond.tank,
        False,
        method=cond.method,
    )
    return _first_members(res, m) if f < 1.0 else res


def _first_members(res: EnsembleResult, m: int) -> EnsembleResult:
    snaps = {t: _first_members(s, m) for t, s in res.snapshots.items()}
    return replace(
        res,
        silver=res.silver[:m],
        initiated=res.initiated[:m],
        silver_m2=res.silver_m2[:m],
        silver_physical=res.silver_physical[:m],
        halide=res.halide[:m],
        growth_rate=res.growth_rate[:m],
        snapshots=freeze_mapping(snaps) if snaps else res.snapshots,
    )
