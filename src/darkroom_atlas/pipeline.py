"""The three operations of the simulation: ``expose``, ``process`` and ``readout``.

* :func:`expose` turns light on the film (any :class:`~darkroom_atlas.light.LightField`)
  into an :class:`~darkroom_atlas.state.EmulsionState` with the exposure model that the
  material's ``model`` names.
* :func:`process` runs the steps of a process in order, each with the handler registered
  for ``(material model, step kind)``, tiling the frame when the run asks for it.
* :func:`readout` turns a state into a :class:`Readout` (an image, spectra or densities)
  with a named readout device.

The physics lives in the registered models (:mod:`darkroom_atlas.registry`); this module
only dispatches, checks what the models return and keeps the history. A spectral
transmittance read out of a processed negative can be applied to a light with
:meth:`Light.filtered <darkroom_atlas.light.Light.filtered>` and exposed onto the next
material, which is how printing reuses the same engine.

Segments and tiling (spec D1 section 4.5): :func:`segments` groups adjacent steps whose
handlers declare a finite ``interaction_range_um`` into one segment; a step with
``interaction_range_um=None`` is a segment of its own and always receives the whole frame.
When ``RunContext.tile_px`` is set, a segment of finite-range steps is run tile by tile:
each tile is cut from the segment's input with :meth:`EmulsionState.crop` (views, no
copies), grown by the sum over the segment's steps of
``ceil(interaction_range_um / pixel_pitch_um)`` pixels, passed through all of the segment's
steps in turn, and only the segment's final fields are written into the frame with a
:class:`~darkroom_atlas.state.FrameAssembler` (O(pixels) in all). A final field that is
still the input's own array in every tile (no step touched it) is passed through without
being copied. Without ``tile_px`` every step runs on the whole frame. Assembled fields are
allocated in the run's :class:`~darkroom_atlas.fieldstore.FieldStore`
(``RunContext.store``) under the store names ``f"s{index}:{key}"``, ``index`` being the
segment's last step (spec §4.3).

A tiled handler must return the same field names, each with the same dtype and trailing
dimensions, and the same ``meta`` for every tile (differing fields raise ``ValueError``):
the assembled state takes its ``meta`` from the first tile, and the handler must not touch
the history (the pipeline appends each step's record, also between the steps of a segment
inside a tile, so a later step sees the earlier ones as it would untiled).
"""

import math
import os
import uuid
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from numbers import Integral, Real
from pathlib import Path
from typing import Any, Literal, get_args

import numpy as np

from darkroom_atlas._frozen import freeze_mapping, readonly_array
from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import Bath, Material, Process, ProcessStep
from darkroom_atlas.fieldstore import FieldStore
from darkroom_atlas.light import LightField
from darkroom_atlas.registry import RegisteredStep, get_exposure, get_readout, get_step
from darkroom_atlas.spectral import DEFAULT_GRID, SpectralGrid
from darkroom_atlas.state import EmulsionState, FieldSpec, FrameAssembler, FrameGeometry, StepRecord
from darkroom_atlas.tiling import Tile, iter_tiles

__all__ = [
    "Readout",
    "ReadoutKind",
    "RunContext",
    "StepContext",
    "expose",
    "process",
    "readout",
    "segments",
    "unexposed_state",
]

ReadoutKind = Literal[
    "image_linear_srgb", "spectral_transmittance", "spectral_reflectance", "density"
]
_READOUT_KINDS: tuple[str, ...] = get_args(ReadoutKind)
_SPECTRAL_KINDS = frozenset({"spectral_transmittance", "spectral_reflectance"})

Profile = Literal["image", "calibration", "reference"]
_PROFILES: tuple[str, ...] = get_args(Profile)


# --------------------------------------------------------------------------------------
# Context and results
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class RunContext:
    """Settings shared by every operation of one run.

    ``library`` supplies data (spectral curves, baths), ``grid`` is the spectral grid every
    light and state must use, ``seed`` keys all random numbers, ``tile_px`` (``None``: no
    tiling) is the tile size for processing steps and ``long_edge_mm`` the length of the
    frame's long edge on the film.

    D1 additions (spec section 4.7): ``film_format`` names the film format (``None``: the
    model's default), ``profile`` the discretisation profile (section 3.1), ``workers`` the
    number of worker processes, ``memory_budget_gb`` the memory budget that bounds them
    (section 4.6) and ``scratch_dir`` the directory for memory-mapped fields (``None``:
    every field stays in memory). The runner owns ``scratch_dir``: it creates and removes it.

    ``store`` is the run's one :class:`~darkroom_atlas.fieldstore.FieldStore` (section
    4.3): created and closed by the runner, it takes no part in comparison, hashing or
    ``repr`` and is not recorded in the manifest or ``recipe_hash``. Steps reach it as
    ``StepContext.run.store``. Without a store, a ``process`` call with ``scratch_dir`` set
    maps its fields into a store of its own under ``scratch_dir``; with neither, every
    field stays in memory.
    """

    library: Library
    grid: SpectralGrid = DEFAULT_GRID
    seed: int = 0
    tile_px: int | None = None
    long_edge_mm: float = 36.0
    film_format: str | None = None
    profile: Profile = "image"
    workers: int = 1
    memory_budget_gb: float = 6.0
    scratch_dir: Path | None = None
    store: FieldStore | None = field(default=None, compare=False, hash=False, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.library, Library):
            raise TypeError(f"library must be a Library, got {type(self.library).__name__}")
        if not isinstance(self.grid, SpectralGrid):
            raise TypeError(f"grid must be a SpectralGrid, got {type(self.grid).__name__}")
        object.__setattr__(self, "seed", _int_at_least("seed", self.seed, 0))
        if self.tile_px is not None:
            object.__setattr__(self, "tile_px", _int_at_least("tile_px", self.tile_px, 1))
        long_edge_mm = _positive_finite("long_edge_mm", self.long_edge_mm)
        object.__setattr__(self, "long_edge_mm", long_edge_mm)
        if self.film_format is not None and (
            not isinstance(self.film_format, str) or not self.film_format
        ):
            raise ValueError(
                f"film_format must be None or a non-empty string, got {self.film_format!r}"
            )
        if self.profile not in _PROFILES:
            raise ValueError(f"profile must be one of {', '.join(_PROFILES)}, got {self.profile!r}")
        object.__setattr__(self, "workers", _int_at_least("workers", self.workers, 1))
        budget = _positive_finite("memory_budget_gb", self.memory_budget_gb)
        object.__setattr__(self, "memory_budget_gb", budget)
        if self.scratch_dir is not None:
            if not isinstance(self.scratch_dir, (str, os.PathLike)):
                raise TypeError(
                    f"scratch_dir must be None or a path, got {type(self.scratch_dir).__name__}"
                )
            object.__setattr__(self, "scratch_dir", Path(self.scratch_dir))
        if self.store is not None and not isinstance(self.store, FieldStore):
            raise TypeError(f"store must be None or a FieldStore, got {type(self.store).__name__}")


@dataclass(frozen=True)
class StepContext:
    """What a step handler gets besides the state: the run, the step's index and its bath.

    ``bath`` is the bath named by the step, loaded from ``run.library`` (``None`` for steps
    without one).
    """

    run: RunContext
    step_index: int
    bath: Bath | None

    def __post_init__(self) -> None:
        if not isinstance(self.run, RunContext):
            raise TypeError(f"run must be a RunContext, got {type(self.run).__name__}")
        object.__setattr__(self, "step_index", _int_at_least("step_index", self.step_index, 0))
        if self.bath is not None and not isinstance(self.bath, Bath):
            raise TypeError(f"bath must be a Bath or None, got {type(self.bath).__name__}")


@dataclass(frozen=True, eq=False)
class Readout:
    """What a readout device produced.

    ``kind`` says what ``data`` holds; its first two dimensions are the frame ``(H, W)``:

    * ``"image_linear_srgb"``: linear sRGB, shape ``(H, W, 3)``.
    * ``"spectral_transmittance"`` / ``"spectral_reflectance"``: per band, shape
      ``(H, W, B)``; ``grid`` is required and ``B == grid.n_bands``.
    * ``"density"``: optical densities (base 10), shape ``(H, W)`` or ``(H, W, C)``.

    ``grid``, when given, describes the last axis of ``data``. ``data`` is stored read-only
    and ``meta`` as a read-only mapping.
    """

    kind: ReadoutKind
    data: np.ndarray
    grid: SpectralGrid | None = None
    meta: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.kind not in _READOUT_KINDS:
            raise ValueError(
                f"unknown readout kind {self.kind!r}; expected one of: {', '.join(_READOUT_KINDS)}"
            )
        try:
            data = readonly_array(self.data)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"readout data is not an array: {exc}") from exc
        if data.dtype.kind not in "iuf":
            raise ValueError(f"readout data must be real numbers, got dtype {data.dtype}")
        if data.ndim < 2 or data.shape[0] < 1 or data.shape[1] < 1:
            raise ValueError(f"readout data must have shape (H, W, ...), got {data.shape}")
        if self.grid is not None and not isinstance(self.grid, SpectralGrid):
            raise TypeError(f"grid must be a SpectralGrid or None, got {type(self.grid).__name__}")
        if self.kind in _SPECTRAL_KINDS and self.grid is None:
            raise ValueError(f"a {self.kind!r} readout needs the spectral grid of its data")
        if self.grid is not None and (data.ndim != 3 or data.shape[2] != self.grid.n_bands):
            raise ValueError(
                f"readout data of shape {data.shape} does not match the grid's "
                f"{self.grid.n_bands} bands; expected (H, W, {self.grid.n_bands})"
            )
        if self.kind == "image_linear_srgb" and (data.ndim != 3 or data.shape[2] != 3):
            raise ValueError(f"an image readout must have shape (H, W, 3), got {data.shape}")
        object.__setattr__(self, "data", data)
        object.__setattr__(self, "meta", freeze_mapping(self.meta))


# --------------------------------------------------------------------------------------
# Operations
# --------------------------------------------------------------------------------------


def unexposed_state(material: Material, light: LightField, ctx: RunContext) -> EmulsionState:
    """An empty state of ``material`` for a frame lit by ``light`` in the run ``ctx``.

    Its geometry is ``FrameGeometry.from_long_edge(light.shape_hw, ctx.long_edge_mm)`` and
    its grid and seed are those of ``ctx``, as :func:`expose` requires; exposure models
    start from it and add their fields.
    """
    geometry = FrameGeometry.from_long_edge(light.shape_hw, ctx.long_edge_mm)
    return EmulsionState(material, geometry, ctx.grid, ctx.seed)


def expose(material: Material, light: LightField, time_s: float, ctx: RunContext) -> EmulsionState:
    """Expose ``material`` to ``light`` for ``time_s`` seconds.

    ``light`` is any :class:`~darkroom_atlas.light.LightField` (a full
    :class:`~darkroom_atlas.light.Light` or a
    :class:`~darkroom_atlas.light.LowRankLight`); it must be on ``ctx.grid`` and ``time_s``
    positive and finite, otherwise ``ValueError``. Dispatches to the exposure model
    registered as ``material.model``
    (:class:`~darkroom_atlas.registry.RegistryError` if there is none) and checks that the
    state it returns has the geometry of :func:`unexposed_state`, the run's grid and seed.
    ``light.meta`` is stored as ``state.meta["light"]``.
    """
    _check_type("material", material, Material)
    _check_type("light", light, LightField)
    _check_type("ctx", ctx, RunContext)
    _check_grid("light", light.grid, ctx)
    time_s = _positive_finite("time_s", time_s)
    model = get_exposure(material.model)
    state = model(material, light, time_s, ctx)

    what = f"exposure model {material.model!r}"
    if not isinstance(state, EmulsionState):
        raise TypeError(f"{what} returned a {type(state).__name__}, not an EmulsionState")
    expected = unexposed_state(material, light, ctx)
    for name in ("geometry", "grid", "seed"):
        got, want = getattr(state, name), getattr(expected, name)
        if got != want:
            raise ValueError(
                f"{what} returned a state with {name} {got}, expected {want}; build the "
                "state from pipeline.unexposed_state"
            )
    return state.with_meta(light=light.meta)


def segments(proc: Process, model: str) -> list[tuple[int, ...]]:
    """Step indices of ``proc`` grouped into the segments :func:`process` runs (spec 4.5).

    Adjacent steps whose handlers (for material model ``model``) declare a finite
    ``interaction_range_um`` form one segment; a step with ``interaction_range_um=None``
    is a segment of its own. A step without a handler raises
    :class:`~darkroom_atlas.registry.RegistryError`.
    """
    _check_type("proc", proc, Process)
    handlers = [get_step(model, step.kind) for step in proc.steps]
    return _group(handlers)


def process(state: EmulsionState, proc: Process, ctx: RunContext) -> EmulsionState:
    """Run the steps of ``proc`` on ``state`` in order and return the processed state.

    Every step's handler, ``(state.material.model, step.kind)`` in the registry, and every
    named bath are looked up before the first step runs, so a process that cannot run fails
    at once (:class:`~darkroom_atlas.registry.RegistryError` or
    :class:`~darkroom_atlas.data.library.DataNotFound`). ``state`` must be on the run's grid
    and carry the run's seed (``ValueError``). The steps are grouped by :func:`segments`
    and each segment is run (tiled or not) as described in the module docs; each step gets
    a :class:`StepContext` with its index and bath, and a
    :class:`~darkroom_atlas.state.StepRecord` is appended to the history after it.
    """
    _check_type("state", state, EmulsionState)
    _check_type("proc", proc, Process)
    _check_type("ctx", ctx, RunContext)
    _check_grid("state", state.grid, ctx)
    _check_seed(state, ctx)
    model = state.material.model
    handlers = [get_step(model, step.kind) for step in proc.steps]
    baths = [None if step.bath is None else ctx.library.bath(step.bath) for step in proc.steps]
    steps = [
        _Step(
            handler=handler,
            step=step,
            ctx=StepContext(ctx, index, bath),
            record=StepRecord(index, step.kind, step.bath, step.time_s, step.note or ""),
            what=f"the {step.kind.value!r} handler of model {model!r}",
        )
        for index, (step, handler, bath) in enumerate(zip(proc.steps, handlers, baths, strict=True))
    ]

    store: FieldStore | None = ctx.store
    for indices in _group(handlers):
        segment = [steps[i] for i in indices]
        whole_frame = ctx.tile_px is None or segment[0].handler.interaction_range_um is None
        if whole_frame:
            for item in segment:
                state = _run_whole(item, state)
        else:
            if store is None:
                store = _field_store(ctx)
            state = _run_tiled(segment, state, ctx.tile_px, store)
    return state


def readout(state: EmulsionState, device: str, ctx: RunContext, **params: Any) -> Readout:
    """Read ``state`` out with the readout device registered as ``device``.

    ``state`` must be on the run's grid and carry the run's seed (``ValueError``).
    ``params`` are passed on to the device. The result must be a :class:`Readout` whose data
    cover the frame.
    """
    _check_type("state", state, EmulsionState)
    _check_type("ctx", ctx, RunContext)
    _check_grid("state", state.grid, ctx)
    _check_seed(state, ctx)
    device_fn = get_readout(device)
    result = device_fn(state, ctx, **params)
    what = f"readout device {device!r}"
    if not isinstance(result, Readout):
        raise TypeError(f"{what} returned a {type(result).__name__}, not a Readout")
    if result.data.shape[:2] != state.geometry.shape_hw:
        raise ValueError(
            f"{what} returned data of shape {result.data.shape} for a frame of shape "
            f"{state.geometry.shape_hw}"
        )
    return result


# --------------------------------------------------------------------------------------
# Segments and tiling
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class _Step:
    """One step of a process, ready to run."""

    handler: RegisteredStep
    step: ProcessStep
    ctx: StepContext
    record: StepRecord
    what: str


def _group(handlers: list[RegisteredStep]) -> list[tuple[int, ...]]:
    groups: list[tuple[int, ...]] = []
    current: list[int] = []
    for index, handler in enumerate(handlers):
        if handler.interaction_range_um is None:
            if current:
                groups.append(tuple(current))
                current = []
            groups.append((index,))
        else:
            current.append(index)
    if current:
        groups.append(tuple(current))
    return groups


def _field_store(ctx: RunContext) -> FieldStore | None:
    """A store for the fields one ``process`` call assembles when the run has none.

    ``None`` (in memory) without ``ctx.scratch_dir``. Otherwise each call gets its own
    subdirectory of ``ctx.scratch_dir``. The store is not closed: the returned state refers
    to its files, which the runner deletes with the scratch directory.
    """
    if ctx.scratch_dir is None:
        return None
    return FieldStore(ctx.scratch_dir / f"process-{uuid.uuid4().hex[:12]}")


def _run_whole(item: _Step, state: EmulsionState) -> EmulsionState:
    result = _checked_result(item.handler.fn(state, item.step, item.ctx), state, item.what)
    return result.with_history(item.record)


def _run_tiled(
    segment: list[_Step], state: EmulsionState, tile_px: int, store: FieldStore | None
) -> EmulsionState:
    """Run a segment of finite-range steps tile by tile; assemble its final fields."""
    pitch = state.geometry.pixel_pitch_um
    # Each step needs its own reach of valid input around a tile, so the halos add up.
    overlap_px = sum(
        math.ceil(item.handler.interaction_range_um / pitch)  # type: ignore[operator]
        for item in segment
    )
    # Store names of the segment's final fields (spec §4.3): the last step's index.
    name_prefix = f"s{segment[-1].ctx.step_index}:"
    first: list[tuple[frozenset[str], dict[str, str]] | None] = [None] * len(segment)
    assemblers: dict[str, FrameAssembler] = {}
    done: list[Tile] = []
    final_keys: frozenset[str] = frozenset()
    meta: Mapping[str, Any] = state.meta

    for tile in iter_tiles(state.geometry.shape_hw, tile_px, overlap_px):
        window = state.crop(tile.outer)
        current = window
        for position, item in enumerate(segment):
            result = _checked_result(
                item.handler.fn(current, item.step, item.ctx), current, item.what
            )
            first[position] = _check_tile_fields(result, tile, first[position], item.what)
            current = result.with_history(item.record)
        if not done:
            final_keys, meta = frozenset(current.fields), current.meta

        for key, value in current.fields.items():
            if key not in assemblers:
                if _same_view(value, window.fields.get(key)):
                    continue  # still the input's array: pass it through
                spec = FieldSpec(key, value.shape[2:], value.dtype.name)
                assembler = FrameAssembler(state.geometry, [spec], store, name_prefix=name_prefix)
                # Earlier tiles passed this field through unchanged.
                for earlier in done:
                    assembler.write(key, earlier.inner, state.fields[key][earlier.inner])
                assemblers[key] = assembler
            assemblers[key].write(key, tile.inner, value[tile.inner_in_outer])
        done.append(tile)

    fields = {
        key: assemblers[key].finish()[key] if key in assemblers else state.fields[key]
        for key in sorted(final_keys)
    }
    history = state.history + tuple(item.record for item in segment)
    return replace(state, fields=fields, meta=meta, history=history)


def _check_tile_fields(
    result: EmulsionState,
    tile: Tile,
    first: tuple[frozenset[str], dict[str, str]] | None,
    what: str,
) -> tuple[frozenset[str], dict[str, str]]:
    """Check one tile's result against the first tile's fields and layouts."""
    if first is None:
        keys = frozenset(result.fields)
        return keys, {key: _layout(value) for key, value in result.fields.items()}
    keys, layouts = first
    if frozenset(result.fields) != keys:
        raise ValueError(
            f"{what} returned fields {sorted(result.fields)} for tile {tile.inner} but "
            f"{sorted(keys)} for the first tile; every tile must give the same fields"
        )
    for key, value in result.fields.items():
        if _layout(value) != layouts[key]:
            raise ValueError(
                f"{what} returned field {key!r} as {_layout(value)} for tile "
                f"{tile.inner} but as {layouts[key]} for the first tile; every tile "
                "must give each field the same dtype and trailing dimensions"
            )
    return first


def _same_view(value: np.ndarray, given: np.ndarray | None) -> bool:
    """Whether ``value`` is the very array view ``given`` (same memory, shape and strides)."""
    return (
        given is not None
        and value.dtype == given.dtype
        and value.shape == given.shape
        and value.strides == given.strides
        and value.__array_interface__["data"][0] == given.__array_interface__["data"][0]
    )


def _layout(field_value: np.ndarray) -> str:
    """A field's dtype and trailing dimensions, which every tile must agree on."""
    return f"{field_value.dtype} with trailing dimensions {field_value.shape[2:]}"


def _checked_result(result: object, given: EmulsionState, what: str) -> EmulsionState:
    if not isinstance(result, EmulsionState):
        raise TypeError(f"{what} returned a {type(result).__name__}, not an EmulsionState")
    for name in ("geometry", "grid"):
        got, want = getattr(result, name), getattr(given, name)
        if got != want:
            raise ValueError(f"{what} changed the state's {name} from {want} to {got}")
    return result


# --------------------------------------------------------------------------------------
# Argument checks
# --------------------------------------------------------------------------------------


def _check_type(name: str, value: object, expected: type) -> None:
    if not isinstance(value, expected):
        raise TypeError(f"{name} must be a {expected.__name__}, got {type(value).__name__}")


def _check_seed(state: EmulsionState, ctx: RunContext) -> None:
    if state.seed != ctx.seed:
        raise ValueError(
            f"the state was made with seed {state.seed}, but the run uses seed {ctx.seed}; "
            "use the same seed for both"
        )


def _check_grid(what: str, grid: SpectralGrid, ctx: RunContext) -> None:
    if grid != ctx.grid:
        raise ValueError(
            f"the {what} is on the spectral grid {grid}, but the run uses the grid "
            f"{ctx.grid}; use the same grid for both"
        )


def _int_at_least(name: str, value: object, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}, got {value!r}")
    return int(value)


def _positive_finite(name: str, value: object) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not math.isfinite(value)
        or value <= 0.0
    ):
        raise ValueError(f"{name} must be a positive, finite number, got {value!r}")
    return float(value)


# Importing the built-in models registers them. It comes last because they import the
# classes defined above.
import darkroom_atlas.models  # noqa: F401
