"""The three operations of the simulation: ``expose``, ``process`` and ``readout``.

* :func:`expose` turns light on the film into an :class:`~darkroom_atlas.state.EmulsionState`
  with the exposure model that the material's ``model`` names.
* :func:`process` runs the steps of a process one by one, each with the handler registered
  for ``(material model, step kind)``, tiling the frame when the run asks for it.
* :func:`readout` turns a state into a :class:`Readout` (an image, spectra or densities)
  with a named readout device.

The physics lives in the registered models (:mod:`darkroom_atlas.registry`); this module
only dispatches, checks what the models return and keeps the history. A spectral
transmittance read out of a processed negative can be applied to a light with
:meth:`Light.filtered <darkroom_atlas.light.Light.filtered>` and exposed onto the next
material, which is how printing reuses the same engine.

Tiling: when ``RunContext.tile_px`` is set, a step whose handler declares a finite
``interaction_range_um`` is run on overlapping tiles (overlap
``ceil(interaction_range_um / pixel_pitch_um)`` pixels) cut with
:meth:`EmulsionState.crop` and stitched back with :meth:`EmulsionState.paste`; a step with
``interaction_range_um=None`` always receives the whole frame. A tiled handler must return
the same field names, each with the same dtype and trailing dimensions, and the same
``meta`` for every tile (differing fields raise ``ValueError``): the stitched state takes its
``meta`` from the first tile, and the handler must not touch the history (the pipeline
appends the step's record).
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from numbers import Integral, Real
from typing import Any, Literal, get_args

import numpy as np

from darkroom_atlas._frozen import freeze_mapping, readonly_array
from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import Bath, Material, Process, ProcessStep
from darkroom_atlas.light import Light
from darkroom_atlas.registry import RegisteredStep, get_exposure, get_readout, get_step
from darkroom_atlas.spectral import DEFAULT_GRID, SpectralGrid
from darkroom_atlas.state import EmulsionState, FrameGeometry, StepRecord
from darkroom_atlas.tiling import iter_tiles

__all__ = [
    "Readout",
    "ReadoutKind",
    "RunContext",
    "StepContext",
    "expose",
    "process",
    "readout",
    "unexposed_state",
]

ReadoutKind = Literal[
    "image_linear_srgb", "spectral_transmittance", "spectral_reflectance", "density"
]
_READOUT_KINDS: tuple[str, ...] = get_args(ReadoutKind)
_SPECTRAL_KINDS = frozenset({"spectral_transmittance", "spectral_reflectance"})


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
    """

    library: Library
    grid: SpectralGrid = DEFAULT_GRID
    seed: int = 0
    tile_px: int | None = None
    long_edge_mm: float = 36.0

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


def unexposed_state(material: Material, light: Light, ctx: RunContext) -> EmulsionState:
    """An empty state of ``material`` for a frame lit by ``light`` in the run ``ctx``.

    Its geometry is ``FrameGeometry.from_long_edge(light.shape_hw, ctx.long_edge_mm)`` and
    its grid and seed are those of ``ctx``, as :func:`expose` requires; exposure models
    start from it and add their fields.
    """
    geometry = FrameGeometry.from_long_edge(light.shape_hw, ctx.long_edge_mm)
    return EmulsionState(material, geometry, ctx.grid, ctx.seed)


def expose(material: Material, light: Light, time_s: float, ctx: RunContext) -> EmulsionState:
    """Expose ``material`` to ``light`` for ``time_s`` seconds.

    ``light`` must be on ``ctx.grid`` and ``time_s`` positive and finite, otherwise
    ``ValueError``. Dispatches to the exposure model registered as ``material.model``
    (:class:`~darkroom_atlas.registry.RegistryError` if there is none) and checks that the
    state it returns has the geometry of :func:`unexposed_state`, the run's grid and seed.
    ``light.meta`` is stored as ``state.meta["light"]``.
    """
    _check_type("material", material, Material)
    _check_type("light", light, Light)
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


def process(state: EmulsionState, proc: Process, ctx: RunContext) -> EmulsionState:
    """Run the steps of ``proc`` on ``state`` in order and return the processed state.

    Every step's handler, ``(state.material.model, step.kind)`` in the registry, and every
    named bath are looked up before the first step runs, so a process that cannot run fails
    at once (:class:`~darkroom_atlas.registry.RegistryError` or
    :class:`~darkroom_atlas.data.library.DataNotFound`). Each step gets a
    :class:`StepContext` with its index and bath; it is tiled as described in the module
    docs; and a :class:`~darkroom_atlas.state.StepRecord` is appended to the history.
    """
    _check_type("state", state, EmulsionState)
    _check_type("proc", proc, Process)
    _check_type("ctx", ctx, RunContext)
    _check_grid("state", state.grid, ctx)
    model = state.material.model
    handlers = [get_step(model, step.kind) for step in proc.steps]
    baths = [None if step.bath is None else ctx.library.bath(step.bath) for step in proc.steps]

    for index, (step, handler, bath) in enumerate(zip(proc.steps, handlers, baths, strict=True)):
        step_ctx = StepContext(ctx, index, bath)
        state = _run_step(handler, state, step, step_ctx)
        record = StepRecord(index, step.kind, step.bath, step.time_s, step.note or "")
        state = state.with_history(record)
    return state


def readout(state: EmulsionState, device: str, ctx: RunContext, **params: Any) -> Readout:
    """Read ``state`` out with the readout device registered as ``device``.

    ``params`` are passed on to the device. The result must be a :class:`Readout` whose data
    cover the frame.
    """
    _check_type("state", state, EmulsionState)
    _check_type("ctx", ctx, RunContext)
    _check_grid("state", state.grid, ctx)
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
# Steps and tiling
# --------------------------------------------------------------------------------------


def _run_step(
    handler: RegisteredStep, state: EmulsionState, step: ProcessStep, ctx: StepContext
) -> EmulsionState:
    what = f"the {step.kind.value!r} handler of model {state.material.model!r}"
    tile_px = ctx.run.tile_px
    if tile_px is None or handler.interaction_range_um is None:
        return _checked_result(handler.fn(state, step, ctx), state, what)

    overlap_px = math.ceil(handler.interaction_range_um / state.geometry.pixel_pitch_um)
    # Start from no fields, so that fields the step removes stay removed, as untiled.
    assembled = replace(state, fields={})
    keys: frozenset[str] | None = None
    layouts: dict[str, str] = {}
    meta: Mapping[str, Any] = state.meta
    for tile in iter_tiles(state.geometry.shape_hw, tile_px, overlap_px):
        window = state.crop(tile.outer)
        result = _checked_result(handler.fn(window, step, ctx), window, what)
        if keys is None:
            keys, meta = frozenset(result.fields), result.meta
            layouts = {key: _layout(value) for key, value in result.fields.items()}
        elif frozenset(result.fields) != keys:
            raise ValueError(
                f"{what} returned fields {sorted(result.fields)} for tile {tile.inner} but "
                f"{sorted(keys)} for the first tile; every tile must give the same fields"
            )
        else:
            for key, value in result.fields.items():
                if _layout(value) != layouts[key]:
                    raise ValueError(
                        f"{what} returned field {key!r} as {_layout(value)} for tile "
                        f"{tile.inner} but as {layouts[key]} for the first tile; every tile "
                        "must give each field the same dtype and trailing dimensions"
                    )
        assembled = assembled.paste(result, tile)
    return replace(assembled, meta=meta)


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
