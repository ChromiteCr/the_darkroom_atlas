"""The state of an exposed and (partly) processed material, and its frame geometry.

An :class:`EmulsionState` is immutable: every operation returns a new state, and every field
array is stored read-only. Fields are keyed ``"<layer_id>/<name>"`` (``"x/exposure"``);
fields that belong to no layer use the reserved prefix ``"frame/<name>"``. The first two
dimensions of every field are the frame ``(H, W)``; trailing dimensions and dtype are free.

A state may describe a window of a larger frame (a tile). Its :class:`FrameGeometry` then
records where the window sits (``origin_px``) and how big the whole frame is
(``frame_shape_hw``), so that position-dependent computations such as
:func:`darkroom_atlas.rng.pixel_uniform` give the same result tiled or untiled.
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from numbers import Integral, Real
from typing import Any

import numpy as np

from darkroom_atlas._frozen import FrozenDict, freeze_mapping, readonly_array
from darkroom_atlas.data.models import NAME_PATTERN, RESERVED_LAYER_IDS, Material
from darkroom_atlas.spectral import SpectralGrid
from darkroom_atlas.tiling import Tile
from darkroom_atlas.units import pixel_pitch_um

__all__ = ["EmulsionState", "FrameGeometry", "StepRecord"]

Window = tuple[slice, slice]


@dataclass(frozen=True)
class FrameGeometry:
    """Pixel grid of a frame (or of a window of one) on the film plane.

    ``shape_hw`` is the size of this (possibly cropped) array, ``pixel_pitch_um`` the pixel
    spacing on the film in um, ``origin_px`` the position of its top-left pixel in the whole
    frame and ``frame_shape_hw`` the size of the whole frame (``shape_hw`` when omitted).
    """

    shape_hw: tuple[int, int]
    pixel_pitch_um: float
    origin_px: tuple[int, int] = (0, 0)
    frame_shape_hw: tuple[int, int] | None = None

    def __post_init__(self) -> None:
        shape = _int_pair("shape_hw", self.shape_hw, minimum=1)
        origin = _int_pair("origin_px", self.origin_px, minimum=0)
        frame = shape if self.frame_shape_hw is None else self.frame_shape_hw
        frame = _int_pair("frame_shape_hw", frame, minimum=1)
        pitch = self.pixel_pitch_um
        if isinstance(pitch, bool) or not isinstance(pitch, Real) or not math.isfinite(pitch):
            raise ValueError(f"pixel_pitch_um must be a finite number, got {pitch!r}")
        if pitch <= 0.0:
            raise ValueError(f"pixel_pitch_um must be positive, got {pitch!r}")
        if any(o + s > f for o, s, f in zip(origin, shape, frame, strict=True)):
            raise ValueError(f"a {shape} window at origin {origin} does not fit in a {frame} frame")
        object.__setattr__(self, "shape_hw", shape)
        object.__setattr__(self, "origin_px", origin)
        object.__setattr__(self, "frame_shape_hw", frame)
        object.__setattr__(self, "pixel_pitch_um", float(pitch))

    @classmethod
    def from_long_edge(
        cls, shape_hw: tuple[int, int], long_edge_mm: float = 36.0
    ) -> "FrameGeometry":
        """Geometry of an uncropped frame whose long edge spans ``long_edge_mm`` on the film."""
        shape = _int_pair("shape_hw", shape_hw, minimum=1)
        return cls(shape, pixel_pitch_um(long_edge_mm, shape))

    def um_to_px(self, um: float) -> float:
        """Convert a length on the film from um to pixels."""
        return float(um) / self.pixel_pitch_um

    def crop(self, window: Window) -> "FrameGeometry":
        """Geometry of ``window`` (a ``(rows, cols)`` pair of slices into this array)."""
        rows, cols = _window(window, self.shape_hw)
        return replace(
            self,
            shape_hw=(rows.stop - rows.start, cols.stop - cols.start),
            origin_px=(self.origin_px[0] + rows.start, self.origin_px[1] + cols.start),
        )


@dataclass(frozen=True)
class StepRecord:
    """What one processing step did: its position in the process, kind, bath and time."""

    index: int
    kind: str
    bath: str | None
    time_s: float | None
    note: str = ""

    def __post_init__(self) -> None:
        if isinstance(self.index, bool) or not isinstance(self.index, Integral) or self.index < 0:
            raise ValueError(f"index must be a non-negative integer, got {self.index!r}")
        if not isinstance(self.kind, str) or not self.kind:
            raise ValueError(f"kind must be a non-empty string, got {self.kind!r}")
        if self.bath is not None and not isinstance(self.bath, str):
            raise ValueError(f"bath must be a string or None, got {self.bath!r}")
        if self.time_s is not None:
            t = self.time_s
            if isinstance(t, bool) or not isinstance(t, Real) or not math.isfinite(t) or t < 0:
                raise ValueError(f"time_s must be a non-negative number or None, got {t!r}")
            object.__setattr__(self, "time_s", float(t))
        if not isinstance(self.note, str):
            raise TypeError(f"note must be a string, got {self.note!r}")
        object.__setattr__(self, "index", int(self.index))
        object.__setattr__(self, "kind", str(self.kind))  # a StepKind becomes its value


@dataclass(frozen=True, eq=False)
class EmulsionState:
    """Immutable state of a material after exposure and any processing steps so far.

    ``fields`` maps ``"<layer_id>/<name>"`` (or ``"frame/<name>"``) to arrays whose first two
    dimensions equal ``geometry.shape_hw``; ``layer_id`` must be a layer of ``material``.
    Arrays are stored read-only (copied unless already read-only), ``fields`` and ``meta``
    as read-only mappings. ``history`` lists the processing steps applied so far.
    """

    material: Material
    geometry: FrameGeometry
    grid: SpectralGrid
    seed: int
    fields: Mapping[str, np.ndarray] = field(default_factory=dict)
    history: tuple[StepRecord, ...] = ()
    meta: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.geometry, FrameGeometry):
            raise TypeError(f"geometry must be a FrameGeometry, got {type(self.geometry)}")
        if not isinstance(self.grid, SpectralGrid):
            raise TypeError(f"grid must be a SpectralGrid, got {type(self.grid).__name__}")
        if isinstance(self.seed, bool) or not isinstance(self.seed, Integral) or self.seed < 0:
            raise ValueError(f"seed must be a non-negative integer, got {self.seed!r}")
        if not isinstance(self.fields, Mapping):
            raise TypeError(f"fields must be a mapping, got {type(self.fields).__name__}")
        prefixes = _field_prefixes(self.material)
        fields = FrozenDict(
            (key, self._checked_field(key, value, prefixes)) for key, value in self.fields.items()
        )
        history = tuple(self.history)
        for record in history:
            if not isinstance(record, StepRecord):
                raise TypeError(f"history entries must be StepRecord, got {type(record)}")
        object.__setattr__(self, "seed", int(self.seed))
        object.__setattr__(self, "fields", fields)
        object.__setattr__(self, "history", history)
        object.__setattr__(self, "meta", freeze_mapping(self.meta))

    def _checked_field(self, key: object, value: Any, prefixes: frozenset[str] | None):
        if not isinstance(key, str):
            raise TypeError(f"field keys must be strings, got {key!r}")
        layer_id, sep, name = key.partition("/")
        if not sep or NAME_PATTERN.fullmatch(layer_id) is None:
            raise ValueError(f"field key {key!r} must look like '<layer_id>/<name>'")
        if NAME_PATTERN.fullmatch(name) is None:
            raise ValueError(
                f"field key {key!r}: the name must be one segment of letters, digits and "
                "'. _ + -', starting with a letter or digit"
            )
        if prefixes is not None and layer_id not in prefixes:
            raise ValueError(
                f"field key {key!r}: {layer_id!r} is not a layer of material "
                f"{getattr(self.material, 'id', '?')!r}; use one of "
                f"{', '.join(sorted(prefixes))}"
            )
        try:
            array = readonly_array(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"field {key!r} is not an array: {exc}") from exc
        if array.dtype.kind not in "biufc":
            raise ValueError(f"field {key!r} must be numeric or boolean, got dtype {array.dtype}")
        if array.ndim < 2 or array.shape[:2] != self.geometry.shape_hw:
            raise ValueError(
                f"field {key!r} has shape {array.shape}; its first two dimensions must be "
                f"the frame shape {self.geometry.shape_hw}"
            )
        return array

    # ----------------------------------------------------------------------------------
    # Fields, history and metadata
    # ----------------------------------------------------------------------------------

    def get(self, layer_id: str, name: str) -> np.ndarray:
        """Return field ``"<layer_id>/<name>"``; ``KeyError`` listing the keys if absent."""
        key = f"{layer_id}/{name}"
        try:
            return self.fields[key]
        except KeyError:
            available = ", ".join(sorted(self.fields)) or "(none)"
            raise KeyError(f"no field {key!r}; available fields: {available}") from None

    def has(self, layer_id: str, name: str) -> bool:
        """Whether field ``"<layer_id>/<name>"`` exists."""
        return f"{layer_id}/{name}" in self.fields

    def with_fields(self, updates: Mapping[str, np.ndarray]) -> "EmulsionState":
        """Return a new state with ``updates`` added to (or replacing) the fields."""
        return replace(self, fields={**self.fields, **updates})

    def with_history(self, record: StepRecord) -> "EmulsionState":
        """Return a new state with ``record`` appended to the history."""
        if not isinstance(record, StepRecord):
            raise TypeError(f"record must be a StepRecord, got {type(record).__name__}")
        return replace(self, history=(*self.history, record))

    def with_meta(self, **kw: Any) -> "EmulsionState":
        """Return a new state with the keyword arguments added to (or replacing) ``meta``."""
        return replace(self, meta={**self.meta, **kw})

    # ----------------------------------------------------------------------------------
    # Tiling
    # ----------------------------------------------------------------------------------

    def crop(self, window: Window) -> "EmulsionState":
        """Return the state of ``window``, a ``(rows, cols)`` pair of slices into this frame.

        Slices must have a step of 1, bounds within the frame and a non-empty extent. Fields
        are cut to the window (as read-only views) and ``geometry`` gets the window's shape
        and its origin in the whole frame (crops of crops accumulate the offsets).
        """
        rows, cols = _window(window, self.geometry.shape_hw)
        return replace(
            self,
            geometry=self.geometry.crop((rows, cols)),
            fields={key: value[rows, cols] for key, value in self.fields.items()},
        )

    def paste(self, tile_state: "EmulsionState", tile: Tile) -> "EmulsionState":
        """Return a copy of this state with the inner part of ``tile_state`` written in.

        ``tile_state`` must be the state of window ``tile.outer`` of this frame (for example
        a processed ``self.crop(tile.outer)``). For each of its fields, the
        ``tile.inner_in_outer`` part is written at ``tile.inner``; fields of this state that
        ``tile_state`` lacks are kept as they are. A field this state does not have yet
        starts as zeros with the trailing dimensions and dtype of ``tile_state``'s field, so
        pasting a full set of tiles defines it everywhere. A field this state has with other
        trailing dimensions or another dtype raises ``ValueError`` (rather than silently
        discarding what was pasted before).
        """
        if not isinstance(tile_state, EmulsionState):
            raise TypeError(f"tile_state must be an EmulsionState, got {type(tile_state)}")
        outer = _window(tile.outer, self.geometry.shape_hw)
        inner = _window(tile.inner, self.geometry.shape_hw)
        inner_in_outer = _window(tile.inner_in_outer, tile_state.geometry.shape_hw)
        expected = self.geometry.crop(outer)
        if tile_state.geometry != expected:
            raise ValueError(
                f"tile_state has geometry {tile_state.geometry}, but window {tile.outer} of "
                f"this state has {expected}"
            )
        if tile_state.grid != self.grid:
            raise ValueError(f"tile_state grid {tile_state.grid} differs from {self.grid}")
        for axis, (o, i, io) in enumerate(zip(outer, inner, inner_in_outer, strict=True)):
            if io.start != i.start - o.start or io.stop != i.stop - o.start:
                raise ValueError(
                    f"tile is inconsistent along axis {axis}: inner {i} inside outer {o} "
                    f"is not inner_in_outer {io}"
                )

        fields: dict[str, np.ndarray] = dict(self.fields)
        for key, value in tile_state.fields.items():
            source = value[inner_in_outer]
            current = self.fields.get(key)
            if current is None:
                target = np.zeros((*self.geometry.shape_hw, *source.shape[2:]), source.dtype)
            elif current.shape[2:] != source.shape[2:] or current.dtype != source.dtype:
                raise ValueError(
                    f"cannot paste field {key!r}: the tile has {source.dtype} with trailing "
                    f"dimensions {source.shape[2:]}, but this state has {current.dtype} with "
                    f"trailing dimensions {current.shape[2:]}"
                )
            else:
                target = current.copy()
            target[inner] = source
            target.flags.writeable = False
            fields[key] = target
        return replace(self, fields=fields)


def _field_prefixes(material: Any) -> frozenset[str] | None:
    layers = getattr(material, "layers", None)
    if layers is None:
        return None
    return frozenset(layer.id for layer in layers) | RESERVED_LAYER_IDS


def _int_pair(name: str, value: object, *, minimum: int) -> tuple[int, int]:
    try:
        a, b = value  # type: ignore[misc]
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a pair of integers, got {value!r}") from None
    for v in (a, b):
        if isinstance(v, bool) or not isinstance(v, Integral) or v < minimum:
            raise ValueError(f"{name} must be a pair of integers >= {minimum}, got {value!r}")
    return (int(a), int(b))


def _window(window: object, shape_hw: tuple[int, int]) -> Window:
    """Check a ``(rows, cols)`` pair of slices against ``shape_hw``; resolve ``None`` bounds."""
    try:
        rows, cols = window  # type: ignore[misc]
    except (TypeError, ValueError):
        raise ValueError(f"window must be a (rows, cols) pair of slices, got {window!r}") from None
    return (_bounded_slice("rows", rows, shape_hw[0]), _bounded_slice("cols", cols, shape_hw[1]))


def _bounded_slice(name: str, s: object, size: int) -> slice:
    if not isinstance(s, slice):
        raise TypeError(f"{name} must be a slice, got {s!r}")
    if s.step not in (None, 1):
        raise ValueError(f"{name} slice must have a step of 1, got {s!r}")
    start = 0 if s.start is None else s.start
    stop = size if s.stop is None else s.stop
    for bound in (start, stop):
        if isinstance(bound, bool) or not isinstance(bound, Integral) or not 0 <= bound <= size:
            raise ValueError(f"{name} slice {s!r} must lie within 0..{size}")
    if stop <= start:
        raise ValueError(f"{name} slice {s!r} is empty")
    return slice(int(start), int(stop))
