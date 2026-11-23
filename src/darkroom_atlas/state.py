"""The state of an exposed and (partly) processed material, and its frame geometry.

An :class:`EmulsionState` is immutable: every operation returns a new state, and every field
array is stored read-only. Fields are keyed ``"<layer_id>/<name>"`` (``"x/exposure"``);
fields that belong to no layer use the reserved prefix ``"frame/<name>"``. The first two
dimensions of every field are the frame ``(H, W)``; trailing dimensions and dtype are free.

A state may describe a window of a larger frame (a tile). Its :class:`FrameGeometry` then
records where the window sits (``origin_px``) and how big the whole frame is
(``frame_shape_hw``), so that position-dependent computations such as
:func:`darkroom_atlas.rng.pixel_uniform` give the same result tiled or untiled.

A :class:`FrameAssembler` builds whole-frame output fields from per-block results: each
field is allocated once (through a :class:`darkroom_atlas.fieldstore.FieldStore`) and each
block writes only its inner part, so assembling costs O(pixels) (spec §4.4).
"""

import itertools
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from numbers import Integral, Real
from typing import Any

import numpy as np

from darkroom_atlas._frozen import FrozenDict, freeze_mapping, readonly_array
from darkroom_atlas.data.models import NAME_PATTERN, RESERVED_LAYER_IDS, Material
from darkroom_atlas.fieldstore import FieldDescriptor, FieldStore
from darkroom_atlas.spectral import SpectralGrid
from darkroom_atlas.tiling import Tile
from darkroom_atlas.units import pixel_pitch_um

__all__ = ["EmulsionState", "FieldSpec", "FrameAssembler", "FrameGeometry", "StepRecord"]

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

    def drop(self, *keys: str) -> "EmulsionState":
        """Return a new state without fields ``keys``; the other fields are shared as is.

        ``KeyError`` listing the available fields if a key is absent.
        """
        missing = [key for key in keys if key not in self.fields]
        if missing:
            available = ", ".join(sorted(self.fields)) or "(none)"
            raise KeyError(f"no field {missing[0]!r}; available fields: {available}")
        drop = set(keys)
        return replace(
            self, fields={key: value for key, value in self.fields.items() if key not in drop}
        )

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

        Each call copies every pasted field of the whole frame, so it costs O(frame pixels);
        assembling many tiles this way grows with tiles x pixels. Use
        :class:`FrameAssembler` to assemble a frame from tiles.
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


@dataclass(frozen=True)
class FieldSpec:
    """An output field of a :class:`FrameAssembler`: key, trailing dimensions and dtype."""

    key: str
    trailing: tuple[int, ...] = ()
    dtype: str = "float64"

    def __post_init__(self) -> None:
        if not isinstance(self.key, str) or not self.key:
            raise ValueError(f"key must be a non-empty string, got {self.key!r}")
        try:
            trailing = tuple(self.trailing)
        except TypeError:
            raise ValueError(
                f"trailing must be a tuple of integers, got {self.trailing!r}"
            ) from None
        for d in trailing:
            if isinstance(d, bool) or not isinstance(d, Integral) or d < 1:
                raise ValueError(f"trailing must be a tuple of integers >= 1, got {trailing!r}")
        dt = np.dtype(self.dtype)
        if dt.kind not in "biufc":
            raise ValueError(f"field {self.key!r}: dtype must be numeric or boolean, got {dt}")
        object.__setattr__(self, "trailing", tuple(int(d) for d in trailing))
        object.__setattr__(self, "dtype", dt.name)


# Side of the square buckets that index written regions for the overlap check: each write
# is compared only with the regions sharing a bucket, so checking costs O(1) per block.
_BUCKET_PX = 128

Rect = tuple[int, int, int, int]  # (r0, r1, c0, c1)


class _Coverage:
    """Regions written into one field since the last reset, indexed by bucket."""

    def __init__(self) -> None:
        self.rects: list[Rect] = []
        self.area = 0
        self._buckets: dict[tuple[int, int], list[int]] = {}

    @staticmethod
    def _cells(rect: Rect) -> itertools.product:
        r0, r1, c0, c1 = rect
        return itertools.product(
            range(r0 // _BUCKET_PX, (r1 - 1) // _BUCKET_PX + 1),
            range(c0 // _BUCKET_PX, (c1 - 1) // _BUCKET_PX + 1),
        )

    def overlap(self, rect: Rect) -> Rect | None:
        """A region written before that overlaps ``rect``, or None."""
        r0, r1, c0, c1 = rect
        for cell in self._cells(rect):
            for index in self._buckets.get(cell, ()):
                a0, a1, b0, b1 = self.rects[index]
                if r0 < a1 and a0 < r1 and c0 < b1 and b0 < c1:
                    return self.rects[index]
        return None

    def add(self, rect: Rect) -> None:
        index = len(self.rects)
        self.rects.append(rect)
        self.area += (rect[1] - rect[0]) * (rect[3] - rect[2])
        for cell in self._cells(rect):
            self._buckets.setdefault(cell, []).append(index)


class _Assembly:
    """Book-keeping for one field of a FrameAssembler."""

    def __init__(self, spec: FieldSpec, name: str) -> None:
        self.spec = spec
        self.name = name  # store name
        self.array: np.ndarray | None = None  # allocated on the first non-uniform write
        self.allocated = False
        self.mapped = False
        self.issued = False  # a descriptor was handed out for workers to write
        self.uniform: np.ndarray | None = None  # (1, 1, *trailing) while all writes match
        self.coverage = _Coverage()


class FrameAssembler:
    """Assemble whole-frame fields from per-block (or per-tile) results (spec §4.4).

    Each field of ``specs`` is allocated once, through ``store`` (RAM only when None), under
    the store name ``f"{name_prefix}{spec.key}"`` (steps pass ``name_prefix=f"s{index}:"``,
    spec §4.3), with shape ``(*geometry.shape_hw, *spec.trailing)`` and ``spec.dtype``.
    :meth:`write` copies a block's values straight into its inner region. Worker processes
    can instead write their inner block themselves into the outputs :meth:`descriptors`
    hands out (opened ``"r+"``); the main process then registers the region with
    :meth:`mark`. Every pixel of every field must be written (or marked) exactly once since
    construction or the last :meth:`reset`: overlapping writes raise ``ValueError`` at once,
    gaps raise in :meth:`finish`. Values are cast to the field's dtype only when the cast is
    safe. Copying and book-keeping cost O(pixels) in all, whatever the number of blocks.

    A field whose every write is a broadcast view (all strides of the two frame axes 0) of
    one bit-identical value is never allocated: :meth:`finish` returns it as a read-only
    broadcast view of a ``(1, 1, *trailing)`` array. Allocation is therefore deferred to the
    field's first non-uniform write (or to :meth:`descriptors`), which also fills in the
    regions written before. Releasing the store's files is up to the store's owner.
    """

    def __init__(
        self,
        geometry: FrameGeometry,
        specs: Sequence[FieldSpec],
        store: FieldStore | None = None,
        *,
        name_prefix: str = "",
    ) -> None:
        if not isinstance(geometry, FrameGeometry):
            raise TypeError(f"geometry must be a FrameGeometry, got {type(geometry).__name__}")
        if store is not None and not isinstance(store, FieldStore):
            raise TypeError(f"store must be a FieldStore or None, got {type(store).__name__}")
        if not isinstance(name_prefix, str):
            raise TypeError(f"name_prefix must be a string, got {type(name_prefix).__name__}")
        self.geometry = geometry
        self._store = FieldStore(None) if store is None else store
        self._fields: dict[str, _Assembly] = {}
        for spec in specs:
            if not isinstance(spec, FieldSpec):
                raise TypeError(f"specs must be FieldSpec instances, got {type(spec).__name__}")
            if spec.key in self._fields:
                raise ValueError(f"duplicate field spec {spec.key!r}")
            self._fields[spec.key] = _Assembly(spec, f"{name_prefix}{spec.key}")
        if not self._fields:
            raise ValueError("a FrameAssembler needs at least one FieldSpec")
        self._finished = False

    @property
    def store_names(self) -> dict[str, str]:
        """Field key -> store name of its array (for fields that were allocated)."""
        return {k: a.name for k, a in self._fields.items() if a.allocated}

    def write(self, key: str, inner: tuple[slice, slice], values: np.ndarray) -> None:
        """Write ``values`` (shape ``(rows, cols, *trailing)``) into region ``inner``."""
        self._check_open()
        item = self._item(key)
        rows, cols = _window(inner, self.geometry.shape_hw)
        spec = item.spec
        values = np.asarray(values)
        expected = (rows.stop - rows.start, cols.stop - cols.start, *spec.trailing)
        if values.shape != expected:
            raise ValueError(
                f"field {key!r}: values for region {(rows, cols)} must have shape {expected}, "
                f"got {values.shape}"
            )
        if not np.can_cast(values.dtype, spec.dtype, "safe"):
            raise ValueError(
                f"field {key!r}: cannot safely cast values of dtype {values.dtype} to {spec.dtype}"
            )
        rect = (rows.start, rows.stop, cols.start, cols.stop)
        self._check_overlap(key, item, rect)

        if item.array is None and values.strides[:2] == (0, 0):
            value = np.array(values[:1, :1], dtype=spec.dtype)
            if item.uniform is None and not item.coverage.rects:
                item.uniform = value
                item.coverage.add(rect)
                return
            if item.uniform is not None and value.tobytes() == item.uniform.tobytes():
                item.coverage.add(rect)
                return
        if item.array is None:
            self._allocate(item)
        assert item.array is not None
        item.array[rows, cols] = values
        item.coverage.add(rect)

    def write_state(self, tile_state: "EmulsionState", tile: Tile) -> None:
        """Write the inner part of ``tile_state`` (window ``tile.outer`` of the frame).

        Every spec field is taken from ``tile_state.fields[key][tile.inner_in_outer]`` and
        written at ``tile.inner``; fields of ``tile_state`` without a spec are ignored. A
        spec field missing from ``tile_state`` raises ``KeyError``.
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
                f"this frame has {expected}"
            )
        for axis, (o, i, io) in enumerate(zip(outer, inner, inner_in_outer, strict=True)):
            if io.start != i.start - o.start or io.stop != i.stop - o.start:
                raise ValueError(
                    f"tile is inconsistent along axis {axis}: inner {i} inside outer {o} "
                    f"is not inner_in_outer {io}"
                )
        missing = [key for key in self._fields if key not in tile_state.fields]
        if missing:
            raise KeyError(f"tile_state lacks field(s) {', '.join(missing)}")
        for key in self._fields:
            self.write(key, inner, tile_state.fields[key][inner_in_outer])

    def descriptors(self) -> dict[str, FieldDescriptor]:
        """Descriptors of the memory-mapped output fields, for workers to open ``"r+"``.

        Allocates every field not allocated yet (filling in uniform writes made so far). A
        worker writes its inner block straight into these maps (and flushes them); the main
        process then calls :meth:`mark` with that block's region. Fields that stay in RAM
        (store without ``root``, or below its threshold) are absent and must be written
        with :meth:`write`.
        """
        self._check_open()
        out: dict[str, FieldDescriptor] = {}
        for key, item in self._fields.items():
            if not item.allocated:
                self._allocate(item)
            if item.mapped:
                item.issued = True
                out[key] = self._store.descriptor(item.name)
        return out

    def mark(self, inner: tuple[slice, slice]) -> None:
        """Register region ``inner`` as written by a worker in every field of :meth:`descriptors`.

        ``ValueError`` if no descriptors were handed out or the region overlaps one written
        (or marked) before in any of those fields.
        """
        self._check_open()
        issued = [(key, item) for key, item in self._fields.items() if item.issued]
        if not issued:
            raise ValueError(
                "mark registers regions written through descriptors(); no descriptors were "
                "handed out (fields in RAM must be written with write())"
            )
        rows, cols = _window(inner, self.geometry.shape_hw)
        rect = (rows.start, rows.stop, cols.start, cols.stop)
        for key, item in issued:
            self._check_overlap(key, item, rect)
        for _, item in issued:
            item.coverage.add(rect)

    def reset(self) -> None:
        """Forget every write and mark, keeping the allocated outputs (re-pass, spec §3.7).

        The next pass overwrites the outputs in place; a field that was only ever written
        uniformly starts over.
        """
        self._check_open()
        for item in self._fields.values():
            item.coverage = _Coverage()
            if not item.allocated:
                item.uniform = None

    def finish(self) -> dict[str, np.ndarray]:
        """Check that every pixel of every field was written once; return read-only fields.

        Allocated fields are frozen through the store (read-only RAM arrays or ``mode="r"``
        memory maps); uniform fields are read-only broadcast views.
        """
        if self._finished:
            raise ValueError("this FrameAssembler is already finished")
        h, w = self.geometry.shape_hw
        for key, item in self._fields.items():
            if item.coverage.area != h * w:
                raise ValueError(
                    f"field {key!r}: {item.coverage.area} of {h * w} pixels written; every "
                    "pixel must be written (or marked) exactly once"
                )
        out: dict[str, np.ndarray] = {}
        for key, item in self._fields.items():
            if not item.allocated:
                assert item.uniform is not None
                item.uniform.flags.writeable = False
                out[key] = np.broadcast_to(item.uniform, (h, w, *item.spec.trailing))
            else:
                item.array = None
                out[key] = self._store.freeze(item.name)
        self._finished = True
        return out

    # ----------------------------------------------------------------------------------

    def _item(self, key: str) -> _Assembly:
        try:
            return self._fields[key]
        except KeyError:
            known = ", ".join(self._fields)
            raise KeyError(f"no field spec {key!r}; specs: {known}") from None

    def _check_open(self) -> None:
        if self._finished:
            raise ValueError("this FrameAssembler is finished")

    @staticmethod
    def _check_overlap(key: str, item: _Assembly, rect: Rect) -> None:
        hit = item.coverage.overlap(rect)
        if hit is not None:
            raise ValueError(
                f"field {key!r}: region {(slice(rect[0], rect[1]), slice(rect[2], rect[3]))} "
                f"overlaps region {(slice(hit[0], hit[1]), slice(hit[2], hit[3]))} written "
                "before"
            )

    def _allocate(self, item: _Assembly) -> None:
        spec = item.spec
        shape = (*self.geometry.shape_hw, *spec.trailing)
        array = self._store.allocate(item.name, shape, spec.dtype)
        if item.uniform is not None:
            for r0, r1, c0, c1 in item.coverage.rects:
                array[r0:r1, c0:c1] = item.uniform
            item.uniform = None
        item.array = array
        item.allocated = True
        item.mapped = self._store.is_mapped(item.name)


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
