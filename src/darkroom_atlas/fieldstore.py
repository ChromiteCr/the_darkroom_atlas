"""Allocation of frame-sized fields in memory or in memory-mapped scratch files (spec §4.3).

A :class:`FieldStore` hands out writable arrays under a store name such as
``"s3:emulsion/silver"`` (step index and field key; spec §4.3). Names are unique for the
lifetime of a store: a name used once raises ``ValueError`` when allocated again, even after
it was released. Fields larger than ``ram_threshold_bytes`` go to an
``np.memmap(mode="w+")`` file under ``root`` (when ``root`` is not None); smaller ones stay
in RAM. Every backing file gets a fresh name
``<name with / and : replaced by _>-<uuid4 hex[:12]>.bin``, so a new allocation never
truncates a file that an earlier, still mapped field lives in.

Once a field is written, :meth:`FieldStore.freeze` returns it read-only: a RAM field is
marked read-only in place, a mapped field is flushed and reopened with ``mode="r"``, which
:func:`darkroom_atlas._frozen.readonly_array` shares without copying.
:meth:`FieldStore.release` deletes a field's backing file; :meth:`FieldStore.name_of`
finds the store name of a field handed on (for example a state field a later step drops).

Arrays never cross process boundaries: a mapped field is described by a picklable
:class:`FieldDescriptor` (:meth:`FieldStore.descriptor`, :func:`descriptor_of`) and a worker
opens it itself with :func:`open_field` (inputs ``"r"``, preallocated outputs ``"r+"``).

``root`` is the run's scratch directory; the runner owns it (creates it, and removes it
at the end), so the store deletes only the files it created.
"""

import math
import mmap
import re
import uuid
from dataclasses import dataclass
from numbers import Integral
from pathlib import Path
from typing import Literal, NoReturn, Self

import numpy as np

__all__ = ["FieldDescriptor", "FieldStore", "descriptor_of", "open_field"]

# Default size above which a field is memory-mapped (spec §4.3).
_DEFAULT_THRESHOLD = 256 * 2**20

# Spec §4.3 replaces "/" and ":"; anything else outside a portable file name too.
_UNSAFE = re.compile(r"[^A-Za-z0-9._+#-]")


@dataclass(frozen=True)
class FieldDescriptor:
    """Where a memory-mapped field lives: what crosses process boundaries (spec §4.3).

    ``path`` is the backing file, ``dtype`` the NumPy dtype name, ``shape`` the array shape
    and ``offset`` the byte offset of element 0 in the file (C order). A dtype name carries
    no byte order, so non-native byte orders raise ``ValueError``.
    """

    path: Path
    dtype: str
    shape: tuple[int, ...]
    offset: int = 0

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", Path(self.path))
        object.__setattr__(self, "dtype", _native(np.dtype(self.dtype)).name)
        object.__setattr__(self, "shape", _shape(self.shape))
        offset = self.offset
        if isinstance(offset, bool) or not isinstance(offset, Integral) or offset < 0:
            raise ValueError(f"offset must be an integer >= 0, got {offset!r}")
        object.__setattr__(self, "offset", int(offset))


def open_field(desc: FieldDescriptor, mode: Literal["r", "r+"]) -> np.ndarray:
    """Open the field ``desc`` describes as an ``np.memmap`` (``"r"`` or ``"r+"``).

    A worker writing through an ``"r+"`` map should ``flush()`` it before returning.
    ``FileNotFoundError`` if the file is gone.
    """
    if not isinstance(desc, FieldDescriptor):
        raise TypeError(f"desc must be a FieldDescriptor, got {type(desc).__name__}")
    if mode not in ("r", "r+"):
        raise ValueError(f"mode must be 'r' or 'r+', got {mode!r}")
    return np.memmap(desc.path, dtype=desc.dtype, mode=mode, shape=desc.shape, offset=desc.offset)


def descriptor_of(array: np.ndarray) -> FieldDescriptor:
    """Descriptor of a whole memory-mapped array (for example a frozen state field).

    ``array`` must be an ``np.memmap`` that owns its mapping (not a view of one) and is
    C-contiguous; anything else raises ``ValueError``, since a descriptor cannot express
    strides. A non-native byte order raises ``ValueError`` too.
    """
    if (
        not isinstance(array, np.memmap)
        or not isinstance(array.base, mmap.mmap)
        or array.filename is None
        or not array.flags.c_contiguous
    ):
        raise ValueError(
            "descriptor_of needs a whole, C-contiguous np.memmap (not a view or RAM array)"
        )
    dtype = _native(array.dtype).name
    return FieldDescriptor(Path(array.filename), dtype, array.shape, array.offset)


@dataclass
class _Entry:
    array: np.ndarray
    path: Path | None
    frozen: bool = False


class FieldStore:
    """Allocate fields in RAM or as memory-mapped files under ``root``, by store name.

    Writers must not keep writing through the array returned by :meth:`allocate` (or views
    of it) after :meth:`freeze`; the frozen array is the one to hand on. Usable as a context
    manager: leaving the ``with`` block calls :meth:`close`. A store is never pickled (it
    would carry whole mapped fields into a worker); pass descriptors instead.
    """

    def __init__(self, root: Path | None, *, ram_threshold_bytes: int = _DEFAULT_THRESHOLD):
        if isinstance(ram_threshold_bytes, bool) or not isinstance(ram_threshold_bytes, Integral):
            raise TypeError(f"ram_threshold_bytes must be an integer, got {ram_threshold_bytes!r}")
        if ram_threshold_bytes < 0:
            raise ValueError(f"ram_threshold_bytes must be >= 0, got {ram_threshold_bytes!r}")
        self.root = None if root is None else Path(root)
        self.ram_threshold_bytes = int(ram_threshold_bytes)
        self._entries: dict[str, _Entry] = {}
        self._used: set[str] = set()
        self._closed = False

    # ----------------------------------------------------------------------------------

    def allocate(self, name: str, shape: tuple[int, ...], dtype: str = "float64") -> np.ndarray:
        """Return a new zero-filled writable array under store name ``name``.

        It is an ``np.memmap(mode="w+")`` file under ``root`` when its size is above the
        threshold and ``root`` is set, otherwise an in-memory array. ``ValueError`` if
        ``name`` was allocated before in this store (released or not), the shape is
        invalid or ``dtype`` has a non-native byte order (a descriptor could not express it).
        """
        self._check_open()
        if not isinstance(name, str) or not name:
            raise ValueError(f"name must be a non-empty string, got {name!r}")
        if name in self._used:
            raise ValueError(
                f"store name {name!r} was already allocated in this store; names are unique "
                "per run (spec §4.3)"
            )
        shape = _shape(shape)
        dt = _native(np.dtype(dtype))
        nbytes = math.prod(shape) * dt.itemsize
        if self.root is not None and nbytes > self.ram_threshold_bytes:
            self.root.mkdir(parents=True, exist_ok=True)
            path = self.root / f"{_UNSAFE.sub('_', name)}-{uuid.uuid4().hex[:12]}.bin"
            array: np.ndarray = np.memmap(path, dtype=dt, mode="w+", shape=shape)
        else:
            path = None
            array = np.zeros(shape, dtype=dt)
        self._used.add(name)
        self._entries[name] = _Entry(array, path)
        return array

    def freeze(self, name: str) -> np.ndarray:
        """Return field ``name`` read-only (RAM array or ``np.memmap`` with ``mode="r"``).

        Freezing again returns the same frozen array.
        """
        entry = self._entry(name)
        if not entry.frozen:
            if entry.path is None:
                entry.array.flags.writeable = False
            else:
                writable = entry.array
                assert isinstance(writable, np.memmap)
                writable.flush()
                shape, dtype = writable.shape, writable.dtype
                del writable
                entry.array = np.memmap(entry.path, dtype=dtype, mode="r", shape=shape)
            entry.frozen = True
        return entry.array

    def descriptor(self, name: str) -> FieldDescriptor:
        """Descriptor of mapped field ``name`` for a worker to open; RAM fields raise."""
        entry = self._entry(name)
        if entry.path is None:
            raise ValueError(
                f"field {name!r} is in RAM, not memory-mapped; only mapped fields have a descriptor"
            )
        return FieldDescriptor(entry.path, entry.array.dtype.name, entry.array.shape, 0)

    def name_of(self, array: np.ndarray) -> str:
        """Store name of the field ``array`` belongs to, for :meth:`release`.

        ``array`` is a field this store handed out (writable or frozen), a view of one, or a
        map of its backing file opened elsewhere (:func:`open_field`). A later step that
        drops a state field releases it with ``store.release(store.name_of(field))``.
        ``KeyError`` if the array is not a field of this store (or was released).
        """
        if not isinstance(array, np.ndarray):
            raise TypeError(f"array must be a NumPy array, got {type(array).__name__}")
        filename = getattr(array, "filename", None)  # np.memmap and its views
        path = None if filename is None else Path(filename).resolve()
        for name, entry in self._entries.items():
            current: object = array
            while isinstance(current, np.ndarray):
                if current is entry.array:
                    return name
                current = current.base
            if path is not None and entry.path is not None and entry.path.resolve() == path:
                return name
        raise KeyError("array is not a field of this store (or it was released)")

    def release(self, name: str) -> None:
        """Forget field ``name`` and delete its backing file, if any.

        On POSIX systems arrays still referring to a deleted file stay readable until they
        are garbage-collected. The name stays used.
        """
        entry = self._entry(name)
        del self._entries[name]
        path, array = entry.path, entry.array
        del entry
        if path is not None:
            if isinstance(array, np.memmap) and array.mode != "r":
                array.flush()
            del array
            path.unlink(missing_ok=True)

    def close(self) -> None:
        """Release every field of this store; later allocations raise ``RuntimeError``."""
        for name in list(self._entries):
            self.release(name)
        self._closed = True

    def __contains__(self, name: object) -> bool:
        return name in self._entries

    def is_mapped(self, name: str) -> bool:
        """Whether field ``name`` is backed by a file (rather than RAM)."""
        return self._entry(name).path is not None

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def __reduce__(self) -> NoReturn:
        raise TypeError(
            "a FieldStore cannot be pickled; pass FieldDescriptor objects to worker processes"
        )

    # ----------------------------------------------------------------------------------

    def _entry(self, name: str) -> _Entry:
        try:
            return self._entries[name]
        except KeyError:
            known = ", ".join(sorted(self._entries)) or "(none)"
            raise KeyError(f"no field {name!r} in this store; allocated: {known}") from None

    def _check_open(self) -> None:
        if self._closed:
            raise RuntimeError("this FieldStore is closed")


def _native(dtype: np.dtype) -> np.dtype:
    if not dtype.isnative:
        raise ValueError(
            f"dtype {dtype.str!r} has a non-native byte order; a FieldDescriptor names the "
            "dtype without one, so fields must use the native byte order"
        )
    return dtype


def _shape(shape: object) -> tuple[int, ...]:
    try:
        dims = tuple(shape)  # type: ignore[arg-type]
    except TypeError:
        raise ValueError(f"shape must be a tuple of integers, got {shape!r}") from None
    for d in dims:
        if isinstance(d, bool) or not isinstance(d, Integral) or d < 0:
            raise ValueError(f"shape must be a tuple of integers >= 0, got {shape!r}")
    return tuple(int(d) for d in dims)
