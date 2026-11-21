"""Allocation of frame-sized fields in memory or in memory-mapped scratch files (spec §4.3).

A :class:`FieldStore` hands out writable arrays. Fields larger than ``ram_threshold_bytes``
go to an ``np.memmap(mode="w+")`` file under ``root`` (when ``root`` is not None); smaller
ones stay in RAM. Once a field is written, :meth:`FieldStore.freeze` returns it read-only:
a RAM field is marked read-only in place, a mapped field is flushed and reopened with
``mode="r"``, which :func:`darkroom_atlas._frozen.readonly_array` shares without copying.
:meth:`FieldStore.release` deletes a field's backing file.

``root`` is the run's scratch directory; the runner owns it (creates it, and removes it
at the end), so the store deletes only the files it created.
"""

import itertools
import math
import re
from dataclasses import dataclass
from numbers import Integral
from pathlib import Path
from typing import Self

import numpy as np

__all__ = ["FieldStore"]

# Default size above which a field is memory-mapped (spec §4.3).
_DEFAULT_THRESHOLD = 256 * 2**20

_UNSAFE = re.compile(r"[^A-Za-z0-9._+-]+")


@dataclass
class _Entry:
    array: np.ndarray
    path: Path | None
    frozen: bool = False


class FieldStore:
    """Allocate fields in RAM or as memory-mapped files under ``root``.

    Writers must not keep writing through the array returned by :meth:`allocate` (or views
    of it) after :meth:`freeze`; the frozen array is the one to hand on. Usable as a context
    manager: leaving the ``with`` block calls :meth:`close`.
    """

    def __init__(self, root: Path | None, *, ram_threshold_bytes: int = _DEFAULT_THRESHOLD):
        if isinstance(ram_threshold_bytes, bool) or not isinstance(ram_threshold_bytes, Integral):
            raise TypeError(f"ram_threshold_bytes must be an integer, got {ram_threshold_bytes!r}")
        if ram_threshold_bytes < 0:
            raise ValueError(f"ram_threshold_bytes must be >= 0, got {ram_threshold_bytes!r}")
        self.root = None if root is None else Path(root)
        self.ram_threshold_bytes = int(ram_threshold_bytes)
        self._entries: dict[str, _Entry] = {}
        self._serial = itertools.count()
        self._closed = False

    # ----------------------------------------------------------------------------------

    def allocate(self, key: str, shape: tuple[int, ...], dtype: str = "float64") -> np.ndarray:
        """Return a new zero-filled writable array for ``key``.

        It is an ``np.memmap(mode="w+")`` file under ``root`` when its size is above the
        threshold and ``root`` is set, otherwise an in-memory array. ``ValueError`` if
        ``key`` is already allocated (and not released) or the shape is invalid.
        """
        self._check_open()
        if not isinstance(key, str) or not key:
            raise ValueError(f"key must be a non-empty string, got {key!r}")
        if key in self._entries:
            raise ValueError(f"field {key!r} is already allocated in this store")
        shape = _shape(shape)
        dt = np.dtype(dtype)
        nbytes = math.prod(shape) * dt.itemsize
        if self.root is not None and nbytes > self.ram_threshold_bytes:
            self.root.mkdir(parents=True, exist_ok=True)
            path = self.root / f"{next(self._serial):05d}-{_UNSAFE.sub('_', key)}.dat"
            array: np.ndarray = np.memmap(path, dtype=dt, mode="w+", shape=shape)
        else:
            path = None
            array = np.zeros(shape, dtype=dt)
        self._entries[key] = _Entry(array, path)
        return array

    def freeze(self, key: str) -> np.ndarray:
        """Return field ``key`` read-only (RAM array or ``np.memmap`` with ``mode="r"``).

        Freezing again returns the same frozen array.
        """
        entry = self._entry(key)
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

    def release(self, key: str) -> None:
        """Forget field ``key`` and delete its backing file, if any.

        On POSIX systems arrays still referring to a deleted file stay readable until they
        are garbage-collected.
        """
        entry = self._entry(key)
        del self._entries[key]
        path, array = entry.path, entry.array
        del entry
        if path is not None:
            if isinstance(array, np.memmap) and array.mode != "r":
                array.flush()
            del array
            path.unlink(missing_ok=True)

    def close(self) -> None:
        """Release every field of this store; later allocations raise ``RuntimeError``."""
        for key in list(self._entries):
            self.release(key)
        self._closed = True

    def __contains__(self, key: object) -> bool:
        return key in self._entries

    def is_mapped(self, key: str) -> bool:
        """Whether field ``key`` is backed by a file (rather than RAM)."""
        return self._entry(key).path is not None

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    # ----------------------------------------------------------------------------------

    def _entry(self, key: str) -> _Entry:
        try:
            return self._entries[key]
        except KeyError:
            known = ", ".join(sorted(self._entries)) or "(none)"
            raise KeyError(f"no field {key!r} in this store; allocated: {known}") from None

    def _check_open(self) -> None:
        if self._closed:
            raise RuntimeError("this FieldStore is closed")


def _shape(shape: object) -> tuple[int, ...]:
    try:
        dims = tuple(shape)  # type: ignore[arg-type]
    except TypeError:
        raise ValueError(f"shape must be a tuple of integers, got {shape!r}") from None
    for d in dims:
        if isinstance(d, bool) or not isinstance(d, Integral) or d < 0:
            raise ValueError(f"shape must be a tuple of integers >= 0, got {shape!r}")
    return tuple(int(d) for d in dims)
