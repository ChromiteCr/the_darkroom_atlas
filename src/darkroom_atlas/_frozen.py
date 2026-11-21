"""Read-only containers for the immutable state objects (``Light``, ``EmulsionState``).

:class:`FrozenDict` is a ``dict`` that refuses changes, so metadata stays JSON- and
pydantic-serialisable (unlike ``types.MappingProxyType``) while it cannot be altered in
place. :func:`readonly_array` stores arrays so that nobody, including the caller who passed
them in, can change them afterwards. Read-only memory maps (``np.memmap`` opened with
``mode="r"``, as :class:`darkroom_atlas.fieldstore.FieldStore` freezes large fields) and
their views are shared without copying (spec §4.3).
"""

from collections.abc import Mapping
from typing import Any, NoReturn

import numpy as np

__all__ = ["FrozenDict", "freeze_mapping", "readonly_array"]


class FrozenDict(dict):
    """A ``dict`` that cannot be changed after construction.

    Reading, iterating, comparing, copying, pickling and JSON encoding work as for ``dict``;
    every mutating method raises ``TypeError``. ``a | b`` returns a new plain ``dict``.
    """

    __slots__ = ()

    def _read_only(self, *args: Any, **kwargs: Any) -> NoReturn:
        raise TypeError(f"{type(self).__name__} is read-only; build a new mapping instead")

    __setitem__ = _read_only
    __delitem__ = _read_only
    __ior__ = _read_only
    clear = _read_only
    pop = _read_only
    popitem = _read_only
    setdefault = _read_only
    update = _read_only

    def __reduce__(self):
        return (type(self), (dict(self),))


def freeze_mapping(mapping: Mapping[str, Any]) -> FrozenDict:
    """Return a :class:`FrozenDict` copy of ``mapping``, freezing nested mappings too."""
    if not isinstance(mapping, Mapping):
        raise TypeError(f"expected a mapping, got {type(mapping).__name__}")
    return FrozenDict(
        (key, value if isinstance(value, FrozenDict) else _freeze_value(value))
        for key, value in mapping.items()
    )


def _freeze_value(value: Any) -> Any:
    return freeze_mapping(value) if isinstance(value, Mapping) else value


def readonly_array(value: Any, dtype: Any = None) -> np.ndarray:
    """Return ``value`` as a NumPy array that cannot be written through any reference.

    An ndarray (or ``np.memmap``) that is already read-only all the way down its ``base``
    chain is returned as is, so read-only fields and their views are shared, not copied. A
    chain may end in a ``mode="r"`` memory map: such a mapping cannot be written through
    NumPy. Anything else, including writable maps (``"r+"``, ``"w+"``, ``"c"``), is copied
    into a new array, optionally converted to ``dtype``, which is then marked read-only; the
    caller's array is never modified.
    """
    shareable = type(value) in (np.ndarray, np.memmap) and _read_only_chain(value)
    if shareable and (dtype is None or value.dtype == np.dtype(dtype)):
        return value
    array = np.array(value, dtype=dtype, copy=True)
    array.flags.writeable = False
    return array


def _read_only_chain(array: np.ndarray) -> bool:
    current: Any = array
    while isinstance(current, np.ndarray):
        if current.flags.writeable:
            return False
        if isinstance(current, np.memmap) and current.mode == "r":
            # A read-only mapping (mmap ACCESS_READ) cannot be written through any view.
            return True
        current = current.base
    # Memory owned by something other than an ndarray (a buffer, a writable memory map)
    # may be writable through that owner, so only a chain that ends in an owning array
    # counts.
    return current is None
