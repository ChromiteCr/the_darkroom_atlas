"""Deterministic random numbers keyed by seed, stream name and absolute pixel position.

Every random number in the simulation comes from this module. A value depends only on
``(seed, stream, row, column, trailing index)``, never on how the frame was tiled or in
which order tiles were computed, so a tiled run reproduces an untiled run bit for bit.

Each absolute row ``r`` gets its own PCG64 bit generator seeded with
``SeedSequence([seed, stream_key(seed, stream), r])``. Only the raw 64-bit output of the
bit generator is used (never the distribution methods of ``numpy.random.Generator``), so
the values do not change between NumPy versions.
"""

import hashlib
from math import prod
from numbers import Integral

import numpy as np

__all__ = ["pixel_uniform", "stream_key"]

_MANTISSA_SHIFT = np.uint64(11)
_INV_2_53 = 2.0**-53


def stream_key(seed: int, stream: str) -> int:
    """Return the 64-bit integer key of ``stream`` under ``seed``.

    The key is the integer value of ``blake2b(f"{seed}:{stream}".encode(), digest_size=8)``,
    which is stable across processes and platforms (unlike the salted built-in ``hash``).
    """
    digest = hashlib.blake2b(f"{seed}:{stream}".encode(), digest_size=8)
    return int(digest.hexdigest(), 16)


def pixel_uniform(
    seed: int,
    stream: str,
    rows: tuple[int, int],
    cols: tuple[int, int],
    frame_width: int,
    trailing: tuple[int, ...] = (),
) -> np.ndarray:
    """Return uniform ``float64`` samples in ``[0, 1)`` for a window of the frame.

    ``rows`` and ``cols`` are half-open ranges in absolute frame coordinates, and
    ``frame_width`` is the width of the whole frame. The result has shape
    ``(rows[1] - rows[0], cols[1] - cols[0], *trailing)``; the sample for a given absolute
    pixel and trailing index is the same whatever window it is requested through.

    Raises ``ValueError`` if ``seed`` is not a non-negative integer, if a range is reversed
    or negative, or if ``cols`` reaches beyond ``frame_width``.
    """
    seed = _non_negative_int("seed", seed)
    frame_width = _non_negative_int("frame_width", frame_width)
    if frame_width < 1:
        raise ValueError(f"frame_width must be at least 1, got {frame_width}")
    r0, r1 = _half_open("rows", rows)
    c0, c1 = _half_open("cols", cols)
    if c1 > frame_width:
        raise ValueError(f"cols {(c0, c1)} reach beyond frame_width {frame_width}")
    trailing = tuple(_non_negative_int("trailing", n) for n in trailing)

    key = stream_key(seed, stream)
    per_pixel = prod(trailing)
    out = np.empty((r1 - r0, c1 - c0, *trailing), dtype=np.float64)
    for i, r in enumerate(range(r0, r1)):
        bitgen = np.random.PCG64(np.random.SeedSequence([seed, key, r]))
        raw = bitgen.random_raw(frame_width * per_pixel)
        row = ((raw >> _MANTISSA_SHIFT) * _INV_2_53).reshape(frame_width, *trailing)
        out[i] = row[c0:c1]
    return out


def _non_negative_int(name: str, value: object) -> int:
    is_int = isinstance(value, Integral) and not isinstance(value, bool)
    if not is_int or value < 0:
        raise ValueError(f"{name} must be a non-negative integer, got {value!r}")
    return int(value)


def _half_open(name: str, bounds: tuple[int, int]) -> tuple[int, int]:
    try:
        start, stop = bounds
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a (start, stop) pair, got {bounds!r}") from None
    start = _non_negative_int(f"{name} start", start)
    stop = _non_negative_int(f"{name} stop", stop)
    if stop < start:
        raise ValueError(f"{name} range {(start, stop)} is reversed")
    return start, stop
