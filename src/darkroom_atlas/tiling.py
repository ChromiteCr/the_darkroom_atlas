"""Tile iteration and seamless stitching for frame-sized arrays.

A frame of shape ``(H, W, ...)`` is cut into *inner* tiles on a ``tile_px`` grid; the inner
tiles partition the frame. Each inner tile is grown by ``overlap_px`` on every side and
clipped to the frame to give its *outer* tile, which is what a tiled computation actually
sees. As long as ``overlap_px`` covers the lateral reach of the computation, the inner part
of every outer result equals the corresponding part of the untiled result, so stitching
the inner parts back together reproduces the untiled computation without seams.
"""

from collections.abc import Callable, Iterator
from dataclasses import dataclass
from numbers import Integral

import numpy as np

__all__ = ["Tile", "iter_tiles", "map_tiles"]


@dataclass(frozen=True)
class Tile:
    """One tile of a frame, as ``(row_slice, col_slice)`` pairs.

    ``outer`` and ``inner`` are in absolute frame coordinates; ``inner_in_outer`` locates the
    inner tile inside an array cut out with ``outer``, so that
    ``array[outer][inner_in_outer]`` is ``array[inner]``. Every slice has concrete integer
    ``start`` and ``stop`` and a step of ``None``.
    """

    outer: tuple[slice, slice]
    inner: tuple[slice, slice]
    inner_in_outer: tuple[slice, slice]


def iter_tiles(shape_hw: tuple[int, int], tile_px: int, overlap_px: int) -> Iterator[Tile]:
    """Iterate over the tiles of a frame of shape ``shape_hw``, row-major.

    Inner tiles are ``tile_px`` square (smaller along the bottom and right edges) and
    partition the frame; outer tiles extend them by ``overlap_px`` on each side, clipped to
    the frame. Arguments are checked immediately: ``ValueError`` if ``tile_px < 1``,
    ``overlap_px < 0`` or a frame dimension is below 1.
    """
    h, w = _frame_shape(shape_hw)
    tile_px = _int_at_least("tile_px", tile_px, 1)
    overlap_px = _int_at_least("overlap_px", overlap_px, 0)
    return _generate_tiles(h, w, tile_px, overlap_px)


def map_tiles(
    fn: Callable[[np.ndarray, Tile], np.ndarray],
    array: np.ndarray,
    tile_px: int,
    overlap_px: int,
) -> np.ndarray:
    """Apply ``fn`` tile by tile to ``array`` (shape ``(H, W, ...)``) and stitch the results.

    ``fn(outer_block, tile)`` receives a read-only view ``array[tile.outer]`` and must return
    an array whose first two dimensions match that block; its trailing dimensions may differ
    from the input's. Only the ``tile.inner_in_outer`` part of each result is kept. The
    trailing dimensions and dtype of the output are taken from the first tile's result.

    Raises ``ValueError`` for bad tiling arguments, an input with fewer than two dimensions,
    or a result whose shape does not fit.
    """
    array = np.asarray(array)
    if array.ndim < 2:
        raise ValueError(f"array must have shape (H, W, ...), got shape {array.shape}")
    h, w = array.shape[:2]
    out: np.ndarray | None = None
    for tile in iter_tiles((h, w), tile_px, overlap_px):
        block = array[tile.outer]
        block.flags.writeable = False
        result = np.asarray(fn(block, tile))
        if result.shape[:2] != block.shape[:2]:
            raise ValueError(
                f"fn returned shape {result.shape} for a tile of shape {block.shape}; "
                "the first two dimensions must match"
            )
        if out is None:
            out = np.empty((h, w, *result.shape[2:]), dtype=result.dtype)
        elif result.shape[2:] != out.shape[2:]:
            raise ValueError(
                f"fn returned trailing shape {result.shape[2:]} for tile {tile.inner}, "
                f"but the first tile gave {out.shape[2:]}"
            )
        out[tile.inner] = result[tile.inner_in_outer]
    assert out is not None  # iter_tiles yields at least one tile for a non-empty frame
    return out


def _generate_tiles(h: int, w: int, tile_px: int, overlap_px: int) -> Iterator[Tile]:
    for r0 in range(0, h, tile_px):
        rows = _axis_slices(r0, min(r0 + tile_px, h), h, overlap_px)
        for c0 in range(0, w, tile_px):
            cols = _axis_slices(c0, min(c0 + tile_px, w), w, overlap_px)
            yield Tile(
                outer=(rows[0], cols[0]),
                inner=(rows[1], cols[1]),
                inner_in_outer=(rows[2], cols[2]),
            )


def _axis_slices(start: int, stop: int, size: int, overlap_px: int) -> tuple[slice, slice, slice]:
    outer_start = max(0, start - overlap_px)
    outer_stop = min(size, stop + overlap_px)
    return (
        slice(outer_start, outer_stop),
        slice(start, stop),
        slice(start - outer_start, stop - outer_start),
    )


def _frame_shape(shape_hw: tuple[int, int]) -> tuple[int, int]:
    try:
        h, w = shape_hw
    except (TypeError, ValueError):
        raise ValueError(f"shape_hw must be an (H, W) pair, got {shape_hw!r}") from None
    return _int_at_least("frame height", h, 1), _int_at_least("frame width", w, 1)


def _int_at_least(name: str, value: object, minimum: int) -> int:
    is_int = isinstance(value, Integral) and not isinstance(value, bool)
    if not is_int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}, got {value!r}")
    return int(value)
