"""Canonical block grid and order-preserving parallel execution (spec §3.6).

Computations with lateral reach (optical convolution, development) run on a fixed grid of
square blocks anchored at the frame origin. The grid depends only on the frame shape, the
block size and the halo, never on ``RunContext.tile_px``, the worker count or scheduling.
Every block, edge blocks included, has the same padded shape ``block + 2 * halo`` on each
axis, so FFT sizes and in-block reductions have a fixed shape and every block's result is
a pure function of its padded input. :func:`run_blocks` returns results in row-major block
order whatever the number of workers, so cross-block sums made in that order are
bit-identical for any worker count.
"""

import math
import multiprocessing
import os
import threading
from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from dataclasses import dataclass
from numbers import Integral

import numpy as np

__all__ = ["Block", "block_grid", "padded_window", "run_blocks", "worker_environment"]

# Worker processes run single-threaded BLAS so reductions do not depend on thread count
# (spec §3.6).
_BLAS_THREAD_VARIABLES = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)

_ENV_LOCK = threading.Lock()


@dataclass(frozen=True)
class Block:
    """One block of a canonical grid.

    ``inner`` is the part of the frame the block is responsible for, in absolute grid units
    (pixels or chemistry cells); inner parts partition the frame. ``padded_origin`` is the
    absolute position of ``padded[0, 0]`` (negative near the top and left edges) and
    ``padded_shape`` the padded array's shape, the same for every block of a grid.
    """

    index: tuple[int, int]
    inner: tuple[slice, slice]
    padded_origin: tuple[int, int]
    padded_shape: tuple[int, int]

    @property
    def inner_in_padded(self) -> tuple[slice, slice]:
        """``inner`` relative to the padded array: ``padded[inner_in_padded]``."""
        return tuple(  # type: ignore[return-value]
            slice(s.start - o, s.stop - o)
            for s, o in zip(self.inner, self.padded_origin, strict=True)
        )


def block_grid(shape: tuple[int, int], block: int, halo: int) -> tuple[Block, ...]:
    """Return the canonical blocks of a frame of ``shape``, row-major.

    Inner blocks are ``block`` square (cut by the bottom and right frame edges) on a grid
    anchored at ``(0, 0)``; every padded block is ``(block + 2 * halo,) * 2``, even when the
    frame is smaller than one block. ``ValueError`` for ``block < 1``, ``halo < 0`` or a
    frame dimension below 1.
    """
    try:
        h, w = shape
    except (TypeError, ValueError):
        raise ValueError(f"shape must be an (H, W) pair, got {shape!r}") from None
    h, w = _int_at_least("frame height", h, 1), _int_at_least("frame width", w, 1)
    block = _int_at_least("block", block, 1)
    halo = _int_at_least("halo", halo, 0)
    padded = (block + 2 * halo, block + 2 * halo)
    blocks = []
    for i in range(math.ceil(h / block)):
        r0 = i * block
        for j in range(math.ceil(w / block)):
            c0 = j * block
            blocks.append(
                Block(
                    index=(i, j),
                    inner=(slice(r0, min(r0 + block, h)), slice(c0, min(c0 + block, w))),
                    padded_origin=(r0 - halo, c0 - halo),
                    padded_shape=padded,
                )
            )
    return tuple(blocks)


def padded_window(source: np.ndarray, block: Block, fill: float = 0.0) -> np.ndarray:
    """Return a new array of shape ``block.padded_shape + source.shape[2:]`` for ``block``.

    Positions inside ``source`` (shape ``(H, W, ...)``) get its values; positions outside
    the frame get ``fill``. The dtype is ``source``'s.
    """
    source = np.asarray(source)
    if source.ndim < 2:
        raise ValueError(f"source must have shape (H, W, ...), got shape {source.shape}")
    if not isinstance(block, Block):
        raise TypeError(f"block must be a Block, got {type(block).__name__}")
    out = np.full((*block.padded_shape, *source.shape[2:]), fill, dtype=source.dtype)
    src: list[slice] = []
    dst: list[slice] = []
    for origin, size, extent in zip(
        block.padded_origin, block.padded_shape, source.shape[:2], strict=True
    ):
        lo, hi = max(origin, 0), min(origin + size, extent)
        if hi <= lo:
            return out
        src.append(slice(lo, hi))
        dst.append(slice(lo - origin, hi - origin))
    out[tuple(dst)] = source[tuple(src)]
    return out


def worker_environment() -> dict[str, str]:
    """Environment variables that make BLAS single-threaded in worker processes (§3.6)."""
    return dict.fromkeys(_BLAS_THREAD_VARIABLES, "1")


def run_blocks[T](
    fn: Callable[[Block], T], blocks: Sequence[Block], *, workers: int = 1
) -> list[T]:
    """Apply ``fn`` to each block and return the results in the order of ``blocks``.

    With ``workers == 1`` (or a single block) the blocks run in a plain loop in this
    process. Otherwise they run in a pool of ``min(workers, len(blocks))`` freshly spawned
    processes started with :func:`worker_environment`. ``fn`` and its results must then be
    picklable: a module-level function, or a ``functools.partial`` of one with picklable
    arguments. Everything bound into ``fn`` is pickled into every task, so pass large
    arrays as file-backed data (for example a path to reopen) rather than in-memory
    arrays. Exceptions raised by ``fn`` propagate.
    """
    workers = _int_at_least("workers", workers, 1)
    blocks = list(blocks)
    for b in blocks:
        if not isinstance(b, Block):
            raise TypeError(f"blocks must be Block instances, got {type(b).__name__}")
    workers = min(workers, len(blocks))
    if workers <= 1:
        return [fn(b) for b in blocks]
    context = multiprocessing.get_context("spawn")
    with (
        _patched_environment(worker_environment()),
        ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool,
    ):
        return list(pool.map(fn, blocks))


@contextmanager
def _patched_environment(values: dict[str, str]) -> Iterator[None]:
    # Spawned workers inherit os.environ when they start; ProcessPoolExecutor starts them
    # lazily on submit, so the variables stay set until the pool is shut down.
    with _ENV_LOCK:
        saved = {k: os.environ.get(k) for k in values}
        os.environ.update(values)
        try:
            yield
        finally:
            for k, v in saved.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v


def _int_at_least(name: str, value: object, minimum: int) -> int:
    is_int = isinstance(value, Integral) and not isinstance(value, bool)
    if not is_int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}, got {value!r}")
    return int(value)
