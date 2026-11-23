"""Canonical block grid and order-preserving parallel execution (spec §3.6).

Computations with lateral reach (optical convolution, development) run on a fixed grid of
square blocks anchored at the frame origin. The grid depends only on the frame shape, the
block size and the halo, never on ``RunContext.tile_px``, the worker count or scheduling.
Every block, edge blocks included, has the same padded shape ``block + 2 * halo`` on each
axis, so FFT sizes and in-block reductions have a fixed shape and every block's result is
a pure function of its padded input. :func:`run_blocks` returns results in row-major block
order whatever the number of workers, so cross-block sums made in that order are
bit-identical for any worker count.

With several workers, block outputs never travel back through the main process (spec §4.4):
a worker opens its inputs and the preallocated outputs from
:class:`~darkroom_atlas.fieldstore.FieldDescriptor` objects, writes its inner block into the
outputs itself and returns only a small summary, which the main process registers with
:meth:`~darkroom_atlas.state.FrameAssembler.mark`.

The development grid is made of fixed chemistry cells of :data:`CHEMISTRY_CELL_UM` anchored
at the frame origin, independent of the pixel pitch (spec §3.1, §3.6):
:func:`cell_count` gives the cells along one frame axis, :func:`owned_pixels` the pixels a
block of cells owns (those whose centres lie in it, so every pixel belongs to exactly one
block) and :func:`overlapping_pixels` the pixels that overlap a range of cells.
:func:`workers_for_budget` turns the memory budget into a worker count (spec §4.6).
"""

import math
import multiprocessing
import os
import threading
import warnings
from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from dataclasses import dataclass
from numbers import Integral

import numpy as np

__all__ = [
    "CHEMISTRY_CELL_UM",
    "MAIN_PROCESS_GB",
    "WORKER_GB",
    "Block",
    "block_grid",
    "cell_count",
    "overlapping_pixels",
    "owned_pixels",
    "padded_window",
    "run_blocks",
    "worker_environment",
    "workers_for_budget",
]

CHEMISTRY_CELL_UM = 24.0
"""Side of a chemistry cell on the film, um, whatever the pixel pitch (spec §3.1, §3.6)."""

MAIN_PROCESS_GB = 1.5
"""Memory kept back for the main process in the worker formula, GB (spec §4.6)."""

WORKER_GB = 0.5
"""Memory of one development worker process, GB (spec §4.6)."""

# Relative slack when counting cells, so a frame length that is a whole number of cells up
# to rounding (pitch from long edge / pixels) does not get an extra cell.
_CELL_COUNT_RTOL = 1e-9

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
    process. Otherwise they run in one pool of ``min(workers, len(blocks))`` freshly spawned
    processes started with :func:`worker_environment`. ``fn`` and its results must then be
    small and picklable: ``fn`` is a module-level function, or a ``functools.partial`` of
    one, that receives its large inputs and outputs as
    :class:`~darkroom_atlas.fieldstore.FieldDescriptor` objects, opens them with
    :func:`~darkroom_atlas.fieldstore.open_field` (inputs ``"r"``, outputs ``"r+"``),
    writes its inner block into the outputs itself and returns a summary (the inner
    region, block fluxes, counts), never a block array (spec §4.4). Everything bound into
    ``fn`` is pickled into every task. Exceptions raised by ``fn`` propagate.
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


def workers_for_budget(
    memory_budget_gb: float, requested: int | None = None, *, cpu_count: int | None = None
) -> int:
    """Number of worker processes for a memory budget (spec §4.6).

    ``min(cpu_count, max(1, floor((budget - MAIN_PROCESS_GB) / WORKER_GB)))``, 9 for the
    default 6 GB. A ``requested`` count is used when smaller, with a ``RuntimeWarning`` when
    it exceeds the budget's count. ``cpu_count`` defaults to ``os.cpu_count()``.
    """
    if (
        isinstance(memory_budget_gb, bool)
        or not isinstance(memory_budget_gb, (int, float))
        or not math.isfinite(memory_budget_gb)
        or memory_budget_gb <= 0.0
    ):
        raise ValueError(
            f"memory_budget_gb must be a positive, finite number, got {memory_budget_gb!r}"
        )
    if cpu_count is None:
        cpus = os.cpu_count() or 1
    else:
        cpus = _int_at_least("cpu_count", cpu_count, 1)
    # Small epsilon so an exact multiple (6 GB -> 9) is not lost to binary rounding.
    budget = max(1, math.floor((memory_budget_gb - MAIN_PROCESS_GB) / WORKER_GB + 1e-9))
    allowed = min(cpus, budget)
    if requested is None:
        return allowed
    requested = _int_at_least("requested", requested, 1)
    if requested > budget:
        warnings.warn(
            f"{requested} workers exceed the {memory_budget_gb:g} GB memory budget "
            f"({budget} workers of {WORKER_GB} GB each); using {min(requested, allowed)}",
            RuntimeWarning,
            stacklevel=2,
        )
    return min(requested, allowed)


def cell_count(n_px: int, pitch_um: float, cell_um: float = CHEMISTRY_CELL_UM) -> int:
    """Chemistry cells covering ``n_px`` pixels of pitch ``pitch_um`` along one axis.

    ``ceil(n_px * pitch / cell)``: the grid starts at the frame origin, and a last cell
    reaching past the frame edge counts (its outside area is unexposed film, spec §3.6).
    """
    n_px = _int_at_least("n_px", n_px, 1)
    pitch_um, cell_um = _length("pitch_um", pitch_um), _length("cell_um", cell_um)
    cells = n_px * pitch_um / cell_um
    return max(1, math.ceil(cells - _CELL_COUNT_RTOL * cells))


def owned_pixels(
    cells: slice, pitch_um: float, n_px: int, cell_um: float = CHEMISTRY_CELL_UM
) -> slice:
    """Pixels along one axis owned by the cells ``[i0, i1)`` (spec §3.6).

    A block owns the pixels whose centres lie in ``[i0 cell, i1 cell)``, i.e. pixels
    ``ceil(i0 cell / p - 1/2)`` up to ``ceil(i1 cell / p - 1/2) - 1``, clipped to the frame
    of ``n_px`` pixels. Adjacent blocks evaluate the same expression at their shared
    boundary, so the owned ranges of a row of blocks partition the frame; the block ending
    at the last cell owns up to the frame edge. The result may be empty when a block lies
    beyond the last pixel centre.
    """
    i0, i1 = _cell_slice(cells)
    pitch_um, cell_um = _length("pitch_um", pitch_um), _length("cell_um", cell_um)
    n_px = _int_at_least("n_px", n_px, 1)

    def first_centre(i: int) -> int:
        return min(max(math.ceil(i * cell_um / pitch_um - 0.5), 0), n_px)

    start = first_centre(i0)
    stop = n_px if i1 >= cell_count(n_px, pitch_um, cell_um) else first_centre(i1)
    return slice(start, max(start, stop))


def overlapping_pixels(
    cells: slice, pitch_um: float, cell_um: float = CHEMISTRY_CELL_UM
) -> tuple[int, int]:
    """Pixel range ``[j0, j1)`` overlapping cells ``[i0, i1)`` along one axis, unclipped.

    ``j0 = floor(i0 cell / p)``, ``j1 = ceil(i1 cell / p)``. Cells may be negative (a padded
    block near the frame edge), and so may the result; positions outside the frame are
    unexposed film for the caller to fill (spec §3.6).
    """
    if not isinstance(cells, slice) or cells.step not in (None, 1):
        raise ValueError(f"cells must be a slice with step 1, got {cells!r}")
    i0, i1 = cells.start, cells.stop
    for bound in (i0, i1):
        if isinstance(bound, bool) or not isinstance(bound, Integral):
            raise TypeError(f"cells must have integer bounds, got {cells!r}")
    if i1 <= i0:
        raise ValueError(f"cells {cells!r} is empty")
    pitch_um, cell_um = _length("pitch_um", pitch_um), _length("cell_um", cell_um)
    return math.floor(i0 * cell_um / pitch_um), math.ceil(i1 * cell_um / pitch_um)


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


def _length(name: str, value: object) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0.0
    ):
        raise ValueError(f"{name} must be a positive, finite number, got {value!r}")
    return float(value)


def _cell_slice(cells: object) -> tuple[int, int]:
    if not isinstance(cells, slice) or cells.step not in (None, 1):
        raise ValueError(f"cells must be a slice with step 1, got {cells!r}")
    i0, i1 = cells.start, cells.stop
    for bound in (i0, i1):
        if isinstance(bound, bool) or not isinstance(bound, Integral) or bound < 0:
            raise ValueError(f"cells must have integer bounds >= 0, got {cells!r}")
    if i1 <= i0:
        raise ValueError(f"cells {cells!r} is empty")
    return int(i0), int(i1)
