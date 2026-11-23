import functools
import math
import os
import pickle

import numpy as np
import pytest

from darkroom_atlas.blocks import (
    CHEMISTRY_CELL_UM,
    Block,
    block_grid,
    cell_count,
    overlapping_pixels,
    owned_pixels,
    padded_window,
    run_blocks,
    worker_environment,
    workers_for_budget,
)
from darkroom_atlas.fieldstore import FieldDescriptor, FieldStore, open_field
from darkroom_atlas.state import FieldSpec, FrameAssembler, FrameGeometry

BLAS_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")


@pytest.mark.parametrize(
    ("shape", "block", "halo"), [((50, 70), 16, 3), ((16, 32), 16, 0), ((5, 9), 16, 4)]
)
def test_block_grid_tiles_frame_once(shape, block, halo):
    blocks = block_grid(shape, block, halo)
    count = np.zeros(shape, dtype=int)
    for b in blocks:
        count[b.inner] += 1
    assert (count == 1).all()
    assert {b.padded_shape for b in blocks} == {(block + 2 * halo, block + 2 * halo)}
    for b in blocks:
        assert b.padded_origin == (b.inner[0].start - halo, b.inner[1].start - halo)
        assert b.inner[0].start == b.index[0] * block and b.inner[1].start == b.index[1] * block
    # Row-major order.
    assert [b.index for b in blocks] == sorted(b.index for b in blocks)


def test_block_grid_anchored_at_origin():
    # The grid depends only on frame shape, block and halo: no other tiling parameter.
    a = block_grid((100, 130), 32, 5)
    assert a == block_grid((100, 130), 32, 5)
    assert [b.inner for b in a[:5]] == [
        (slice(0, 32), slice(c, min(c + 32, 130))) for c in (0, 32, 64, 96, 128)
    ]


def test_block_grid_rejects_bad_arguments():
    for args in [((0, 5), 4, 0), ((5, 5), 0, 0), ((5, 5), 4, -1), ((5,), 4, 0), ((5, 5), 2.0, 0)]:
        with pytest.raises(ValueError):
            block_grid(*args)


def test_padded_window_fills_outside():
    source = np.arange(7 * 9 * 2, dtype=float).reshape(7, 9, 2)
    for b in block_grid((7, 9), 4, 2):
        window = padded_window(source, b, fill=-1.0)
        assert window.shape == (8, 8, 2) and window.dtype == source.dtype
        np.testing.assert_array_equal(window[b.inner_in_padded], source[b.inner])
        r0, c0 = b.padded_origin
        for i in range(8):
            for j in range(8):
                r, c = r0 + i, c0 + j
                expected = source[r, c] if 0 <= r < 7 and 0 <= c < 9 else -1.0
                np.testing.assert_array_equal(window[i, j], expected)
        window[...] = 0.0  # a new array: writing it leaves source alone
    assert source[0, 0, 1] == 1.0


def test_padded_window_halo_larger_than_frame():
    source = np.ones((2, 3))
    (b,) = block_grid((2, 3), 4, 6)
    window = padded_window(source, b)
    assert window.shape == (16, 16)
    assert window.sum() == 6.0


def test_padded_window_rejects_bad_input():
    (b,) = block_grid((2, 2), 2, 0)
    with pytest.raises(ValueError):
        padded_window(np.ones(4), b)
    with pytest.raises(TypeError):
        padded_window(np.ones((2, 2)), "block")


def _spectral_block(source: np.ndarray, block: Block) -> np.ndarray:
    """A fixed-size FFT filter of a padded block, cropped to the inner part."""
    window = padded_window(source, block)
    filtered = np.fft.irfft2(np.fft.rfft2(window) * 0.5, s=window.shape)
    weights = np.linspace(0.0, 1.0, block.inner[1].stop - block.inner[1].start)
    # Element-wise product and np.sum, not "@" (spec §3.6 coding rule).
    return np.sum(filtered[block.inner_in_padded] * weights, axis=1)


def _blas_env(_block: Block) -> dict[str, str | None]:
    return {k: os.environ.get(k) for k in BLAS_VARS}


def test_worker_environment_single_threaded():
    assert worker_environment() == dict.fromkeys(BLAS_VARS, "1")


def test_run_blocks_order_and_identity_across_workers():
    source = np.random.default_rng(1).random((90, 70))
    blocks = block_grid(source.shape, 32, 4)
    fn = functools.partial(_spectral_block, source)
    reference = run_blocks(fn, blocks, workers=1)
    assert len(reference) == len(blocks)
    for workers in (2, 3):
        results = run_blocks(fn, blocks, workers=workers)
        assert [r.tobytes() for r in results] == [r.tobytes() for r in reference]


def test_run_blocks_workers_get_single_threaded_blas():
    before = dict(os.environ)
    blocks = block_grid((4, 4), 2, 0)
    envs = run_blocks(_blas_env, blocks, workers=2)
    assert envs == [dict.fromkeys(BLAS_VARS, "1")] * len(blocks)
    assert dict(os.environ) == before  # the parent environment is restored


def test_run_blocks_edge_cases():
    blocks = block_grid((4, 4), 2, 0)
    assert run_blocks(_blas_env, [], workers=4) == []
    assert run_blocks(lambda b: b.index, blocks) == [b.index for b in blocks]
    with pytest.raises(ValueError):
        run_blocks(_blas_env, blocks, workers=0)
    with pytest.raises(TypeError):
        run_blocks(_blas_env, [object()])  # type: ignore[list-item]


def _write_block_directly(
    source: FieldDescriptor, outputs: dict[str, FieldDescriptor], block: Block
) -> dict:
    """Worker side of spec §4.4: read the input, write the inner block, return a summary."""
    data = open_field(source, "r")
    window = padded_window(data, block)
    smoothed = np.fft.irfft2(np.fft.rfft2(window, axes=(0, 1)) * 0.5, s=window.shape[:2],
                             axes=(0, 1))  # fmt: skip
    inner = smoothed[block.inner_in_padded]
    out = open_field(outputs["x/out"], "r+")
    out[block.inner] = inner
    out.flush()
    del out
    return {"index": block.index, "inner": block.inner, "total": float(np.sum(inner))}


def test_run_blocks_workers_write_directly(tmp_path):
    # Spec §10 D1a: with workers = 3 the workers write the outputs themselves through
    # FieldDescriptor ("r+"); the main process gets only small summaries; bit-identical to
    # workers = 1.
    source = np.random.default_rng(7).random((90, 70, 3))
    geometry = FrameGeometry((90, 70), 6.0)
    blocks = block_grid(source.shape[:2], 32, 5)
    results = {}
    with FieldStore(tmp_path, ram_threshold_bytes=0) as store:
        src = store.allocate("in:x/src", source.shape)
        src[...] = source
        store.freeze("in:x/src")
        for workers in (1, 3):
            asm = FrameAssembler(
                geometry, [FieldSpec("x/out", (3,))], store, name_prefix=f"w{workers}:"
            )
            outputs = asm.descriptors()
            fn = functools.partial(_write_block_directly, store.descriptor("in:x/src"), outputs)
            summaries = run_blocks(fn, blocks, workers=workers)
            assert [s["index"] for s in summaries] == [b.index for b in blocks]
            for summary in summaries:
                assert len(pickle.dumps(summary)) < 64 * 1024
                assert not any(isinstance(v, np.ndarray) for v in summary.values())
                asm.mark(summary["inner"])
            results[workers] = (np.array(asm.finish()["x/out"]), summaries)
    reference, ref_summaries = results[1]
    assert np.array_equal(results[3][0], reference)
    assert [s["total"] for s in results[3][1]] == [s["total"] for s in ref_summaries]
    assert np.any(reference != 0.0)


@pytest.mark.parametrize("pitch", [6.0, 9.0, 12.0, 30.0, 48.0, 36000.0 / 7000])
def test_owned_pixels_partition_frame(pitch):
    # Spec §3.6: every pixel belongs to exactly one block of cells, the one its centre lies
    # in; the cell grid is anchored at the frame origin with 24 um cells (spec §3.1).
    n_px = 731
    n_cells = cell_count(n_px, pitch)
    assert n_cells == math.ceil(round(n_px * pitch / CHEMISTRY_CELL_UM, 6))
    for block in (1, 7, 256):
        count = np.zeros(n_px, dtype=int)
        for start in range(0, n_cells, block):
            owned = owned_pixels(slice(start, min(start + block, n_cells)), pitch, n_px)
            count[owned] += 1
            for j in range(owned.start, owned.stop):
                centre = (j + 0.5) * pitch
                lo, hi = start * CHEMISTRY_CELL_UM, min(start + block, n_cells) * 24.0
                # Centres exactly on a face belong to the block that starts there.
                assert lo - 1e-9 <= centre < hi + 1e-9 or start + block >= n_cells
        assert (count == 1).all()


def test_cell_grid_counts_and_overlaps():
    # 135 at 6000 px (p = 6 um): 36 mm = 1500 cells; 120 at 56 mm: 2333.3 -> 2334 cells.
    assert cell_count(6000, 36000.0 / 6000) == 1500
    assert cell_count(7000, 36000.0 / 7000) == 1500
    assert cell_count(4000, 56000.0 / 4000) == 2334
    assert cell_count(1, 6.0) == 1
    # p = 48 um: one pixel spans two cells; centres land on faces (pixel 0 at 24 um).
    assert owned_pixels(slice(0, 1), 48.0, 10) == slice(0, 0)
    assert owned_pixels(slice(1, 3), 48.0, 10) == slice(0, 1)
    assert overlapping_pixels(slice(0, 1), 48.0) == (0, 1)
    assert overlapping_pixels(slice(-3, 2), 9.0) == (-8, 6)
    assert overlapping_pixels(slice(0, 1), 6.0) == (0, 4)
    with pytest.raises(ValueError):
        owned_pixels(slice(2, 2), 6.0, 10)
    with pytest.raises(ValueError):
        owned_pixels(slice(-1, 2), 6.0, 10)
    with pytest.raises(ValueError):
        cell_count(10, 0.0)
    with pytest.raises(TypeError):
        overlapping_pixels(slice(0.5, 2), 6.0)
    with pytest.raises(ValueError):
        overlapping_pixels(slice(0, 4, 2), 6.0)


def test_workers_for_budget():
    # Spec §4.6: min(cpu, max(1, floor((budget - 1.5) / 0.5))), 9 for 6 GB.
    assert workers_for_budget(6.0, cpu_count=64) == 9
    assert workers_for_budget(6.0, cpu_count=4) == 4
    assert workers_for_budget(1.0, cpu_count=8) == 1
    assert workers_for_budget(3.4, cpu_count=8) == 3
    assert workers_for_budget(6.0, 2, cpu_count=64) == 2
    with pytest.warns(RuntimeWarning, match="exceed"):
        assert workers_for_budget(6.0, 12, cpu_count=64) == 9
    assert 1 <= workers_for_budget(6.0) <= max(1, os.cpu_count() or 1)
    for bad in (0.0, -1.0, float("nan"), True):
        with pytest.raises(ValueError):
            workers_for_budget(bad)
    with pytest.raises(ValueError):
        workers_for_budget(6.0, 0)
