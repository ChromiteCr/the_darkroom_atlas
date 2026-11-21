import functools
import os

import numpy as np
import pytest

from darkroom_atlas.blocks import (
    Block,
    block_grid,
    padded_window,
    run_blocks,
    worker_environment,
)

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
    return filtered[block.inner_in_padded] @ np.linspace(
        0.0, 1.0, block.inner[1].stop - block.inner[1].start
    )


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
