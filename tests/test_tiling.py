import numpy as np
import pytest
from scipy.ndimage import gaussian_filter

from darkroom_atlas.tiling import Tile, iter_tiles, map_tiles


def _coverage(shape_hw, tile_px, overlap_px):
    count = np.zeros(shape_hw, dtype=np.int64)
    for tile in iter_tiles(shape_hw, tile_px, overlap_px):
        count[tile.inner] += 1
    return count


def _blur(a, tile):
    return gaussian_filter(a, sigma=2, truncate=3, mode="nearest")


def _random_frame(shape_hw):
    return np.random.default_rng(20261120).random(shape_hw)


def test_inner_tiles_partition():
    count = _coverage((50, 70), 16, 5)
    assert np.all(count == 1)


def test_outer_clipped_and_contains_inner():
    # Rows 50 = 3 x 16 + 2 and cols 70 = 4 x 16 + 6: the 2-row bottom strip clips the
    # overlap below the tile above it, every other interior edge has room for all of it.
    h, w, overlap = 50, 70, 5
    full_extensions = 0
    for tile in iter_tiles((h, w), 16, overlap):
        for outer, inner, size in zip(tile.outer, tile.inner, (h, w), strict=True):
            assert 0 <= outer.start <= inner.start < inner.stop <= outer.stop <= size
            if inner.start > 0:
                assert inner.start - outer.start == overlap
                full_extensions += 1
            else:
                assert outer.start == 0
            if inner.stop < size:
                assert outer.stop - inner.stop == min(overlap, size - inner.stop)
                full_extensions += outer.stop - inner.stop == overlap
            else:
                assert outer.stop == size
    assert full_extensions > 0


def test_inner_in_outer_consistent():
    x = np.arange(50 * 70 * 2, dtype=np.float64).reshape(50, 70, 2)
    for tile in iter_tiles((50, 70), 16, 5):
        assert np.array_equal(x[tile.outer][tile.inner_in_outer], x[tile.inner])


def test_frame_smaller_than_tile():
    tiles = list(iter_tiles((10, 12), 64, 3))
    assert len(tiles) == 1
    (tile,) = tiles
    full = (slice(0, 10), slice(0, 12))
    assert tile.outer == full
    assert tile.inner == full
    assert tile.inner_in_outer == full
    assert tile == Tile(outer=full, inner=full, inner_in_outer=full)


def test_overlap_larger_than_tile():
    h, w = 30, 41
    assert np.all(_coverage((h, w), 8, 20) == 1)
    x = np.arange(h * w, dtype=np.float64).reshape(h, w)
    for tile in iter_tiles((h, w), 8, 20):
        for outer, inner, size in zip(tile.outer, tile.inner, (h, w), strict=True):
            assert outer.start == max(0, inner.start - 20)
            assert outer.stop == min(size, inner.stop + 20)
        assert np.array_equal(x[tile.outer][tile.inner_in_outer], x[tile.inner])
    assert np.max(np.abs(map_tiles(_blur, x, 8, 20) - _blur(x, None))) < 1e-12


def test_rejects_bad_args():
    x = np.zeros((10, 10))
    for tile_px, overlap_px in ((0, 0), (-1, 2), (16, -1)):
        with pytest.raises(ValueError):
            iter_tiles((10, 10), tile_px, overlap_px)
        with pytest.raises(ValueError):
            map_tiles(lambda a, t: a, x, tile_px, overlap_px)
    with pytest.raises(ValueError):
        iter_tiles((10, 10), 4.0, 1)
    with pytest.raises(ValueError):
        iter_tiles((0, 10), 4, 1)
    with pytest.raises(ValueError):
        map_tiles(lambda a, t: a, np.zeros(10), 4, 1)


def test_map_identity():
    x = _random_frame((45, 61))
    assert np.array_equal(map_tiles(lambda a, t: a, x, 16, 3), x)


def test_map_seamless_blur():
    x = _random_frame((60, 90))
    diff = np.max(np.abs(map_tiles(_blur, x, 16, 6) - _blur(x, None)))
    assert diff < 1e-12


def test_map_insufficient_overlap_differs():
    x = _random_frame((60, 90))
    diff = np.max(np.abs(map_tiles(_blur, x, 16, 1) - _blur(x, None)))
    assert diff > 1e-6


def test_map_trailing_dims():
    x = _random_frame((37, 50, 3))

    def widen(a, tile):
        return np.concatenate([a, a[..., :2] * 2.0], axis=-1)

    out = map_tiles(widen, x, 16, 2)
    assert out.shape == (37, 50, 5)
    assert np.array_equal(out, widen(x, None))


def test_map_passes_read_only_tiles():
    x = _random_frame((20, 20))

    def mutate(a, tile):
        a += 1.0
        return a

    with pytest.raises(ValueError):
        map_tiles(mutate, x, 8, 2)
    assert np.array_equal(x, _random_frame((20, 20)))


def test_map_rejects_wrong_output_shape():
    x = _random_frame((20, 20))
    with pytest.raises(ValueError):
        map_tiles(lambda a, t: a[1:], x, 8, 2)
