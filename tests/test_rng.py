import hashlib

import numpy as np
import pytest

from darkroom_atlas.rng import pixel_uniform, stream_key


def test_stream_key_is_blake2b_integer():
    expected = int(hashlib.blake2b(b"7:grain", digest_size=8).hexdigest(), 16)
    assert stream_key(7, "grain") == expected
    assert 0 <= stream_key(7, "grain") < 2**64


def test_deterministic():
    a = pixel_uniform(3, "grain", (5, 25), (10, 40), 64, (2,))
    b = pixel_uniform(3, "grain", (5, 25), (10, 40), 64, (2,))
    assert a.shape == (20, 30, 2)
    assert a.dtype == np.float64
    assert np.array_equal(a, b)


def test_range_and_mean():
    u = pixel_uniform(0, "s", (0, 200), (0, 300), 300)
    assert u.shape == (200, 300)
    assert u.dtype == np.float64
    assert np.all(u >= 0.0) and np.all(u < 1.0)
    assert abs(u.mean() - 0.5) < 0.01


def test_tile_equals_crop_of_full():
    full = pixel_uniform(7, "s", (0, 50), (0, 80), 80, (3,))
    part = pixel_uniform(7, "s", (10, 30), (20, 60), 80, (3,))
    assert part.shape == (20, 40, 3)
    assert np.array_equal(part, full[10:30, 20:60])


def test_streams_and_seeds_differ():
    base = pixel_uniform(1, "a", (0, 8), (0, 8), 8)
    other_stream = pixel_uniform(1, "b", (0, 8), (0, 8), 8)
    other_seed = pixel_uniform(2, "a", (0, 8), (0, 8), 8)
    assert not np.array_equal(base, other_stream)
    assert not np.array_equal(base, other_seed)


def test_rows_independent():
    u = pixel_uniform(0, "s", (0, 2), (0, 32), 32)
    assert not np.array_equal(u[0], u[1])


def test_rejects_bad_ranges():
    with pytest.raises(ValueError):
        pixel_uniform(0, "s", (0, 4), (0, 81), 80)  # cols beyond frame_width
    with pytest.raises(ValueError):
        pixel_uniform(0, "s", (0, 4), (60, 70), 64)  # cols stop beyond frame_width
    with pytest.raises(ValueError):
        pixel_uniform(0, "s", (4, 0), (0, 8), 8)  # reversed rows
    with pytest.raises(ValueError):
        pixel_uniform(0, "s", (0, 4), (8, 0), 8)  # reversed cols
    with pytest.raises(ValueError):
        pixel_uniform(0, "s", (-1, 4), (0, 8), 8)  # negative row
    with pytest.raises(ValueError):
        pixel_uniform(0, "s", (0, 4), (-2, 8), 8)  # negative col
    with pytest.raises(ValueError):
        pixel_uniform(0, "s", (0, 4), (0, 0), -1)  # negative frame_width


def test_rejects_bad_seed():
    with pytest.raises(ValueError):
        pixel_uniform(-1, "s", (0, 2), (0, 2), 2)
    with pytest.raises(ValueError):
        pixel_uniform(1.5, "s", (0, 2), (0, 2), 2)
    with pytest.raises(ValueError):
        pixel_uniform(True, "s", (0, 2), (0, 2), 2)


def test_numpy_integer_arguments_accepted():
    a = pixel_uniform(np.int64(4), "s", (np.int64(1), 3), (0, np.int32(5)), np.int64(8))
    b = pixel_uniform(4, "s", (1, 3), (0, 5), 8)
    assert np.array_equal(a, b)
