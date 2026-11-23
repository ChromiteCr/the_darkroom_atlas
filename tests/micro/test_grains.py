"""Microscope-mode grain placement (spec D1 2.7.2, step 1; 10 D1f)."""

import dataclasses
import math

import numpy as np
import pytest

from darkroom_atlas.emulsion.classes import SizeClasses
from darkroom_atlas.micro import GrainSet, place_grains
from darkroom_atlas.micro.grains import MICRO_CELL_UM


def _ro(values) -> np.ndarray:
    array = np.array(values, dtype=np.float64)
    array.flags.writeable = False
    return array


def make_classes(number_per_um2, ecd_um=None) -> SizeClasses:
    """Size classes with given surface densities; only ``number_per_um2`` matters here."""
    nu = np.atleast_1d(np.asarray(number_per_um2, dtype=np.float64))
    d = np.full(nu.shape, 0.5) if ecd_um is None else np.asarray(ecd_um, dtype=np.float64)
    return SizeClasses(
        morphology="conventional",
        ecd_um=_ro(d),
        thickness_um=_ro(d),
        volume_um3=_ro(np.pi * d**3 / 6.0),
        surface_um2=_ro(np.pi * d**2),
        projected_um2=_ro(np.pi * d**2 / 4.0),
        number_fraction=_ro(nu / nu.sum()),
        number_per_um2=_ro(nu),
        silver_g_per_grain=_ro(np.full(nu.shape, 1e-12)),
        iodide_fraction=0.0,
        dry_thickness_um=5.0,
        silver_coverage_g_m2=1.0,
    )


CLASSES = make_classes([0.4, 1.5, 0.05])


def _as_tuples(g: GrainSet, mask: np.ndarray) -> list[tuple]:
    return list(
        zip(
            g.x_um[mask].tolist(),
            g.y_um[mask].tolist(),
            g.size_class[mask].tolist(),
            g.depth_index[mask].tolist(),
            strict=True,
        )
    )


def test_micro_deterministic():
    a = place_grains(CLASSES, (120.0, 90.0), 8, seed=5, origin_um=(36.0, 20.0))
    b = place_grains(CLASSES, (120.0, 90.0), 8, seed=5, origin_um=(36.0, 20.0))
    for field in dataclasses.fields(GrainSet):
        np.testing.assert_array_equal(getattr(a, field.name), getattr(b, field.name))
    assert a.count > 1000

    # A different origin and size: the shared micro-cells hold the same grains, bit for bit.
    c = place_grains(CLASSES, (101.0, 77.0), 8, seed=5, origin_um=(81.0, 51.0))
    lo_x, hi_x = 81.0, 156.0
    lo_y, hi_y = 51.0, 110.0

    def overlap(g):
        return (g.x_um >= lo_x) & (g.x_um < hi_x) & (g.y_um >= lo_y) & (g.y_um < hi_y)

    shared_a, shared_c = _as_tuples(a, overlap(a)), _as_tuples(c, overlap(c))
    assert len(shared_a) > 500
    assert shared_a == shared_c

    # Another seed gives other grains.
    d = place_grains(CLASSES, (120.0, 90.0), 8, seed=6, origin_um=(36.0, 20.0))
    assert _as_tuples(d, overlap(d)) != shared_a


def test_patch_is_a_union_of_sub_patches():
    whole = place_grains(CLASSES, (80.0, 80.0), 4, seed=1, origin_um=(10.0, 10.0))
    parts = [
        place_grains(CLASSES, (40.0, 80.0), 4, seed=1, origin_um=(10.0, 10.0)),
        place_grains(CLASSES, (40.0, 80.0), 4, seed=1, origin_um=(50.0, 10.0)),
    ]
    joined = sorted(t for p in parts for t in _as_tuples(p, np.ones(p.count, bool)))
    assert joined == sorted(_as_tuples(whole, np.ones(whole.count, bool)))


def test_grains_lie_in_patch_with_empty_development():
    w, h, x0, y0, n_z = 50.5, 33.3, 7.25, 101.0, 6
    g = place_grains(CLASSES, (w, h), n_z, seed=3, origin_um=(x0, y0))
    assert np.all((g.x_um >= x0) & (g.x_um < x0 + w))
    assert np.all((g.y_um >= y0) & (g.y_um < y0 + h))
    assert g.depth_index.min() >= 0 and g.depth_index.max() <= n_z - 1
    assert set(np.unique(g.size_class)) <= {0, 1, 2}
    np.testing.assert_array_equal(g.cluster_size, 0)
    assert np.all(np.isinf(g.init_time_s))
    np.testing.assert_array_equal(g.developed_fraction, 0.0)
    assert g.depth_index.dtype == np.int64 and g.size_class.dtype == np.int64
    assert not g.x_um.flags.writeable


def test_counts_positions_and_depths_are_poisson_uniform():
    side, n_z = 400.0, 8
    g = place_grains(CLASSES, (side, side), n_z, seed=11)
    area = side * side
    for c, nu in enumerate(CLASSES.number_per_um2):
        n = np.count_nonzero(g.size_class == c)
        mean = nu * area
        assert abs(n - mean) <= 4.0 * math.sqrt(mean)
    # Uniform within cells and over depth.
    frac_x = np.mod(g.x_um, MICRO_CELL_UM) / MICRO_CELL_UM
    counts, _ = np.histogram(frac_x, bins=8, range=(0.0, 1.0))
    expected = g.count / 8
    assert np.all(np.abs(counts - expected) <= 5.0 * math.sqrt(expected))
    depth_counts = np.bincount(g.depth_index, minlength=n_z)
    expected = g.count / n_z
    assert np.all(np.abs(depth_counts - expected) <= 5.0 * math.sqrt(expected))
    # Counts per cell are Poisson: variance equals mean.
    cells = (np.floor(g.y_um / MICRO_CELL_UM) * 1000 + np.floor(g.x_um / MICRO_CELL_UM))[
        g.size_class == 1
    ]
    per_cell = np.bincount(np.unique(cells, return_inverse=True)[1])
    per_cell = np.concatenate([per_cell, np.zeros(10_000 - per_cell.size)])
    mean = CLASSES.number_per_um2[1] * MICRO_CELL_UM**2
    assert per_cell.var() == pytest.approx(mean, rel=0.05)


def test_large_cell_means_use_cornish_fisher():
    # 40 grains per um2 gives a mean of 640 per cell, above the exact inverse-CDF range.
    g = place_grains(make_classes([40.0]), (40.0, 40.0), 2, seed=0)
    assert abs(g.count - 64_000) <= 4.0 * math.sqrt(64_000)


def test_empty_class_and_order():
    g = place_grains(make_classes([0.0, 0.3]), (60.0, 60.0), 3, seed=2)
    assert g.count > 0 and np.all(g.size_class == 1)
    # Grains come in order of cell row, then cell column.
    key = np.floor(g.y_um / MICRO_CELL_UM) * 1e6 + np.floor(g.x_um / MICRO_CELL_UM)
    assert np.all(np.diff(key) >= 0)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"patch_um": (0.0, 10.0)}, "patch_um"),
        ({"patch_um": (10.0, -1.0)}, "patch_um"),
        ({"patch_um": (2001.0, 10.0)}, "exceed"),
        ({"patch_um": (10.0, float("nan"))}, "patch_um"),
        ({"patch_um": 10.0}, "patch_um"),
        ({"origin_um": (-1.0, 0.0)}, "origin_um"),
        ({"origin_um": (0.0, float("inf"))}, "origin_um"),
        ({"n_z": 0}, "n_z"),
        ({"n_z": 2.0}, "n_z"),
        ({"n_z": True}, "n_z"),
        ({"seed": -1}, "seed"),
        ({"seed": 1.5}, "seed"),
    ],
)
def test_place_grains_rejects_bad_input(kwargs, match):
    args = {"patch_um": (10.0, 10.0), "n_z": 4, "seed": 0, "origin_um": (0.0, 0.0)} | kwargs
    with pytest.raises(ValueError, match=match):
        place_grains(CLASSES, **args)


def test_place_grains_rejects_bad_densities():
    with pytest.raises(ValueError, match="number_per_um2"):
        place_grains(make_classes([-0.1, 1.0]), (10.0, 10.0), 4, seed=0)


def _grainset(n=3, **over):
    fields = {
        "x_um": np.zeros(n),
        "y_um": np.zeros(n),
        "depth_index": np.zeros(n, dtype=np.int64),
        "size_class": np.zeros(n, dtype=np.int64),
        "cluster_size": np.zeros(n, dtype=np.int64),
        "init_time_s": np.full(n, np.inf),
        "developed_fraction": np.zeros(n),
    }
    return GrainSet(**(fields | over))


@pytest.mark.parametrize(
    ("over", "match"),
    [
        ({"x_um": np.zeros(2)}, "entries"),
        ({"y_um": np.zeros((3, 1))}, "1-D"),
        ({"x_um": np.array([0.0, np.nan, 0.0])}, "finite"),
        ({"size_class": np.array([0, -1, 0])}, "negative"),
        ({"depth_index": np.array([0.0, 1.0, 0.0])}, "integer"),
        ({"developed_fraction": np.array([0.0, 1.5, 0.0])}, r"\[0, 1\]"),
        ({"developed_fraction": np.array([0.0, np.nan, 0.0])}, r"\[0, 1\]"),
        ({"init_time_s": np.array([0.0, -1.0, 0.0])}, "init_time_s"),
    ],
)
def test_grainset_rejects_bad_fields(over, match):
    with pytest.raises(ValueError, match=match):
        _grainset(**over)


def test_grainset_copies_and_freezes():
    x = np.zeros(3)
    g = _grainset(x_um=x)
    x[0] = 5.0
    assert g.x_um[0] == 0.0
    with pytest.raises(ValueError):
        g.x_um[0] = 1.0
    assert g.count == 3
