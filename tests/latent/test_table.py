"""Latent-image lookup table (spec D1 2.2.5, 2.2.9, D1b)."""

import math
import time
from pathlib import Path

import numpy as np
import pytest

from darkroom_atlas.data.library import Library
from darkroom_atlas.emulsion import size_classes
from darkroom_atlas.latent import ChainParams, LatentTable, cluster_distribution, latent_table

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "library"


@pytest.fixture(scope="module")
def layer():
    material = Library(FIXTURES).material("test/bw-mechanistic")
    (emulsion,) = material.emulsion_layers
    return emulsion


@pytest.fixture(scope="module")
def classes(layer):
    return size_classes(layer)


def _params(**overrides) -> ChainParams:
    values = {
        "efficiency": 0.5,
        "efficiency_size_exponent": 0.3,
        "tau0_median_s": 0.074,
        "tau0_sigma_ln": 1.12,
        "ionic_time_s": 7.8e-6,
        "competing_sites": 3.0,
    }
    values.update(overrides)
    return ChainParams(**values)


_CASES = [
    (_params(), 1e-2),
    (_params(), 1e-5),
    (_params(), 100.0),
    (_params(tau0_median_s=1.92, tau0_sigma_ln=3.47, ionic_time_s=1.8e-5), 1.0),
    (_params(efficiency=1.0, tau0_median_s=math.inf, ionic_time_s=0.0, competing_sites=2.0), 1.0),
]


@pytest.mark.parametrize(("params", "t"), _CASES)
def test_monotone_in_exposure(classes, params, t):
    table = latent_table(classes, params, t)
    tails = np.cumsum(table.probabilities[..., ::-1], axis=-1)[..., ::-1]  # P(n >= r)
    steps = np.diff(tails[..., 2:], axis=1)
    # Non-decreasing up to floating-point jitter on the saturated plateaus.
    assert steps.min() >= -1e-12


@pytest.mark.parametrize(("params", "t"), _CASES)
def test_latent_table_interpolation(classes, params, t):
    table = latent_table(classes, params, t)
    mid = 0.5 * (table.log10_m[1:] + table.log10_m[:-1])
    for c in (0, classes.count // 2, classes.count - 1):
        looked_up = table.lookup(10.0**mid, c)
        direct = np.array(
            [cluster_distribution(10.0**u, t, params, ecd_um=float(classes.ecd_um[c])) for u in mid]
        )
        assert np.abs(looked_up - direct).max() <= 1e-4


def test_latent_table_deterministic(classes, layer):
    params = ChainParams.from_layer(layer)
    a = latent_table(classes, params, 1e-2)
    b = latent_table(classes, params, 1e-2)
    assert a.probabilities.tobytes() == b.probabilities.tobytes()
    assert a.log10_m.tobytes() == b.log10_m.tobytes()
    m = np.logspace(-7, 5, 333)
    assert a.lookup(m, 3).tobytes() == b.lookup(m, 3).tobytes()


def test_latent_table_build_time(classes):
    """K = 6, 201 nodes, 9 quadrature nodes in under 1 s on one core (D1b); best of three."""
    params = _params(tau0_sigma_ln=3.0)
    assert classes.count == 6
    best = math.inf
    for _ in range(3):
        start = time.perf_counter()
        latent_table(classes, params, 100.0, nodes=201, quadrature=9)
        best = min(best, time.perf_counter() - start)
    assert best < 1.0


def test_table_layout_and_nodes(classes):
    params = _params()
    table = latent_table(classes, params, 1e-2)
    assert isinstance(table, LatentTable)
    assert table.exposure_time_s == 1e-2
    assert table.log10_m.shape == (201,)
    np.testing.assert_allclose(table.log10_m[[0, 100, -1]], [-6.0, -1.0, 4.0], atol=1e-12)
    assert table.probabilities.shape == (6, 201, 6)
    assert not table.probabilities.flags.writeable
    assert not table.log10_m.flags.writeable
    assert np.all(table.probabilities[..., 1] == 0.0)
    np.testing.assert_allclose(table.probabilities.sum(axis=-1), 1.0, rtol=0, atol=1e-12)
    # Below the range p(n >= 2) is negligible (spec D1 2.2.5).
    assert table.probabilities[:, 0, 2:].sum(axis=-1).max() < 1e-12
    # Nodes are the direct chain values.
    for c in (0, 5):
        for i in (0, 77, 200):
            direct = cluster_distribution(
                10.0 ** table.log10_m[i], 1e-2, params, ecd_um=float(classes.ecd_um[c])
            )
            np.testing.assert_allclose(table.probabilities[c, i], direct, rtol=0, atol=1e-14)


def test_lookup_nodes_clamping_and_shape(classes):
    table = latent_table(classes, _params(), 1.0)
    nodes = 10.0**table.log10_m
    np.testing.assert_allclose(table.lookup(nodes, 2), table.probabilities[2], atol=1e-13)
    np.testing.assert_allclose(table.lookup(0.0, 2), table.probabilities[2, 0], atol=1e-15)
    np.testing.assert_allclose(table.lookup(1e-9, 2), table.probabilities[2, 0], atol=1e-15)
    np.testing.assert_allclose(table.lookup(1e7, 2), table.probabilities[2, -1], atol=1e-13)
    np.testing.assert_allclose(table.lookup(math.inf, 2), table.probabilities[2, -1], atol=1e-13)
    m = np.logspace(-8, 6, 24).reshape(2, 3, 4)
    p = table.lookup(m, 4)
    assert p.shape == (2, 3, 4, 6)
    assert np.all(p >= 0.0)
    assert np.all(p[..., 1] == 0.0)
    np.testing.assert_allclose(p.sum(axis=-1), 1.0, rtol=0, atol=1e-12)


@pytest.mark.parametrize(
    ("m", "size_class", "match"),
    [
        (-1.0, 0, "mean_photons"),
        (np.array([1.0, np.nan]), 0, "mean_photons"),
        (1.0, 6, "size_class"),
        (1.0, -1, "size_class"),
        (1.0, 1.0, "size_class"),
    ],
)
def test_lookup_rejects_bad_input(classes, m, size_class, match):
    table = latent_table(classes, _params(), 1.0, nodes=11)
    with pytest.raises(ValueError, match=match):
        table.lookup(m, size_class)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"exposure_time_s": 0.0}, "exposure_time_s"),
        ({"exposure_time_s": 1.0, "nodes": 1}, "nodes"),
        ({"exposure_time_s": 1.0, "log10_range": (4.0, -6.0)}, "log10_range"),
        ({"exposure_time_s": 1.0, "log10_range": (-6.0, math.inf)}, "log10_range"),
        ({"exposure_time_s": 1.0, "quadrature": 0}, "quadrature"),
    ],
)
def test_latent_table_rejects_bad_input(classes, kwargs, match):
    with pytest.raises(ValueError, match=match):
        latent_table(classes, _params(), **kwargs)


def test_latent_table_rate_overflow_raises(classes):
    """10**400 photons overflow to inf: the build raises instead of hanging."""
    with pytest.raises(FloatingPointError, match="overflow"):
        latent_table(classes, _params(), 1.0, nodes=5, log10_range=(-6.0, 400.0))


def test_latent_table_rejects_bad_layout():
    with pytest.raises(ValueError, match="evenly spaced"):
        LatentTable(np.array([0.0, 1.0, 3.0]), np.full((1, 3, 3), 1.0 / 3.0), 1.0)
    with pytest.raises(ValueError, match="shape"):
        LatentTable(np.array([0.0, 1.0, 2.0]), np.full((1, 4, 3), 1.0 / 3.0), 1.0)
