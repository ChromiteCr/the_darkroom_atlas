import numpy as np
import pytest
from scipy.stats import poisson

from darkroom_atlas.grain.poisson import poisson_inverse
from darkroom_atlas.rng import pixel_uniform


@pytest.mark.parametrize("mu", [0.1, 3.0, 30.0, 300.0, 3000.0])
def test_poisson_sampler_moments(mu):
    n = 1_000_000
    u = pixel_uniform(7, f"test/poisson/{mu}", (0, 1000), (0, 1000), 1000)
    k = poisson_inverse(u, np.full(u.shape, mu)).astype(np.float64)
    mean = k.mean()
    var = k.var(ddof=1)
    # Standard errors of the sample mean and variance for a Poisson variable.
    se_mean = np.sqrt(mu / n)
    se_var = np.sqrt((mu + 2.0 * mu * mu) / n)
    assert abs(mean - mu) <= 4.0 * se_mean
    assert abs(var - mu) <= 4.0 * se_var


def test_poisson_sampler_deterministic():
    seed, stream, width = 3, "grain/emulsion/0", 64
    mean = np.linspace(0.0, 400.0, 48 * width).reshape(48, width)
    whole = poisson_inverse(pixel_uniform(seed, stream, (0, 48), (0, width), width), mean)
    again = poisson_inverse(pixel_uniform(seed, stream, (0, 48), (0, width), width), mean)
    np.testing.assert_array_equal(whole, again)
    tiled = np.empty_like(whole)
    for r0 in range(0, 48, 17):
        for c0 in range(0, width, 23):
            r1, c1 = min(r0 + 17, 48), min(c0 + 23, width)
            u = pixel_uniform(seed, stream, (r0, r1), (c0, c1), width)
            tiled[r0:r1, c0:c1] = poisson_inverse(u, mean[r0:r1, c0:c1])
    np.testing.assert_array_equal(tiled, whole)


@pytest.mark.parametrize("mu", [0.0, 0.1, 1.0, 3.0, 12.5, 30.0])
def test_exact_branch_matches_scipy_ppf(mu):
    rng = np.random.default_rng(11)
    u = rng.uniform(1e-12, 1.0 - 1e-12, 20_000)
    k = poisson_inverse(u, np.full(u.shape, mu))
    expected = np.minimum(poisson.ppf(u, mu), np.floor(mu + 10 * np.sqrt(mu) + 10))
    np.testing.assert_array_equal(k, expected.astype(np.int64))


def test_cornish_fisher_close_to_ppf_above_threshold():
    u = np.linspace(0.001, 0.999, 999)
    for mu in (31.0, 300.0, 3000.0):
        k = poisson_inverse(u, np.full(u.shape, mu))
        assert np.max(np.abs(k - poisson.ppf(u, mu))) <= 1


def test_threshold_selects_branch():
    u = np.linspace(0.01, 0.99, 99)
    mu = np.full(u.shape, 20.0)
    exact = poisson_inverse(u, mu)
    approx = poisson_inverse(u, mu, exact_below=10.0)
    np.testing.assert_array_equal(exact, poisson.ppf(u, 20.0).astype(np.int64))
    assert np.max(np.abs(approx - exact)) <= 1


def test_broadcast_dtype_and_extremes():
    k = poisson_inverse(np.array([[0.0], [1.0]]), np.array([0.0, 2.0, 100.0]))
    assert k.shape == (2, 3)
    assert k.dtype == np.int64
    # u = 0: the exact branch gives 0; Cornish-Fisher clips u to 2^-53 (z about -8.2).
    np.testing.assert_array_equal(k[0, :2], [0, 0])
    assert 0 <= k[0, 2] < 40
    # u = 1: finite counts, never above floor(mu + 10 sqrt(mu) + 10).
    cap = np.floor(np.array([0.0, 2.0, 100.0]) + 10 * np.sqrt([0.0, 2.0, 100.0]) + 10)
    assert k[1, 0] == 0
    assert np.all(k[1, 1:] > k[0, 1:])
    assert np.all(k[1] <= cap)


@pytest.mark.parametrize(
    ("u", "mean"),
    [(-0.1, 1.0), (1.5, 1.0), (np.nan, 1.0), (0.5, -1.0), (0.5, np.inf), ([0.1, 0.2], [1, 2, 3])],
)
def test_bad_input_raises(u, mean):
    with pytest.raises(ValueError):
        poisson_inverse(np.asarray(u), np.asarray(mean))


def test_bad_threshold_raises():
    with pytest.raises(ValueError):
        poisson_inverse(np.array([0.5]), np.array([1.0]), exact_below=-1.0)
