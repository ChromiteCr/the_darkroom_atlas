"""Latent-image chain: limits and design-time checks (spec D1 2.2.9, 12.1-12.3)."""

import math
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.optimize import brentq
from scipy.special import ndtr
from scipy.stats import poisson

from darkroom_atlas.data.library import Library
from darkroom_atlas.latent import ChainParams, cluster_distribution, generator

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "library"


def _params(**overrides) -> ChainParams:
    values = {
        "efficiency": 0.5,
        "efficiency_size_exponent": 0.0,
        "tau0_median_s": 0.3,
        "tau0_sigma_ln": 1.0,
        "ionic_time_s": 1e-5,
        "competing_sites": 3.0,
        "hole_attack": 0.0,
        "cluster_cap": 5,
    }
    values.update(overrides)
    return ChainParams(**values)


def _fold_states(p_states: np.ndarray, cap: int) -> np.ndarray:
    """Fold a 36-state vector to n = max(k1, k2) with n = 1 merged into n = 0."""
    size = np.maximum.outer(np.arange(cap + 1), np.arange(cap + 1)).ravel()
    p = np.bincount(size, weights=p_states, minlength=cap + 1)
    p[0] += p[1]
    p[1] = 0.0
    return p


# Toy population of spec D1 12.2: log-normal ECD classes, absorption cross-section
# proportional to d**1.5, criterion number-weighted P(n >= 3), eta = 0.5. The spec does not
# give the size distribution; d50 = 1 um, sigma_d = 0.4, K = 6 (+-3 sigma, equal ln d widths)
# reproduces its tables to about 0.01 stop.
_T_EDGES = np.linspace(-3.0, 3.0, 7)
_CDF = ndtr(_T_EDGES)
_TOY_FRACTION = np.diff(_CDF) / (_CDF[-1] - _CDF[0])
_TOY_ECD = np.exp(0.4 * 0.5 * (_T_EDGES[:-1] + _T_EDGES[1:]))
_REFERENCE_TIME_S = 1e-2  # stops are quoted relative to 1/100 s (spec D1 12.2)


def _toy_response(log_h: float, t: float, params: ChainParams) -> float:
    h = math.exp(log_h)
    total = 0.0
    for f, d in zip(_TOY_FRACTION, _TOY_ECD, strict=True):
        p = cluster_distribution(h * d**1.5, t, params, ecd_um=float(d))
        total += f * p[3:].sum()
    return total


def _toy_exposure(level: float, t: float, params: ChainParams) -> float:
    return brentq(
        lambda lh: _toy_response(lh, t, params) - level,
        math.log(1e-3),
        math.log(1e6),
        xtol=1e-12,
    )


def _toy_stops(level: float, t: float, params: ChainParams) -> float:
    reference = _toy_exposure(level, _REFERENCE_TIME_S, params)
    return (_toy_exposure(level, t, params) - reference) / math.log(2.0)


# --- generator -----------------------------------------------------------------------------


def test_generator_shape_and_columns():
    params = _params(hole_attack=0.2)
    q = generator(1e3, 0.05, params)
    assert q.shape == (36, 36)
    np.testing.assert_allclose(q.sum(axis=0), 0.0, atol=1e-9)
    off = q - np.diag(np.diag(q))
    assert np.all(off >= 0.0)


def test_generator_rates():
    params = _params(competing_sites=4.0, ionic_time_s=1e-4, hole_attack=0.1)
    nu, tau0 = 2e4, 0.02
    q = generator(nu, tau0, params)
    q1 = 1.0 / (1.0 + nu * 1e-4)
    q2 = (1.0 - q1) / 3.0
    n = 6

    def idx(k1, k2):
        return k1 * n + k2

    assert q[idx(1, 0), idx(0, 0)] == pytest.approx(nu * q1, rel=1e-14)
    assert q[idx(0, 1), idx(0, 0)] == pytest.approx(nu * q2, rel=1e-14)
    assert q[idx(0, 0), idx(1, 0)] == pytest.approx(1 / tau0 + 0.1 * nu, rel=1e-14)
    assert q[idx(2, 0), idx(1, 0)] == pytest.approx(nu * q1, rel=1e-14)
    # The capped site does not grow; lost electrons change nothing.
    assert q[idx(5, 0), idx(5, 0)] == pytest.approx(-nu * q2, rel=1e-14)
    assert q[idx(5, 5), idx(5, 5)] == 0.0


def test_generator_without_decay():
    q = generator(10.0, math.inf, _params(tau0_median_s=math.inf))
    assert q[0, 6] == 0.0  # (1, 0) -> (0, 0)
    assert np.all(np.isfinite(q))


@pytest.mark.parametrize(("nu", "tau0"), [(-1.0, 1.0), (math.nan, 1.0), (1.0, 0.0), (1.0, -1)])
def test_generator_rejects_bad_input(nu, tau0):
    with pytest.raises(ValueError, match="must be"):
        generator(nu, tau0, _params())


def test_matches_36_state_expm():
    """The factorised solution equals expm of the 36-state generator (spec D1 2.2.2)."""
    rng = np.random.default_rng(7)
    worst = 0.0
    for _ in range(100):
        params = _params(
            efficiency=rng.uniform(0.05, 1.0),
            tau0_median_s=10 ** rng.uniform(-1, 1),
            tau0_sigma_ln=0.0,
            ionic_time_s=10 ** rng.uniform(-7, -4),
            competing_sites=rng.uniform(2.0, 6.0),
            hole_attack=rng.uniform(0.0, 1.0),
        )
        t = 10 ** rng.uniform(-4, 1)
        m = 10 ** rng.uniform(-1, 1.5)
        nu = params.efficiency * m / t
        p_states = expm(generator(nu, params.tau0_median_s, params) * t)[:, 0]
        expected = _fold_states(p_states, 5)
        worst = max(worst, np.abs(cluster_distribution(m, t, params) - expected).max())
    assert worst <= 1e-12


# --- spec D1 2.2.9 -------------------------------------------------------------------------


def test_poisson_r_hit_limit():
    """Pure birth at rate eta nu: P(n >= r) is the Poisson tail of mean eta m (spec 12.1)."""
    params = _params(efficiency=0.7, tau0_median_s=math.inf, ionic_time_s=0.0, competing_sites=2.0)
    for m in (0.5, 3.0, 8.0, 20.0):
        for t in (1e-4, 1.0, 100.0):
            p = cluster_distribution(m, t, params)
            for r in range(2, 6):
                expected = poisson.sf(r - 1, 0.7 * m)
                assert abs(p[r:].sum() - expected) <= 1e-10, (m, t, r)


def test_probability_conservation():
    rng = np.random.default_rng(2026)
    for _ in range(200):
        params = _params(
            efficiency=rng.uniform(0.05, 1.0),
            efficiency_size_exponent=rng.uniform(-1.0, 1.0),
            tau0_median_s=10 ** rng.uniform(-2, 1),
            tau0_sigma_ln=rng.uniform(0.0, 4.0),
            ionic_time_s=10 ** rng.uniform(-7, -4),
            competing_sites=rng.uniform(2.0, 10.0),
            hole_attack=rng.uniform(0.0, 1.0),
            cluster_cap=int(rng.integers(2, 8)),
        )
        m = 10 ** rng.uniform(-6, 4)
        t = 10 ** rng.uniform(-5, 3)
        p = cluster_distribution(m, t, params, ecd_um=rng.uniform(0.2, 3.0))
        assert p.shape == (params.cluster_cap + 1,)
        assert abs(p.sum() - 1.0) <= 1e-12
        assert p[1] == 0.0
        assert np.all(p >= 0.0)


@pytest.mark.parametrize(("competing_sites", "hole_attack"), [(2.0, 0.0), (3.5, 0.0), (3.0, 0.4)])
def test_no_reciprocity_failure_without_decay_or_blocking(competing_sites, hole_attack):
    params = _params(
        tau0_median_s=math.inf,
        ionic_time_s=0.0,
        competing_sites=competing_sites,
        hole_attack=hole_attack,
    )
    for m in (0.1, 2.0, 15.0, 300.0):
        reference = cluster_distribution(m, 1.0, params)
        for t in np.logspace(-5, 3, 9):
            assert np.abs(cluster_distribution(m, t, params) - reference).max() <= 1e-10


def test_single_site_nucleation_formula():
    """Three-state chain 0 -> 1 -> 2 with 1 -> 0 decay; includes very stiff decay rates."""
    for tau0 in (1e-12, 1e-3, 0.3, 50.0):
        params = _params(
            efficiency=1.0,
            tau0_median_s=tau0,
            tau0_sigma_ln=0.0,
            ionic_time_s=0.0,
            competing_sites=2.0,
            cluster_cap=2,
        )
        for t in (1e-4, 1.0, 100.0):
            for m in (0.1, 3.0, 1e3):
                nu, d = m / t, 1.0 / tau0
                lam_minus = -(2 * nu + d + math.sqrt(d * d + 4 * nu * d)) / 2
                lam_plus = nu * nu / lam_minus
                survival = (
                    lam_plus * math.exp(lam_minus * t) - lam_minus * math.exp(lam_plus * t)
                ) / (lam_plus - lam_minus)
                p = cluster_distribution(m, t, params)
                assert p.shape == (3,)
                assert abs(p[2] - (1.0 - survival)) <= 1e-10, (tau0, t, m)


def test_lirf_raises_contrast():
    """Long exposures need more compensation at higher criteria (spec 12.2: 2.46 vs 1.89)."""
    params = _params(tau0_median_s=0.3, tau0_sigma_ln=2.0, ionic_time_s=0.0, competing_sites=2.0)
    assert _toy_stops(0.4, 100.0, params) > _toy_stops(0.05, 100.0, params)


def test_hirf_asymptote():
    """For nu tau_i -> 1e4 the rate saturates at 1/(S_eff - 1): 1 stop for S_eff = 3."""
    t, eta = 1.0, 0.5
    reciprocal = _params(
        efficiency=eta, tau0_median_s=math.inf, ionic_time_s=0.0, competing_sites=3.0
    )

    def solve(params):
        return math.exp(
            brentq(
                lambda lm: cluster_distribution(math.exp(lm), t, params)[2:].sum() - 0.5,
                math.log(1e-2),
                math.log(1e3),
                xtol=1e-13,
            )
        )

    m_reciprocal = solve(reciprocal)
    nu_at_solution = eta * 2.0 * m_reciprocal / t
    blocked = _params(
        efficiency=eta,
        tau0_median_s=math.inf,
        ionic_time_s=1e4 / nu_at_solution,
        competing_sites=3.0,
    )
    stops = math.log2(solve(blocked) / m_reciprocal)
    assert stops == pytest.approx(1.0, abs=0.02)


def test_hirf_lowers_contrast():
    """Short exposures need more compensation at higher criteria (spec 12.3: 0.37 vs 0.21)."""
    params = _params(tau0_median_s=0.074, tau0_sigma_ln=1.12, ionic_time_s=1e-5)
    assert _toy_stops(0.4, 1e-4, params) > _toy_stops(0.1, 1e-4, params)


@pytest.mark.xfail(
    strict=True,
    reason=(
        "spec D1 2.2.9: 9-node Gauss-Hermite is not within 1e-5 of 21 nodes for sigma_tau <= 4."
        " Measured max |dp| over m in [1e-6, 1e4], tau0 in {0.01, 0.3, 3} s, t in 1e-5..1e3 s:"
        " 1.4e-4 (sigma 1), 5.9e-3 (sigma 2), 2.0e-2 (sigma 3), 3.5e-2 (sigma 4); only"
        " sigma = 0.5 passes. 21 nodes is itself 9.4e-3 from converged at sigma 4. The decay"
        " transition is ~1/sigma wide in x."
    ),
)
def test_gauss_hermite_converged():
    grid = np.logspace(-6, 4, 201)
    worst = 0.0
    for sigma in (0.5, 1.0, 2.0, 3.0, 4.0):
        for tau0 in (0.01, 0.3, 3.0):
            for t in (1e-5, 1e-2, 1.0, 100.0, 1e3):
                params = _params(tau0_median_s=tau0, tau0_sigma_ln=sigma)
                for m in grid[::5]:
                    a = cluster_distribution(m, t, params, quadrature=9)
                    b = cluster_distribution(m, t, params, quadrature=21)
                    worst = max(worst, np.abs(a - b).max())
    assert worst <= 1e-5


def test_gauss_hermite_converged_small_spread():
    """The part of the convergence criterion that holds: sigma_tau = 0.5."""
    grid = np.logspace(-6, 4, 41)
    for tau0 in (0.01, 0.3, 3.0):
        for t in (1e-5, 1e-2, 1.0, 100.0, 1e3):
            params = _params(tau0_median_s=tau0, tau0_sigma_ln=0.5)
            for m in grid:
                a = cluster_distribution(m, t, params, quadrature=9)
                b = cluster_distribution(m, t, params, quadrature=21)
                assert np.abs(a - b).max() <= 1e-5


# --- design-time tables (spec D1 12.2, 12.3) -----------------------------------------------


@pytest.mark.parametrize(
    ("sigma", "level", "expected"),
    [
        (0.0, 0.05, (0.31, 1.22, 2.56, 4.10)),
        (2.0, 0.05, (0.35, 0.93, 1.89, 3.21)),
        (2.0, 0.4, (0.45, 1.22, 2.46, 3.98)),
        (3.0, 0.05, (0.35, 0.80, 1.47, 2.53)),
    ],
)
def test_lirf_design_table(sigma, level, expected):
    """Spec D1 12.2 with the toy population above (two sites, no loss, no blocking)."""
    params = _params(tau0_median_s=0.3, tau0_sigma_ln=sigma, ionic_time_s=0.0, competing_sites=2.0)
    stops = [_toy_stops(level, t, params) for t in (1.0, 10.0, 100.0, 1000.0)]
    np.testing.assert_allclose(stops, expected, atol=0.05)


def test_trap_depth_spread_flattens_lirf():
    """A single tau0 gives ~1.3-1.5 stops per decade over 10-1000 s; sigma_tau flattens it."""

    def slope(sigma):
        params = _params(
            tau0_median_s=0.3, tau0_sigma_ln=sigma, ionic_time_s=0.0, competing_sites=2.0
        )
        return (_toy_stops(0.05, 1000.0, params) - _toy_stops(0.05, 10.0, params)) / 2.0

    single, spread2, spread3 = slope(0.0), slope(2.0), slope(3.0)
    assert 1.3 <= single <= 1.5
    assert spread3 < spread2 < single


@pytest.mark.parametrize(
    ("competing_sites", "ionic_time_s", "expected"),
    [
        (2.0, 3e-6, ((0.43, 0.45), (0.04, 0.09), (-0.02, None))),
        (2.0, 1e-5, ((0.43, 0.32), (0.19, 0.29), (0.0, None))),
        (3.0, 1e-5, ((1.14, 1.13), (0.21, 0.37), (0.0, None))),
        (3.0, 3e-5, ((1.10, 1.05), (0.66, 0.90), (0.04, None))),
        (5.0, 1e-5, ((1.88, 1.99), (0.21, 0.40), (0.0, None))),
    ],
)
def test_hirf_design_table(competing_sites, ionic_time_s, expected):
    """Spec D1 12.3: criterion 0.1 (and 0.4) at 1/100 000, 1/10 000 and 1/1000 s."""
    params = _params(
        tau0_median_s=0.074,
        tau0_sigma_ln=1.12,
        ionic_time_s=ionic_time_s,
        competing_sites=competing_sites,
    )
    for t, (low, high) in zip((1e-5, 1e-4, 1e-3), expected, strict=True):
        assert _toy_stops(0.1, t, params) == pytest.approx(low, abs=0.05)
        if high is not None:
            assert _toy_stops(0.4, t, params) == pytest.approx(high, abs=0.05)


# --- parameters and input checks -----------------------------------------------------------


def test_from_layer_fixture_film():
    material = Library(FIXTURES).material("test/bw-mechanistic")
    (layer,) = material.emulsion_layers
    params = ChainParams.from_layer(layer)
    assert params == ChainParams(
        efficiency=0.3,
        efficiency_size_exponent=0.0,
        tau0_median_s=0.3,
        tau0_sigma_ln=1.0,
        ionic_time_s=1e-5,
        competing_sites=3.0,
        hole_attack=0.0,
        cluster_cap=5,
    )
    assert params.states == 36


def test_from_layer_requires_latent():
    material = Library(FIXTURES).material("test/bw-mechanistic")
    (layer,) = material.emulsion_layers
    with pytest.raises(ValueError, match="latent parameters are required"):
        ChainParams.from_layer(layer.model_copy(update={"latent": None}))


@pytest.mark.parametrize(
    ("overrides", "match"),
    [
        ({"efficiency": 0.0}, "efficiency"),
        ({"efficiency": math.nan}, "efficiency"),
        ({"efficiency_size_exponent": math.inf}, "efficiency_size_exponent"),
        ({"tau0_median_s": 0.0}, "tau0_median_s"),
        ({"tau0_median_s": math.nan}, "tau0_median_s"),
        ({"tau0_sigma_ln": -0.1}, "tau0_sigma_ln"),
        ({"ionic_time_s": -1e-6}, "ionic_time_s"),
        ({"competing_sites": 1.5}, "competing_sites"),
        ({"hole_attack": -0.1}, "hole_attack"),
        ({"cluster_cap": 1}, "cluster_cap"),
        ({"cluster_cap": 2.0}, "cluster_cap"),
        ({"cluster_cap": True}, "cluster_cap"),
    ],
)
def test_chain_params_validation(overrides, match):
    with pytest.raises(ValueError, match=match):
        _params(**overrides)


@pytest.mark.parametrize(
    ("args", "kwargs", "match"),
    [
        ((-1.0, 1.0), {}, "mean_photons"),
        ((math.nan, 1.0), {}, "mean_photons"),
        ((1.0, 0.0), {}, "exposure_time_s"),
        ((1.0, math.inf), {}, "exposure_time_s"),
        ((1.0, 1.0), {"ecd_um": 0.0}, "ecd_um"),
        ((1.0, 1.0), {"quadrature": 0}, "quadrature"),
        ((1.0, 1.0), {"quadrature": 9.0}, "quadrature"),
    ],
)
def test_cluster_distribution_rejects_bad_input(args, kwargs, match):
    with pytest.raises(ValueError, match=match):
        cluster_distribution(*args, _params(), **kwargs)


def test_zero_exposure_is_empty():
    p = cluster_distribution(0.0, 1.0, _params())
    np.testing.assert_array_equal(p, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0])


def test_efficiency_size_exponent():
    """eta_c = eta_0 (d_c / 1 um)**xi: doubling d with xi = 1 equals doubling m."""
    params = _params(efficiency_size_exponent=1.0, tau0_median_s=math.inf, ionic_time_s=0.0)
    np.testing.assert_allclose(
        cluster_distribution(2.0, 1.0, params, ecd_um=2.0),
        cluster_distribution(4.0, 1.0, params, ecd_um=1.0),
        rtol=0,
        atol=1e-14,
    )
