import itertools
import math

import numpy as np
import pytest
from scipy.integrate import quad

from darkroom_atlas.grain.overlap import boolean_excess_variance, disk_intersection_area, kappa

# OPT 5.4 / spec 12.7 table of kappa(D) for constant-radius disks.
_KAPPA_TABLE = {0.1: 1.055, 0.3: 1.182, 0.5: 1.336, 1.0: 1.891, 1.5: 2.861, 2.0: 4.655}


@pytest.mark.parametrize(("d", "expected"), sorted(_KAPPA_TABLE.items()))
def test_overlap_kappa_table(d, expected):
    assert abs(kappa(d) - expected) <= 0.005


def test_kappa_gauss_legendre_digits():
    # Spec 12.7: 48-point Gauss-Legendre gives 1.8913 (D = 1) and 4.6548 (D = 2).
    assert round(kappa(1.0), 3) == 1.891
    assert round(kappa(2.0), 3) == 4.655
    assert abs(kappa(1.0) - 1.8913) <= 1e-4
    assert abs(kappa(2.0) - 4.6548) <= 1.5e-4


def _var_count(lam, r, area):
    # Spec 2.7.1: Var_count / T^2 = exp(mu <q^2>) - 1 with mu = lam A, q = pi r^2 / A.
    return math.expm1(lam * area * (math.pi * r * r / area) ** 2)


def test_siedentopf_limit():
    r, area = 0.5, 1.0e4
    ratios = []
    for d in (1e-1, 1e-2, 1e-3, 1e-4):
        lam = d * math.log(10) / (math.pi * r * r)
        var_bool = boolean_excess_variance(np.array([lam]), np.array([r]), 1.0, area)
        ratios.append(float(var_bool) / _var_count(lam, r, area))
    # The ratio tends to 1 from above, roughly linearly in D.
    assert all(a > b > 1.0 for a, b in itertools.pairwise(ratios))
    assert abs(ratios[-1] - 1.0) <= 1e-3


@pytest.mark.xfail(
    strict=True,
    reason="spec 2.7.4 asks for <= 0.3 % at D = 0.01, but kappa(0.01) = 1.0053 for "
    "constant radius, so Var_bool / Var_count - 1 is about 0.53 % at large pixel area",
)
def test_siedentopf_limit_d001_tolerance():
    r, area = 0.5, 1.0e4
    lam = 0.01 * math.log(10) / (math.pi * r * r)
    var_bool = boolean_excess_variance(np.array([lam]), np.array([r]), 1.0, area)
    assert abs(float(var_bool) / _var_count(lam, r, area) - 1.0) <= 0.003


def test_disk_intersection_limits():
    r = np.array([0.3, 1.0, 2.5])
    np.testing.assert_allclose(disk_intersection_area(0.0, r), np.pi * r * r, rtol=1e-14)
    np.testing.assert_array_equal(disk_intersection_area(2.0 * r, r), 0.0)
    np.testing.assert_array_equal(disk_intersection_area(3.0 * r, r), 0.0)
    assert disk_intersection_area(1.0, 0.0) == 0.0
    # Integral over the plane equals (pi r^2)^2.
    val, _ = quad(lambda t: float(disk_intersection_area(t, 1.0)) * 2 * np.pi * t, 0, 2)
    assert abs(val - np.pi**2) <= 1e-9
    # Decreasing in tau.
    tau = np.linspace(0.0, 2.0, 101)
    assert np.all(np.diff(disk_intersection_area(tau, 1.0)) < 0.0)


def test_disk_intersection_half_overlap():
    # tau = r: 2 r^2 arccos(1/2) - (r/2) sqrt(3) r = r^2 (2 pi / 3 - sqrt(3) / 2).
    assert math.isclose(
        float(disk_intersection_area(1.0, 1.0)), 2 * math.pi / 3 - math.sqrt(3) / 2, rel_tol=1e-14
    )


def test_boolean_variance_matches_quadrature_for_mixed_sizes():
    lam = np.array([0.4, 0.1])
    r = np.array([0.3, 0.7])
    t_mean, area = 0.35, 900.0

    def integrand(t):
        g = sum(lam[c] * float(disk_intersection_area(t, r[c])) for c in range(2))
        return math.expm1(g) * 2 * math.pi * t

    ref = 0.0
    for a, b in ((0.0, 0.6), (0.6, 1.4)):
        ref += quad(integrand, a, b, epsabs=0, epsrel=1e-12)[0]
    ref *= t_mean**2 / area
    got = float(boolean_excess_variance(lam, r, t_mean, area))
    # 48 nodes across the kink at 2 r_1 = 0.6 um: a loose check.
    assert abs(got / ref - 1.0) <= 2e-3


def test_boolean_variance_shapes_and_zero():
    lam = np.zeros((4, 5, 3))
    r = np.array([0.2, 0.4, 0.6])
    out = boolean_excess_variance(lam, r, np.full((4, 5), 0.5), 100.0)
    assert out.shape == (4, 5)
    np.testing.assert_array_equal(out, 0.0)
    lam[1, 2] = [0.5, 0.2, 0.1]
    out = boolean_excess_variance(lam, r, 0.5, 100.0)
    assert out.shape == (4, 5)
    assert out[1, 2] > 0.0
    single = boolean_excess_variance(lam[1, 2], r, 0.5, 100.0)
    assert out[1, 2] == pytest.approx(float(single), rel=1e-14)
    # Scaling with T^2 / A.
    assert float(boolean_excess_variance(lam[1, 2], r, 0.25, 50.0)) == pytest.approx(
        float(single) / 2.0, rel=1e-14
    )


@pytest.mark.parametrize(
    ("lam", "r", "t", "a"),
    [
        ([-0.1], [1.0], 0.5, 1.0),
        ([0.1], [-1.0], 0.5, 1.0),
        ([0.1], [1.0], 1.5, 1.0),
        ([0.1], [1.0], 0.5, 0.0),
        (0.1, 1.0, 0.5, 1.0),
        ([0.1, 0.2], [1.0, 2.0, 3.0], 0.5, 1.0),
    ],
)
def test_boolean_variance_bad_input(lam, r, t, a):
    with pytest.raises(ValueError):
        boolean_excess_variance(np.asarray(lam), np.asarray(r), t, a)


@pytest.mark.parametrize("bad", [0.0, -1.0, np.nan, np.inf])
def test_kappa_bad_density(bad):
    with pytest.raises(ValueError):
        kappa(bad)


def test_bad_nodes_and_tau():
    with pytest.raises(ValueError):
        kappa(1.0, nodes=0)
    with pytest.raises(ValueError):
        disk_intersection_area(-1.0, 1.0)
