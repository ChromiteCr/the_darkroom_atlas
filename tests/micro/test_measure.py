"""Microscope-mode aperture measurement (spec D1 2.7.2, step 5; 10 D1f)."""

import dataclasses
import math

import numpy as np
import pytest
from scipy import integrate

from darkroom_atlas.emulsion.classes import SizeClasses
from darkroom_atlas.grain.overlap import disk_intersection_area, kappa
from darkroom_atlas.micro import GrainSet, aperture_transmittances, place_grains
from darkroom_atlas.micro.measure import default_pitch


def _ro(values) -> np.ndarray:
    array = np.array(values, dtype=np.float64)
    array.flags.writeable = False
    return array


def make_classes(number_per_um2) -> SizeClasses:
    """Size classes with given surface densities; only ``number_per_um2`` matters here."""
    nu = np.atleast_1d(np.asarray(number_per_um2, dtype=np.float64))
    d = np.full(nu.shape, 0.5)
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


def grains_at(x, y, size_class=None, developed=None) -> GrainSet:
    x = np.atleast_1d(np.asarray(x, dtype=np.float64))
    y = np.atleast_1d(np.asarray(y, dtype=np.float64))
    n = x.size
    return GrainSet(
        x_um=x,
        y_um=y,
        depth_index=np.zeros(n, dtype=np.int64),
        size_class=np.zeros(n, dtype=np.int64) if size_class is None else np.asarray(size_class),
        cluster_size=np.zeros(n, dtype=np.int64),
        init_time_s=np.zeros(n),
        developed_fraction=np.ones(n) if developed is None else np.asarray(developed, float),
    )


def lens_area(d: float, r1: float, r2: float) -> float:
    """Intersection area of two disks of radii r1, r2 whose centres are d apart."""
    if d >= r1 + r2:
        return 0.0
    if d <= abs(r1 - r2):
        return math.pi * min(r1, r2) ** 2
    a1 = r1 * r1 * math.acos((d * d + r1 * r1 - r2 * r2) / (2 * d * r1))
    a2 = r2 * r2 * math.acos((d * d + r2 * r2 - r1 * r1) / (2 * d * r2))
    k = math.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
    return a1 + a2 - 0.5 * k


def transmit_one(x, y, r, aperture_radius=5.0, centre=(0.0, 0.0)) -> float:
    g = grains_at(x, y)
    t = aperture_transmittances(
        g, np.array([math.pi * r * r]), np.array([centre]), aperture_um=2.0 * aperture_radius
    )
    return float(t[0])


@pytest.mark.parametrize(
    ("d", "r"),
    [(0.0, 1.0), (2.0, 1.0), (4.0, 1.0), (4.5, 1.0), (5.0, 1.0), (5.8, 1.0), (5.3, 2.5)],
)
def test_single_disk_is_the_lens_area(d, r):
    for angle in (0.0, 0.7, 2.0, math.pi, 4.4):
        x, y = d * math.cos(angle), d * math.sin(angle)
        t = transmit_one([3.0 + x], [-2.0 + y], r, centre=(3.0, -2.0))
        assert 1.0 - t == pytest.approx(lens_area(d, 5.0, r) / (25.0 * math.pi), abs=1e-12)


@pytest.mark.parametrize(("d", "r"), [(0.0, 5.0), (0.0, 6.0), (1.0, 7.0), (1.0, 6.0)])
def test_disk_containing_the_aperture_is_opaque(d, r):
    assert transmit_one([d], [0.0], r) == pytest.approx(0.0, abs=1e-12)


def test_disks_outside_or_tangent_outside_are_clear():
    assert transmit_one([6.0], [0.0], 1.0) == 1.0
    assert transmit_one([0.0], [-9.0], 1.0) == 1.0


def test_no_developed_grain_is_clear():
    g = grains_at([0.0, 1.0], [0.0, 1.0], size_class=[0, 1], developed=[0.0, 1.0])
    t = aperture_transmittances(g, np.array([1.0, 0.0]), np.array([[0.0, 0.0], [1.0, 1.0]]))
    np.testing.assert_array_equal(t, 1.0)
    empty = grains_at([], [])
    np.testing.assert_array_equal(aperture_transmittances(empty, np.array([1.0]), [[0, 0]]), 1.0)
    assert aperture_transmittances(g, np.array([1.0, 1.0]), np.zeros((0, 2))).shape == (0,)


def test_coincident_and_nested_disks_count_once():
    a = math.pi
    # Two identical disks on top of each other, inside and straddling the edge.
    for cx in (0.0, 4.5):
        two = aperture_transmittances(
            grains_at([cx, cx], [0.0, 0.0]), np.array([a]), [[0.0, 0.0]], aperture_um=10.0
        )
        one = aperture_transmittances(
            grains_at([cx], [0.0]), np.array([a]), [[0.0, 0.0]], aperture_um=10.0
        )
        assert two[0] == pytest.approx(one[0], abs=1e-12)
    # A small disk internally tangent to a large one adds nothing.
    g = grains_at([0.0, 1.5], [0.0, 0.0], size_class=[0, 1])
    t = aperture_transmittances(g, np.array([4.0 * a, 0.25 * a]), [[0.0, 0.0]], aperture_um=10.0)
    assert 1.0 - t[0] == pytest.approx(4.0 / 25.0, abs=1e-12)
    # A disk internally tangent to the aperture edge lies wholly inside it.
    assert 1.0 - transmit_one([4.0], [0.0], 1.0) == pytest.approx(1.0 / 25.0, abs=1e-12)


def test_two_overlapping_disks_are_their_union():
    r, d = 1.0, 1.2
    t = aperture_transmittances(
        grains_at([-0.6, 0.6], [0.3, 0.3]), np.array([math.pi]), [[0.0, 0.0]], aperture_um=10.0
    )
    union = 2.0 * math.pi - lens_area(d, r, r)
    assert 1.0 - t[0] == pytest.approx(union / (25.0 * math.pi), abs=1e-12)


def test_exact_agrees_with_fine_raster():
    rng = np.random.default_rng(3)
    n = 3000
    g = grains_at(
        rng.uniform(-2.0, 32.0, n),
        rng.uniform(-2.0, 32.0, n),
        size_class=rng.integers(0, 3, n),
        developed=rng.uniform(0.0, 1.0, n),
    )
    areas = np.array([0.3, 0.8, 1.5])
    centres = rng.uniform(8.0, 22.0, (12, 2))
    exact = aperture_transmittances(g, areas, centres, aperture_um=12.0)
    assert 0.2 < exact.mean() < 0.6
    coarse = aperture_transmittances(g, areas, centres, aperture_um=12.0, pitch_um=0.05)
    fine = aperture_transmittances(g, areas, centres, aperture_um=12.0, pitch_um=0.025)
    assert np.abs(fine - exact).max() <= 1e-4
    # The raster converges to the exact value.
    assert np.abs(fine - exact).max() < np.abs(coarse - exact).max()


def test_order_and_batching_do_not_matter():
    classes = make_classes([1.5])
    g = place_grains(classes, (300.0, 300.0), 2, seed=4)
    g = dataclasses.replace(g, developed_fraction=np.full(g.count, 0.8))
    areas = np.array([0.5])
    c = np.arange(20.0, 281.0, 10.0)
    cx, cy = np.meshgrid(c, c)
    centres = np.stack([cx.ravel(), cy.ravel()], axis=1)  # 729 apertures, several batches
    t = aperture_transmittances(g, areas, centres, aperture_um=30.0)
    perm = np.random.default_rng(0).permutation(centres.shape[0])
    t_perm = aperture_transmittances(g, areas, centres[perm], aperture_um=30.0)
    np.testing.assert_allclose(t_perm, t[perm], rtol=0.0, atol=1e-12)
    single = [aperture_transmittances(g, areas, centres[i : i + 1], aperture_um=30.0)[0]
              for i in (0, 300, 728)]  # fmt: skip
    np.testing.assert_allclose(single, t[[0, 300, 728]], rtol=0.0, atol=1e-12)


def test_default_pitch_is_spec_p_r():
    g = grains_at([0.0, 1.0, 2.0], [0.0, 0.0, 0.0], size_class=[0, 0, 1], developed=[1.0, 0.0, 1.0])
    # Class 0: only the developed grain counts, r = sqrt(4 pi * 1 / pi) = 2 -> 0.5 > 0.05.
    assert default_pitch(g, np.array([4.0 * math.pi, 0.01 * math.pi])) == pytest.approx(0.025)
    assert default_pitch(g, np.array([4.0 * math.pi, 4.0 * math.pi])) == pytest.approx(0.05)
    assert default_pitch(grains_at([0.0], [0.0], developed=[0.0]), np.array([1.0])) is None


def _aperture_kappa(net_density: float, r: float, aperture_radius: float) -> float:
    """Expected aperture variance over Siedentopf for constant-radius disks (exact integral)."""
    a = math.pi * r * r
    lam = net_density * math.log(10.0) / a
    area = math.pi * aperture_radius**2

    def integrand(tau: float) -> float:
        gamma = float(disk_intersection_area(np.array(tau), np.array(r)))
        overlap = float(disk_intersection_area(np.array(tau), np.array(aperture_radius)))
        return math.expm1(lam * gamma) * overlap / area * 2.0 * math.pi * tau

    value, _ = integrate.quad(integrand, 0.0, 2.0 * r, limit=200)
    return value / (lam * a * a)


_KAPPA_CASES = [
    pytest.param(0.3, id="D0.3"),
    pytest.param(1.0, id="D1.0", marks=pytest.mark.slow),
    pytest.param(2.0, id="D2.0", marks=pytest.mark.slow),
]


@pytest.mark.parametrize("density", _KAPPA_CASES)
def test_micro_constant_radius_kappa(density):
    # One class of equal, fully developed disks of radius r with lambda a = D ln 10 (Nutting,
    # spec 2.5.1); apertures of diameter 60 r on a grid of spacing 30 r over six 1 mm patches.
    r = 0.5
    a = math.pi * r * r
    nu = density * math.log(10.0) / a
    classes = make_classes([nu])
    radius = 30.0 * r
    margin = radius + 2.0 * r
    side = 1000.0
    c = np.arange(margin, side - margin + 1e-9, radius)
    cx, cy = np.meshgrid(c, c)
    local = np.stack([cx.ravel(), cy.ravel()], axis=1)
    values = []
    for k in range(6):
        origin = (1000.0 * (k % 3), 1000.0 * (k // 3))
        g = place_grains(classes, (side, side), 1, seed=17, origin_um=origin)
        g = dataclasses.replace(g, developed_fraction=np.ones(g.count))
        values.append(
            aperture_transmittances(g, np.array([a]), local + origin, aperture_um=2.0 * radius)
        )
    t = np.concatenate(values)
    t_mean = 10.0**-density
    siedentopf = t_mean**2 * nu * a * a / (math.pi * radius * radius)
    ratio = float(t.var(ddof=1)) / siedentopf
    assert abs(ratio / kappa(density) - 1.0) <= 0.05

    # Tighter: the exact aperture value (b_ap below kappa, spec 2.7.1) within four standard
    # errors; overlapping apertures at spacing R carry 0.5735 independent samples per R^2.
    expected = _aperture_kappa(density, r, radius)
    m_eff = 0.5735 * 6 * (c.size * radius) ** 2 / radius**2
    assert abs(ratio / expected - 1.0) <= 4.0 * math.sqrt(2.0 / m_eff)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"centres_um": np.zeros(2)}, "centres_um"),
        ({"centres_um": np.zeros((2, 3))}, "centres_um"),
        ({"centres_um": np.array([[0.0, np.nan]])}, "centres_um"),
        ({"aperture_um": 0.0}, "aperture_um"),
        ({"aperture_um": -48.0}, "aperture_um"),
        ({"aperture_um": float("inf")}, "aperture_um"),
        ({"aperture_um": True}, "aperture_um"),
        ({"pitch_um": 0.0}, "pitch_um"),
        ({"pitch_um": float("nan")}, "pitch_um"),
        ({"pitch_um": True}, "pitch_um"),
        ({"areas_um2": np.array([-1.0])}, "areas_um2"),
        ({"areas_um2": np.array([])}, "areas_um2"),
    ],
)
def test_aperture_rejects_bad_input(kwargs, match):
    args = {"areas_um2": np.array([1.0]), "centres_um": np.zeros((1, 2))} | kwargs
    with pytest.raises(ValueError, match=match):
        aperture_transmittances(grains_at([0.0], [0.0]), **args)
