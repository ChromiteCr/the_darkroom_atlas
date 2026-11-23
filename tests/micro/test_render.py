"""Microscope-mode coverage rendering (spec D1 2.7.2, steps 4 and 5; 10 D1f)."""

import dataclasses
import math

import numpy as np
import pytest

from darkroom_atlas.data.models import EmulsionLayer
from darkroom_atlas.emulsion import size_classes
from darkroom_atlas.grain.overlap import boolean_excess_variance
from darkroom_atlas.micro import GrainSet, place_grains, render_coverage


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


@pytest.mark.parametrize("radius", [0.05, 0.3, 1.0, 2.7])
@pytest.mark.parametrize("ratio", [4.0, 4.5, 7.3])
def test_render_single_disk_area(radius, ratio):
    pitch = radius / ratio  # p_r <= r / 4
    rng = np.random.default_rng(int(radius * 1000 + ratio * 10))
    side = math.ceil(4.0 * radius / pitch) * pitch
    for _ in range(5):
        cx, cy = rng.uniform(1.5 * radius, side - 1.5 * radius, size=2)
        g = grains_at(cx, cy)
        cov = render_coverage(g, np.array([math.pi * radius**2]), (0.0, 0.0, side, side), pitch)
        area = cov.sum() * pitch * pitch
        assert abs(area / (math.pi * radius**2) - 1.0) <= 1e-4
        assert cov.min() >= 0.0 and cov.max() <= 1.0


def test_developed_fraction_scales_the_disk_area():
    pitch, side = 0.05, 4.0
    g = grains_at([2.0], [2.0], developed=[0.36])
    cov = render_coverage(g, np.array([3.0]), (0.0, 0.0, side, side), pitch)
    assert cov.sum() * pitch**2 == pytest.approx(0.36 * 3.0, rel=1e-10)


def test_overlapping_disks_cover_once():
    pitch, side, r = 0.05, 6.0, 1.0
    area = np.array([math.pi * r * r])
    same = render_coverage(grains_at([3.0, 3.0], [3.0, 3.0]), area, (0, 0, side, side), pitch)
    one = render_coverage(grains_at([3.0], [3.0]), area, (0, 0, side, side), pitch)
    # Pixels a disk covers fully stay at 1; only edge pixels combine as 1 - (1 - c)^2, so the
    # excess is at most a quarter pixel along the perimeter.
    np.testing.assert_array_equal(same[one == 1.0], 1.0)
    excess = (same - one).sum() * pitch**2
    assert 0.0 <= excess <= 0.25 * 2.0 * math.pi * r * pitch * math.sqrt(2.0)
    # Disjoint disks add.
    apart = render_coverage(grains_at([1.5, 4.5], [3.0, 3.0]), area, (0, 0, side, side), pitch)
    assert apart.sum() * pitch**2 == pytest.approx(2.0 * math.pi, rel=1e-12)
    # Partly overlapping disks: covered area lies between one disk and the exact union.
    d = 1.2
    lens = 2 * r * r * math.acos(d / (2 * r)) - 0.5 * d * math.sqrt(4 * r * r - d * d)
    union = 2.0 * math.pi * r * r - lens
    part = render_coverage(grains_at([2.4, 3.6], [3.0, 3.0]), area, (0, 0, side, side), pitch)
    assert part.sum() * pitch**2 == pytest.approx(union, rel=2e-3)
    assert part.max() <= 1.0


def test_disk_reaching_into_the_window_counts():
    pitch, r = 0.05, 1.0
    # Centred on the left edge of the window: half the disk is inside.
    cov = render_coverage(grains_at([0.0], [2.0]), np.array([math.pi]), (0, 0, 4.0, 4.0), pitch)
    assert cov.sum() * pitch**2 == pytest.approx(0.5 * math.pi * r * r, rel=1e-12)
    # Entirely outside: nothing.
    cov = render_coverage(grains_at([-1.5], [2.0]), np.array([math.pi]), (0, 0, 4.0, 4.0), pitch)
    assert not cov.any()


def test_undeveloped_and_zero_area_grains_are_clear():
    g = grains_at([1.0, 2.0], [1.0, 2.0], size_class=[0, 1], developed=[0.0, 1.0])
    cov = render_coverage(g, np.array([1.0, 0.0]), (0, 0, 3.0, 3.0), 0.1)
    assert cov.shape == (30, 30)
    assert not cov.any()


def test_raster_layout_rows_are_y():
    g = grains_at([0.5], [2.5])
    cov = render_coverage(g, np.array([0.04]), (0.0, 0.0, 1.0, 3.0), 0.1)
    assert cov.shape == (30, 10)
    j, i = np.unravel_index(np.argmax(cov), cov.shape)
    assert (i, j) in {(4, 24), (5, 24), (4, 25), (5, 25)}


def test_render_is_deterministic():
    classes = size_classes(_layer(0.25))
    g = place_grains(classes, (60.0, 60.0), 4, seed=9)
    g = dataclasses.replace(g, developed_fraction=np.linspace(0.0, 1.0, g.count))
    areas = 2.0 * np.asarray(classes.projected_um2)
    a = render_coverage(g, areas, (5.0, 5.0, 55.0, 55.0), 0.1)
    b = render_coverage(g, areas, (5.0, 5.0, 55.0, 55.0), 0.1)
    np.testing.assert_array_equal(a, b)


def _q(value, unit="1"):
    return {"value": value, "unit": unit, "source": "assumed"}


def _layer(sigma_log: float) -> EmulsionLayer:
    return EmulsionLayer.model_validate(
        {
            "id": "emulsion",
            "gelatin_coverage": _q(5.0, "g/m2"),
            "grain": {
                "morphology": "conventional",
                "mean_diameter": _q(0.9, "um"),
                "sigma_log": _q(sigma_log),
                "silver_coverage": _q(4.0, "g/m2"),
                "n_classes": 6,
            },
        }
    )


# (net density, patch side in um, number of patches, raster pitch in um)
_NUTTING_CASES = [(0.3, 600.0, 1, 0.25), (1.0, 1200.0, 1, 0.5), (2.0, 1900.0, 2, 1.0)]


@pytest.mark.parametrize(("density", "side", "patches", "pitch"), _NUTTING_CASES)
def test_micro_mean_density_matches_nutting(density, side, patches, pitch):
    # Six classes (spec 2.0.2), fully developed, disk areas proportional to the projected
    # areas and scaled so that the Nutting density 0.4343 sum_c nu_c a_c is `density` (2.5.1).
    classes = size_classes(_layer(0.25))
    nu = np.asarray(classes.number_per_um2)
    areas = density * math.log(10.0) / np.sum(nu * classes.projected_um2) * classes.projected_um2
    radii = np.sqrt(areas / math.pi)
    margin = (math.ceil(radii.max() / pitch) + 1) * pitch
    inner = side - 2.0 * margin

    clear = []
    for k in range(patches):
        x0 = k * side
        g = place_grains(classes, (side, side), 8, seed=21, origin_um=(x0, 0.0))
        g = dataclasses.replace(g, developed_fraction=np.ones(g.count))
        window = (x0 + margin, margin, x0 + margin + inner, margin + inner)
        clear.append(1.0 - render_coverage(g, areas, window, pitch).mean())
    measured = -math.log10(float(np.mean(clear)))

    # The patch is big enough that 0.005 is at least four standard errors (Boolean model).
    t_mean = 10.0**-density
    var = boolean_excess_variance(nu, radii, t_mean, patches * inner * inner)
    sigma_d = math.sqrt(float(var)) / (t_mean * math.log(10.0))
    assert 4.0 * sigma_d <= 0.005
    assert abs(measured - density) <= 0.005


@pytest.mark.parametrize(
    ("window", "pitch", "match"),
    [
        ((0.0, 0.0, 1.05, 1.0), 0.1, "multiple"),
        ((0.0, 0.0, 0.0, 1.0), 0.1, "x1 > x0"),
        ((0.0, 1.0, 1.0, 0.0), 0.1, "x1 > x0"),
        ((0.0, 0.0, float("nan"), 1.0), 0.1, "finite"),
        ((0.0, 0.0, 1.0), 0.1, "window_um"),
        ((0.0, 0.0, 1.0, 1.0), 0.0, "pitch_um"),
        ((0.0, 0.0, 1.0, 1.0), -0.1, "pitch_um"),
        ((0.0, 0.0, 1.0, 1.0), float("inf"), "pitch_um"),
    ],
)
def test_render_rejects_bad_window(window, pitch, match):
    with pytest.raises(ValueError, match=match):
        render_coverage(grains_at([0.5], [0.5]), np.array([0.1]), window, pitch)


@pytest.mark.parametrize(
    "areas", [np.array([0.1]), np.array([-0.1, 0.1]), np.array([np.nan, 0.1]), np.ones((2, 1))]
)
def test_render_rejects_bad_areas(areas):
    g = grains_at([0.5, 0.6], [0.5, 0.6], size_class=[0, 1])
    with pytest.raises(ValueError, match="areas_um2"):
        render_coverage(g, areas, (0.0, 0.0, 1.0, 1.0), 0.1)
