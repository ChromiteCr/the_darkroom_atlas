"""Grain size classes (spec D1 2.0.2)."""

from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import norm

from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import EmulsionLayer
from darkroom_atlas.emulsion import SizeClasses, size_classes
from darkroom_atlas.emulsion.classes import M_AG, M_AGBR, M_AGI, RHO_AGX, RHO_GELATIN

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "library"


def _q(value, unit="1"):
    return {"value": value, "unit": unit, "source": "assumed"}


def _layer(
    *,
    morphology="conventional",
    d50=0.9,
    sigma=0.45,
    silver=4.0,
    gelatin=5.0,
    iodide=0.04,
    thickness=None,
    n_classes=6,
) -> EmulsionLayer:
    grain = {
        "morphology": morphology,
        "mean_diameter": _q(d50, "um"),
        "sigma_log": _q(sigma),
        "silver_coverage": _q(silver, "g/m2"),
        "n_classes": n_classes,
    }
    if iodide is not None:
        grain["iodide_fraction"] = _q(iodide)
    if thickness is not None:
        grain["thickness"] = _q(thickness, "um")
    data = {"id": "emulsion", "grain": grain}
    if gelatin is not None:
        data["gelatin_coverage"] = _q(gelatin, "g/m2")
    return EmulsionLayer.model_validate(data)


def test_size_classes_conserve_silver():
    rng = np.random.default_rng(20261121)
    for _ in range(200):
        tabular = bool(rng.integers(2))
        layer = _layer(
            morphology="tabular" if tabular else "conventional",
            d50=float(rng.uniform(0.2, 3.0)),
            sigma=float(rng.uniform(0.05, 0.9)),
            silver=float(rng.uniform(0.5, 8.0)),
            gelatin=float(rng.uniform(1.0, 10.0)),
            iodide=float(rng.uniform(0.0, 0.2)),
            thickness=float(rng.uniform(0.05, 0.3)) if tabular else None,
            n_classes=int(rng.integers(2, 65)),
        )
        classes = size_classes(layer)
        coated = 1e12 * np.sum(classes.number_per_um2 * classes.silver_g_per_grain)
        assert coated == pytest.approx(classes.silver_coverage_g_m2, rel=1e-12, abs=0.0)


def test_layer_thickness_formula():
    # LAT 2.3: 1.75 g/m2 Ag + 4.84 g/m2 gelatin gives 4.1 um (pure AgBr).
    classes = size_classes(_layer(silver=1.75, gelatin=4.84, iodide=0.0))
    assert classes.dry_thickness_um == pytest.approx(4.1, abs=0.05)
    expected = 1.75 * (M_AGBR / M_AG) / RHO_AGX + 4.84 / RHO_GELATIN
    assert classes.dry_thickness_um == pytest.approx(expected, rel=1e-15)


def test_dry_thickness_with_iodide():
    x = 0.04
    classes = size_classes(_layer(silver=4.0, gelatin=5.0, iodide=x))
    m_agx = (1 - x) * M_AGBR + x * M_AGI
    assert classes.dry_thickness_um == pytest.approx(
        4.0 * (m_agx / M_AG) / RHO_AGX + 5.0 / RHO_GELATIN, rel=1e-15
    )


@pytest.mark.parametrize("k", [2, 6, 12, 64])
def test_number_fractions_sum_to_one_and_match_truncated_lognormal(k):
    classes = size_classes(_layer(n_classes=k))
    assert classes.count == k
    assert np.sum(classes.number_fraction) == pytest.approx(1.0, abs=1e-14)
    assert np.all(classes.number_fraction > 0.0)
    # symmetric in t, so symmetric across classes
    np.testing.assert_allclose(classes.number_fraction, classes.number_fraction[::-1], rtol=1e-12)
    edges = np.linspace(-3.0, 3.0, k + 1)
    expected = np.diff(norm.cdf(edges)) / (norm.cdf(3.0) - norm.cdf(-3.0))
    np.testing.assert_allclose(classes.number_fraction, expected, rtol=1e-12)


def test_truncated_mass_is_0_27_percent():
    assert 1.0 - (norm.cdf(3.0) - norm.cdf(-3.0)) == pytest.approx(0.0027, abs=5e-5)


def test_class_centres_are_geometric_means_of_bounds():
    d50, sigma, k = 0.9, 0.45, 6
    classes = size_classes(_layer(d50=d50, sigma=sigma, n_classes=k))
    c = np.arange(1, k + 1)
    ln_lo = np.log(d50) + sigma * (-3.0 + 6.0 * (c - 1) / k)
    ln_hi = np.log(d50) + sigma * (-3.0 + 6.0 * c / k)
    np.testing.assert_allclose(classes.ecd_um, np.exp((ln_lo + ln_hi) / 2.0), rtol=1e-14)
    assert np.all(np.diff(classes.ecd_um) > 0.0)
    # with an even K the distribution median is a class bound: the middle centres straddle it
    assert np.sqrt(classes.ecd_um[2] * classes.ecd_um[3]) == pytest.approx(d50, rel=1e-14)


def test_conventional_geometry():
    classes = size_classes(_layer(morphology="conventional"))
    d = classes.ecd_um
    np.testing.assert_allclose(classes.volume_um3, np.pi * d**3 / 6.0, rtol=1e-15)
    np.testing.assert_allclose(classes.surface_um2, np.pi * d**2, rtol=1e-15)
    np.testing.assert_allclose(classes.projected_um2, np.pi * d**2 / 4.0, rtol=1e-15)
    np.testing.assert_array_equal(classes.thickness_um, d)


def test_cubic_is_treated_as_conventional():
    cubic = size_classes(_layer(morphology="cubic"))
    sphere = size_classes(_layer(morphology="conventional"))
    assert cubic.morphology == "cubic"
    np.testing.assert_array_equal(cubic.volume_um3, sphere.volume_um3)
    np.testing.assert_array_equal(cubic.surface_um2, sphere.surface_um2)


def test_tabular_geometry():
    t = 0.12
    classes = size_classes(_layer(morphology="tabular", d50=1.0, thickness=t, iodide=0.03))
    d = classes.ecd_um
    np.testing.assert_allclose(classes.volume_um3, np.pi * d**2 * t / 4.0, rtol=1e-15)
    np.testing.assert_allclose(classes.surface_um2, np.pi * d**2 / 2.0 + np.pi * d * t, rtol=1e-15)
    np.testing.assert_allclose(classes.projected_um2, np.pi * d**2 / 4.0, rtol=1e-15)
    np.testing.assert_array_equal(classes.thickness_um, np.full(6, t))
    # thin tabular grains carry more surface per volume than spheres of the same ECD
    sphere = size_classes(_layer(d50=1.0, iodide=0.03))
    assert np.all(classes.surface_um2 / classes.volume_um3 > sphere.surface_um2 / sphere.volume_um3)


def test_silver_per_grain_and_counts():
    x = 0.04
    classes = size_classes(_layer(iodide=x))
    m_agx = (1 - x) * M_AGBR + x * M_AGI
    np.testing.assert_allclose(
        classes.silver_g_per_grain, 1e-12 * RHO_AGX * classes.volume_um3 * M_AG / m_agx, rtol=1e-15
    )
    n_tot = np.sum(classes.number_per_um2)
    np.testing.assert_allclose(classes.number_per_um2, classes.number_fraction * n_tot, rtol=1e-14)
    assert classes.iodide_fraction == x
    assert classes.silver_coverage_g_m2 == 4.0


def test_fine_classes_approach_the_continuous_mean_volume():
    # The class sum is a quadrature of the truncated log-normal mean volume.
    d50, sigma = 0.9, 0.45
    classes = size_classes(_layer(d50=d50, sigma=sigma, n_classes=64))
    mean_volume = np.sum(classes.number_fraction * classes.volume_um3)

    def integrand(t):
        return np.pi / 6.0 * np.exp(3.0 * (np.log(d50) + sigma * t)) * norm.pdf(t)

    exact = quad(integrand, -3.0, 3.0)[0] / (norm.cdf(3.0) - norm.cdf(-3.0))
    assert mean_volume == pytest.approx(exact, rel=2e-3)


def test_n_classes_override():
    layer = _layer(n_classes=6)
    assert size_classes(layer).count == 6
    assert size_classes(layer, 12).count == 12


def test_arrays_read_only_and_float64():
    classes = size_classes(_layer())
    assert isinstance(classes, SizeClasses)
    for name in (
        "ecd_um", "thickness_um", "volume_um3", "surface_um2", "projected_um2",
        "number_fraction", "number_per_um2", "silver_g_per_grain",
    ):  # fmt: skip
        array = getattr(classes, name)
        assert array.dtype == np.float64 and array.shape == (6,)
        with pytest.raises(ValueError, match="read-only"):
            array[0] = 1.0
    with pytest.raises(AttributeError):
        classes.iodide_fraction = 0.5  # frozen dataclass


def test_missing_iodide_counts_as_pure_bromide():
    assert size_classes(_layer(iodide=None)).iodide_fraction == 0.0


def test_fixture_film_classes():
    material = Library(FIXTURES).material("test/bw-mechanistic")
    (layer,) = material.emulsion_layers
    classes = size_classes(layer)
    assert classes.morphology == "tabular"
    assert classes.count == 6
    assert 2.0 < classes.dry_thickness_um < 10.0


@pytest.mark.parametrize(
    ("layer_kwargs", "n_classes", "match"),
    [
        ({"gelatin": None}, None, "gelatin_coverage is required"),
        ({}, 1, "n_classes must be an integer from 2 to 64"),
        ({}, 65, "n_classes"),
        ({}, 2.5, "n_classes"),
        ({}, True, "n_classes"),
        ({"d50": 0.0}, None, "mean_diameter must be positive"),
        ({"sigma": 0.0}, None, "sigma_log must be positive"),
        ({"silver": 0.0}, None, "silver_coverage must be positive"),
    ],
)
def test_size_classes_rejects(layer_kwargs, n_classes, match):
    with pytest.raises(ValueError, match=match):
        size_classes(_layer(**layer_kwargs), n_classes)


def test_size_classes_needs_grain():
    with pytest.raises(ValueError, match="grain is required"):
        size_classes(EmulsionLayer(id="emulsion"))
