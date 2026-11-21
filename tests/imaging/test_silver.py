"""Silver image: Nutting, covering power and spectral transmittance (spec D1 2.5.4)."""

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

from darkroom_atlas.emulsion.classes import size_classes
from darkroom_atlas.imaging import (
    TransmittanceFactors,
    covering_powers,
    developed_areas_um2,
    imaging_constants,
    product_density,
    silver_density,
    transmittance_factors,
)
from darkroom_atlas.readout.products import type2_product, visual_product
from darkroom_atlas.spectral import DEFAULT_GRID


def _load_synthetic():
    name = "d1d_synthetic"
    if name not in sys.modules:
        path = Path(__file__).resolve().parents[1] / "imaging" / "synthetic.py"
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


syn = _load_synthetic()
K = 6
B = DEFAULT_GRID.n_bands


@pytest.fixture(scope="module")
def film():
    return syn.material()


def _random_state(film, shape=(7, 9), seed=0, **kw):
    rng = np.random.default_rng(seed)
    f = rng.uniform(0.0, 1.0, (*shape, K))
    x = f * rng.uniform(0.0, 1.0, (*shape, K))
    y = x * rng.uniform(0.5, 1.0, (*shape, K))
    phys = rng.uniform(0.0, 0.2, shape)
    halide = rng.uniform(0.0, 2.0, shape)
    return syn.state(
        film,
        shape,
        5.0,
        k=K,
        initiated=f,
        developed=x,
        developed_m2=y,
        physical=phys,
        halide=halide,
        **kw,
    )


def test_covering_power_law(film):
    layer = film.emulsion_layers[0]
    classes = size_classes(layer, K)
    cp = covering_powers(classes, layer.development)
    # CP_c = CP_ref (d_c / 1 um)^-0.5 with the fixture's CP_ref = 0.6 and beta = 0.5.
    np.testing.assert_allclose(cp, 0.6 * classes.ecd_um**-0.5, rtol=1e-15)
    assert not cp.flags.writeable
    a = developed_areas_um2(classes, cp)
    np.testing.assert_allclose(
        a, math.log(10.0) * 1e12 * cp * classes.silver_g_per_grain, rtol=1e-15
    )
    assert not a.flags.writeable


def test_nutting_linear_in_silver(film):
    s1 = _random_state(film)
    d1 = silver_density(s1)
    alpha = 1.7
    s2 = s1.with_fields(
        {
            "emulsion/silver": alpha * s1.get("emulsion", "silver"),
            "emulsion/silver_physical": alpha * s1.get("emulsion", "silver_physical"),
        }
    )
    d2 = silver_density(s2)
    assert np.max(np.abs(d2 - alpha * d1) / (alpha * d1)) <= 1e-12
    # Additivity of two silver images.
    s3 = s1.with_fields(
        {
            "emulsion/silver": s1.get("emulsion", "silver") + s2.get("emulsion", "silver"),
            "emulsion/silver_physical": s1.get("emulsion", "silver_physical")
            + s2.get("emulsion", "silver_physical"),
        }
    )
    d3 = silver_density(s3)
    assert np.max(np.abs(d3 - (d1 + d2)) / d3) <= 1e-12


def test_covering_power_area_identity(film):
    s = _random_state(film)
    const = imaging_constants(s)
    nu = const.classes.number_per_um2
    silver = s.get("emulsion", "silver")
    x = silver / const.class_silver_g_m2
    lhs = np.zeros(silver.shape[:2])
    rhs = np.zeros(silver.shape[:2])
    for c in range(K):
        lhs = lhs + nu[c] * const.developed_areas_um2[c] * x[..., c]
        rhs = rhs + const.covering_power[c] * silver[..., c]
    lhs = lhs / math.log(10.0)  # 0.4343 = 1 / ln 10
    assert np.max(np.abs(lhs - rhs) / rhs) <= 1e-12


def test_silver_density_window(film):
    s = _random_state(film)
    full = silver_density(s)
    part = silver_density(s, slice(2, 5), slice(1, 8))
    assert np.array_equal(part, full[2:5, 1:8])
    assert not full.flags.writeable
    # The density adds r_CP CP_ref silver_physical (0.5 * 0.6 in the fixture).
    no_phys = s.with_fields({"emulsion/silver_physical": np.zeros((7, 9))})
    np.testing.assert_allclose(
        full - silver_density(no_phys), 0.3 * s.get("emulsion", "silver_physical"), rtol=1e-12
    )


def test_transmittance_factorisation(film):
    neutral = syn.material(development={"dye_density": 0.0})
    s = _random_state(neutral)
    rows, cols = slice(0, 7), slice(0, 9)
    factors = transmittance_factors(s, rows, cols)
    t_ag = 10.0 ** -silver_density(s)
    d_visual = product_density(t_ag, factors, visual_product(DEFAULT_GRID))
    d_type2 = product_density(t_ag, factors, type2_product(DEFAULT_GRID))
    assert np.max(np.abs(d_visual - d_type2)) <= 1e-12
    # T_b does not depend on b: a one-band product at any band gives the same density.
    expected = silver_density(s) + factors.halide_density + 0.10
    for b in (0, 17, 40, 80):
        delta = np.zeros(B)
        delta[b] = 1.0
        d_b = product_density(t_ag, factors, delta)
        assert np.max(np.abs(d_b - expected)) <= 1e-12
    assert np.max(np.abs(d_visual - expected)) <= 1e-12


def test_transmittance_factors_fields(film):
    s = _random_state(film, film_format="120")
    factors = transmittance_factors(s, slice(1, 4), slice(2, 6))
    assert isinstance(factors, TransmittanceFactors)
    assert factors.base_density.shape == (B,)
    np.testing.assert_array_equal(factors.base_density, 0.0)  # 120 grey density 0, backing removed
    np.testing.assert_array_equal(factors.dye_density, 0.02)
    np.testing.assert_allclose(
        factors.halide_density, 0.3 * s.get("emulsion", "halide")[1:4, 2:6], rtol=1e-15
    )
    assert factors.dye_residual.shape == (3, 4)
    s135 = _random_state(film)
    np.testing.assert_array_equal(
        transmittance_factors(s135, slice(0, 1), slice(0, 1)).base_density, 0.1
    )


def test_product_density_with_dye(film):
    # A non-uniform residual dye with a band-dependent dye density, against the direct sum.
    s = _random_state(film)
    rng = np.random.default_rng(3)
    residual = rng.uniform(0.0, 1.0, (7, 9))
    s = s.with_fields({"emulsion/residual_dye": residual})
    factors = transmittance_factors(s, slice(0, 7), slice(0, 9))
    dye = np.linspace(0.0, 0.1, B)
    factors = TransmittanceFactors(
        factors.base_density, dye, factors.halide_density, factors.dye_residual
    )
    product = visual_product(DEFAULT_GRID)
    t_ag = 10.0 ** -silver_density(s)
    got = product_density(t_ag, factors, product)
    t_b = t_ag[..., None] * 10.0 ** -(
        factors.base_density + factors.halide_density[..., None] + dye * residual[..., None]
    )
    expected = -np.log10(np.sum(product * t_b, axis=-1) / np.sum(product))
    np.testing.assert_allclose(got, expected, rtol=1e-12, atol=1e-13)
    # Uniform residual dye: the same through the one-number path.
    s_u = s.with_fields({"emulsion/residual_dye": np.full((7, 9), 0.4)})
    f_u = transmittance_factors(s_u, slice(0, 7), slice(0, 9))
    f_u = TransmittanceFactors(f_u.base_density, dye, f_u.halide_density, f_u.dye_residual)
    t_b = t_ag[..., None] * 10.0 ** -(f_u.base_density + f_u.halide_density[..., None] + dye * 0.4)
    expected = -np.log10(np.sum(product * t_b, axis=-1) / np.sum(product))
    np.testing.assert_allclose(product_density(t_ag, f_u, product), expected, rtol=1e-12)


def test_product_density_zero_transmittance_is_inf(film):
    s = _random_state(film)
    factors = transmittance_factors(s, slice(0, 7), slice(0, 9))
    t = np.full((7, 9), 0.5)
    t[0, 0] = 0.0
    d = product_density(t, factors, visual_product(DEFAULT_GRID))
    assert np.isinf(d[0, 0]) and np.all(np.isfinite(d.ravel()[1:]))


def test_bad_input(film):
    s = _random_state(film)
    factors = transmittance_factors(s, slice(0, 7), slice(0, 9))
    product = visual_product(DEFAULT_GRID)
    with pytest.raises(ValueError, match="shape"):
        product_density(np.ones((7, 8)), factors, product)
    with pytest.raises(ValueError, match="product"):
        product_density(np.ones((7, 9)), factors, product[:-1])
    with pytest.raises(ValueError, match="weight"):
        product_density(np.ones((7, 9)), factors, np.zeros(B))
    with pytest.raises(ValueError, match="lie in"):
        product_density(np.full((7, 9), 1.5), factors, product)
    with pytest.raises(TypeError):
        product_density(np.ones((7, 9)), object(), product)
    bad = s.with_fields({"emulsion/silver": -s.get("emulsion", "silver")})
    with pytest.raises(ValueError, match="negative"):
        silver_density(bad)
    with pytest.raises(ValueError, match="empty"):
        silver_density(s, slice(3, 3), None)
    with pytest.raises(ValueError, match="step"):
        silver_density(s, slice(0, 7, 2), None)
    with pytest.raises(KeyError):
        silver_density(s.drop("emulsion/silver"))
    classes = size_classes(film.emulsion_layers[0], K)
    with pytest.raises(ValueError, match="shape"):
        developed_areas_um2(classes, np.ones(K + 1))
    with pytest.raises(ValueError, match="non-negative"):
        developed_areas_um2(classes, -np.ones(K))
    with pytest.raises(ValueError, match="film format"):
        transmittance_factors(_random_state(film, film_format="220"), slice(0, 1), slice(0, 1))
