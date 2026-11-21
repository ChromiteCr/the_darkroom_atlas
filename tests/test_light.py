import array_api_strict as xps
import numpy as np
import pytest

from darkroom_atlas.light import Light, LightField, LowRankLight, apply_transmittance
from darkroom_atlas.spectral import DEFAULT_GRID, SpectralGrid

B = DEFAULT_GRID.n_bands


def _band(wavelength_nm, grid=DEFAULT_GRID):
    return int(np.flatnonzero(grid.wavelengths == wavelength_nm)[0])


def test_luminous_exposure_anchor():
    flux = np.zeros((1, 1, B))
    flux[0, 0, _band(555.0)] = 4090.7
    light = Light(flux, DEFAULT_GRID)
    exposure = light.luminous_exposure_lux_s(1.0)
    assert exposure.shape == (1, 1)
    assert exposure[0, 0] == pytest.approx(1.000, rel=2e-3)


def test_luminous_exposure_scales_with_time():
    flux = np.full((2, 3, B), 10.0)
    light = Light(flux, DEFAULT_GRID)
    np.testing.assert_allclose(
        light.luminous_exposure_lux_s(0.5), 0.5 * light.luminous_exposure_lux_s(1.0), rtol=1e-15
    )
    with pytest.raises(ValueError):
        light.luminous_exposure_lux_s(0.0)


def test_filtered_scalar_band():
    flux = np.random.default_rng(1).random((4, 5, B)) * 100.0
    light = Light(flux, DEFAULT_GRID, meta={"iso": 100.0})
    half = light.filtered(np.full(B, 0.5))
    np.testing.assert_array_equal(half.photon_flux, 0.5 * flux)
    assert half.meta == {"iso": 100.0}
    assert half.grid == light.grid
    per_pixel = light.filtered(np.full((4, 5, B), 0.5))
    np.testing.assert_array_equal(per_pixel.photon_flux, 0.5 * flux)
    assert per_pixel.meta == {"iso": 100.0}


def test_filtered_rejects_out_of_range():
    light = Light(np.ones((2, 2, B)), DEFAULT_GRID)
    for bad in (1.2, -0.1, np.nan):
        t = np.full(B, 0.5)
        t[10] = bad
        with pytest.raises(ValueError):
            light.filtered(t)
    with pytest.raises(ValueError):
        light.filtered(np.full((2, 2, B), 1.2))
    with pytest.raises(ValueError):
        light.filtered(np.full(B - 1, 0.5))


def test_rejects_bad_flux():
    good = np.ones((2, 3, B))
    negative = good.copy()
    negative[0, 0, 0] = -1.0
    with pytest.raises(ValueError):
        Light(negative, DEFAULT_GRID)
    nan = good.copy()
    nan[1, 2, 3] = np.nan
    with pytest.raises(ValueError):
        Light(nan, DEFAULT_GRID)
    inf = good.copy()
    inf[1, 2, 3] = np.inf
    with pytest.raises(ValueError):
        Light(inf, DEFAULT_GRID)
    with pytest.raises(ValueError):
        Light(np.ones((2, B)), DEFAULT_GRID)
    with pytest.raises(ValueError):
        Light(np.ones((1, 2, 3, B)), DEFAULT_GRID)
    with pytest.raises(ValueError):
        Light(good, SpectralGrid(380.0, 780.0, 10.0))


def test_light_is_immutable():
    flux = np.ones((2, 3, B))
    light = Light(flux, DEFAULT_GRID, meta={"iso": 100.0})
    assert light.shape_hw == (2, 3)
    assert not light.photon_flux.flags.writeable
    flux[0, 0, 0] = 5.0  # the caller's array is copied, not shared
    assert light.photon_flux[0, 0, 0] == 1.0
    with pytest.raises(TypeError):
        light.meta["iso"] = 200.0


def test_filtered_array_api_strict():
    flux = xps.ones((2, 3, B), dtype=xps.float64)
    t = xps.full((B,), 0.5, dtype=xps.float64)
    out = apply_transmittance(flux, t)
    assert type(out) is type(flux)
    assert bool(xps.all(out == 0.5))
    per_pixel = apply_transmittance(flux, xps.full((2, 3, B), 0.25, dtype=xps.float64))
    assert bool(xps.all(per_pixel == 0.25))
    with pytest.raises(ValueError):
        apply_transmittance(flux, xps.full((B,), 1.5, dtype=xps.float64))


# --------------------------------------------------------------------------------------
# D1a: window / project and LowRankLight (spec D1 section 4.2, section 10)
# --------------------------------------------------------------------------------------


def _low_rank(h=6, w=5, rank=3, seed=3, meta=None):
    rng = np.random.default_rng(seed)
    coefficients = rng.random((h, w, rank)) * 2.0
    basis = rng.random((rank, B)) * 1e3
    return LowRankLight(coefficients, basis, DEFAULT_GRID, meta or {"iso": 100.0})


def _assert_rel(actual, expected, rtol=1e-12):
    actual, expected = np.asarray(actual), np.asarray(expected)
    assert actual.shape == expected.shape
    assert np.all(np.abs(actual - expected) <= rtol * np.abs(expected))


def test_light_window_shares_and_project():
    flux = np.random.default_rng(2).random((5, 7, B))
    light = Light(flux, DEFAULT_GRID, meta={"iso": 100.0})
    block = light.window(slice(1, 4), slice(2, 7))
    assert block.shape_hw == (3, 5)
    assert np.shares_memory(block.photon_flux, light.photon_flux)
    np.testing.assert_array_equal(block.photon_flux, flux[1:4, 2:7])
    assert block.meta == {"iso": 100.0}
    weights = np.random.default_rng(4).random((B, 2))
    np.testing.assert_array_equal(light.project(weights), flux @ weights)
    np.testing.assert_array_equal(
        light.project(weights, rows=slice(1, 4), cols=slice(2, 7)), flux[1:4, 2:7] @ weights
    )
    assert isinstance(light, LightField)


def test_low_rank_light_window_and_project():
    low = _low_rank()
    full = low.to_light()
    assert isinstance(low, LightField)
    assert full.shape_hw == low.shape_hw == (6, 5)
    _assert_rel(full.photon_flux, low.coefficients @ low.basis)
    weights = np.random.default_rng(5).random((B, 4))
    _assert_rel(low.project(weights), full.project(weights))
    rows, cols = slice(2, 5), slice(1, 4)
    block = low.window(rows, cols)
    _assert_rel(block.photon_flux, full.window(rows, cols).photon_flux)
    _assert_rel(low.project(weights, rows, cols), full.project(weights, rows, cols))
    _assert_rel(low.project(weights, rows=rows), full.project(weights, rows=rows))
    _assert_rel(low.luminous_exposure_lux_s(0.01), full.luminous_exposure_lux_s(0.01))
    assert block.meta == low.meta


def test_low_rank_light_memory():
    shape = (6000, 4000, 3)
    coefficients = np.zeros(shape)
    coefficients.flags.writeable = False
    low = LowRankLight(coefficients, np.ones((3, B)), DEFAULT_GRID)
    assert low.coefficients is coefficients  # read-only input is shared, not copied
    assert low.nbytes == coefficients.nbytes + 3 * B * 8
    assert low.nbytes <= 0.6e9
    assert low.rank == 3


def test_low_rank_light_filtered():
    low = _low_rank()
    t = np.linspace(0.0, 1.0, B)
    out = low.filtered(t)
    assert isinstance(out, LowRankLight)
    assert out.coefficients is low.coefficients
    np.testing.assert_array_equal(out.basis, low.basis * t)
    _assert_rel(out.to_light().photon_flux, low.to_light().filtered(t).photon_flux)
    assert out.meta == low.meta
    with pytest.raises(ValueError):
        low.filtered(np.full((6, 5, B), 0.5))
    with pytest.raises(ValueError):
        low.filtered(np.full(B, 1.5))
    bad = np.full(B, 0.5)
    bad[3] = np.nan
    with pytest.raises(ValueError):
        low.filtered(bad)


def test_low_rank_light_is_immutable():
    coefficients = np.ones((2, 3, 3))
    low = LowRankLight(coefficients, np.ones((3, B)), DEFAULT_GRID, meta={"iso": 100.0})
    assert not low.coefficients.flags.writeable
    assert not low.basis.flags.writeable
    coefficients[0, 0, 0] = 5.0  # a writable input is copied
    assert low.coefficients[0, 0, 0] == 1.0
    with pytest.raises(TypeError):
        low.meta["iso"] = 200.0


def test_low_rank_light_rejects_bad_input():
    good_c, good_b = np.ones((2, 3, 3)), np.ones((3, B))
    with pytest.raises(ValueError):
        LowRankLight(np.ones((2, 3)), good_b, DEFAULT_GRID)
    with pytest.raises(ValueError):
        LowRankLight(np.ones((0, 3, 3)), good_b, DEFAULT_GRID)
    with pytest.raises(ValueError):
        LowRankLight(good_c, np.ones((2, B)), DEFAULT_GRID)
    with pytest.raises(ValueError):
        LowRankLight(good_c, np.ones((3, B - 1)), DEFAULT_GRID)
    for bad in (-1.0, np.nan, np.inf):
        c = good_c.copy()
        c[0, 0, 0] = bad
        with pytest.raises(ValueError):
            LowRankLight(c, good_b, DEFAULT_GRID)
        b = good_b.copy()
        b[0, 0] = bad
        with pytest.raises(ValueError):
            LowRankLight(good_c, b, DEFAULT_GRID)
    with pytest.raises(TypeError):
        LowRankLight(good_c, good_b, "grid")


def test_window_and_project_reject_bad_arguments():
    low = _low_rank()
    for light in (low, low.to_light()):
        with pytest.raises(ValueError):
            light.window(slice(3, 3), slice(0, 2))
        with pytest.raises(TypeError):
            light.window(0, slice(0, 2))
        with pytest.raises(ValueError):
            light.project(np.ones(B))
        with pytest.raises(ValueError):
            light.project(np.ones((B - 1, 2)))
        with pytest.raises(ValueError):
            light.project(np.full((B, 2), np.nan))
        with pytest.raises(ValueError):
            light.project(np.ones((B, 2)), rows=slice(10, 12))
        with pytest.raises(TypeError):
            light.project(np.ones((B, 2)), rows=1)
        with pytest.raises(ValueError):
            light.luminous_exposure_lux_s(0.0)
