import array_api_strict as xps
import numpy as np
import pytest

from darkroom_atlas.light import Light, apply_transmittance
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
