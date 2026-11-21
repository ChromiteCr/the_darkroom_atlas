import numpy as np
import pytest

from darkroom_atlas.readout.products import (
    TYPE2_ANCHORS,
    scan_product,
    type2_log_product,
    type2_product,
    visual_product,
)
from darkroom_atlas.spectral import DEFAULT_GRID, SpectralGrid, cmf_xyz, illuminant

_PRODUCTS = (visual_product, type2_product, scan_product)


def test_visual_product_matches_iso_anchors():
    # CAL 1.4: ISO 5-3:1995 Table 2 log10 of the visual product, peak-normalised at 570 nm.
    iso = {450.0: 3.090, 500.0: 4.276, 560.0: 4.989, 600.0: 4.902, 650.0: 4.238}
    wl = DEFAULT_GRID.wavelengths
    log_p = np.log10(visual_product(DEFAULT_GRID))
    log_p = log_p - log_p[wl == 570.0][0] + 5.0
    for nm, value in iso.items():
        assert abs(log_p[wl == nm][0] - value) <= 0.002


def test_type2_anchors():
    for nm, value in TYPE2_ANCHORS:
        assert abs(float(type2_log_product(nm)) - value) <= 1e-12
    # On the grid the anchors inside [380, 540] nm are reproduced relative to the 430 nm peak.
    wl = DEFAULT_GRID.wavelengths
    p = type2_product(DEFAULT_GRID)
    for nm, value in TYPE2_ANCHORS[1:]:
        rel = np.log10(p[wl == nm][0] / p[wl == 430.0][0]) + 5.0
        assert abs(rel - value) <= 1e-12
    assert p[wl == 430.0][0] == p.max()


def test_type2_support_and_monotone():
    wl = DEFAULT_GRID.wavelengths
    p = type2_product(DEFAULT_GRID)
    assert np.all(p[(wl < 380.0) | (wl > 540.0)] == 0.0)
    inside = (wl >= 380.0) & (wl <= 540.0)
    assert np.all(p[inside] > 0.0)
    rise = inside & (wl <= 430.0)
    fall = inside & (wl >= 430.0)
    assert np.all(np.diff(p[rise]) > 0.0)
    assert np.all(np.diff(p[fall]) < 0.0)
    fine = SpectralGrid(360.0, 540.0, 1.0)
    p_fine = type2_product(fine)
    assert p_fine[fine.wavelengths < 380.0].sum() == 0.0


@pytest.mark.parametrize("product", _PRODUCTS)
def test_products_normalised_and_read_only(product):
    p = product(DEFAULT_GRID)
    assert p.shape == (DEFAULT_GRID.n_bands,)
    assert p.dtype == np.float64
    assert np.all(p >= 0.0)
    assert abs(p.sum() - 1.0) <= 1e-12
    assert not p.flags.writeable


def test_visual_and_scan_definitions():
    ybar = cmf_xyz(DEFAULT_GRID)[:, 1]
    for product, name in ((visual_product, "A"), (scan_product, "D50")):
        raw = illuminant(name, DEFAULT_GRID) * ybar
        np.testing.assert_allclose(product(DEFAULT_GRID), raw / raw.sum(), rtol=1e-14)


def test_neutral_film_equal_densities():
    t = np.array([1.0, 0.5, 0.1, 1e-3])
    expected = -np.log10(t)
    for product in _PRODUCTS:
        p = product(DEFAULT_GRID)
        d = -np.log10((t[:, None] * p).sum(axis=-1) / p.sum())
        np.testing.assert_allclose(d, expected, rtol=0, atol=1e-12)


def test_type2_grid_without_support_raises():
    with pytest.raises(ValueError):
        type2_product(SpectralGrid(600.0, 780.0, 5.0))
