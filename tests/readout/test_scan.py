"""Standard scan readout (spec D1 2.8): fixed inversion, grain on and off."""

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

from darkroom_atlas.pipeline import readout
from darkroom_atlas.readout import scan
from darkroom_atlas.readout.densitometer import density_frame
from darkroom_atlas.readout.products import scan_product


def _load_e2e():
    name = "d1d_end_to_end_bw"
    if name not in sys.modules:
        path = Path(__file__).resolve().parents[1] / "test_end_to_end_bw.py"
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


_e2e = _load_e2e()


@pytest.fixture(scope="module")
def developed(tmp_path_factory):
    library = _e2e.mechanistic_library(tmp_path_factory.mktemp("lib") / "library")
    return _e2e.developed_state(library)


def test_scan_constants():
    assert scan.SCAN_G_REF == 0.56
    assert scan.SCAN_Y_REF == 0.18
    assert scan.SCAN_D_REF == pytest.approx(0.10 + 0.56 * math.log10(12.5), abs=1e-15)
    assert round(scan.SCAN_D_REF, 3) == 0.714


def test_scan_inversion_fixed_point():
    for base in (0.0, 0.10, 0.25):
        y = scan.scan_inversion(np.array([base + scan.SCAN_D_REF]), base)
        assert abs(y[0] - 0.18) <= 1e-12


def test_scan_inversion_gradient_and_opaque():
    base = 0.1
    d = np.array([base + scan.SCAN_D_REF + scan.SCAN_G_REF, math.inf, 1e6])
    y = scan.scan_inversion(d, base)
    assert y[0] == pytest.approx(1.8, rel=1e-12)  # one G_ref of density is a decade
    assert y[1] == 1.0 and y[2] == 1.0


def test_read_scan_follows_the_rule(developed):
    state, ctx = developed
    for grain in (True, False):
        result = readout(state, "scan_standard", ctx, grain=grain)
        assert result.kind == "image_linear_srgb"
        assert result.data.shape == (*state.geometry.shape_hw, 3)
        d_scan, _ = density_frame(state, scan_product(state.grid), grain=grain)
        expected = scan.scan_inversion(d_scan, scan.base_density(state))
        for channel in range(3):
            assert np.array_equal(result.data[..., channel], expected)
        assert result.meta["base_density"] == pytest.approx(0.10)  # fixture 135 base
        assert ("grain" in result.meta) is grain


def test_scan_default_has_grain(developed):
    state, ctx = developed
    default = readout(state, "scan_standard", ctx)
    sampled = readout(state, "scan_standard", ctx, grain=True)
    mean = readout(state, "scan_standard", ctx, grain=False)
    assert np.array_equal(default.data, sampled.data)
    assert not np.array_equal(sampled.data, mean.data)
    assert sampled.meta["grain"]["mean_only_classes"] is not None


def test_scan_is_seeded(developed, tmp_path):
    state, ctx = developed
    a = readout(state, "scan_standard", ctx).data
    library = _e2e.mechanistic_library(tmp_path / "library")
    other, other_ctx = _e2e.developed_state(library, seed=1)
    b = readout(other, "scan_standard", other_ctx).data
    assert not np.array_equal(a, b)


def test_scan_opaque_pixels_counts_every_replaced_pixel(developed, monkeypatch):
    state, ctx = developed
    h, w = state.geometry.shape_hw
    d = np.full((h, w), 0.5)
    d[0, 0] = math.inf  # non-finite density
    d[0, 1] = 1e6  # finite density whose Y_out overflows
    monkeypatch.setattr(scan, "density_frame", lambda *a, **k: (d, None))
    result = readout(state, "scan_standard", ctx, grain=False)
    assert result.meta["opaque_pixels"] == 2
    assert result.data[0, 0, 0] == 1.0 and result.data[0, 1, 0] == 1.0


def test_bad_scan_parameters(developed):
    state, ctx = developed
    with pytest.raises(TypeError, match="grain"):
        readout(state, "scan_standard", ctx, grain="yes")
    with pytest.raises(TypeError):
        readout(state, "scan_standard", ctx, aperture_um=48.0)
