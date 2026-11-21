"""Exposure fields and the uniform one-column path (spec §2.1.4-§2.1.7, tests of §2.1.10)."""

import math
from pathlib import Path

import numpy as np
import pytest
from scipy import signal

from darkroom_atlas.data.library import Library
from darkroom_atlas.emulsion import size_classes
from darkroom_atlas.light import Light, LowRankLight
from darkroom_atlas.optics.expose import (
    _uniform,
    depth_weights,
    exposure_fields,
    uniform_exposure,
)
from darkroom_atlas.optics.halation import halation_fraction, halation_kernel
from darkroom_atlas.optics.kubelka import depth_response
from darkroom_atlas.optics.psf import scatter_kernel
from darkroom_atlas.optics.stack import optical_stack
from darkroom_atlas.spectral import DEFAULT_GRID, SpectralGrid
from darkroom_atlas.state import FrameGeometry

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "library"
TIME_S = 0.01
B = DEFAULT_GRID.n_bands


@pytest.fixture(scope="module")
def film():
    material = Library(FIXTURES).material("test/bw-mechanistic")
    classes = size_classes(material.emulsion_layers[0])
    stack = optical_stack(material, "135", classes, DEFAULT_GRID)
    return material, classes, stack


@pytest.fixture(scope="module")
def stack(film):
    return film[2]


def _response(stack, n_z=2):
    return depth_response(stack, n_z)


def _basis(seed=1, rank=3):
    rng = np.random.default_rng(seed)
    return rng.uniform(0.2, 1.0, (rank, B)) * 50.0


def _random_light(shape, seed=0, rank=3) -> LowRankLight:
    rng = np.random.default_rng(seed)
    coefficients = rng.uniform(0.1, 1.0, (*shape, rank))
    return LowRankLight(coefficients, _basis(seed + 1, rank), DEFAULT_GRID)


def _rel(a, b):
    return float(np.max(np.abs(a - b) / np.abs(b)))


def _expected_global(light: Light, stack, response, pitch):
    """Whole-frame linear convolution (zero fill) with one FFT per channel."""
    plus, minus = depth_weights(stack, response)
    flux = light.photon_flux * TIME_S
    psf = scatter_kernel(pitch, stack.psf_weight, stack.psf_length_um)
    hal = halation_kernel(pitch, stack.base).kernel
    out = np.zeros((*light.shape_hw, plus.shape[1]))
    for m in range(plus.shape[1]):
        a = np.einsum("hwb,b->hw", flux, plus[:, m])
        c = np.einsum("hwb,b->hw", flux, minus[:, m])
        out[..., m] = signal.fftconvolve(a, psf, mode="same") + signal.fftconvolve(
            c, hal, mode="same"
        )
    n_z = response.n_z
    return out[..., :n_z], out[..., n_z:]


# ------------------------------------------------------------------------------------------
# Tests listed in spec §2.1.10
# ------------------------------------------------------------------------------------------


@pytest.mark.parametrize("n_z", [1, 2, 4])
def test_low_rank_projection_matches_full(stack, n_z):
    # 2R < 4 n_z convolves the coefficients (n_z = 2, 4); n_z = 1 projects first.
    response = _response(stack, n_z)
    light = _random_light((64, 48))
    geometry = FrameGeometry((64, 48), 15.0)
    e_v, e_s, _ = exposure_fields(light, TIME_S, stack, response, geometry)
    f_v, f_s, _ = exposure_fields(light.to_light(), TIME_S, stack, response, geometry)
    assert e_v.shape == e_s.shape == (64, 48, n_z)
    assert _rel(e_v, f_v) <= 1e-12
    assert _rel(e_s, f_s) <= 1e-12


def test_expose_linear_in_light(stack):
    response = _response(stack)
    light = _random_light((60, 50))
    doubled = LowRankLight(light.coefficients * 2.0, light.basis, DEFAULT_GRID)
    geometry = FrameGeometry((60, 50), 15.0)
    e_v, e_s, _ = exposure_fields(light, TIME_S, stack, response, geometry)
    d_v, d_s, _ = exposure_fields(doubled, TIME_S, stack, response, geometry)
    assert _rel(d_v, 2.0 * e_v) <= 1e-12
    assert _rel(d_s, 2.0 * e_s) <= 1e-12


def _uniform_frame(stack, shape=(400, 300), pitch=20.0, block_px=128, n_z=2):
    response = _response(stack, n_z)
    basis = _basis()
    light = LowRankLight(np.ones((*shape, 3)), basis, DEFAULT_GRID)
    geometry = FrameGeometry(shape, pitch)
    e_v, e_s, meta = exposure_fields(light, TIME_S, stack, response, geometry, block_px=block_px)
    h = meta["halo_px"]
    interior = (slice(h + 1, shape[0] - h - 1), slice(h + 1, shape[1] - h - 1))
    photons = basis.sum(axis=0) * TIME_S
    return response, e_v, e_s, meta, interior, photons


def test_expose_uniform_interior(stack):
    _, e_v, e_s, meta, interior, _ = _uniform_frame(stack)
    assert meta["halo_px"] == math.ceil(meta["halation_radius_um"] / 20.0)
    for field in (e_v, e_s):
        inner = field[interior]
        assert inner.shape[0] > 128 and inner.shape[1] > 0  # spans several blocks
        assert _rel(inner, np.broadcast_to(inner[:1, :1], inner.shape)) <= 1e-12


def test_uniform_frame_equals_uniform_exposure(stack):
    response, e_v, e_s, meta, interior, photons = _uniform_frame(stack)
    kernel_sum = meta["halation_fraction"] - meta["halation_truncated"]
    # The image path's halation kernel holds f_H - truncated; with that weight the one-column
    # path reproduces the interior exactly.
    u_v, u_s = _uniform(photons, stack, response, kernel_sum)
    assert _rel(e_v[interior], np.broadcast_to(u_v, e_v[interior].shape)) <= 1e-12
    assert _rel(e_s[interior], np.broadcast_to(u_s, e_s[interior].shape)) <= 1e-12
    # uniform_exposure uses the untruncated f_H (spec §6.1): it exceeds the interior by
    # exactly the truncated tail of the halation.
    v, s = uniform_exposure(photons, stack, response)
    _, minus = depth_weights(stack, response)
    tail = meta["halation_truncated"] * (photons @ minus)
    n_z = response.n_z
    np.testing.assert_allclose(v - e_v[interior][0, 0], tail[:n_z], rtol=1e-6, atol=1e-12 * v.max())
    np.testing.assert_allclose(s - e_s[interior][0, 0], tail[n_z:], rtol=1e-6, atol=1e-12 * s.max())
    assert _rel(e_v[interior][0, 0], v) <= 1e-3
    assert _rel(e_s[interior][0, 0], s) <= 1e-3


@pytest.mark.parametrize("block_px", [64, 512])
def test_block_convolution_equals_global_fft(stack, block_px):
    pitch = 20.0
    response = _response(stack)
    light = _random_light((300, 200), seed=4)
    geometry = FrameGeometry((300, 200), pitch)
    e_v, e_s, meta = exposure_fields(light, TIME_S, stack, response, geometry, block_px=block_px)
    radius_px = halation_kernel(pitch, stack.base).kernel.shape[0] // 2
    assert meta["halo_px"] >= radius_px
    assert meta["block_px"] == block_px
    g_v, g_s = _expected_global(light.to_light(), stack, response, pitch)
    assert _rel(e_v, g_v) <= 1e-12
    assert _rel(e_s, g_s) <= 1e-12


# ------------------------------------------------------------------------------------------
# Further limits
# ------------------------------------------------------------------------------------------


def test_photon_conservation(film):
    # A lit square well inside the frame: every scattered and returned photon lands in the
    # frame, so the photons absorbed by the grains equal sum_b N_b [A+_b + tau_b sum(H) A-_b].
    _, classes, stack = film
    pitch, n_z = 20.0, 2
    response = _response(stack, n_z)
    shape = (300, 300)
    coefficients = np.zeros((*shape, 3))
    coefficients[140:160, 145:155] = np.random.default_rng(7).uniform(0.1, 1.0, (20, 10, 3))
    basis = _basis()
    light = LowRankLight(coefficients, basis, DEFAULT_GRID)
    e_v, e_s, meta = exposure_fields(light, TIME_S, stack, response, FrameGeometry(shape, pitch))
    assert 140 - meta["halo_px"] > 0 and 160 + meta["halo_px"] < 300
    cell = classes.dry_thickness_um / n_z
    nu = classes.number_per_um2
    # Photons per um^2 absorbed in cell z: (nu_c / n_z)(V_c e_V + S_c e_S) summed over c.
    per_area = np.sum(nu * classes.volume_um3) * e_v + np.sum(nu * classes.surface_um2) * e_s
    absorbed = np.sum(per_area) / n_z * pitch**2

    photons = coefficients.reshape(-1, 3).sum(axis=0) @ basis * TIME_S * pitch**2  # (B,)
    up_absorbed = np.sum(stack.mu_a_per_um * response.minus, axis=0) * cell
    down_absorbed = np.sum(stack.mu_a_per_um * response.plus, axis=0) * cell
    np.testing.assert_allclose(down_absorbed, response.absorbed, rtol=1e-12)
    kernel_sum = meta["halation_fraction"] - meta["halation_truncated"]
    expected = np.sum(photons * (response.absorbed + response.into_base * kernel_sum * up_absorbed))
    assert abs(absorbed - expected) / expected <= 1e-12


def test_halation_ring_position(stack):
    pitch = 10.0
    response = _response(stack)
    n = 401
    coefficients = np.zeros((n, n, 1))
    coefficients[n // 2, n // 2, 0] = 1.0
    light = LowRankLight(coefficients, _basis(rank=1), DEFAULT_GRID)
    e_v, e_s, meta = exposure_fields(light, TIME_S, stack, response, FrameGeometry((n, n), pitch))
    assert n // 2 > meta["halo_px"]
    y, x = np.mgrid[:n, :n] - n // 2
    radius = np.hypot(y, x) * pitch
    r_c = stack.base.critical_radius_um
    assert r_c == pytest.approx(2 * 127 * math.tan(math.asin(1 / 1.48)), rel=1e-14)
    outside_psf = radius > 3.0 * meta["psf_radius_um"]
    for field in (e_v[..., -1], e_s[..., -1]):
        ring = np.where(outside_psf, field, 0.0)
        peak = radius.flat[np.argmax(ring)]
        assert abs(peak - r_c) <= pitch
        # Inside r_c only the Fresnel part returns (a few per cent).
        inner = field[(radius > 0.5 * r_c) & (radius < r_c - 3 * pitch)]
        assert inner.max() < 0.2 * ring.max()
        # Nothing beyond the halation radius.
        assert np.all(field[radius > meta["halation_radius_um"] + 2 * pitch] < 1e-12 * ring.max())


def test_edge_spread(stack):
    # A straight edge: away from the top and bottom frame edges each row is the 1-D
    # convolution of the step with the column marginals of the two kernels.
    pitch = 20.0
    response = _response(stack)
    shape, x0 = (220, 400), 200
    coefficients = np.zeros((*shape, 1))
    coefficients[:, x0:, 0] = 1.0
    basis = _basis(rank=1)
    light = LowRankLight(coefficients, basis, DEFAULT_GRID)
    e_v, e_s, meta = exposure_fields(
        light, TIME_S, stack, response, FrameGeometry(shape, pitch), block_px=96
    )
    h = meta["halo_px"]
    rows = slice(h + 1, shape[0] - h - 1)
    assert rows.stop > rows.start
    step = np.zeros(shape[1])
    step[x0:] = 1.0
    psf = scatter_kernel(pitch, stack.psf_weight, stack.psf_length_um).sum(axis=0)
    hal = halation_kernel(pitch, stack.base).kernel.sum(axis=0)
    scattered = np.convolve(step, psf, mode="same")
    returned = np.convolve(step, hal, mode="same")
    plus, minus = depth_weights(stack, response)
    q_plus = TIME_S * (basis @ plus)[0]
    q_minus = TIME_S * (basis @ minus)[0]
    expected = scattered[:, None] * q_plus + returned[:, None] * q_minus  # (W, 2 n_z)
    n_z = response.n_z
    for field, cols in ((e_v, slice(0, n_z)), (e_s, slice(n_z, 2 * n_z))):
        want = np.broadcast_to(expected[:, cols], field[rows].shape)
        assert np.max(np.abs(field[rows] - want)) <= 1e-12 * want.max()
        # Dark side beyond the halo: no light at all; the spread is symmetric about the edge.
        assert np.all(field[rows, : x0 - h - 1] <= 1e-12 * want.max())
        far_col = (x0 + h + 1 + shape[1] - h - 1) // 2  # lit, farther than the halo from edges
        assert x0 + h < far_col < shape[1] - h - 1
        far = field[rows, far_col]
        np.testing.assert_allclose(field[rows, x0 - 5] + field[rows, x0 + 4], far, rtol=1e-12)


# ------------------------------------------------------------------------------------------
# Bookkeeping and bad input
# ------------------------------------------------------------------------------------------


def test_exposure_fields_meta_and_read_only(stack):
    light = _random_light((20, 30))
    e_v, e_s, meta = exposure_fields(
        light, TIME_S, stack, _response(stack), FrameGeometry((20, 30), 15.0)
    )
    assert set(meta) == {
        "halation_fraction",
        "halation_radius_um",
        "halation_truncated",
        "psf_radius_um",
        "block_px",
        "halo_px",
    }
    assert meta["halation_fraction"] == pytest.approx(halation_fraction(stack.base), rel=1e-15)
    # About energy_tol * f_H; pixel sampling at R_H adds a little.
    assert 0.0 < meta["halation_truncated"] < 2e-3 * meta["halation_fraction"]
    assert meta["psf_radius_um"] == pytest.approx(10.65, abs=0.01)
    assert meta["block_px"] == 512
    for field in (e_v, e_s):
        assert field.dtype == np.float64
        assert not field.flags.writeable
        assert np.all(field > 0.0)


def test_uniform_exposure_batches_and_validates(stack):
    response = _response(stack, 4)
    photons = _basis()  # (3, B)
    v, s = uniform_exposure(photons, stack, response)
    assert v.shape == s.shape == (3, 4)
    v0, s0 = uniform_exposure(photons[0], stack, response)
    np.testing.assert_array_equal(v[0], v0)
    np.testing.assert_array_equal(s[0], s0)
    with pytest.raises(ValueError, match="photons_per_band must have shape"):
        uniform_exposure(np.ones(B - 1), stack, response)
    with pytest.raises(ValueError, match="non-negative"):
        uniform_exposure(-np.ones(B), stack, response)
    with pytest.raises(TypeError, match="OpticalStack"):
        uniform_exposure(np.ones(B), None, response)


def test_exposure_fields_rejects_bad_input(stack):
    response = _response(stack)
    light = _random_light((20, 30))
    geometry = FrameGeometry((20, 30), 15.0)
    with pytest.raises(ValueError, match="geometry has shape"):
        exposure_fields(light, TIME_S, stack, response, FrameGeometry((30, 20), 15.0))
    with pytest.raises(ValueError, match="whole frame"):
        window = FrameGeometry((20, 30), 15.0, origin_px=(1, 0), frame_shape_hw=(21, 30))
        exposure_fields(light, TIME_S, stack, response, window)
    with pytest.raises(ValueError, match="time_s"):
        exposure_fields(light, 0.0, stack, response, geometry)
    with pytest.raises(ValueError, match="block_px"):
        exposure_fields(light, TIME_S, stack, response, geometry, block_px=0)
    with pytest.raises(ValueError, match="workers"):
        exposure_fields(light, TIME_S, stack, response, geometry, workers=0)
    with pytest.raises(TypeError, match="LightField"):
        exposure_fields(np.ones((20, 30, B)), TIME_S, stack, response, geometry)
    other_grid = SpectralGrid(380.0, 780.0, 10.0)
    other = LowRankLight(np.ones((20, 30, 1)), np.ones((1, other_grid.n_bands)), other_grid)
    with pytest.raises(ValueError, match="grid"):
        exposure_fields(other, TIME_S, stack, response, geometry)
