"""Scatter PSF kernel and transfer function (spec §2.1.6, tests of §2.1.10)."""

import math

import numpy as np
import pytest
from scipy import special

from darkroom_atlas.optics.psf import (
    scatter_half_width,
    scatter_kernel,
    scatter_otf,
    tail_radius_ratio,
)

PITCH = 2.0
LENGTHS = [0.3 * PITCH, PITCH, 4.0 * PITCH]


def _line_transfer(kernel: np.ndarray, nu: np.ndarray, pitch: float) -> np.ndarray:
    """Transfer function of a pixel kernel along x at nu_y = 0 (the kernel is symmetric)."""
    half = kernel.shape[0] // 2
    lsf = kernel.sum(axis=0)
    m = np.arange(-half, half + 1) * pitch
    return (lsf[None, :] * np.cos(2 * math.pi * nu[:, None] * m[None, :])).sum(axis=1)


def _pixel_sampled_otf(nu, pitch, weight, length, n_alias=8):
    """Box-filtered, aliased OTF of spec §2.1.10: sum_{k=-8..8} P^(nu + k/p) sinc(nu p + k).

    The transfer function of the pixel-integrated, pixel-sampled continuous PSF along x
    (the box's y factor vanishes at k' != 0).
    """
    k = np.arange(-n_alias, n_alias + 1)
    nu_k = nu[:, None] + k[None, :] / pitch
    box = np.sinc(nu[:, None] * pitch + k[None, :])
    return (scatter_otf(nu_k, weight, length) * box).sum(axis=1)


@pytest.mark.parametrize("length", LENGTHS)
@pytest.mark.parametrize("weight", [1.0, 0.4])
def test_scatter_kernel_normalised(length, weight):
    kernel = scatter_kernel(PITCH, weight, length)
    assert kernel.ndim == 2 and kernel.shape[0] == kernel.shape[1] and kernel.shape[0] % 2 == 1
    assert abs(kernel.sum() - 1.0) <= 1e-12
    assert np.all(kernel >= 0)
    np.testing.assert_array_equal(kernel, kernel.T)
    np.testing.assert_array_equal(kernel, kernel[::-1, ::-1])
    assert not kernel.flags.writeable


@pytest.mark.parametrize("length", LENGTHS)
@pytest.mark.parametrize("weight", [1.0, 0.4])
def test_scatter_kernel_otf(length, weight):
    # Spec §2.1.10 (revised): the kernel holds pixel integrals, so its DTFT is compared with
    # the box-filtered, aliased OTF, not with P^(nu) itself.
    nu = np.linspace(0.0, 0.25 / PITCH, 41)
    kernel = scatter_kernel(PITCH, weight, length)
    expected = _pixel_sampled_otf(nu, PITCH, weight, length)
    assert np.max(np.abs(_line_transfer(kernel, nu, PITCH) - expected)) <= 2e-3


@pytest.mark.xfail(
    strict=True,
    reason="spec §2.1.10 (revised) notes that a pixel-integrated kernel cannot match the "
    "continuous P^(nu) to 2e-3 (it differs by 0.002-0.018 w_P); kept as the record of why "
    "test_scatter_kernel_otf compares with the box-filtered, aliased OTF",
)
def test_scatter_kernel_otf_vs_continuous():
    nu = np.linspace(0.0, 0.25 / PITCH, 41)
    for length in LENGTHS:
        kernel = scatter_kernel(PITCH, 1.0, length)
        diff = np.abs(_line_transfer(kernel, nu, PITCH) - scatter_otf(nu, 1.0, length))
        assert np.max(diff) <= 2e-3


def test_scatter_otf_values():
    nu = np.array([0.0, 0.1, 1.0])
    np.testing.assert_allclose(scatter_otf(nu, 0.0, 3.0), 1.0)
    np.testing.assert_allclose(scatter_otf(nu, 1.0, 3.0), 1 / (1 + (2 * math.pi * 3.0 * nu) ** 2))
    assert scatter_otf(np.array(0.0), 0.7, 1.0) == 1.0


def test_tail_radius_matches_spec():
    # §12.7: (R/a) K1(R/a) < 1e-3, 1e-4, 1e-6 needs R = 8.24a, 10.66a, 15.44a (rounded up).
    for tail, ratio in ((1e-3, 8.24), (1e-4, 10.66), (1e-6, 15.44)):
        x = tail_radius_ratio(tail)
        assert x * special.k1(x) == pytest.approx(tail, rel=1e-9)
        assert ratio - 0.011 <= x <= ratio


def test_scatter_kernel_size_follows_tail():
    length = 1.5
    radius = tail_radius_ratio(1e-4) * length
    half = scatter_half_width(PITCH, length)
    assert (half - 0.5) * PITCH <= radius < (half + 0.5) * PITCH
    assert scatter_kernel(PITCH, 0.5, length).shape == (2 * half + 1, 2 * half + 1)
    assert scatter_kernel(PITCH, 0.5, length, tail=1e-6).shape[0] > 2 * half + 1


def test_scatter_kernel_centre_pixel_disc():
    # With a tiny tail the centre holds the equal-area disc mass, up to the renormalisation
    # that absorbs the midpoint error of the other pixels.
    length = 4.0 * PITCH
    kernel = scatter_kernel(PITCH, 1.0, length, tail=1e-9)
    half = kernel.shape[0] // 2
    rho = PITCH / math.sqrt(math.pi)
    disc = 1 - (rho / length) * special.k1(rho / length)
    assert kernel[half, half] == pytest.approx(disc, rel=1e-3)


def test_scatter_kernel_zero_weight_is_delta():
    kernel = scatter_kernel(PITCH, 0.0, 3.0)
    half = kernel.shape[0] // 2
    assert kernel[half, half] == 1.0 and kernel.sum() == 1.0


def test_scatter_kernel_delta_part():
    full = scatter_kernel(PITCH, 1.0, PITCH)
    mixed = scatter_kernel(PITCH, 0.3, PITCH)
    delta = np.zeros_like(full)
    delta[full.shape[0] // 2, full.shape[0] // 2] = 1.0
    np.testing.assert_allclose(mixed, 0.3 * full + 0.7 * delta, rtol=0, atol=1e-15)


@pytest.mark.parametrize(
    ("args", "kwargs", "match"),
    [
        ((0.0, 0.5, 1.0), {}, "pitch_um"),
        ((-1.0, 0.5, 1.0), {}, "pitch_um"),
        ((1.0, 1.5, 1.0), {}, "weight"),
        ((1.0, -0.1, 1.0), {}, "weight"),
        ((1.0, 0.5, 0.0), {}, "length_um"),
        ((1.0, 0.5, float("inf")), {}, "length_um"),
        ((1.0, 0.5, 1.0), {"tail": 0.0}, "tail"),
        ((1.0, 0.5, 1.0), {"tail": 1.0}, "tail"),
        ((1.0, 0.5, 1.0), {"supersample": 0}, "supersample"),
    ],
)
def test_scatter_kernel_rejects_bad_input(args, kwargs, match):
    with pytest.raises(ValueError, match=match):
        scatter_kernel(*args, **kwargs)


def test_scatter_otf_rejects_bad_input():
    with pytest.raises(ValueError, match="weight"):
        scatter_otf(np.array([0.1]), 2.0, 1.0)
    with pytest.raises(ValueError, match="length_um"):
        scatter_otf(np.array([0.1]), 0.5, -1.0)
