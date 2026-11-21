"""Halation fraction, kernel and transfer function (spec §2.1.5, §12.7, tests of §2.1.10)."""

import itertools
import math

import numpy as np
import pytest

from darkroom_atlas.optics.halation import (
    HalationKernel,
    _energy_radius,
    _fraction_parts,
    _irradiance,
    halation_fraction,
    halation_kernel,
    halation_otf,
)
from darkroom_atlas.optics.stack import BaseOptics

TAC_135 = BaseOptics(
    thickness_um=127.0, refractive_index=1.475, grey_density=0.10, backing_density=0
)
CLEAR = BaseOptics(thickness_um=127.0, refractive_index=1.475, grey_density=0.0, backing_density=0)
PITCH = 6.0


def _radial_integral(base: BaseOptics, r_max: float, n: int = 400_001) -> float:
    """Independent check: trapezoid of H(r) 2 pi r dr in r, split at r_c and r = 1e4 um."""
    r_c = base.critical_radius_um
    total = 0.0
    edges = [e for e in (0.0, r_c, 4.0 * r_c, 1e4) if e < r_max] + [r_max]
    for lo, hi in itertools.pairwise(edges):
        r = np.linspace(lo, hi, n)
        total += np.trapezoid(_irradiance(r, base) * 2 * math.pi * r, r)
    return total


def test_halation_fraction_clear_base():
    f_h = halation_fraction(CLEAR)
    # H ~ (2 t_b / r)^4 / (4 pi t_b^2) far out, so the tail beyond r_max is 4 t_b^2 / r_max^2.
    r_max = 2e6
    numeric = _radial_integral(CLEAR, r_max) + 4 * CLEAR.thickness_um**2 / r_max**2
    assert abs(f_h - numeric) <= 1e-4
    fresnel, total = _fraction_parts(CLEAR)
    assert abs(total - (1 - 1 / CLEAR.refractive_index**2)) <= 1e-4
    assert fresnel + total == pytest.approx(f_h, rel=1e-14)
    # §12.7: 0.581, of which 0.540 is total reflection.
    assert f_h == pytest.approx(0.581, abs=5e-4)
    assert total == pytest.approx(0.540, abs=5e-4)


def test_halation_grey_base_numbers():
    f_h = halation_fraction(TAC_135)
    assert abs(f_h - 0.222) <= 0.002
    hk = halation_kernel(PITCH, TAC_135)
    assert abs(hk.critical_radius_um - 234.0) <= PITCH
    # The steepest radial rise of the kernel (the onset of total reflection) sits at r_c.
    half = hk.kernel.shape[0] // 2
    row = hk.kernel[half, half:]
    jump = int(np.argmax(np.diff(row)))
    assert abs((jump + 0.5) * PITCH - 234.0) <= PITCH
    assert math.degrees(TAC_135.critical_angle_rad) == pytest.approx(42.7, abs=0.05)


def test_halation_grey_base_parts_and_radii():
    # §12.7: total-reflection part 0.199; 99 % and 99.9 % radii 1167 and 1831 um.
    _, total = _fraction_parts(TAC_135)
    assert total == pytest.approx(0.199, abs=5e-4)
    assert _energy_radius(TAC_135, 0.99) == pytest.approx(1167.0, rel=1e-3)
    assert _energy_radius(TAC_135, 0.999) == pytest.approx(1831.0, rel=1e-3)


def test_halation_clear_base_radii():
    # §12.7: 3324 and 10532 um (this solver gives 10536 um, 4e-4 relative).
    assert _energy_radius(CLEAR, 0.99) == pytest.approx(3324.0, rel=1e-3)
    assert _energy_radius(CLEAR, 0.999) == pytest.approx(10532.0, rel=1e-3)


def test_halation_critical_radius_120():
    base = BaseOptics(
        thickness_um=119.0, refractive_index=1.475, grey_density=0.0, backing_density=0.5
    )
    assert base.critical_radius_um == pytest.approx(220.0, abs=0.5)


def test_halation_kernel_truncation():
    hk = halation_kernel(PITCH, TAC_135)
    assert isinstance(hk, HalationKernel)
    assert hk.kernel.shape[0] == hk.kernel.shape[1] and hk.kernel.shape[0] % 2 == 1
    assert not hk.kernel.flags.writeable
    np.testing.assert_array_equal(hk.kernel, hk.kernel.T)
    np.testing.assert_array_equal(hk.kernel, hk.kernel[::-1, ::-1])
    assert hk.fraction == pytest.approx(halation_fraction(TAC_135), rel=1e-15)
    assert hk.truncated == pytest.approx(hk.fraction - hk.kernel.sum(), abs=1e-15)
    # Cut at 99.9 %: the truncated share is 1e-3 f_H up to the pixel quadrature error.
    assert hk.truncated == pytest.approx(1e-3 * hk.fraction, abs=1e-5)
    assert hk.radius_um == pytest.approx(1831.0, rel=1e-3)
    half = hk.kernel.shape[0] // 2
    assert (half - 0.5) * PITCH <= hk.radius_um < (half + 0.5) * PITCH


def test_halation_kernel_radius_cap():
    hk = halation_kernel(50.0, CLEAR)
    assert hk.radius_um == 3000.0
    # Everything returned beyond 3 mm is reported as truncated.
    beyond = hk.fraction - _radial_integral(CLEAR, 3000.0, n=40_001)
    assert hk.truncated == pytest.approx(beyond, abs=2e-4)


def test_halation_otf_zero_frequency():
    assert halation_otf(np.array([0.0]), TAC_135)[0] == pytest.approx(
        halation_fraction(TAC_135), rel=1e-10
    )


def test_halation_otf_matches_kernel():
    hk = halation_kernel(PITCH, TAC_135)
    half = hk.kernel.shape[0] // 2
    lsf = hk.kernel.sum(axis=0)
    m = np.arange(-half, half + 1) * PITCH
    nu = np.linspace(0.0, 0.25 / PITCH, 21)
    discrete = (lsf[None, :] * np.cos(2 * math.pi * nu[:, None] * m[None, :])).sum(axis=1)
    # Pixel aperture; aliases are negligible because H varies on the scale of r_c >> pitch.
    expected = halation_otf(nu, TAC_135) * np.sinc(nu * PITCH)
    assert np.max(np.abs(discrete - expected)) <= hk.truncated + 2e-5


def test_halation_otf_is_even_and_bounded():
    nu = np.array([-0.002, 0.002, 0.01, 0.05])
    values = halation_otf(nu, TAC_135)
    assert values[0] == values[1]
    assert np.all(np.abs(values) < halation_fraction(TAC_135))
    assert abs(values[3]) < abs(values[2]) < abs(values[1])


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"pitch_um": 0.0}, "pitch_um"),
        ({"pitch_um": float("nan")}, "pitch_um"),
        ({"energy_tol": 0.0}, "energy_tol"),
        ({"energy_tol": 1.0}, "energy_tol"),
        ({"max_radius_um": -1.0}, "max_radius_um"),
    ],
)
def test_halation_kernel_rejects_bad_input(kwargs, match):
    args = {"pitch_um": PITCH, "energy_tol": 1e-3, "max_radius_um": 3000.0} | kwargs
    pitch = args.pop("pitch_um")
    with pytest.raises(ValueError, match=match):
        halation_kernel(pitch, TAC_135, **args)


def test_halation_rejects_non_base():
    with pytest.raises(TypeError, match="BaseOptics"):
        halation_fraction(None)
    with pytest.raises(TypeError, match="BaseOptics"):
        halation_kernel(PITCH, "base")
    with pytest.raises(ValueError, match="finite"):
        halation_otf(np.array([np.inf]), TAC_135)
