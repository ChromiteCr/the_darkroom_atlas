import numpy as np
import pytest

from darkroom_atlas.backend import xp_of
from darkroom_atlas.units import (
    C_LIGHT,
    H_PLANCK,
    K_M,
    N_AVOGADRO,
    R_GAS,
    convert,
    density_to_transmittance,
    parse_dilution,
    parse_duration_s,
    photon_energy_j,
    pixel_pitch_um,
    transmittance_to_density,
)


def test_constants():
    assert H_PLANCK == 6.62607015e-34
    assert C_LIGHT == 299792458.0
    assert N_AVOGADRO == 6.02214076e23
    assert R_GAS == 8.314462618
    assert K_M == 683.0


def test_photon_energy_555nm():
    assert photon_energy_j(555.0) == pytest.approx(3.5792e-19, rel=1e-4)


def test_photon_energy_across_the_grid():
    # The ideal round trip cannot catch a wrong photon energy (camera and model share it),
    # so anchor it at every band against h*c = 1239.841984 eV nm and 1 eV = 1.602176634e-19 J.
    wavelengths = np.arange(380.0, 781.0, 5.0)
    expected = 1239.841984 / wavelengths * 1.602176634e-19
    np.testing.assert_allclose(photon_energy_j(wavelengths), expected, rtol=1e-9, atol=0.0)


def test_density_roundtrip():
    d = np.array([0.0, 0.3, 1.0, 3.0])
    back = transmittance_to_density(density_to_transmittance(d))
    np.testing.assert_allclose(back, d, rtol=0, atol=1e-12)


def test_density_floor():
    assert transmittance_to_density(0.0) == pytest.approx(12.0)


def test_density_array_api_strict():
    import array_api_strict as xps

    def same_namespace(arr):
        ns = xp_of(arr)
        return ns is xps or "strict" in ns.__name__

    t = xps.asarray([1.0, 0.5, 0.0])
    d = transmittance_to_density(t)
    assert same_namespace(d)
    assert np.allclose(np.asarray(d), [0.0, -np.log10(0.5), 12.0], rtol=0, atol=1e-12)

    t_back = density_to_transmittance(xps.asarray([0.0, 1.0]))
    assert same_namespace(t_back)
    assert np.allclose(np.asarray(t_back), [1.0, 0.1], rtol=0, atol=1e-15)

    e = photon_energy_j(xps.asarray([555.0]))
    assert same_namespace(e)
    assert float(np.asarray(e)[0]) == pytest.approx(3.5792e-19, rel=1e-4)


def test_convert():
    assert convert(1.0, "mm", "um") == 1000.0
    assert convert(20.0, "C", "K") == 293.15
    assert convert(2.0, "min", "s") == 120.0
    assert convert(5, "g/L", "g/L") == 5


def test_convert_rejects():
    with pytest.raises(ValueError):
        convert(1, "mm", "s")
    with pytest.raises(ValueError):
        convert(1, "furlong", "m")


def test_parse_duration():
    assert parse_duration_s("9:45") == 585.0
    assert parse_duration_s("1:02:03") == 3723.0
    assert parse_duration_s("30s") == 30.0
    assert parse_duration_s("8min") == 480.0
    assert parse_duration_s(90) == 90.0
    assert parse_duration_s("1.5min") == 90.0


@pytest.mark.parametrize("text", ["", "abc", "-5", "9:75"])
def test_parse_duration_rejects(text):
    with pytest.raises(ValueError):
        parse_duration_s(text)


def test_parse_dilution():
    assert parse_dilution("stock") == 1.0
    assert parse_dilution("1+1") == 0.5
    assert parse_dilution("1+50") == 1 / 51
    assert parse_dilution("1+0") == 1.0


@pytest.mark.parametrize("text", ["1:1", "0+1", "x"])
def test_parse_dilution_rejects(text):
    with pytest.raises(ValueError):
        parse_dilution(text)


def test_pixel_pitch():
    assert pixel_pitch_um(36.0, (4000, 6000)) == 6.0
    assert pixel_pitch_um(36.0, (6000, 4000)) == 6.0
