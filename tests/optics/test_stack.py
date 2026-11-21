"""Optical stack builder (spec §2.1.2, §2.1.3, §2.1.5, §2.1.8)."""

import copy
import math
from pathlib import Path

import numpy as np
import pytest
import yaml

from darkroom_atlas.data.models import Material
from darkroom_atlas.emulsion import size_classes
from darkroom_atlas.optics.stack import (
    BaseOptics,
    base_optics,
    interface_reflectances,
    optical_stack,
)
from darkroom_atlas.spectral import DEFAULT_GRID, SpectralGrid
from darkroom_atlas.units import FARADAY

FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "library"
    / "films"
    / "test"
    / "bw-mechanistic.yaml"
)


@pytest.fixture(scope="module")
def raw() -> dict:
    return yaml.safe_load(FIXTURE.read_text())


def _material(raw: dict, edit=None) -> Material:
    data = copy.deepcopy(raw)
    if edit is not None:
        edit(data)
    return Material.model_validate(data)


@pytest.fixture(scope="module")
def material(raw) -> Material:
    return _material(raw)


@pytest.fixture(scope="module")
def classes(material):
    return size_classes(material.emulsion_layers[0])


def test_interface_reflectances_numbers():
    t_top, r_i, r_bot = interface_reflectances(1.54, 1.48)
    r_top = (0.54 / 2.54) ** 2
    assert t_top == pytest.approx(1.0 - r_top, rel=1e-15)
    assert t_top == pytest.approx(0.9548019, abs=1e-7)
    assert r_i == pytest.approx(1.0 - (1.0 - r_top) / 1.54**2, rel=1e-15)
    assert r_i == pytest.approx(0.5974018, abs=1e-7)
    assert r_bot == pytest.approx(1.0 - (1.48 / 1.54) ** 2, rel=1e-15)
    assert r_bot == pytest.approx(0.0764041, abs=1e-7)


def test_interface_no_total_reflection_into_denser_base():
    assert interface_reflectances(1.5, 1.5)[2] == 0.0
    assert interface_reflectances(1.5, 1.6)[2] == 0.0


@pytest.mark.parametrize(("n_gel", "n_b"), [(1.0, 1.48), (0.9, 1.48), (1.54, 0.0), (1.54, -1)])
def test_interface_reflectances_rejects(n_gel, n_b):
    with pytest.raises(ValueError, match="must be"):
        interface_reflectances(n_gel, n_b)


def test_stack_fixture_values(material, classes):
    stack = optical_stack(material, "135", classes, DEFAULT_GRID)
    lam = DEFAULT_GRID.wavelengths
    e_ph = 1239.84 / lam
    expected_alpha = np.where(e_ph >= 2.7, 1.0, np.exp((e_ph - 2.7) / 0.06))
    np.testing.assert_allclose(stack.alpha_per_um, expected_alpha, rtol=1e-12, atol=0)
    nu = 1000.0 / lam
    expected_dye = 0.02 * np.exp(-0.5 * ((nu - 1.8) / 0.10) ** 2) + 0.01 * np.exp(
        -0.5 * ((nu - 1.6) / 0.08) ** 2
    )
    np.testing.assert_allclose(stack.dye_cross_section, expected_dye, rtol=1e-12, atol=1e-300)
    # mu_a,b = sum_c (nu_c/L_e) sigma_c,b, sigma_c,b = V_c alpha_b + S_c s_b (spec §2.1.2).
    sigma = (
        classes.volume_um3[:, None] * expected_alpha[None, :]
        + classes.surface_um2[:, None] * expected_dye[None, :]
    )
    mu_a = np.sum(classes.number_per_um2[:, None] * sigma, axis=0) / classes.dry_thickness_um
    np.testing.assert_allclose(stack.mu_a_per_um, mu_a, rtol=1e-12, atol=0)
    np.testing.assert_array_equal(stack.km_k_per_um, 2.0 * stack.mu_a_per_um)
    np.testing.assert_allclose(stack.km_s_per_um, 0.2 * 550.0 / lam, rtol=1e-14)
    assert stack.km_s_per_um[DEFAULT_GRID.wavelengths == 550.0][0] == pytest.approx(0.2, rel=1e-15)
    assert stack.emulsion_thickness_um == classes.dry_thickness_um
    t_top, r_i, r_bot = interface_reflectances(1.54, 1.48)
    assert (stack.t_top, stack.r_internal_top, stack.r_bottom) == (t_top, r_i, r_bot)
    assert stack.psf_weight == 0.5
    assert stack.psf_length_um == 1.0
    assert stack.grid == DEFAULT_GRID


def test_urbach_edge_is_continuous(material, classes):
    # A fine grid around the edge (459.2 nm for E_0 = 2.7 eV): no jump across it.
    grid = SpectralGrid(455.0, 465.0, 0.01)
    alpha = optical_stack(material, "135", classes, grid).alpha_per_um
    assert np.max(np.abs(np.diff(alpha))) < 2e-3
    assert np.all(np.diff(alpha) <= 0.0)
    assert alpha[0] == 1.0
    assert alpha[-1] < 1.0


def test_energies_in_other_units(raw, classes, material):
    # 2.7 eV and 0.06 eV written in kJ/mol give the same edge.
    def edit(data):
        optics = data["layers"][1]["optics"]["intrinsic_absorption"]
        optics["edge_energy"] = {
            "value": 2.7 * FARADAY / 1e3,
            "unit": "kJ/mol",
            "source": "assumed",
        }
        optics["urbach_energy"] = {
            "value": 0.06 * FARADAY / 1e3,
            "unit": "kJ/mol",
            "source": "assumed",
        }

    other = optical_stack(_material(raw, edit), "135", classes, DEFAULT_GRID)
    ref = optical_stack(material, "135", classes, DEFAULT_GRID)
    np.testing.assert_allclose(other.alpha_per_um, ref.alpha_per_um, rtol=1e-12)


def test_base_optics_per_format(material):
    assert base_optics(material, "135") == BaseOptics(127.0, 1.48, 0.10, 0.0)
    assert base_optics(material, "120") == BaseOptics(119.0, 1.48, 0.0, 0.5)


def test_stack_120_uses_its_base(material, classes):
    stack = optical_stack(material, "120", classes, DEFAULT_GRID)
    assert stack.base == BaseOptics(119.0, 1.48, 0.0, 0.5)
    assert stack.base.critical_radius_um == pytest.approx(
        2 * 119.0 * math.tan(math.asin(1 / 1.48)), rel=1e-14
    )


def test_stack_rejects_unknown_format(material, classes):
    with pytest.raises(ValueError, match=r"unknown film format '220'.*135, 120"):
        optical_stack(material, "220", classes, DEFAULT_GRID)


def test_stack_rejects_missing_optics(raw, classes):
    def edit(data):
        del data["layers"][1]["optics"]

    with pytest.raises(ValueError, match="has no optics"):
        optical_stack(_material(raw, edit), "135", classes, DEFAULT_GRID)


def test_stack_rejects_missing_refractive_index(raw, classes):
    def edit(data):
        del data["layers"][2]["refractive_index"]

    with pytest.raises(ValueError, match="no refractive_index"):
        optical_stack(_material(raw, edit), "135", classes, DEFAULT_GRID)


def test_stack_rejects_bad_arguments(material, classes):
    with pytest.raises(TypeError, match="material must be a Material"):
        optical_stack("film", "135", classes, DEFAULT_GRID)
    with pytest.raises(TypeError, match="classes must be SizeClasses"):
        optical_stack(material, "135", None, DEFAULT_GRID)
    with pytest.raises(TypeError, match="grid must be a SpectralGrid"):
        optical_stack(material, "135", classes, None)


def test_stack_rejects_zero_dye_width(raw, classes):
    def edit(data):
        data["layers"][1]["optics"]["dye_bands"][0]["width"] = {
            "value": 0.0,
            "unit": "um-1",
            "source": "assumed",
        }

    with pytest.raises(ValueError, match=r"dye_bands\[0\] needs amplitude >= 0 and width > 0"):
        optical_stack(_material(raw, edit), "135", classes, DEFAULT_GRID)
