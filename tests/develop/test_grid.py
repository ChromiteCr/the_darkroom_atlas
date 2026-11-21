import numpy as np
import pytest

from darkroom_atlas.develop.grid import PROFILE_LEVELS, column_grid

PROFILES = ("reference", "calibration", "image")


def test_profile_levels_match_spec_table():
    # Spec D1 3.1.
    assert PROFILE_LEVELS == {
        "reference": (8, 16, 1.25),
        "calibration": (4, 8, 1.4),
        "image": (2, 4, 2.0),
    }


@pytest.mark.parametrize("profile", PROFILES)
def test_cell_counts_and_kinds(profile):
    n_z, n_l, _ = PROFILE_LEVELS[profile]
    grid = column_grid(10.0, 1.0, 3.0, 300.0, profile)
    assert (grid.n_z, grid.n_l) == (n_z, n_l)
    assert grid.kinds == ("emulsion",) * n_z + ("overcoat",) + ("liquid",) * n_l
    assert grid.thickness_um.shape == grid.height_um.shape == (n_z + 1 + n_l,)
    assert grid.delta_max_um == 300.0


@pytest.mark.parametrize("profile", PROFILES)
def test_cell_thicknesses_sum_to_swollen_thickness(profile):
    grid = column_grid(12.0, 1.5, 3.0, 300.0, profile)
    n_z = grid.n_z
    np.testing.assert_allclose(grid.thickness_um[:n_z], 3.0 * 12.0 / n_z, rtol=1e-15)
    assert grid.thickness_um[n_z] == pytest.approx(4.5, rel=1e-15)
    assert grid.thickness_um[: n_z + 1].sum() == pytest.approx(3.0 * 13.5, rel=1e-12)
    assert grid.gel_thickness_um == pytest.approx(40.5, rel=1e-12)
    assert grid.thickness_um[n_z + 1 :].sum() == pytest.approx(300.0, rel=1e-12)


@pytest.mark.parametrize("profile", PROFILES)
def test_liquid_cells_are_geometric(profile):
    _, n_l, r_l = PROFILE_LEVELS[profile]
    grid = column_grid(10.0, 1.0, 3.0, 300.0, profile)
    liquid = grid.thickness_um[grid.n_z + 1 :]
    h1 = 300.0 * (r_l - 1.0) / (r_l**n_l - 1.0)
    np.testing.assert_allclose(liquid, h1 * r_l ** np.arange(n_l), rtol=1e-12)


@pytest.mark.parametrize("profile", PROFILES)
def test_heights_are_cell_centres_above_base(profile):
    grid = column_grid(10.0, 1.0, 3.0, 300.0, profile)
    tops = np.cumsum(grid.thickness_um)
    np.testing.assert_allclose(grid.height_um, tops - 0.5 * grid.thickness_um, rtol=1e-15)
    assert np.all(np.diff(grid.height_um) > 0)
    # Liquid centres measured from the overcoat surface follow the geometric series.
    _, n_l, r_l = PROFILE_LEVELS[profile]
    h1 = 300.0 * (r_l - 1.0) / (r_l**n_l - 1.0)
    cells = h1 * r_l ** np.arange(n_l)
    expected = np.cumsum(cells) - 0.5 * cells
    above_surface = grid.height_um[grid.n_z + 1 :] - grid.gel_thickness_um
    np.testing.assert_allclose(above_surface, expected, rtol=1e-12)


def test_arrays_are_read_only():
    grid = column_grid(10.0, 1.0, 3.0, 300.0, "image")
    with pytest.raises(ValueError):
        grid.thickness_um[0] = 1.0
    with pytest.raises(ValueError):
        grid.height_um[0] = 1.0


@pytest.mark.parametrize(
    "args",
    [
        (0.0, 1.0, 3.0, 300.0, "image"),
        (-1.0, 1.0, 3.0, 300.0, "image"),
        (10.0, 0.0, 3.0, 300.0, "image"),
        (10.0, 1.0, 0.5, 300.0, "image"),
        (10.0, 1.0, 3.0, 0.0, "image"),
        (10.0, 1.0, 3.0, float("inf"), "image"),
        (10.0, 1.0, float("nan"), 300.0, "image"),
        (10.0, 1.0, 3.0, 300.0, "fine"),
    ],
)
def test_bad_input_raises(args):
    with pytest.raises(ValueError):
        column_grid(*args)


def test_non_number_raises_type_error():
    with pytest.raises(TypeError):
        column_grid("10", 1.0, 3.0, 300.0, "image")
