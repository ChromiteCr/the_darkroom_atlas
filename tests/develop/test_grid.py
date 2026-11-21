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
    grid = column_grid(10.0, 1.0, 3.0, 50.0, 300.0, profile)
    assert (grid.n_z, grid.n_l) == (n_z, n_l)
    assert grid.kinds == ("emulsion",) * n_z + ("overcoat",) + ("liquid",) * n_l
    assert grid.thickness_um.shape == grid.height_um.shape == (n_z + 1 + n_l,)
    assert grid.delta_max_um == 300.0


@pytest.mark.parametrize("profile", PROFILES)
def test_cell_thicknesses_sum_to_swollen_thickness(profile):
    grid = column_grid(12.0, 1.5, 3.0, 50.0, 300.0, profile)
    n_z = grid.n_z
    np.testing.assert_allclose(grid.thickness_um[:n_z], 3.0 * 12.0 / n_z, rtol=1e-15)
    assert grid.thickness_um[n_z] == pytest.approx(4.5, rel=1e-15)
    assert grid.thickness_um[: n_z + 1].sum() == pytest.approx(3.0 * 13.5, rel=1e-12)
    assert grid.gel_thickness_um == pytest.approx(40.5, rel=1e-12)
    assert grid.thickness_um[n_z + 1 :].sum() == pytest.approx(300.0, rel=1e-12)


def _sections(n_l, r_l, delta_agit, delta_max):
    def section(length, n):
        h1 = length * (r_l - 1.0) / (r_l**n - 1.0)
        return h1 * r_l ** np.arange(n)

    n_a = n_l // 2
    return np.concatenate([section(delta_agit, n_a), section(delta_max - delta_agit, n_l - n_a)])


@pytest.mark.parametrize(("profile", "n_a"), [("reference", 8), ("calibration", 4), ("image", 2)])
def test_liquid_split_matches_spec_table(profile, n_a):
    # Spec D1 3.1: liquid sections n_a + (n_l - n_a) = 8 + 8 / 4 + 4 / 2 + 2.
    grid = column_grid(10.0, 1.0, 3.0, 50.0, 300.0, profile)
    assert grid.n_a == n_a
    assert grid.n_l - grid.n_a == n_a
    assert grid.agit_cell == grid.n_z + n_a
    assert grid.delta_agit_um == 50.0


@pytest.mark.parametrize("profile", PROFILES)
@pytest.mark.parametrize("delta_agit", [20.0, 39.0, 100.0, 299.0])
def test_liquid_sections_are_geometric(profile, delta_agit):
    _, n_l, r_l = PROFILE_LEVELS[profile]
    grid = column_grid(10.0, 1.0, 3.0, delta_agit, 300.0, profile)
    liquid = grid.thickness_um[grid.n_z + 1 :]
    np.testing.assert_allclose(liquid, _sections(n_l, r_l, delta_agit, 300.0), rtol=1e-12)
    lower, upper = liquid[: grid.n_a], liquid[grid.n_a :]
    np.testing.assert_allclose(lower[1:] / lower[:-1], r_l, rtol=1e-12)
    np.testing.assert_allclose(upper[1:] / upper[:-1], r_l, rtol=1e-12)
    assert lower.sum() == pytest.approx(delta_agit, rel=1e-14)
    assert upper.sum() == pytest.approx(300.0 - delta_agit, rel=1e-14)


@pytest.mark.parametrize("profile", PROFILES)
def test_agitation_face_is_a_cell_face(profile):
    # Spec D1 2.4.1 / 2.4.8: h = delta_agit is the top face of cell agit_cell.
    grid = column_grid(10.0, 1.0, 3.0, 41.0, 300.0, profile)
    tops = np.cumsum(grid.thickness_um)
    assert tops[grid.agit_cell] - grid.gel_thickness_um == pytest.approx(41.0, rel=1e-13)
    assert tops[-1] - grid.gel_thickness_um == pytest.approx(300.0, rel=1e-13)
    assert grid.kinds[grid.agit_cell] == "liquid"
    assert grid.kinds[grid.agit_cell + 1] == "liquid"


@pytest.mark.parametrize("profile", PROFILES)
def test_heights_are_cell_centres_above_base(profile):
    grid = column_grid(10.0, 1.0, 3.0, 50.0, 300.0, profile)
    tops = np.cumsum(grid.thickness_um)
    np.testing.assert_allclose(grid.height_um, tops - 0.5 * grid.thickness_um, rtol=1e-15)
    assert np.all(np.diff(grid.height_um) > 0)
    # Liquid centres measured from the overcoat surface follow the two geometric sections.
    _, n_l, r_l = PROFILE_LEVELS[profile]
    cells = _sections(n_l, r_l, 50.0, 300.0)
    expected = np.cumsum(cells) - 0.5 * cells
    above_surface = grid.height_um[grid.n_z + 1 :] - grid.gel_thickness_um
    np.testing.assert_allclose(above_surface, expected, rtol=1e-12)


def test_arrays_are_read_only():
    grid = column_grid(10.0, 1.0, 3.0, 50.0, 300.0, "image")
    with pytest.raises(ValueError):
        grid.thickness_um[0] = 1.0
    with pytest.raises(ValueError):
        grid.height_um[0] = 1.0


@pytest.mark.parametrize(
    "args",
    [
        (0.0, 1.0, 3.0, 50.0, 300.0, "image"),
        (-1.0, 1.0, 3.0, 50.0, 300.0, "image"),
        (10.0, 0.0, 3.0, 50.0, 300.0, "image"),
        (10.0, 1.0, 0.5, 50.0, 300.0, "image"),
        (10.0, 1.0, 3.0, 50.0, 0.0, "image"),
        (10.0, 1.0, 3.0, 50.0, float("inf"), "image"),
        (10.0, 1.0, float("nan"), 50.0, 300.0, "image"),
        (10.0, 1.0, 3.0, 50.0, 300.0, "fine"),
        (10.0, 1.0, 3.0, 0.0, 300.0, "image"),
        (10.0, 1.0, 3.0, -5.0, 300.0, "image"),
        (10.0, 1.0, 3.0, float("nan"), 300.0, "image"),
        (10.0, 1.0, 3.0, 300.0, 300.0, "image"),
        (10.0, 1.0, 3.0, 400.0, 300.0, "image"),
    ],
)
def test_bad_input_raises(args):
    with pytest.raises(ValueError):
        column_grid(*args)


def test_non_number_raises_type_error():
    with pytest.raises(TypeError):
        column_grid("10", 1.0, 3.0, 50.0, 300.0, "image")
