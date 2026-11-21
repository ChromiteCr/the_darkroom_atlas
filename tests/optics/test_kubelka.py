"""Two-stream depth response (spec §2.1.3, tests of §2.1.10)."""

import numpy as np
import pytest

from darkroom_atlas.optics.kubelka import DepthResponse, _two_stream, depth_response
from darkroom_atlas.optics.stack import BaseOptics, OpticalStack
from darkroom_atlas.spectral import SpectralGrid

BASE = BaseOptics(thickness_um=127.0, refractive_index=1.475, grey_density=0.1, backing_density=0)


def _grid(n_bands: int) -> SpectralGrid:
    return SpectralGrid(380.0, 380.0 + 5.0 * (n_bands - 1), 5.0)


def _stack(k, s, thickness=5.0, t_top=0.955, r_i=0.6, r_bot=0.08) -> OpticalStack:
    k = np.atleast_1d(np.asarray(k, dtype=float))
    if k.size == 1:  # a SpectralGrid has at least two bands
        k = np.repeat(k, 2)
    s = np.broadcast_to(np.asarray(s, dtype=float), k.shape)
    return OpticalStack(
        grid=_grid(k.size),
        emulsion_thickness_um=thickness,
        alpha_per_um=np.full(k.shape, 0.5),
        dye_cross_section=np.full(k.shape, 0.01),
        mu_a_per_um=k / 2.0,
        km_k_per_um=k,
        km_s_per_um=s,
        t_top=t_top,
        r_internal_top=r_i,
        r_bottom=r_bot,
        base=BASE,
        psf_weight=0.5,
        psf_length_um=1.0,
    )


def _solve(k, s, length, t_top, r_i, r_bot, n_in, j_in, n_z=8, branch="auto"):
    return _two_stream(
        k=np.atleast_1d(np.asarray(k, dtype=float)),
        s=np.atleast_1d(np.asarray(s, dtype=float)),
        thickness_um=length,
        t_top=t_top,
        r_top=r_i,
        r_bottom=r_bot,
        n_in=n_in,
        j_in=j_in,
        n_z=n_z,
        branch=branch,
    )


def _balance_error(sol, t_top, r_i, r_bot, n_in, j_in):
    source = t_top * n_in + j_in
    out = sol["absorbed"] + (1 - r_i) * sol["up_top"] + (1 - r_bot) * sol["down_bottom"]
    return np.abs(out - source) / source


def test_km_energy_balance():
    rng = np.random.default_rng(20261121)
    worst = 0.0
    for _ in range(100):
        k = 10.0 ** rng.uniform(-4.0, 0.5)
        s = 10.0 ** rng.uniform(-4.0, 1.0)
        length = rng.uniform(0.5, 30.0)
        r_i = rng.uniform(0.0, 0.95)
        r_bot = rng.uniform(0.0, 0.95)
        t_top = rng.uniform(0.5, 1.0)
        for j_in in (0.0, 1.0):
            sol = _solve(k, s, length, t_top, r_i, r_bot, 1.0, j_in)
            worst = max(worst, float(_balance_error(sol, t_top, r_i, r_bot, 1.0, j_in)[0]))
    assert worst <= 1e-12


@pytest.mark.parametrize("k", [0.0, 1e-30, 1e-20, 1e-12, 1e-6])
def test_km_energy_balance_vanishing_absorption(k):
    # K/S -> 0 (red bands with an Urbach tail far below the edge): the hyperbolic branch.
    for j_in in (0.0, 1.0):
        sol = _solve(k, 0.2, 5.0, 0.955, 0.6, 0.08, 1.0, j_in)
        assert np.all(np.isfinite(sol["cell_mean_u"]))
        assert _balance_error(sol, 0.955, 0.6, 0.08, 1.0, j_in)[0] <= 1e-12
        assert 0.0 <= sol["absorbed"][0] <= 10.0 * k * 5.0


def test_km_branches_agree_near_switch():
    rng = np.random.default_rng(7)
    for _ in range(50):
        length = rng.uniform(1.0, 20.0)
        beta_l = rng.uniform(0.3, 3.0)
        s = 10.0 ** rng.uniform(-3.0, 1.0)
        beta = beta_l / length
        k = -s + np.sqrt(s * s + beta * beta)  # solves K (K + 2S) = beta^2
        args = (k, s, length, 0.95, rng.uniform(0, 0.9), rng.uniform(0, 0.9), 1.0, 0.7)
        hyp = _solve(*args, branch="hyperbolic")
        exp = _solve(*args, branch="exponential")
        for key in ("cell_mean_u", "up_top", "down_bottom", "absorbed"):
            np.testing.assert_allclose(hyp[key], exp[key], rtol=1e-11, atol=1e-14)


def test_km_beer_lambert_limit():
    k = np.array([1e-3, 0.01, 0.1, 0.5, 1.0, 3.0])
    length = 5.0
    resp = depth_response(_stack(k, 1e-12, thickness=length, r_i=0.0, r_bot=0.0), 4)
    t_top = 0.955
    np.testing.assert_allclose(resp.absorbed / t_top, 1 - np.exp(-k * length), rtol=0, atol=1e-10)


def test_km_thick_layer_reflectance():
    k, s = 0.01, 1.0
    beta = np.sqrt(k * (k + 2 * s))
    r_inf = s / (k + s + beta)
    # Escaped / entering with no top interface.
    resp = depth_response(_stack(k, s, thickness=1e4, t_top=1.0, r_i=0.0, r_bot=0.3), 16)
    assert abs(resp.escaped_top[0] - r_inf) <= 1e-9
    # The internal ratio I-(0)/I+(0) tends to R_inf whatever the top interface.
    sol = _solve(k, s, 1e4, 0.955, 0.6, 0.3, 1.0, 0.0)
    assert abs(sol["up_top"][0] / sol["down_top"][0] - r_inf) <= 1e-9


def test_depth_response_balance_and_cells():
    k = np.array([0.0, 1e-25, 1e-4, 0.03, 0.3, 2.0])
    s = np.array([0.2, 0.2, 0.05, 0.3, 0.1, 0.0])
    stack = _stack(k, s, thickness=6.0)
    resp = depth_response(stack, 7)
    assert resp.plus.shape == (7, 6) and resp.minus.shape == (7, 6) and resp.n_z == 7
    total = resp.absorbed + resp.escaped_top + resp.into_base
    np.testing.assert_allclose(total, stack.t_top, rtol=1e-12)
    # mu_a * psi = K (I+ + I-): absorption from the cell means of psi.
    dz = stack.emulsion_thickness_um / 7
    absorbed = (stack.mu_a_per_um * resp.plus).sum(axis=0) * dz
    np.testing.assert_allclose(absorbed, resp.absorbed, rtol=1e-12, atol=1e-300)
    assert np.all(resp.plus > 0) and np.all(resp.minus > 0)
    assert np.all(resp.into_base >= 0) and np.all(resp.escaped_top >= 0)


def test_depth_response_cells_are_additive():
    stack = _stack([1e-20, 0.02, 0.4, 1.5], [0.3, 0.1, 0.2, 0.05])
    one = depth_response(stack, 3)
    two = depth_response(stack, 6)
    np.testing.assert_allclose(one.plus, two.plus.reshape(3, 2, -1).mean(axis=1), rtol=1e-12)
    np.testing.assert_allclose(one.minus, two.minus.reshape(3, 2, -1).mean(axis=1), rtol=1e-12)
    np.testing.assert_array_equal(one.into_base, two.into_base)


def test_depth_response_light_from_below_peaks_at_bottom():
    resp = depth_response(_stack([0.5], [0.1]), 5)
    assert np.all(np.diff(resp.plus[:, 0]) < 0)
    assert np.all(np.diff(resp.minus[:, 0]) > 0)


def test_depth_response_is_read_only():
    resp = depth_response(_stack([0.1, 0.2], [0.1, 0.1]), 2)
    assert isinstance(resp, DepthResponse)
    for name in ("plus", "minus", "into_base", "escaped_top", "absorbed"):
        with pytest.raises(ValueError):
            getattr(resp, name)[0] = 1.0


@pytest.mark.parametrize("n_z", [0, -1, 2.5, True, "3"])
def test_depth_response_rejects_bad_n_z(n_z):
    with pytest.raises(ValueError, match="n_z"):
        depth_response(_stack([0.1], [0.1]), n_z)


def test_depth_response_rejects_non_stack():
    with pytest.raises(TypeError, match="OpticalStack"):
        depth_response(object(), 4)


def test_optical_stack_validation():
    good = _stack([0.1, 0.2], [0.1, 0.1])
    assert not good.km_k_per_um.flags.writeable
    with pytest.raises(ValueError, match="2 \\* mu_a"):
        _stack_with(good, km_k_per_um=np.array([0.1, 0.3]))
    with pytest.raises(ValueError, match="shape"):
        _stack_with(good, km_s_per_um=np.array([0.1, 0.1, 0.1]))
    with pytest.raises(ValueError, match="non-negative"):
        _stack_with(good, km_s_per_um=np.array([0.1, -0.1]))
    with pytest.raises(ValueError, match="finite"):
        _stack_with(good, alpha_per_um=np.array([np.nan, 0.1]))
    with pytest.raises(ValueError, match="t_top"):
        _stack_with(good, t_top=0.0)
    with pytest.raises(ValueError, match="r_bottom"):
        _stack_with(good, r_bottom=1.0)
    with pytest.raises(ValueError, match="r_internal_top"):
        _stack_with(good, r_internal_top=-0.1)
    with pytest.raises(ValueError, match="emulsion_thickness_um"):
        _stack_with(good, emulsion_thickness_um=0.0)
    with pytest.raises(ValueError, match="psf_weight"):
        _stack_with(good, psf_weight=1.5)
    with pytest.raises(ValueError, match="psf_length_um"):
        _stack_with(good, psf_length_um=0.0)
    with pytest.raises(TypeError, match="BaseOptics"):
        _stack_with(good, base=None)


def _stack_with(stack: OpticalStack, **changes) -> OpticalStack:
    fields = {name: getattr(stack, name) for name in stack.__dataclass_fields__}
    fields.update(changes)
    return OpticalStack(**fields)


def test_base_optics_validation():
    assert BASE.critical_radius_um == pytest.approx(234.26, abs=0.01)
    with pytest.raises(ValueError, match="thickness_um"):
        BaseOptics(0.0, 1.475, 0.1, 0.0)
    with pytest.raises(ValueError, match="refractive_index"):
        BaseOptics(127.0, 1.0, 0.1, 0.0)
    with pytest.raises(ValueError, match="grey_density"):
        BaseOptics(127.0, 1.475, -0.1, 0.0)
    with pytest.raises(ValueError, match="backing_density"):
        BaseOptics(127.0, 1.475, 0.1, float("nan"))
