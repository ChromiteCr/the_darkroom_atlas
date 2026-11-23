"""Tank balance, agitation mixing and prepass norms (spec D1 2.4.8, 2.4.9, 3.7).

These check the tank kernels. The chain-level ``test_tank_mass_balance`` and
``test_tank_infinite_volume`` of spec D1 2.4.13 (film + liquid + tank inventories of a whole
strip) belong with the ensemble integrator.
"""

import math

import numpy as np
import pytest

from darkroom_atlas.develop.tank import (
    CONTRACTION_WARN,
    DEPLETION_FLOOR_MOL_L,
    MAX_EXTRA_PASSES,
    SILVER_FLOOR_G_M2,
    TANK_TOLERANCE,
    UM3_TO_L,
    advance_bulk,
    area_mean,
    bulk_at,
    combine_flux,
    contraction,
    depletion_change,
    depletion_scale,
    integrate_bulk,
    mixing_flux,
    mixing_inventory,
    prepass_delta,
    prepass_passed,
    tank_coefficient,
)

# One 135-36 roll, 1.5 x 0.035 m, in a 0.25 L small tank.
AREA_UM2 = 1.5e6 * 3.5e4
VOLUME_L = 0.25


def _schedule(rng, n=200):
    h = rng.uniform(0.05, 5.0, size=n)
    return np.concatenate([[0.0], np.cumsum(h)])


def test_constants_match_spec():
    assert UM3_TO_L == 1e-15
    assert TANK_TOLERANCE == 1e-3
    assert SILVER_FLOOR_G_M2 == 0.3
    assert DEPLETION_FLOOR_MOL_L == 1e-6
    assert CONTRACTION_WARN == 0.5
    assert MAX_EXTRA_PASSES == 2


def test_tank_coefficient():
    assert tank_coefficient(film_area_um2=AREA_UM2, volume_l=VOLUME_L) == pytest.approx(
        1e-15 * AREA_UM2 / VOLUME_L
    )


def test_tank_mass_balance_kernel():
    # Spec D1 2.4.13 test_tank_mass_balance at the kernel level: the tank change equals
    # -integral(Jbar) A / V to 1e-12 of the tank inventory, and tank + film moles are conserved.
    rng = np.random.default_rng(1)
    t = _schedule(rng)
    flux = rng.normal(scale=0.05, size=(t.size - 1, 6))  # mol L^-1 um s^-1 (D-76 order)
    c0 = np.array([0.79, 0.021, 0.045, 0.0, 0.0, 0.08])
    coef = tank_coefficient(film_area_um2=AREA_UM2, volume_l=VOLUME_L)
    traj = integrate_bulk(c0, flux, t, coef)
    integral = (np.diff(t)[:, None] * flux).sum(axis=0)
    expect = -coef * integral
    change = traj[-1] - c0
    assert np.max(np.abs(change - expect)) <= 1e-12 * np.max(np.abs(c0))
    # Moles: V c (tank, mol) + 1e-15 A * integral (into the film, mol) is constant.
    moles = VOLUME_L * traj[-1] + UM3_TO_L * AREA_UM2 * integral
    np.testing.assert_allclose(moles, VOLUME_L * c0, rtol=0, atol=1e-12 * VOLUME_L * c0.max())


def test_tank_infinite_volume_kernel():
    # Spec D1 2.4.13 test_tank_infinite_volume at the kernel level: V = 1e9 L.
    rng = np.random.default_rng(2)
    t = _schedule(rng)
    flux = rng.normal(scale=0.05, size=(t.size - 1, 3))
    c0 = np.array([0.79, 0.0, 0.08])
    traj = integrate_bulk(c0, flux, t, tank_coefficient(film_area_um2=AREA_UM2, volume_l=1e9))
    assert np.max(np.abs(traj - c0)) <= 1e-12


def test_integrate_bulk_reproduces_stepwise_updates():
    # Re-integrating recorded fluxes (spec D1 3.7 step 4) is bit-identical to the coupled run.
    rng = np.random.default_rng(3)
    t = _schedule(rng, 50)
    flux = rng.normal(size=(t.size - 1, 4))
    c0 = rng.uniform(size=4)
    coef = tank_coefficient(film_area_um2=AREA_UM2, volume_l=VOLUME_L)
    c = c0
    stepwise = [c0]
    for n in range(t.size - 1):
        c = advance_bulk(c, flux[n], t[n + 1] - t[n], coef)
        stepwise.append(c)
    assert np.array_equal(integrate_bulk(c0, flux, t, coef), np.array(stepwise))


def test_mixing_impulse_conserves_and_matches_jump():
    # Spec D1 2.4.8/2.4.9: the upper-section excess joins the tank; booked as -dI/h in the
    # step's flux it gives the jump +coef dI, and film + tank moles are unchanged.
    thick = np.array([20.0, 35.0, 61.0, 106.0])  # um, upper section
    bulk = np.array([0.79, 0.002, 0.08])
    rng = np.random.default_rng(4)
    conc = bulk + rng.normal(scale=1e-3, size=(5, 4, 3))  # 5 columns
    weights = np.array([1.0, 2.0, 1.0, 3.0, 1.0])
    inv = mixing_inventory(conc, thick, bulk)
    expect = ((conc - bulk) * thick[:, None]).sum(axis=-2)
    np.testing.assert_allclose(inv, expect, rtol=1e-12, atol=1e-18)
    inv_mean = area_mean(inv, weights)
    coef = tank_coefficient(film_area_um2=AREA_UM2, volume_l=VOLUME_L)
    h = 0.1
    after = advance_bulk(bulk, mixing_flux(inv_mean, h), h, coef)
    np.testing.assert_allclose(after, bulk + coef * inv_mean, rtol=1e-15)
    # Moles: before, the film's upper section holds sum A_m Delta (C - C_bulk) extra; after,
    # the section is at bath values and the tank holds it.
    film_area = AREA_UM2 * weights / weights.sum()
    extra = UM3_TO_L * (film_area[:, None] * inv).sum(axis=0)
    np.testing.assert_allclose(VOLUME_L * after, VOLUME_L * bulk + extra, rtol=1e-15)
    # A section at bath values has no excess.
    assert np.array_equal(mixing_inventory(np.broadcast_to(bulk, (4, 3)), thick, bulk), 0 * bulk)


def test_area_mean_and_combine():
    v = np.array([[1.0, 2.0], [3.0, 6.0]])
    np.testing.assert_allclose(area_mean(v, np.array([1.0, 3.0])), [2.5, 5.0])
    jb, jz = np.array([1.0, -2.0]), np.array([0.5, 0.0])
    np.testing.assert_allclose(combine_flux(jb, jz, 0.65), 0.65 * jb + 0.35 * jz)
    np.testing.assert_array_equal(combine_flux(jb, jz, 1.0), jb)


def test_bulk_at_is_linear_and_exact_at_nodes():
    t = np.array([0.0, 0.1, 0.26, 1.0])
    v = np.array([[1.0, 0.0], [0.9, 0.01], [0.85, 0.02], [0.8, 0.05]])
    for n in range(t.size):
        assert np.array_equal(bulk_at(t, v, t[n]), v[n])
    np.testing.assert_allclose(bulk_at(t, v, 0.18), 0.5 * (v[1] + v[2]), rtol=1e-14)


def test_depletion_norm_handles_zero_initial_values():
    # Spec D1 3.7 (revised): bromide, iodide and silver start at 0; s_k uses the change.
    cur = np.array([[0.79, 0.0, 0.0], [0.78, 2e-3, 0.0], [0.77, 4e-3, 0.0]])
    new = cur + np.array([[0.0, 0.0, 0.0], [1e-4, 1e-5, 0.0], [2e-4, 4e-5, 5e-9]])
    s = depletion_scale(cur)
    np.testing.assert_allclose(s, [0.79, 4e-3, 1e-6])
    assert depletion_change(new, cur) == pytest.approx(max(2e-4 / 0.79, 4e-5 / 4e-3, 5e-9 / 1e-6))
    assert depletion_change(cur, cur) == 0.0


def test_prepass_delta_and_pass():
    cur = np.array([0.0, 0.05, 1.2, 2.5])
    new = cur + np.array([2e-4, -1e-4, 1.4e-3, -2e-3])
    delta = prepass_delta(new, cur)
    assert delta == pytest.approx(max(2e-4 / 0.3, 1e-4 / 0.35, 1.4e-3 / 1.5, 2e-3 / 2.8))
    assert prepass_passed(delta)
    assert prepass_passed(TANK_TOLERANCE)
    assert not prepass_passed(1.0001e-3)
    assert prepass_delta(cur, cur) == 0.0


def test_contraction():
    assert contraction([4e-3, 1e-3, 2e-4]) == pytest.approx((0.25, 0.2))
    assert contraction([1e-3]) == ()
    assert contraction([0.0, 0.0]) == (0.0,)
    assert contraction([0.0, 1e-4]) == (math.inf,)


# ---- bad input -------------------------------------------------------------------------------


def test_tank_rejects_bad_input():
    with pytest.raises(ValueError, match="volume"):
        tank_coefficient(film_area_um2=AREA_UM2, volume_l=0.0)
    with pytest.raises(ValueError, match="film_area"):
        tank_coefficient(film_area_um2=-1.0, volume_l=1.0)
    with pytest.raises(TypeError):
        tank_coefficient(AREA_UM2, VOLUME_L)  # positional order would mirror TankSpec's
    with pytest.raises(ValueError, match="h_s"):
        advance_bulk(np.ones(2), np.ones(2), 0.0, 1.0)
    with pytest.raises(ValueError, match="shape"):
        advance_bulk(np.ones(2), np.ones(3), 1.0, 1.0)
    with pytest.raises(ValueError, match="coefficient"):
        advance_bulk(np.ones(2), np.ones(2), 1.0, -1.0)
    with pytest.raises(ValueError, match="finite"):
        advance_bulk(np.array([1.0, np.nan]), np.ones(2), 1.0, 1.0)
    with pytest.raises(ValueError, match="flux_mean"):
        integrate_bulk(np.ones(2), np.ones((3, 2)), np.arange(3.0), 1.0)
    with pytest.raises(ValueError, match="increase"):
        integrate_bulk(np.ones(2), np.ones((2, 2)), np.array([0.0, 1.0, 1.0]), 1.0)
    with pytest.raises(ValueError, match="conc"):
        mixing_inventory(np.ones((3, 2)), np.ones(2), np.ones(2))
    with pytest.raises(ValueError, match="thickness"):
        mixing_inventory(np.ones((2, 2)), np.array([1.0, 0.0]), np.ones(2))
    with pytest.raises(ValueError, match="weights"):
        area_mean(np.ones((2, 1)), np.zeros(2))
    with pytest.raises(ValueError, match="rows"):
        area_mean(np.ones((3, 1)), np.ones(2))
    with pytest.raises(ValueError, match="fraction"):
        combine_flux(np.ones(2), np.ones(2), 1.5)
    with pytest.raises(ValueError, match="outside"):
        bulk_at(np.array([0.0, 1.0]), np.ones((2, 1)), 1.5)
    with pytest.raises(ValueError, match="shape"):
        depletion_change(np.ones((2, 2)), np.ones((3, 2)))
    with pytest.raises(ValueError, match=">= 0"):
        prepass_delta(np.array([-1.0]), np.array([1.0]))
    with pytest.raises(ValueError):
        prepass_passed(float("nan"))
    with pytest.raises(ValueError):
        contraction([1.0, -1.0])
