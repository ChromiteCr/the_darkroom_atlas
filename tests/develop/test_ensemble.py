"""Development engine: film column model and ensemble integrator (spec D1 2.4, 3.1-3.5, 3.7,
2.4.13, 3.9). Developed silver is asserted in g/m2 (spec D1 3.8, change num11)."""

import dataclasses
import importlib.util
import itertools
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.special import erf

from darkroom_atlas.chemistry.kinetics import (
    arrhenius_factor,
    cluster_hazards,
    positive_part,
    resolve_kinetics,
    surface_state,
)
from darkroom_atlas.chemistry.solution import working_solution
from darkroom_atlas.chemistry.speciation import mean_protons, water_ions
from darkroom_atlas.data.models import Agitation, BathComponent
from darkroom_atlas.data.provenance import Quantity
from darkroom_atlas.develop import ensemble as ens
from darkroom_atlas.develop import rosenbrock as ros
from darkroom_atlas.develop import tank as tk
from darkroom_atlas.develop.ensemble import (
    BulkTrajectory,
    Lateral,
    Members,
    TankSpec,
    film_column_model,
    run_ensemble,
)
from darkroom_atlas.develop.grid import column_grid
from darkroom_atlas.develop.schedule import Schedule, agitation_periods, step_schedule
from darkroom_atlas.emulsion.classes import M_AG


def _load_helpers():
    # The helpers live in this directory's conftest.py; load it by path under its own name
    # (another tests/*/conftest.py would shadow a plain "import conftest").
    name = "d1c_develop_conftest"
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name("conftest.py"))
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


_h = _load_helpers()
EXPOSURE_TIME_S = _h.EXPOSURE_TIME_S
INFINITE_TANK = _h.INFINITE_TANK
LARGE_TANK = _h.LARGE_TANK
ROLL_TANK = _h.ROLL_TANK
SMALL_TANK = _h.SMALL_TANK
columns = _h.columns
developed_silver = _h.developed_silver
make_schedule = _h.make_schedule
make_solution = _h.make_solution
run_strip = _h.run_strip
strip_exposures = _h.strip_exposures
with_solvent = _h.with_solvent
with_test_rates = _h.with_test_rates
with_transport = _h.with_transport

CONTINUOUS = Agitation(method="continuous")
G_PER_MOL_UM = M_AG * 1e-3  # 1 mol L^-1 um of Ag = 0.10787 g/m2


def _integrate(film, sol, kin, schedule, ev, es, *, tank=ROLL_TANK, **kw):
    n = ev.shape[0]
    return ens._integrate(
        film,
        sol,
        kin,
        schedule,
        columns(ev, es),
        n_cells=n,
        cell_weights=np.ones(n),
        tank=tank,
        **kw,
    )


def _constant_bulk(sol, schedule):
    s = schedule.times_s.size
    return BulkTrajectory(
        times_s=schedule.times_s,
        totals=np.tile(np.asarray(sol.totals), (s, 1)),
        toth=np.full(s, sol.toth),
    )


@pytest.fixture(scope="module")
def strip_run(film_cal, d76, kin_test):
    """12 steps, small tank, 7.5 min, calibration profile, a one-roll tank."""
    schedule = make_schedule("calibration", 450.0, SMALL_TANK, record_s=(120.0, 300.0))
    ev, es = strip_exposures(12, film_cal.grid.n_z)
    return _integrate(film_cal, d76, kin_test, schedule, ev, es, record_cells=16)


# ---------------------------------------------------------------------------------------
# Film column model
# ---------------------------------------------------------------------------------------


def test_film_column_model_arrays(film_cal, film_material, kin_d76, d76):
    classes = film_cal.classes
    grid = film_cal.grid
    assert classes.count == 6
    assert film_cal.latent.log10_m.size == 201
    ref = column_grid(classes.dry_thickness_um, 1.0, 3.0, 50.0, 300.0, "calibration")
    np.testing.assert_array_equal(grid.thickness_um, ref.thickness_um)
    # Spec D1 2.0.2: sum_c C_c over the column is the coated silver.
    coat = np.sum(film_cal.silver_per_cell) * grid.n_z * grid.thickness_um[0] * G_PER_MOL_UM
    assert coat == pytest.approx(classes.silver_coverage_g_m2, rel=1e-12)
    expected = (
        1e15
        * classes.number_per_um2
        * classes.silver_g_per_grain
        / (M_AG * grid.n_z * grid.thickness_um[0])
    )
    np.testing.assert_allclose(film_cal.silver_per_cell, expected, rtol=1e-14)
    # Spec D1 2.4.7: T_gel = 1e3 m_gel / (alpha_sw L_e) g_OH; overcoat 1e3 rho_gel / alpha_sw.
    gel = film_cal.gelatin_total
    assert gel[0] == pytest.approx(1e3 * 5.0 / (3.0 * classes.dry_thickness_um) * 4.8e-4)
    assert gel[grid.n_z] == pytest.approx(1e3 * 1.35 / 3.0 * 4.8e-4)
    assert np.all(gel[grid.n_z + 1 :] == 0.0)
    # Spec D1 2.4.2: hazards with the bath's free sulfite.
    free = ens._free_sulfite_bath(d76)
    assert 0.0 < free < d76.total("sulfite")
    np.testing.assert_array_equal(film_cal.hazards, cluster_hazards(kin_d76, 1e-3, 5, free))
    assert film_cal.fog_ratio == 1e-3
    # Spec D1 2.3.7: D(20 C, I) = D_w / (1 + b_eta I).
    k = d76.component_ids.index("hydroquinone")
    assert film_cal.diffusion_um2_s[k] == pytest.approx(855.0 / (1 + 0.125 * d76.ionic_strength))
    assert film_cal.gel_ratio[k] == 0.10
    for name in ("silver_per_cell", "gelatin_total", "diffusion_um2_s", "hazards"):
        assert not getattr(film_cal, name).flags.writeable


def test_film_column_model_profiles_and_overrides(film_material, kin_d76, d76, species):
    image = film_column_model(film_material, "135", 0.01, kin_d76, d76, "image", species=species)
    assert image.grid.n_z == 2 and image.grid.n_l == 4 and image.classes.count == 6
    custom = film_column_model(
        film_material,
        "135",
        0.01,
        kin_d76,
        d76,
        "image",
        species=species,
        n_classes=4,
        quadrature=5,
        latent_nodes=51,
        levels=(3, 6, 1.5),
    )
    assert custom.classes.count == 4
    assert custom.latent.log10_m.size == 51
    assert custom.grid.n_z == 3 and custom.grid.n_l == 6 and custom.grid.n_a == 3
    liquid = custom.grid.thickness_um[4:]
    assert liquid[:3].sum() == pytest.approx(50.0) and liquid.sum() == pytest.approx(300.0)


def test_diffusion_factor(kin_d76):
    assert ens.diffusion_factor(kin_d76.transport, 20.0, 0.0) == pytest.approx(1.0, rel=1e-14)
    expected = (298.15 / 293.15) * (1.002 / 0.890) / (1.0 + 0.125 * 2.0)
    assert ens.diffusion_factor(kin_d76.transport, 25.0, 2.0) == pytest.approx(expected)
    # ln eta linear in 1/T: half way (in 1/T) between 20 and 24 C.
    t = 1.0 / (0.5 * (1 / 293.15 + 1 / 297.15))
    eta = math.sqrt(1.002 * 0.911)
    assert ens._viscosity(kin_d76.transport.water_viscosity, t) == pytest.approx(eta, rel=1e-12)


def test_film_column_model_rejects_bad_input(film_material, kin_d76, d76, species):
    with pytest.raises(ValueError, match="profile"):
        film_column_model(film_material, "135", 0.01, kin_d76, d76, "fine", species=species)
    with pytest.raises(ValueError, match="format"):
        film_column_model(film_material, "220", 0.01, kin_d76, d76, "image", species=species)
    with pytest.raises(ValueError, match="exposure_time_s"):
        film_column_model(film_material, "135", 0.0, kin_d76, d76, "image", species=species)
    with pytest.raises(ValueError, match="levels"):
        film_column_model(
            film_material, "135", 0.01, kin_d76, d76, "image", species=species, levels=(0, 4, 2)
        )


# ---------------------------------------------------------------------------------------
# Rate law, derivatives and the Schur solve
# ---------------------------------------------------------------------------------------


def _setup(film, sol, kin, presoaked=False):
    return ens._Setup(film, sol, kin, presoaked)


def test_surface_matches_kinetics_and_derivatives(film_cal, d76, kin_test):
    st = _setup(film_cal, d76, kin_test)
    rng = np.random.default_rng(3)
    a = rng.uniform(1e-6, 1e-3, (2, 3, 4))
    br = rng.uniform(0.0, 5e-3, (3, 4))
    io = rng.uniform(0.0, 1e-4, (3, 4))
    s = ens._surface(st, a, br, io)
    ref = surface_state(a.reshape(2, -1), br.ravel(), io.ravel(), kin_test, d76.ionic_strength)
    np.testing.assert_allclose(s.kappa.ravel(), ref.kappa_tilde, rtol=1e-13)
    np.testing.assert_allclose(s.lam.ravel(), ref.lambda_tilde, rtol=1e-13)
    np.testing.assert_allclose(s.partition.reshape(2, -1), ref.partition, rtol=1e-13)
    inputs = [a[0], a[1], br, io]
    for u in range(4):
        step = 1e-6 * np.maximum(inputs[u], 1e-6)
        args = [x.copy() for x in inputs]
        args[u] = args[u] + step
        plus = ens._surface(st, np.stack(args[:2]), args[2], args[3])
        args[u] = args[u] - 2 * step
        minus = ens._surface(st, np.stack(args[:2]), args[2], args[3])
        for name, d in (("kappa", s.dkappa[u]), ("lam", s.dlam[u])):
            fd = (getattr(plus, name) - getattr(minus, name)) / (2 * step)
            np.testing.assert_allclose(d, fd, rtol=1e-6, atol=1e-12 * np.abs(fd).max())


def _mid_state(strip_run):
    return strip_run.state, strip_run.setup


def test_rate_gradients_match_finite_difference(strip_run, d76):
    # The rates R_dev, R_diss, R_pd behind the reaction Jacobian (spec D1 3.3): analytic
    # gradients over the components, TOTH and Lambda against central differences.
    state, st = _mid_state(strip_run)
    kern = ens._Kernels(omega=ens._latent(st, *strip_exposures(12, st.n_z)))
    bulk = np.concatenate([d76.totals, [d76.toth]])
    car = ens._bath_carriers(st, bulk)

    def ev(conc, loc):
        return ens._evaluate(
            st,
            st.swelling(450.0),
            conc,
            loc,
            state.g,
            state.h,
            kern,
            bulk,
            car,
            st.L,
            state.ph,
            True,
        )

    ev0 = ev(state.conc, state.loc)
    lev = st.n_z - 1
    for v in range(st.nc + 1):
        conc_p, conc_m = state.conc.copy(), state.conc.copy()
        loc_p, loc_m = state.loc.copy(), state.loc.copy()
        if v < st.nc:
            x = state.conc[:, lev, v]
            step = 1e-6 * np.maximum(np.abs(x), 1e-7)
            conc_p[:, lev, v] += step
            conc_m[:, lev, v] -= step
            col = v
        else:
            step = 1e-6 * np.maximum(state.loc[:, lev, 0], 1e-3)
            loc_p[:, lev, 0] += step
            loc_m[:, lev, 0] -= step
            col = st.iL
        rp, rm = ev(conc_p, loc_p).rates, ev(conc_m, loc_m).rates
        for t_ in range(3):
            fd = (rp[t_][:, lev] - rm[t_][:, lev]) / (2 * step)
            grad = ev0.grad_rates[t_, :, lev, col]
            tol = 1e-5 * np.abs(fd).max() + 1e-12 * np.abs(ev0.rates[t_][:, lev]).max() / 1e-7
            np.testing.assert_allclose(grad, fd, rtol=1e-5, atol=tol)


def _full_matrix(st, jac, kern, gh):
    """Dense ``I - gh J`` of column plus kernels for one batch (test only)."""
    n = jac.coef.shape[0]
    lay = st.layout
    jcc = jac.column.to_dense()  # (N, size, size)
    nz, k, nb = st.n_z, st.K, kern.omega.shape[-1]
    nu = nz * k * nb
    size = lay.size + 3 * nu
    full = np.zeros((n, size, size))
    full[:, : lay.size, : lay.size] = jcc

    def row(lev, r):  # dense row index of level lev, unknown r
        return lev * lay.block + r

    def ku(part, lev, c, b):
        return lay.size + part * nu + (lev * k + c) * nb + b

    for lev in range(nz):
        for c in range(k):
            for b in range(nb):
                ig, ih, iq = ku(0, lev, c, b), ku(1, lev, c, b), ku(2, lev, c, b)
                kap = jac.kappa_c[:, lev, c, 0]
                for r in range(st.b):
                    full[:, row(lev, r), ig] = sum(
                        jac.coef[:, lev, r, t] * jac.g_term[t, :, lev, c] * kern.omega[:, lev, c, b]
                        for t in range(3)
                    )
                for v in range(st.nc):
                    full[:, ig, row(lev, v)] = jac.s_eg[:, lev, c, b] * jac.grad_kappa[:, lev, v]
                    full[:, ih, row(lev, v)] = jac.s_gh[:, lev, c, b] * jac.grad_kappa[:, lev, v]
                full[:, ig, row(lev, st.iL)] = (
                    kap * jac.de[:, lev, 0 if jac.de.shape[2] == 1 else c, b]
                )
                full[:, iq, row(lev, st.i_ag)] = jac.q_ag[:, lev, c, b]
                full[:, ig, ig] = -kap
                full[:, ih, ig] = 2 * kap
                full[:, ih, ih] = -2 * kap
                full[:, iq, ig] = jac.q_g[:, lev, c, 0]
    return np.eye(size)[None] - gh * full


def test_member_schur_elimination_exact(strip_run, d76):
    # Spec D1 2.4.13: the Schur elimination of the kernels equals the fully coupled linear
    # solve (one member per column) to <= 1e-12.
    state, st = _mid_state(strip_run)
    kern = ens._Kernels(omega=ens._latent(st, *strip_exposures(12, st.n_z)))
    bulk = np.concatenate([d76.totals, [d76.toth]])
    car = ens._bath_carriers(st, bulk)
    sl = slice(6, 9)
    sub = ens._Kernels(omega=kern.omega[sl])
    ev = ens._evaluate(
        st,
        st.swelling(450.0),
        state.conc[sl],
        state.loc[sl],
        state.g[sl],
        state.h[sl],
        sub,
        bulk,
        car,
        st.L,
        state.ph[sl],
        True,
    )
    gh = ros.GAMMA * 3.0
    solve = ens._make_solver(st, ev.jac, sub, gh)
    k_sol = solve(ev.f_conc, ev.f_loc, ev.f_g, ev.f_h, ev.f_q)
    full = _full_matrix(st, ev.jac, sub, gh)
    rhs = np.concatenate(
        [
            ens._pack(st, ev.f_conc, ev.f_loc),
            ev.f_g.reshape(3, -1),
            ev.f_h.reshape(3, -1),
            ev.f_q.reshape(3, -1),
        ],
        axis=1,
    )
    x = np.stack([np.linalg.solve(full[i], rhs[i]) for i in range(3)])
    mine = np.concatenate(
        [
            ens._pack(st, k_sol[0], k_sol[1]),
            k_sol[2].reshape(3, -1),
            k_sol[3].reshape(3, -1),
            k_sol[4].reshape(3, -1),
        ],
        axis=1,
    )
    scale = np.abs(x).max(axis=1, keepdims=True)
    assert np.max(np.abs(mine - x) / scale) <= 1e-12


# ---------------------------------------------------------------------------------------
# Conservation and the tank (spec D1 2.4.9, 2.4.13)
# ---------------------------------------------------------------------------------------


def _inventory(st, conc):
    """Film-side inventory per column (mol L^-1 um): every level, held cells included."""
    return np.einsum("nlk,l->nk", conc, st.delta, optimize=False)


def _initial_inventory(run, n):
    st = run.setup
    bulk0 = run.bulk[0]
    s0 = ens._initial_state(st, n, bulk0, run.state.g.shape[-1])
    return _inventory(st, s0.conc)


def _balance(run, result, coef):
    """(change of film side + tank, per area) and the sources from the outputs, per
    component group, all in mol L^-1 um."""
    st = run.setup
    n = result.silver.shape[0]
    film_change = (_inventory(st, run.state.conc) - _initial_inventory(run, n)).mean(axis=0)
    tank_change = (run.bulk[-1] - run.bulk[0]) / coef
    change = film_change + tank_change
    m_coat = st.film.classes.silver_coverage_g_m2
    chem = result.silver.sum(axis=1).mean() / G_PER_MOL_UM
    phys = result.silver_physical.mean() / G_PER_MOL_UM
    dissolved = (run.state.loc[:, :, 1] * st.delta[: st.n_z]).sum(axis=1).mean()
    # The halide output books the same D_z (spec D1 2.4.10), to the rounding of m_coat.
    halide_dissolved = (m_coat - result.silver.sum(axis=1) - result.halide).mean() / G_PER_MOL_UM
    assert halide_dissolved == pytest.approx(dissolved, abs=1e-13 * m_coat / G_PER_MOL_UM)
    return change, chem, phys, dissolved


def test_silver_conservation(strip_run):
    # Spec D1 2.4.13: chemical + physical + halide + dissolved silver (in the film, the liquid
    # and the tank) = coating, <= 1e-12.
    run = strip_run
    res = run.result
    coef = tk.tank_coefficient(film_area_um2=ROLL_TANK.film_area_um2, volume_l=ROLL_TANK.volume_l)
    change, _, phys, dissolved = _balance(run, res, coef)
    in_solution = change[run.setup.i_ag] * G_PER_MOL_UM  # g/m2; the bath starts silver-free
    total = (
        res.silver.sum(axis=1).mean() + res.silver_physical.mean() + res.halide.mean() + in_solution
    )
    assert total == pytest.approx(run.setup.film.classes.silver_coverage_g_m2, rel=1e-12)
    assert dissolved > 0.0 and phys > 0.0 and in_solution > 0.0


def _check_components(run, coef):
    res = run.result
    st = run.setup
    change, chem, phys, dissolved = _balance(run, res, coef)
    ids = st.ids
    br, io, ag = change[st.i_br], change[st.i_io], change[st.i_ag]
    x_i = st.x_i
    assert br == pytest.approx((1 - x_i) * chem + dissolved, rel=1e-12)
    assert io == pytest.approx(x_i * chem, rel=1e-12)
    assert ag == pytest.approx(dissolved - phys, rel=1e-12, abs=1e-12 * dissolved)
    # Developer electrons: each developed Ag takes half a developer molecule (spec D1 2.4.6).
    dev = sum(change[ids.index(a)] for a in st.kin.agent_ids)
    assert dev == pytest.approx(-0.5 * (chem + phys), rel=1e-12)
    # Sulfite and TOTH: d(sulfite + n_lig Ag - TOTH) = -R_e.
    combo = change[st.i_sulf] + st.n_lig * ag - change[st.it]
    assert combo == pytest.approx(-(chem + phys), rel=1e-11)
    assert chem > 0.1


@pytest.mark.parametrize("agitation", [SMALL_TANK, LARGE_TANK], ids=["small", "large"])
def test_component_conservation(film_cal, d76, kin_test, agitation):
    # Spec D1 2.4.13: film + liquid + tank inventory change = integrated sources, <= 1e-12,
    # with the held upper section counted on the film side (module docstring).
    schedule = make_schedule("calibration", 200.0, agitation)
    ev, es = strip_exposures(5, film_cal.grid.n_z, lo=0.0, hi=3.0)
    run = _integrate(film_cal, d76, kin_test, schedule, ev, es)
    coef = tk.tank_coefficient(film_area_um2=ROLL_TANK.film_area_um2, volume_l=ROLL_TANK.volume_l)
    _check_components(run, coef)


def test_component_conservation_ends_agitated(film_cal, d76, kin_test):
    schedule = make_schedule("calibration", 92.0, SMALL_TANK)  # ends inside 90-95 s
    assert schedule.agitated[-1]
    ev, es = strip_exposures(3, film_cal.grid.n_z, lo=1.0, hi=3.0)
    run = _integrate(film_cal, d76, kin_test, schedule, ev, es)
    assert run.held_final
    coef = tk.tank_coefficient(film_area_um2=ROLL_TANK.film_area_um2, volume_l=ROLL_TANK.volume_l)
    _check_components(run, coef)


def test_tank_mass_balance(strip_run):
    # Spec D1 2.4.13 / 3.7: the bulk trajectory is the tank ODE re-integrated from flux_mean
    # (mixing impulses included), bit for bit; and it changes.
    res = strip_run.result
    coef = tk.tank_coefficient(film_area_um2=ROLL_TANK.film_area_um2, volume_l=ROLL_TANK.volume_l)
    traj = np.concatenate([res.bulk.totals, res.bulk.toth[:, None]], axis=1)
    again = tk.integrate_bulk(traj[0], res.flux_mean, res.bulk.times_s, coef)
    np.testing.assert_array_equal(again, traj)
    st = strip_run.setup
    assert traj[-1, st.i_br] > 1e-4  # released bromide reaches the bath
    assert traj[-1, st.ids.index("metol")] < traj[0, st.ids.index("metol")]


def test_tank_infinite_volume(film_cal, d76, kin_test):
    # Spec D1 2.4.13: V = 1e9 L leaves the bath unchanged to <= 1e-12 mol/L.
    schedule = make_schedule("calibration", 120.0)
    ev, es = strip_exposures(3, film_cal.grid.n_z, lo=1.0, hi=3.0)
    tank = TankSpec(volume_l=1e9, film_area_um2=ROLL_TANK.film_area_um2)
    res = run_strip(film_cal, d76, kin_test, schedule, ev, es, tank=tank)
    traj = np.concatenate([res.bulk.totals, res.bulk.toth[:, None]], axis=1)
    # Absolute change in mol/L (the bromide-free bath has no relative scale).
    assert np.max(np.abs(traj - traj[0])) <= 1e-12
    assert np.max(np.abs(traj - traj[0])) > 0.0


def test_bulk_mode_books_the_same_tank(film_cal, d76, kin_test):
    # A given trajectory: outputs are per column and flux_mean re-integrates to the tank the
    # columns would have produced (spec D1 3.7 step 4 relies on this).
    schedule = make_schedule("calibration", 100.0)
    ev, es = strip_exposures(3, film_cal.grid.n_z, lo=1.0, hi=3.0)
    coupled = run_strip(film_cal, d76, kin_test, schedule, ev, es)
    given = run_strip(film_cal, d76, kin_test, schedule, ev, es, tank=None, bulk=coupled.bulk)
    # The coupled run predicts the end-of-step bulk for the second stage; a given trajectory
    # interpolates it: the two agree to second order in the step.
    np.testing.assert_allclose(developed_silver(given), developed_silver(coupled), rtol=5e-4)
    coef = tk.tank_coefficient(film_area_um2=ROLL_TANK.film_area_um2, volume_l=ROLL_TANK.volume_l)
    traj = np.concatenate([coupled.bulk.totals, coupled.bulk.toth[:, None]], axis=1)
    again = tk.integrate_bulk(traj[0], given.flux_mean, given.bulk.times_s, coef)
    change = np.abs(traj[-1] - traj[0]) + 1e-12
    assert np.max(np.abs(again[-1] - traj[-1]) / change) < 1e-2


# ---------------------------------------------------------------------------------------
# Members, kernels and vectorisation
# ---------------------------------------------------------------------------------------


def test_ensemble_vectorisation(film_image, d76, kin_test):
    # Spec D1 2.4.13: N independent columns (given bulk) equal the one-by-one runs bit for bit.
    schedule = make_schedule("image", 60.0)
    bulk = _constant_bulk(d76, schedule)
    ev, es = strip_exposures(4, film_image.grid.n_z, lo=0.0, hi=3.0)
    es = es + 0.1 * ev
    batch = run_strip(film_image, d76, kin_test, schedule, ev, es, tank=None, bulk=bulk)
    for i in range(4):
        one = run_strip(
            film_image, d76, kin_test, schedule, ev[i : i + 1], es[i : i + 1], tank=None, bulk=bulk
        )
        for name in (
            "silver",
            "initiated",
            "silver_m2",
            "silver_physical",
            "halide",
            "growth_rate",
        ):
            np.testing.assert_array_equal(getattr(one, name)[0], getattr(batch, name)[i])


def test_kernel_reduction_exact(film_cal, d76, kin_test):
    # Spec D1 2.4.13: X, Y, P from the kernels equal the directly integrated member equations
    # (one member per column) to <= 1e-12, and the cell chemistry is the same.
    schedule = make_schedule("calibration", 150.0)
    ev, es = strip_exposures(4, film_cal.grid.n_z, lo=0.0, hi=3.0)
    kern_run = _integrate(film_cal, d76, kin_test, schedule, ev, es)
    st = kern_run.setup
    p = ens._latent(st, ev, es)  # (N, nz, K, Nn)
    direct = ens._Kernels(omega=np.ones((*p.shape[:3], 1)), basis=p[:, :, :, None, :])
    direct_run = _integrate(film_cal, d76, kin_test, schedule, ev, es, kernels=direct)
    s_k, s_d = kern_run.state, direct_run.state
    x_k = (p * s_k.g).sum(axis=-1)
    y_k = (p * s_k.h).sum(axis=-1)
    q_k = (p * s_k.q).sum(axis=-1)
    assert np.max(np.abs(x_k - s_d.g[..., 0])) <= 1e-12
    assert np.max(np.abs(y_k - s_d.h[..., 0])) <= 1e-12
    assert np.max(np.abs(q_k - s_d.q[..., 0])) <= 1e-12 * max(np.abs(q_k).max(), 1e-30) + 1e-18
    scale = np.maximum(np.abs(s_k.conc), 1e-6)
    assert np.max(np.abs(s_k.conc - s_d.conc) / scale) <= 1e-12


# ---------------------------------------------------------------------------------------
# Limits (spec D1 2.4.13)
# ---------------------------------------------------------------------------------------


def test_zero_exposure_fog_only(film_cal, d76, kin_test, film_material, species):
    schedule = make_schedule("calibration", 300.0)
    zero = np.zeros((2, film_cal.grid.n_z))
    res = run_strip(film_cal, d76, kin_test, schedule, zero, zero)
    # Unexposed grains develop only through the fog hazard h_0 = h_fog: identical columns.
    np.testing.assert_array_equal(res.silver[0], res.silver[1])
    assert res.silver.sum() > 0.0
    no_fog_material = _with_development(film_material, fog_hazard_ratio=0.0)
    film0 = film_column_model(
        no_fog_material, "135", EXPOSURE_TIME_S, kin_test, d76, "calibration", species=species
    )
    res0 = run_strip(film0, d76, kin_test, schedule, zero, zero)
    # p(n >= 2 | m = 0) is the first table node (< 1e-12, spec D1 2.2.5), not exactly 0.
    assert np.all(res0.silver <= 1e-12 * film_cal.classes.silver_coverage_g_m2)
    assert np.all(res0.initiated <= 1e-11)


def _with_development(material, **values):

    layers = []
    for layer in material.layers:
        if layer.type == "emulsion":
            dev = layer.development.model_copy(
                update={k: Quantity(value=v, unit="1", source="assumed") for k, v in values.items()}
            )
            layer = layer.model_copy(update={"development": dev})
        layers.append(layer)
    return material.model_copy(update={"layers": layers})


def test_long_development_limit(film_cal, d76, kin_test):
    # Spec D1 2.4.13: strong exposure, ample developer, 5000 s: developable classes X -> 1
    # within 1e-3 (X = silver_c / coated silver of class c).
    schedule = make_schedule("calibration", 5000.0, CONTINUOUS)
    ev = np.full((1, film_cal.grid.n_z), 1e4)
    res = run_strip(film_cal, d76, kin_test, schedule, ev, np.zeros_like(ev), tank=INFINITE_TANK)
    classes = film_cal.classes
    coated = 1e12 * classes.silver_g_per_grain * classes.number_per_um2
    x = res.silver[0] / coated
    assert np.all(np.abs(1.0 - x) <= 1e-3), x
    assert np.all(res.initiated[0] >= 1.0 - 1e-3)


def test_gel_equilibrium_no_toth_flux(film_cal, d76, kin_test):
    # Spec D1 2.4.13: unexposed, h_fog = 0, k_diss = k_pd = 0, s = 1; gel cells hold the bath's
    # mobile totals and TOTH = TOTH_mob,bath + T_gel nbar_gel(pH_bath); after 300 s every level
    # keeps the bath pH (<= 1e-10) and every face carries |J_toth| <= 1e-14.
    kin = with_solvent(kin_test, k_diss=(0.0, "um/s"), k_pd=(0.0, "um/s"))
    film = dataclasses.replace(film_cal, hazards=np.zeros_like(film_cal.hazards))
    st = _setup(film, d76, kin, presoaked=True)
    n = 1
    conc = np.broadcast_to(np.asarray(d76.totals), (n, st.L, st.kc)).copy()
    gel_pk = st.spec.pka[st.nf, 0]
    nbar_gel = 1.0 / (1.0 + 10.0 ** (d76.ph - gel_pk))
    toth = np.full((n, st.L), d76.toth) + st.gelatin[None] * nbar_gel
    schedule = make_schedule("calibration", 300.0, SMALL_TANK)
    zero = np.zeros((n, st.n_z))
    run = ens._integrate(
        film,
        d76,
        kin,
        schedule,
        columns(zero, zero),
        n_cells=n,
        cell_weights=np.ones(n),
        tank=INFINITE_TANK,
        presoaked=True,
        initial=(conc, toth),
    )
    state = run.state
    sp = ens._speciate(
        st,
        positive_part(state.conc[:, :, : st.kc]),
        state.conc[:, :, st.it],
        st.gelatin,
        state.ph,
    )
    assert np.max(np.abs(sp.ph - d76.ph)) <= 1e-10
    carriers = [state.conc[:, :, st.fam_comp[f]] * sp.nbar[f] for f in st.mobile]
    carriers += [sp.h, -sp.oh]
    d_car = st.carrier_diffusion(1.0)
    flux = np.zeros(st.L - 1)
    for c, car in enumerate(carriers):
        face = ens._face_mean(d_car[c], st.delta)
        flux = flux - face * (car[0, 1:] - car[0, :-1]) / st.hf
    assert np.max(np.abs(flux)) <= 1e-14
    assert np.all(run.result.silver == 0.0)


def _material_with(material, *, gelatin=None, overcoat=None):

    layers = []
    for layer in material.layers:
        if layer.type == "emulsion" and gelatin is not None:
            q = Quantity(value=gelatin, unit="g/m2", source="assumed")
            layer = layer.model_copy(update={"gelatin_coverage": q})
        if layer.type == "overcoat" and overcoat is not None:
            layer = layer.model_copy(
                update={"thickness": Quantity(value=overcoat, unit="um", source="assumed")}
            )
        layers.append(layer)
    return material.model_copy(update={"layers": layers})


def _borate_bath(library):

    comp = BathComponent(
        species="borax-decahydrate", amount=Quantity(value=20.0, unit="g/L", source="assumed")
    )
    return library.bath("d-76").model_copy(
        update={"id": "test/borate", "components": [comp], "nuisance": []}
    )


def test_immobile_buffer_retards_ph_front(library, species, kfile, film_material):
    # Spec D1 2.4.13: one mobile buffer (borate) plus gelatin, no reaction, s = 1, a thick gel;
    # a small pH step at the surface spreads as erf(x / 2 sqrt(D_app t)) with
    # D_app = sum_c D_c beta_c / beta (= D / (1 + beta_gel / beta_mob) up to the water ions),
    # within 1e-2 of the step.

    sol = working_solution(
        _borate_bath(library), species, kfile, dilution="stock", temperature_c=20.0
    )
    kin = resolve_kinetics(kfile, sol, species)
    assert kin.agent_ids == ()
    kin = with_solvent(kin, k_diss=(0.0, "um/s"), k_pd=(0.0, "um/s"))
    kin = with_transport(kin, delta_agit=(1e-3, "um"), delta_max=(2e-3, "um"))
    material = _material_with(film_material, gelatin=100.0, overcoat=1e-3)
    film = film_column_model(
        material, "135", 0.01, kin, sol, "image", species=species, levels=(100, 4, 2.0)
    )
    film = dataclasses.replace(film, hazards=np.zeros_like(film.hazards))
    st = _setup(film, sol, kin, presoaked=True)
    # Gel equilibrated with the bath at pH_bath + a: TOTH lower by beta_gel-cell * a.
    amp = 1e-3
    ph_gel = sol.ph + amp
    totals = np.zeros((st.nf + 1, 1, st.L))
    for f in range(st.nf):
        totals[f] = sol.totals[st.fam_comp[f]]
    totals[st.nf] = st.gelatin[None]

    nbar = mean_protons(st.spec, np.full((1, st.L), ph_gel))
    hh, oh = water_ions(st.spec, np.full((1, st.L), ph_gel))
    toth = (totals * nbar).sum(axis=0) + hh - oh
    conc = np.broadcast_to(np.asarray(sol.totals), (1, st.L, st.kc)).copy()
    conc[:, ~st.gel] = sol.totals
    toth[:, ~st.gel] = sol.toth
    t_end = 4.0
    schedule = step_schedule(
        agitation_periods(CONTINUOUS, t_end), t_end, h0_s=0.002, ratio=1.1, h_max_s=0.02
    )
    zero = np.zeros((1, st.n_z))
    run = ens._integrate(
        film,
        sol,
        kin,
        schedule,
        columns(zero, zero),
        n_cells=1,
        cell_weights=np.ones(1),
        tank=INFINITE_TANK,
        presoaked=True,
        initial=(conc, toth),
    )
    state = run.state
    sp = ens._speciate(
        st,
        positive_part(state.conc[:, :, : st.kc]),
        state.conc[:, :, st.it],
        st.gelatin,
        state.ph,
    )
    # D_app from the carriers at the gel's pH (emulsion levels).
    ph0 = np.full((1, st.L), ph_gel)
    sp0 = ens._speciate(st, positive_part(conc), toth, st.gelatin, ph0)
    d_car = st.carrier_diffusion(1.0)[:, 0]
    beta_c = [ens._LN10 * sol.totals[st.fam_comp[f]] * sp0.var[f, 0, 0] for f in st.mobile]
    beta_c += [ens._LN10 * sp0.h[0, 0], ens._LN10 * sp0.oh[0, 0]]
    d_app = sum(d * b for d, b in zip(d_car, beta_c, strict=True)) / sp0.beta[0, 0]
    beta_gel = sp0.beta[0, 0] - sum(beta_c)
    assert 0.1 < beta_gel / sum(beta_c) < 10.0  # the gelatin matters
    top = st.grid.thickness_um[: st.n_z + 1].sum()
    x = top - st.grid.height_um[: st.n_z]
    expected = sol.ph + amp * erf(x / (2.0 * math.sqrt(d_app * t_end)))
    err = np.abs(sp.ph[0, : st.n_z] - expected)
    assert err.max() <= 1e-2 * amp, err.max() / amp
    # And it is slower than the front without gelatin.
    d_mob = sum(d * b for d, b in zip(d_car, beta_c, strict=True)) / sum(beta_c)
    fast = sol.ph + amp * erf(x / (2.0 * math.sqrt(d_mob * t_end)))
    assert np.abs(sp.ph[0, : st.n_z] - fast).max() > 3.0 * err.max()


def test_presoak_removes_delay(film_cal, d76, kin_test):
    # Spec D1 2.4.13: the dry film starts slower; as T_half -> 0 dry equals presoaked (<= 1e-6).
    schedule = make_schedule("calibration", 120.0, SMALL_TANK, record_s=(20.0,))
    ev, es = strip_exposures(3, film_cal.grid.n_z, lo=1.0, hi=3.0)
    dry = run_strip(film_cal, d76, kin_test, schedule, ev, es)
    wet = run_strip(film_cal, d76, kin_test, schedule, ev, es, presoaked=True)
    assert np.all(developed_silver(dry.snapshots[20.0]) < developed_silver(wet.snapshots[20.0]))
    fast = with_transport(kin_test, swelling_half_time=(1e-9, "s"))
    dry0 = run_strip(film_cal, d76, fast, schedule, ev, es)
    np.testing.assert_allclose(dry0.silver, wet.silver, rtol=0, atol=1e-6)
    np.testing.assert_allclose(dry0.initiated, wet.initiated, rtol=0, atol=1e-6)


def _with_bromide(library, grams):

    bath = library.bath("d-76")
    comp = BathComponent(
        species="potassium-bromide", amount=Quantity(value=grams, unit="g/L", source="assumed")
    )
    return bath.model_copy(update={"components": [*bath.components, comp]})


def test_bromide_restrains(library, species, kfile, film_cal, kin_test):
    # Spec D1 2.4.13: with initial bromide, X at a fixed time falls monotonically and the
    # initiation clock Lambda grows more slowly.

    schedule = make_schedule("calibration", 240.0)
    ev, es = strip_exposures(2, film_cal.grid.n_z, lo=1.5, hi=2.5)
    silver, clocks = [], []
    for grams in (1e-6, 0.5, 2.0):
        sol = working_solution(
            _with_bromide(library, grams), species, kfile, dilution="stock", temperature_c=20.0
        )
        kin = with_test_rates(resolve_kinetics(kfile, sol, species))
        film = dataclasses.replace(film_cal, component_ids=sol.component_ids)
        run = _integrate(film, sol, kin, schedule, ev, es, record_cells=4)
        silver.append(developed_silver(run.result))
        clocks.append(run.result.cell_history.lambda_clock[-1])
    assert np.all(np.diff(np.array(silver), axis=0) < 0.0)
    assert np.all(np.diff(np.array(clocks), axis=0) < 0.0)


@pytest.mark.parametrize("profile", ["image", "calibration"])
def test_dry_film_bromide_free_start(profile, film_material, kin_test, d76, species):
    # Spec D1 2.4.13: dry film, bromide-free D-76, image (h0 = 0.1 s) and calibration
    # schedules: no NaN, min C/s_k >= -1e-3 over the run, silver balance holds.
    film = film_column_model(film_material, "135", 0.01, kin_test, d76, profile, species=species)
    schedule = make_schedule(profile, 450.0)
    ev, es = strip_exposures(4, film.grid.n_z)
    run = _integrate(film, d76, kin_test, schedule, ev, es)
    res = run.result
    for name in ("silver", "initiated", "silver_m2", "silver_physical", "halide", "growth_rate"):
        assert np.all(np.isfinite(getattr(res, name)))
    assert np.all(np.isfinite(run.state.conc))
    assert -1e-3 <= res.min_scaled_concentration < 0.0  # the dry start dips below zero
    coef = tk.tank_coefficient(film_area_um2=ROLL_TANK.film_area_um2, volume_l=ROLL_TANK.volume_l)
    change, _, _, _ = _balance(run, res, coef)
    total = res.silver.sum(axis=1).mean() + res.silver_physical.mean() + res.halide.mean()
    total += change[run.setup.i_ag] * G_PER_MOL_UM
    assert total == pytest.approx(film.classes.silver_coverage_g_m2, rel=1e-12)


def test_arrhenius_scaling(film_cal, d76, kfile, species, library):
    # Spec D1 2.4.13 at the engine level: one developing agent, temperature-independent
    # speciation (the 20 C solution), Delta H_Br = 0, no bromide, no dissolution or physical
    # development, very fast transport: the time to reach a given X scales by
    # exp(Ea/R (1/T2 - 1/T1)), to <= 1e-3. Running the warm bath on the cold schedule scaled
    # by that factor must reproduce the cold run.

    sol = working_solution(
        library.bath("d-76"),
        species,
        kfile,
        dilution="stock",
        temperature_c=20.0,
        overrides={"hydroquinone": 0.0},
    )
    zero_dh = Quantity(value=0.0, unit="kJ/mol", source="assumed")
    halides = kfile.halides.model_copy(update={"adsorption_enthalpy": zero_dh})
    kf = kfile.model_copy(update={"halides": halides})
    cold = with_test_rates(resolve_kinetics(kf, sol, species))
    assert cold.agent_ids == ("metol",)
    warm_sol = dataclasses.replace(sol, temperature_c=24.0)
    warm = with_test_rates(resolve_kinetics(kf, warm_sol, species))
    a = arrhenius_factor(61.8, 297.15)
    warm = dataclasses.replace(
        warm, k_et=cold.k_et * a, lambda0=cold.lambda0 * a, k_des=cold.k_des * a
    )
    assert warm.K_br == cold.K_br
    cold = with_solvent(cold, k_diss=(0.0, "um/s"), k_pd=(0.0, "um/s"))
    warm = with_solvent(warm, k_diss=(0.0, "um/s"), k_pd=(0.0, "um/s"))
    film = dataclasses.replace(
        film_cal,
        component_ids=sol.component_ids,
        diffusion_um2_s=film_cal.diffusion_um2_s * 1e3,
        water_ion_diffusion_um2_s=tuple(1e3 * d for d in film_cal.water_ion_diffusion_um2_s),
    )
    t_cold = 300.0
    base = make_schedule("calibration", t_cold, CONTINUOUS)
    scaled = Schedule(times_s=base.times_s / a, agitated=base.agitated, record_steps={})
    ev, es = strip_exposures(3, film.grid.n_z, lo=1.0, hi=3.0)
    r_cold = run_strip(film, sol, cold, base, ev, es, tank=INFINITE_TANK, presoaked=True)
    r_warm = run_strip(film, sol, warm, scaled, ev, es, tank=INFINITE_TANK, presoaked=True)
    m_cold, m_warm = developed_silver(r_cold), developed_silver(r_warm)
    assert np.all(m_cold > 0.05)
    # Equal X at the scaled time <=> the time to a given X scales by 1/a; dX/dlnt converts.
    np.testing.assert_allclose(m_warm, m_cold, rtol=1e-3)


def test_snapshots_equal_shorter_runs(strip_run, film_cal, d76, kin_test):
    res = strip_run.result
    assert sorted(res.snapshots) == [120.0, 300.0]
    schedule = make_schedule("calibration", 300.0, SMALL_TANK)
    ev, es = strip_exposures(12, film_cal.grid.n_z)
    short = run_strip(film_cal, d76, kin_test, schedule, ev, es)
    snap = res.snapshots[300.0]
    for name in ("silver", "initiated", "silver_m2", "silver_physical", "halide", "growth_rate"):
        np.testing.assert_array_equal(getattr(snap, name), getattr(short, name))
    np.testing.assert_array_equal(snap.flux_mean, short.flux_mean)
    np.testing.assert_array_equal(snap.bulk.totals, short.bulk.totals)
    assert np.all(developed_silver(res.snapshots[120.0]) < developed_silver(snap))


def test_cell_histories(strip_run):
    res = strip_run.result
    hist = res.cell_history
    st = strip_run.setup
    times = strip_run.result.bulk.times_s
    assert hist.times_s[0] == 0.0 and hist.times_s[-1] == times[-1]
    assert np.all(np.isin(hist.times_s, times))
    assert 8 <= hist.times_s.size <= 16
    assert hist.lambda_clock.shape == (hist.times_s.size, 12, st.n_z)
    assert np.all(np.diff(hist.lambda_clock, axis=0) >= 0.0)
    assert np.all(np.diff(hist.kappa_integral, axis=0) >= 0.0)
    np.testing.assert_array_equal(hist.lambda_clock[-1], strip_run.state.loc[:, ::-1, 0])
    assert hist.cell_spacing_um == ens.CHEM_CELL_UM
    assert np.all(hist.kappa_integral[-1] > 0.0)


def test_strip_outputs_are_read_only_and_shaped(strip_run):
    res = strip_run.result
    assert res.silver.shape == (12, 6) and res.halide.shape == (12,)
    assert res.flux_mean.shape == (res.bulk.times_s.size - 1, len(strip_run.setup.ids) + 1)
    for name in ("silver", "initiated", "silver_m2", "silver_physical", "halide", "flux_mean"):
        assert not getattr(res, name).flags.writeable
    m = developed_silver(res)
    assert np.all(np.diff(m) >= 0.0)  # more exposure, more silver
    assert np.all(res.initiated <= 1.0) and np.all(res.silver_m2 <= res.initiated + 1e-12)


def test_run_ensemble_rejects_bad_input(film_cal, d76, kin_test):
    schedule = make_schedule("calibration", 30.0)
    ev, es = strip_exposures(2, film_cal.grid.n_z)
    mem = columns(ev, es)
    ok = {"n_cells": 2, "cell_weights": np.ones(2), "tank": ROLL_TANK}

    def call(**kw):
        args = {**ok, **kw}
        members = args.pop("members", mem)
        return run_ensemble(film_cal, d76, kin_test, schedule, members, **args)

    with pytest.raises(ValueError, match="exactly one"):
        call(tank=None)
    with pytest.raises(ValueError, match="exactly one"):
        call(bulk=_constant_bulk(d76, schedule))
    with pytest.raises(ValueError, match="method"):
        call(method="rk4")
    with pytest.raises(ValueError, match="rtol"):
        call(rtol=0.0)
    with pytest.raises(ValueError, match="record_cells"):
        call(record_cells=-1)
    with pytest.raises(ValueError, match="cell_weights"):
        call(cell_weights=np.ones(3))
    with pytest.raises(ValueError, match="positive sum"):
        call(cell_weights=np.zeros(2))
    with pytest.raises(ValueError, match="columns layout"):
        call(n_cells=3, cell_weights=np.ones(3))
    with pytest.raises(ValueError, match=">= 0"):
        call(members=columns(-ev, es))
    with pytest.raises(ValueError, match="layout"):
        call(members=Members(ev, es, layout="grid"))
    with pytest.raises(ValueError, match="share a shape"):
        call(members=columns(ev, es[:, :1]))
    other = make_schedule("calibration", 40.0)
    with pytest.raises(ValueError, match="times_s"):
        call(tank=None, bulk=_constant_bulk(d76, other))
    with pytest.raises(ValueError, match="shape_cells"):
        call(lateral=Lateral((1, 3), 10.0, "periodic"))
    with pytest.raises(ValueError, match="raster"):
        call(members=Members(ev[None], es[None], layout="raster", pixel_um=6.0))
    bad_film = dataclasses.replace(film_cal, component_ids=film_cal.component_ids[::-1])
    with pytest.raises(ValueError, match="components"):
        run_ensemble(bad_film, d76, kin_test, schedule, mem, **ok)


# ---------------------------------------------------------------------------------------
# Lateral coupling, compensation and the split (spec D1 2.4.13, 3.2, 3.5, 3.9)
# ---------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def edge_run(film_image, d76, kin_test):
    """A periodic strip, 128 cells of 12 um: a high-exposure band over half the period."""
    nx = 128
    i = np.arange(nx)
    level = np.where((i >= nx // 4) & (i < 3 * nx // 4), 300.0, 3.0)
    ev = np.repeat(level[:, None], film_image.grid.n_z, axis=1)
    schedule = make_schedule("image", 450.0)
    res = run_strip(
        film_image,
        d76,
        kin_test,
        schedule,
        ev,
        np.zeros_like(ev),
        lateral=Lateral((1, nx), 12.0, "periodic"),
    )
    return nx, developed_silver(res)


def test_edge_shows_adjacency(edge_run):
    # Spec D1 2.4.13 (adjacency): near the edge the dense side overshoots its plateau (fresh
    # developer from the thin side) and the thin side undershoots its plateau (bromide from
    # the dense side). Plateaus are taken in the middle of each band, 384 um from the edges.
    nx, m = edge_run
    edge = nx // 4
    high_plateau, low_plateau = m[nx // 2], m[0]
    high_near = m[edge : edge + 8].max()
    low_near = m[edge - 8 : edge].min()
    assert high_near > 1.02 * high_plateau
    assert low_near < 0.98 * low_plateau
    # The pattern is mirror-symmetric about the band centre, and so is the result.
    np.testing.assert_allclose(m, m[::-1], rtol=1e-9)


def test_lateral_uniform_strip_matches_columns(film_image, d76, kin_test):
    # Lateral diffusion leaves a uniform strip uniform: it equals the column without lateral.
    schedule = make_schedule("image", 120.0)
    ev = np.full((8, film_image.grid.n_z), 100.0)
    lat = run_strip(
        film_image,
        d76,
        kin_test,
        schedule,
        ev,
        np.zeros_like(ev),
        lateral=Lateral((2, 4), 24.0, "neumann"),
    )
    col = run_strip(film_image, d76, kin_test, schedule, ev[:1], np.zeros_like(ev[:1]))
    np.testing.assert_allclose(
        lat.silver, np.broadcast_to(col.silver, lat.silver.shape), rtol=1e-10
    )


def test_dilute_developer_shows_compensation(library, species, kfile, film_material):
    # Compensation: in D-76 1+1 developed longer, the net shoulder silver relative to the net
    # toe silver is lower than in stock (local exhaustion restrains the dense steps).
    ratios = []
    for dilution, minutes in (("stock", 7.5), ("1+1", 11.0)):
        sol = make_solution(library, species, kfile, dilution=dilution)
        kin = with_test_rates(resolve_kinetics(kfile, sol, species))
        film = film_column_model(film_material, "135", 0.01, kin, sol, "image", species=species)
        ev, es = strip_exposures(10, film.grid.n_z)
        res = run_strip(film, sol, kin, make_schedule("image", 60.0 * minutes), ev, es)
        net = developed_silver(res) - developed_silver(res)[0]
        assert net[5] > 0.05
        ratios.append(net[-1] / net[5])
    assert ratios[1] < 0.95 * ratios[0], ratios


def test_splitting_second_order(film_image, d76, kin_test):
    # Spec D1 3.9: periodic strip, fixed steps halved: successive differences fall by 3.5-4.5.
    # A lateral bromide wave over a presoaked, bath-equilibrated column (no reactions), so the
    # only lateral-vertical coupling is the split itself.
    kin = with_solvent(kin_test, k_diss=(0.0, "um/s"), k_pd=(0.0, "um/s"))
    film = dataclasses.replace(film_image, hazards=np.zeros_like(film_image.hazards))
    st = _setup(film, d76, kin, presoaked=True)
    nx = 16
    conc = np.broadcast_to(np.asarray(d76.totals), (nx, st.L, st.kc)).copy()
    conc[:, :, st.i_br] = 5e-3 * (1.0 + np.cos(2.0 * np.pi * np.arange(nx) / nx))[:, None]
    nbar = 1.0 / (1.0 + 10.0 ** (d76.ph - st.spec.pka[st.nf, 0]))
    toth = np.full((nx, st.L), d76.toth) + st.gelatin[None] * nbar
    zero = np.zeros((nx, st.n_z))
    fields = []
    for steps in (40, 80, 160, 320):
        schedule = Schedule(
            times_s=np.linspace(0.0, 20.0, steps + 1),
            agitated=np.zeros(steps, bool),
            record_steps={},
        )
        run = ens._integrate(
            film,
            d76,
            kin,
            schedule,
            columns(zero, zero),
            n_cells=nx,
            cell_weights=np.ones(nx),
            tank=INFINITE_TANK,
            presoaked=True,
            initial=(conc, toth),
            lateral=Lateral((1, nx), 10.0, "periodic"),
        )
        fields.append(run.state.conc[:, :, st.i_br].copy())
    d = [np.max(np.abs(a - b)) for a, b in itertools.pairwise(fields)]
    ratios = [d[i] / d[i + 1] for i in range(len(d) - 1)]
    assert all(3.5 <= r <= 4.5 for r in ratios), ratios


# ---------------------------------------------------------------------------------------
# Raster members (spec D1 2.4.4, 3.6, 3.8)
# ---------------------------------------------------------------------------------------


def test_raster_pbar_is_the_area_mean(film_image, d76, kin_test):
    st = _setup(film_image, d76, kin_test)
    rng = np.random.default_rng(5)
    ev = 10.0 ** rng.uniform(-1, 3, (8, 12, st.n_z))
    es = np.zeros_like(ev)
    ras = ens._Raster(
        2, 3, 24.0, 6.0, (0.0, 0.0), 8, 12, (slice(0, 2), slice(0, 3)), (0, 8), (0, 12)
    )
    pbar = ens._raster_pbar(st, ras, ev, es).reshape(
        2, 3, *ens._latent(st, ev[:1, :1], es[:1, :1]).shape[2:]
    )
    p = ens._latent(st, ev, es)
    for a in range(2):
        for b in range(3):
            block = p[4 * a : 4 * a + 4, 4 * b : 4 * b + 4].mean(axis=(0, 1))
            np.testing.assert_allclose(pbar[a, b], block, rtol=0, atol=1e-13)
    # Splitting every pixel into four identical 3 um pixels leaves p-bar unchanged.
    fine = np.repeat(np.repeat(ev, 2, axis=0), 2, axis=1)
    ras2 = dataclasses.replace(ras, pixel=3.0, h=16, w=24, rows=(0, 16), cols=(0, 24))
    pbar2 = ens._raster_pbar(st, ras2, fine, np.zeros_like(fine)).reshape(pbar.shape)
    assert np.max(np.abs(pbar2 - pbar)) <= 1e-13
    # Cell area outside the pixel window counts as unexposed film (spec D1 3.6).
    ras3 = dataclasses.replace(ras, h=6, rows=(0, 6))
    pbar3 = ens._raster_pbar(st, ras3, ev[:6], es[:6]).reshape(pbar.shape)
    zero = ens._zero_latent(st)
    expected = 0.5 * p[4:6, 4:8].mean(axis=(0, 1)) + 0.5 * zero
    np.testing.assert_allclose(pbar3[1, 1], expected, rtol=0, atol=1e-13)


def test_kernel_reduction_multiple_members(film_cal, d76, kin_test):
    # Spec D1 2.4.13: several members per cell, nearest-cell coupling: X, Y, P synthesised
    # from the kernels equal the directly integrated member equations (<= 1e-12), and the cell
    # chemistry depends on p-bar only (reordering the members changes nothing, <= 1e-13).
    schedule = make_schedule("calibration", 150.0)
    rng = np.random.default_rng(11)
    ev = 10.0 ** rng.uniform(0.0, 3.0, (3, 4, film_cal.grid.n_z))  # 3 cells x 4 members
    es = np.zeros_like(ev)
    st = _setup(film_cal, d76, kin_test)
    p = ens._latent(st, ev, es)  # (3, 4, nz, K, Nn)
    pbar = p.mean(axis=1)
    ev_cell = ev.mean(axis=1)
    kern = ens._Kernels(omega=np.ascontiguousarray(pbar))
    k_run = _integrate(film_cal, d76, kin_test, schedule, ev_cell, ev_cell * 0, kernels=kern)
    basis = np.ascontiguousarray(np.moveaxis(p, 1, 3))  # (3, nz, K, 4, Nn)
    direct = ens._Kernels(omega=np.full(basis.shape[:4], 0.25), basis=basis)
    d_run = _integrate(film_cal, d76, kin_test, schedule, ev_cell, ev_cell * 0, kernels=direct)
    for m in range(4):
        pm = p[:, m]
        x = (pm * k_run.state.g).sum(axis=-1)
        assert np.max(np.abs(x - d_run.state.g[..., m])) <= 1e-12
        y = (pm * k_run.state.h).sum(axis=-1)
        assert np.max(np.abs(y - d_run.state.h[..., m])) <= 1e-12
    scale = np.maximum(np.abs(k_run.state.conc), 1e-6)
    assert np.max(np.abs(k_run.state.conc - d_run.state.conc) / scale) <= 1e-12
    perm = basis[:, :, :, ::-1]
    d2 = _integrate(
        film_cal,
        d76,
        kin_test,
        schedule,
        ev_cell,
        ev_cell * 0,
        kernels=ens._Kernels(omega=direct.omega, basis=np.ascontiguousarray(perm)),
    )
    assert np.max(np.abs(d2.state.conc - d_run.state.conc) / scale) <= 1e-13


def test_raster_uniform_matches_column(film_image, d76, kin_test):
    # A uniform raster (6 um pixels over 24 um cells, a halo of one cell) gives every inner
    # pixel the column result.
    schedule = make_schedule("image", 60.0)
    nz = film_image.grid.n_z
    ev = np.full((16, 16, nz), 100.0)
    weights = np.zeros((4, 4))
    weights[1:3, 1:3] = 1.0
    members = Members(ev, np.zeros_like(ev), layout="raster", pixel_um=6.0)
    res = run_ensemble(
        film_image,
        d76,
        kin_test,
        schedule,
        members,
        n_cells=16,
        cell_weights=weights.ravel(),
        lateral=Lateral((4, 4), 24.0, "neumann"),
        bulk=_constant_bulk(d76, schedule),
    )
    assert res.silver.shape == (64, 6)  # inner cells 2x2 own 8x8 pixels
    col = run_strip(
        film_image,
        d76,
        kin_test,
        schedule,
        ev[0, :1],
        np.zeros((1, nz)),
        tank=None,
        bulk=_constant_bulk(d76, schedule),
    )
    np.testing.assert_allclose(res.silver, np.broadcast_to(col.silver, res.silver.shape), rtol=1e-9)
    np.testing.assert_allclose(res.halide, col.halide[0], rtol=1e-12)


def test_raster_edge_interpolates_between_cells(film_image, d76, kin_test):
    # A step between two cells: pixel outputs are bilinear in the cell kernels at the pixel
    # centres, so they vary monotonically across the step.
    schedule = make_schedule("image", 60.0)
    nz = film_image.grid.n_z
    ev = np.full((4, 24, nz), 1.0)
    ev[:, 12:] = 1000.0
    members = Members(ev, np.zeros_like(ev), layout="raster", pixel_um=6.0)
    res = run_ensemble(
        film_image,
        d76,
        kin_test,
        schedule,
        members,
        n_cells=6,
        cell_weights=np.ones(6),
        lateral=Lateral((1, 6), 24.0, "neumann"),
        bulk=_constant_bulk(d76, schedule),
    )
    assert res.silver.shape == (96, 6)
    m = developed_silver(res).reshape(4, 24)
    # Within each half the member's own p dominates; across the cells the kernels interpolate.
    assert np.all(m[:, 12:].min() > m[:, :12].max())
    np.testing.assert_array_equal(m[0], m[3])
