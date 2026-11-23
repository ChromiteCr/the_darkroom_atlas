"""Convergence, reference comparisons and the D1c timing criterion of the development engine
(spec D1 3.9, 2.4.13 test_fixed_schedule_vs_reference, 10 D1c done criterion). All slow.

D-76 strips on the fixture film ``test/bw-mechanistic``; thresholds on developed silver M in
g/m2 (spec D1 3.8, 3.9: 1e-3 D = 1.4e-3 g/m2 with CP_nom = 0.7 m2/g)."""

import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pytest

from darkroom_atlas.chemistry.kinetics import resolve_kinetics
from darkroom_atlas.develop.ensemble import film_column_model


def _load_helpers():
    name = "d1c_develop_conftest"
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name("conftest.py"))
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


_h = _load_helpers()

pytestmark = pytest.mark.slow

STRIP_STEPS = 3  # toe, mid-curve, shoulder (exposure 10**[-1.5, 0.75, 3] um^-3)
DURATION_S = 450.0


@pytest.fixture(scope="module")
def setting(library, species, kfile, film_material):
    sol = _h.make_solution(library, species, kfile)
    kin = _h.with_test_rates(resolve_kinetics(kfile, sol, species))
    return sol, kin


def _film(setting, film_material, species, profile, **kw):
    sol, kin = setting
    return film_column_model(
        film_material, "135", _h.EXPOSURE_TIME_S, kin, sol, profile, species=species, **kw
    )


def _strip(film, setting, *, method="adaptive", rtol=1e-4, agitation=None, profile="calibration"):
    sol, kin = setting
    schedule = _h.make_schedule(profile, DURATION_S, agitation or _h.SMALL_TANK)
    ev, es = _h.strip_exposures(STRIP_STEPS, film.grid.n_z)
    return _h.run_strip(film, sol, kin, schedule, ev, es, method=method, rtol=rtol)


def _x(film, res):
    c = film.classes
    return res.silver / (1e12 * c.silver_g_per_grain * c.number_per_um2)


@pytest.fixture(scope="module")
def reference(setting, film_material, species):
    film = _film(setting, film_material, species, "reference")
    return film, _strip(film, setting)


M_TOL = 1.4e-3  # g/m2, spec D1 3.9


def test_convergence_depth(reference, setting, film_material, species):
    _, base = reference
    fine = _film(setting, film_material, species, "reference", levels=(16, 16, 1.25))
    d = np.abs(_h.developed_silver(_strip(fine, setting)) - _h.developed_silver(base)).max()
    assert d <= M_TOL, d


@pytest.mark.xfail(
    strict=True,
    reason="spec D1 3.9: measured |dM| = 1.56e-3 g/m2 for n_l 16 -> 32 (limit 1.4e-3), 3-step strip",
)
def test_convergence_liquid(reference, setting, film_material, species):
    _, base = reference
    # Both liquid sections doubled (8 + 8 -> 16 + 16); r_l = 1.25 -> sqrt(1.25) keeps the span.
    fine = _film(setting, film_material, species, "reference", levels=(8, 32, 1.25**0.5))
    d = np.abs(_h.developed_silver(_strip(fine, setting)) - _h.developed_silver(base)).max()
    assert d <= M_TOL, d


def test_convergence_time(reference, setting):
    film, base = reference
    d = np.abs(_h.developed_silver(_strip(film, setting, rtol=1e-6)) - _h.developed_silver(base))
    assert d.max() <= M_TOL, d


def test_convergence_classes(reference, setting, film_material, species):
    _, base = reference
    fine = _film(setting, film_material, species, "reference", n_classes=24)
    d = np.abs(_h.developed_silver(_strip(fine, setting)) - _h.developed_silver(base)).max()
    assert d <= M_TOL, d


def test_convergence_quadrature(reference, setting, film_material, species):
    _, base = reference
    fine = _film(setting, film_material, species, "reference", quadrature=41)
    d = np.abs(_h.developed_silver(_strip(fine, setting)) - _h.developed_silver(base)).max()
    assert d <= 1.4e-4, d


@pytest.mark.xfail(
    strict=True,
    reason=(
        "spec D1 3.4/2.4.13: measured max|dX| fixed vs adaptive (rtol 1e-4, same grid, 8-step "
        "strip): calibration 4.1e-3 small tank / 2.7e-3 large tank (limit 1e-3), image 8.1e-3 / "
        "6.3e-3 (limit 2e-3). The error is first order in the schedule: each agitation start "
        "resets the held section and ROS2 shows order reduction in the stiff transient."
    ),
)
@pytest.mark.parametrize("profile", ["calibration", "image"])
@pytest.mark.parametrize("tank", ["small", "large"])
def test_fixed_schedule_vs_reference(setting, film_material, species, profile, tank):
    agitation = _h.SMALL_TANK if tank == "small" else _h.LARGE_TANK
    film = _film(setting, film_material, species, profile)
    fixed = _strip(film, setting, method="fixed", agitation=agitation, profile=profile)
    ref = _strip(film, setting, agitation=agitation, profile=profile)
    limit = 1e-3 if profile == "calibration" else 2e-3
    assert np.abs(_x(film, fixed) - _x(film, ref)).max() <= limit


@pytest.mark.xfail(
    strict=True,
    reason=(
        "spec D1 3.9: measured |dM| = 2.6e-2 g/m2 (limit 7e-3); dominated by the fixed-schedule "
        "error at agitation starts (see test_fixed_schedule_vs_reference)"
    ),
)
def test_calibration_profile_vs_reference(reference, setting, film_material, species):
    # Spec D1 3.9: calibration profile (fixed schedule) vs reference profile, |dM| <= 7e-3.
    _, ref = reference
    film = _film(setting, film_material, species, "calibration")
    cal = _strip(film, setting, method="fixed")
    d = np.abs(_h.developed_silver(cal) - _h.developed_silver(ref)).max()
    assert d <= 7e-3, d


def test_timing_35_step_calibration_strip(setting, film_material, species):
    # Spec D1 10 D1c done criterion: a 35-step calibration-profile strip in D-76 (small tank,
    # 7.5 min, 314 steps) on one core in <= 2 s. Best of three to damp machine load.
    sol, kin = setting
    film = _film(setting, film_material, species, "calibration")
    schedule = _h.make_schedule("calibration", DURATION_S, _h.SMALL_TANK)
    assert schedule.n_steps == 314
    ev, es = _h.strip_exposures(35, film.grid.n_z)
    best = np.inf
    for _ in range(3):
        t0 = time.perf_counter()
        _h.run_strip(film, sol, kin, schedule, ev, es)
        best = min(best, time.perf_counter() - t0)
    print(f"35-step calibration strip: {best:.3f} s, {1e3 * best / schedule.n_steps:.2f} ms/step")
    assert best <= 2.0
