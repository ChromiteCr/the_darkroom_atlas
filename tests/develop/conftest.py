"""Shared fixtures of the development-engine tests (spec D1 2.4.13, 3.9).

The film is the fixture ``test/bw-mechanistic`` (spec D1 3.9: D1c uses its parameters); the
developer is the packaged ``d-76`` of the data library at 20 C.
"""

import dataclasses
from pathlib import Path

import numpy as np
import pytest

from darkroom_atlas.chemistry.kinetics import resolve_kinetics
from darkroom_atlas.chemistry.solution import working_solution
from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import Agitation
from darkroom_atlas.data.provenance import Quantity
from darkroom_atlas.develop.ensemble import Members, TankSpec, film_column_model, run_ensemble
from darkroom_atlas.develop.schedule import agitation_periods, step_schedule

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "library"
FILM_ID = "test/bw-mechanistic"
EXPOSURE_TIME_S = 0.01

# Spec D1 3.1 fixed schedules (h0, ratio, h_max).
SCHEDULES = {"image": (0.1, 1.6, 5.0), "calibration": (0.05, 1.4, 3.0)}
# Spec D1 3.4 table conditions (initial_s = 15 s reproduces its step counts, see test_schedule).
SMALL_TANK = Agitation(method="inversion", initial_s=15.0, interval_s=30.0, duration_s=5.0)
LARGE_TANK = Agitation(method="lift", initial_s=22.5, interval_s=60.0, duration_s=7.5)

# A 135-36 roll (35 mm x 1.6 m, both sides of the reel wetted on one side only) in 0.3 L: a
# typical one-shot small-tank load; it only sets the depletion scale.
ROLL_AREA_UM2 = 0.035 * 1.6 * 1e12
ROLL_TANK = TankSpec(volume_l=0.3, film_area_um2=ROLL_AREA_UM2)
INFINITE_TANK = TankSpec(volume_l=1e9, film_area_um2=1e12)


@pytest.fixture(scope="session")
def library():
    return Library()


@pytest.fixture(scope="session")
def species(library):
    return library.species()


@pytest.fixture(scope="session")
def kfile(library):
    return library.kinetics()


@pytest.fixture(scope="session")
def film_material():
    return Library(FIXTURES).material(FILM_ID)


def make_solution(library, species, kfile, bath="d-76", dilution="stock", t=20.0):
    return working_solution(library.bath(bath), species, kfile, dilution=dilution, temperature_c=t)


@pytest.fixture(scope="session")
def d76(library, species, kfile):
    return make_solution(library, species, kfile)


@pytest.fixture(scope="session")
def kin_d76(kfile, d76, species):
    return resolve_kinetics(kfile, d76, species)


@pytest.fixture(scope="session")
def film_cal(film_material, kin_d76, d76, species):
    return film_column_model(
        film_material, "135", EXPOSURE_TIME_S, kin_d76, d76, "calibration", species=species
    )


@pytest.fixture(scope="session")
def film_image(film_material, kin_d76, d76, species):
    return film_column_model(
        film_material, "135", EXPOSURE_TIME_S, kin_d76, d76, "image", species=species
    )


def make_schedule(profile, duration_s, agitation=SMALL_TANK, record_s=()):
    h0, ratio, h_max = SCHEDULES.get(profile, SCHEDULES["calibration"])
    periods = agitation_periods(agitation, duration_s)
    return step_schedule(
        periods, duration_s, h0_s=h0, ratio=ratio, h_max_s=h_max, record_s=record_s
    )


def strip_exposures(n_steps, n_z, lo=-1.5, hi=3.0):
    """Exposure steps: ``e_V`` log-spaced in ``[10**lo, 10**hi]`` um^-3, uniform in depth,
    ``e_S = 0`` (a grain of volume V absorbs ``m = V e_V`` photons)."""
    ev = np.repeat(np.logspace(lo, hi, n_steps)[:, None], n_z, axis=1)
    return ev, np.zeros_like(ev)


def columns(ev, es):
    return Members(exposure_volume=ev, exposure_surface=es, layout="columns")


def developed_silver(result):
    """Total developed silver M = sum_c silver_c + silver_physical, g/m2 (spec D1 3.8)."""
    return result.silver.sum(axis=1) + result.silver_physical


def run_strip(film, solution, kin, schedule, ev, es, *, tank=ROLL_TANK, **kw):
    n = ev.shape[0]
    return run_ensemble(
        film,
        solution,
        kin,
        schedule,
        columns(ev, es),
        n_cells=n,
        cell_weights=np.ones(n),
        tank=tank,
        **kw,
    )


# Kinetics for the engine tests: the kinetics.yaml values are uncalibrated placeholders under
# which D-76 hardly develops in 7.5 min. These values sit inside the prior bounds of spec D1
# 2.4.11 and give a characteristic curve of about 0.1-1.8 g/m2 on the fixture film.
TEST_RATES = {"k_et": 0.02, "lambda0": 10.0, "K_ads": 1000.0}


def with_test_rates(kin):
    a = len(kin.agent_ids)
    return dataclasses.replace(kin, **{k: np.full(a, v) for k, v in TEST_RATES.items()})


@pytest.fixture(scope="session")
def kin_test(kin_d76):
    return with_test_rates(kin_d76)


def with_solvent(kin, **values):
    """Copy of ``kin`` with solvent constants replaced (``name=(value, unit)``)."""
    update = {k: Quantity(value=v, unit=u, source="assumed") for k, (v, u) in values.items()}
    return dataclasses.replace(kin, solvent=kin.solvent.model_copy(update=update))


def with_transport(kin, **values):
    update = {k: Quantity(value=v, unit=u, source="assumed") for k, (v, u) in values.items()}
    return dataclasses.replace(kin, transport=kin.transport.model_copy(update=update))
