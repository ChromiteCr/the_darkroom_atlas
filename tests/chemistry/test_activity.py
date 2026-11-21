"""SIT-type activity corrections and water constants (spec D1 2.3.3, 2.3.6)."""

import itertools
import math

import numpy as np
import pytest

from darkroom_atlas.chemistry.activity import (
    WATER_IONISATION_ENTHALPY_KJ_MOL,
    conditional_pka,
    debye_huckel_a,
    debye_huckel_d,
    delta_eps,
    pkw,
    water_constants,
)
from darkroom_atlas.chemistry.species import AcidBaseFamily
from darkroom_atlas.data.models import ActivityParams


def _q(value, unit="1"):
    return {"value": value, "unit": unit, "source": "assumed", "note": "test"}


def params(**eps) -> ActivityParams:
    # Spec D1 2.3.3, 2.3.6: A(20 C) 0.5050, A(25 C) 0.5091, b = 1.5.
    return ActivityParams(
        debye_huckel_A_20C=_q(0.5050),
        debye_huckel_A_25C=_q(0.5091),
        sit_b=_q(1.5),
        delta_eps={k: _q(v, "kg/mol") for k, v in eps.items()},
    )


def family(charge0: int, pka: tuple[float, ...], cls: str) -> AcidBaseFamily:
    n = len(pka)
    return AcidBaseFamily(
        component="x",
        charges=tuple(charge0 - j for j in range(n + 1)),
        protons_rel=tuple(n - j for j in range(n + 1)),
        pka25=pka,
        temperature_law="van_t_hoff",
        enthalpy_kj_mol=(0.0,) * n,
        dpka_dt=(),
        activity_class=cls,
    )


SULFITE = family(0, (1.86, 7.19), "sulfite")
BORATE = family(0, (9.235,), "borate")
DEA = family(1, (8.883,), "amine")


def test_debye_huckel_a_interpolates_and_extrapolates():
    p = params()
    assert debye_huckel_a(293.15, p) == pytest.approx(0.5050, abs=1e-15)
    assert debye_huckel_a(298.15, p) == pytest.approx(0.5091, abs=1e-15)
    assert debye_huckel_a(303.15, p) == pytest.approx(0.5132, abs=1e-12)


def test_debye_huckel_d_at_d76_ionic_strength():
    # Spec D1 12.4: D(I) = 0.237 at I = 2.39 M (25 C).
    d = debye_huckel_d(2.39, 298.15, params())
    assert d == pytest.approx(0.237, abs=5e-4)
    root = math.sqrt(2.39)
    assert d == pytest.approx(0.5091 * root / (1 + 1.5 * root), rel=1e-14)


def test_debye_huckel_d_zero_and_monotone():
    p = params()
    assert debye_huckel_d(0.0, 298.15, p) == 0.0
    values = [debye_huckel_d(i, 298.15, p) for i in (0.01, 0.1, 1.0, 3.0)]
    assert all(b > a for a, b in itertools.pairwise(values))


@pytest.mark.parametrize("ionic_strength", [-0.1, float("nan"), float("inf")])
def test_debye_huckel_d_rejects_bad_ionic_strength(ionic_strength):
    with pytest.raises(ValueError, match="ionic strength"):
        debye_huckel_d(ionic_strength, 298.15, params())


def test_conditional_pka_charge_factors():
    """(2 z_HA - 1) D: -3D for HSO3-/SO3 2-, -D for boric acid, +D for DEAH+ (spec D1 2.3.3)."""
    p = params()
    i, t = 2.39, 298.15
    d = debye_huckel_d(i, t, p)
    np.testing.assert_allclose(
        conditional_pka(SULFITE, t, i, 0.0, p), [1.86 - d, 7.19 - 3 * d], rtol=0, atol=1e-14
    )
    np.testing.assert_allclose(conditional_pka(BORATE, t, i, 0.0, p), [9.235 - d], atol=1e-14)
    np.testing.assert_allclose(conditional_pka(DEA, t, i, 0.0, p), [8.883 + d], atol=1e-14)


def test_conditional_pka_sit_term():
    p = params(sulfite=0.2, amine=-0.1)
    d = debye_huckel_d(1.0, 298.15, p)
    m = 1.6
    np.testing.assert_allclose(
        conditional_pka(SULFITE, 298.15, 1.0, m, p), [1.86 - d + 0.32, 7.19 - 3 * d + 0.32]
    )
    np.testing.assert_allclose(conditional_pka(DEA, 298.15, 1.0, m, p), [8.883 + d - 0.16])
    # A class without a coefficient, and 'none', get 0.
    np.testing.assert_allclose(conditional_pka(BORATE, 298.15, 1.0, m, p), [9.235 - d])
    assert delta_eps("none", p) == 0.0
    assert delta_eps("organic", p) == 0.0


def test_conditional_pka_ideal_limit_is_thermodynamic():
    p = params(sulfite=0.3)
    np.testing.assert_array_equal(
        conditional_pka(SULFITE, 293.15, 0.0, 0.0, p), SULFITE.pka(293.15)
    )


def test_conditional_pka_rejects_negative_molality():
    with pytest.raises(ValueError, match="molality"):
        conditional_pka(SULFITE, 298.15, 1.0, -1.0, params())


def test_pkw_through_both_points():
    assert pkw(293.15) == pytest.approx(14.17, abs=1e-12)
    assert pkw(298.15) == pytest.approx(14.00, abs=1e-12)
    # Spec D1 2.3.3: Delta H_w = 56.9 kJ/mol, derived from the two values.
    assert WATER_IONISATION_ENTHALPY_KJ_MOL == pytest.approx(56.9, abs=0.05)
    assert pkw(303.15) < 14.0


def test_water_constants():
    p = params()
    log_kw, log_gh, log_goh = water_constants(298.15, 0.0, p)
    assert (log_kw, log_gh, log_goh) == (pytest.approx(-14.0, abs=1e-12), 0.0, 0.0)
    d = debye_huckel_d(2.0, 293.15, p)
    log_kw, log_gh, log_goh = water_constants(293.15, 2.0, p)
    assert log_kw == pytest.approx(-14.17, abs=1e-12)
    assert log_gh == pytest.approx(-d, abs=1e-15)
    assert log_goh == pytest.approx(-d, abs=1e-15)
