"""Speciation and pH (spec D1 2.3.4, 2.3.8)."""

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from darkroom_atlas.chemistry.activity import debye_huckel_d
from darkroom_atlas.chemistry.speciation import (
    Speciation,
    buffer_capacity,
    fractions,
    make_speciation,
    mean_protons,
    solve_ph,
    water_ions,
)
from darkroom_atlas.chemistry.species import AcidBaseFamily, component_table
from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import ActivityParams

T25 = 298.15
T20 = 293.15


def _q(value, unit="1"):
    return {"value": value, "unit": unit, "source": "assumed", "note": "test"}


def params(a20=0.5050, a25=0.5091, **eps) -> ActivityParams:
    # Spec D1 2.3.3, 2.3.6.
    return ActivityParams(
        debye_huckel_A_20C=_q(a20),
        debye_huckel_A_25C=_q(a25),
        sit_b=_q(1.5),
        delta_eps={k: _q(v, "kg/mol") for k, v in eps.items()},
    )


IDEAL = params(0.0, 0.0)


def weak_acid(pka: float, name="ha") -> AcidBaseFamily:
    """HA / A-, reference state A- (so TOTH counts the protons on HA)."""
    return AcidBaseFamily(name, (0, -1), (1, 0), (pka,), "van_t_hoff", (0.0,), (), "organic")


@pytest.fixture(scope="module")
def table():
    return component_table(Library().species())


@pytest.fixture(scope="module")
def d76(table):
    """Families and totals of the published D-76 (2/100/5/2 g/L), assembled by hand."""
    fams = [table[k].family for k in ("sulfite", "borate", "hydroquinone", "metol")]
    metol = 2.0 / 344.39 * 2
    borax = 2.0 / 381.37
    totals = np.array([100.0 / 126.04, 4 * borax, 5.0 / 110.11, metol])
    toth = 2 * 2.0 / 344.39 - 2 * borax
    return fams, totals, toth


def _proton_excess(spec, totals, ph):
    """TOTH recomputed from the state fractions (independent of the solver's moments)."""
    alpha = fractions(spec, ph)
    h, oh = water_ions(spec, ph)
    out = h - oh
    for f, fam in enumerate(spec.families):
        n = np.asarray(fam.protons_rel, dtype=float)
        out = out + totals[f] * np.tensordot(n, alpha[f, : n.size], axes=(0, 0))
    return out


# ---------------------------------------------------------------------------------------
# Spec D1 2.3.8
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("pka", "c_acid", "c_base"),
    [(4.76, 0.05, 0.1), (7.0, 0.02, 0.02), (9.5, 1e-3, 3e-3), (11.0, 0.2, 0.01)],
)
def test_speciation_ideal_limit(pka, c_acid, c_base):
    """I -> 0, d_eps = 0: pH equals Henderson-Hasselbalch within 1e-10.

    The water ions are kept in the exact form pH = pK + log10((C_b + [H] - [OH]) /
    (C_a - [H] + [OH])), which reduces to the textbook form when they are negligible.
    """
    spec = make_speciation([weak_acid(pka)], T25, 0.0, 0.0, params(sulfite=0.5))
    ph = solve_ph(spec, np.array([[c_acid + c_base]]), np.array([c_acid]))[0]
    h, oh = 10.0**-ph, 10.0 ** (ph - 14.0)
    hh = pka + math.log10((c_base + h - oh) / (c_acid - h + oh))
    assert abs(ph - hh) <= 1e-10


def test_speciation_ideal_limit_textbook_form():
    """At pK = pK_w / 2 and a 1:1 ratio the water terms cancel: pH = pK exactly."""
    spec = make_speciation([weak_acid(7.0)], T25, 0.0, 0.0, IDEAL)
    ph = solve_ph(spec, np.array([[0.02]]), np.array([0.01]))[0]
    assert abs(ph - 7.0) <= 1e-10


def test_buffer_titration(d76):
    """pH after adding delta of TOTH equals the integral of -dTOTH / beta within 1e-9."""
    fams, totals, toth = d76
    spec = make_speciation(fams, T20, 2.39, 1.6, params(sulfite=-0.1, borate=0.05))
    t = totals[:, None]
    ph_a = solve_ph(spec, t, np.array([toth]))[0]
    for delta in (5e-3, -5e-3, 2e-2):
        ph_b = solve_ph(spec, t, np.array([toth + delta]))[0]

        def rhs(_, y):
            return [-1.0 / buffer_capacity(spec, t, np.array([y[0]]))[0]]

        sol = solve_ivp(rhs, (0.0, delta), [ph_a], method="DOP853", rtol=1e-13, atol=1e-14)
        assert abs(sol.y[0, -1] - ph_b) <= 1e-9


def test_toth_round_trip(d76):
    """TOTH recomputed from the speciation equals the input within 1e-13 (relative)."""
    fams, totals, _ = d76
    spec = make_speciation(fams, T25, 2.39, 1.6, params())
    toth = np.array([0.05, 0.2, -0.03, -0.1, 0.5])
    t = np.repeat(totals[:, None], toth.size, axis=1)
    ph = solve_ph(spec, t, toth)
    back = _proton_excess(spec, t, ph)
    np.testing.assert_array_less(np.abs(back - toth) / np.abs(toth), 1e-13)


# ---------------------------------------------------------------------------------------
# Balances
# ---------------------------------------------------------------------------------------


def test_charge_balance_residual(table):
    """A solution made of electroneutral salts is electroneutral at the solved pH."""
    fams = [table[k].family for k in ("sulfite", "borate", "diethanolamine", "dtpa")]
    na2so3, nahso3, borax, dea, dtpa_na5 = 0.4, 0.05, 0.01, 0.1, 0.002
    totals = np.array([na2so3 + nahso3, 4 * borax, dea, dtpa_na5])
    toth = nahso3 - 2 * borax - 5 * dtpa_na5
    sodium = 2 * na2so3 + nahso3 + 2 * borax + 5 * dtpa_na5
    spec = make_speciation(fams, T25, 1.4, sodium, params(sulfite=0.1, amine=-0.05))
    ph = solve_ph(spec, totals[:, None], np.array([toth]))
    alpha = fractions(spec, ph)[:, :, 0]
    h, oh = water_ions(spec, ph)
    charge = sodium + h[0] - oh[0]
    for f, fam in enumerate(fams):
        z = np.asarray(fam.charges, dtype=float)
        charge += totals[f] * np.dot(z, alpha[f, : z.size])
    assert abs(charge) <= 1e-14 * sodium


def test_mass_balance(d76):
    fams, totals, _ = d76
    spec = make_speciation(fams, T25, 1.0, 1.0, params())
    ph = np.linspace(0.0, 14.5, 59)
    alpha = fractions(spec, ph)
    assert alpha.shape == (4, 3, 59)
    np.testing.assert_allclose(alpha.sum(axis=1), 1.0, rtol=0, atol=1e-15)
    assert np.all(alpha >= 0.0)
    # Borate has 2 states: its padded third state is exactly 0.
    assert np.all(alpha[1, 2] == 0.0)
    species_amounts = totals[:, None, None] * alpha
    np.testing.assert_allclose(species_amounts.sum(axis=1), totals[:, None] + 0 * ph, rtol=1e-15)


def test_fractions_follow_conditional_pka(table):
    spec = make_speciation([table["sulfite"].family], T25, 2.0, 1.0, params(sulfite=0.2))
    pk2 = spec.pka[0, 1]
    alpha = fractions(spec, np.array([pk2]))[0, :, 0]
    assert alpha[1] == pytest.approx(alpha[2], rel=1e-13)


# ---------------------------------------------------------------------------------------
# Known buffers
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(("temperature", "ph_ref"), [(T25, 9.180), (T20, 9.225)])
def test_borax_standard_buffer(table, temperature, ph_ref):
    """0.01 mol/kg borax against the NIST pH(S) values 9.180 (25 C) and 9.225 (20 C).

    Reference values recalled from the NIST borax standard (SRM 187; Bates, Determination of
    pH, 1973); the SIT D(I) term at I = 0.02 should land within 0.01 of them.
    """
    borax = 0.01
    totals = np.array([[4 * borax]])
    spec = make_speciation([table["borate"].family], temperature, 0.02, 0.02, params())
    ph = solve_ph(spec, totals, np.array([-2 * borax]))[0]
    assert ph == pytest.approx(ph_ref, abs=0.01)


def test_pure_water_is_neutral():
    for t, ph_ref in ((T25, 7.0), (T20, 7.085)):
        spec = make_speciation([], t, 0.0, 0.0, IDEAL)
        ph = solve_ph(spec, np.zeros((0, 1)), np.zeros(1))[0]
        assert ph == pytest.approx(ph_ref, abs=1e-12)


def test_strong_base_and_strong_acid():
    spec = make_speciation([], T25, 0.0, 0.0, IDEAL)
    ph = solve_ph(spec, np.zeros((0, 2)), np.array([-0.1, 0.1]))
    np.testing.assert_allclose(ph, [13.0, 1.0], atol=1e-6)


def test_half_neutralised_amine(table):
    """Equimolar DEA / DEAH+: pH = pK'(T) up to the water ions."""
    dea = table["diethanolamine"].family
    for t in (T20, T25):
        spec = make_speciation([dea], t, 0.5, 0.0, params())
        ph = solve_ph(spec, np.array([[0.2]]), np.array([0.1]))[0]
        assert ph == pytest.approx(spec.pka[0, 0], abs=5e-4)


# ---------------------------------------------------------------------------------------
# Buffer capacity and derivatives
# ---------------------------------------------------------------------------------------


def test_buffer_capacity_matches_numerical_derivative(d76):
    fams, totals, _ = d76
    spec = make_speciation(fams, T25, 2.39, 1.6, params(sulfite=0.1))
    ph = np.linspace(2.0, 13.0, 23)
    t = np.repeat(totals[:, None], ph.size, axis=1)
    h = 1e-5
    numeric = -(_proton_excess(spec, t, ph + h) - _proton_excess(spec, t, ph - h)) / (2 * h)
    beta = buffer_capacity(spec, t, ph)
    np.testing.assert_allclose(beta, numeric, rtol=1e-7)
    assert np.all(beta > 0.0)


def test_ph_derivatives(d76):
    """dpH/dTOTH = -1/beta and dpH/dT_f = nbar_f / beta (spec D1 2.3.4)."""
    fams, totals, toth = d76
    spec = make_speciation(fams, T25, 2.39, 1.6, params())
    t = totals[:, None]
    ph = solve_ph(spec, t, np.array([toth]))
    beta = buffer_capacity(spec, t, ph)[0]
    nbar = mean_protons(spec, ph)[:, 0]
    d = 1e-7
    up = solve_ph(spec, t, np.array([toth + d]))[0]
    dn = solve_ph(spec, t, np.array([toth - d]))[0]
    assert (up - dn) / (2 * d) == pytest.approx(-1.0 / beta, rel=1e-6)
    for f in range(len(fams)):
        tp, tm = t.copy(), t.copy()
        tp[f] += d
        tm[f] -= d
        slope = solve_ph(spec, tp, np.array([toth]))[0] - solve_ph(spec, tm, np.array([toth]))[0]
        assert slope / (2 * d) == pytest.approx(nbar[f] / beta, rel=1e-5, abs=1e-9)


# ---------------------------------------------------------------------------------------
# Temperature and ionic strength
# ---------------------------------------------------------------------------------------


def test_temperature_dependence(d76, table):
    fams, totals, toth = d76
    t = totals[:, None]
    ph20 = solve_ph(make_speciation(fams, T20, 0.0, 0.0, IDEAL), t, np.array([toth]))[0]
    ph25 = solve_ph(make_speciation(fams, T25, 0.0, 0.0, IDEAL), t, np.array([toth]))[0]
    assert ph20 > ph25  # borate pK and pK_w rise on cooling
    dea = table["diethanolamine"].family
    s20 = make_speciation([dea], T20, 0.0, 0.0, IDEAL)
    s25 = make_speciation([dea], T25, 0.0, 0.0, IDEAL)
    p20 = solve_ph(s20, np.array([[0.2]]), np.array([0.1]))[0]
    p25 = solve_ph(s25, np.array([[0.2]]), np.array([0.1]))[0]
    assert p20 - p25 == pytest.approx(0.125, abs=1e-4)


def test_ionic_strength_lowers_sulfite_buffer_ph(table):
    """At high I the bisulfite pK' falls by 3D, so a sulfite buffer reads lower."""
    sulfite = table["sulfite"].family
    t, q = np.array([[0.8]]), np.array([0.4])
    ideal = solve_ph(make_speciation([sulfite], T25, 0.0, 0.0, IDEAL), t, q)[0]
    p = params()
    real = solve_ph(make_speciation([sulfite], T25, 2.0, 1.6, p), t, q)[0]
    assert real - ideal == pytest.approx(-3 * debye_huckel_d(2.0, T25, p), abs=1e-3)


# ---------------------------------------------------------------------------------------
# Solver behaviour
# ---------------------------------------------------------------------------------------


def _random_cells(table, n, seed=7):
    rng = np.random.default_rng(seed)
    ids = ("sulfite", "borate", "hydroquinone", "metol", "diethanolamine", "dtpa", "ascorbate")
    fams = [table[k].family for k in ids]
    totals = rng.uniform(0.0, 0.3, size=(len(fams), n)) * (rng.uniform(size=(len(fams), n)) < 0.8)
    toth = rng.uniform(-0.1, 0.4, size=n)
    return fams, totals, toth


def test_vectorised_equals_scalar(table):
    fams, totals, toth = _random_cells(table, 64)
    spec = make_speciation(fams, T20, 1.2, 0.9, params(sulfite=0.1, organic=-0.05))
    batch = solve_ph(spec, totals, toth)
    single = np.array(
        [solve_ph(spec, totals[:, i : i + 1], toth[i : i + 1])[0] for i in range(toth.size)]
    )
    np.testing.assert_array_equal(batch, single)
    shaped = solve_ph(spec, totals.reshape(len(fams), 8, 8), toth.reshape(8, 8))
    np.testing.assert_array_equal(shaped.ravel(), batch)
    np.testing.assert_array_equal(
        buffer_capacity(spec, totals, batch)[3],
        buffer_capacity(spec, totals[:, 3:4], batch[3:4])[0],
    )


@pytest.mark.parametrize("ph0", [0.0, 14.5, 7.0, 3.3, 12.9])
def test_convergence_from_poor_starting_values(table, ph0):
    fams, totals, toth = _random_cells(table, 200, seed=11)
    spec = make_speciation(fams, T25, 0.8, 0.5, params())
    ref = solve_ph(spec, totals, toth)
    got = solve_ph(spec, totals, toth, ph0=np.full(toth.size, ph0), max_iter=30)
    np.testing.assert_allclose(got, ref, rtol=0, atol=1e-11)


@pytest.mark.parametrize("ph0", [0.0, 14.5])
def test_convergence_from_bracket_ends_for_extreme_solutions(ph0):
    spec = make_speciation([weak_acid(4.0), weak_acid(10.0, "hb")], T25, 0.0, 0.0, IDEAL)
    totals = np.array([[1e-4, 0.0, 0.5, 0.0], [0.0, 0.5, 0.5, 0.0]])
    toth = np.array([-0.2, 0.3, 0.9, 0.0])  # strong base, strong acid, acids, water
    got = solve_ph(spec, totals, toth, ph0=np.full(4, ph0))
    ref = solve_ph(spec, totals, toth)
    np.testing.assert_allclose(got, ref, atol=1e-11)
    assert got[3] == pytest.approx(7.0, abs=1e-12)


def test_root_outside_bracket_raises():
    spec = make_speciation([], T25, 0.0, 0.0, IDEAL)
    with pytest.raises(ValueError, match="above 14.5"):
        solve_ph(spec, np.zeros((0, 2)), np.array([0.0, -5.0]))
    with pytest.raises(ValueError, match="below 0"):
        solve_ph(spec, np.zeros((0, 1)), np.array([5.0]))


def test_no_convergence_raises(d76):
    fams, totals, toth = d76
    spec = make_speciation(fams, T25, 0.0, 0.0, IDEAL)
    with pytest.raises(RuntimeError, match="did not converge"):
        solve_ph(spec, totals[:, None], np.array([toth]), max_iter=2)


def test_input_validation(d76):
    fams, totals, toth = d76
    spec = make_speciation(fams, T25, 0.0, 0.0, IDEAL)
    with pytest.raises(ValueError, match="shape"):
        solve_ph(spec, totals[:3, None], np.array([toth]))
    with pytest.raises(ValueError, match=">= 0"):
        solve_ph(spec, -totals[:, None], np.array([toth]))
    with pytest.raises(ValueError, match="toth shape"):
        solve_ph(spec, totals[:, None], np.array([toth, toth, toth]))
    with pytest.raises(ValueError, match="finite"):
        solve_ph(spec, totals[:, None], np.array([np.nan]))
    with pytest.raises(ValueError, match="tol"):
        solve_ph(spec, totals[:, None], np.array([toth]), tol=0.0)


def test_speciation_shape_and_padding(table):
    fams = [table["borate"].family, table["dtpa"].family]
    spec = make_speciation(fams, T25, 0.5, 0.5, params())
    assert spec.pka.shape == (2, 5)
    assert np.isnan(spec.pka[0, 1:]).all() and np.isfinite(spec.pka[1]).all()
    assert not spec.pka.flags.writeable
    with pytest.raises(ValueError, match="shape"):
        Speciation(tuple(fams), np.zeros((2, 4)), -14.0, 0.0, 0.0)
    with pytest.raises(ValueError, match="NaN padding"):
        Speciation(tuple(fams), np.zeros((2, 5)), -14.0, 0.0, 0.0)
