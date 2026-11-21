"""Working solutions from bath formulas (spec D1 2.3.1, 2.3.2, 2.3.8, 12.4)."""

import numpy as np
import pytest

from darkroom_atlas.chemistry.solution import (
    WorkingSolution,
    component_charges,
    working_solution,
)
from darkroom_atlas.chemistry.speciation import fractions, solve_ph, water_ions
from darkroom_atlas.chemistry.species import BUILTIN_COMPONENTS, component_table
from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import (
    ActivityParams,
    AgentKinetics,
    Bath,
    GelatinParams,
    KineticsFile,
)
from darkroom_atlas.data.provenance import Quantity


def _q(value, unit="1"):
    return {"value": value, "unit": unit, "source": "assumed", "note": "test"}


def make_kinetics(a20=0.5050, a25=0.5091, pyrazolidone_pka=None, **eps) -> KineticsFile:
    """Kinetics with only the blocks a working solution reads (activity, agents, gelatin).

    The activity block is validated; the rest of the file is not needed here, so it is
    built with ``model_construct`` to keep these tests independent of the other blocks.
    A and b are the spec D1 2.3.6 / 5.8 values.
    """
    activity = ActivityParams.model_validate(
        {
            "debye_huckel_A_20C": _q(a20),
            "debye_huckel_A_25C": _q(a25),
            "sit_b": _q(1.5),
            "delta_eps": {k: _q(v, "kg/mol") for k, v in eps.items()},
        }
    )
    agents = {}
    if pyrazolidone_pka is not None:
        agents["pyrazolidone"] = AgentKinetics.model_construct(
            pka=Quantity.model_validate(_q(pyrazolidone_pka))
        )
    gelatin = GelatinParams.model_construct(pka=Quantity.model_validate(_q(10.5)))
    return KineticsFile.model_construct(activity=activity, agents=agents, gelatin=gelatin)


SIT = make_kinetics()  # D(I) only, delta eps = 0
IDEAL = make_kinetics(0.0, 0.0)  # D = 0, delta eps = 0

# XTOL Part A in 850 mL (US 5,756,271 Example 1; DEV 2.1), as overrides of the xtol file.
XTOL_PART_A = {
    "sodium-sulfite": 10.0 / 0.85,
    "dtpa-pentasodium": 1.0 / 0.85,
    "sodium-metaborate-tetrahydrate": 4.0 / 0.85,
    "dimezone-s": 0.2 / 0.85,
    "sodium-metabisulfite": 0.0,
    "sodium-erythorbate-monohydrate": 0.0,
}

# Spec D1 12.4 computes the XTOL values with Dimezone S non-ionising (as DEV 5.1, which had
# no pKa for it). Its only effect on pH is through its acid step, so an amount of 0 is the
# same reading; species.yaml carries a placeholder pKa 8.5 that lowers the pH (spec issue).
NO_DIMEZONE_ACID = {"dimezone-s": 0.0}


@pytest.fixture(scope="module")
def library():
    return Library()


@pytest.fixture(scope="module")
def species(library):
    return library.species()


def ws(library, species, bath_id, kinetics=SIT, dilution="stock", t=25.0, **overrides):
    return working_solution(
        library.bath(bath_id),
        species,
        kinetics,
        dilution=dilution,
        temperature_c=t,
        overrides=overrides or None,
    )


def _with_component(bath: Bath, species_id: str, grams_per_l: float) -> Bath:
    comps = [c.model_dump(mode="json") for c in bath.components]
    comps.append({"species": species_id, "amount": _q(grams_per_l, "g/L")})
    return Bath.model_validate({**bath.model_dump(mode="json"), "components": comps})


# ---------------------------------------------------------------------------------------
# Spec D1 2.3.8 / 12.4 pH anchors
# ---------------------------------------------------------------------------------------


def test_d76_published_ideal_ph(library, species):
    sol = ws(library, species, "d-76-published", IDEAL)
    assert sol.ph == pytest.approx(8.95, abs=0.02)  # spec 12.4; DEV 5.1 gives 8.94


def test_xtol_ideal_ph(library, species):
    """Spec 12.4 values assume a non-ionising Dimezone S (DEV 5.1); its amount is set to 0,
    which is equivalent for pH."""
    work = ws(library, species, "xtol", IDEAL, **NO_DIMEZONE_ACID)
    assert work.ph == pytest.approx(8.89, abs=0.02)
    part_a = ws(library, species, "xtol", IDEAL, **{**XTOL_PART_A, **NO_DIMEZONE_ACID})
    assert part_a.ph == pytest.approx(11.04, abs=0.03)


def test_sit_debye_huckel_only(library, species):
    """Delta eps = 0 with D(I). XTOL as in test_xtol_ideal_ph (non-ionising Dimezone S)."""
    d76 = ws(library, species, "d-76-published", SIT)
    assert d76.ph == pytest.approx(8.44, abs=0.02)
    xtol = ws(library, species, "xtol", SIT, **NO_DIMEZONE_ACID)
    assert xtol.ph == pytest.approx(8.29, abs=0.02)


def test_ionic_strength_d76(library, species):
    sol = ws(library, species, "d-76-published", SIT)
    assert sol.ionic_strength == pytest.approx(2.39, abs=0.03)


@pytest.mark.parametrize(
    ("b2o3", "ideal", "sit"),
    [(1.0, 8.82, 8.32), (2.0, 8.73, 8.23), (3.0, 8.67, 8.17)],
)
def test_d76_published_plus_boric_oxide(library, species, b2o3, ideal, sit):
    """Spec 12.4 table row 'published + B2O3 1 / 2 / 3 g/L'."""
    bath = _with_component(library.bath("d-76-published"), "boric-oxide", b2o3)
    for kin, ref in ((IDEAL, ideal), (SIT, sit)):
        sol = working_solution(bath, species, kin, dilution="stock", temperature_c=25.0)
        assert sol.ph == pytest.approx(ref, abs=0.02)


def test_xtol_part_a_sit(library, species):
    """Spec 12.4: Part A with D(I) 10.76 at I = 0.32 M (non-ionising Dimezone S)."""
    part_a = ws(library, species, "xtol", SIT, **{**XTOL_PART_A, **NO_DIMEZONE_ACID})
    assert part_a.ph == pytest.approx(10.76, abs=0.03)
    assert part_a.ionic_strength == pytest.approx(0.32, abs=0.03)


def test_placeholder_dimezone_pka_lowers_xtol_ph(library, species):
    """With the species-file pKa (8.5 placeholder) the acid step lowers the pH."""
    with_acid = ws(library, species, "xtol", SIT)
    without = ws(library, species, "xtol", SIT, **NO_DIMEZONE_ACID)
    assert with_acid.ph < without.ph


# ---------------------------------------------------------------------------------------
# Dilution (spec D1 2.3.2, 2.3.8; DEV 5.2)
# ---------------------------------------------------------------------------------------


def test_dilution_raises_ph_in_sulfite_buffer(library, species):
    stock = ws(library, species, "xtol", SIT)
    diluted = ws(library, species, "xtol", SIT, dilution="1+1")
    assert diluted.ph > stock.ph  # DEV 5.2 direction
    assert diluted.ionic_strength < stock.ionic_strength


def test_ideal_ph_unchanged_by_dilution(library, species):
    """DEV 5.2: the ideal model keeps the buffer ratio, so pH does not change."""
    for bath_id in ("d-76-published", "xtol"):
        stock = ws(library, species, bath_id, IDEAL)
        for dilution in ("1+1", "1+3"):
            assert ws(library, species, bath_id, IDEAL, dilution=dilution).ph == pytest.approx(
                stock.ph, abs=1e-3
            )


def test_dilution_scales_totals(library, species):
    stock = ws(library, species, "xtol", SIT)
    for dilution, f in (("1+1", 0.5), ("1+2", 1 / 3), ("1+3", 0.25)):
        sol = ws(library, species, "xtol", SIT, dilution=dilution)
        assert sol.dilution == dilution
        assert sol.dilution_fraction == pytest.approx(f, rel=1e-15)
        np.testing.assert_allclose(sol.totals, f * stock.totals, rtol=1e-14, atol=0)
        assert sol.toth == pytest.approx(f * stock.toth, rel=1e-14)
        assert sol.cation_molality == pytest.approx(f * stock.cation_molality, rel=1e-14)


def test_concentrate_dilution(library, species):
    """HC-110 B is 1+31 from the concentrate (spec D1 5.5): f = 1/32."""
    conc = ws(library, species, "hc-110-2019", dilution="stock")
    b = ws(library, species, "hc-110-2019", dilution="1+31")
    assert b.dilution_fraction == pytest.approx(1 / 32, rel=1e-15)
    np.testing.assert_allclose(b.totals, conc.totals / 32, rtol=1e-14, atol=0)


def test_rodinal_ph_falls_with_dilution(library, species):
    """DEV 5.2: the hydroxide/phenolate developer loses pH on dilution."""
    phs = [ws(library, species, "rodinal", SIT, dilution=d, t=20.0).ph for d in ("1+25", "1+50")]
    assert phs[1] < phs[0]


# ---------------------------------------------------------------------------------------
# Bookkeeping (spec D1 2.3.1, 2.3.2)
# ---------------------------------------------------------------------------------------


def test_working_solution_conserves_mass(library, species):
    """Element totals (S, B, Br, agents) equal the formula, relative difference <= 1e-13."""
    xtol = ws(library, species, "xtol")
    expected = {
        "sulfite": 85.0 / 126.04 + 2 * 3.5 / 190.11,  # S
        "borate": 4.0 / 137.86,  # B
        "dimezone-s": 0.2 / 206.24,
        "ascorbate": 12.0 / 216.12,
        "dtpa": 1.0 / 503.26,
        "bromide": 0.0,
    }
    for cid, ref in expected.items():
        assert xtol.total(cid) == pytest.approx(ref, rel=1e-13, abs=0.0)
    d76 = ws(library, species, "d-76")
    expected = {
        "sulfite": 100.0 / 126.04,
        "borate": 4 * 1.53 / 291.29 + 2 * 2.0 / 69.62,
        "hydroquinone": 5.0 / 110.11,
        "metol": 2 * 2.0 / 344.39,
        "dtpa": 0.5 / 393.35,
    }
    for cid, ref in expected.items():
        assert d76.total(cid) == pytest.approx(ref, rel=1e-13, abs=0.0)
    hc = ws(library, species, "hc-110-pre2019", dilution="1+31")
    f = 1 / 32
    expected = {
        "bromide": f * 23.0 / 119.00,
        "sulfite": f * 201.25 / 64.07,
        "hydroquinone": f * 86.25 / 110.11,
        "methylphenidone": f * 8.625 / 176.21,
        "diethanolamine": f * 373.75 / 105.14,
        "ethanolamine": f * 86.25 / 61.08,
    }
    for cid, ref in expected.items():
        assert hc.total(cid) == pytest.approx(ref, rel=1e-13, abs=0.0)


def test_toth_and_cations(library, species):
    sol = ws(library, species, "d-76-published")
    metol, borax = 2.0 / 344.39, 2.0 / 381.37
    assert sol.toth == pytest.approx(2 * metol - 2 * borax, rel=1e-13)  # spec 2.3.1 examples
    assert sol.cation_molality == pytest.approx(2 * 100.0 / 126.04 + 2 * borax, rel=1e-13)


def test_component_ids_and_speciation(library, species):
    sol = ws(library, species, "xtol")
    table = component_table(species)
    assert sol.component_ids == tuple(c for c in table if c not in BUILTIN_COMPONENTS)
    assert sol.totals.shape == (len(sol.component_ids),)
    assert not sol.totals.flags.writeable
    fam = [f.component for f in sol.speciation.families]
    assert fam == [c for c in sol.component_ids if table[c].family is not None]
    assert [sol.component_ids[i] for i in sol.family_index] == fam
    # The stored pH solves the proton condition at the stored speciation.
    ph = solve_ph(sol.speciation, sol.family_totals, np.asarray(sol.toth))
    assert float(ph) == pytest.approx(sol.ph, abs=1e-10)


def test_ionic_strength_is_a_fixed_point(library, species):
    sol = ws(library, species, "d-76-published")
    charges = component_charges(species)
    alpha = fractions(sol.speciation, np.asarray(sol.ph))
    h, oh = water_ions(sol.speciation, np.asarray(sol.ph))
    acc = sol.cation_molality + float(h) + float(oh)
    acc += 4 * 2.0 / 344.39  # SO4 2- of metol sulfate
    for f, fam in enumerate(sol.speciation.families):
        z2 = np.asarray(fam.charges, dtype=float) ** 2
        acc += sol.family_totals[f] * float(z2 @ alpha[f, : z2.size])
    for cid, total in zip(sol.component_ids, sol.totals, strict=True):
        if cid not in fam_ids(sol):
            acc += charges.get(cid, 0.0) ** 2 * total
    assert 0.5 * acc == pytest.approx(sol.ionic_strength, abs=2e-6)


def fam_ids(sol: WorkingSolution) -> set[str]:
    return {f.component for f in sol.speciation.families}


def test_component_charges(species):
    charges = component_charges(species)
    assert charges["bromide"] == -1.0
    assert charges["sulfite"] == -2.0
    assert charges["borate"] == 0.0
    assert charges["hydroquinone"] == 0.0


def test_temperature_changes_ph(library, species):
    """Borate and sulfite pK rise on cooling, so a borax developer is more alkaline at 20 C."""
    t20 = ws(library, species, "d-76-published", t=20.0)
    t25 = ws(library, species, "d-76-published", t=25.0)
    assert t20.temperature_c == 20.0
    assert t20.ph > t25.ph


def test_class_pka_from_kinetics(library, species):
    """The pyrazolidone class pK of kinetics.yaml reaches the XTOL Dimezone S family."""
    low = ws(library, species, "xtol", make_kinetics(pyrazolidone_pka=7.6))
    high = ws(library, species, "xtol", make_kinetics(pyrazolidone_pka=9.4))
    assert low.ph < high.ph


def test_delta_eps_shifts_ph(library, species):
    """A positive delta eps for the sulfite class raises the bisulfite pK' and the pH."""
    base = ws(library, species, "xtol", SIT)
    shifted = ws(library, species, "xtol", make_kinetics(sulfite=0.1))
    assert shifted.ph > base.ph


# ---------------------------------------------------------------------------------------
# Overrides (nuisance parameters, spec D1 2.3.2, 5.5)
# ---------------------------------------------------------------------------------------


def test_overrides_replace_amounts(library, species):
    base = ws(library, species, "d-76")
    zero = ws(library, species, "d-76", **{"boric-oxide": 0.0})
    more = ws(library, species, "d-76", **{"boric-oxide": 3.8})
    assert zero.total("borate") == pytest.approx(4 * 1.53 / 291.29, rel=1e-13)
    assert more.total("borate") == pytest.approx(4 * 1.53 / 291.29 + 2 * 3.8 / 69.62, rel=1e-13)
    assert zero.ph > base.ph > more.ph  # boric acid lowers the pH
    # Other components are untouched.
    assert zero.total("sulfite") == base.total("sulfite")


def test_override_equals_edited_file(library, species):
    bath = library.bath("d-76")
    comps = [c.model_dump(mode="json") for c in bath.components]
    for c in comps:
        if c["species"] == "boric-oxide":
            c["amount"]["value"] = 1.0
    edited = Bath.model_validate({**bath.model_dump(mode="json"), "components": comps})
    a = working_solution(edited, species, SIT, dilution="1+1", temperature_c=20.0)
    b = ws(library, species, "d-76", dilution="1+1", t=20.0, **{"boric-oxide": 1.0})
    np.testing.assert_array_equal(a.totals, b.totals)
    assert a.ph == b.ph


def test_override_unknown_ingredient_raises(library, species):
    with pytest.raises(ValueError, match="not an ingredient"):
        ws(library, species, "d-76", **{"sodium-bromide": 1.0})


@pytest.mark.parametrize("value", [-1.0, float("nan"), float("inf")])
def test_override_bad_amount_raises(library, species, value):
    with pytest.raises(ValueError, match="overrides"):
        ws(library, species, "d-76", **{"boric-oxide": value})


# ---------------------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize("dilution", ["1:1", "half", "0+1"])
def test_bad_dilution_raises(library, species, dilution):
    with pytest.raises(ValueError, match="dilution"):
        ws(library, species, "d-76", dilution=dilution)


@pytest.mark.parametrize("t", [float("nan"), -300.0])
def test_bad_temperature_raises(library, species, t):
    with pytest.raises(ValueError, match="temperature"):
        ws(library, species, "d-76", t=t)


def test_unknown_species_raises(library, species):
    bath = _with_component(library.bath("d-76"), "unobtainium", 1.0)
    with pytest.raises(ValueError, match="unknown species"):
        working_solution(bath, species, SIT, dilution="stock", temperature_c=20.0)


def test_charged_component_weighed_directly_raises(library, species):
    bath = _with_component(library.bath("d-76-published"), "sulfite", 1.0)
    with pytest.raises(ValueError, match="weighed directly"):
        working_solution(bath, species, SIT, dilution="stock", temperature_c=20.0)


def test_non_neutral_ingredient_raises(species):
    bad = dict(species)
    bad["sodium-sulfite"] = species["sodium-sulfite"].model_copy(update={"toth": 1.0})
    with pytest.raises(ValueError, match="electroneutral"):
        component_charges(bad)


def test_total_unknown_component_raises(library, species):
    with pytest.raises(KeyError, match="not a component"):
        ws(library, species, "d-76").total("nothing")
