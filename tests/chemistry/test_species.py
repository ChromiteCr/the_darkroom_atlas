"""Species data file and acid-base families (spec D1 2.3.1, 2.3.3, 2.3.7)."""

import math

import numpy as np
import pytest

from darkroom_atlas.chemistry.species import (
    GAS_CONSTANT,
    AcidBaseFamily,
    ComponentInfo,
    component_table,
    manov_boric_pka,
)
from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import iter_quantities
from darkroom_atlas.data.provenance import Source

# Spec D1 2.3.7: D_w (um2/s, 20 C) and r_gel per component.
TABLE_237 = {
    "sulfite": (1030, 0.5),
    "borate": (1030, 0.5),
    "hydroquinone": (855, 0.10),
    "metol": (740, 0.13),
    "p-aminophenol": (825, 0.12),
    "dimezone-s": (555, 0.2),
    "phenidone": (650, 0.2),
    "methylphenidone": (605, 0.2),
    "ascorbate": (710, 0.2),
    "diethanolamine": (785, 0.4),
    "ethanolamine": (1090, 0.4),
    "methylethanolamine": (1090, 0.4),
    "dtpa": (555, 0.2),
    "hedta": (555, 0.2),
    "bromide": (1816, 0.5),
    "iodide": (1786, 0.47),
    "silver": (1030, 0.5),
}

# Spec D1 2.3.7: pK values at 25 C, most protonated first (placeholders for pyrazolidones).
PKA_237 = {
    "sulfite": (1.86, 7.19),
    "borate": (9.235,),
    "hydroquinone": (9.85, 11.45),
    "metol": (5.9, 10.25),
    "p-aminophenol": (5.48, 10.46),
    "ascorbate": (4.04, 11.34),
    "diethanolamine": (8.883,),
    "ethanolamine": (9.498,),
    "methylethanolamine": (9.82,),
    "dtpa": (1.8, 2.6, 4.3, 8.6, 10.5),
    "hedta": (2.6, 5.4, 9.8),
    "gelatin": (10.5,),
}

# Spec D1 2.3.5: kinetic class, role and charge of the active form.
AGENTS_235 = {
    "hydroquinone": ("hydroquinone", "regenerator", -2),
    "ascorbate": ("ascorbate", "regenerator", -1),
    "metol": ("metol", "eta", -1),
    "dimezone-s": ("pyrazolidone", "eta", -1),
    "phenidone": ("pyrazolidone", "eta", -1),
    "methylphenidone": ("pyrazolidone", "eta", -1),
    "p-aminophenol": ("p-aminophenol", "single", -1),
}

# Charges of the components without an acid-base family (for the electroneutrality check).
FIXED_CHARGE = {"bromide": -1, "iodide": -1}


@pytest.fixture(scope="module")
def species():
    return Library().species()


@pytest.fixture(scope="module")
def table(species):
    return component_table(species)


def test_species_file_has_every_component_of_237(table):
    assert set(TABLE_237) | {"gelatin"} <= set(table)
    for cid, (d_w, r_gel) in TABLE_237.items():
        info = table[cid]
        assert isinstance(info, ComponentInfo)
        assert info.d_water_um2_s == pytest.approx(d_w, rel=1e-12), cid
        assert info.gel_ratio == pytest.approx(r_gel, rel=1e-12), cid
        assert not info.immobile


def test_gelatin_is_immobile(table):
    gel = table["gelatin"]
    assert gel.immobile
    assert gel.d_water_um2_s == 0.0 and gel.gel_ratio == 0.0
    assert gel.family.charges == (1, 0)
    assert gel.family.protons_rel == (1, 0)


def test_component_table_skips_ingredients(species, table):
    assert "sodium-sulfite" in species
    assert "sodium-sulfite" not in table
    assert all(species[k].kind == "component" for k in table)


def test_pka_values_of_237(table):
    for cid, pka in PKA_237.items():
        assert table[cid].family.pka25 == pytest.approx(pka, abs=1e-12), cid


def test_pyrazolidone_pka_is_an_assumed_placeholder_in_the_prior_box(species):
    for cid in ("dimezone-s", "phenidone", "methylphenidone"):
        (q,) = species[cid].acid_base.pka
        assert q.source is Source.ASSUMED
        assert q.bounds == (7.5, 9.5)
        assert "PLACEHOLDER" in q.note


def test_family_states(table):
    sulfite = table["sulfite"].family
    assert sulfite.charges == (0, -1, -2)
    assert sulfite.protons_rel == (2, 1, 0)
    borate = table["borate"].family
    assert borate.charges == (0, -1)
    assert borate.protons_rel == (0, -1)
    metol = table["metol"].family
    assert metol.charges == (1, 0, -1)
    assert metol.protons_rel == (1, 0, -1)
    dtpa = table["dtpa"].family
    assert dtpa.charges == (0, -1, -2, -3, -4, -5)
    assert dtpa.protons_rel == (0, -1, -2, -3, -4, -5)
    assert table["bromide"].family is None


def test_developing_agents_of_235(table):
    for cid, (kinetic_class, role, charge) in AGENTS_235.items():
        agent = table[cid].agent
        assert agent.kinetic_class == kinetic_class, cid
        assert agent.role == role, cid
        assert table[cid].family.charges[agent.active_state] == charge, cid
    assert table["hydroquinone"].agent.product_sulfonated
    assert not table["ascorbate"].agent.product_sulfonated


def _component_reference_charge(table, cid):
    family = table[cid].family
    if family is None:
        return FIXED_CHARGE[cid]
    return family.charges[family.protons_rel.index(0)]


def test_ingredients_are_electroneutral(species, table):
    """sum_k s_k z_ref,k + toth + cations + spectators = 0 (spec D1 2.3.1)."""
    ingredients = [s for s in species.values() if s.kind == "ingredient"]
    assert len(ingredients) >= 20
    for s in ingredients:
        charge = sum(n * _component_reference_charge(table, k) for k, n in s.dissolves_to.items())
        charge += s.toth + sum(s.cations.values()) + sum(s.spectator_charge.values())
        assert charge == pytest.approx(0.0, abs=1e-12), s.id


def test_toth_examples_of_231(species):
    assert species["sodium-metabisulfite"].toth == 2
    assert species["sodium-metabisulfite"].dissolves_to == {"sulfite": 2}
    assert species["borax-decahydrate"].toth == -2
    assert species["borax-decahydrate"].dissolves_to == {"borate": 4}
    assert species["metol-sulfate"].toth == 2
    assert species["metol-sulfate"].dissolves_to == {"metol": 2}


def test_species_file_provenance(species):
    for s in species.values():
        for loc, q in iter_quantities(s):
            assert q.source is not Source.FITTED, f"{s.id}.{loc}"
            if q.source is Source.ASSUMED:
                assert q.note, f"{s.id}.{loc}: an assumed value needs a note"


def test_species_file_validates_except_for_missing_kinetics():
    """The only species-file problems allowed are references to kinetics.yaml."""
    problems = [p for p in Library().validate() if "chemistry/species.yaml" in p]
    other = [p for p in problems if "kinetic_class" not in p]
    assert other == []


# ---------------------------------------------------------------------------------------
# Temperature laws (spec D1 2.3.3)
# ---------------------------------------------------------------------------------------


def test_boric_acid_pka_temperature(table):
    borate = table["borate"].family
    assert borate.pka(293.15)[0] == pytest.approx(9.278, abs=1e-3)
    assert borate.pka(298.15)[0] == pytest.approx(9.235, abs=1e-3)
    # The Manov law ignores the stored 25 C value; the two must agree.
    assert borate.pka25[0] == pytest.approx(manov_boric_pka(298.15), abs=1e-3)


def test_van_t_hoff_sulfite(table):
    sulfite = table["sulfite"].family
    t = 293.15
    expected = 7.19 + (-3.65e3) / (GAS_CONSTANT * math.log(10)) * (1 / t - 1 / 298.15)
    pk = sulfite.pka(t)
    assert pk[1] == pytest.approx(expected, abs=1e-12)
    assert pk[0] == pytest.approx(1.86, abs=1e-12)  # zero enthalpy
    assert pk[1] < 7.19  # exothermic step: pK falls on cooling
    np.testing.assert_allclose(sulfite.pka(298.15), (1.86, 7.19), atol=1e-12)


def test_linear_law_amines(table):
    dea = table["diethanolamine"].family
    assert dea.pka(293.15)[0] == pytest.approx(8.883 + 0.025 * 5, abs=1e-12)
    # Close to the measured 9.005 at 20 C (Bower, Robinson & Bates 1962; DEV 6.2).
    assert dea.pka(293.15)[0] == pytest.approx(9.005, abs=0.005)
    mea = table["ethanolamine"].family
    assert mea.pka(293.15)[0] == pytest.approx(9.647, abs=0.005)


def test_organic_agents_have_no_temperature_dependence(table):
    hq = table["hydroquinone"].family
    np.testing.assert_array_equal(hq.pka(288.15), hq.pka(298.15))


def test_pka_is_read_only(table):
    pk = table["sulfite"].family.pka(293.15)
    with pytest.raises(ValueError):
        pk[0] = 1.0


@pytest.mark.parametrize("temperature", [0.0, -5.0, float("nan"), float("inf")])
def test_pka_rejects_bad_temperature(table, temperature):
    with pytest.raises(ValueError, match="temperature"):
        table["sulfite"].family.pka(temperature)


def test_family_validation():
    with pytest.raises(ValueError, match="charges and protons_rel"):
        AcidBaseFamily("x", (0,), (0,), (7.0,), "van_t_hoff", (0.0,), (), "organic")
    with pytest.raises(ValueError, match="enthalpy"):
        AcidBaseFamily("x", (0, -1), (0, -1), (7.0,), "van_t_hoff", (), (), "organic")
    with pytest.raises(ValueError, match="dpka_dt"):
        AcidBaseFamily("x", (1, 0), (1, 0), (9.0,), "linear", (), (), "amine")


def test_from_species_without_acid_base_raises(species):
    with pytest.raises(ValueError, match="no acid_base"):
        AcidBaseFamily.from_species(species["bromide"])


def test_component_table_key_mismatch(species):
    with pytest.raises(ValueError, match="differs"):
        component_table({"not-sulfite": species["sulfite"]})
