"""The seven D1 developer files (spec D1 1.1, 1.3 decision 1, 5.5, 5.8, 7.1, 7.2)."""

import math

import pytest

from darkroom_atlas.chemistry.solution import working_solution
from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import (
    ActivityParams,
    GelatinParams,
    KineticsFile,
    iter_quantities,
)
from darkroom_atlas.data.provenance import FormulaConfidence, Quantity, Source
from darkroom_atlas.units import parse_dilution

# id -> (formula_confidence, variant, base_solution), spec D1 1.1 and 5.5.
DEVELOPERS = {
    "d-76": ("msds-estimated", "packaged", "stock"),
    "d-76-published": ("published", "published", "stock"),
    "xtol": ("patent", "patent-example-1", "stock"),
    "hc-110-pre2019": ("msds-estimated", "pre2019", "concentrate"),
    "hc-110-2019": ("msds-estimated", "2019", "concentrate"),
    "t-max-dev": ("msds-estimated", "generation-I", "concentrate"),
    "rodinal": ("msds-estimated", "adox-2024", "concentrate"),
}

# Spec D1 7.2 anchor table, one per file: (dilution, temperature C, pH, uncertainty, method).
ANCHORS = {
    "d-76": [("stock", 25.0, 8.5, 0.05, "typical")],
    "d-76-published": [("stock", 25.0, 8.5, 0.1, "colorimetric")],
    "xtol": [("stock", 25.0, 8.20, 0.05, "glass_electrode")],
    "hc-110-pre2019": [("1+15", 25.0, 8.9, 0.10, "typical")],
    "hc-110-2019": [],
    "t-max-dev": [("1+4", 25.0, 8.37, 0.05, "typical")],
    "rodinal": [],
}

# Nuisance ingredients fitted within bounds (spec D1 7.1 stages A and D).
NUISANCE = {
    "d-76": {"boric-oxide"},
    "t-max-dev": {"sulfur-dioxide", "hydroquinone", "dimezone-s", "sodium-bromide"},
    "hc-110-pre2019": {
        "potassium-hydroxide",
        "hydroquinone",
        "methylphenidone",
        "potassium-bromide",
    },
}

# Dilutions as explicit factors (spec D1 1.3 decision 1, 5.5).
DILUTIONS = {
    "d-76": ["stock", "1+1"],
    "d-76-published": ["stock", "1+1"],
    "xtol": ["stock", "1+1", "1+2", "1+3"],
    "hc-110-pre2019": ["1+15", "1+31", "1+19", "1+39", "1+47", "1+79"],
    "hc-110-2019": ["1+15", "1+31", "1+19", "1+39", "1+47", "1+79"],
    "t-max-dev": ["1+4"],
    "rodinal": ["1+25", "1+50"],
}

# Initial bromide of the base solution in g/L (KBr or NaBr), from the reports.
BROMIDE = {
    "d-76": None,
    "d-76-published": None,
    "xtol": None,
    "hc-110-pre2019": ("potassium-bromide", 23.0),
    "hc-110-2019": ("sodium-bromide", 21.4),
    "t-max-dev": ("sodium-bromide", 0.0),
    "rodinal": ("potassium-bromide", 24.0),
}


def _q(value, unit="1"):
    return {"value": value, "unit": unit, "source": "assumed", "note": "test"}


SIT = KineticsFile.model_construct(
    activity=ActivityParams.model_validate(
        {
            "debye_huckel_A_20C": _q(0.5050),
            "debye_huckel_A_25C": _q(0.5091),
            "sit_b": _q(1.5),
            "delta_eps": {},
        }
    ),
    agents={},
    gelatin=GelatinParams.model_construct(pka=Quantity.model_validate(_q(10.5))),
)


@pytest.fixture(scope="module")
def library():
    return Library()


@pytest.fixture(scope="module")
def species(library):
    return library.species()


@pytest.mark.parametrize("bath_id", sorted(DEVELOPERS))
def test_developer_file_loads(library, bath_id):
    bath = library.bath(bath_id)
    confidence, variant, base = DEVELOPERS[bath_id]
    assert bath.id == bath_id
    assert bath.kind.value == "developer"
    assert bath.formula_confidence is FormulaConfidence(confidence)
    assert bath.variant == variant
    assert bath.base_solution == base
    assert bath.refs


@pytest.mark.parametrize("bath_id", sorted(DEVELOPERS))
def test_developer_components_are_known_and_in_g_per_l(library, species, bath_id):
    bath = library.bath(bath_id)
    for c in bath.components:
        assert c.species in species, c.species
        assert c.amount.unit == "g/L"
        assert c.amount.value >= 0.0


@pytest.mark.parametrize("bath_id", sorted(DEVELOPERS))
def test_developer_provenance(library, bath_id):
    """assumed values carry a note, literature/datasheet values a ref; nothing is fitted
    before the D1e calibration; msds-estimated amounts all have bounds."""
    bath = library.bath(bath_id)
    for loc, q in iter_quantities(bath):
        assert q.source is not Source.FITTED, loc
        if q.source is Source.ASSUMED:
            assert q.note, f"{bath_id}: {loc}: an assumed value needs a note"
        else:
            assert q.ref, f"{bath_id}: {loc}: a {q.source.value} value needs a ref"
    if bath.formula_confidence is FormulaConfidence.MSDS_ESTIMATED:
        for c in bath.components:
            assert c.amount.bounds is not None, f"{bath_id}: {c.species} needs SDS bounds"


@pytest.mark.parametrize("bath_id", sorted(DEVELOPERS))
def test_developer_nuisance(library, bath_id):
    bath = library.bath(bath_id)
    assert set(bath.nuisance) == NUISANCE.get(bath_id, set())
    amounts = {c.species: c.amount for c in bath.components}
    for name in bath.nuisance:
        assert amounts[name].bounds is not None


def test_packaged_d76_boric_oxide_nuisance(library):
    """Spec D1 5.8: B2O3 in packaged D-76 is a bounded nuisance, 0-3.8 g/L."""
    amounts = {c.species: c.amount for c in library.bath("d-76").components}
    assert amounts["boric-oxide"].bounds == (0.0, 3.8)
    assert "boric-oxide" not in {c.species for c in library.bath("d-76-published").components}


@pytest.mark.parametrize("bath_id", sorted(DEVELOPERS))
def test_developer_ph_anchors(library, bath_id):
    """Exactly the spec D1 7.2 anchors, with their published values and sources."""
    bath = library.bath(bath_id)
    got = [
        (
            a.dilution,
            float(a.temperature.to("C")),
            float(a.ph.value),
            a.ph.uncertainty,
            a.method,
        )
        for a in bath.ph_anchors
    ]
    assert got == ANCHORS[bath_id]
    for a in bath.ph_anchors:
        assert a.ph.source in (Source.DATASHEET, Source.LITERATURE)
        assert a.ph.ref


@pytest.mark.parametrize("bath_id", sorted(DEVELOPERS))
def test_developer_dilutions_are_explicit(library, bath_id):
    bath = library.bath(bath_id)
    assert bath.standard_dilutions == DILUTIONS[bath_id]
    for d in bath.standard_dilutions:
        assert d == "stock" or "+" in d
        parse_dilution(d)


def test_hc110_dilution_b_is_1_plus_31(library):
    """DEV 3.3: B is 1+31 from the concentrate (factor 32)."""
    for bath_id in ("hc-110-pre2019", "hc-110-2019"):
        assert "1+31" in library.bath(bath_id).standard_dilutions
    assert parse_dilution("1+31") == pytest.approx(1 / 32)


@pytest.mark.parametrize("bath_id", sorted(DEVELOPERS))
def test_developer_initial_bromide(library, bath_id):
    amounts = {c.species: c.amount for c in library.bath(bath_id).components}
    halides = {k for k in amounts if k in ("potassium-bromide", "sodium-bromide")}
    if BROMIDE[bath_id] is None:
        assert halides == set()
    else:
        name, value = BROMIDE[bath_id]
        assert halides == {name}
        assert amounts[name].value == pytest.approx(value)


@pytest.mark.parametrize("bath_id", sorted(DEVELOPERS))
def test_every_standard_dilution_has_a_working_solution(library, species, bath_id):
    bath = library.bath(bath_id)
    for d in bath.standard_dilutions:
        sol = working_solution(bath, species, SIT, dilution=d, temperature_c=20.0)
        assert math.isfinite(sol.ph) and 7.0 < sol.ph < 13.0, (d, sol.ph)
        assert sol.ionic_strength > 0.0


def test_hc110_pre2019_alkali_meets_the_anchor(library, species):
    """The KOH centre was chosen to give the CIS-61 dilution-A pH with D(I), delta eps = 0."""
    sol = working_solution(
        library.bath("hc-110-pre2019"), species, SIT, dilution="1+15", temperature_c=25.0
    )
    assert sol.ph == pytest.approx(8.9, abs=0.01)


def test_developer_files_pass_check_data(library):
    """Library.validate reports nothing about the developer files."""
    problems = [p for p in library.validate() if "/developers/" in p]
    assert problems == []
