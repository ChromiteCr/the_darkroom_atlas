import numpy as np
import pytest
from pydantic import TypeAdapter, ValidationError

from darkroom_atlas.data.models import (
    BaseLayer,
    Bath,
    BathKind,
    EmulsionLayer,
    InlineCurve,
    Layer,
    Material,
    MaterialKind,
    Process,
    ProcessStep,
    ReferenceCurve,
    SpeciesFile,
    SpeciesRole,
    SpectralCurve,
    StepKind,
)
from darkroom_atlas.data.provenance import CurveOrigin

ASSUMED = {"source": "assumed"}


def _material(**overrides):
    data = {
        "schema_version": 1,
        "id": "m",
        "maker": "test",
        "name": "test-m",
        "kind": "ideal",
        "model": "ideal",
        "layers": [{"type": "emulsion", "id": "e"}],
    }
    data.update(overrides)
    return Material.model_validate(data)


def _bath(**overrides):
    data = {
        "schema_version": 1,
        "id": "b",
        "name": "test-b",
        "kind": "developer",
        "formula_confidence": "published",
        "components": [],
    }
    data.update(overrides)
    return Bath.model_validate(data)


def test_enum_values():
    assert {k.value for k in StepKind} == {
        "identity", "presoak", "develop", "stop", "fix", "bleach", "bleach_fix",
        "reversal", "wash", "stabilize", "tone", "dry",
    }  # fmt: skip
    assert StepKind.identity is StepKind.IDENTITY
    assert MaterialKind.ideal == "ideal"
    assert BathKind.bleach_fix == "bleach_fix"
    assert SpeciesRole.color_developing_agent == "color_developing_agent"


def test_layer_union_dispatch():
    adapter = TypeAdapter(Layer)
    layer = adapter.validate_python({"type": "base", "id": "base", "material": "test-polymer"})
    assert isinstance(layer, BaseLayer)
    assert layer.thickness is None
    emulsion = adapter.validate_python({"type": "emulsion", "id": "e"})
    assert isinstance(emulsion, EmulsionLayer)
    assert emulsion.couplers == [] and emulsion.parameters == {}
    with pytest.raises(ValidationError):
        adapter.validate_python({"type": "nope", "id": "a"})
    with pytest.raises(ValidationError):
        adapter.validate_python({"id": "a"})


def test_material_layer_ids_unique():
    with pytest.raises(ValidationError, match="duplicate layer id"):
        _material(layers=[{"type": "emulsion", "id": "e"}, {"type": "emulsion", "id": "e"}])


@pytest.mark.parametrize("bad_id", ["frame", "a/b", "", "has space"])
def test_layer_id_rules(bad_id):
    with pytest.raises(ValidationError):
        _material(layers=[{"type": "emulsion", "id": bad_id}])


def test_material_requires_emulsion():
    with pytest.raises(ValidationError, match="emulsion"):
        _material(layers=[{"type": "overcoat", "id": "o"}])


def test_real_material_requires_base():
    with pytest.raises(ValidationError, match="base"):
        _material(kind="bw_negative", layers=[{"type": "emulsion", "id": "e"}])
    with pytest.raises(ValidationError, match="base"):
        _material(
            kind="bw_negative",
            layers=[
                {"type": "emulsion", "id": "e"},
                {"type": "base", "id": "b1", "material": "p"},
                {"type": "base", "id": "b2", "material": "p"},
            ],
        )
    ideal = _material(kind="ideal", layers=[{"type": "emulsion", "id": "e"}])
    assert ideal.kind is MaterialKind.IDEAL
    real = _material(
        kind="bw_negative",
        layers=[{"type": "emulsion", "id": "e"}, {"type": "base", "id": "b", "material": "p"}],
    )
    assert [layer.id for layer in real.emulsion_layers] == ["e"]


def test_material_quantity_dimension_checked():
    with pytest.raises(ValidationError, match="thickness"):
        _material(
            layers=[
                {"type": "emulsion", "id": "e", "thickness": {"value": 1, "unit": "s", **ASSUMED}}
            ]
        )


def test_spectral_curve_exactly_one():
    inline = {"wavelengths_nm": [400, 500], "values": [1, 2]}
    with pytest.raises(ValidationError, match="exactly one"):
        SpectralCurve(quantity="q", source="assumed")
    with pytest.raises(ValidationError, match="exactly one"):
        SpectralCurve(quantity="q", source="assumed", inline=inline, builtin="cie1931_y")
    with pytest.raises(ValidationError, match="exactly one"):
        SpectralCurve(quantity="q", source="assumed", file="c.csv", builtin="cie1931_y")
    assert SpectralCurve(quantity="q", source="assumed", builtin="cie1931_y").builtin == "cie1931_y"
    assert SpectralCurve(quantity="q", source="assumed", inline=inline).inline is not None
    assert SpectralCurve(quantity="q", source="assumed", file="curves/c.csv").file == "curves/c.csv"


def test_spectral_curve_rules():
    with pytest.raises(ValidationError, match="builtin"):
        SpectralCurve(quantity="q", source="assumed", builtin="nope")
    with pytest.raises(ValidationError, match="ref"):
        SpectralCurve(quantity="q", source="datasheet", builtin="cie1931_y")
    for bad in ["/abs/c.csv", "../c.csv", "a/../../c.csv", ""]:
        with pytest.raises(ValidationError):
            SpectralCurve(quantity="q", source="assumed", file=bad)


def test_inline_curve_monotonic():
    with pytest.raises(ValidationError, match="increasing"):
        InlineCurve(wavelengths_nm=[400, 400, 500], values=[1, 2, 3])
    with pytest.raises(ValidationError, match="increasing"):
        InlineCurve(wavelengths_nm=[500, 400], values=[1, 2])
    with pytest.raises(ValidationError, match="same length"):
        InlineCurve(wavelengths_nm=[400, 500, 600], values=[1, 2])
    with pytest.raises(ValidationError):
        InlineCurve(wavelengths_nm=[400], values=[1])
    with pytest.raises(ValidationError):
        InlineCurve(wavelengths_nm=[400, 500], values=[1, float("nan")])
    ok = InlineCurve(wavelengths_nm=[400, 500], values=[1, 2])
    assert ok.wavelengths_nm == [400.0, 500.0]


@pytest.mark.parametrize("kind", ["develop", "fix", "bleach", "bleach_fix", "tone"])
def test_step_requires_bath(kind):
    with pytest.raises(ValidationError, match="bath"):
        ProcessStep(kind=kind, time="5:00")
    assert ProcessStep(kind=kind, bath="b", time="5:00").bath == "b"


def test_step_requires_time():
    with pytest.raises(ValidationError, match="time"):
        ProcessStep(kind="wash")
    with pytest.raises(ValidationError, match="time"):
        ProcessStep(kind="wash", time=0)
    with pytest.raises(ValidationError, match="time"):
        ProcessStep(kind="wash", time="0:00")
    with pytest.raises(ValidationError):
        ProcessStep(kind="wash", time="9:75")
    identity = ProcessStep(kind="identity")
    assert identity.time_s is None
    assert identity.dilution_fraction == 1.0


def test_step_time_and_dilution():
    step = ProcessStep(kind="develop", bath="b", time="9:45", dilution="1+1")
    assert step.time_s == 585.0
    assert step.dilution_fraction == 0.5
    assert ProcessStep(kind="wash", time=90).time_s == 90.0
    assert ProcessStep(kind="wash", time="1.5min").time_s == 90.0
    with pytest.raises(ValidationError, match="dilution"):
        ProcessStep(kind="wash", time=60, dilution="1:1")


def test_step_agitation_and_tank():
    step = ProcessStep(
        kind="develop",
        bath="b",
        time="5:00",
        temperature={"value": 20, "unit": "C", **ASSUMED},
        agitation={"method": "inversion", "initial_s": 30, "interval_s": 60, "duration_s": 10},
        tank={"volume_ml": 500, "film_area_cm2": 80},
    )
    assert step.agitation.method == "inversion"
    assert step.tank.volume_ml == 500.0
    with pytest.raises(ValidationError):
        ProcessStep(kind="wash", time=60, agitation={"method": "shake"})
    with pytest.raises(ValidationError):
        ProcessStep(kind="wash", time=60, tank={"volume_ml": 0})
    with pytest.raises(ValidationError, match="temperature"):
        ProcessStep(kind="wash", time=60, temperature={"value": 20, "unit": "um", **ASSUMED})


def test_process_requires_steps():
    base = {"schema_version": 1, "id": "p", "name": "test-p", "family": "bw"}
    with pytest.raises(ValidationError):
        Process.model_validate({**base, "steps": []})
    with pytest.raises(ValidationError):
        Process.model_validate({**base, "schema_version": 2, "steps": [{"kind": "identity"}]})
    proc = Process.model_validate({**base, "steps": [{"kind": "identity"}]})
    assert proc.steps[0].kind is StepKind.identity


def test_bad_dilution_in_bath():
    with pytest.raises(ValidationError, match="1:1"):
        _bath(standard_dilutions=["1:1"])
    assert _bath().standard_dilutions == ["stock"]
    assert _bath(standard_dilutions=["stock", "1+1"]).standard_dilutions == ["stock", "1+1"]


def test_bath_component_species_unique():
    amount = {"value": 1, "unit": "g/L", **ASSUMED}
    with pytest.raises(ValidationError, match="duplicate"):
        _bath(components=[{"species": "a", "amount": amount}, {"species": "a", "amount": amount}])


def test_species_ids_unique():
    species = {
        "id": "test-a",
        "name": "test-a",
        "molar_mass": {"value": 100, "unit": "g/mol", **ASSUMED},
        "roles": ["developing_agent"],
    }
    with pytest.raises(ValidationError, match="duplicate species id"):
        SpeciesFile.model_validate({"schema_version": 1, "species": [species, species]})
    ok = SpeciesFile.model_validate(
        {"schema_version": 1, "species": [species, {**species, "id": "test-b"}]}
    )
    assert [s.id for s in ok.species] == ["test-a", "test-b"]
    assert ok.species[0].roles == [SpeciesRole.developing_agent]


def test_species_rules():
    species = {
        "id": "test-a",
        "name": "test-a",
        "molar_mass": {"value": 100, "unit": "g/mol", **ASSUMED},
        "roles": ["developing_agent"],
    }
    with pytest.raises(ValidationError, match="molar_mass"):
        SpeciesFile.model_validate(
            {
                "schema_version": 1,
                "species": [{**species, "molar_mass": {"value": 1, "unit": "um", **ASSUMED}}],
            }
        )
    with pytest.raises(ValidationError, match="CAS"):
        SpeciesFile.model_validate({"schema_version": 1, "species": [{**species, "cas": "1-2-3"}]})
    # Synthetic numbers: 1234-56-6 has a correct check digit, 1234-56-7 does not.
    with pytest.raises(ValidationError, match="check digit"):
        SpeciesFile.model_validate(
            {"schema_version": 1, "species": [{**species, "cas": "1234-56-7"}]}
        )
    ok = SpeciesFile.model_validate(
        {"schema_version": 1, "species": [{**species, "cas": "1234-56-6"}]}
    )
    assert ok.species[0].cas == "1234-56-6"


def test_reference_curve_model():
    curve = ReferenceCurve(
        material="m",
        quantity="q",
        x_unit="1",
        y_unit="1",
        origin="datasheet",
        ref="test ref",
        conditions=None,
        x=[0.0, 1.0, 2.0],
        y=[1.0, 2.0, 3.0],
    )
    assert curve.origin is CurveOrigin.DATASHEET
    assert isinstance(curve.x, np.ndarray) and curve.x.dtype == np.float64
    assert not curve.x.flags.writeable and not curve.y.flags.writeable
    with pytest.raises(ValidationError):
        ReferenceCurve(
            material="m", quantity="q", x_unit="1", y_unit="1", origin="datasheet",
            ref="r", conditions=None, x=[0.0, 1.0], y=[1.0],
        )  # fmt: skip
    with pytest.raises(ValidationError):
        ReferenceCurve(
            material="m", quantity="q", x_unit="1", y_unit="1", origin="datasheet",
            ref="", conditions=None, x=[0.0, 1.0], y=[1.0, 2.0],
        )  # fmt: skip
