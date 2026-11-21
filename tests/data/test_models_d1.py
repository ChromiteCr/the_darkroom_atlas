"""D1 schema additions (spec D1 section 5): every new field optional, new file models."""

import copy
from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError

from darkroom_atlas.data.models import (
    Agitation,
    Bath,
    CalibrationRecord,
    EmulsionLayer,
    Grain,
    KineticsFile,
    Material,
    PhAnchor,
    Species,
    SpeciesFile,
    Tank,
    bw_mechanistic_problems,
    iter_quantities,
)

SAMPLES = Path(__file__).resolve().parent / "samples"
FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "library"


def _yaml(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _q(value, unit="1", **extra):
    return {"value": value, "unit": unit, "source": "assumed", **extra}


@pytest.fixture
def film() -> dict:
    return _yaml(FIXTURES / "films" / "test" / "bw-mechanistic.yaml")


# --------------------------------------------------------------------------------------
# D0 files are unchanged
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["bw-neg", "color-neg"])
def test_d0_materials_load_with_d1_defaults(name):
    material = Material.model_validate(_yaml(FIXTURES / "films" / "test" / f"{name}.yaml"))
    assert material.default_format is None
    assert material.formats == {}
    for layer in material.emulsion_layers:
        assert layer.optics is None and layer.latent is None and layer.development is None
        assert layer.gelatin_coverage is None and layer.swelling_ratio is None
        if layer.grain is not None:
            assert layer.grain.n_classes == 6
            assert layer.grain.thickness is None and layer.grain.iodide_fraction is None


def test_d0_species_default_to_ingredients():
    species = SpeciesFile.model_validate(_yaml(FIXTURES / "chemistry" / "species.yaml"))
    for s in species.species:
        assert s.kind == "ingredient"
        assert s.dissolves_to == {} and s.toth == 0.0 and s.acid_base is None


# --------------------------------------------------------------------------------------
# Films (spec D1 5.4)
# --------------------------------------------------------------------------------------


def test_fixture_film_has_every_d1_field(film):
    material = Material.model_validate(film)
    assert material.model == "bw-mechanistic"
    assert material.default_format == "135"
    assert sorted(material.formats) == ["120", "135"]
    assert bw_mechanistic_problems(material) == []
    (emulsion,) = material.emulsion_layers
    for name in ("gelatin_coverage", "swelling_ratio", "optics", "latent", "development"):
        assert getattr(emulsion, name) is not None
    grain = emulsion.grain
    assert grain.morphology == "tabular"
    assert grain.thickness.to("um") == 0.12
    assert grain.iodide_fraction.to("1") == 0.03
    assert emulsion.optics.intrinsic_absorption.edge_energy.to("eV") == 2.7
    assert len(emulsion.optics.dye_bands) == 2
    assert emulsion.latent.cluster_cap == 5
    assert emulsion.development.covering_power.unit == "1/(g/m2)"
    # every set field of 5.4 appears in the fixture (none left at a default)
    for model in (grain, emulsion.optics, emulsion.latent, emulsion.development):
        assert model.model_fields_set == set(type(model).model_fields)
    # every parameter is tagged assumed
    sources = {q.source.value for _, q in iter_quantities(material)}
    assert sources == {"assumed"}


def test_tabular_grain_needs_thickness():
    grain = {
        "morphology": "tabular",
        "mean_diameter": _q(1.0, "um"),
        "sigma_log": _q(0.3),
        "silver_coverage": _q(4.0, "g/m2"),
    }
    with pytest.raises(ValidationError, match="thickness: required for tabular"):
        Grain.model_validate(grain)
    Grain.model_validate({**grain, "morphology": "conventional"})


@pytest.mark.parametrize(
    ("change", "match"),
    [
        ({"n_classes": 1}, "greater than or equal to 2"),
        ({"n_classes": 65}, "less than or equal to 64"),
        ({"iodide_fraction": _q(1.5)}, r"iodide_fraction must be a number in \[0, 1\]"),
        ({"thickness": _q(0.1, "s")}, "not a length unit"),
    ],
)
def test_grain_rejects(change, match):
    grain = {
        "morphology": "conventional",
        "mean_diameter": _q(1.0, "um"),
        "sigma_log": _q(0.3),
        "silver_coverage": _q(4.0, "g/m2"),
        **change,
    }
    with pytest.raises(ValidationError, match=match):
        Grain.model_validate(grain)


def test_default_format_rules(film):
    del film["default_format"]
    with pytest.raises(ValidationError, match="default_format is required"):
        Material.model_validate(film)
    film["default_format"] = "220"
    with pytest.raises(ValidationError, match="not one of the base formats"):
        Material.model_validate(film)
    film["default_format"] = 135
    with pytest.raises(ValidationError, match="quote it"):
        Material.model_validate(film)


def test_format_ids_must_be_quoted(film):
    formats = film["layers"][2]["formats"]
    formats[135] = formats.pop("135")
    with pytest.raises(ValidationError, match="quote it in YAML"):
        Material.model_validate(film)


def test_default_format_without_formats_rejected():
    data = _yaml(FIXTURES / "films" / "test" / "bw-neg.yaml")
    data["default_format"] = "135"
    with pytest.raises(ValidationError, match="not one of the base formats"):
        Material.model_validate(data)


@pytest.mark.parametrize(
    ("path", "value", "match"),
    [
        (("optics", "dye_bands"), [], "at least 1"),
        (("optics", "psf_length"), _q(1.0, "um-1"), "not a length unit"),
        (("optics", "intrinsic_absorption", "edge_energy"), _q(2.7, "um"), "molar energy"),
        (("latent", "tau0_median"), _q(0.3, "um"), "not a time unit"),
        (("development", "covering_power"), _q(0.6, "g/m2"), "inverse areal mass"),
    ],
)
def test_emulsion_d1_fields_check_units(film, path, value, match):
    target = film["layers"][1]
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    with pytest.raises(ValidationError, match=match):
        Material.model_validate(film)


def test_dye_bands_at_most_four(film):
    bands = film["layers"][1]["optics"]["dye_bands"]
    bands.extend(copy.deepcopy(bands) + copy.deepcopy(bands[:1]))
    with pytest.raises(ValidationError, match="at most 4"):
        Material.model_validate(film)


def test_bw_mechanistic_problems_lists_missing_fields(film):
    emulsion = film["layers"][1]
    for name in ("optics", "latent", "development", "gelatin_coverage", "swelling_ratio"):
        del emulsion[name]
    del emulsion["grain"]["iodide_fraction"]
    del film["layers"][0]["thickness"]
    del film["layers"][2]["refractive_index"]
    del film["layers"][2]["formats"]
    del film["default_format"]
    problems = bw_mechanistic_problems(Material.model_validate(film))
    assert problems == [
        "layers[0].thickness: required by the bw-mechanistic model",
        "layers[1].gelatin_coverage: required by the bw-mechanistic model",
        "layers[1].swelling_ratio: required by the bw-mechanistic model",
        "layers[1].optics: required by the bw-mechanistic model",
        "layers[1].latent: required by the bw-mechanistic model",
        "layers[1].development: required by the bw-mechanistic model",
        "layers[1].grain.iodide_fraction: required by the bw-mechanistic model",
        "layers[2].refractive_index: required by the bw-mechanistic model",
        "layers[2].formats: required by the bw-mechanistic model",
    ]


def test_bw_mechanistic_problems_layer_structure():
    material = Material.model_validate(_yaml(FIXTURES / "films" / "test" / "bw-neg.yaml"))
    (problem,) = bw_mechanistic_problems(material)
    assert problem.startswith("layers: the bw-mechanistic model needs exactly the layers")
    assert "antihalation 'antihalation'" in problem


def test_emulsion_layer_unknown_d1_field_rejected():
    with pytest.raises(ValidationError, match="Extra inputs"):
        EmulsionLayer.model_validate({"id": "e", "optic": {}})


# --------------------------------------------------------------------------------------
# Species (spec D1 5.2)
# --------------------------------------------------------------------------------------


@pytest.fixture
def species_data() -> dict:
    return _yaml(SAMPLES / "species.yaml")


def test_d1_species_sample_loads(species_data):
    by_id = {s.id: s for s in SpeciesFile.model_validate(species_data).species}
    salt = by_id["test-sulfite-salt"]
    assert salt.kind == "ingredient" and salt.dissolves_to == {"test-sulfite": 1.0}
    assert salt.cations == {"Na": 2.0}
    sulfite = by_id["test-sulfite"]
    assert sulfite.acid_base.n_states == 3
    assert sulfite.acid_base.temperature_law == "van_t_hoff"
    agent = by_id["test-agent"].developing_agent
    assert (agent.role, agent.kinetic_class, agent.active_state) == (
        "regenerator",
        "hydroquinone",
        2,
    )
    assert by_id["test-gelatin"].immobile


def _species(species_data, sid):
    return next(s for s in species_data["species"] if s["id"] == sid)


@pytest.mark.parametrize(
    ("sid", "change", "match"),
    [
        ("test-sulfite", {"gel_diffusion_ratio": None}, "needs gel_diffusion_ratio"),
        ("test-sulfite", {"diffusion_coefficient": None}, "needs diffusion_coefficient"),
        ("test-sulfite", {"diffusion_coefficient": _q(1.0, "um")}, "not a diffusion coefficient"),
        ("test-sulfite", {"cations": {"Na": 1}}, "only allowed on species of kind 'ingredient'"),
        ("test-sulfite", {"toth": 1.0}, "toth: only allowed"),
        ("test-sulfite-salt", {"immobile": True}, "only allowed on species of kind 'component'"),
        ("test-sulfite-salt", {"gel_diffusion_ratio": _q(0.5)}, "only allowed"),
        ("test-sulfite-salt", {"cations": {"Li": 1}}, "Na"),
        ("test-sulfite-salt", {"dissolves_to": {"x": 0}}, "greater than 0"),
        ("test-agent", {"acid_base": None}, "needs acid_base"),
        ("test-agent", {"developing_agent": {"role": "single", "kinetic_class": "x",
                                             "active_state": 3, "product_sulfonated": False}},
         "out of range"),
        ("test-agent", {"developing_agent": {"role": "boss", "kinetic_class": "x",
                                             "active_state": 0, "product_sulfonated": False}},
         "role"),
    ],
)  # fmt: skip
def test_species_kind_rules(species_data, sid, change, match):
    entry = _species(species_data, sid)
    for key, value in change.items():
        if value is None:
            entry.pop(key, None)
        else:
            entry[key] = value
    with pytest.raises(ValidationError, match=match):
        Species.model_validate(entry)


@pytest.mark.parametrize(
    ("change", "match"),
    [
        ({"reference_protons": 3}, "must not exceed"),
        ({"ionisation_enthalpy": [_q(1.0, "kJ/mol")]}, "one entry per pka"),
        ({"temperature_law": "linear"}, "one dpka_dt per pka"),
        ({"dpka_dt": [_q(-0.01, "1/K"), _q(-0.01, "1/K")]}, "only used when"),
        ({"activity_class": "salt"}, "activity_class"),
    ],
)
def test_acid_base_rules(species_data, change, match):
    entry = _species(species_data, "test-sulfite")
    entry["acid_base"].update(change)
    with pytest.raises(ValidationError, match=match):
        Species.model_validate(entry)


def test_acid_base_linear_law(species_data):
    entry = _species(species_data, "test-sulfite")
    entry["acid_base"].update(
        temperature_law="linear", dpka_dt=[_q(-0.01, "1/K"), _q(-0.012, "1/K")]
    )
    assert len(Species.model_validate(entry).acid_base.dpka_dt) == 2


# --------------------------------------------------------------------------------------
# Kinetics (spec D1 5.3)
# --------------------------------------------------------------------------------------


@pytest.fixture
def kinetics_data() -> dict:
    return _yaml(SAMPLES / "kinetics.yaml")


def test_kinetics_sample_loads(kinetics_data):
    kin = KineticsFile.model_validate(kinetics_data)
    assert kin.agents["hydroquinone"].k_reg.unit == "L/(mol*s)"
    assert kin.agents["hydroquinone"].k_des is None
    assert kin.agents["metol"].k_des.to("1/s") == 0.1
    assert kin.halides.K_Br.to("L/mol") == 300.0
    assert kin.activity.delta_eps["sulfite"].bounds == (-0.45, 0.45)
    assert kin.transport.water_viscosity.viscosity.value == [1.053, 1.002, 0.911, 0.890]
    assert kin.gelatin.dry_density.to("g/cm3") == 1.35
    assert kin.solvent.sulfonation_switch.to("mmol/L") == pytest.approx(0.1)
    assert kin.calibration is None


@pytest.mark.parametrize(
    ("section", "key", "value", "match"),
    [
        ("halides", "K_Br", _q(1.0, "1/s"), "inverse concentration"),
        ("solvent", "k_diss", _q(1.0, "1/s"), "velocity"),
        ("transport", "swelling_half_time", _q(1.0, "um"), "time"),
        ("gelatin", "dry_density", _q(1.35, "g/L"), "mass density"),
        ("activity", "delta_eps", {"none": _q(0.0, "kg/mol")}, "sulfite"),
        ("barriers", "ionic_strength_ref", None, "Field required"),
        ("transport", "water_viscosity", {"temperature": _q([20, 18], "C"),
                                          "viscosity": _q([1.0, 1.05], "cP")}, "increasing"),
        ("transport", "water_viscosity", {"temperature": _q([18, 20], "C"),
                                          "viscosity": _q([1.0], "cP")}, "same length"),
        ("transport", "water_viscosity", {"temperature": _q(20, "C"),
                                          "viscosity": _q(1.0, "cP")}, "lists"),
    ],
)  # fmt: skip
def test_kinetics_rejects(kinetics_data, section, key, value, match):
    if value is None:
        del kinetics_data[section][key]
    else:
        kinetics_data[section][key] = value
    with pytest.raises(ValidationError, match=match):
        KineticsFile.model_validate(kinetics_data)


def test_kinetics_agent_units(kinetics_data):
    kinetics_data["agents"]["metol"]["activation_energy"] = _q(0.7, "eV")
    kin = KineticsFile.model_validate(kinetics_data)
    assert kin.agents["metol"].activation_energy.to("kJ/mol") == pytest.approx(67.54, abs=0.01)
    kinetics_data["agents"]["metol"]["k_reg"] = _q(1.0, "1/s")
    with pytest.raises(ValidationError, match="second-order rate"):
        KineticsFile.model_validate(kinetics_data)


# --------------------------------------------------------------------------------------
# Baths and processes (spec D1 5.5, 5.6)
# --------------------------------------------------------------------------------------


def _bath(**extra) -> dict:
    return {
        "schema_version": 1,
        "id": "test-dev",
        "name": "test developer",
        "kind": "developer",
        "formula_confidence": "msds-estimated",
        "components": [
            {"species": "a", "amount": _q(100.0, "g/L", bounds=[84.0, 119.0])},
            {"species": "b", "amount": _q(2.0, "g/L")},
        ],
        **extra,
    }


def test_bath_d1_fields():
    anchor = {
        "dilution": "1+1",
        "temperature": _q(25, "C"),
        "ph": _q(8.5, uncertainty=0.05),
        "method": "typical",
    }
    bath = Bath.model_validate(
        _bath(
            variant="packaged",
            base_solution="concentrate",
            ph_anchors=[anchor],
            nuisance=["a"],
            constants={"neutralisation_time": _q(2.0, "s")},
        )
    )
    assert bath.base_solution == "concentrate"
    assert bath.ph_anchors[0].dilution == "1+1"
    assert bath.constants.neutralisation_time.to("s") == 2.0
    assert bath.constants.wash_time_constant is None
    plain = Bath.model_validate(_bath())
    assert (plain.variant, plain.base_solution, plain.ph_anchors, plain.nuisance) == (
        None, "stock", [], []
    )  # fmt: skip
    assert plain.constants is None


@pytest.mark.parametrize(
    ("extra", "match"),
    [
        ({"nuisance": ["b"]}, "needs bounds"),
        ({"nuisance": ["c"]}, "not a component"),
        ({"nuisance": ["a", "a"]}, "duplicate"),
        ({"base_solution": "working"}, "base_solution"),
        ({"constants": {"wash_time_constant": _q(1.0, "um")}}, "not a time unit"),
        ({"constants": {"wash_time": _q(1.0, "s")}}, "Extra inputs"),
    ],
)
def test_bath_rejects(extra, match):
    with pytest.raises(ValidationError, match=match):
        Bath.model_validate(_bath(**extra))


@pytest.mark.parametrize(
    ("change", "match"),
    [
        ({"ph": _q(8.5)}, "needs an uncertainty"),
        ({"dilution": "1:1"}, "ambiguous"),
        ({"method": "guess"}, "method"),
    ],
)
def test_ph_anchor_rejects(change, match):
    anchor = {
        "dilution": "stock",
        "temperature": _q(25, "C"),
        "ph": _q(8.5, uncertainty=0.05),
        "method": "glass_electrode",
        **change,
    }
    with pytest.raises(ValidationError, match=match):
        PhAnchor.model_validate(anchor)


def test_process_d1_fields():
    assert Agitation(method="lift", initial_s=22.5, interval_s=60.0, duration_s=7.5).method == (
        "lift"
    )
    tank = Tank.model_validate(
        {"volume_ml": 3785.0, "film_area_cm2": 35.0, "image_area_fraction": _q(1.0)}
    )
    assert tank.image_area_fraction.to("1") == 1.0 and tank.rolls_before == 0
    assert Tank(volume_ml=500.0, rolls_before=3).rolls_before == 3
    for bad, match in (
        ({"image_area_fraction": _q(0.0)}, r"\(0, 1\]"),
        ({"image_area_fraction": _q(1.2)}, r"\(0, 1\]"),
        ({"rolls_before": -1}, "greater than or equal to 0"),
    ):
        with pytest.raises(ValidationError, match=match):
            Tank.model_validate({"volume_ml": 500.0, **bad})


# --------------------------------------------------------------------------------------
# Calibration records (spec D1 5.7)
# --------------------------------------------------------------------------------------


@pytest.fixture
def calibration_data() -> dict:
    return _yaml(SAMPLES / "calibration.yaml")


def test_calibration_record_loads(calibration_data):
    record = CalibrationRecord.model_validate(calibration_data)
    assert record.created_utc.tzinfo is not None  # PyYAML reads the trailing Z as UTC
    assert record.stages[0].name == "C_film"
    assert record.stages[0].optimizer["starts"] == 8
    assert record.parameters[0].interval is None
    assert record.singular_values == {"C_film": [3.0, 0.5]}
    assert not record.frozen and not record.held_out_evaluated and record.notes == []


@pytest.mark.parametrize(
    ("change", "match"),
    [
        ({"created_utc": "2026-12-01T00:00:00"}, "timezone"),
        ({"profile": "image"}, "calibration"),
        ({"training_set_sha256": "PLACEHOLDER" + "0" * 55}, "hexadecimal"),
        ({"held_out_evaluated": True}, "requires frozen"),
        ({"stages": [{"name": "G_extra"}]}, "A_activity"),
    ],
)
def test_calibration_record_rejects(calibration_data, change, match):
    calibration_data.update(change)
    with pytest.raises(ValidationError, match=match):
        CalibrationRecord.model_validate(calibration_data)


def test_calibration_parameter_interval_required_key(calibration_data):
    del calibration_data["parameters"][0]["interval"]
    with pytest.raises(ValidationError, match="Field required"):
        CalibrationRecord.model_validate(calibration_data)
    calibration_data["parameters"][0]["interval"] = [0.7, 0.9]
    calibration_data["parameters"].append(dict(calibration_data["parameters"][0]))
    with pytest.raises(ValidationError, match="duplicate parameter"):
        CalibrationRecord.model_validate(calibration_data)


# --------------------------------------------------------------------------------------
# Walking quantities
# --------------------------------------------------------------------------------------


def test_iter_quantities_paths(film):
    material = Material.model_validate(film)
    paths = dict(iter_quantities(material))
    assert "layers[1].grain.mean_diameter" in paths
    assert "layers[1].optics.dye_bands[1].width" in paths
    assert "layers[2].formats.120.backing_density" in paths
    assert paths["iso_speed"].value == 100.0
