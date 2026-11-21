"""D1 library accessors and checks (spec D1 5.2, 5.3, 5.7, 5.9)."""

import shutil
from pathlib import Path

import pytest
import yaml
from typer.testing import CliRunner

from darkroom_atlas.cli import app
from darkroom_atlas.data.library import (
    DEVELOPER_FACTS_FILES,
    DataError,
    DataNotFound,
    Library,
)
from darkroom_atlas.data.models import CalibrationRecord, KineticsFile

TESTS = Path(__file__).resolve().parents[1]
FIXTURES = TESTS / "fixtures" / "library"
SAMPLES = Path(__file__).resolve().parent / "samples"


@pytest.fixture
def root(tmp_path) -> Path:
    """A copy of the fixture library with the D1 sample files added."""
    lib = tmp_path / "library"
    shutil.copytree(FIXTURES, lib)
    shutil.copy(SAMPLES / "kinetics.yaml", lib / "chemistry" / "kinetics.yaml")
    species = yaml.safe_load((lib / "chemistry" / "species.yaml").read_text())
    species["species"] += yaml.safe_load((SAMPLES / "species.yaml").read_text())["species"]
    _write(lib / "chemistry" / "species.yaml", species)
    (lib / "calibration").mkdir()
    shutil.copy(SAMPLES / "calibration.yaml", lib / "calibration" / "d1e-test.yaml")
    _write(
        lib / "reference" / "test" / "bw-neg" / "tables.yaml",
        {"schema_version": 1, "material": "test/bw-neg", "facts": [{"id": "x", "value": 1}]},
    )
    _write(
        lib / "reference" / "developers" / "hc-110.yaml",
        {"schema_version": 1, "developer": "hc-110", "facts": [{"id": "ph", "value": None}]},
    )
    return lib


def _write(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")


def _edit(path: Path, fn) -> None:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    fn(data)
    _write(path, data)


def test_d1_library_valid(root):
    assert Library(root).validate() == []


def test_fixture_library_with_bw_mechanistic_passes_check_data():
    result = CliRunner().invoke(app, ["check-data", "--data", str(FIXTURES)])
    assert result.exit_code == 0, result.output
    assert "OK:" in result.output and "3 materials" in result.output


def test_fixture_film_loads_through_library():
    material = Library(FIXTURES).material("test/bw-mechanistic")
    assert material.model == "bw-mechanistic"
    assert material.formats["120"].backing_density.to("1") == 0.5


# ---- kinetics -------------------------------------------------------------------------


def test_kinetics_accessor(root):
    kinetics = Library(root).kinetics()
    assert isinstance(kinetics, KineticsFile)
    assert sorted(kinetics.agents) == ["hydroquinone", "metol"]


def test_kinetics_missing():
    with pytest.raises(DataNotFound, match="kinetics file does not exist"):
        Library(FIXTURES).kinetics()


def test_kinetics_invalid_reported(root):
    _edit(root / "chemistry" / "kinetics.yaml", lambda d: d["halides"].pop("K_Br"))
    with pytest.raises(DataError, match=r"kinetics.yaml: halides.K_Br: Field required"):
        Library(root).kinetics()
    assert any("halides.K_Br: Field required" in p for p in Library(root).validate())


def test_validate_unknown_kinetic_class(root):
    def change(data):
        agent = next(s for s in data["species"] if s["id"] == "test-agent")
        agent["developing_agent"]["kinetic_class"] = "ascorbate"

    _edit(root / "chemistry" / "species.yaml", change)
    problems = Library(root).validate()
    assert len(problems) == 1
    assert "developing_agent.kinetic_class: unknown agent class 'ascorbate'" in problems[0]
    assert "hydroquinone, metol" in problems[0]


def test_validate_developing_agent_needs_kinetics_file(root):
    (root / "chemistry" / "kinetics.yaml").unlink()
    problems = Library(root).validate()
    assert len(problems) == 1
    assert "kinetic_class: needs" in problems[0] and "does not exist" in problems[0]


def test_validate_dissolves_to_must_name_components(root):
    def change(data):
        salt = next(s for s in data["species"] if s["id"] == "test-sulfite-salt")
        salt["dissolves_to"] = {"test-preservative": 1, "nothing": 2}

    _edit(root / "chemistry" / "species.yaml", change)
    problems = Library(root).validate()
    assert len(problems) == 2
    assert all("is not a component species" in p for p in problems)
    assert "dissolves_to.test-preservative" in problems[0]
    assert "dissolves_to.nothing" in problems[1]
    assert "species[3]" in problems[0]  # index in the file's species list


def test_validate_component_without_diffusion(root):
    def change(data):
        sulfite = next(s for s in data["species"] if s["id"] == "test-sulfite")
        del sulfite["gel_diffusion_ratio"]

    _edit(root / "chemistry" / "species.yaml", change)
    problems = Library(root).validate()
    assert any("a mobile component needs gel_diffusion_ratio" in p for p in problems)


# ---- bw-mechanistic films --------------------------------------------------------------


def test_validate_incomplete_bw_mechanistic_film(root):
    path = root / "films" / "test" / "bw-mechanistic.yaml"
    _edit(path, lambda d: d["layers"][1].pop("latent"))
    assert Library(root).validate() == [
        f"{path}: layers[1].latent: required by the bw-mechanistic model"
    ]


def test_other_models_not_checked_for_bw_fields():
    # test/bw-neg uses model test-bw and lacks the D1 fields; it stays valid.
    assert Library(FIXTURES).validate() == []


# ---- calibration records ---------------------------------------------------------------


def test_calibration_record_accessor(root):
    record = Library(root).calibration_record("d1e-test")
    assert isinstance(record, CalibrationRecord)
    assert record.stages[0].data_files[0].path == "reference/test/bw-neg/characteristic.csv"


@pytest.mark.parametrize("bad", ["missing", "../d1e-test", "", "a/b"])
def test_calibration_record_not_found(root, bad):
    with pytest.raises(DataNotFound, match="available: d1e-test"):
        Library(root).calibration_record(bad)


def test_calibration_record_run_id_must_match_name(root):
    shutil.copy(root / "calibration" / "d1e-test.yaml", root / "calibration" / "other.yaml")
    with pytest.raises(DataError, match="declares 'd1e-test' but its name implies 'other'"):
        Library(root).calibration_record("other")
    assert any("its name implies 'other'" in p for p in Library(root).validate())


def _fit(run_id: str) -> dict:
    return {
        "run_id": run_id,
        "model_version": "D1e",
        "training_set_sha256": "0" * 64,
        "method": "covariance",
        "transform": "log10",
        "sigma_transformed": 0.02,
    }


def test_validate_fitted_values_need_calibration_record(root):
    film = root / "films" / "test" / "bw-mechanistic.yaml"

    def fitted(run_id):
        def change(data):
            q = data["layers"][1]["grain"]["mean_diameter"]
            q.update(source="fitted", fit=_fit(run_id))

        return change

    _edit(film, fitted("d1e-test"))
    assert Library(root).validate() == []
    _edit(film, fitted("d1e-gone"))
    expected = (
        f"{film}: layers[1].grain.mean_diameter.fit.run_id: no calibration record 'd1e-gone' "
        f"in {root / 'calibration'}"
    )
    assert Library(root).validate() == [expected]


def test_validate_fitted_kinetics_need_calibration_record(root):
    def change(data):
        data["halides"]["K_Br"].update(source="fitted", fit=_fit("nope"))

    _edit(root / "chemistry" / "kinetics.yaml", change)
    (problem,) = Library(root).validate()
    assert "kinetics.yaml: halides.K_Br.fit.run_id: no calibration record 'nope'" in problem


# ---- reference tables and developer facts -----------------------------------------------


def test_reference_tables_accessor(root):
    tables = Library(root).reference_tables("test/bw-neg")
    assert tables["facts"] == [{"id": "x", "value": 1}]


def test_reference_tables_not_found(root):
    with pytest.raises(DataNotFound, match="available: test/bw-neg"):
        Library(root).reference_tables("test/color-neg")
    with pytest.raises(DataNotFound):
        Library(root).reference_tables("../chemistry")


def test_reference_tables_material_mismatch(root):
    path = root / "reference" / "test" / "bw-neg" / "tables.yaml"
    _edit(path, lambda d: d.update(material="test/color-neg"))
    with pytest.raises(DataError, match="declares 'test/color-neg' but its path implies"):
        Library(root).reference_tables("test/bw-neg")
    assert len(Library(root).validate()) == 1


def test_reference_tables_schema_version(root):
    path = root / "reference" / "test" / "bw-neg" / "tables.yaml"
    _edit(path, lambda d: d.update(schema_version=2))
    with pytest.raises(DataError, match="schema_version: expected 1"):
        Library(root).reference_tables("test/bw-neg")


def test_developer_facts_accessor(root):
    facts = Library(root).developer_facts("hc-110-pre2019")
    assert facts["developer"] == "hc-110"


def test_developer_facts_mapping():
    assert DEVELOPER_FACTS_FILES == {
        "d-76": "d-76",
        "xtol": "xtol",
        "hc-110-pre2019": "hc-110",
        "t-max-dev": "t-max-dev",
    }


@pytest.mark.parametrize("bad", ["hc-110-2019", "d-76", "hc-110", "t-max-rs", "rodinal"])
def test_developer_facts_not_found(root, bad):
    with pytest.raises(DataNotFound, match="available: hc-110-pre2019"):
        Library(root).developer_facts(bad)


def test_developer_facts_name_mismatch(root):
    _edit(root / "reference" / "developers" / "hc-110.yaml", lambda d: d.update(developer="xtol"))
    with pytest.raises(DataError, match="declares 'xtol'"):
        Library(root).developer_facts("hc-110-pre2019")
    assert len(Library(root).validate()) == 1


def test_validate_checks_unmapped_developer_facts(root):
    path = root / "reference" / "developers" / "t-max-rs.yaml"
    _write(path, {"developer": "t-max-rs"})
    assert Library(root).validate() == []
    with pytest.raises(DataNotFound):
        Library(root).developer_facts("t-max-rs")
    _write(path, {"developer": "t-max"})
    (problem,) = Library(root).validate()
    assert "t-max-rs.yaml: developer: the file declares 't-max'" in problem


def test_developer_facts_do_not_disturb_reference_curves(root):
    lib = Library(root)
    assert lib.reference_curve("test/bw-neg/characteristic.csv").material == "test/bw-neg"


@pytest.mark.skipif(
    not (TESTS.parent / "tools" / "digitize" / "staged" / "developers").is_dir(),
    reason="staged digitising data not present",
)
def test_staged_reference_files_readable(tmp_path):
    """The digitising team's staged files read through the accessors once merged (5.9)."""
    staged = TESTS.parent / "tools" / "digitize" / "staged"
    lib = tmp_path / "library"
    shutil.copytree(FIXTURES, lib)
    shutil.copytree(staged / "developers", lib / "reference" / "developers")
    (lib / "reference" / "developers" / "README.md").unlink()
    for film in ("tri-x-400", "t-max-100"):
        target = lib / "reference" / "kodak" / film
        target.mkdir(parents=True)
        shutil.copy(staged / "kodak" / film / "tables.yaml", target / "tables.yaml")
    library = Library(lib)
    assert library.validate() == []
    for bath in DEVELOPER_FACTS_FILES:
        assert library.developer_facts(bath)["developer"] == DEVELOPER_FACTS_FILES[bath]
    for film in ("tri-x-400", "t-max-100"):
        assert library.reference_tables(f"kodak/{film}")["material"] == f"kodak/{film}"


def test_validate_fitted_species_values_use_species_paths(root):
    def change(data):
        sulfite = next(s for s in data["species"] if s["id"] == "test-sulfite")
        sulfite["acid_base"]["pka"][1].update(source="fitted", fit=_fit("nope"))

    _edit(root / "chemistry" / "species.yaml", change)
    (problem,) = Library(root).validate()
    assert "species.yaml: species[4].acid_base.pka[1].fit.run_id: no calibration record" in problem


# ---- revised checks (spec D1 5.2, 5.7) -------------------------------------------------


@pytest.mark.parametrize("builtin", ["toth", "gelatin"])
def test_validate_rejects_builtin_pseudo_components_in_species_file(root, builtin):
    def add(data):
        entry = {
            "id": builtin,
            "name": builtin,
            "molar_mass": {"value": 1.0, "unit": "g/mol", "source": "assumed"},
            "roles": ["buffer"],
            "kind": "component",
            "immobile": True,
        }
        data["species"].append(entry)

    _edit(root / "chemistry" / "species.yaml", add)
    (problem,) = Library(root).validate()
    assert f"species[{_species_count(root) - 1}].id: {builtin!r} is a built-in" in problem


def _species_count(root: Path) -> int:
    return len(yaml.safe_load((root / "chemistry" / "species.yaml").read_text())["species"])


def _pair(material: str, developer: str) -> dict:
    return {
        "material_id": material,
        "developer_id": developer,
        "tier": "training",
        "status": "validated",
        "thresholds_sha256": "a" * 64,
    }


def test_validate_pair_statuses_name_known_files(root):
    path = root / "calibration" / "d1e-test.yaml"
    _edit(path, lambda d: d.update(pair_statuses=[_pair("test/bw-neg", "test-bw-dev")]))
    assert Library(root).validate() == []
    _edit(path, lambda d: d.update(pair_statuses=[_pair("test/nope", "nope-dev")]))
    problems = Library(root).validate()
    assert problems == [
        f"{path}: pair_statuses[0].material_id: unknown material 'test/nope'; available: "
        + ", ".join(Library(root).ids("material")),
        f"{path}: pair_statuses[0].developer_id: unknown bath 'nope-dev'; available: "
        + ", ".join(Library(root).ids("bath")),
    ]
