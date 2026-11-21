import hashlib
import re
from pathlib import Path

import numpy as np
import pytest
import yaml

from darkroom_atlas.data.library import DataError, DataNotFound, Library, default_data_root
from darkroom_atlas.data.models import (
    BaseLayer,
    Bath,
    FilterLayer,
    Material,
    Process,
    ReferenceCurve,
    SpectralCurve,
)
from darkroom_atlas.data.provenance import CurveOrigin, Source
from darkroom_atlas.spectral import DEFAULT_GRID, SpectralGrid, cmf_xyz

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "library"
REPO_DATA = Path(__file__).resolve().parents[2] / "data"


@pytest.fixture
def fx():
    return Library(FIXTURES)


@pytest.fixture
def repo(monkeypatch):
    monkeypatch.delenv("DARKROOM_ATLAS_DATA", raising=False)
    return Library()


def _write_yaml(path: Path, data) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return path


def _material_data(id_="a", **overrides):
    data = {
        "schema_version": 1,
        "id": id_,
        "maker": "test",
        "name": "test-a",
        "kind": "ideal",
        "model": "ideal",
        "layers": [{"type": "emulsion", "id": "e"}],
    }
    data.update(overrides)
    return data


def _process_data(id_="p", steps=None):
    return {
        "schema_version": 1,
        "id": id_,
        "name": "test-p",
        "family": "bw",
        "steps": steps if steps is not None else [{"kind": "identity"}],
    }


def _bath_data(id_="d", components=()):
    return {
        "schema_version": 1,
        "id": id_,
        "name": "test-d",
        "kind": "developer",
        "formula_confidence": "published",
        "components": list(components),
    }


# --------------------------------------------------------------------------------------
# Whole libraries
# --------------------------------------------------------------------------------------


def test_repo_library_valid(repo):
    assert repo.validate() == []


def test_fixture_library_valid(fx):
    assert fx.validate() == []


def test_load_ideal(repo):
    material = repo.material("ideal")
    assert isinstance(material, Material)
    assert material.model == "ideal"
    assert material.default_readout == "ideal_colorimetric"
    assert [layer.id for layer in material.emulsion_layers] == ["x", "y", "z"]
    assert len(material.emulsion_layers) == 3
    for layer in material.emulsion_layers:
        assert layer.spectral_sensitivity.source is Source.ASSUMED
        assert layer.spectral_sensitivity.builtin == f"cie1931_{layer.id}"


def test_load_repo_identity(repo):
    proc = repo.process("identity")
    assert isinstance(proc, Process)
    assert [step.kind.value for step in proc.steps] == ["identity"]
    # The repository species file grows with D1c; it must load whatever it holds.
    assert all(s.id == key for key, s in repo.species().items())


def test_load_color_negative(fx):
    material = fx.material("test/color-neg")
    assert material.id == "test/color-neg"
    emulsions = material.emulsion_layers
    assert [layer.id for layer in emulsions] == ["blue", "green", "red"]
    assert len([layer for layer in material.layers if isinstance(layer, FilterLayer)]) == 1
    assert len([layer for layer in material.layers if isinstance(layer, BaseLayer)]) == 1
    assert [layer.type for layer in material.layers] == [
        "overcoat", "emulsion", "filter", "emulsion", "interlayer", "emulsion",
        "antihalation", "base",
    ]  # fmt: skip
    image_dyes = [[c.dye for c in layer.couplers if c.kind == "image"] for layer in emulsions]
    assert image_dyes == [["yellow"], ["magenta"], ["cyan"]]
    kinds = {c.kind for layer in emulsions for c in layer.couplers}
    assert kinds == {"image", "dir", "colored_mask"}
    green = emulsions[1]
    assert {c.kind for c in green.couplers} == {"image", "dir", "colored_mask"}
    for layer in emulsions:
        assert layer.spectral_sensitivity.inline is not None


def test_load_bw_negative_and_baths(fx):
    material = fx.material("test/bw-neg")
    (emulsion,) = material.emulsion_layers
    assert emulsion.grain.morphology == "conventional"
    assert emulsion.grain.mean_diameter.to("nm") == pytest.approx(500.0)
    bath = fx.bath("test-bw-dev")
    assert isinstance(bath, Bath)
    assert [c.species for c in bath.components] == ["test-dev-agent", "test-preservative"]
    assert fx.bath("test-bleach").kind.value == "bleach"
    species = fx.species()
    assert sorted(species) == ["test-cd-agent", "test-dev-agent", "test-preservative"]
    proc = fx.process("test-bw")
    develop = proc.steps[0]
    assert develop.time_s == 585.0 and develop.dilution_fraction == 0.5
    assert develop.agitation.method == "inversion" and develop.tank.volume_ml == 500.0


def test_ids_listing(fx):
    assert fx.ids("material") == ["test/bw-mechanistic", "test/bw-neg", "test/color-neg"]
    assert fx.ids("bath") == ["test-bleach", "test-bw-dev", "test-color-dev", "test-fix"]
    assert fx.ids("process") == ["test-bw", "test-c41"]
    with pytest.raises(ValueError, match="material"):
        fx.ids("film")


def test_path_of(fx):
    assert fx.path_of("material", "test/bw-neg") == FIXTURES / "films" / "test" / "bw-neg.yaml"
    assert fx.path_of("bath", "test-fix") == FIXTURES / "baths" / "test-fix.yaml"
    assert fx.path_of("bath", "test-bw-dev") == FIXTURES / "developers" / "test-bw-dev.yaml"
    assert fx.path_of("process", "test-c41") == FIXTURES / "processes" / "test-c41.yaml"
    with pytest.raises(DataNotFound):
        fx.path_of("process", "nope")


def test_not_found_lists_available(fx):
    with pytest.raises(DataNotFound) as info:
        fx.material("nope")
    message = str(info.value)
    assert "test/bw-neg" in message and "test/color-neg" in message
    assert isinstance(info.value, DataError)
    with pytest.raises(DataNotFound, match="test-bw-dev"):
        fx.bath("nope")
    with pytest.raises(DataNotFound, match="test-c41"):
        fx.process("nope")
    with pytest.raises(DataNotFound, match="test/bw-neg/characteristic.csv"):
        fx.reference_curve("nope.csv")


@pytest.mark.parametrize("bad_id", ["../films/test/bw-neg", "/etc/passwd", "", "a//b"])
def test_not_found_for_malformed_ids(fx, bad_id):
    with pytest.raises(DataNotFound):
        fx.material(bad_id)


def test_papers_fallback(tmp_path):
    _write_yaml(tmp_path / "papers" / "test" / "p.yaml", _material_data("test/p"))
    lib = Library(tmp_path)
    assert lib.material("test/p").id == "test/p"
    assert lib.ids("material") == ["test/p"]
    assert lib.validate() == []


def test_id_must_match_path(tmp_path):
    path = _write_yaml(tmp_path / "films" / "a.yaml", _material_data("b"))
    lib = Library(tmp_path)
    with pytest.raises(DataError) as info:
        lib.material("a")
    message = str(info.value)
    assert message.startswith(str(path))
    assert "'b'" in message and "'a'" in message
    problems = lib.validate()
    assert len(problems) == 1 and problems[0].startswith(str(path))


def test_unknown_field_reports_path(tmp_path):
    data = _material_data()
    data["layers"][0]["foo"] = 1
    path = _write_yaml(tmp_path / "films" / "a.yaml", data)
    with pytest.raises(DataError) as info:
        Library(tmp_path).material("a")
    message = str(info.value)
    assert message.startswith(str(path))
    assert "foo" in message
    assert "Traceback" not in message
    assert info.value.__cause__ is None


def test_missing_ref_reports_path(tmp_path):
    data = _material_data()
    data["layers"][0]["thickness"] = {"value": 10.0, "unit": "um", "source": "datasheet"}
    path = _write_yaml(tmp_path / "films" / "a.yaml", data)
    with pytest.raises(DataError) as info:
        Library(tmp_path).material("a")
    message = str(info.value)
    assert message.startswith(str(path))
    assert "ref" in message and "thickness" in message


def test_each_problem_on_its_own_line(tmp_path):
    data = _material_data(kind="nope", extra_field=1)
    path = _write_yaml(tmp_path / "films" / "a.yaml", data)
    with pytest.raises(DataError) as info:
        Library(tmp_path).material("a")
    lines = str(info.value).splitlines()
    assert len(lines) >= 2
    assert all(line.startswith(f"{path}: ") for line in lines)


@pytest.mark.parametrize(
    ("text", "fragment"),
    [
        ("", "empty"),
        ("- 1\n- 2\n", "mapping"),
        ("id: [unclosed\n", "YAML"),
        ("id: a\nid: a\n", "duplicate key"),
    ],
)
def test_malformed_yaml_reports_path(tmp_path, text, fragment):
    path = tmp_path / "films" / "a.yaml"
    path.parent.mkdir(parents=True)
    path.write_text(text, encoding="utf-8")
    with pytest.raises(DataError) as info:
        Library(tmp_path).material("a")
    assert str(info.value).startswith(str(path))
    assert fragment in str(info.value)


def test_missing_root_rejected(tmp_path):
    with pytest.raises(DataNotFound, match="does not exist"):
        Library(tmp_path / "nope")


# --------------------------------------------------------------------------------------
# validate(): collects everything, never stops at the first error
# --------------------------------------------------------------------------------------


def test_validate_dangling_bath(tmp_path):
    steps = [{"kind": "develop", "bath": "nope-bath", "time": "5:00"}]
    path = _write_yaml(tmp_path / "processes" / "p.yaml", _process_data(steps=steps))
    problems = Library(tmp_path).validate()
    assert len(problems) == 1
    assert "nope-bath" in problems[0]
    assert problems[0].startswith(str(path))


def test_validate_dangling_species(tmp_path):
    amount = {"value": 1.0, "unit": "g/L", "source": "assumed"}
    components = [{"species": "nope-species", "amount": amount}]
    path = _write_yaml(tmp_path / "developers" / "d.yaml", _bath_data(components=components))
    problems = Library(tmp_path).validate()
    assert len(problems) == 1
    assert "nope-species" in problems[0]
    assert problems[0].startswith(str(path))


def test_validate_missing_curve_file(tmp_path):
    data = _material_data()
    data["layers"][0]["spectral_sensitivity"] = {
        "quantity": "relative_energy_sensitivity",
        "source": "assumed",
        "file": "curves/missing.csv",
    }
    path = _write_yaml(tmp_path / "films" / "a.yaml", data)
    problems = Library(tmp_path).validate()
    assert len(problems) == 1
    assert "curves/missing.csv" in problems[0]
    assert problems[0].startswith(str(path))


def test_validate_collects_all_problems(tmp_path):
    bad_material = _material_data(foo=1)
    _write_yaml(tmp_path / "films" / "a.yaml", bad_material)
    _write_yaml(tmp_path / "films" / "b.yaml", _material_data("c"))
    steps = [{"kind": "fix", "bath": "nope-bath", "time": "5:00"}]
    _write_yaml(tmp_path / "processes" / "p.yaml", _process_data(steps=steps))
    (tmp_path / "chemistry").mkdir()
    (tmp_path / "chemistry" / "species.yaml").write_text("species: [\n", encoding="utf-8")
    problems = Library(tmp_path).validate()
    joined = "\n".join(problems)
    assert "a.yaml" in joined and "foo" in joined
    assert "b.yaml" in joined
    assert "nope-bath" in joined
    assert "species.yaml" in joined
    assert len(problems) == 4


def test_validate_reports_ambiguous_ids_and_yml(tmp_path):
    _write_yaml(tmp_path / "films" / "a.yaml", _material_data("a"))
    _write_yaml(tmp_path / "papers" / "a.yaml", _material_data("a"))
    _write_yaml(tmp_path / "processes" / "q.yml", _process_data("q"))
    problems = Library(tmp_path).validate()
    assert len(problems) == 2
    assert any("papers" in p and "films" in p for p in problems)
    assert any("q.yml" in p and ".yaml" in p for p in problems)


def test_validate_reference_material_must_exist(tmp_path):
    csv = tmp_path / "reference" / "test" / "gone" / "c.csv"
    csv.parent.mkdir(parents=True)
    csv.write_text(
        "# material: test/gone\n# quantity: q\n# x_unit: 1\n# y_unit: 1\n"
        "# origin: datasheet\n# ref: synthetic\n0,1\n1,2\n",
        encoding="utf-8",
    )
    problems = Library(tmp_path).validate()
    assert len(problems) == 1
    assert "test/gone" in problems[0] and problems[0].startswith(str(csv))


# --------------------------------------------------------------------------------------
# Curves
# --------------------------------------------------------------------------------------


def test_resolve_builtin(fx):
    for index, name in enumerate(["cie1931_x", "cie1931_y", "cie1931_z"]):
        curve = SpectralCurve(quantity="q", source="assumed", builtin=name)
        out = fx.resolve_curve(curve, DEFAULT_GRID)
        np.testing.assert_array_equal(out, cmf_xyz(DEFAULT_GRID)[:, index])
    grid = SpectralGrid(400.0, 700.0, 10.0)
    curve = SpectralCurve(quantity="q", source="assumed", builtin="cie1931_y")
    assert fx.resolve_curve(curve, grid).shape == (grid.n_bands,)


def test_resolve_inline_zero_outside(fx):
    curve = SpectralCurve(
        quantity="q",
        source="assumed",
        inline={"wavelengths_nm": [400, 550, 700], "values": [1.0, 2.0, 3.0]},
    )
    out = fx.resolve_curve(curve, DEFAULT_GRID)
    wl = DEFAULT_GRID.wavelengths
    assert out.shape == (81,)
    assert out[0] == 0.0  # 380 nm
    assert out[wl == 395.0][0] == 0.0
    assert out[wl == 400.0][0] == 1.0
    assert out[wl == 475.0][0] == pytest.approx(1.5, abs=1e-12)
    assert out[wl == 700.0][0] == 3.0
    assert out[-1] == 0.0  # 780 nm


def test_curve_files_lists_every_named_file(fx):
    def curve(**given):
        return {"quantity": "q", "source": "assumed", **given}

    material = Material.model_validate(
        _material_data(
            layers=[
                {
                    "type": "emulsion",
                    "id": "e",
                    "spectral_sensitivity": curve(file="./curves/s.csv"),
                    "couplers": [
                        {
                            "id": "c",
                            "dye": "cyan",
                            "kind": "image",
                            "dye_spectral_density": curve(file="curves/dye.csv"),
                        }
                    ],
                },
                {"type": "emulsion", "id": "f", "spectral_sensitivity": curve(file="curves/s.csv")},
                {"type": "filter", "id": "y", "spectral_density": curve(builtin="cie1931_y")},
                {
                    "type": "base",
                    "id": "b",
                    "material": "PET",
                    "spectral_density": curve(file="curves//base.csv"),
                },
            ]
        )
    )
    assert Library.curve_files(material) == ["curves/base.csv", "curves/dye.csv", "curves/s.csv"]
    assert fx.curve_files(fx.material("test/bw-neg")) == []


def test_resolve_file_curve(tmp_path):
    csv = tmp_path / "curves" / "s.csv"
    csv.parent.mkdir(parents=True)
    csv.write_text("# quantity: q\n400,1\n500,2\n600,3\n", encoding="utf-8")
    lib = Library(tmp_path)
    curve = SpectralCurve(quantity="q", source="assumed", file="curves/s.csv")
    out = lib.resolve_curve(curve, DEFAULT_GRID)
    wl = DEFAULT_GRID.wavelengths
    assert out[0] == 0.0
    assert out[wl == 450.0][0] == pytest.approx(1.5, abs=1e-12)
    missing = SpectralCurve(quantity="q", source="assumed", file="curves/none.csv")
    with pytest.raises(DataError, match="none.csv"):
        lib.resolve_curve(missing, DEFAULT_GRID)


def test_reference_curve(fx):
    curve = fx.reference_curve("test/bw-neg/characteristic.csv")
    assert isinstance(curve, ReferenceCurve)
    assert curve.origin == CurveOrigin.DATASHEET
    assert curve.material == "test/bw-neg"
    assert curve.x_unit == "log10(lux*s)"
    assert curve.conditions == "test-bw-dev 1+1, 9:45 at 20 C"
    assert len(curve.x) == len(curve.y) == 7
    assert curve.x[0] == -3.0 and curve.y[-1] == 1.8


@pytest.mark.parametrize(
    ("text", "fragment"),
    [
        ("# material: test/a\n0,1\n1,2\n", "Field required"),
        (
            (
                "# material: test/a\n# quantity: q\n# x_unit: 1\n# y_unit: 1\n"
                "# origin: datasheet\n# ref: r\n# colour: red\n0,1\n1,2\n"
            ),
            "colour",
        ),
        (
            (
                "# material: test/a\n# quantity: q\n# x_unit: 1\n# y_unit: 1\n"
                "# origin: datasheet\n# ref: r\n0,1\n1,x\n"
            ),
            "line 8",
        ),
        (
            (
                "# material: test/b\n# quantity: q\n# x_unit: 1\n# y_unit: 1\n"
                "# origin: datasheet\n# ref: r\n0,1\n1,2\n"
            ),
            "'test/a'",
        ),
    ],
)
def test_reference_curve_errors_report_path(tmp_path, text, fragment):
    csv = tmp_path / "reference" / "test" / "a" / "c.csv"
    csv.parent.mkdir(parents=True)
    csv.write_text(text, encoding="utf-8")
    with pytest.raises(DataError) as info:
        Library(tmp_path).reference_curve("test/a/c.csv")
    assert str(info.value).startswith(str(csv))
    assert fragment in str(info.value)


# --------------------------------------------------------------------------------------
# Hashing and the data root
# --------------------------------------------------------------------------------------


def test_sha256_stable(fx):
    path = fx.path_of("material", "test/bw-neg")
    first = fx.sha256(path)
    second = fx.sha256(path)
    assert first == second
    assert re.fullmatch(r"[0-9a-f]{64}", first)
    assert first == hashlib.sha256(path.read_bytes()).hexdigest()


def test_env_override(monkeypatch):
    monkeypatch.setenv("DARKROOM_ATLAS_DATA", str(FIXTURES))
    assert default_data_root() == FIXTURES
    assert Library().root == FIXTURES
    monkeypatch.delenv("DARKROOM_ATLAS_DATA")
    assert default_data_root() == REPO_DATA
