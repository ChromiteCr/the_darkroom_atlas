import os
import shutil
import subprocess
from pathlib import Path

import pytest
from typer.testing import CliRunner

from darkroom_atlas.cli import app
from darkroom_atlas.data.library import ENV_DATA_ROOT
from darkroom_atlas.scene.io import EIGHT_BIT_WARNING

TESTS = Path(__file__).resolve().parent
REPO = TESTS.parent
SAMPLE = TESTS / "assets" / "sample.png"
OUTPUT_NAME = "sample__ideal__identity.png"


@pytest.fixture
def cli(monkeypatch):
    monkeypatch.delenv(ENV_DATA_ROOT, raising=False)
    return CliRunner()


def _run_args(out_dir, *extra, material="ideal", process="identity", image=SAMPLE):
    return ["run", str(image), "-m", material, "-p", process, "-o", str(out_dir), *extra]


def _assert_one_line_error(result, code=2):
    """Exit code ``code`` through a clean exit, one message line on stderr, no traceback."""
    assert result.exit_code == code, result.output
    assert isinstance(result.exception, SystemExit), repr(result.exception)
    assert "Traceback" not in result.output
    lines = result.stderr.strip().splitlines()
    assert len(lines) == 1, result.stderr
    assert lines[0].startswith("error: ")
    return lines[0]


def test_run_ok(tmp_path, cli):
    out_dir = tmp_path / "out"
    result = cli.invoke(app, _run_args(out_dir))
    assert result.exit_code == 0, result.output
    assert OUTPUT_NAME in result.stdout
    assert "sample__ideal__identity.run.json" in result.stdout
    assert EIGHT_BIT_WARNING in result.stderr
    assert (out_dir / OUTPUT_NAME).is_file()
    assert (out_dir / "sample__ideal__identity.run.json").is_file()


def test_run_options_reach_manifest(tmp_path, cli):
    from darkroom_atlas.manifest import RunManifest

    out_dir = tmp_path / "out"
    extra = ["--iso", "400", "--ev", "1", "--time", "0.5", "--seed", "7", "--tile", "16",
             "--long-edge-mm", "24", "--illuminant", "A", "-r", "ideal_colorimetric"]  # fmt: skip
    result = cli.invoke(app, _run_args(out_dir, *extra))
    assert result.exit_code == 0, result.output
    manifest = RunManifest.read(out_dir / "sample__ideal__identity.run.json")
    assert manifest.exposure["iso"] == 400
    assert manifest.exposure["ev"] == 1
    assert manifest.exposure["time_s"] == 0.5
    assert manifest.exposure["illuminant"] == "A"
    assert manifest.seed == 7
    assert manifest.tile_px == 16
    assert manifest.geometry["long_edge_mm"] == 24
    assert manifest.readout.device == "ideal_colorimetric"


def test_run_unknown_material(tmp_path, cli):
    result = cli.invoke(app, _run_args(tmp_path / "out", material="nope"))
    line = _assert_one_line_error(result)
    assert "ideal" in line
    assert "'nope'" in line
    assert not (tmp_path / "out").exists()


def test_run_unknown_process(tmp_path, cli):
    result = cli.invoke(app, _run_args(tmp_path / "out", process="nope"))
    line = _assert_one_line_error(result)
    assert "identity" in line


def test_run_missing_input(tmp_path, cli):
    result = cli.invoke(app, _run_args(tmp_path / "out", image=tmp_path / "missing.png"))
    line = _assert_one_line_error(result)
    assert "missing.png" in line


def test_run_unknown_readout(tmp_path, cli):
    result = cli.invoke(app, _run_args(tmp_path / "out", "-r", "nope"))
    line = _assert_one_line_error(result)
    assert "ideal_colorimetric" in line


def test_run_bad_data_path(tmp_path, cli):
    result = cli.invoke(app, _run_args(tmp_path / "out", "--data", str(tmp_path / "no-data")))
    line = _assert_one_line_error(result)
    assert "no-data" in line


def test_run_output_is_a_file(tmp_path, cli):
    taken = tmp_path / "taken"
    taken.write_text("occupied", encoding="utf-8")
    result = cli.invoke(app, _run_args(taken))
    line = _assert_one_line_error(result)
    assert "taken" in line
    assert taken.read_text(encoding="utf-8") == "occupied"


@pytest.mark.parametrize("extra", [["--iso", "-5"], ["--tile", "0"], ["--illuminant", "X"]])
def test_run_bad_values(tmp_path, cli, extra):
    result = cli.invoke(app, _run_args(tmp_path / "out", *extra))
    _assert_one_line_error(result)


@pytest.mark.parametrize("ev", ["1024", "1100", "1010", "-1100"])
def test_run_ev_out_of_range(tmp_path, cli, ev):
    result = cli.invoke(app, _run_args(tmp_path / "out", f"--ev={ev}"))
    line = _assert_one_line_error(result)
    assert "ev" in line
    assert "energy_scale_j_per_um2" not in line
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize(("name", "content"), [("empty.png", b""), ("tiny.png", b"ab")])
def test_run_unreadable_input(tmp_path, cli, name, content):
    image = tmp_path / name
    image.write_bytes(content)
    result = cli.invoke(app, _run_args(tmp_path / "out", image=image))
    line = _assert_one_line_error(result)
    assert name in line
    assert "pip" not in line


def test_check_data_ok(cli):
    result = cli.invoke(app, ["check-data"])
    assert result.exit_code == 0, result.output
    assert "OK" in result.stdout
    assert "material" in result.stdout


def test_check_data_reports(tmp_path, cli):
    films = tmp_path / "library" / "films"
    films.mkdir(parents=True)
    (films / "broken-film.yaml").write_text("schema_version: 1\nid: broken-film\n", "utf-8")
    result = cli.invoke(app, ["check-data", "--data", str(tmp_path / "library")])
    assert result.exit_code == 1, result.output
    assert "broken-film.yaml" in result.output
    assert "OK" not in result.stdout


def test_check_data_bad_path(tmp_path, cli):
    result = cli.invoke(app, ["check-data", "--data", str(tmp_path / "no-data")])
    _assert_one_line_error(result)


def _uv():
    uv = shutil.which("uv")
    if uv is None:
        pytest.skip("uv is not on PATH")
    return uv


def _subprocess_env():
    return {key: value for key, value in os.environ.items() if key != ENV_DATA_ROOT}


def test_entry_point_subprocess(tmp_path):
    command = [_uv(), "run", "atlas", "run", "tests/assets/sample.png",
               "--material", "ideal", "--process", "identity", "-o", str(tmp_path)]  # fmt: skip
    result = subprocess.run(
        command,
        cwd=REPO,
        env=_subprocess_env(),
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / OUTPUT_NAME).is_file()
    assert OUTPUT_NAME in result.stdout


def test_entry_point_error_has_no_traceback(tmp_path):
    command = [_uv(), "run", "atlas", "run", "tests/assets/sample.png",
               "--material", "nope", "--process", "identity", "-o", str(tmp_path)]  # fmt: skip
    result = subprocess.run(
        command,
        cwd=REPO,
        env=_subprocess_env(),
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    assert result.returncode == 2
    assert "Traceback" not in result.stderr
    assert "ideal" in result.stderr
