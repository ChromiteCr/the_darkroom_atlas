"""The Kodak reference processes (spec D1 2.4.8, 5.6, 5.8, 6.1, 7.6; data/processes/kodak/README.md).

The conditions are rebuilt here from the staged reference data, independently of the files:
characteristic curves from their ``# conditions:`` header, time-table cells from
``tables.yaml`` (classified with the split registry before their minutes are read, so no
held-out row is ever read).
"""

import re
import shutil
from pathlib import Path

import pytest
import yaml

from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import BathKind, Process, StepKind
from darkroom_atlas.data.provenance import Source
from darkroom_atlas.units import parse_duration_s
from darkroom_atlas.validate.conditions import parse_conditions_header, parse_temperature_c
from darkroom_atlas.validate.splits import classify, default_reference_root

FILMS = ("tri-x-400", "t-max-100")
BATH = {"d-76": "d-76", "xtol": "xtol", "hc-110": "hc-110-pre2019", "t-max-dev": "t-max-dev"}
FIXED_DILUTION = {"hc-110": "1+31", "t-max-dev": "1+4"}
# Developers of the training pairs (spec D1 1.3 decision 2); the other rows are not read.
TRAINING_DEVELOPERS = {"tri-x-400": {"d-76", "xtol", "hc-110", "t-max-dev"}, "t-max-100": {"d-76"}}
# Tank -> (agitation fields, volume mL) (spec D1 2.4.8, 6.1).
TANKS = {
    "small-tank": (
        {"method": "inversion", "initial_s": 5.0, "interval_s": 30.0, "duration_s": 5.0},
        473.0,
    ),
    "large-tank": (
        {"method": "lift", "initial_s": 22.5, "interval_s": 60.0, "duration_s": 7.5},
        3785.0,
    ),
    "rotary-tube": (
        {"method": "rotary", "initial_s": 0.0, "interval_s": None, "duration_s": None},
        250.0,
    ),
    "tray": ({"method": "tray", "initial_s": 0.0, "interval_s": None, "duration_s": None}, 1000.0),
}
KODAK_F = (65, 68, 70, 72, 75)
ID_RULE = re.compile(
    r"kodak/(?P<film>[a-z0-9-]+)/(?P<bath>[a-z0-9-]+?)-(?P<dilution>stock|\d+\+\d+)-"
    r"(?P<temp>\d+(?:\.\d)?)c-(?P<tank>small-tank|large-tank|rotary-tube|tray)-"
    r"(?P<time>\d+(?:\.\d+)?)min"
)


@pytest.fixture(scope="module")
def lib():
    return Library()


@pytest.fixture(scope="module")
def processes(lib) -> dict[str, Process]:
    return {pid: lib.process(pid) for pid in lib.ids("process") if pid.startswith("kodak/")}


@pytest.fixture(scope="module")
def staged() -> Path:
    root = default_reference_root() / "kodak"
    if not root.is_dir():
        pytest.skip(f"reference data not found at {root}")
    return root


def _temp_label(c: float) -> str:
    text = f"{c:.1f}"
    return text.removesuffix(".0")


def _time_label(minutes: float) -> str:
    return f"{round(minutes, 4):.4f}".rstrip("0").rstrip(".")


def process_id(film, bath, dilution, temperature_c, tank, minutes) -> str:
    """The naming rule of data/processes/kodak/README.md."""
    return (
        f"kodak/{film}/{bath}-{dilution}-{_temp_label(temperature_c)}c-{tank}-"
        f"{_time_label(minutes)}min"
    )


def _f_to_c(f: float) -> float:
    return (f - 32.0) * 5.0 / 9.0


def _tank_of(text: str) -> str:
    for key in ("small tank", "large tank", "rotary tube", "tray"):
        if key in text.lower():
            return key.replace(" ", "-")
    raise AssertionError(f"no tank in {text!r}")


def _dilution(developer: str, printed: str) -> str:
    if developer in FIXED_DILUTION:
        return FIXED_DILUTION[developer]
    return "1+1" if str(printed).startswith("1+1") else "stock"


def _usable(item: str) -> bool:
    c = classify(item)
    return c.split in ("training", "within_checks") or (
        c.split == "unused" and (c.reason or "").startswith("duplicates")
    )


def curve_conditions(staged: Path) -> list[tuple[str, str]]:
    """(curve item id, process id) for every characteristic curve D1 uses."""
    found = []
    for film in FILMS:
        for path in sorted((staged / film).glob("char-*.csv")):
            fields = parse_conditions_header(path.read_text())
            if fields["developer"] == "t-max-rs":
                continue
            developer = fields["developer"]
            pid = process_id(
                film,
                BATH[developer],
                _dilution(developer, fields["dilution"]),
                parse_temperature_c(fields["temperature"]),
                _tank_of(fields["agitation"]),
                float(fields["time"].removesuffix("min")),
            )
            found.append((f"kodak/{film}/{path.name}", pid))
    return found


def table_conditions(staged: Path) -> list[tuple[str, str]]:
    """(cell item id, process id) for every time-table cell the README says has a file."""
    found = []
    for film in FILMS:
        tables = yaml.safe_load((staged / film / "tables.yaml").read_text())
        for block in tables["development_times"]:
            for row in block["rows"]:
                developer = row["developer"]
                if developer not in TRAINING_DEVELOPERS[film]:
                    continue
                dilution = _dilution(developer, row["dilution"])
                suffix = "/1+1" if dilution == "1+1" else ""
                for i, temp_f in enumerate(block["temperature_f"]):
                    item = (
                        f"kodak/{film}/tables.yaml#development_times/{block['id']}/"
                        f"{developer}{suffix}@{temp_f}"
                    )
                    if not _usable(item) or row["minutes"][i] is None:
                        continue
                    pid = process_id(
                        film,
                        BATH[developer],
                        dilution,
                        _f_to_c(temp_f),
                        _tank_of(block["tank"]),
                        row["minutes"][i],
                    )
                    found.append((item, pid))
    block = yaml.safe_load((staged / "tri-x-400/tables.yaml").read_text())[
        "xtol_contrast_index_times"
    ]
    for row in block["rows"]:
        for dilution, temps, minutes in (
            ("stock", "full_strength_temperature_f", "full_strength_minutes"),
            ("1+1", "diluted_1_1_temperature_f", "diluted_1_1_minutes"),
        ):
            suffix = "/1+1" if dilution == "1+1" else ""
            for i, temp_f in enumerate(block[temps]):
                item = (
                    f"kodak/tri-x-400/tables.yaml#xtol_contrast_index_times/"
                    f"ci-{row['contrast_index']}/xtol{suffix}@{temp_f}"
                )
                try:
                    usable = _usable(item)
                except KeyError:
                    continue
                if usable and row[minutes][i] is not None:
                    pid = process_id(
                        "tri-x-400",
                        "xtol",
                        dilution,
                        _f_to_c(temp_f),
                        "small-tank",
                        row[minutes][i],
                    )
                    found.append((item, pid))
    tables = yaml.safe_load((staged / "t-max-100/tables.yaml").read_text())
    factors = tables["development_time_adjustment_factors"]
    normal = factors["columns"].index("Normal Contrast")
    (block,) = [b for b in factors["blocks"] if b["developers"] == ["d-76"]]
    for row in block["rows"]:
        item = (
            "kodak/t-max-100/tables.yaml#development_time_adjustment_factors/d-76/"
            f"normal@{row['temperature_f']}"
        )
        pid = process_id(
            "t-max-100",
            "d-76",
            "stock",
            _f_to_c(row["temperature_f"]),
            "small-tank",
            row["factor"][normal] * 6.5,  # t_std = 6 1/2 min (spec D1 7.6)
        )
        found.append((item, pid))
    return found


def _develop(process: Process):
    (step,) = [s for s in process.steps if s.kind is StepKind.DEVELOP]
    return step


def test_all_kodak_processes_load_and_check_data_is_clean(lib, processes):
    assert len(processes) > 50
    assert [p for p in lib.validate() if "processes/kodak/" in p] == []


def test_spec_example_process(processes):
    # Spec D1 5.8 example, field by field.
    process = processes["kodak/tri-x-400/d-76-stock-20c-large-tank-8min"]
    assert process.family == "bw"
    assert [s.kind for s in process.steps] == [
        StepKind.DEVELOP,
        StepKind.STOP,
        StepKind.FIX,
        StepKind.WASH,
    ]
    develop = _develop(process)
    assert (develop.bath, develop.dilution, develop.time_s) == ("d-76", "stock", 480.0)
    assert develop.temperature.to("C") == 20.0
    assert develop.temperature.source is Source.DATASHEET
    assert develop.agitation.model_dump() == TANKS["large-tank"][0]
    assert develop.tank.volume_ml == 3785.0
    assert develop.tank.film_area_cm2 == 35.0
    assert develop.tank.image_area_fraction.value == 1.0
    assert [(s.bath, s.time_s) for s in process.steps[1:]] == [
        ("stop-acetic", 30.0),
        ("fixer-thiosulfate", 300.0),
        ("water", 600.0),
    ]


def test_ids_follow_the_naming_rule(processes):
    for pid, process in processes.items():
        match = ID_RULE.fullmatch(pid)
        assert match, pid
        develop = _develop(process)
        tank = match["tank"]
        agitation, volume = TANKS[tank]
        assert develop.agitation.model_dump() == agitation, pid
        assert develop.tank.volume_ml == volume, pid
        assert pid == process_id(
            match["film"],
            develop.bath,
            develop.dilution,
            develop.temperature.to("C"),
            tank,
            develop.time_s / 60.0,
        )
        assert parse_duration_s(develop.time) == develop.time_s


def test_every_process_names_an_existing_developer_and_dilution(lib, processes):
    for pid, process in processes.items():
        develop = _develop(process)
        assert develop.bath in BATH.values(), pid
        bath = lib.bath(develop.bath)
        assert bath.kind is BathKind.DEVELOPER, pid
        assert lib.path_of("bath", develop.bath).parent.name == "developers", pid
        assert develop.dilution in bath.standard_dilutions, pid
        for step in process.steps[1:]:
            assert lib.bath(step.bath).kind in (BathKind.STOP, BathKind.FIXER, BathKind.WASH)


def test_temperatures_are_exact_fahrenheit_conversions(processes):
    # Spec D1 6.1: 65/68/70/72/75 F converted exactly, not the rounded C labels.
    for pid, process in processes.items():
        temperature = _develop(process).temperature
        assert temperature.source is Source.DATASHEET and temperature.ref, pid
        fahrenheit = temperature.to("C") * 9.0 / 5.0 + 32.0
        assert min(abs(fahrenheit - f) for f in KODAK_F) < 1e-9, pid


def test_every_characteristic_curve_condition_has_its_process(processes, staged):
    pairs = curve_conditions(staged)
    assert len(pairs) == 16 + 3 + 4  # Tri-X 35 mm + 120; T-Max D-76; T-Max T-MAX (held out)
    for item, pid in pairs:
        assert pid in processes, (item, pid)
        assert item in processes[pid].refs, (item, pid)


def test_every_training_and_within_check_table_cell_has_its_process(processes, staged):
    pairs = table_conditions(staged)
    assert len(pairs) > 100
    for item, pid in pairs:
        assert pid in processes, (item, pid)
        assert item in processes[pid].refs, (item, pid)
    # The spec D1 6.6 granularity conditions.
    assert "kodak/t-max-100/d-76-stock-20c-small-tank-6.5min" in processes
    assert "kodak/tri-x-400/hc-110-pre2019-1+31-20c-small-tank-3.75min" in processes


def test_every_training_ci_curve_has_a_process_at_its_condition(processes, staged):
    # A CI(t) curve is a time series; its condition (developer, dilution, 20 C, tank) has the
    # Kodak 20 C table time as a process.
    checked = 0
    for film in FILMS:
        for path in sorted((staged / film).glob("ci-*.csv")):
            item = f"kodak/{film}/{path.name}"
            if classify(item).split not in ("training", "within_checks"):
                continue
            fields = parse_conditions_header(path.read_text())
            developer = fields["developer"]
            prefix = (
                f"kodak/{film}/{BATH[developer]}-{_dilution(developer, fields['dilution'])}-"
                f"{_temp_label(parse_temperature_c(fields['temperature']))}c-"
                f"{_tank_of(fields['agitation'])}-"
            )
            assert any(pid.startswith(prefix) for pid in processes), (item, prefix)
            checked += 1
    # Tri-X small tank: D-76, D-76 1+1, XTOL, XTOL 1+1, HC-110 B, T-MAX; T-Max: D-76 small
    # tank stock and 1+1, D-76 large tank.
    assert checked == 6 + 3


def test_processes_cover_only_the_listed_conditions(processes, staged):
    # Every file is used by at least one listed item, and refs hold no held-out table row.
    listed = {pid for _, pid in curve_conditions(staged) + table_conditions(staged)}
    assert set(processes) == listed
    for pid, process in processes.items():
        assert process.refs, pid
        for item in process.refs:
            split = classify(item).split
            if split == "held_out":
                assert re.fullmatch(r"kodak/t-max-100/char-tmax-dev-\d+min\.csv", item), item


def test_unknown_bath_is_reported(tmp_path):
    root = tmp_path / "data"
    shutil.copytree(Library().root, root)
    path = root / "processes/kodak/tri-x-400/d-76-stock-20c-large-tank-8min.yaml"
    path.write_text(path.read_text().replace("bath: d-76", "bath: d-77"))
    problems = Library(root).validate()
    assert any("steps[0].bath: unknown bath 'd-77'" in p for p in problems)
