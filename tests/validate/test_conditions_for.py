"""Reference items to conditions and Kodak processes (spec D1 6.1, 7.6, 10 D1e conditions_for)."""

import dataclasses

import pytest

from darkroom_atlas.data.library import Library
from darkroom_atlas.validate import conditions as cond_mod
from darkroom_atlas.validate.conditions import (
    Condition,
    NoProcessError,
    bath_for_developer,
    conditions_for,
    dilution_for,
    fahrenheit_to_c,
    film_format_for,
    process_id_for,
    reference_conditions,
    tank_kind,
)
from darkroom_atlas.validate.splits import (
    HELD_OUT_EVALUATED,
    TRAINING,
    WITHIN_CHECKS,
    LeakageError,
    TrainingLibrary,
)

_TX, _TM = "kodak/tri-x-400", "kodak/t-max-100"


@pytest.fixture(scope="module")
def library() -> Library:
    return Library()


@pytest.fixture(scope="module")
def items(library):
    return {
        (film, split): reference_conditions(film, library, split=split)
        for film in (_TX, _TM)
        for split in ("train", "within")
    }


def _by_id(items, film, split="train"):
    return {i.item_id: i for i in items[(film, split)]}


def test_every_condition_names_an_existing_process_and_developer(library, items):
    processes = set(library.ids("process"))
    for (film, split), group in items.items():
        assert group, (film, split)
        for item in group:
            c = item.condition
            assert isinstance(c, Condition)
            assert c.material_id == film and c.source == item.item_id
            assert item.process_id in processes, item.item_id
            library.bath(c.developer)  # raises DataNotFound if there is no such developer
            step = next(s for s in library.process(item.process_id).steps if s.bath == c.developer)
            assert step.dilution == c.dilution, item.item_id
            assert step.agitation == c.agitation and step.tank == c.tank, item.item_id
            assert step.temperature.value == c.temperature_c, item.item_id
            if item.kind in ("density", "table_ci"):
                assert c.times_s == (step.time_s,), item.item_id


def test_conditions_cover_the_split(items):
    for film in (_TX, _TM):
        for split, registry in (("train", TRAINING), ("within", WITHIN_CHECKS)):
            ids = {i.item_id for i in items[(film, split)]}
            assert ids <= registry
            curves = {
                i
                for i in registry
                if i.startswith(film + "/")
                and "#" not in i
                and i.rsplit("/", 1)[1].startswith(("char-", "ci-", "mtf-", "spectral-"))
            }
            assert curves <= ids, sorted(curves - ids)


def test_conditions_for_returns_the_conditions(library, items):
    got = conditions_for(_TM, library, split="train")
    assert got == [i.condition for i in items[(_TM, "train")]]


def test_characteristic_curve_conditions(items):
    tx = _by_id(items, _TX)
    c35 = tx[f"{_TX}/char-35mm-d76-08min.csv"]
    assert c35.process_id == f"{_TX}/d-76-stock-20c-large-tank-8min"
    assert c35.condition.film_format == "135" and c35.condition.times_s == (480.0,)
    assert tx[f"{_TX}/char-120-tmax-dev-07min.csv"].condition.film_format == "120"
    tmd = tx[f"{_TX}/char-120-tmax-dev-07min.csv"]
    assert (tmd.condition.developer, tmd.condition.dilution) == ("t-max-dev", "1+4")
    tm = _by_id(items, _TM)[f"{_TM}/char-d76-07.5min.csv"]
    assert tm.process_id == f"{_TM}/d-76-stock-20c-small-tank-7.5min"
    assert tm.condition.film_format == "135"  # format not stated (spec 6.7 step 7, 7.6)


def test_ci_curve_conditions(library, items):
    hc = _by_id(items, _TX)[f"{_TX}/ci-small-tank-hc110-b.csv"]
    assert (hc.condition.developer, hc.condition.dilution) == ("hc-110-pre2019", "1+31")
    assert hc.process_id == f"{_TX}/hc-110-pre2019-1+31-20c-small-tank-3.75min"
    curve = library.reference_curve(hc.item_id)
    times = hc.condition.times_s
    assert len(times) == 5
    assert times[0] == pytest.approx(curve.x[0] * 60.0) and times[-1] == pytest.approx(
        curve.x[-1] * 60.0
    )
    large = _by_id(items, _TM)[f"{_TM}/ci-large-tank-d76.csv"]
    assert large.process_id == f"{_TM}/d-76-stock-20c-large-tank-7.25min"


def test_film_level_conditions(items):
    tm = _by_id(items, _TM)
    assert tm[f"{_TM}/mtf-d76-small-tank.csv"].condition.illuminant == "A"
    assert tm[f"{_TM}/spectral-sensitivity-d0.3.csv"].condition.exposure_time_s == 1.4
    tx = _by_id(items, _TX)
    assert tx[f"{_TX}/spectral-sensitivity-d1.0.csv"].condition.exposure_time_s == 0.5
    mtf = tx[f"{_TX}/mtf-d76-large-tank.csv"]
    assert mtf.condition.illuminant == "iso_daylight"
    assert mtf.process_id == f"{_TX}/d-76-stock-20c-large-tank-7.75min"
    gran = tm[f"{_TM}/tables.yaml#facts/rms-granularity"]
    assert gran.process_id == f"{_TM}/d-76-stock-20c-small-tank-6.5min"  # spec 6.6 step 5
    rec = tm[f"{_TM}/tables.yaml#reciprocity/10s"]
    assert rec.condition.exposure_time_s == 10.0


def test_table_cell_conditions(items):
    tx = _by_id(items, _TX)
    cell = tx[f"{_TX}/tables.yaml#development_times/normal-small-tank/d-76@65"]
    assert cell.condition.temperature_c == fahrenheit_to_c(65)  # exact (spec 6.1)
    assert cell.ci_aim == 0.56 and cell.temperature_f == 65 and cell.tank == "small-tank"
    assert cell.condition.times_s == (8.0 * 60.0,)
    j109 = tx[f"{_TX}/tables.yaml#xtol_contrast_index_times/ci-0.82/xtol@70"]
    assert j109.ci_aim == 0.82
    assert f"{_TX}/tables.yaml#xtol_contrast_index_times/ci-0.82/xtol@65" not in tx  # null
    tm = _by_id(items, _TM)
    sheets = tm[f"{_TM}/tables.yaml#development_times/normal-tray-sheets/d-76@70"]
    assert sheets.condition.film_format == "135" and sheets.tank == "tray"
    assert sheets.ci_aim == 0.56  # F-4016 prints no CI aim: 0.56 assumed (spec 7.2)
    # Within check: F-4016 [d-76] Normal column, t = f x 6.5 min (spec 7.6).
    factor = _by_id(items, _TM, "within")[
        f"{_TM}/tables.yaml#development_time_adjustment_factors/d-76/normal@65"
    ]
    assert factor.condition.times_s == pytest.approx((1.2 * 6.5 * 60.0,))
    assert factor.process_id == f"{_TM}/d-76-stock-18.3c-small-tank-7.8min"


def test_process_id_for_round_trip(library, items):
    for key in ((_TX, "train"), (_TM, "within")):
        for item in items[key][:12]:
            assert process_id_for(item.condition, library) == item.process_id


def test_missing_process_names_the_item(library, monkeypatch):
    monkeypatch.setattr(cond_mod, "_develop_steps", lambda *a: [])
    with pytest.raises(NoProcessError, match=r"kodak/tri-x-400/\S+: no process"):
        conditions_for(_TX, library, split="train")


def test_held_out_conditions_refuse_a_training_library(library):
    with pytest.raises(LeakageError):
        conditions_for(_TM, TrainingLibrary(library), split="held_out")


def test_conditions_for_bad_input(library):
    with pytest.raises(ValueError, match="split"):
        conditions_for(_TX, library, split="test")  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        conditions_for(_TX, "data", split="train")  # type: ignore[arg-type]
    # A condition without a reference item has no process to name.
    orphan = dataclasses.replace(conditions_for(_TX, library, split="train")[0], source="")
    with pytest.raises(NoProcessError, match="no reference item"):
        process_id_for(orphan, library)


@pytest.mark.parametrize(
    ("bath", "text", "expected"),
    [
        ("d-76", "stock (the figure names no dilution)", "stock"),
        ("d-76", "1+1 (printed 1:1)", "1+1"),
        ("xtol", "not stated in the figure", "stock"),
        ("hc-110-pre2019", "B (Kodak's dilution letter)", "1+31"),
        ("xtol", "full strength", "stock"),
        ("xtol", "1:1 (developer:water)", "1+1"),
        ("t-max-dev", "not stated in the figure", "1+4"),
        ("t-max-dev", "1+7 (printed 1:7)", "1+7"),
    ],
)
def test_dilution_for(bath, text, expected):
    assert dilution_for(bath, text) == expected


def test_mapping_helpers_reject_unknown_values():
    assert bath_for_developer("hc-110") == "hc-110-pre2019"  # spec 5.9
    with pytest.raises(ValueError, match="t-max-rs"):
        bath_for_developer("t-max-rs")
    with pytest.raises(ValueError, match="dilution B"):
        dilution_for("hc-110-pre2019", "A")
    with pytest.raises(ValueError):
        dilution_for("d-76", "1:3")
    assert tank_kind("small tank (8- or 16-ounce)") == "small-tank"
    assert tank_kind("tray, sheets") == "tray"
    with pytest.raises(ValueError):
        tank_kind("bucket")
    assert film_format_for("35 mm") == "135"
    assert film_format_for("120-size") == "120"
    assert film_format_for("35 mm and 120-size") == "135"
    assert film_format_for(None) == "135"


def test_held_out_conditions_cover_every_evaluated_curve_and_cell(library):
    # Structure only: no held-out time is compared or printed (spec 7.7).
    nulls: set[str] = set()
    items = reference_conditions(_TM, library, split="held_out", require_process=False, nulls=nulls)
    conditions = conditions_for(_TM, library, split="held_out")
    assert len(conditions) == len(items)
    assert all(a == b.condition for a, b in zip(conditions, items, strict=True))
    evaluated = {i for i in HELD_OUT_EVALUATED if i.startswith(_TM + "/") or "/t-max-100/" in i}
    sources = {c.source for c in conditions}
    missing = sorted(evaluated - sources - nulls)
    assert not missing, "held-out items without a condition:\n" + "\n".join(missing)
    assert sources <= HELD_OUT_EVALUATED
    tanks = cond_mod._templates(library)
    processes = set(library.ids("process"))
    for item in items:
        # Locals keep a failure message free of the Condition repr (it shows times_s).
        is_condition = isinstance(item.condition, Condition)
        material, developer = item.condition.material_id, item.condition.developer
        same_tank = (item.condition.agitation, item.condition.tank) == tanks[item.tank]
        known = item.process_id is None or item.process_id in processes
        assert is_condition and material == _TM, item.item_id
        assert developer in ("xtol", "hc-110-pre2019", "t-max-dev"), item.item_id
        library.bath(developer)
        assert same_tank and known, item.item_id
    assert conditions_for(_TX, library, split="held_out") == []


def test_process_id_for_held_out(library):
    items = reference_conditions(_TM, library, split="held_out", require_process=False)
    char = next(i for i in items if i.item_id == f"{_TM}/char-tmax-dev-06min.csv")
    assert process_id_for(char.condition, library) == f"{_TM}/t-max-dev-1+4-20c-small-tank-6min"
    cell = next(i for i in items if i.process_id is None)
    with pytest.raises(NoProcessError, match="no process file"):
        process_id_for(cell.condition, library)
