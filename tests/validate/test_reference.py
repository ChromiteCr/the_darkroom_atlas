"""Reference observations, sigma and weights, and the leakage rules (spec D1 7.2, 7.6, 7.7)."""

import math
import shutil
from collections import Counter

import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator

from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import CalibrationRecord
from darkroom_atlas.validate import reference as ref_mod
from darkroom_atlas.validate.conditions import ConditionedItem, reference_conditions
from darkroom_atlas.validate.metrics import time_to_ci
from darkroom_atlas.validate.reference import (
    Observation,
    held_out_observations,
    observation_counts,
    observation_sources,
    training_observations,
    within_check_observations,
)
from darkroom_atlas.validate.splits import (
    HELD_OUT,
    TRAINING,
    UNUSED,
    WITHIN_CHECKS,
    LeakageError,
    TrainingLibrary,
    classify,
    enumerate_items,
)

_TX, _TM = "kodak/tri-x-400", "kodak/t-max-100"

# Training items for which spec 7.2 defines no observable of their own (reasons for the
# coverage test); every other training item must give at least one observation.
_NO_OBSERVABLE = {
    "developer_definitions/t-max-developers-contrast-aim": "J-86 CI aim of T-MAX cells (0.56)",
    "facts/base-plus-fog/": "left end of a curve: compared in the 8.3 report, no 7.2 residual",
    "facts/processing-aim": "CI aim of the table cells (the y of table_ci), not a residual",
    "developer_definitions/": "dilution definitions: they fix conditions, not observations",
}


@pytest.fixture(scope="module")
def library() -> Library:
    return Library()


@pytest.fixture(scope="module")
def training(library) -> list[Observation]:
    return training_observations(TrainingLibrary(library))


@pytest.fixture(scope="module")
def within(library) -> list[Observation]:
    return within_check_observations(library)


def _one(observations, item_id):
    hits = [o for o in observations if o.item_id == item_id]
    assert len(hits) == 1, item_id
    return hits[0]


def _null_cell(library, item_id) -> bool:
    """True for a registered table cell whose reference value is null (skipped by 7.2)."""
    film, _, rest = item_id.partition("/tables.yaml#")
    if not rest:
        return False
    return not any(i.item_id == item_id for i in reference_conditions(film, library, split="train"))


# ---- leakage (spec 7.7) ----------------------------------------------------------------------


def test_training_observations_contain_only_training_ids(training):
    sources = {s for o in training for s in observation_sources(o)}
    assert sources <= TRAINING
    assert not sources & (HELD_OUT | UNUSED | WITHIN_CHECKS)
    for o in training:
        if not observation_sources(o):
            assert o.kind == "ph" and o.item_id.startswith(ref_mod.PH_PREFIX)


def test_training_observations_require_training_library(library):
    with pytest.raises(TypeError, match="TrainingLibrary"):
        training_observations(library)  # type: ignore[arg-type]


def test_calibration_loader_refuses_held_out(library, monkeypatch):
    # Spec 7.7: the calibration data loader (training_observations) raises on a held-out id.
    real = ref_mod.reference_conditions
    held = f"{_TM}/tables.yaml#development_times/normal-small-tank-rolls/xtol@68"
    assert held in HELD_OUT

    def leaky(material_id, lib, *, split, require_process=True):
        items = real(material_id, lib, split=split, require_process=require_process)
        if material_id == _TM:
            gran = next(i for i in items if i.kind == "granularity")
            items.append(ConditionedItem(held, "granularity", gran.condition, gran.process_id))
        return items

    monkeypatch.setattr(ref_mod, "reference_conditions", leaky)
    with pytest.raises(LeakageError, match="held out"):
        training_observations(TrainingLibrary(library))


def test_held_out_requires_frozen_record(library):
    with pytest.raises(LeakageError, match="frozen"):
        held_out_observations(library, CalibrationRecord.model_construct(frozen=False))
    with pytest.raises(LeakageError, match="frozen"):
        held_out_observations(library, object())  # type: ignore[arg-type]
    with pytest.raises(LeakageError, match="TrainingLibrary"):
        held_out_observations(
            TrainingLibrary(library), CalibrationRecord.model_construct(frozen=True)
        )


_SYNTHETIC_CI = """# material: kodak/t-max-100
# quantity: contrast_index
# x_unit: min
# y_unit: 1
# origin: datasheet
# ref: synthetic test curve, not Kodak data
# conditions: developer=xtol; dilution=stock; temperature=20 C; agitation=small tank; format=not stated; densitometry=diffuse visual
5.0,0.40
7.0,0.55
9.0,0.66
11.0,0.74
"""
_SYNTHETIC_TABLES = """schema_version: 1
material: kodak/t-max-100
development_times:
- id: normal-small-tank-rolls
  tank: small tank (8- or 16-ounce), rolls
  temperature_c: [18, 20, 21, 22, 24]
  temperature_f: [65, 68, 70, 72, 75]
  rows:
  - developer: xtol
    dilution: stock
    minutes: [9.0, 8.0, 7.0, null, 5.0]
"""


def _synthetic_char(minutes: int, gamma: float) -> str:
    """A synthetic T-MAX Developer characteristic curve (not Kodak data)."""
    log_h = np.linspace(-3.0, 0.0, 31)
    density = 0.20 + gamma * np.log1p(np.exp(4.0 * (log_h + 2.2))) / 4.0
    rows = "".join(f"{x:.4f},{y:.4f}\n" for x, y in zip(log_h, density, strict=True))
    return (
        "# material: kodak/t-max-100\n# quantity: characteristic_curve\n"
        "# x_unit: log10(lux*s)\n# y_unit: 1\n# origin: datasheet\n"
        "# ref: synthetic test curve, not Kodak data\n"
        f"# conditions: developer=t-max-dev; dilution=not stated; time={minutes}min; "
        "temperature=20 C; agitation=small tank; format=not stated; densitometry=diffuse visual\n"
        + rows
    )


# Table ids the split registry classifies as evaluated held-out (7.6).
_TABLE_ID = {
    "xtol": "table1-small-tank-full-strength",
    "hc-110": "roll-dilution-b-small-tank",
    "t-max-dev": "small-tank-rolls",
}


def _synthetic_developer(name: str, ci: float | None, minutes: list) -> str:
    """A developer-publication file with one T-Max row (and, for T-MAX Dev, a push row)."""
    ci_field = f"    ci: {ci}\n" if ci is not None else ""
    ei_field = "" if name == "hc-110" else "    ei: '100/200'\n"  # J-24 rows print no EI
    push = (
        "  - {film_as_printed: T-MAX 100, material: kodak/t-max-100, format: roll, ei: '400',"
        " minutes: [6.0, 5.0, null]}\n"
        if name == "t-max-dev"
        else ""
    )
    return (
        f"schema_version: 1\ndeveloper: {name}\ndevelopment_times:\n"
        f"- id: {_TABLE_ID[name]}\n  tank: small tank\n"
        f"  dilution: {'B' if name == 'hc-110' else 'full strength'}\n"
        "  temperature_c: [20, 21, 24]\n  temperature_f: [68, 70, 75]\n  rows:\n"
        "  - film_as_printed: T-MAX 100\n    material: kodak/t-max-100\n"
        f"    format: {'135/120' if name == 'xtol' else 'roll'}\n"
        + ei_field
        + ci_field
        + f"    minutes: {minutes}\n".replace("None", "null")
        + push
    )


def _synthetic_held_out(library, tmp_path) -> Library:
    # Spec 8.3: the held-out pipeline is exercised on synthetic data, not on Kodak's files.
    ref = tmp_path / "reference" / _TM
    ref.mkdir(parents=True)
    (ref / "ci-small-tank-xtol.csv").write_text(_SYNTHETIC_CI, encoding="utf-8")
    (ref / "tables.yaml").write_text(_SYNTHETIC_TABLES, encoding="utf-8")
    for minutes, gamma in ((6, 0.45), (7, 0.55), (10, 0.75), (12, 0.90)):
        text = _synthetic_char(minutes, gamma)
        (ref / f"char-tmax-dev-{minutes:02d}min.csv").write_text(text, encoding="utf-8")
    dev = tmp_path / "reference" / "developers"
    dev.mkdir()
    for name, ci, minutes in (
        ("xtol", 0.62, [8.0, 7.0, 5.0]),
        ("hc-110", None, [6.0, None, 4.0]),
        ("t-max-dev", None, [None, None, None]),
    ):
        text = _synthetic_developer(name, ci, minutes)
        (dev / f"{name}.yaml").write_text(text, encoding="utf-8")
    proc = tmp_path / "processes" / _TM
    proc.mkdir(parents=True)
    src = library.path_of("process", f"{_TM}/d-76-stock-20c-small-tank-6.5min")
    shutil.copy(src, proc / src.name)
    return Library(tmp_path)


_DT = "#development_times"
_X_ROW = f"developers/xtol.yaml{_DT}/{_TABLE_ID['xtol']}/t-max-100/135-120/ei100-200/ci0.62"
_HC_ROW = f"developers/hc-110.yaml{_DT}/{_TABLE_ID['hc-110']}/t-max-100/roll"
_TMD_ROW = f"developers/t-max-dev.yaml{_DT}/{_TABLE_ID['t-max-dev']}/t-max-100/roll/ei100-200"


def test_held_out_observations_on_synthetic_reference(library, tmp_path):
    synthetic = _synthetic_held_out(library, tmp_path)
    obs = held_out_observations(synthetic, CalibrationRecord.model_construct(frozen=True))
    kinds = Counter(o.kind for o in obs)
    # 4 tables.yaml cells (the null 72 F cell skipped), 3 J-109 and 2 J-24 cells; the
    # all-null T-MAX Dev row and its push row (not evaluated) give nothing.
    assert kinds == {"ci": 1, "table_ci": 9, "density": 4, "speed": 1}
    for o in obs:
        for source in observation_sources(o):
            assert classify(source).split == "held_out" and classify(source).evaluated
        assert not o.in_objective and np.all(np.isnan(o.sigma))
        assert o.condition is not None and o.condition.material_id == _TM
    ci = _one(obs, f"{_TM}/ci-small-tank-xtol.csv")
    assert ci.x[0] == 300.0 and ci.x[-1] == 660.0 and ci.y[0] == pytest.approx(0.40)
    cell = _one(obs, f"{_TM}/tables.yaml#development_times/normal-small-tank-rolls/xtol@65")
    assert cell.x[0] == 540.0 and cell.y[0] == 0.56
    assert cell.condition.temperature_c == (65 - 32) * 5 / 9
    j109 = [o for o in obs if o.item_id == _X_ROW]
    assert [o.x[0] for o in j109] == [480.0, 420.0, 300.0]
    assert all(o.y[0] == 0.62 for o in j109)  # J-109: the row's own CI (spec 7.2)
    assert j109[1].condition.temperature_c == (70 - 32) * 5 / 9
    assert j109[0].condition.dilution == "stock"
    j24 = [o for o in obs if o.item_id == _HC_ROW]
    assert [o.y[0] for o in j24] == [0.56, 0.56]  # J-24: CI_aim 0.56 (spec 7.2)
    assert j24[0].condition.developer == "hc-110-pre2019"
    assert j24[0].condition.dilution == "1+31"
    speed = _one(
        [o for o in obs if o.kind == "speed"], next(o.item_id for o in obs if o.kind == "speed")
    )
    assert len(observation_sources(speed)) == 4 and math.isfinite(speed.y[0])
    assert not speed.in_objective  # 8.3: log H_m* of the held-out T-MAX Dev series


def test_every_evaluated_held_out_item_yields_an_observation(library, tmp_path):
    synthetic = _synthetic_held_out(library, tmp_path)
    obs = held_out_observations(synthetic, CalibrationRecord.model_construct(frozen=True))
    covered = {s for o in obs for s in observation_sources(o)}
    evaluated = {
        i
        for i in enumerate_items(synthetic.root / "reference")
        if classify(i).split == "held_out" and classify(i).evaluated
    }
    # Reasoned exclusions: cells or rows whose reference value is null (spec 7.2).
    null_ids = {
        f"{_TM}/tables.yaml#development_times/normal-small-tank-rolls/xtol@72",
        _TMD_ROW,
    }
    assert null_ids <= evaluated
    assert evaluated - null_ids == covered, sorted(evaluated ^ covered)


def test_held_out_observations_refuse_to_skip_an_item(library, tmp_path, monkeypatch):
    synthetic = _synthetic_held_out(library, tmp_path)
    real = ref_mod._film_observations

    def lossy(material_id, lib, items, *, mode):
        return [o for o in real(material_id, lib, items, mode=mode) if o.item_id != _HC_ROW]

    monkeypatch.setattr(ref_mod, "_film_observations", lossy)
    with pytest.raises(ValueError, match="hc-110.yaml"):
        held_out_observations(synthetic, CalibrationRecord.model_construct(frozen=True))


# ---- coverage ----------------------------------------------------------------------------------


def test_every_training_item_with_an_observable_yields_an_observation(library, training):
    covered = {s for o in training for s in observation_sources(o)}
    missing = []
    for item_id in sorted(TRAINING - covered):
        reason = next((r for k, r in _NO_OBSERVABLE.items() if k in item_id), None)
        if reason is None and _null_cell(library, item_id):
            reason = "reference value is null"
        if reason is None:
            missing.append(item_id)
    assert not missing, "training items without an observation:\n" + "\n".join(missing)


def test_observation_counts(training, within, capsys):
    counts = observation_counts(training)
    with capsys.disabled():
        print("\ntraining observations per (film, kind):")
        for key, n in counts.items():
            print(f"  {key[0]:18s} {key[1]:12s} {n}")
        print("within-check observations per (film, kind):")
        for key, n in observation_counts(within).items():
            print(f"  {key[0]:18s} {key[1]:12s} {n}")
    assert counts == {
        ("bath", "ph"): 5,
        (_TM, "ci"): 3,
        (_TM, "density"): 3,
        (_TM, "granularity"): 1,
        (_TM, "mtf"): 1,
        (_TM, "reciprocity"): 6,
        (_TM, "spectral"): 2,
        (_TM, "speed"): 1,
        (_TM, "table_ci"): 28,
        (_TX, "ci"): 6,
        (_TX, "density"): 16,
        (_TX, "mtf"): 1,
        (_TX, "reciprocity"): 6,
        (_TX, "spectral"): 2,
        (_TX, "speed"): 4,
        (_TX, "table_ci"): 38,
    }


# ---- sigma and weights (spec 7.2) --------------------------------------------------------------


def test_density_sigma_and_weights(library, training):
    o = _one(training, f"{_TX}/char-35mm-d76-08min.csv")
    curve = library.reference_curve(o.item_id)
    d_min = float(curve.y[0])
    assert o.d_min_ref == d_min
    assert o.x.size == 15 and o.n_effective == 6.0
    assert np.all(o.sigma == 0.025 * math.sqrt(15 / 6))
    assert o.sigma[0] == pytest.approx(0.040, abs=5e-4)  # spec 7.2 writes the rounded 0.040
    assert o.x[-1] == curve.x[-1]
    assert o.y[0] == pytest.approx(d_min + 0.05, abs=1e-9)
    assert np.allclose(np.diff(o.x), np.diff(o.x)[0])
    assert np.array_equal(o.weight, np.where(o.y < d_min + 0.1, 0.5, 1.0))
    assert 0.5 in o.weight and 1.0 in o.weight
    assert o.in_objective


def test_ci_sigma(library, training):
    tx = _one(training, f"{_TX}/ci-small-tank-d76.csv")
    assert np.all(tx.sigma == 0.025 * math.sqrt(5 / 3))  # sigma_floor 0 for Tri-X
    assert tx.x.size == 5 and tx.n_effective == 3.0
    curve = library.reference_curve(tx.item_id)
    assert tx.x[0] == pytest.approx(curve.x[0] * 60) and tx.x[-1] == pytest.approx(curve.x[-1] * 60)
    floor = ref_mod.sigma_floor(_TM, TrainingLibrary(library))
    assert floor == pytest.approx(0.036 / 2, abs=1e-3)  # spec 7.2: 0.694 vs 0.658 at 7.5 min
    tm = _one(training, f"{_TM}/ci-small-tank-d76.csv")
    assert np.all(tm.sigma == math.sqrt(0.025**2 + floor**2) * math.sqrt(5 / 3))


def _slope(library, rel, aim):
    curve = library.reference_curve(rel)
    p = PchipInterpolator(curve.x, curve.y)
    t = time_to_ci(curve.x, curve.y, aim)
    return float(p.derivative()(t))


def test_table_ci_sigma(library, training):
    rounding = 0.25 / math.sqrt(12)
    s = _slope(library, f"{_TX}/ci-small-tank-d76.csv", 0.56)
    cell = _one(training, f"{_TX}/tables.yaml#development_times/normal-small-tank/d-76@65")
    assert cell.y[0] == 0.56 and cell.x[0] == 480.0 and cell.weight[0] == 1.0
    n_class = 4  # 65/70/72/75 F of D-76 stock small tank at CI 0.56
    assert cell.sigma[0] == pytest.approx(math.sqrt(0.02**2 + (s * rounding) ** 2) * 2, rel=1e-9)
    assert cell.n_effective == pytest.approx(1 / n_class)
    at_20 = _one(training, f"{_TX}/tables.yaml#development_times/normal-large-tank/d-76@68")
    assert at_20.sigma[0] == pytest.approx(math.sqrt(0.02**2 + (s * rounding) ** 2), rel=1e-9)
    sx = _slope(library, f"{_TX}/ci-small-tank-xtol.csv", 0.56)
    xtol = _one(training, f"{_TX}/tables.yaml#development_times/normal-small-tank/xtol@65")
    expected = math.sqrt(0.02**2 + (sx * rounding) ** 2) * 2 * 2  # sqrt(4), XTOL stock x 2
    assert xtol.sigma[0] == pytest.approx(expected, rel=1e-9)
    s82 = _slope(library, f"{_TX}/ci-small-tank-xtol.csv", 0.82)
    j109 = _one(training, f"{_TX}/tables.yaml#xtol_contrast_index_times/ci-0.82/xtol@70")
    assert j109.y[0] == 0.82
    expected = math.sqrt(0.02**2 + (s82 * rounding) ** 2) * math.sqrt(2) * 2  # 65 F is "—"
    assert j109.sigma[0] == pytest.approx(expected, rel=1e-9)


def test_ph_sigma(training):
    xtol = _one(training, "bath:xtol#ph_anchors/stock")
    assert xtol.y[0] == 8.20 and xtol.sigma[0] == 0.05 and xtol.condition is None
    cc = _one(training, "bath:d-76-published#ph_anchors/stock#delta_cc")
    assert cc.sigma[0] == 0.10  # Carlton & Crabtree colorimetric (spec 7.2)
    assert ref_mod.takes_delta_cc(cc) and not ref_mod.takes_delta_cc(xtol)
    assert [o.item_id for o in training if ref_mod.takes_delta_cc(o)] == [cc.item_id]
    hc = _one(training, "bath:hc-110-pre2019#ph_anchors/1+15")
    assert hc.sigma[0] == 0.10 and not ref_mod.takes_delta_cc(hc)  # CIS-61, its own sigma


def test_spectral_sigma_and_shape(library, training):
    o = _one(training, f"{_TX}/spectral-sensitivity-d1.0.csv")
    curve = library.reference_curve(o.item_id)
    assert np.all(o.sigma == 0.03) and np.all(o.weight == 1.0)
    assert np.all(o.x % 20 == 0) and o.x[0] >= 400 and o.x[-1] <= min(700, curve.x[-1])
    assert o.x[-1] < 660  # Tri-X D = 1.0 ends at 655.3 nm (spec 6.9 step 4)
    window = (o.x >= 450) & (o.x <= 600)
    assert np.all(o.y[window] <= 1e-12)
    tm = _one(training, f"{_TM}/spectral-sensitivity-d0.3.csv")
    span = library.reference_curve(tm.item_id).x
    expected = [b for b in range(400, 701, 20) if span[0] <= b <= span[-1]]
    assert tm.x.tolist() == expected  # 400-700 nm every 20 nm inside the span (6.9 step 4)


def test_mtf_tail_sigma(library, training, within):
    tx = _one(training, f"{_TX}/mtf-d76-large-tank.csv")
    # Pre-registered tail bands (6.7 step 8): any change to them must fail here.
    assert tx.x.tolist() == [30.0, 40.0, 50.0, 60.0, 70.0]
    assert tx.sigma.tolist() == [2.0, 2.0, 3.0, 3.0, 3.0]
    assert np.array_equal(tx.sigma, np.where(tx.x >= 50, 3.0, 2.0))
    curve = library.reference_curve(tx.item_id)
    last_full = curve.x[np.flatnonzero(curve.y >= 100)[-1]]
    assert np.all(tx.x > last_full) and tx.in_objective
    tm = _one(training, f"{_TM}/mtf-d76-small-tank.csv")
    assert tm.x.tolist() == [50.0, 60.0, 70.0, 80.0, 100.0, 120.0, 150.0]
    assert np.all(tm.sigma == 2.0)
    peak = _one(within, f"{_TX}/mtf-d76-large-tank.csv")
    assert not peak.in_objective and np.all(peak.x < tx.x[0])


def test_granularity_and_reciprocity_sigma(training):
    g = _one(training, f"{_TM}/tables.yaml#facts/rms-granularity")
    assert (g.x[0], g.y[0], g.sigma[0]) == (1.0, 8.0, 0.6)
    r = _one(training, f"{_TX}/tables.yaml#reciprocity/10s")
    assert (r.x[0], r.y[0], r.sigma[0]) == (10.0, 2.0, 0.15)
    assert r.condition.exposure_time_s == 10.0
    assert _one(training, f"{_TM}/tables.yaml#reciprocity/0.0001s").y[0] == pytest.approx(0.3333)


def test_speed_is_reported_not_fitted(training):
    speeds = {o.item_id.split(" + ")[0]: o for o in training if o.kind == "speed"}
    s35 = speeds[f"{_TX}/char-35mm-d76-06min.csv"]
    s120 = speeds[f"{_TX}/char-120-d76-06min.csv"]
    assert not s35.in_objective and s35.weight[0] == 0.0
    assert len(observation_sources(s35)) == 4
    # Spec 6.3 test values (DS 5.3): log H_m about -2.73 (35 mm) and -2.69 (120), +- 0.03.
    assert s35.y[0] == pytest.approx(-2.73, abs=0.03)
    assert s120.y[0] == pytest.approx(-2.69, abs=0.03)


def test_within_checks_are_never_in_the_objective(within):
    assert within
    for o in within:
        assert not o.in_objective
        for s in observation_sources(o):
            assert s in WITHIN_CHECKS or s in TRAINING
            assert s not in HELD_OUT
    ids = {o.item_id for o in within}
    assert f"{_TM}/tables.yaml#development_time_adjustment_factors/d-76/normal@65" in ids
    assert f"{_TX}/tables.yaml#reciprocity/100s" in ids


# ---- Observation itself ------------------------------------------------------------------------


def test_observation_validation():
    ok = {"kind": "ph", "item_id": "bath:x#ph_anchors/stock", "condition": None, "x": [25.0],
          "y": [8.2], "sigma": [0.05], "weight": [1.0], "n_effective": 1.0}  # fmt: skip
    o = Observation(**ok)
    assert o.x.dtype == np.float64 and not o.x.flags.writeable
    with pytest.raises(ValueError, match="kind"):
        Observation(**ok | {"kind": "fog"})
    with pytest.raises(ValueError, match="1-D"):
        Observation(**ok | {"y": [8.2, 8.3]})
    with pytest.raises(ValueError, match="sigma"):
        Observation(**ok | {"sigma": [math.nan]})
    Observation(**ok | {"sigma": [math.nan], "in_objective": False})
    with pytest.raises(ValueError, match="weight"):
        Observation(**ok | {"weight": [-1.0]})
    with pytest.raises(ValueError, match="finite"):
        Observation(**ok | {"y": [math.inf]})
    with pytest.raises(ValueError, match="n_effective"):
        Observation(**ok | {"n_effective": 0.0})
    with pytest.raises(ValueError, match="item_id"):
        Observation(**ok | {"item_id": ""})


def test_table_ci_sigma_rejects_bad_class():
    with pytest.raises(ValueError, match="n_class"):
        ref_mod.table_ci_sigma(0.1, 0, xtol_stock_off_20c=False)
