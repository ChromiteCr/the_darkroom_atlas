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
    HELD_OUT_EVALUATED,
    TRAINING,
    UNUSED,
    WITHIN_CHECKS,
    LeakageError,
    TrainingLibrary,
)

_TX, _TM = "kodak/tri-x-400", "kodak/t-max-100"

# Training items for which spec 7.2 defines no observable of their own (reasons for the
# coverage test); every other training item must give at least one observation.
_NO_OBSERVABLE = {
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


def test_training_observations_refuse_a_held_out_item(library, monkeypatch):
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


def test_held_out_observations_on_synthetic_reference(library, tmp_path):
    # Spec 8.3: the held-out pipeline is exercised on synthetic data, not on Kodak's files.
    ref = tmp_path / "reference" / _TM
    ref.mkdir(parents=True)
    (ref / "ci-small-tank-xtol.csv").write_text(_SYNTHETIC_CI, encoding="utf-8")
    (ref / "tables.yaml").write_text(_SYNTHETIC_TABLES, encoding="utf-8")
    proc = tmp_path / "processes" / _TM
    proc.mkdir(parents=True)
    src = library.path_of("process", f"{_TM}/d-76-stock-20c-small-tank-6.5min")
    shutil.copy(src, proc / src.name)
    synthetic = Library(tmp_path)
    obs = held_out_observations(synthetic, CalibrationRecord.model_construct(frozen=True))
    kinds = Counter(o.kind for o in obs)
    assert kinds == {"ci": 1, "table_ci": 4}  # the null 72 F cell is skipped
    for o in obs:
        assert set(observation_sources(o)) <= HELD_OUT_EVALUATED
        assert not o.in_objective and np.all(np.isnan(o.sigma))
        assert o.condition is not None and o.condition.developer == "xtol"
    ci = _one(obs, f"{_TM}/ci-small-tank-xtol.csv")
    assert ci.x[0] == 300.0 and ci.x[-1] == 660.0 and ci.y[0] == pytest.approx(0.40)
    cell = _one(obs, f"{_TM}/tables.yaml#development_times/normal-small-tank-rolls/xtol@65")
    assert cell.x[0] == 540.0 and cell.y[0] == 0.56
    assert cell.condition.temperature_c == (65 - 32) * 5 / 9


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
    cc = _one(training, "bath:d-76-published#ph_anchors/stock")
    assert cc.sigma[0] == 0.10  # Carlton & Crabtree colorimetric (spec 7.2)
    assert _one(training, "bath:hc-110-pre2019#ph_anchors/1+15").sigma[0] == 0.10


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
    assert np.array_equal(tx.sigma, np.where(tx.x >= 50, 3.0, 2.0))
    curve = library.reference_curve(tx.item_id)
    last_full = curve.x[np.flatnonzero(curve.y >= 100)[-1]]
    assert np.all(tx.x > last_full) and tx.in_objective
    tm = _one(training, f"{_TM}/mtf-d76-small-tank.csv")
    assert np.all(tm.sigma == 2.0) and 150.0 in tm.x
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
