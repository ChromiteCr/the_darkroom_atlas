"""Train / within-check / held-out / unused split registry (spec D1 7.6, 7.7, 8.3, 8.4)."""

import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import replace
from pathlib import Path

import pytest
import yaml

from darkroom_atlas.data.library import Library
from darkroom_atlas.validate import splits
from darkroom_atlas.validate.splits import (
    LeakageError,
    TrainingLibrary,
    assert_trainable,
    build_registry,
    classify,
    enumerate_items,
)

_STAGED = Path(__file__).resolve().parents[2] / "data" / "reference"


@pytest.fixture(scope="module")
def root() -> Path:
    if not _STAGED.is_dir():
        pytest.skip(f"reference data not found at {_STAGED}")
    return _STAGED


@pytest.fixture(scope="module")
def registry(root: Path) -> splits.Registry:
    return build_registry(root)


@pytest.fixture(scope="module")
def items(root: Path) -> list[str]:
    return enumerate_items(root)


def _counts(registry: splits.Registry) -> str:
    return (
        f"training={len(registry.training)} within_checks={len(registry.within_checks)} "
        f"held_out={len(registry.held_out)} (evaluated {len(registry.held_out_evaluated)}) "
        f"unused={len(registry.unused)}"
    )


# --------------------------------------------------------------------------------------
# Completeness (spec D1 10, D1e)
# --------------------------------------------------------------------------------------


def test_split_registry_complete(root: Path, registry: splits.Registry, items: list[str]) -> None:
    counts = _counts(registry)
    classes = [registry.training, registry.within_checks, registry.held_out, registry.unused]
    membership = Counter(i for c in classes for i in c)
    assert set(membership) == set(items), counts
    multiple = sorted(i for i, n in membership.items() if n != 1)
    assert not multiple, f"{counts}; in more than one class: {multiple[:10]}"
    assert len(items) == len(set(items)), counts

    assert registry.held_out_evaluated <= registry.held_out, counts
    needs_reason = (registry.held_out - registry.held_out_evaluated) | registry.unused
    assert set(registry.not_evaluated_reason) == needs_reason, counts
    assert all(r.strip() for r in registry.not_evaluated_reason.values()), counts

    # Every file under the root is classified, as a whole or through its rows.
    for path in root.rglob("*"):
        if path.is_file() and not any(p.startswith(".") for p in path.relative_to(root).parts):
            rel = path.relative_to(root).as_posix()
            assert any(i == rel or i.startswith(rel + "#") for i in items), f"{counts}; {rel}"

    # J-78: every row with material: null is unused (spec 7.6).
    j78 = yaml.safe_load((root / "developers" / "d-76.yaml").read_text(encoding="utf-8"))
    n_null = sum(
        1 for t in j78["development_times"] for r in t["rows"] if r.get("material") is None
    )
    j78_rows = [
        i
        for i in items
        if i.startswith("developers/d-76.yaml#development_times/") and "/~" not in i
    ]
    assert n_null == len(j78_rows) > 0, counts
    assert set(j78_rows) <= registry.unused, counts

    # The 14 pairs of the pair-metric table (spec 8.3).
    assert len(splits.PAIR_METRICS) == 14, counts
    assert set(splits.PAIR_METRICS) == set(splits.PAIR_TIERS)
    print(counts)


def test_tables_rows_and_cells_complete(root: Path, items: list[str]) -> None:
    """Every printed time cell of every film table has its own id (spec 7.6)."""
    for film in ("tri-x-400", "t-max-100"):
        rel = f"kodak/{film}/tables.yaml"
        data = yaml.safe_load((root / rel).read_text(encoding="utf-8"))
        n_cells = sum(len(t["temperature_f"]) * len(t["rows"]) for t in data["development_times"])
        cells = [i for i in items if "/~" not in i]  # free-text sub-items are counted below
        ids = [i for i in cells if i.startswith(f"{rel}#development_times/")]
        assert len(ids) == n_cells
        assert sum(i.startswith(f"{rel}#reciprocity/") and i.count("/") == 3 for i in cells) == len(
            data["reciprocity"]["rows"]
        )
        xt = data["xtol_contrast_index_times"]
        n_xt = len(xt["rows"]) * (
            len(xt["full_strength_temperature_f"]) + len(xt["diluted_1_1_temperature_f"])
        )
        assert sum(i.startswith(f"{rel}#xtol_contrast_index_times/") for i in cells) == n_xt


def test_module_constants_match_default_root() -> None:
    reg = build_registry()
    assert splits.TRAINING == reg.training
    assert splits.WITHIN_CHECKS == reg.within_checks
    assert splits.HELD_OUT == reg.held_out
    assert splits.UNUSED == reg.unused
    assert splits.HELD_OUT_EVALUATED == reg.held_out_evaluated
    assert dict(splits.NOT_EVALUATED_REASON) == dict(reg.not_evaluated_reason)
    with pytest.raises(AttributeError):
        _ = splits.NO_SUCH_CONSTANT  # type: ignore[attr-defined]


# --------------------------------------------------------------------------------------
# Items the spec names (spec D1 7.6), transcribed by hand
# --------------------------------------------------------------------------------------

_X = "kodak/tri-x-400/"
_M = "kodak/t-max-100/"
_XT = "kodak/tri-x-400/tables.yaml#"
_MT = "kodak/t-max-100/tables.yaml#"

_TRAINING = [
    # Film block, Tri-X.
    *(f"{_X}char-35mm-d76-{m}min.csv" for m in ("06", "08", "10", "12")),
    *(f"{_X}char-120-d76-{m}min.csv" for m in ("06", "08", "10", "12")),
    f"{_X}mtf-d76-large-tank.csv",
    f"{_X}spectral-sensitivity-d0.3.csv",
    f"{_X}spectral-sensitivity-d1.0.csv",
    *(f"{_XT}reciprocity/{t}" for t in ("0.0001s", "0.001s", "0.01s", "0.1s", "1s", "10s")),
    # Film block, T-Max.
    *(f"{_M}char-d76-{m}min.csv" for m in ("06", "07.5", "10")),
    f"{_M}mtf-d76-small-tank.csv",
    f"{_M}spectral-sensitivity-d0.3.csv",
    f"{_M}spectral-sensitivity-d1.0.csv",
    f"{_MT}facts/rms-granularity",
    *(f"{_MT}reciprocity/{t}" for t in ("0.0001s", "0.001s", "0.01s", "0.1s", "1s", "10s")),
    # Tri-X x D-76.
    f"{_X}ci-small-tank-d76.csv",
    f"{_X}ci-small-tank-d76-1-1.csv",
    *(f"{_XT}development_times/normal-small-tank/d-76@{t}" for t in (65, 70, 72, 75)),
    *(f"{_XT}development_times/normal-large-tank/d-76@{t}" for t in (65, 68, 70, 72, 75)),
    f"{_XT}development_times/normal-large-tank/d-76/1+1@68",
    # Tri-X x XTOL.
    f"{_X}ci-small-tank-xtol.csv",
    f"{_X}ci-small-tank-xtol-1-1.csv",
    *(f"{_XT}development_times/normal-small-tank/xtol@{t}" for t in (65, 70, 72, 75)),
    *(f"{_XT}development_times/normal-large-tank/xtol@{t}" for t in (65, 68, 70, 72, 75)),
    f"{_XT}development_times/normal-large-tank/xtol/1+1@68",
    *(f"{_XT}xtol_contrast_index_times/ci-0.72/xtol@{t}" for t in (65, 70, 75)),
    *(f"{_XT}xtol_contrast_index_times/ci-0.82/xtol@{t}" for t in (65, 70, 75)),
    # Tri-X x HC-110.
    f"{_X}ci-small-tank-hc110-b.csv",
    *(f"{_XT}development_times/normal-small-tank/hc-110@{t}" for t in (65, 70, 72, 75)),
    *(f"{_XT}development_times/normal-large-tank/hc-110@{t}" for t in (65, 68, 70, 72, 75)),
    # Tri-X x T-MAX Developer.
    *(f"{_X}char-35mm-tmax-dev-{m}min.csv" for m in ("06", "07", "09", "11")),
    *(f"{_X}char-120-tmax-dev-{m}min.csv" for m in ("06", "07", "09", "11")),
    f"{_X}ci-small-tank-tmax-dev.csv",
    *(f"{_XT}development_times/normal-small-tank/t-max-dev@{t}" for t in (65, 70, 72, 75)),
    # T-Max x D-76.
    f"{_M}ci-small-tank-d76.csv",
    f"{_M}ci-small-tank-d76-1-1.csv",
    f"{_M}ci-large-tank-d76.csv",
    *(f"{_MT}development_times/normal-small-tank-rolls/d-76@{t}" for t in (65, 70, 72, 75)),
    *(f"{_MT}development_times/normal-large-tank-rolls/d-76@{t}" for t in (65, 70, 72, 75)),
    *(
        f"{_MT}development_times/{table}/d-76@{t}"
        for table in (
            "normal-large-tank-sheets",
            "normal-tray-sheets",
            "normal-rotary-tube-rolls",
            "normal-rotary-tube-sheets",
        )
        for t in (65, 68, 70, 72, 75)
    ),
    # Mixed entries (spec 7.6 "held-out granularity").
    f"{_MT}facts/base-plus-fog/char-d76-06min",
]

_WITHIN = [
    *(f"{_XT}development_times/normal-small-tank/d-76/1+1@{t}" for t in (65, 70, 72, 75)),
    *(f"{_XT}development_times/normal-small-tank/xtol/1+1@{t}" for t in (65, 70, 72, 75)),
    *(f"{_XT}development_times/normal-large-tank/d-76/1+1@{t}" for t in (65, 70, 72, 75)),
    *(f"{_XT}development_times/normal-large-tank/xtol/1+1@{t}" for t in (65, 70, 72, 75)),
    *(
        f"{_XT}xtol_contrast_index_times/ci-{c}/xtol/1+1@{t}"
        for c in ("0.72", "0.82")
        for t in (70, 75)
    ),
    *(f"{_MT}development_times/normal-small-tank-rolls/d-76/1+1@{t}" for t in (65, 70, 72, 75)),
    f"{_XT}reciprocity/1e-05s",
    f"{_XT}reciprocity/100s",
    f"{_XT}reciprocity/1e-05s/development",
    f"{_XT}reciprocity/100s/development",
    *(f"{_XT}reciprocity/{t}/development" for t in ("0.0001s", "0.001s", "1s", "10s")),
    f"{_MT}reciprocity/100s",
    *(f"{_MT}development_time_adjustment_factors/d-76/normal@{t}" for t in (65, 68, 70, 72, 75)),
    f"{_XT}facts/rms-granularity",
    "developers/d-76.yaml#facts/capacity-without-replenishment",
    "developers/d-76.yaml#facts/dilution-1-1-two-rolls-16-ounce-tank",
    "developers/xtol.yaml#facts/dilution-1-1-minimum-stock",
    "developers/xtol.yaml#facts/time-compensation",
]

_UNUSED = [
    f"{_XT}development_times/normal-small-tank/d-76@68",
    f"{_XT}development_times/normal-small-tank/d-76/1+1@68",
    f"{_XT}development_times/normal-small-tank/xtol@68",
    f"{_XT}development_times/normal-small-tank/xtol/1+1@68",
    f"{_XT}development_times/normal-small-tank/hc-110@68",
    f"{_XT}development_times/normal-small-tank/t-max-dev@68",
    f"{_MT}development_times/normal-small-tank-rolls/d-76@68",
    f"{_MT}development_times/normal-small-tank-rolls/d-76/1+1@68",
    f"{_MT}development_times/normal-large-tank-rolls/d-76@68",
    *(f"{_XT}xtol_contrast_index_times/ci-0.56/xtol@{t}" for t in (65, 68, 70, 75, 80)),
    *(f"{_XT}xtol_contrast_index_times/ci-0.56/xtol/1+1@{t}" for t in (68, 70, 75, 80)),
    f"{_XT}xtol_contrast_index_times/ci-0.72/xtol@68",
    f"{_XT}xtol_contrast_index_times/ci-0.82/xtol@68",
    # T-MAX RS.
    f"{_M}char-tmax-rs-08min.csv",
    f"{_M}char-tmax-rs-10.5min.csv",
    f"{_M}ci-small-tank-tmax-rs.csv",
    f"{_M}ci-large-tank-tmax-rs.csv",
    f"{_X}ci-small-tank-tmax-rs.csv",
    f"{_MT}facts/base-plus-fog/char-tmax-rs-08min",
    f"{_XT}development_times/normal-small-tank/t-max-rs@68",
    f"{_MT}development_times/normal-small-tank-rolls/t-max-rs/1+7@68",
    f"{_MT}replenishment/t-max-rs",
    "developers/t-max-rs.yaml#facts/form",
    # Push tables of the training pairs.
    f"{_XT}development_times/push-ei1600-small-tank/d-76@68",
    f"{_XT}development_times/push-ei3200-rotary-tube/xtol/1+1@75",
    f"{_MT}push_processing/push-small-tank-rolls/d-76",
]

_HELD_OUT_EVALUATED = [
    # T-Max x XTOL.
    *(
        f"{_M}ci-{tank}-xtol{d}.csv"
        for tank, d in (
            ("small-tank", ""),
            ("small-tank", "-1-1"),
            ("large-tank", ""),
            ("rotary-tube", ""),
            ("rotary-tube", "-1-1"),
            ("tray", ""),
            ("tray", "-1-1"),
        )
    ),
    f"{_MT}development_times/normal-small-tank-rolls/xtol@68",
    f"{_MT}development_times/normal-tray-sheets/xtol/1+1@75",
    f"{_MT}xtol_contrast_index_times/ci-0.52/xtol@65",
    f"{_MT}xtol_contrast_index_times/ci-0.62/xtol/1+1@80",
    # T-Max x HC-110.
    f"{_M}ci-small-tank-hc110-b.csv",
    f"{_M}ci-large-tank-hc110-b.csv",
    f"{_MT}development_times/normal-large-tank-rolls/hc-110@70",
    f"{_MT}development_time_adjustment_factors/hc-110/normal@68",
    # T-Max x T-MAX Developer.
    *(f"{_M}char-tmax-dev-{m}min.csv" for m in ("06", "07", "10", "12")),
    f"{_M}ci-small-tank-tmax-dev.csv",
    f"{_MT}development_times/normal-small-tank-rolls/t-max-dev/1+4@68",
    f"{_MT}development_times/normal-rotary-tube-rolls/t-max-dev/1+9@75",
    f"{_MT}development_time_adjustment_factors/t-max-dev+t-max-rs/normal@75",
]

_HELD_OUT_NOT_EVALUATED = [
    f"{_MT}facts/characteristic-vs-contrast-index-curves",
    f"{_MT}facts/base-plus-fog/char-tmax-dev-06min",
    f"{_MT}push_processing/push-small-tank-rolls/xtol",
    f"{_MT}push_processing/push-small-tank-rolls/hc-110",
    f"{_MT}push_processing/push-rotary-tube-rolls/t-max-dev",
    f"{_MT}rack_and_tank/normal/xtol",
    f"{_MT}replenishment/xtol",
    *(
        f"{_MT}development_time_adjustment_factors/{b}/{c}@75"
        for b in ("hc-110", "t-max-dev+t-max-rs")
        for c in ("less-20", "more-20", "more-40")
    ),
    f"{_M}ci-small-tank-d76.csv#conditions/ci_cross_check",
    f"{_M}char-d76-07.5min.csv#conditions/ci_cross_check",
    # Free text quoting held-out cells (7.6 design-time exposure, 7.7).
    f"{_MT}development_times/normal-rotary-tube-sheets/~note",
    f"{_MT}development_times/normal-small-tank-rolls/~contrast_aim",
    f"{_MT}development_time_adjustment_factors/~note",
    f"{_MT}push_exposure/~note",
    f"{_MT}facts/base-plus-fog/~statement",
    "developers/hc-110.yaml#development_times/roll-dilution-b-small-tank/~note",
]

# Prefixes of developer-publication rows whose class 7.6 fixes as a group.
_DEVELOPER_ROW_GROUPS = [
    ("developers/xtol.yaml#development_times/table1-small-tank-full-strength/t-max-100/", "ho+"),
    ("developers/xtol.yaml#development_times/table4-large-tank-fresh/t-max-100/", "ho+"),
    ("developers/xtol.yaml#development_times/table8-large-tank-fresh-sheets/t-max-100/", "ho+"),
    ("developers/xtol.yaml#development_times/table2-rotary-tube-seasoned/t-max-100/", "ho-"),
    ("developers/xtol.yaml#development_times/table3-large-tank-seasoned/t-max-100/", "ho-"),
    (
        "developers/xtol.yaml#development_times/table6-rotary-tube-sheets-replenished/t-max-100/",
        "ho-",
    ),
    ("developers/xtol.yaml#development_times/table7-large-tank-seasoned-sheets/t-max-100/", "ho-"),
    ("developers/hc-110.yaml#development_times/roll-dilution-b-small-tank/t-max-100/", "ho+"),
    (
        "developers/t-max-dev.yaml#development_times/small-tank-rolls/t-max-100/roll/ei100-200",
        "ho+",
    ),
    ("developers/t-max-dev.yaml#development_times/small-tank-rolls/t-max-100/roll/ei400", "ho-"),
    ("developers/xtol.yaml#development_times/table1-small-tank-full-strength/tri-x-400/", "un"),
    ("developers/xtol.yaml#development_times/table2-rotary-tube-seasoned/tri-x-400/", "un"),
    ("developers/xtol.yaml#development_times/table3-large-tank-seasoned/tri-x-400/", "un"),
    ("developers/xtol.yaml#development_times/table4-large-tank-fresh/tri-x-400/", "un"),
    ("developers/hc-110.yaml#development_times/roll-dilution-b-small-tank/tri-x-400/", "un"),
    ("developers/t-max-dev.yaml#development_times/small-tank-rolls/tri-x-400/", "un"),
    ("developers/t-max-rs.yaml#development_times/small-tank-rolls/t-max-100/", "un"),
]


def test_spec_named_items_land_in_their_class(registry: splits.Registry, items: list[str]) -> None:
    known = set(items)
    expected = (
        [(i, registry.training) for i in _TRAINING]
        + [(i, registry.within_checks) for i in _WITHIN]
        + [(i, registry.unused) for i in _UNUSED]
        + [(i, registry.held_out_evaluated) for i in _HELD_OUT_EVALUATED]
        + [(i, registry.held_out - registry.held_out_evaluated) for i in _HELD_OUT_NOT_EVALUATED]
    )
    missing = [i for i, _ in expected if i not in known]
    assert not missing, f"named by the spec but not enumerated: {missing}"
    wrong = [(i, classify(i).split) for i, cls in expected if i not in cls]
    assert not wrong


def test_developer_publication_rows(registry: splits.Registry, items: list[str]) -> None:
    for prefix, tag in _DEVELOPER_ROW_GROUPS:
        ids = [i for i in items if i.startswith(prefix)]
        assert ids, prefix
        target = {
            "ho+": registry.held_out_evaluated,
            "ho-": registry.held_out - registry.held_out_evaluated,
            "un": registry.unused,
        }[tag]
        assert set(ids) <= target, prefix
    # Every T-Max row of the three held-out developers' publications is held out.
    for name in ("xtol", "hc-110", "t-max-dev"):
        rows = [
            i
            for i in items
            if i.startswith(f"developers/{name}.yaml#development_times/") and "/~" not in i
        ]
        tmax = [i for i in rows if "/t-max-100/" in i]
        assert tmax and set(tmax) <= registry.held_out
        assert set(rows) - set(tmax) <= registry.unused
    # Nothing of T-MAX RS is used anywhere.
    rs = [
        i for i in items if re.search(r"t-?max-rs", i.split("#")[0]) or re.search(r"/t-max-rs\b", i)
    ]
    rs = [i for i in rs if "t-max-dev+t-max-rs" not in i]
    assert rs and set(rs) <= registry.unused


def test_held_out_combinations_are_wholly_held_out(
    registry: splits.Registry, items: list[str]
) -> None:
    """Every T-Max row or cell naming xtol, hc-110 or t-max-dev is held out (spec 7.6)."""
    for i in items:
        if not i.startswith(_MT):
            continue
        sub = i.split("#", 1)[1]
        segments = sub.split("@")[0].split("/")
        if {"xtol", "hc-110", "t-max-dev"} & set(segments):
            assert i in registry.held_out, i
    assert all(i in registry.held_out for i in items if i.startswith(f"{_MT}xtol_contrast_index"))


def test_classify_returns_reasons() -> None:
    c = classify(f"{_MT}push_processing/push-small-tank-rolls/xtol")
    assert c.split == "held_out" and c.evaluated is False and c.reason
    c = classify(f"{_M}ci-tray-xtol.csv")
    assert c.split == "held_out" and c.evaluated is True
    c = classify(f"{_X}char-35mm-d76-06min.csv")
    assert c.split == "training" and c.evaluated is None
    assert classify(f"{_XT}development_times/normal-small-tank/d-76@68").reason


@pytest.mark.parametrize(
    "item_id",
    [
        "",
        "kodak/tri-x-400/unknown.csv",
        "kodak/ilford-hp5/char-d76-06min.csv",
        f"{_XT}development_times/normal-small-tank/d-76@69",
        f"{_XT}no_such_block/x",
        f"{_MT}development_times/normal-small-tank-rolls/rodinal@68",
        "developers/xtol.yaml#development_times/table9-unknown/t-max-100/135/ei100",
        "developers/rodinal.yaml#facts/x",
    ],
)
def test_classify_rejects_unregistered_ids(item_id: str) -> None:
    with pytest.raises(KeyError):
        classify(item_id)
    with pytest.raises(LeakageError):
        assert_trainable(item_id)


# --------------------------------------------------------------------------------------
# Leakage (spec D1 7.7)
# --------------------------------------------------------------------------------------


def test_assert_trainable(registry: splits.Registry) -> None:
    for item in sorted(registry.training)[:50]:
        assert_trainable(item)
    for group in (registry.held_out, registry.unused, registry.within_checks):
        for item in sorted(group):
            with pytest.raises(LeakageError):
                assert_trainable(item)
    with pytest.raises(LeakageError, match="held out"):
        assert_trainable(f"{_M}ci-small-tank-xtol.csv")
    with pytest.raises(TypeError):
        assert_trainable(None)  # type: ignore[arg-type]


def _library_copy(root: Path, tmp_path: Path) -> Library:
    shutil.copytree(root, tmp_path / "reference")
    return Library(tmp_path)


def _walk_ids(rel: str, data: dict) -> set[str]:
    return {item.id for item in splits.walk_yaml(rel, data)}


def test_training_library_filters_held_out(
    root: Path, registry: splits.Registry, tmp_path: Path
) -> None:
    lib = TrainingLibrary(_library_copy(root, tmp_path))

    for rel in sorted(i for i in registry.held_out if i.endswith(".csv")):
        with pytest.raises(LeakageError):
            lib.reference_curve(rel)
    # Training curves load, minus the held-out header field.
    curve = lib.reference_curve(f"{_M}ci-small-tank-d76.csv")
    assert curve.conditions and "ci_cross_check" not in curve.conditions
    assert "developer=d-76" in curve.conditions
    assert lib.reference_curve(f"{_X}char-35mm-d76-08min.csv").x.size > 10

    tables = lib.reference_tables("kodak/t-max-100")
    ids = _walk_ids("kodak/t-max-100/tables.yaml", tables)
    assert ids and not (ids & registry.held_out)
    assert ids >= {i for i in registry.training if i.startswith(_MT)}
    assert "xtol_contrast_index_times" not in tables
    facts = {f["id"] for f in tables["facts"]}
    assert "characteristic-vs-contrast-index-curves" not in facts
    assert "rms-granularity" in facts
    developers = {r["developer"] for t in tables["development_times"] for r in t["rows"]}
    assert developers.isdisjoint({"xtol", "hc-110", "t-max-dev"})
    # The base library is untouched.
    base = Library(tmp_path).reference_tables("kodak/t-max-100")
    assert "xtol_contrast_index_times" in base

    tx = lib.reference_tables("kodak/tri-x-400")
    assert _walk_ids("kodak/tri-x-400/tables.yaml", tx) == _walk_ids(
        "kodak/tri-x-400/tables.yaml", Library(tmp_path).reference_tables("kodak/tri-x-400")
    )

    for bath in ("xtol", "hc-110-pre2019", "t-max-dev"):
        facts = lib.developer_facts(bath)
        materials = {
            r.get("material") for t in facts.get("development_times", []) for r in t["rows"]
        }
        assert "kodak/t-max-100" not in materials
    assert lib.developer_facts("d-76")["facts"]


def _strings(node: object, path: tuple = ()) -> list[tuple[tuple, str]]:
    if isinstance(node, dict):
        return [s for k, v in node.items() for s in _strings(v, (*path, k))]
    if isinstance(node, list):
        return [s for i, v in enumerate(node) for s in _strings(v, (*path, i))]
    return [(path, node)] if isinstance(node, str) else []


# Names of the developers held out for T-Max 100. "T-MAX RS" (unused, kept) does not match.
_HELD_OUT_NAMES = re.compile(r"xtol|hc-?110|t-?max dev(eloper)?\b", re.IGNORECASE)
_TITLE_KEYS = {"table", "push_table", "where"}  # publication table titles, no values


def test_training_library_free_text_has_no_held_out_values(root: Path, tmp_path: Path) -> None:
    """Notes, statements, footnotes and contrast aims quoting held-out cells are removed (7.7)."""
    lib = TrainingLibrary(_library_copy(root, tmp_path))
    tables = lib.reference_tables("kodak/t-max-100")
    leaks = [
        (p, s[:80])
        for p, s in _strings(tables)
        if p[0] != "sources" and p[-1] not in _TITLE_KEYS and _HELD_OUT_NAMES.search(s)
    ]
    assert not leaks
    # The specific quotes found in review are gone; training cells and conditions remain.
    flat = " ".join(s for _, s in _strings(tables))
    for quote in ("8¾, 7¾, 6¾", "1.2 x 6¼", "gives 7¾ min", "CI 0.56", "-3.15 for D-76 and T-MAX"):
        assert quote not in flat, quote
    small = next(t for t in tables["development_times"] if t["id"] == "normal-small-tank-rolls")
    assert "agitation" in small and "contrast_aim" not in small and "footnotes" not in small
    assert {r["developer"] for r in small["rows"]} == {"t-max-rs", "d-76"}
    for name in ("xtol", "hc-110-pre2019", "t-max-dev"):
        facts = lib.developer_facts(name)
        for t in facts.get("development_times", []):
            assert "note" not in t and "footnotes" not in t, (name, t["id"])
    # Tri-X free text is unused, not held out, and passes through.
    tx = lib.reference_tables("kodak/tri-x-400")
    assert tx["development_times"][0]["note"]


def test_free_text_of_units_with_held_out_items_is_held_out(
    root: Path, registry: splits.Registry
) -> None:
    """A note or footnote beside held-out cells is held out with them (7.6, 7.7)."""
    for path in sorted(root.rglob("*.yaml")):
        rel = path.relative_to(root).as_posix()
        walked = splits.walk_yaml(rel, yaml.safe_load(path.read_text(encoding="utf-8")))
        held_units = {u for it in walked if it.id in registry.held_out for u in it.units}
        for it in walked:
            if "/~" in it.id and it.units[-1] in held_units:
                assert it.id in registry.held_out, it.id
                assert it.id not in registry.held_out_evaluated, it.id
                assert registry.not_evaluated_reason[it.id], it.id


def test_free_text_ids() -> None:
    items = splits.walk_yaml(
        "kodak/tri-x-400/tables.yaml",
        {"reciprocity": {"note": "n", "statement": "s", "rows": []}, "push_exposure": {}},
    )
    assert [i.id for i in items] == [
        "kodak/tri-x-400/tables.yaml#reciprocity/~statement",
        "kodak/tri-x-400/tables.yaml#reciprocity/~note",
    ]
    assert items[0].fields == (("reciprocity", "statement"),)
    with pytest.raises(KeyError):
        classify(f"{_MT}development_times/normal-small-tank-rolls/~unknown key!")
    with pytest.raises(LeakageError, match="held out"):
        assert_trainable(f"{_MT}push_exposure/~note")


def test_training_library_refuses_unregistered_curves(root: Path, tmp_path: Path) -> None:
    lib = TrainingLibrary(_library_copy(root, tmp_path))
    for rel in (f"{_M}CI-Small-Tank-XTOL.csv", f"{_M}./ci-small-tank-xtol.csv", "x/y.csv"):
        with pytest.raises(LeakageError):
            lib.reference_curve(rel)


def test_training_library_requires_library() -> None:
    with pytest.raises(TypeError):
        TrainingLibrary("data")  # type: ignore[arg-type]


# --------------------------------------------------------------------------------------
# Hash (spec D1 8.4)
# --------------------------------------------------------------------------------------


def test_registry_hash_stable(root: Path, registry: splits.Registry, tmp_path: Path) -> None:
    digest = registry.sha256()
    assert len(digest) == 64 and int(digest, 16) >= 0
    assert digest == build_registry(root).sha256()
    assert splits.registry_sha256(root) == digest

    code = (
        "from darkroom_atlas.validate.splits import registry_sha256; "
        f"print(registry_sha256({str(root)!r}))"
    )
    seen = set()
    for seed in ("1", "2"):
        env = {**os.environ, "PYTHONHASHSEED": seed}
        out = subprocess.run(
            [sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True
        )
        seen.add(out.stdout.strip())
    assert seen == {digest}

    # The hash does not depend on where the reference data live.
    moved = tmp_path / "data" / "reference"
    shutil.copytree(root, moved)
    assert splits.registry_sha256(moved) == digest

    # It does depend on every class and every reason.
    some = next(iter(sorted(registry.not_evaluated_reason)))
    reasons = {**registry.not_evaluated_reason, some: "changed"}
    assert replace(registry, not_evaluated_reason=reasons).sha256() != digest
    moved_item = next(iter(sorted(registry.training)))
    assert (
        replace(
            registry,
            training=registry.training - {moved_item},
            within_checks=registry.within_checks | {moved_item},
        ).sha256()
        != digest
    )


# --------------------------------------------------------------------------------------
# Pair-metric table (spec D1 8.3)
# --------------------------------------------------------------------------------------


def test_pair_metrics_reference_registered_items(registry: splits.Registry) -> None:
    tiers = Counter(splits.PAIR_TIERS.values())
    assert tiers == {"training": 5, "held_out": 3, "predicted": 6}
    for pair, metrics in splits.PAIR_METRICS.items():
        tier = splits.PAIR_TIERS[pair]
        if tier == "predicted":
            assert metrics == ()
            continue
        assert metrics and len(set(metrics)) == len(metrics)
        pool = registry.training if tier == "training" else registry.held_out_evaluated
        for metric in metrics:
            name, _, target = metric.partition(":")
            assert name in splits.METRIC_NAMES, metric
            if target.endswith(".csv"):
                assert target in pool, metric
    # Counts from the pair-metric table.
    key = ("kodak/tri-x-400", "d-76")
    assert sum(m.startswith("density_rms:") for m in splits.PAIR_METRICS[key]) == 8
    assert sum(m.startswith("ci_at_time:") for m in splits.PAIR_METRICS[key]) == 2
    assert sum(m.startswith("log_h_m:") for m in splits.PAIR_METRICS[key]) == 2
    key = ("kodak/t-max-100", "xtol")
    assert sum(m.startswith("ci_at_time:") for m in splits.PAIR_METRICS[key]) == 7
    key = ("kodak/t-max-100", "t-max-dev")
    assert sum(m.startswith("toe_density:") for m in splits.PAIR_METRICS[key]) == 4
    assert ("kodak/t-max-100", "hc-110-pre2019") in splits.PAIR_METRICS


# --------------------------------------------------------------------------------------
# Bad input
# --------------------------------------------------------------------------------------


def test_enumerate_rejects_missing_root(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        enumerate_items(tmp_path / "missing")


def test_walk_rejects_unknown_blocks_and_duplicates() -> None:
    with pytest.raises(ValueError, match="unknown"):
        splits.walk_yaml("kodak/tri-x-400/tables.yaml", {"schema_version": 1, "surprise": []})
    with pytest.raises(ValueError, match="unknown"):
        splits.walk_yaml("developers/xtol.yaml", {"schema_version": 1, "surprise": {}})
    dup = {"facts": [{"id": "a", "value": 1}, {"id": "a", "value": 2}]}
    with pytest.raises(ValueError, match="duplicate"):
        splits.walk_yaml("kodak/tri-x-400/tables.yaml", dup)
    bad_cells = {
        "development_times": [
            {
                "id": "t",
                "temperature_f": [65, 68],
                "rows": [{"developer": "d-76", "dilution": "stock", "minutes": [1.0]}],
            }
        ]
    }
    with pytest.raises(ValueError, match="cells"):
        splits.walk_yaml("kodak/tri-x-400/tables.yaml", bad_cells)
    with pytest.raises(ValueError):
        splits.walk_yaml("kodak/tri-x-400/other.yaml", {})


def test_enumerate_flags_unclassified_file(root: Path, tmp_path: Path) -> None:
    copy = tmp_path / "reference"
    shutil.copytree(root, copy)
    (copy / "kodak" / "tri-x-400" / "char-35mm-rodinal-08min.csv").write_text(
        "# material: kodak/tri-x-400\n0,0\n1,1\n", encoding="utf-8"
    )
    assert "kodak/tri-x-400/char-35mm-rodinal-08min.csv" in enumerate_items(copy)
    with pytest.raises(KeyError):
        build_registry(copy)
