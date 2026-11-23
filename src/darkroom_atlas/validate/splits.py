"""Train / within-check / held-out / unused registry of the reference data (spec D1 7.6, 7.7).

Every reference file, table row and sub-item has an id: its path relative to the reference
root (``data/reference/``, spec 5.9), plus for
YAML rows ``#<block>/<id>/<developer>[/<dilution>][@<temperature F>]`` (spec 7.6). Examples::

    kodak/tri-x-400/char-35mm-d76-08min.csv
    kodak/tri-x-400/tables.yaml#development_times/normal-small-tank/d-76/1+1@70
    kodak/t-max-100/tables.yaml#facts/base-plus-fog/char-tmax-dev-06min
    kodak/t-max-100/ci-small-tank-d76.csv#conditions/ci_cross_check
    kodak/t-max-100/tables.yaml#development_times/normal-rotary-tube-sheets/~note
    developers/xtol.yaml#development_times/table4-large-tank-fresh/t-max-100/120/ei100/ci0.56

:func:`classify` maps an id to its class by explicit rules over ids; an id that no rule
covers raises ``KeyError``, and two rules covering one id is a bug. :func:`enumerate_items`
walks a reference root, :func:`build_registry` classifies everything it finds. The module
constants ``TRAINING``, ``WITHIN_CHECKS``, ``HELD_OUT``, ``UNUSED``, ``HELD_OUT_EVALUATED``
and ``NOT_EVALUATED_REASON`` are the registry of the default root, built on first access.
"""

import copy
import functools
import hashlib
import json
import math
import re
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from types import MappingProxyType
from typing import Any, Literal

import yaml

from darkroom_atlas.data.library import (
    DEVELOPER_FACTS_FILES,
    Library,
    default_data_root,
)
from darkroom_atlas.data.models import ReferenceCurve

__all__ = [
    "HELD_OUT",
    "HELD_OUT_EVALUATED",
    "METRIC_NAMES",
    "NOT_EVALUATED_REASON",
    "PAIR_METRICS",
    "PAIR_TIERS",
    "TRAINING",
    "UNUSED",
    "WITHIN_CHECKS",
    "Classification",
    "Item",
    "LeakageError",
    "Registry",
    "TrainingLibrary",
    "assert_trainable",
    "build_registry",
    "classify",
    "default_reference_root",
    "enumerate_items",
    "registry_sha256",
    "walk_yaml",
]

Split = Literal["training", "within_checks", "held_out", "unused"]
_SPLITS: tuple[Split, ...] = ("training", "within_checks", "held_out", "unused")


class LeakageError(RuntimeError):
    """A held-out (or otherwise non-training) reference item reached a training path (7.7)."""


def default_reference_root() -> Path:
    """``<data root>/reference`` (spec 5.9)."""
    return default_data_root() / "reference"


# --------------------------------------------------------------------------------------
# Rules (spec 7.6)
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Classification:
    """Class of one item; ``evaluated`` only for held-out items, ``reason`` when required."""

    split: Split
    evaluated: bool | None = None
    reason: str | None = None


@dataclass(frozen=True)
class _Rule:
    pattern: re.Pattern[str]
    result: Classification


def _rule(regex: str, split: Split, reason: str | None = None, *, evaluated: bool | None = None):
    if split == "held_out" and evaluated is None:
        raise ValueError(f"held-out rule {regex!r} must say whether it is evaluated")
    if (split == "unused" or evaluated is False) and not reason:
        raise ValueError(f"rule {regex!r} needs a reason")
    return _Rule(re.compile(regex), Classification(split, evaluated, reason))


_X = r"kodak/tri-x-400/"
_M = r"kodak/t-max-100/"
_XT = r"kodak/tri-x-400/tables\.yaml#"
_MT = r"kodak/t-max-100/tables\.yaml#"
_T4 = r"@(65|70|72|75)"  # the non-20 C columns of the F-4016/F-4017 tables
_ANY_T = r"@\d+"
_TMAX_RS = r"/t-max-rs(/1\+\d+)?"

_R_RS = "T-MAX RS: D1 does not use T-MAX RS data (7.6)"
_R_CI_CURVE = "duplicates the small-tank CI(t) curve at 20 C (7.2, 7.6)"
_R_PUSH_TRAIN = "push table of a training pair: no CI target (7.2, 7.6)"
_R_PUSH_HELD = "push processing: no CI target (7.6)"
_R_CONTRAST_COLS = (
    "contrast columns of the temperature-factor block: no pre-registered CI target (7.2)"
)
_R_NOT_OBS = "not an observation of any D1 metric (7.2)"
_R_PROSE = "free text (note, statement, footnotes): not an observation of any D1 metric (7.2)"
_R_PROSE_HELD = (
    "free text of a unit with held-out items; it may quote held-out cells, so it is removed "
    "with them (7.6 design-time exposure, 7.7)"
)
_R_DEV_DUP = (
    "developer-publication copy of 400TX rows: cells shared with F-4017 are identical (reference "
    "developers/README.md cross-check) and the rest are not selected (7.6)"
)

_RULES: tuple[_Rule, ...] = (
    # ---- documentation -----------------------------------------------------------------
    _rule(r"developers/README\.md", "unused", "documentation of the reference files, not data"),
    # ---- curve files: Tri-X ------------------------------------------------------------
    _rule(_X + r"char-(35mm|120)-d76-(06|08|10|12)min\.csv", "training"),
    _rule(_X + r"char-(35mm|120)-tmax-dev-(06|07|09|11)min\.csv", "training"),
    _rule(_X + r"mtf-d76-large-tank\.csv", "training"),
    _rule(_X + r"spectral-sensitivity-d(0\.3|1\.0)\.csv", "training"),
    _rule(_X + r"ci-small-tank-(d76|d76-1-1|xtol|xtol-1-1|hc110-b|tmax-dev)\.csv", "training"),
    # ---- curve files: T-Max ------------------------------------------------------------
    _rule(_M + r"char-d76-(06|07\.5|10)min\.csv", "training"),
    _rule(_M + r"mtf-d76-small-tank\.csv", "training"),
    _rule(_M + r"spectral-sensitivity-d(0\.3|1\.0)\.csv", "training"),
    _rule(_M + r"ci-(small-tank-d76|small-tank-d76-1-1|large-tank-d76)\.csv", "training"),
    _rule(
        _M + r"ci-(small-tank|large-tank|rotary-tube|tray)-xtol(-1-1)?\.csv",
        "held_out",
        evaluated=True,
    ),
    _rule(_M + r"ci-(small|large)-tank-hc110-b\.csv", "held_out", evaluated=True),
    _rule(_M + r"char-tmax-dev-(06|07|10|12)min\.csv", "held_out", evaluated=True),
    _rule(_M + r"ci-small-tank-tmax-dev\.csv", "held_out", evaluated=True),
    _rule(
        r"kodak/(tri-x-400|t-max-100)/(char-tmax-rs-[0-9.]+min|ci-[a-z-]+-tmax-rs)\.csv",
        "unused",
        _R_RS,
    ),
    # ---- Tri-X tables.yaml: facts --------------------------------------------------------
    _rule(_XT + r"facts/processing-aim-contrast-index", "training"),
    _rule(_XT + r"facts/rms-granularity", "within_checks"),
    _rule(_XT + r"facts/exposure-index", "within_checks"),
    _rule(
        _XT + r"facts/base-plus-fog/char-(35mm|120)-(d76-(06|08|10|12)|tmax-dev-(06|07|09|11))min",
        "training",
    ),
    _rule(
        _XT + r"facts/(short-development-warning|single-use-camera-film-speed|"
        r"base-type-and-thickness|formats)",
        "unused",
        _R_NOT_OBS,
    ),
    # ---- Tri-X tables.yaml: reciprocity (7.6; the corrections are checks, 6.8) ---------------
    _rule(_XT + r"reciprocity/(0\.0001|0\.001|0\.01|0\.1|1|10)s", "training"),
    _rule(_XT + r"reciprocity/(1e-05|100)s", "within_checks"),
    _rule(_XT + r"reciprocity/(1e-05|0\.0001|0\.001|1|10|100)s/development", "within_checks"),
    _rule(
        _XT + r"reciprocity/(0\.01|0\.1)s/development",
        "unused",
        "no development change printed for this row; 7.6 checks only the +20 % ... -30 % rows",
    ),
    # ---- Tri-X tables.yaml: development times ------------------------------------------
    _rule(
        _XT + r"development_times/normal-small-tank/(d-76|xtol|hc-110|t-max-dev)" + _T4, "training"
    ),
    _rule(_XT + r"development_times/normal-small-tank/(d-76|xtol)/1\+1" + _T4, "within_checks"),
    _rule(
        _XT + r"development_times/normal-small-tank/((d-76|xtol)(/1\+1)?|hc-110|t-max-dev)@68",
        "unused",
        _R_CI_CURVE,
    ),
    _rule(
        _XT + r"development_times/normal-large-tank/(d-76|xtol|hc-110)@(65|68|70|72|75)", "training"
    ),
    _rule(_XT + r"development_times/normal-large-tank/(d-76|xtol)/1\+1@68", "training"),
    _rule(_XT + r"development_times/normal-large-tank/(d-76|xtol)/1\+1" + _T4, "within_checks"),
    _rule(
        _XT + r"development_times/normal-large-tank/t-max-dev" + _ANY_T,
        "unused",
        "printed NR (not recommended) in F-4017",
    ),
    _rule(
        _XT + r"development_times/normal-rotary-tube/(d-76|xtol|hc-110|t-max-dev)(/1\+1)?" + _ANY_T,
        "unused",
        "every cell equals the small-tank table (7.2)",
    ),
    _rule(
        _XT + r"development_times/push-[a-z0-9-]+/(d-76|xtol|hc-110|t-max-dev)(/1\+1)?" + _ANY_T,
        "unused",
        _R_PUSH_TRAIN,
    ),
    _rule(_XT + r"development_times/[a-z0-9-]+" + _TMAX_RS + _ANY_T, "unused", _R_RS),
    # ---- Tri-X tables.yaml: J-109 table 1 Tri-X rows (7.6) ---------------------------------
    _rule(
        _XT + r"xtol_contrast_index_times/ci-0\.56/xtol(/1\+1)?" + _ANY_T,
        "unused",
        "same cells as the F-4017 small-tank XTOL rows (reference README cross-check); "
        "counted once (7.6)",
    ),
    _rule(_XT + r"xtol_contrast_index_times/ci-0\.(72|82)/xtol@(65|70|75)", "training"),
    _rule(_XT + r"xtol_contrast_index_times/ci-0\.(72|82)/xtol/1\+1@(70|75)", "within_checks"),
    _rule(
        _XT + r"xtol_contrast_index_times/ci-0\.(72|82)/xtol(/1\+1)?@68",
        "unused",
        "68 F cell readable from the CI(t) curves, which span CI 0.34-0.98 (7.6)",
    ),
    _rule(
        _XT + r"xtol_contrast_index_times/ci-0\.(72|82)/xtol(/1\+1)?@80",
        "unused",
        "80 F cell: outside the F-4017 temperature range, not selected by 7.6",
    ),
    # ---- Tri-X tables.yaml: other blocks --------------------------------------------------
    _rule(
        _XT + r"rack_and_tank/normal/(d-76|xtol)",
        "unused",
        "rack-and-tank machine processing: not fitted (7.2)",
    ),
    _rule(_XT + r"rack_and_tank/normal/t-max-rs", "unused", _R_RS),
    _rule(_XT + r"developer_definitions/[a-z0-9-]+", "training"),
    _rule(
        _XT + r"(filter_factors|daylight_exposure)/[a-z0-9+-]+",
        "unused",
        "filter factors and daylight exposure tables are not used (7.2, 7.6)",
    ),
    # ---- T-Max tables.yaml: facts ---------------------------------------------------------
    _rule(_MT + r"facts/(processing-aim|rms-granularity)", "training"),
    _rule(_MT + r"facts/(exposure-index|iso-speed)", "within_checks"),
    _rule(
        _MT + r"facts/t-max-developers-contrast-aim",
        "held_out",
        "CI aim of the held-out T-MAX Developer cells; used as their target, not scored itself",
        evaluated=False,
    ),
    _rule(
        _MT + r"facts/characteristic-vs-contrast-index-curves",
        "held_out",
        "design-time reference containing T-MAX numbers; sigma_floor is recomputed from the "
        "training CSVs (7.2, 7.6)",
        evaluated=False,
    ),
    _rule(_MT + r"facts/base-plus-fog/char-d76-(06|07\.5|10)min", "training"),
    _rule(
        _MT + r"facts/base-plus-fog/char-tmax-dev-(06|07|10|12)min",
        "held_out",
        "duplicates the left end of the held-out curve, which is scored from the CSV",
        evaluated=False,
    ),
    _rule(_MT + r"facts/base-plus-fog/char-tmax-rs-[0-9.]+min", "unused", _R_RS),
    _rule(
        _MT + r"facts/(short-development-warning|development-adjustment-advice|"
        r"resolving-power-toc-1\.6|resolving-power-toc-1000|base-type|base-thickness-120|"
        r"base-thickness-135-and-sheets|sizes|spectral-sensitivity-blue-note)",
        "unused",
        _R_NOT_OBS,
    ),
    # ---- T-Max tables.yaml: reciprocity ---------------------------------------------------
    _rule(_MT + r"reciprocity/(0\.0001|0\.001|0\.01|0\.1|1|10)s", "training"),
    _rule(_MT + r"reciprocity/100s", "within_checks"),
    # ---- T-Max tables.yaml: development times ---------------------------------------------
    _rule(_MT + r"development_times/normal-small-tank-rolls/d-76" + _T4, "training"),
    _rule(_MT + r"development_times/normal-small-tank-rolls/d-76/1\+1" + _T4, "within_checks"),
    _rule(
        _MT + r"development_times/normal-small-tank-rolls/d-76(/1\+1)?@68", "unused", _R_CI_CURVE
    ),
    _rule(_MT + r"development_times/normal-large-tank-rolls/d-76" + _T4, "training"),
    _rule(
        _MT + r"development_times/normal-large-tank-rolls/d-76@68",
        "unused",
        "duplicates ci-large-tank-d76.csv at 20 C (7.6)",
    ),
    _rule(
        _MT + r"development_times/(normal-large-tank-sheets|normal-tray-sheets|"
        r"normal-rotary-tube-rolls|normal-rotary-tube-sheets)/d-76" + _ANY_T,
        "training",
    ),
    _rule(
        _MT + r"development_times/[a-z0-9-]+/(xtol|hc-110|t-max-dev)(/1\+\d+)?" + _ANY_T,
        "held_out",
        evaluated=True,
    ),
    _rule(_MT + r"development_times/[a-z0-9-]+" + _TMAX_RS + _ANY_T, "unused", _R_RS),
    # ---- T-Max tables.yaml: F-4016 temperature factors (7.2, 7.6) ----------------------------
    _rule(_MT + r"development_time_adjustment_factors/d-76/normal" + _ANY_T, "within_checks"),
    _rule(
        _MT + r"development_time_adjustment_factors/d-76/(less-20|more-20|more-40)" + _ANY_T,
        "unused",
        _R_CONTRAST_COLS,
    ),
    _rule(
        _MT + r"development_time_adjustment_factors/(hc-110|t-max-dev\+t-max-rs)/normal" + _ANY_T,
        "held_out",
        evaluated=True,
    ),
    _rule(
        _MT + r"development_time_adjustment_factors/(hc-110|t-max-dev\+t-max-rs)/"
        r"(less-20|more-20|more-40)" + _ANY_T,
        "held_out",
        _R_CONTRAST_COLS,
        evaluated=False,
    ),
    # ---- T-Max tables.yaml: push, machine, replenishment -----------------------------------
    _rule(
        _MT + r"push_exposure/[a-z0-9+-]+",
        "unused",
        "push exposure indices: push processing has no CI target (7.2, 7.6)",
    ),
    _rule(
        _MT + r"push_processing/[a-z0-9-]+/(xtol|hc-110|t-max-dev)(/1\+\d+)?",
        "held_out",
        _R_PUSH_HELD,
        evaluated=False,
    ),
    _rule(_MT + r"push_processing/[a-z0-9-]+/d-76", "unused", _R_PUSH_TRAIN),
    _rule(_MT + r"push_processing/[a-z0-9-]+" + _TMAX_RS, "unused", _R_RS),
    _rule(
        _MT + r"rack_and_tank/(normal|push-ei\d+)/xtol",
        "held_out",
        "machine-processing (rack-and-tank) time range (7.6)",
        evaluated=False,
    ),
    _rule(_MT + r"rack_and_tank/(normal|push-ei\d+)/t-max-rs", "unused", _R_RS),
    _rule(
        _MT + r"replenishment/xtol",
        "held_out",
        "replenishment (1.2 non-goal, 7.6)",
        evaluated=False,
    ),
    _rule(_MT + r"replenishment/t-max-rs", "unused", _R_RS),
    # ---- T-Max tables.yaml: J-109 table 1 T-Max rows, definitions, filters -------------------
    _rule(
        _MT + r"xtol_contrast_index_times/ci-\d\.\d\d/xtol(/1\+1)?" + _ANY_T,
        "held_out",
        evaluated=True,
    ),
    _rule(
        _MT + r"developer_definitions/(hc-110-dilution-b|t-max-developer-dilution)",
        "held_out",
        "dilution definition of a held-out developer; defines conditions, not scored",
        evaluated=False,
    ),
    _rule(_MT + r"filter_factors/[a-z0-9+-]+", "unused", "filter factors are not used (7.2, 7.6)"),
    # ---- free text of blocks, tables and facts ("<unit>/~<key>", 7.6, 7.7) -------------------
    # T-Max: every unit that holds held-out cells, plus push_exposure, whose note quotes them.
    _rule(
        _MT + r"((development_times|push_processing)/[a-z0-9-]+|facts/base-plus-fog|"
        r"development_time_adjustment_factors|push_exposure|rack_and_tank|replenishment|"
        r"xtol_contrast_index_times)/~[a-z_]+",
        "held_out",
        _R_PROSE_HELD,
        evaluated=False,
    ),
    _rule(_MT + r"(reciprocity|filter_factors)/~[a-z_]+", "unused", _R_PROSE),
    _rule(_XT + r"[a-z0-9_/+.-]+/~[a-z_]+", "unused", _R_PROSE),
    _rule(
        r"developers/(xtol|hc-110|t-max-dev)\.yaml#development_times/[a-z0-9-]+/~[a-z_]+",
        "held_out",
        _R_PROSE_HELD,
        evaluated=False,
    ),
    # d-76.yaml and t-max-rs.yaml free text falls under their whole-block rules below.
    # ---- CSV header sub-items (7.6 design-time exposure) -------------------------------------
    # handled in classify(): "<csv>#conditions/ci_cross_check"
    # ---- developer publications (7.6) --------------------------------------------------------
    _rule(r"developers/(d-76|xtol|hc-110|t-max-dev)\.yaml#facts/[a-z0-9.+-]+", "within_checks"),
    _rule(
        r"developers/(d-76|xtol|hc-110|t-max-dev)\.yaml#agitation/[a-z0-9-]+",
        "unused",
        "agitation definitions describe processes; each reference condition takes its "
        "agitation from the film publication",
    ),
    _rule(r"developers/t-max-rs\.yaml#.+", "unused", _R_RS),
    _rule(
        r"developers/d-76\.yaml#development_times/.+",
        "unused",
        "J-78 rows carry only legacy product names (material: null, 7.6)",
    ),
    _rule(
        r"developers/xtol\.yaml#development_times/(table1-small-tank-full-strength|"
        r"table1-small-tank-1-1|table2-rotary-tube-fresh-full-strength|"
        r"table2-rotary-tube-fresh-1-1|table4-large-tank-fresh|table5-tray-sheets-full-strength|"
        r"table5-tray-sheets-1-1|table6-rotary-tube-sheets-fresh|table6-rotary-tube-sheets-1-1|"
        r"table8-large-tank-fresh-sheets)/t-max-100/.+",
        "held_out",
        evaluated=True,
    ),
    _rule(
        r"developers/xtol\.yaml#development_times/(table2-rotary-tube-seasoned|"
        r"table6-rotary-tube-sheets-replenished)/t-max-100/.+",
        "held_out",
        "seasoned or replenished developer (1.2 non-goal)",
        evaluated=False,
    ),
    _rule(
        r"developers/xtol\.yaml#development_times/(table3-large-tank-seasoned|"
        r"table7-large-tank-seasoned-sheets)/t-max-100/.+",
        "held_out",
        "seasoned developer (1.2 non-goal) with nitrogen-burst agitation, which D1 does not model",
        evaluated=False,
    ),
    _rule(
        r"developers/xtol\.yaml#development_times/(table1-small-tank-full-strength|"
        r"table1-small-tank-1-1)/tri-x-400/.+",
        "unused",
        "J-109 table 1 Tri-X rows are used through tables.yaml#xtol_contrast_index_times",
    ),
    _rule(
        r"developers/xtol\.yaml#development_times/(table2-rotary-tube-fresh-full-strength|"
        r"table2-rotary-tube-fresh-1-1|table4-large-tank-fresh)/tri-x-400/.+",
        "unused",
        _R_DEV_DUP,
    ),
    _rule(
        r"developers/xtol\.yaml#development_times/table2-rotary-tube-seasoned/tri-x-400/.+",
        "unused",
        "seasoned or replenished developer (1.2 non-goal)",
    ),
    _rule(
        r"developers/xtol\.yaml#development_times/table3-large-tank-seasoned/tri-x-400/.+",
        "unused",
        "seasoned developer with nitrogen-burst agitation (1.2 non-goal; no D1 agitation class)",
    ),
    _rule(
        r"developers/hc-110\.yaml#development_times/[a-z0-9-]+/t-max-100/.+",
        "held_out",
        evaluated=True,
    ),
    _rule(
        r"developers/t-max-dev\.yaml#development_times/[a-z0-9-]+/t-max-100/[a-z]+/ei100-200",
        "held_out",
        evaluated=True,
    ),
    _rule(
        r"developers/t-max-dev\.yaml#development_times/[a-z0-9-]+/t-max-100/[a-z]+/ei(400|800)",
        "held_out",
        "push row (EI 400, 800): no CI target (7.6)",
        evaluated=False,
    ),
    _rule(
        r"developers/(hc-110|t-max-dev)\.yaml#development_times/[a-z0-9-]+/tri-x-400/.+",
        "unused",
        _R_DEV_DUP,
    ),
)

_CROSS_CHECK = "#conditions/ci_cross_check"


@functools.cache
def classify(item_id: str) -> Classification:
    """Class of ``item_id`` (spec 7.6); ``KeyError`` if no rule covers it."""
    if not isinstance(item_id, str):
        raise TypeError(f"item id must be a str, got {type(item_id).__name__}")
    if item_id.endswith(_CROSS_CHECK):
        parent = classify(item_id[: -len(_CROSS_CHECK)])
        if parent.split == "unused":
            return Classification("unused", None, f"header of an unused file: {parent.reason}")
        return Classification(
            "held_out",
            False,
            "header summary derived partly from held-out curves (7.6 design-time exposure)",
        )
    hits = {r.result for r in _RULES if r.pattern.fullmatch(item_id)}
    if not hits:
        raise KeyError(f"{item_id!r}: no rule of the split registry covers this item (7.6)")
    if len(hits) > 1:
        raise RuntimeError(f"{item_id!r}: split rules overlap: {sorted(map(str, hits))}")
    return hits.pop()


def assert_trainable(item_id: str) -> None:
    """Raise :class:`LeakageError` unless ``item_id`` is a training item (spec 7.7)."""
    try:
        c = classify(item_id)
    except KeyError as exc:
        raise LeakageError(f"{item_id!r} is not in the split registry") from exc
    if c.split == "held_out":
        raise LeakageError(f"{item_id!r} is held out (7.6); it must not enter any fit")
    if c.split != "training":
        raise LeakageError(f"{item_id!r} is {c.split}, not a training item (7.6)")


# --------------------------------------------------------------------------------------
# Walking the reference data
# --------------------------------------------------------------------------------------

_Path = tuple[str | int, ...]


@dataclass(frozen=True)
class Item:
    """One registered item of a YAML file and where it sits in the loaded data.

    ``units`` are the enclosing removable nodes, outermost first (block, table, row);
    ``fields`` are the nodes to remove (or, for cells of positional lists, set to ``None``)
    when the item alone is filtered out.
    """

    id: str
    units: tuple[_Path, ...]
    fields: tuple[_Path, ...]
    null_fields: bool = False


_TABLES_META = frozenset({"schema_version", "material", "film", "sources"})
_DEVELOPER_META = frozenset({"schema_version", "developer", "name", "films", "sources", "notes"})
_ADJUST_COLUMNS = {
    "20% Less Contrast": "less-20",
    "Normal Contrast": "normal",
    "20% More Contrast": "more-20",
    "40% More Contrast": "more-40",
}
_DILUTION = re.compile(r"^\s*1\+(\d+)\b")


def _slug(text: object) -> str:
    s = re.sub(r"[^a-z0-9+]+", "-", str(text).lower()).strip("-")
    if not s:
        raise ValueError(f"cannot form an id from {text!r}")
    return s


def _dev(row: Mapping[str, Any]) -> str:
    developer = row.get("developer")
    if not isinstance(developer, str) or not developer:
        raise ValueError(f"row without a developer: {row!r:.80}")
    m = _DILUTION.match(str(row.get("dilution") or ""))
    return f"{developer}/1+{m.group(1)}" if m else developer


def _temp(value: object) -> str:
    if isinstance(value, bool) or not isinstance(value, int | float) or value != int(value):
        raise ValueError(f"temperature {value!r} is not a whole number of degrees F")
    return f"@{int(value)}"


def _cells(
    rel_id: str,
    units: tuple[_Path, ...],
    row: Mapping,
    row_path: _Path,
    temps: list,
    keys: tuple[str, ...],
) -> Iterator[Item]:
    for key in keys:
        if key in row and (not isinstance(row[key], list) or len(row[key]) != len(temps)):
            raise ValueError(
                f"{rel_id}: {key} has {len(row[key] or [])} cells for {len(temps)} temperatures"
            )
    for j, t in enumerate(temps):
        fields = tuple((*row_path, k, j) for k in keys if k in row)
        yield Item(f"{rel_id}{_temp(t)}", units, fields, null_fields=True)


# Prose keys of a block, table or fact that holds rows or sub-items. Each is its own sub-item
# "<unit id>/~<key>" so that free text quoting held-out cells is classified and filtered
# like any other sub-item (7.6 design-time exposure, 7.7 "已去掉留出子项").
_PROSE_KEYS = ("statement", "note", "contrast_aim", "footnotes", "push_footnote")


def _prose(unit_id: str, units: tuple[_Path, ...], node: Mapping) -> Iterator[Item]:
    for key in _PROSE_KEYS:
        if key in node:
            yield Item(f"{unit_id}/~{key}", units, ((*units[-1], key),))


def _rows(block: Mapping, name: str) -> list:
    rows = block.get(name) or []
    if not isinstance(rows, list):
        raise TypeError(f"{name} is not a list")
    return rows


def _walk_tables(rel: str, data: Mapping[str, Any]) -> Iterator[Item]:
    for key, block in data.items():
        if key in _TABLES_META:
            continue
        top: _Path = (key,)
        base = f"{rel}#{key}/"
        if isinstance(block, Mapping):
            yield from _prose(f"{rel}#{key}", (top,), block)
        if key in ("facts", "developer_definitions"):
            for i, fact in enumerate(block):
                path = (key, i)
                subs = fact.get("left_end_of_characteristic_curves") if key == "facts" else None
                if subs:
                    yield from _prose(f"{base}{fact['id']}", (top, path), fact)
                    for j, sub in enumerate(subs):
                        p = (*path, "left_end_of_characteristic_curves", j)
                        stem = PurePosixPath(str(sub["curve"])).stem
                        yield Item(f"{base}{fact['id']}/{stem}", (top, path, p), (p,))
                else:
                    yield Item(f"{base}{fact['id']}", (top, path), (path,))
        elif key == "reciprocity":
            for i, row in enumerate(_rows(block, "rows")):
                path = (key, "rows", i)
                t = float(row["indicated_time_s"])
                if not (math.isfinite(t) and t > 0):
                    raise ValueError(f"{rel}: reciprocity time {t!r}")
                rid = f"{base}{t:g}s"
                dev_keys = ("development_adjustment", "development_change_percent")
                if row.get("development_adjustment") is not None:
                    own = tuple((*path, k) for k in row if k not in dev_keys)
                    yield Item(rid, (top, path), own)
                    dev = tuple((*path, k) for k in dev_keys if k in row)
                    yield Item(f"{rid}/development", (top, path), dev)
                else:
                    yield Item(rid, (top, path), (path,))
        elif key == "development_times":
            for ti, table in enumerate(block):
                temps = table["temperature_f"]
                yield from _prose(f"{base}{table['id']}", (top, (key, ti)), table)
                for ri, row in enumerate(_rows(table, "rows")):
                    path = (key, ti, "rows", ri)
                    units = (top, (key, ti), path)
                    yield from _cells(
                        f"{base}{table['id']}/{_dev(row)}",
                        units,
                        row,
                        path,
                        temps,
                        ("minutes", "printed", "bold", "minutes_range"),
                    )
        elif key == "xtol_contrast_index_times":
            for ri, row in enumerate(_rows(block, "rows")):
                path = (key, "rows", ri)
                units = (top, path)
                rid = f"{base}ci-{float(row['contrast_index']):.2f}/xtol"
                yield from _cells(
                    rid,
                    units,
                    row,
                    path,
                    block["full_strength_temperature_f"],
                    ("full_strength_minutes", "full_strength_printed"),
                )
                yield from _cells(
                    f"{rid}/1+1",
                    units,
                    row,
                    path,
                    block["diluted_1_1_temperature_f"],
                    ("diluted_1_1_minutes", "diluted_1_1_printed"),
                )
        elif key == "development_time_adjustment_factors":
            columns = [_ADJUST_COLUMNS[c] for c in block["columns"]]
            for bi, sub in enumerate(_rows(block, "blocks")):
                bkey = "+".join(sub["developers"])
                for ri, row in enumerate(_rows(sub, "rows")):
                    path = (key, "blocks", bi, "rows", ri)
                    units = (top, (key, "blocks", bi), path)
                    for ci, col in enumerate(columns):
                        fields = tuple((*path, k, ci) for k in ("factor", "printed") if k in row)
                        rid = f"{base}{bkey}/{col}{_temp(row['temperature_f'])}"
                        yield Item(rid, units, fields, null_fields=True)
        elif key == "push_processing":
            for ti, table in enumerate(block):
                yield from _prose(f"{base}{table['id']}", (top, (key, ti)), table)
                for ri, row in enumerate(_rows(table, "rows")):
                    path = (key, ti, "rows", ri)
                    yield Item(f"{base}{table['id']}/{_dev(row)}", (top, (key, ti), path), (path,))
        elif key == "rack_and_tank":
            for name in ("rows", "push_rows"):
                for ri, row in enumerate(_rows(block, name)):
                    path = (key, name, ri)
                    group = "normal"
                    if name == "push_rows":
                        ei = re.match(r"\d+", str(row["exposure_index"]))
                        if ei is None:
                            raise ValueError(f"{rel}: rack_and_tank push row without an EI")
                        group = f"push-ei{ei.group(0)}"
                    for dev in re.split(r"\s+or\s+", str(row["developer"])):
                        yield Item(f"{base}{group}/{_slug(dev)}", (top, path), (path,))
        elif key in ("replenishment",):
            for ri, row in enumerate(_rows(block, "rows")):
                path = (key, "rows", ri)
                yield Item(f"{base}{_slug(row['developer'])}", (top, path), (path,))
        elif key in ("push_exposure", "filter_factors", "daylight_exposure"):
            label = {
                "push_exposure": "push",
                "filter_factors": "filter",
                "daylight_exposure": "lighting",
            }[key]
            for ri, row in enumerate(_rows(block, "rows")):
                path = (key, "rows", ri)
                yield Item(f"{base}{_slug(row[label])}", (top, path), (path,))
        else:
            raise ValueError(f"{rel}: unknown block {key!r}; the split registry must classify it")


def _walk_developer(rel: str, data: Mapping[str, Any]) -> Iterator[Item]:
    for key, block in data.items():
        if key in _DEVELOPER_META:
            continue
        top: _Path = (key,)
        if key in ("facts", "agitation"):
            for i, entry in enumerate(block):
                yield Item(f"{rel}#{key}/{entry['id']}", (top, (key, i)), ((key, i),))
        elif key == "development_times":
            for ti, table in enumerate(block):
                yield from _prose(f"{rel}#{key}/{table['id']}", (top, (key, ti)), table)
                for ri, row in enumerate(_rows(table, "rows")):
                    material = row.get("material")
                    if material is None:
                        film = f"legacy-{_slug(row['film_as_printed'])}"
                    else:
                        film = PurePosixPath(str(material)).name
                    parts = [film, _slug(row.get("format", "any"))]
                    if row.get("ei") is not None:
                        parts.append(f"ei{_slug(row['ei'])}")
                    if row.get("ci") is not None:
                        parts.append(f"ci{float(row['ci']):.2f}")
                    path = (key, ti, "rows", ri)
                    rid = f"{rel}#{key}/{table['id']}/{'/'.join(parts)}"
                    yield Item(rid, (top, (key, ti), path), (path,))
        else:
            raise ValueError(f"{rel}: unknown block {key!r}; the split registry must classify it")


def walk_yaml(rel: str, data: Mapping[str, Any]) -> list[Item]:
    """Items of the loaded YAML file ``rel`` (relative to the reference root).

    Raises ``ValueError`` for an unknown block, a malformed row or a duplicate id.
    """
    parts = PurePosixPath(rel).parts
    if len(parts) == 3 and parts[0] == "kodak" and parts[2] == "tables.yaml":
        walker = _walk_tables
    elif len(parts) == 2 and parts[0] == "developers" and rel.endswith(".yaml"):
        walker = _walk_developer
    else:
        raise ValueError(f"{rel!r} is not a reference tables or developer-facts file")
    try:
        items = list(walker(rel, data))
    except (KeyError, TypeError, AttributeError) as exc:
        raise ValueError(f"{rel}: malformed reference table: {exc!r}") from exc
    seen: set[str] = set()
    for item in items:
        if item.id in seen:
            raise ValueError(f"{rel}: duplicate item id {item.id!r}")
        seen.add(item.id)
    return items


def _csv_items(rel: str, path: Path) -> list[str]:
    ids = [rel]
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            key, _, value = line[1:].partition(":")
            if key.strip() == "conditions" and "ci_cross_check=" in value:
                ids.append(rel + _CROSS_CHECK)
    return ids


def enumerate_items(root: Path | str | None = None) -> list[str]:
    """Sorted ids of every reference file, table row and sub-item under ``root``."""
    base = Path(root) if root is not None else default_reference_root()
    if not base.is_dir():
        raise FileNotFoundError(f"{base}: reference directory does not exist")
    ids: list[str] = []
    for path in sorted(base.rglob("*")):
        rel_path = path.relative_to(base)
        if not path.is_file() or any(p.startswith(".") for p in rel_path.parts):
            continue
        rel = rel_path.as_posix()
        if path.suffix == ".csv":
            ids.extend(_csv_items(rel, path))
        elif path.suffix == ".yaml":
            data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
            ids.extend(item.id for item in walk_yaml(rel, data))
        else:
            ids.append(rel)
    if len(ids) != len(set(ids)):
        raise ValueError(f"{base}: duplicate item ids")
    return sorted(ids)


# --------------------------------------------------------------------------------------
# Registry and hash (spec 7.7, 8.4)
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Registry:
    """The four disjoint classes, the evaluated held-out subset and the reasons (7.7)."""

    training: frozenset[str]
    within_checks: frozenset[str]
    held_out: frozenset[str]
    unused: frozenset[str]
    held_out_evaluated: frozenset[str]
    not_evaluated_reason: Mapping[str, str] = field(default_factory=dict)

    def canonical_json(self) -> str:
        """Sorted, whitespace-free JSON of the registry; independent of the root path."""
        doc = {
            "TRAINING": sorted(self.training),
            "WITHIN_CHECKS": sorted(self.within_checks),
            "HELD_OUT": sorted(self.held_out),
            "UNUSED": sorted(self.unused),
            "HELD_OUT_EVALUATED": sorted(self.held_out_evaluated),
            "NOT_EVALUATED_REASON": dict(sorted(self.not_evaluated_reason.items())),
        }
        return json.dumps(doc, sort_keys=True, separators=(",", ":"), ensure_ascii=True)

    def sha256(self) -> str:
        """Hex SHA-256 of :meth:`canonical_json` (``split_registry_sha256``, 8.4)."""
        return hashlib.sha256(self.canonical_json().encode("ascii")).hexdigest()


def build_registry(root: Path | str | None = None) -> Registry:
    """Classify every item found under ``root`` (default :func:`default_reference_root`)."""
    groups: dict[str, set[str]] = {s: set() for s in _SPLITS}
    evaluated: set[str] = set()
    reasons: dict[str, str] = {}
    for item_id in enumerate_items(root):
        c = classify(item_id)
        groups[c.split].add(item_id)
        if c.evaluated:
            evaluated.add(item_id)
        if c.split == "unused" or c.evaluated is False:
            reasons[item_id] = c.reason or ""
    return Registry(
        training=frozenset(groups["training"]),
        within_checks=frozenset(groups["within_checks"]),
        held_out=frozenset(groups["held_out"]),
        unused=frozenset(groups["unused"]),
        held_out_evaluated=frozenset(evaluated),
        not_evaluated_reason=MappingProxyType(dict(sorted(reasons.items()))),
    )


def registry_sha256(root: Path | str | None = None) -> str:
    """SHA-256 of the registry's canonical JSON (spec 8.4 ``split_registry_sha256``)."""
    return build_registry(root).sha256()


@functools.lru_cache(maxsize=1)
def _default_registry() -> Registry:
    return build_registry()


# Built from the default root on first access (module __getattr__ below).
TRAINING: frozenset[str]
WITHIN_CHECKS: frozenset[str]
HELD_OUT: frozenset[str]
UNUSED: frozenset[str]
HELD_OUT_EVALUATED: frozenset[str]
NOT_EVALUATED_REASON: Mapping[str, str]

_CONSTANTS = {
    "TRAINING": "training",
    "WITHIN_CHECKS": "within_checks",
    "HELD_OUT": "held_out",
    "UNUSED": "unused",
    "HELD_OUT_EVALUATED": "held_out_evaluated",
    "NOT_EVALUATED_REASON": "not_evaluated_reason",
}


def __getattr__(name: str) -> Any:
    if name in _CONSTANTS:
        return getattr(_default_registry(), _CONSTANTS[name])
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


# --------------------------------------------------------------------------------------
# Pair-metric table (spec 8.3)
# --------------------------------------------------------------------------------------

METRIC_NAMES = frozenset({
    "density_rms",          # characteristic-curve density RMS
    "density_max_error",    # characteristic-curve maximum error
    "left_end_density",     # |D_sim(x_left) - D_ref(x_left)|
    "fog_upper_bound",      # D_sim(unexposed) <= D_ref(x_left) + 0.02
    "log_h_m",              # log H_m* of a development series
    "toe_density",          # D at log H_m* + 0.3 and + 0.6 (held out only)
    "ci_at_time",           # CI on a CI(t) curve
    "table_ci",             # CI at the tabulated time of every table cell of the pair
    "mtf_tail",
    "spectral_shape",
    "granularity",
    "reciprocity_fitted",
})  # fmt: skip

_TX, _TM = "kodak/tri-x-400", "kodak/t-max-100"


def _curve_metrics(curves: list[str], *, toe: bool = False) -> list[str]:
    names = ["density_rms", "density_max_error", "left_end_density", "fog_upper_bound"]
    if toe:
        names.append("toe_density")
    return [f"{n}:{c}" for c in curves for n in names]


def _pair_metrics() -> tuple[dict[tuple[str, str], tuple[str, ...]], dict[tuple[str, str], str]]:
    x_d76 = [f"{_TX}/char-{f}-d76-{m}min.csv" for f in ("35mm", "120")
             for m in ("06", "08", "10", "12")]  # fmt: skip
    x_tmd = [f"{_TX}/char-{f}-tmax-dev-{m}min.csv" for f in ("35mm", "120")
             for m in ("06", "07", "09", "11")]  # fmt: skip
    m_d76 = [f"{_TM}/char-d76-{m}min.csv" for m in ("06", "07.5", "10")]
    m_tmd = [f"{_TM}/char-tmax-dev-{m}min.csv" for m in ("06", "07", "10", "12")]
    m_xtol = [f"{_TM}/ci-{t}-xtol{d}.csv" for t, d in (
        ("small-tank", ""), ("small-tank", "-1-1"), ("large-tank", ""), ("rotary-tube", ""),
        ("rotary-tube", "-1-1"), ("tray", ""), ("tray", "-1-1"))]  # fmt: skip

    def film_level(film: str, mtf: str) -> list[str]:
        return [
            f"mtf_tail:{film}/{mtf}",
            f"spectral_shape:{film}/spectral-sensitivity-d0.3.csv",
            f"spectral_shape:{film}/spectral-sensitivity-d1.0.csv",
            f"reciprocity_fitted:{film}/tables.yaml#reciprocity",
        ]

    def ci(film: str, names: list[str]) -> list[str]:
        return [f"ci_at_time:{film}/{n}" for n in names]

    table: dict[tuple[str, str], tuple[list[str], str]] = {
        (_TX, "d-76"): (
            _curve_metrics(x_d76)
            + [f"log_h_m:{_TX}/char-35mm-d76", f"log_h_m:{_TX}/char-120-d76"]
            + ci(_TX, ["ci-small-tank-d76.csv", "ci-small-tank-d76-1-1.csv"])
            + [f"table_ci:{_TX}/d-76"]
            + film_level(_TX, "mtf-d76-large-tank.csv"),
            "training",
        ),
        (_TX, "xtol"): (
            ci(_TX, ["ci-small-tank-xtol.csv", "ci-small-tank-xtol-1-1.csv"])
            + [f"table_ci:{_TX}/xtol"],
            "training",
        ),
        (_TX, "hc-110-pre2019"): (
            ci(_TX, ["ci-small-tank-hc110-b.csv"]) + [f"table_ci:{_TX}/hc-110-pre2019"],
            "training",
        ),
        (_TX, "t-max-dev"): (
            _curve_metrics(x_tmd)
            + [f"log_h_m:{_TX}/char-35mm-tmax-dev", f"log_h_m:{_TX}/char-120-tmax-dev"]
            + ci(_TX, ["ci-small-tank-tmax-dev.csv"])
            + [f"table_ci:{_TX}/t-max-dev"],
            "training",
        ),
        (_TM, "d-76"): (
            _curve_metrics(m_d76)
            + [f"log_h_m:{_TM}/char-d76"]
            + ci(
                _TM, ["ci-small-tank-d76.csv", "ci-small-tank-d76-1-1.csv", "ci-large-tank-d76.csv"]
            )
            + [f"table_ci:{_TM}/d-76"]
            + film_level(_TM, "mtf-d76-small-tank.csv")
            + [f"granularity:{_TM}/tables.yaml#facts/rms-granularity"],
            "training",
        ),
        (_TM, "xtol"): (
            [f"ci_at_time:{c}" for c in m_xtol] + [f"table_ci:{_TM}/xtol"],
            "held_out",
        ),
        (_TM, "hc-110-pre2019"): (
            ci(_TM, ["ci-small-tank-hc110-b.csv", "ci-large-tank-hc110-b.csv"])
            + [f"table_ci:{_TM}/hc-110-pre2019"],
            "held_out",
        ),
        (_TM, "t-max-dev"): (
            _curve_metrics(m_tmd, toe=True)
            + [f"log_h_m:{_TM}/char-tmax-dev"]
            + ci(_TM, ["ci-small-tank-tmax-dev.csv"])
            + [f"table_ci:{_TM}/t-max-dev"],
            "held_out",
        ),
    }
    for film in (_TX, _TM):
        for bath in ("rodinal", "hc-110-2019", "d-76-published"):
            table[(film, bath)] = ([], "predicted")
    metrics = {k: tuple(v) for k, (v, _) in table.items()}
    tiers = {k: t for k, (_, t) in table.items()}
    return metrics, tiers


_PM, _PT = _pair_metrics()
PAIR_METRICS: Mapping[tuple[str, str], tuple[str, ...]] = MappingProxyType(_PM)
"""(material, bath) -> applicable metric ids ``<metric>:<item or group>`` (spec 8.3).

Pruning of metrics the reference data cannot give (e.g. log H_m* when a development series
does not bracket Delta D = 0.80) is done before the freeze (8.3); not yet applied here.
"""
PAIR_TIERS: Mapping[tuple[str, str], Literal["training", "held_out", "predicted"]] = (
    MappingProxyType(_PT)
)
"""(material, bath) -> tier of the pair (spec 8.2)."""


# --------------------------------------------------------------------------------------
# Training view of a library (spec 7.7)
# --------------------------------------------------------------------------------------

_REMOVED = object()


def _get(data: Any, path: _Path) -> Any:
    for key in path:
        data = data[key]
    return data


def _set(data: Any, path: _Path, value: Any) -> None:
    _get(data, path[:-1])[path[-1]] = value


def _prune(node: Any) -> Any:
    if isinstance(node, dict):
        return {k: _prune(v) for k, v in node.items() if v is not _REMOVED}
    if isinstance(node, list):
        return [_prune(v) for v in node if v is not _REMOVED]
    return node


def _filter_held_out(rel: str, data: dict) -> dict:
    """Deep copy of ``data`` without its held-out items; a unit wholly held out goes too."""
    items = walk_yaml(rel, data)
    held = {item.id for item in items if classify(item.id).split == "held_out"}
    if not held:
        return copy.deepcopy(data)
    members: dict[_Path, list[str]] = {}
    for item in items:
        for unit in item.units:
            members.setdefault(unit, []).append(item.id)
    out = copy.deepcopy(data)
    for item in items:
        if item.id in held:
            for path in item.fields:
                if item.null_fields:
                    _set(out, path, None)
                else:
                    _set(out, path, _REMOVED)
    # Outermost units first; a removed ancestor makes its descendants unreachable.
    for unit in sorted(members, key=len):
        if all(i in held for i in members[unit]):
            try:
                _set(out, unit, _REMOVED)
            except TypeError:
                continue  # an ancestor was already replaced
    return _prune(out)


_FIELD_SPLIT = re.compile(r";\s*(?=[a-z][a-z0-9_]*=)")


class TrainingLibrary(Library):
    """Filtered view of ``base`` for calibration (spec 7.7).

    Reading a held-out curve raises :class:`LeakageError`; tables and developer facts come
    back without their held-out rows and sub-items (free text beside held-out cells included),
    and held-out header fields are removed.
    """

    def __init__(self, base: Library) -> None:
        if not isinstance(base, Library):
            raise TypeError(f"base must be a Library, got {type(base).__name__}")
        super().__init__(base.root)
        self._base = base

    def reference_curve(self, relpath: str | PurePosixPath) -> ReferenceCurve:
        item_id = PurePosixPath(str(relpath)).as_posix()
        try:
            held = classify(item_id).split == "held_out"
        except KeyError as exc:  # fail closed, e.g. a case variant of a held-out name
            raise LeakageError(f"{item_id!r} is not in the split registry") from exc
        if held:
            raise LeakageError(f"{item_id!r} is held out (7.6)")
        curve = self._base.reference_curve(relpath)
        if curve.conditions and "ci_cross_check=" in curve.conditions:
            kept = [
                part
                for part in _FIELD_SPLIT.split(curve.conditions)
                if not part.strip().startswith("ci_cross_check=")
            ]
            curve = curve.model_copy(update={"conditions": "; ".join(kept)})
        return curve

    def reference_tables(self, material_id: str) -> dict:
        data = self._base.reference_tables(material_id)
        return _filter_held_out(f"{material_id}/tables.yaml", data)

    def developer_facts(self, developer: str) -> dict:
        data = self._base.developer_facts(developer)
        name = DEVELOPER_FACTS_FILES[developer]
        return _filter_held_out(f"developers/{name}.yaml", data)
