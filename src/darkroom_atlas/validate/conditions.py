"""Test conditions of the virtual sensitometer (spec D1 section 6.1) and reference-CSV headers.

A :class:`Condition` names one film, developer, dilution, temperature, agitation, tank and the
development times at which the strip is read. :func:`parse_conditions_header` splits the
``# conditions:`` line of a digitised reference CSV into its ``key=value`` fields.

:func:`conditions_for` maps every reference curve and table cell of one split (spec 7.6) to a
:class:`Condition` built from the Kodak process file it uses (``data/processes/kodak/``,
naming rule in its README): agitation, tank and the exact temperature come from that file,
never from the header's free text. :func:`reference_conditions` returns the same together
with the process id and the table data a residual needs (CI aim, column temperature).
"""

import itertools
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import PurePosixPath
from typing import Any, Literal

from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import Agitation, Process, ProcessStep, Tank
from darkroom_atlas.units import parse_duration_s

__all__ = [
    "BATH_FOR_DEVELOPER",
    "Condition",
    "ConditionedItem",
    "NoProcessError",
    "bath_for_developer",
    "conditions_for",
    "dilution_for",
    "fahrenheit_to_c",
    "film_format_for",
    "parse_conditions_header",
    "parse_temperature_c",
    "process_id_for",
    "reference_conditions",
    "tank_kind",
]

_PREFIX = "# conditions:"
# A field starts after "; " with a lower-case identifier and "=". Values may contain ";",
# "," or "=" elsewhere, so the split looks ahead for the next "key=".
_FIELD_SPLIT = re.compile(r";\s*(?=[a-z][a-z0-9_]*=)")
_KEY = re.compile(r"[a-z][a-z0-9_]*")
_TEMPERATURE = re.compile(r"^\s*([-+]?\d+(?:\.\d+)?)\s*(?:°\s*)?([CF])\b", re.IGNORECASE)


@dataclass(frozen=True)
class Condition:
    """One sensitometric test condition (spec D1 6.1).

    ``times_s`` are the development times (s) at which the strip is read, positive and
    increasing. ``exposure_time_s`` defaults to 1/50 s (spec 6.1, DS 5.1, CAL 1.1).
    ``source`` identifies the reference file or table row the condition comes from (7.6).
    """

    material_id: str
    film_format: str
    developer: str
    dilution: str
    temperature_c: float
    agitation: Agitation
    times_s: tuple[float, ...]
    tank: Tank
    exposure_time_s: float = 0.02
    illuminant: Literal["iso_daylight", "A"] = "iso_daylight"
    source: str = ""

    def __post_init__(self) -> None:
        times = tuple(float(t) for t in self.times_s)
        object.__setattr__(self, "times_s", times)
        if not times:
            raise ValueError("Condition.times_s: at least one development time is required")
        if not all(math.isfinite(t) and t > 0.0 for t in times):
            raise ValueError(f"Condition.times_s: times must be finite and positive, got {times}")
        if any(b <= a for a, b in itertools.pairwise(times)):
            raise ValueError(f"Condition.times_s: times must be strictly increasing, got {times}")
        if not math.isfinite(self.temperature_c):
            raise ValueError(f"Condition.temperature_c must be finite, got {self.temperature_c}")
        if not (math.isfinite(self.exposure_time_s) and self.exposure_time_s > 0.0):
            raise ValueError(
                f"Condition.exposure_time_s must be positive, got {self.exposure_time_s}"
            )
        if self.illuminant not in ("iso_daylight", "A"):
            raise ValueError(
                f"Condition.illuminant must be 'iso_daylight' or 'A', got {self.illuminant!r}"
            )
        if not self.material_id:
            raise ValueError("Condition.material_id must not be empty")


def parse_conditions_header(text: str) -> dict[str, str]:
    """Return the ``key=value`` fields of the ``# conditions:`` line in ``text``.

    ``text`` is either that line alone or a whole CSV (header and data). Keys are lower-case
    identifiers; values are stripped but otherwise kept verbatim (parenthesised remarks
    included). Raises ``ValueError`` when there is no such line, a field has no ``=``, or a key
    repeats.
    """
    line = next((raw.strip() for raw in text.splitlines() if raw.strip().startswith(_PREFIX)), None)
    if line is None:
        raise ValueError("no '# conditions:' line found")
    body = line[len(_PREFIX) :].strip()
    if not body:
        raise ValueError("the '# conditions:' line is empty")
    fields: dict[str, str] = {}
    for part in _FIELD_SPLIT.split(body):
        key, sep, value = part.partition("=")
        key = key.strip()
        if not sep or not _KEY.fullmatch(key):
            raise ValueError(f"conditions field {part!r} is not of the form key=value")
        if key in fields:
            raise ValueError(f"conditions key {key!r} appears twice")
        fields[key] = value.strip().rstrip(";").strip()
    return fields


def fahrenheit_to_c(value_f: float) -> float:
    """Exact Fahrenheit to Celsius conversion (spec 6.1: Kodak's °F, not rounded °C labels)."""
    return (float(value_f) - 32.0) * 5.0 / 9.0


def parse_temperature_c(value: str) -> float | None:
    """Parse a ``temperature`` field such as ``20 C`` or ``68 F`` to °C; ``None`` if not stated.

    Fahrenheit values are converted exactly (spec 6.1). Anything else raises ``ValueError``.
    """
    text = value.strip()
    if text.lower().startswith("not stated"):
        return None
    match = _TEMPERATURE.match(text)
    if match is None:
        raise ValueError(f"cannot read a temperature from {value!r}")
    number, unit = float(match.group(1)), match.group(2).upper()
    return number if unit == "C" else fahrenheit_to_c(number)


# --------------------------------------------------------------------------------------
# From a reference header or table row to a Kodak process (spec 5.9, 6.1, 7.6)
# --------------------------------------------------------------------------------------

Split = Literal["train", "within", "held_out"]
_SPLIT_CLASS: dict[str, str] = {
    "train": "training",
    "within": "within_checks",
    "held_out": "held_out",
}

BATH_FOR_DEVELOPER: Mapping[str, str] = {
    # Spec 5.9: staged developer names to bath ids; t-max-rs is not used in D1 (7.6).
    "d-76": "d-76",
    "xtol": "xtol",
    "hc-110": "hc-110-pre2019",
    "t-max-dev": "t-max-dev",
}

# Fixed dilutions of the concentrates (data/processes/kodak/README.md, spec 5.5, 5.9):
# HC-110 dilution B is 1+31 of the concentrate (J-24 Dec 2017 p. 2), T-MAX Developer 1+4
# (J-86 Dec 2017 p. 1).
_CONCENTRATE_DILUTION = {"hc-110-pre2019": "1+31", "t-max-dev": "1+4"}

TANKS = ("small-tank", "large-tank", "rotary-tube", "tray")
_TANK_TEXT = (
    ("small tank", "small-tank"),
    ("large tank", "large-tank"),
    ("rotary tube", "rotary-tube"),
    ("rotary-tube", "rotary-tube"),
    ("tray", "tray"),
)
# The 68 F column of the Kodak tables is the 20 C reference temperature (spec 6.1).
_REFERENCE_F = 68
_PROCESS_ID = re.compile(
    r"^kodak/(?P<film>[a-z0-9-]+)/(?P<bath>[a-z0-9-]+)-(?P<dilution>stock|\d+\+\d+)-"
    r"(?P<temp>\d+(?:\.\d+)?)c-(?P<tank>small-tank|large-tank|rotary-tube|tray)-"
    r"(?P<time>\d+(?:\.\d+)?)min$"
)
_EXPOSURE_S = re.compile(r"(\d+(?:\.\d+)?)\s*s\b")
_TIME_TOL_S = 1e-6
_TEMP_TOL_C = 1e-6


class NoProcessError(LookupError):
    """No Kodak process file matches a reference item's condition (spec 7.6)."""


def bath_for_developer(name: str) -> str:
    """Bath id of a staged developer name (spec 5.9); ``ValueError`` if D1 has none."""
    key = str(name).strip().lower()
    if key not in BATH_FOR_DEVELOPER:
        raise ValueError(
            f"developer {name!r} has no D1 bath (spec 5.9 maps {sorted(BATH_FOR_DEVELOPER)})"
        )
    return BATH_FOR_DEVELOPER[key]


def dilution_for(bath: str, text: str | None) -> str:
    """Explicit dilution of a header or table ``dilution`` text for ``bath``.

    A printed ``1+n`` is kept. Otherwise HC-110 dilution B is 1+31 and T-MAX Developer 1+4
    (processes README); D-76 and XTOL read as stock (also when the dilution is not stated).
    """
    value = (text or "").strip().lower()
    m = re.match(r"^1\+(\d+)\b", value)
    if m:
        return f"1+{m.group(1)}"
    if bath == "hc-110-pre2019":
        if not value.startswith("b"):
            raise ValueError(f"HC-110 dilution {text!r}: only dilution B is used (spec 7.6)")
        return _CONCENTRATE_DILUTION[bath]
    if bath in _CONCENTRATE_DILUTION:
        return _CONCENTRATE_DILUTION[bath]
    if value.startswith(("stock", "not stated")) or not value:
        return "stock"
    raise ValueError(f"cannot read a dilution of {bath} from {text!r}")


def tank_kind(text: str) -> str:
    """``small-tank``, ``large-tank``, ``rotary-tube`` or ``tray`` from an agitation or tank text."""
    value = str(text).strip().lower()
    for prefix, kind in _TANK_TEXT:
        if value.startswith(prefix):
            return kind
    raise ValueError(f"cannot read a tank type from {text!r}")


def film_format_for(text: str | None) -> str:
    """Film format of a header ``format`` field.

    ``35 mm`` is 135 and ``120-size`` is 120. Anything else (not stated, sheets, or
    "35 mm and 120-size") is simulated as 135 (spec 6.7 step 7, 7.6 sheet rows).
    """
    value = (text or "").strip().lower()
    if value.startswith("120"):
        return "120"
    return "135"


@dataclass(frozen=True)
class _Develop:
    process_id: str
    bath: str
    dilution: str
    temperature_c: float
    tank: str
    time_s: float
    step: ProcessStep
    refs: frozenset[str]


@dataclass(frozen=True)
class ConditionedItem:
    """One reference item of a split with its :class:`Condition` (spec 6.1, 7.6).

    ``kind`` says which observation the item gives (``density``, ``ci``, ``table_ci``,
    ``mtf``, ``spectral``, ``granularity``, ``reciprocity``). ``process_id`` is the Kodak
    process the condition uses, ``None`` only when built with ``require_process=False``.
    Table cells carry their ``ci_aim``, Fahrenheit column and dilution/tank class key.
    """

    item_id: str
    kind: str
    condition: Condition
    process_id: str | None
    ci_aim: float | None = None
    temperature_f: int | None = None
    tank: str | None = None


def _develop_steps(library: Library, material_id: str) -> list[_Develop]:
    film = material_id.split("/", 1)[-1]
    found = []
    for pid in library.ids("process"):
        m = _PROCESS_ID.match(pid)
        if m is None or not pid.startswith("kodak/"):
            continue
        if m["film"] != film:
            continue
        process: Process = library.process(pid)
        steps = [s for s in process.steps if s.kind.value == "develop"]
        if len(steps) != 1 or steps[0].temperature is None:
            raise ValueError(f"process {pid!r}: expected one develop step with a temperature")
        s = steps[0]
        found.append(
            _Develop(
                process_id=pid,
                bath=str(s.bath),
                dilution=str(s.dilution),
                temperature_c=float(s.temperature.value),  # type: ignore[arg-type]
                tank=m["tank"],
                time_s=float(s.time_s or 0.0),
                step=s,
                refs=frozenset(process.refs),
            )
        )
    return found


def _templates(library: Library) -> dict[str, tuple[Agitation, Tank]]:
    """Agitation and tank of each tank type, as every Kodak process file gives them."""
    out: dict[str, tuple[Agitation, Tank]] = {}
    for pid in library.ids("process"):
        m = _PROCESS_ID.match(pid)
        if m is None:
            continue
        step = next(s for s in library.process(pid).steps if s.kind.value == "develop")
        if step.agitation is None or step.tank is None:
            continue
        pair = (step.agitation, step.tank)
        if m["tank"] in out and out[m["tank"]] != pair:
            raise ValueError(f"process {pid!r}: agitation or tank differs from other {m['tank']}")
        out[m["tank"]] = pair
    return out


@dataclass(frozen=True)
class _Want:
    item_id: str
    kind: str
    bath: str
    dilution: str
    temperature_c: float
    tank: str
    times_s: tuple[float, ...]
    process_time_s: float | None
    film_format: str = "135"
    exposure_time_s: float = 0.02
    illuminant: Literal["iso_daylight", "A"] = "iso_daylight"
    ci_aim: float | None = None
    temperature_f: int | None = None


class _Builder:
    """Matches wanted conditions to process files of one film."""

    def __init__(self, library: Library, material_id: str) -> None:
        self.library = library
        self.material_id = material_id
        self.steps = _develop_steps(library, material_id)
        self._templates: dict[str, tuple[Agitation, Tank]] | None = None

    def match(self, want: _Want) -> _Develop | None:
        hits = [
            d
            for d in self.steps
            if d.bath == want.bath
            and d.dilution == want.dilution
            and d.tank == want.tank
            and abs(d.temperature_c - want.temperature_c) <= _TEMP_TOL_C
            and (want.process_time_s is None or abs(d.time_s - want.process_time_s) <= _TIME_TOL_S)
        ]
        listed = [d for d in hits if want.item_id in d.refs]
        if len(listed) == 1:
            return listed[0]
        if len(hits) == 1:
            return hits[0]
        if len(hits) > 1:
            raise NoProcessError(
                f"{want.item_id}: several processes match ({', '.join(d.process_id for d in hits)})"
            )
        return None

    def build(self, want: _Want, *, require_process: bool) -> ConditionedItem:
        dev = self.match(want)
        if dev is None:
            if require_process:
                # The time is left out of the message: it may be a held-out table value.
                raise NoProcessError(
                    f"{want.item_id}: no process in data/processes/kodak/ for "
                    f"{self.material_id} in {want.bath} {want.dilution}, "
                    f"{want.temperature_c:.4g} C, {want.tank} at the item's time"
                )
            if self._templates is None:
                self._templates = _templates(self.library)
            if want.tank not in self._templates:
                raise NoProcessError(f"{want.item_id}: no process uses a {want.tank}")
            agitation, tank = self._templates[want.tank]
            temperature_c, pid = want.temperature_c, None
        else:
            assert dev.step.agitation is not None and dev.step.tank is not None
            agitation, tank = dev.step.agitation, dev.step.tank
            temperature_c, pid = dev.temperature_c, dev.process_id
        condition = Condition(
            material_id=self.material_id,
            film_format=want.film_format,
            developer=want.bath,
            dilution=want.dilution,
            temperature_c=temperature_c,
            agitation=agitation,
            times_s=want.times_s,
            tank=tank,
            exposure_time_s=want.exposure_time_s,
            illuminant=want.illuminant,
            source=want.item_id,
        )
        return ConditionedItem(
            item_id=want.item_id,
            kind=want.kind,
            condition=condition,
            process_id=pid,
            ci_aim=want.ci_aim,
            temperature_f=want.temperature_f,
            tank=want.tank,
        )


# ---- reading the reference data -------------------------------------------------------


def _csv_kind(name: str) -> str | None:
    for prefix, kind in (
        ("char-", "density"),
        ("ci-", "ci"),
        ("mtf-", "mtf"),
        ("spectral-sensitivity-", "spectral"),
    ):
        if name.startswith(prefix):
            return kind
    return None


def _table_minutes(tables: Mapping[str, Any], dev: str, dilution: str, tank: str) -> float | None:
    """Minutes of the 68 F (20 C) cell of the film's normal table for one development."""
    for table in tables.get("development_times") or []:
        tid = str(table.get("id", ""))
        if not tid.startswith(f"normal-{tank}"):
            continue
        temps = list(table.get("temperature_f") or [])
        if _REFERENCE_F not in temps:
            continue
        j = temps.index(_REFERENCE_F)
        for row in table.get("rows") or []:
            if row.get("developer") not in BATH_FOR_DEVELOPER:
                continue
            bath = bath_for_developer(row["developer"])
            if bath != dev or dilution_for(bath, row.get("dilution")) != dilution:
                continue
            minutes = (row.get("minutes") or [None] * len(temps))[j]
            if minutes is not None:
                return float(minutes)
    return None


class _Reader:
    def __init__(self, library: Library, material_id: str) -> None:
        self.library = library
        self.material_id = material_id
        self._tables: dict | None = None

    @property
    def tables(self) -> dict:
        if self._tables is None:
            self._tables = self.library.reference_tables(self.material_id)
        return self._tables

    def aim_time_s(self, item_id: str, bath: str, dilution: str, tank: str) -> float:
        """Kodak's 20 C time of one development (the 68 F table cell), in s."""
        minutes = _table_minutes(self.tables, bath, dilution, tank)
        if minutes is None:
            raise NoProcessError(
                f"{item_id}: no 68 F time for {bath} {dilution} {tank} in "
                f"{self.material_id}/tables.yaml"
            )
        return minutes * 60.0


def _ci_aim_default(tables: Mapping[str, Any]) -> float:
    """CI aim of the F-4016/F-4017 time tables (spec 7.2): the film's processing-aim fact,
    else 0.56 (``assumed`` for F-4016, which prints no CI aim, spec 7.2)."""
    for fact in tables.get("facts") or []:
        if str(fact.get("id", "")).startswith("processing-aim") and isinstance(
            fact.get("value"), int | float
        ):
            return float(fact["value"])
    return 0.56  # spec 7.2: CI_aim of the F-4016 tables, assumed


def _csv_wants(reader: _Reader, rel: str, kind: str) -> _Want:
    curve = reader.library.reference_curve(rel)
    fields = parse_conditions_header(f"# conditions: {curve.conditions or ''}")
    film = reader.material_id
    name = PurePosixPath(rel).name
    developer = fields.get("developer", "")
    if developer.startswith("not stated"):
        developer = "d-76"  # spec 6.9 step 2: Tri-X spectral, D-76 assumed
    bath = bath_for_developer(developer)
    dilution = dilution_for(bath, fields.get("dilution"))
    temperature = parse_temperature_c(fields.get("temperature", "not stated"))
    if temperature is None:
        temperature = 20.0  # spec 6.7 step 7, 6.9 step 2: 20 C assumed where not stated
    if "agitation" in fields:
        tank = tank_kind(fields["agitation"])
    else:
        # Spec 6.9 step 2: Tri-X spectral curves D-76 large tank (assumed); the T-Max
        # spectral figure names no tank, the small tank of its other D-76 figures is used.
        tank = "large-tank" if film.endswith("tri-x-400") else "small-tank"
    film_format = film_format_for(fields.get("format"))
    exposure = 0.02
    illuminant: Literal["iso_daylight", "A"] = "iso_daylight"
    if kind == "spectral":
        m = _EXPOSURE_S.search(fields.get("exposure", ""))
        if m is None:
            raise ValueError(f"{rel}: no exposure time in the header (spec 6.9 step 1)")
        exposure = float(m.group(1))
    if kind == "mtf" and fields.get("exposure", "").startswith("tungsten"):
        illuminant = "A"  # spec 6.1: CIE A for "Tungsten"
    if kind == "density":
        t = parse_duration_s(fields["time"])
        return _Want(rel, kind, bath, dilution, temperature, tank, (t,), t, film_format)
    aim = reader.aim_time_s(rel, bath, dilution, tank)
    if kind == "ci":
        x = curve.x * 60.0
        times = tuple(float(v) for v in _linspace(float(x[0]), float(x[-1]), 5))
    else:
        times = (aim,)
    del name
    return _Want(
        rel, kind, bath, dilution, temperature, tank, times, aim, film_format, exposure, illuminant
    )


def _linspace(a: float, b: float, n: int) -> list[float]:
    return [a + (b - a) * i / (n - 1) for i in range(n)]


def _field_index(fields: tuple, key: str) -> tuple[tuple, int] | None:
    for path in fields:
        if len(path) >= 2 and path[-2] == key:
            return tuple(path[:-2]), int(path[-1])
    return None


def _get(data: Any, path: tuple) -> Any:
    for key in path:
        data = data[key]
    return data


def _table_wants(reader: _Reader, wanted: set[str]) -> list[_Want]:
    from darkroom_atlas.validate.splits import walk_yaml

    film = reader.material_id
    rel = f"{film}/tables.yaml"
    data = reader.tables
    aim_default = _ci_aim_default(data)
    out: list[_Want] = []
    for item in walk_yaml(rel, data):
        if item.id not in wanted:
            continue
        block = item.id.split("#", 1)[1].split("/", 1)[0]
        if block == "development_times":
            hit = _field_index(item.fields, "minutes")
            if hit is None:
                continue
            row_path, j = hit
            row, table = _get(data, row_path), _get(data, row_path[:2])
            minutes = row["minutes"][j]
            if minutes is None:
                continue
            bath = bath_for_developer(row["developer"])
            dilution = dilution_for(bath, row.get("dilution"))
            temp_f = int(table["temperature_f"][j])
            t = float(minutes) * 60.0
            out.append(
                _Want(
                    item.id,
                    "table_ci",
                    bath,
                    dilution,
                    fahrenheit_to_c(temp_f),
                    tank_kind(table["tank"]),
                    (t,),
                    t,
                    ci_aim=aim_default,
                    temperature_f=temp_f,
                )
            )
        elif block == "xtol_contrast_index_times":
            for key, temps_key, dilution in (
                ("full_strength_minutes", "full_strength_temperature_f", "stock"),
                ("diluted_1_1_minutes", "diluted_1_1_temperature_f", "1+1"),
            ):
                hit = _field_index(item.fields, key)
                if hit is None:
                    continue
                row_path, j = hit
                row, blk = _get(data, row_path), data[block]
                minutes = row[key][j]
                if minutes is None:
                    continue
                temp_f = int(blk[temps_key][j])
                t = float(minutes) * 60.0
                out.append(
                    _Want(
                        item.id,
                        "table_ci",
                        "xtol",
                        dilution,
                        fahrenheit_to_c(temp_f),
                        tank_kind(blk["tank"]),
                        (t,),
                        t,
                        ci_aim=float(row["contrast_index"]),  # J-109 table 1: the row's CI
                        temperature_f=temp_f,
                    )
                )
        elif block == "development_time_adjustment_factors":
            if "/normal@" not in item.id:
                continue
            hit = _field_index(item.fields, "factor")
            if hit is None:
                continue
            row_path, j = hit
            row, sub = _get(data, row_path), _get(data, row_path[:3])
            factor = row["factor"][j]
            if factor is None:
                continue
            bath = bath_for_developer(sub["developers"][0])
            # The HC-110 block is read at dilution B, the dilution of the F-4016 tables.
            dilution = dilution_for(bath, "B" if bath == "hc-110-pre2019" else None)
            # Spec 7.6: t = f x t_std, t_std the small-tank time at the block's T_std.
            t_std = _std_minutes(data, bath, dilution, float(sub["standard_temperature_c"]))
            if t_std is None:
                raise NoProcessError(f"{item.id}: no small-tank time at T_std for {bath}")
            temp_f = int(row["temperature_f"])
            t = round(float(factor) * t_std * 60.0, 9)
            out.append(
                _Want(
                    item.id,
                    "table_ci",
                    bath,
                    dilution,
                    fahrenheit_to_c(temp_f),
                    "small-tank",
                    (t,),
                    t,
                    ci_aim=aim_default,
                    temperature_f=temp_f,
                )
            )
        elif block == "reciprocity" and not item.id.endswith("/development"):
            row = _get(data, item.units[-1])
            t_ind = float(row["indicated_time_s"])
            # Spec 6.8: reference development = D-76 stock, small tank, the Kodak 20 C time
            # (developed to CI 0.56); the film files name no reference process.
            aim = reader.aim_time_s(item.id, "d-76", "stock", "small-tank")
            out.append(
                _Want(
                    item.id,
                    "reciprocity",
                    "d-76",
                    "stock",
                    20.0,
                    "small-tank",
                    (aim,),
                    aim,
                    exposure_time_s=t_ind,
                )
            )
        elif block == "facts" and item.id.endswith("#facts/rms-granularity"):
            # Spec 6.6 step 5: T-Max D-76, Tri-X HC-110 B; 20 C, small tank, CI 0.56 time.
            bath = "hc-110-pre2019" if film.endswith("tri-x-400") else "d-76"
            dilution = dilution_for(bath, "B" if bath == "hc-110-pre2019" else "stock")
            aim = reader.aim_time_s(item.id, bath, dilution, "small-tank")
            out.append(
                _Want(item.id, "granularity", bath, dilution, 20.0, "small-tank", (aim,), aim)
            )
    return out


def _std_minutes(data: Mapping[str, Any], bath: str, dilution: str, t_std_c: float) -> float | None:
    for table in data.get("development_times") or []:
        if not str(table.get("id", "")).startswith("normal-small-tank"):
            continue
        temps_c = [float(v) for v in table.get("temperature_c") or []]
        if t_std_c not in temps_c:
            continue
        j = temps_c.index(t_std_c)
        for row in table.get("rows") or []:
            if row.get("developer") not in BATH_FOR_DEVELOPER:
                continue
            b = bath_for_developer(row["developer"])
            if b == bath and dilution_for(b, row.get("dilution")) == dilution:
                minutes = (row.get("minutes") or [None] * len(temps_c))[j]
                if minutes is not None:
                    return float(minutes)
    return None


def reference_conditions(
    material_id: str,
    library: Library,
    *,
    split: Split,
    require_process: bool = True,
) -> list[ConditionedItem]:
    """Every curve and table cell of ``material_id`` in ``split`` with its condition.

    ``split`` is ``train`` (registry ``TRAINING``), ``within`` (``WITHIN_CHECKS``) or
    ``held_out`` (``HELD_OUT_EVALUATED``). Train and within read through a
    :class:`~darkroom_atlas.validate.splits.TrainingLibrary`; held out refuses one. Cells whose
    time is null are skipped. Raises :class:`NoProcessError` naming the item when no process
    file matches (unless ``require_process`` is false, then agitation and tank come from the
    process files of the same tank type and ``process_id`` is ``None``).
    """
    from darkroom_atlas.validate.splits import (
        LeakageError,
        TrainingLibrary,
        classify,
        enumerate_items,
    )

    if split not in _SPLIT_CLASS:
        raise ValueError(f"split must be one of {sorted(_SPLIT_CLASS)}, got {split!r}")
    if not isinstance(library, Library):
        raise TypeError(f"library must be a Library, got {type(library).__name__}")
    if split == "held_out":
        if isinstance(library, TrainingLibrary):
            raise LeakageError("held-out conditions cannot be read through a TrainingLibrary")
        view = library
    else:
        view = library if isinstance(library, TrainingLibrary) else TrainingLibrary(library)
    target = _SPLIT_CLASS[split]
    prefix = f"{material_id}/"
    ids = []
    for item_id in enumerate_items(view.root / "reference"):
        if not item_id.startswith(prefix):
            continue
        c = classify(item_id)
        if c.split != target or (split == "held_out" and not c.evaluated):
            continue
        ids.append(item_id)
    reader = _Reader(view, material_id)
    builder = _Builder(view, material_id)
    wants: list[_Want] = []
    for item_id in ids:
        rel = PurePosixPath(item_id)
        if "#" in item_id or rel.suffix != ".csv":
            continue
        kind = _csv_kind(rel.name)
        if kind is not None:
            wants.append(_csv_wants(reader, item_id, kind))
    if any("/tables.yaml#" in i for i in ids):
        wants.extend(_table_wants(reader, {i for i in ids if "/tables.yaml#" in i}))
    return [builder.build(w, require_process=require_process) for w in wants]


def conditions_for(material_id: str, library: Library, *, split: Split) -> list[Condition]:
    """One :class:`Condition` per reference curve or table cell of ``split`` (spec 10 D1e).

    ``Condition.source`` is the item id (spec 7.6). Raises :class:`NoProcessError` naming
    the item when no Kodak process file matches its condition.
    """
    return [c.condition for c in reference_conditions(material_id, library, split=split)]


def process_id_for(condition: Condition, library: Library) -> str:
    """Id of the Kodak process file the :class:`Condition` of item ``condition.source`` uses."""
    from darkroom_atlas.validate.splits import classify

    if not condition.source:
        raise NoProcessError("the condition names no reference item (Condition.source)")
    c = classify(condition.source)
    split = {v: k for k, v in _SPLIT_CLASS.items()}.get(c.split)
    if split is None:
        raise NoProcessError(f"{condition.source}: {c.split} items have no condition")
    for item in reference_conditions(condition.material_id, library, split=split):  # type: ignore[arg-type]
        if item.item_id == condition.source and item.condition == condition:
            assert item.process_id is not None
            return item.process_id
    raise NoProcessError(f"{condition.source}: no process for this condition")
