"""Test conditions of the virtual sensitometer (spec D1 section 6.1) and reference-CSV headers.

A :class:`Condition` names one film, developer, dilution, temperature, agitation, tank and the
development times at which the strip is read. :func:`parse_conditions_header` splits the
``# conditions:`` line of a digitised reference CSV into its ``key=value`` fields.
"""

import itertools
import math
import re
from dataclasses import dataclass
from typing import Literal

from darkroom_atlas.data.models import Agitation, Tank

__all__ = ["Condition", "fahrenheit_to_c", "parse_conditions_header", "parse_temperature_c"]

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
