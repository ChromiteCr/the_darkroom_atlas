"""Provenance tags and the :class:`Quantity` used for every physical parameter in data files.

Every number that describes a real material or chemical is a :class:`Quantity` that says
where it came from. Values taken from a manufacturer's datasheet or from the literature must
cite the document in ``ref``; fitted and assumed values may omit it.

Enumerations are :class:`~enum.StrEnum` s whose values are the strings written in data
files. Members can be reached by either an upper-case or a lower-case name
(``CurveOrigin.DATASHEET`` and ``CurveOrigin.datasheet`` are the same member); a hyphen in a
value becomes an underscore in the name (``FormulaConfidence.MSDS_ESTIMATED``).
"""

import math
from enum import EnumType, StrEnum
from numbers import Real
from typing import Any, Self

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from darkroom_atlas.units import convert

__all__ = [
    "CurveOrigin",
    "DataEnum",
    "FormulaConfidence",
    "FrozenModel",
    "Quantity",
    "Source",
    "VerificationStatus",
]


class _CaseInsensitiveNames(EnumType):
    """Enum metaclass that also resolves member names written in the other letter case.

    The lookup is a fallback only: it does not create aliases, so iteration, ``__members__``
    and the lists of allowed values in validation messages are unaffected.
    """

    def __getattr__(cls, name: str) -> Any:
        if not name.startswith("_"):
            members = cls.__dict__.get("_member_map_", {})
            for candidate in (name.upper(), name.lower()):
                if candidate in members:
                    return members[candidate]
        raise AttributeError(f"type object {cls.__name__!r} has no attribute {name!r}")


class DataEnum(StrEnum, metaclass=_CaseInsensitiveNames):
    """Base for the string enumerations used in data files."""


class Source(DataEnum):
    """Where a parameter value comes from."""

    DATASHEET = "datasheet"
    LITERATURE = "literature"
    FITTED = "fitted"
    ASSUMED = "assumed"


class FormulaConfidence(DataEnum):
    """How reliable a bath formula is."""

    PUBLISHED = "published"
    PATENT = "patent"
    MSDS_ESTIMATED = "msds-estimated"


class VerificationStatus(DataEnum):
    """Outcome of checking a material model against reference data."""

    VALIDATED = "validated"
    PREDICTED = "predicted"
    FAILED = "failed"


class CurveOrigin(DataEnum):
    """Origin of a reference curve: digitised from a datasheet, or measured."""

    DATASHEET = "datasheet"
    MEASURED = "measured"


REF_REQUIRED = frozenset({Source.DATASHEET, Source.LITERATURE})
"""Sources whose values must cite a reference."""


class FrozenModel(BaseModel):
    """Base for data-file models: immutable, and unknown fields are errors."""

    model_config = ConfigDict(frozen=True, extra="forbid")


def check_ref(source: Source, ref: str | None) -> None:
    """Raise ``ValueError`` if ``source`` requires a reference and ``ref`` is blank."""
    if source in REF_REQUIRED and (ref is None or not ref.strip()):
        raise ValueError(f"ref is required when source is {source.value!r}")


def _finite_float(value: Any) -> float:
    # YAML 1.1 reads exponent notation without a decimal point (``1e-3``) as a string, so
    # numeric strings are accepted. Booleans are rejected rather than read as 0 or 1.
    # ValueError (not TypeError) throughout, so that pydantic reports it as a validation error.
    if isinstance(value, str):
        try:
            result = float(value.strip())
        except ValueError:
            raise ValueError(f"expected a number, got {value!r}") from None
    elif isinstance(value, Real) and not isinstance(value, bool):
        result = float(value)
    else:
        raise ValueError(f"expected a number, got {value!r}")  # noqa: TRY004
    if not math.isfinite(result):
        raise ValueError(f"value must be finite, got {value!r}")
    return result


class Quantity(FrozenModel):
    """A physical value with its unit and provenance.

    ``value`` is a number or a list of numbers, all finite. ``unit`` is a plain ASCII unit
    string (see :mod:`darkroom_atlas.units`). ``uncertainty`` is one standard uncertainty in
    the same unit.
    """

    value: float | list[float]
    unit: str
    source: Source
    ref: str | None = None
    uncertainty: float | None = Field(default=None, ge=0.0, allow_inf_nan=False)
    note: str | None = None

    @field_validator("value", mode="before")
    @classmethod
    def _check_value(cls, value: Any) -> Any:
        if isinstance(value, (list, tuple)):
            if not value:
                raise ValueError("value list must not be empty")
            return [_finite_float(v) for v in value]
        return _finite_float(value)

    @field_validator("unit")
    @classmethod
    def _check_unit(cls, unit: str) -> str:
        if not unit or not unit.isascii() or not unit.isprintable() or " " in unit:
            raise ValueError(f"unit must be a non-empty ASCII string without spaces, got {unit!r}")
        return unit

    @model_validator(mode="after")
    def _check_ref(self) -> Self:
        check_ref(self.source, self.ref)
        return self

    def to(self, unit: str) -> float | np.ndarray:
        """Return the value converted to ``unit``: a float, or an array for list values."""
        if isinstance(self.value, list):
            return np.asarray(convert(np.asarray(self.value, dtype=np.float64), self.unit, unit))
        return float(convert(self.value, self.unit, unit))
