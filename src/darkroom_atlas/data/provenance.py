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
from typing import Annotated, Any, Literal, Self

import numpy as np
from pydantic import (
    AfterValidator,
    BaseModel,
    ConfigDict,
    Field,
    field_validator,
    model_validator,
)

from darkroom_atlas.manifest import Sha256
from darkroom_atlas.units import convert

__all__ = [
    "CurveOrigin",
    "DataEnum",
    "FitRecord",
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


def _nonblank(value: str) -> str:
    if not value.strip():
        raise ValueError("must not be empty")
    return value


_Text = Annotated[str, AfterValidator(_nonblank)]


class FitRecord(FrozenModel):
    """How a ``fitted`` value was obtained (spec D1 5.1, 7.8).

    ``sigma_transformed`` is one standard deviation in the ``transform`` space; ``interval``
    is a 95 % interval in the quantity's own unit and is required for the profile method.
    """

    run_id: _Text
    model_version: _Text
    training_set_sha256: Sha256
    method: Literal["covariance", "profile", "prior"]
    transform: Literal["linear", "log10", "logit"]
    sigma_transformed: float = Field(ge=0.0, allow_inf_nan=False)
    interval: tuple[float, float] | None = None
    at_bound: bool = False
    identifiable: bool = True

    @field_validator("interval", mode="before")
    @classmethod
    def _check_interval(cls, value: Any) -> Any:
        if value is None:
            return None
        lo, hi = _pair(value, "interval")
        if lo > hi:
            raise ValueError(f"interval must have lo <= hi, got [{lo}, {hi}]")
        return (lo, hi)

    @model_validator(mode="after")
    def _check_method(self) -> Self:
        if self.method == "profile" and self.interval is None:
            raise ValueError("interval is required when method is 'profile'")
        return self


def _pair(value: Any, name: str) -> tuple[float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError(f"{name} must be a pair [lo, hi], got {value!r}")
    return _finite_float(value[0]), _finite_float(value[1])


class Quantity(FrozenModel):
    """A physical value with its unit and provenance.

    ``value`` is a number or a list of numbers, all finite. ``unit`` is a plain ASCII unit
    string (see :mod:`darkroom_atlas.units`). ``uncertainty`` is one standard uncertainty in
    the same unit.

    Calibration fields (spec D1 5.1): ``bounds`` is the prior box ``[lo, hi]`` in the same
    unit, and every value must lie inside it; ``prior`` names the prior shape (a
    ``lognormal`` prior needs positive bounds); ``fit`` is present exactly when ``source``
    is ``fitted``.
    """

    value: float | list[float]
    unit: str
    source: Source
    ref: str | None = None
    uncertainty: float | None = Field(default=None, ge=0.0, allow_inf_nan=False)
    note: str | None = None
    bounds: tuple[float, float] | None = None
    prior: Literal["uniform", "normal", "lognormal"] | None = None
    fit: FitRecord | None = None

    @field_validator("bounds", mode="before")
    @classmethod
    def _check_bounds(cls, value: Any) -> Any:
        if value is None:
            return None
        lo, hi = _pair(value, "bounds")
        if lo > hi:
            raise ValueError(f"bounds must have lo <= hi, got [{lo}, {hi}]")
        return (lo, hi)

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
        if self.source is Source.FITTED and self.fit is None:
            raise ValueError("a fitted value needs a fit record (fit)")
        if self.source is not Source.FITTED and self.fit is not None:
            raise ValueError(
                f"fit is only allowed when source is 'fitted', not {self.source.value!r}"
            )
        if self.bounds is not None:
            lo, hi = self.bounds
            values = self.value if isinstance(self.value, list) else [self.value]
            outside = [v for v in values if not lo <= v <= hi]
            if outside:
                raise ValueError(f"value {outside[0]} is outside bounds [{lo}, {hi}]")
            if self.prior == "lognormal" and not lo > 0.0:
                raise ValueError(f"a lognormal prior needs positive bounds, got [{lo}, {hi}]")
        return self

    def to(self, unit: str) -> float | np.ndarray:
        """Return the value converted to ``unit``: a float, or an array for list values."""
        if isinstance(self.value, list):
            return np.asarray(convert(np.asarray(self.value, dtype=np.float64), self.unit, unit))
        return float(convert(self.value, self.unit, unit))
