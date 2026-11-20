"""Pydantic models for the data files: chemical species, baths, materials and processes.

Every top-level file model carries ``schema_version: 1``. All models are immutable and reject
unknown fields, so a misspelt key is an error rather than a silently ignored value.

Identifiers:

* File ids (materials, baths, processes) are ``/``-separated segments such as
  ``kodak/tri-x-400``; each segment starts with a letter or digit and may contain
  ``. _ + -``. The library requires the id in a file to equal the id implied by its path.
* Species, layer and coupler ids are single segments. Layer ids name state fields
  (``"<layer_id>/<name>"``), so ``frame`` is reserved for fields that belong to no layer.

YAML note: PyYAML reads unquoted ``9:45`` as the base-60 integer 585 (which happens to equal
the duration in seconds) and unquoted ``1:1`` as 61. Quote durations and dilutions.
"""

import re
from pathlib import PurePosixPath
from typing import Annotated, Any, Literal, Self

import numpy as np
from pydantic import (
    AfterValidator,
    BeforeValidator,
    ConfigDict,
    Field,
    FiniteFloat,
    field_validator,
    model_validator,
)

from darkroom_atlas.data.provenance import (
    CurveOrigin,
    DataEnum,
    FormulaConfidence,
    FrozenModel,
    Quantity,
    Source,
    check_ref,
)
from darkroom_atlas.units import convert, parse_dilution, parse_duration_s

__all__ = [
    "BUILTIN_CURVES",
    "Agitation",
    "AntihalationLayer",
    "BaseLayer",
    "Bath",
    "BathComponent",
    "BathKind",
    "Coupler",
    "EmulsionLayer",
    "FilterLayer",
    "Grain",
    "InlineCurve",
    "InterlayerLayer",
    "Layer",
    "Material",
    "MaterialKind",
    "OvercoatLayer",
    "Process",
    "ProcessStep",
    "ReferenceCurve",
    "Species",
    "SpeciesFile",
    "SpeciesRole",
    "SpectralCurve",
    "StepKind",
    "Tank",
]


# --------------------------------------------------------------------------------------
# Enumerations
# --------------------------------------------------------------------------------------


class SpeciesRole(DataEnum):
    """What a chemical species does in a bath."""

    DEVELOPING_AGENT = "developing_agent"
    COLOR_DEVELOPING_AGENT = "color_developing_agent"
    PRESERVATIVE = "preservative"
    ALKALI = "alkali"
    BUFFER = "buffer"
    RESTRAINER = "restrainer"
    SOLVENT = "solvent"
    ACID = "acid"
    FIXER = "fixer"
    BLEACH = "bleach"
    COUPLER = "coupler"
    OTHER = "other"


class BathKind(DataEnum):
    """Kind of processing bath."""

    DEVELOPER = "developer"
    FIRST_DEVELOPER = "first_developer"
    COLOR_DEVELOPER = "color_developer"
    STOP = "stop"
    FIXER = "fixer"
    BLEACH = "bleach"
    BLEACH_FIX = "bleach_fix"
    REVERSAL = "reversal"
    WASH = "wash"
    STABILIZER = "stabilizer"
    TONER = "toner"


class MaterialKind(DataEnum):
    """Kind of light-sensitive material. ``ideal`` marks non-physical test materials."""

    BW_NEGATIVE = "bw_negative"
    COLOR_NEGATIVE = "color_negative"
    COLOR_REVERSAL = "color_reversal"
    BW_PAPER = "bw_paper"
    COLOR_PAPER = "color_paper"
    IDEAL = "ideal"


class StepKind(DataEnum):
    """Kind of process step; selects the step handler registered for a material model."""

    IDENTITY = "identity"
    PRESOAK = "presoak"
    DEVELOP = "develop"
    STOP = "stop"
    FIX = "fix"
    BLEACH = "bleach"
    BLEACH_FIX = "bleach_fix"
    REVERSAL = "reversal"
    WASH = "wash"
    STABILIZE = "stabilize"
    TONE = "tone"
    DRY = "dry"


STEPS_REQUIRING_BATH = frozenset(
    {StepKind.DEVELOP, StepKind.FIX, StepKind.BLEACH, StepKind.BLEACH_FIX, StepKind.TONE}
)
"""Step kinds that cannot run without a bath."""


# --------------------------------------------------------------------------------------
# Field types
# --------------------------------------------------------------------------------------

_SEGMENT = r"[A-Za-z0-9][A-Za-z0-9._+-]*"
DATA_ID_PATTERN = re.compile(rf"{_SEGMENT}(?:/{_SEGMENT})*")
NAME_PATTERN = re.compile(_SEGMENT)
RESERVED_LAYER_IDS = frozenset({"frame"})


def _data_id(value: str) -> str:
    if DATA_ID_PATTERN.fullmatch(value) is None:
        raise ValueError(
            f"invalid id {value!r}: use '/'-separated segments of letters, digits and "
            "'. _ + -', each starting with a letter or digit"
        )
    return value


def _name_id(value: str) -> str:
    if NAME_PATTERN.fullmatch(value) is None:
        raise ValueError(
            f"invalid id {value!r}: use letters, digits and '. _ + -', starting with a "
            "letter or digit (no '/')"
        )
    return value


def _layer_id(value: str) -> str:
    _name_id(value)
    if value in RESERVED_LAYER_IDS:
        raise ValueError(f"layer id {value!r} is reserved for fields that belong to no layer")
    return value


def _text(value: str) -> str:
    if not value.strip():
        raise ValueError("must not be empty")
    return value


DataId = Annotated[str, AfterValidator(_data_id)]
NameId = Annotated[str, AfterValidator(_name_id)]
LayerId = Annotated[str, AfterValidator(_layer_id)]
Text = Annotated[str, AfterValidator(_text)]
PositiveFloat = Annotated[float, Field(gt=0.0, allow_inf_nan=False)]
NonNegativeFloat = Annotated[float, Field(ge=0.0, allow_inf_nan=False)]


def _dimension(reference_unit: str, what: str) -> AfterValidator:
    """Validator requiring a :class:`Quantity` whose unit converts to ``reference_unit``."""

    def check(quantity: Quantity) -> Quantity:
        try:
            convert(1.0, quantity.unit, reference_unit)
        except ValueError:
            raise ValueError(
                f"unit {quantity.unit!r} is not a {what} unit (expected e.g. {reference_unit!r})"
            ) from None
        return quantity

    return AfterValidator(check)


Length = Annotated[Quantity, _dimension("um", "length")]
Temperature = Annotated[Quantity, _dimension("C", "temperature")]
MolarMass = Annotated[Quantity, _dimension("g/mol", "molar mass")]
ArealMass = Annotated[Quantity, _dimension("g/m2", "areal mass")]
Dimensionless = Annotated[Quantity, _dimension("1", "dimensionless")]


def _dilution(value: Any) -> Any:
    # ValueError, not TypeError, so that pydantic reports it as a validation error.
    if isinstance(value, int) and not isinstance(value, bool):
        raise ValueError(  # noqa: TRY004
            f"dilution must be a string such as 'stock' or '1+1', got the number {value}; "
            "quote it (YAML reads an unquoted a:b as a base-60 number)"
        )
    parse_dilution(value)
    return value


# A before-validator, so that the hint above is given instead of a bare type error.
Dilution = Annotated[str, BeforeValidator(_dilution)]


def _duplicates(values) -> list:
    seen: set = set()
    duplicates: list = []
    for value in values:
        if value in seen and value not in duplicates:
            duplicates.append(value)
        seen.add(value)
    return duplicates


# --------------------------------------------------------------------------------------
# Chemistry
# --------------------------------------------------------------------------------------

_CAS_PATTERN = re.compile(r"(\d{2,7})-(\d{2})-(\d)")


def _cas(value: str) -> str:
    match = _CAS_PATTERN.fullmatch(value)
    if match is None:
        raise ValueError(f"invalid CAS registry number {value!r}; expected the form NNNNNNN-NN-N")
    digits = match.group(1) + match.group(2)
    check = sum(i * int(d) for i, d in enumerate(reversed(digits), start=1)) % 10
    if check != int(match.group(3)):
        raise ValueError(
            f"CAS registry number {value!r} has a wrong check digit (computed {check})"
        )
    return value


class Species(FrozenModel):
    """A chemical species and its physical constants."""

    id: NameId
    name: Text
    formula: str | None = None
    cas: Annotated[str, AfterValidator(_cas)] | None = None
    molar_mass: MolarMass
    roles: list[SpeciesRole] = Field(min_length=1)
    pka: list[Dimensionless] = Field(default_factory=list)
    diffusion_coefficient: Quantity | None = None
    parameters: dict[str, Quantity] = Field(default_factory=dict)


class SpeciesFile(FrozenModel):
    """``chemistry/species.yaml``: every species the baths may refer to."""

    schema_version: Literal[1]
    species: list[Species]

    @model_validator(mode="after")
    def _unique_ids(self) -> Self:
        duplicates = _duplicates(s.id for s in self.species)
        if duplicates:
            raise ValueError(f"duplicate species id(s): {', '.join(map(repr, duplicates))}")
        return self


class BathComponent(FrozenModel):
    """One ingredient of a bath: a species id and its amount (e.g. g/L of working solution)."""

    species: Text
    amount: Quantity


class Bath(FrozenModel):
    """A processing bath (developer, bleach, fixer, ...) and its formula."""

    schema_version: Literal[1]
    id: DataId
    name: Text
    kind: BathKind
    formula_confidence: FormulaConfidence
    components: list[BathComponent]
    ph: Dimensionless | None = None
    standard_dilutions: list[Dilution] = Field(default_factory=lambda: ["stock"], min_length=1)
    refs: list[str] = Field(default_factory=list)
    note: str | None = None

    @model_validator(mode="after")
    def _unique_components(self) -> Self:
        duplicates = _duplicates(c.species for c in self.components)
        if duplicates:
            raise ValueError(f"duplicate species in components: {', '.join(map(repr, duplicates))}")
        return self


# --------------------------------------------------------------------------------------
# Spectral curves
# --------------------------------------------------------------------------------------

BUILTIN_CURVES: tuple[str, ...] = ("cie1931_x", "cie1931_y", "cie1931_z")
"""Names accepted by :attr:`SpectralCurve.builtin`."""


class InlineCurve(FrozenModel):
    """A curve given in the data file: values at strictly increasing wavelengths (nm)."""

    wavelengths_nm: list[FiniteFloat] = Field(min_length=2)
    values: list[FiniteFloat]

    @model_validator(mode="after")
    def _check(self) -> Self:
        if len(self.values) != len(self.wavelengths_nm):
            raise ValueError(
                f"wavelengths_nm and values must have the same length "
                f"({len(self.wavelengths_nm)} != {len(self.values)})"
            )
        if any(b <= a for a, b in zip(self.wavelengths_nm, self.wavelengths_nm[1:], strict=False)):
            raise ValueError("wavelengths_nm must be strictly increasing")
        if self.wavelengths_nm[0] <= 0.0:
            raise ValueError("wavelengths_nm must be positive")
        return self


def _curve_file(value: str) -> str:
    path = PurePosixPath(value)
    if not value.strip() or "\\" in value or path.is_absolute() or ".." in path.parts:
        raise ValueError(
            f"invalid curve file {value!r}: give a path relative to the data root, using '/' "
            "and without '..'"
        )
    return value


def _builtin(value: str) -> str:
    if value not in BUILTIN_CURVES:
        raise ValueError(f"unknown builtin curve {value!r}; available: {', '.join(BUILTIN_CURVES)}")
    return value


class SpectralCurve(FrozenModel):
    """A spectral quantity, given in exactly one of three ways.

    * ``inline``: wavelengths and values in the data file.
    * ``file``: a CSV file, relative to the data root, in the reference-curve format
      (optional ``# key: value`` header lines, then ``wavelength,value`` rows; wavelengths in
      nm unless the header gives another length unit as ``x_unit``).
    * ``builtin``: a named standard curve, one of :data:`BUILTIN_CURVES`.
    """

    quantity: Text
    source: Source
    ref: str | None = None
    inline: InlineCurve | None = None
    file: Annotated[str, AfterValidator(_curve_file)] | None = None
    builtin: Annotated[str, AfterValidator(_builtin)] | None = None

    @model_validator(mode="after")
    def _check(self) -> Self:
        given = [name for name in ("inline", "file", "builtin") if getattr(self, name) is not None]
        if len(given) != 1:
            got = ", ".join(given) if given else "none"
            raise ValueError(f"exactly one of inline, file or builtin must be given (got {got})")
        check_ref(self.source, self.ref)
        return self


# --------------------------------------------------------------------------------------
# Layers
# --------------------------------------------------------------------------------------


class Grain(FrozenModel):
    """Silver-halide grain population of an emulsion layer (log-normal size distribution)."""

    morphology: Literal["conventional", "tabular", "cubic"]
    mean_diameter: Length
    sigma_log: Dimensionless
    silver_coverage: ArealMass
    halide: str | None = None


class Coupler(FrozenModel):
    """A colour coupler incorporated in an emulsion layer."""

    id: NameId
    dye: Literal["cyan", "magenta", "yellow"]
    kind: Literal["image", "dir", "colored_mask"]
    dye_spectral_density: SpectralCurve | None = None
    parameters: dict[str, Quantity] = Field(default_factory=dict)


class OvercoatLayer(FrozenModel):
    """Protective top coat."""

    type: Literal["overcoat"] = "overcoat"
    id: LayerId
    thickness: Length | None = None


class EmulsionLayer(FrozenModel):
    """A light-sensitive layer."""

    type: Literal["emulsion"] = "emulsion"
    id: LayerId
    thickness: Length | None = None
    sensitized_to: str | None = None
    grain: Grain | None = None
    spectral_sensitivity: SpectralCurve | None = None
    couplers: list[Coupler] = Field(default_factory=list)
    parameters: dict[str, Quantity] = Field(default_factory=dict)

    @model_validator(mode="after")
    def _unique_couplers(self) -> Self:
        duplicates = _duplicates(c.id for c in self.couplers)
        if duplicates:
            raise ValueError(f"duplicate coupler id(s): {', '.join(map(repr, duplicates))}")
        return self


class FilterLayer(FrozenModel):
    """A filter layer, such as the yellow filter below the blue-sensitive layer."""

    type: Literal["filter"] = "filter"
    id: LayerId
    spectral_density: SpectralCurve | None = None
    removed_in_processing: bool = True


class InterlayerLayer(FrozenModel):
    """A gelatin interlayer, optionally containing a scavenger for oxidised developer."""

    type: Literal["interlayer"] = "interlayer"
    id: LayerId
    thickness: Length | None = None
    scavenger: bool = False


class AntihalationLayer(FrozenModel):
    """Antihalation protection."""

    type: Literal["antihalation"] = "antihalation"
    id: LayerId
    kind: Literal["undercoat", "dyed_base", "remjet", "backing"]
    spectral_density: SpectralCurve | None = None
    removed_in_processing: bool = True


class BaseLayer(FrozenModel):
    """Film base or paper support."""

    type: Literal["base"] = "base"
    id: LayerId
    material: Text
    thickness: Length | None = None
    refractive_index: Dimensionless | None = None
    spectral_density: SpectralCurve | None = None


Layer = Annotated[
    OvercoatLayer | EmulsionLayer | FilterLayer | InterlayerLayer | AntihalationLayer | BaseLayer,
    Field(discriminator="type"),
]
"""Any layer, selected by its ``type`` field."""


# --------------------------------------------------------------------------------------
# Materials
# --------------------------------------------------------------------------------------


class Material(FrozenModel):
    """A light-sensitive material: its layers from top to bottom and the model that runs it."""

    schema_version: Literal[1]
    id: DataId
    maker: Text
    name: Text
    kind: MaterialKind
    model: Text
    iso_speed: Dimensionless | None = None
    default_readout: str | None = None
    layers: list[Layer]
    refs: list[str] = Field(default_factory=list)
    note: str | None = None

    @model_validator(mode="after")
    def _check_layers(self) -> Self:
        duplicates = _duplicates(layer.id for layer in self.layers)
        if duplicates:
            raise ValueError(f"duplicate layer id(s): {', '.join(map(repr, duplicates))}")
        if not self.emulsion_layers:
            raise ValueError("a material needs at least one emulsion layer")
        if self.kind is not MaterialKind.IDEAL:
            n_base = sum(isinstance(layer, BaseLayer) for layer in self.layers)
            if n_base != 1:
                raise ValueError(
                    f"a {self.kind.value!r} material needs exactly one base layer, found {n_base}"
                )
        return self

    @property
    def emulsion_layers(self) -> tuple[EmulsionLayer, ...]:
        """The emulsion layers, from top to bottom."""
        return tuple(layer for layer in self.layers if isinstance(layer, EmulsionLayer))


# --------------------------------------------------------------------------------------
# Processes
# --------------------------------------------------------------------------------------


class Agitation(FrozenModel):
    """Agitation scheme: an initial period, then ``duration_s`` every ``interval_s``."""

    method: Literal["inversion", "rotary", "stand", "continuous", "tray"]
    initial_s: NonNegativeFloat = 0.0
    interval_s: PositiveFloat | None = None
    duration_s: NonNegativeFloat | None = None


class Tank(FrozenModel):
    """Processing tank: working-solution volume and the film area it serves."""

    volume_ml: PositiveFloat
    film_area_cm2: PositiveFloat | None = None


class ProcessStep(FrozenModel):
    """One step of a process.

    ``time`` is anything :func:`~darkroom_atlas.units.parse_duration_s` accepts (``"9:45"``,
    ``"30s"``, ``"8min"``, or a number of seconds) and must be positive. ``dilution`` is
    anything :func:`~darkroom_atlas.units.parse_dilution` accepts. Steps of the kinds
    ``develop fix bleach bleach_fix tone`` need a ``bath``; every step except ``identity``
    needs a ``time``.
    """

    kind: StepKind
    bath: Text | None = None
    dilution: Dilution = "stock"
    time: str | float | None = None
    temperature: Temperature | None = None
    agitation: Agitation | None = None
    tank: Tank | None = None
    note: str | None = None

    @field_validator("time", mode="before")
    @classmethod
    def _check_time(cls, value: Any) -> Any:
        if value is not None and not parse_duration_s(value) > 0.0:
            raise ValueError(f"time must be positive, got {value!r}")
        return value

    @model_validator(mode="after")
    def _check_requirements(self) -> Self:
        if self.kind in STEPS_REQUIRING_BATH and self.bath is None:
            raise ValueError(f"a {self.kind.value!r} step needs a bath")
        if self.kind is not StepKind.IDENTITY and self.time is None:
            raise ValueError(f"a {self.kind.value!r} step needs a time")
        return self

    @property
    def time_s(self) -> float | None:
        """Step duration in seconds, or ``None`` if no time is given."""
        return None if self.time is None else parse_duration_s(self.time)

    @property
    def dilution_fraction(self) -> float:
        """Fraction of stock solution in the working solution."""
        return parse_dilution(self.dilution)


class Process(FrozenModel):
    """A process: an ordered list of steps."""

    schema_version: Literal[1]
    id: DataId
    name: Text
    family: Text
    steps: list[ProcessStep] = Field(min_length=1)
    refs: list[str] = Field(default_factory=list)
    note: str | None = None


# --------------------------------------------------------------------------------------
# Reference curves
# --------------------------------------------------------------------------------------


def _curve_array(value: Any) -> np.ndarray:
    try:
        array = np.array(value, dtype=np.float64)
    except (TypeError, ValueError):
        raise ValueError("expected a sequence of numbers") from None
    if array.ndim != 1 or array.size < 2:
        raise ValueError(
            f"expected a one-dimensional sequence of at least 2 numbers, got {value!r}"
        )
    if not np.all(np.isfinite(array)):
        raise ValueError("values must be finite")
    array.flags.writeable = False
    return array


class ReferenceCurve(FrozenModel):
    """A curve digitised from a datasheet (or measured), used to check a material model.

    ``x`` and ``y`` are read-only ``float64`` arrays of equal length with ``x`` strictly
    increasing.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    material: Text
    quantity: Text
    x_unit: Text
    y_unit: Text
    origin: CurveOrigin
    ref: Text
    conditions: str | None = None
    x: np.ndarray
    y: np.ndarray

    @field_validator("x", "y", mode="before")
    @classmethod
    def _as_array(cls, value: Any) -> np.ndarray:
        return _curve_array(value)

    @model_validator(mode="after")
    def _check(self) -> Self:
        if self.x.shape != self.y.shape:
            raise ValueError(f"x and y must have the same length ({self.x.size} != {self.y.size})")
        if np.any(np.diff(self.x) <= 0.0):
            raise ValueError("x must be strictly increasing")
        return self
