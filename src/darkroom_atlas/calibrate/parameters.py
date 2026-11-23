"""Calibration parameter table, transforms and the in-memory library overlay (spec D1 7.1,
7.3, 7.5, 5.1).

:func:`parameter_specs` lists, for one fitting stage of spec D1 7.3, exactly the parameters
that the 7.1 table assigns to it. The assignment of parameters to blocks and stages and the
transform of each one are declared here, from the spec. Each data file supplies the value, unit,
``bounds`` and ``prior`` of its own fields (spec D1 5.1). A listed field without bounds or
without a prior is a :class:`~darkroom_atlas.data.DataError` that names the file and the field.
There are two exceptions. A normal prior with an ``uncertainty`` and no bounds gets infinite
bounds, because spec D1 2.4.11 gives the metol and hydroquinone activation energies no box. A bath
``nuisance`` amount without a prior is uniform inside its bounds, because spec D1 5.5 defines a
nuisance as fitted within its bounds.

Names are ``"<file id>:<field path>"``. The file id is the file's path under the data root
without ``.yaml`` (``films/kodak/tri-x-400``, ``chemistry/kinetics``, ``developers/d-76``). In a
material the field path starts with the layer id (``emulsion.latent.tau0_median``). In a bath a
component is selected by its species id (``components.boric-oxide.amount``). List items are
indexed (``emulsion.optics.dye_bands[0].amplitude``). Parameters that exist only in the
calibration (spec D1 7.1) are ``"calibration:<name>"``. The Carlton & Crabtree offset is
``calibration:delta_cc``. The two per-film parameters carry the material id:
``calibration:<material id>/s_abs`` and, for tabular grains,
``calibration:<material id>/aspect_ratio``.

Transforms (spec D1 7.1): ``log10`` for the scale parameters with a log prior (``lognormal`` in
the data files) and for the aspect ratio, ``logit`` for the fractions in (0, 1) (latent
efficiency and PSF weight), ``linear`` for the rest. ``s_abs`` is already a log10 factor with box
[-1, 1], so it is linear. A box prior (``uniform`` or ``lognormal`` without a width) has
``prior_tau = inf`` and adds no prior residual. A normal prior has ``prior_mu`` equal to the
transformed file value and ``prior_tau`` equal to its ``uncertainty``.

:class:`LibraryOverlay` is a :class:`~darkroom_atlas.data.Library` whose ``material``, ``bath``
and ``kinetics`` return the base library's models with the overridden values and nothing written
to disk. It also applies the derived parametrisations of spec D1 7.1 and 7.5. With
``calibration:<id>/aspect_ratio`` it sets ``grain.mean_diameter = AR * grain.thickness``. With
``calibration:<id>/s_abs`` it multiplies ``alpha0`` and every dye-band ``amplitude`` by
``10**s_abs``. The derived values are applied once, on top of the data overrides of the whole
overlay chain, and never feed back: :func:`parameter_value` reads each parameter's own value
(the unscaled ``alpha0``), and :func:`parameter_specs` reads bounds and priors from the root
library's files, so ``prior_mu`` stays at the data-file value. Derived values are not checked against the file bounds: the ``[0.6, 1.5] um`` box
of the T-Max ``mean_diameter`` is a one-sided penalty, not a hard bound.
"""

import math
import re
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from numbers import Real
from pathlib import Path, PurePosixPath
from types import MappingProxyType
from typing import Any, Literal, get_args

import numpy as np
from pydantic import BaseModel

from darkroom_atlas.data.library import DataError, DataNotFound, Kind, Library
from darkroom_atlas.data.models import (
    Bath,
    CalibrationRecord,
    EmulsionLayer,
    KineticsFile,
    Material,
    Process,
    ReferenceCurve,
    Species,
    SpectralCurve,
)
from darkroom_atlas.data.provenance import Quantity
from darkroom_atlas.spectral import SpectralGrid

__all__ = [
    "CALIBRATION_FILMS",
    "CALIBRATION_PREFIX",
    "KINETICS_FILE_ID",
    "STAGES",
    "LibraryOverlay",
    "ParameterError",
    "ParameterSpec",
    "current_values",
    "from_theta",
    "parameter_specs",
    "parameter_value",
    "theta_bounds",
    "to_theta",
]

Transform = Literal["linear", "log10", "logit"]
PriorShape = Literal["uniform", "normal", "lognormal"]
Block = Literal["A", "B", "F", "C"]

STAGES: tuple[str, ...] = ("A", "B", "C", "D", "E")
"""Fitting stages that have a parameter list (spec D1 7.3). E is C and D together."""

_STAGE_ALIASES = {
    "A_activity": "A",
    "B_spectral": "B",
    "C_film": "C",
    "D_chemistry": "D",
    "E_polish": "E",
}

CALIBRATION_FILMS: tuple[str, ...] = ("kodak/tri-x-400", "kodak/t-max-100")
"""Materials with a film block (spec D1 7.1, 7.6)."""

CALIBRATION_PREFIX = "calibration:"
"""Name prefix of the parameters that exist only in the calibration (spec D1 7.1)."""

KINETICS_FILE_ID = "chemistry/kinetics"
"""File id of ``chemistry/kinetics.yaml``."""

_MATERIAL_DIRS = ("films", "papers")
_BATH_DIRS = ("developers", "baths")
_KINETICS_PATH = PurePosixPath("chemistry/kinetics.yaml")

# Tolerance for values that the inverse transform puts a few ulp outside their bounds.
_BOUND_SLACK = 1e-12


class ParameterError(ValueError):
    """An unknown parameter name, or a value that is not finite or lies outside its bounds."""


@dataclass(frozen=True)
class ParameterSpec:
    """One fitted parameter (spec D1 10, D1e ``calibrate/parameters.py``).

    ``lower`` and ``upper`` are the prior box in ``unit`` (infinite when the prior is a normal
    distribution without a box). ``prior_mu`` and ``prior_tau`` are in transform space. The
    prior residual is ``(g(value) - prior_mu) / prior_tau``, and ``prior_tau`` is ``inf`` for a box
    prior. ``block`` is the spec D1 7.1 block (A activity and pH nuisance, B spectral, F film,
    C chemistry). ``stage`` is the spec D1 7.3 stage that fits it.
    """

    name: str
    unit: str
    transform: Transform
    lower: float
    upper: float
    prior: PriorShape
    prior_mu: float
    prior_tau: float
    block: Block
    stage: str

    def __post_init__(self) -> None:
        if self.transform not in get_args(Transform):
            raise ValueError(f"{self.name}: unknown transform {self.transform!r}")
        if self.prior not in get_args(PriorShape):
            raise ValueError(f"{self.name}: unknown prior {self.prior!r}")
        if self.block not in get_args(Block):
            raise ValueError(f"{self.name}: unknown block {self.block!r}")
        if not self.lower <= self.upper:
            raise ValueError(f"{self.name}: bounds must have lower <= upper")
        if not self.prior_tau > 0.0:
            raise ValueError(f"{self.name}: prior_tau must be positive")
        if self.transform == "log10" and not self.lower >= 0.0:
            raise ValueError(f"{self.name}: a log10 parameter needs non-negative bounds")
        if self.transform == "logit" and not (self.lower >= 0.0 and self.upper <= 1.0):
            raise ValueError(f"{self.name}: a logit parameter needs bounds inside [0, 1]")


# --------------------------------------------------------------------------------------
# The parameter table of spec D1 7.1
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class _CalibrationOnly:
    """A parameter that has no data field (spec D1 7.1), with the spec's values."""

    key: str
    unit: str
    transform: Transform
    lower: float
    upper: float
    prior: PriorShape
    block: Block
    stage: str
    per_film: bool


# Carlton & Crabtree pH offset: linear, bounds +-0.5, uniform (spec D1 7.1 block A).
_DELTA_CC = _CalibrationOnly("delta_cc", "1", "linear", -0.5, 0.5, "uniform", "A", "A", False)
# Absolute absorption scale: log10 factor in [-1, 1], uniform (spec D1 7.1 block F).
_S_ABS = _CalibrationOnly("s_abs", "1", "linear", -1.0, 1.0, "uniform", "F", "C", True)
# T-Max aspect ratio AR = d_50 / t_g: log10, bounds [8, 19] (spec D1 7.1 block F, 2.0.2).
_ASPECT_RATIO = _CalibrationOnly("aspect_ratio", "1", "log10", 8.0, 19.0, "uniform", "F", "C", True)
_CALIBRATION_ONLY = {p.key: p for p in (_DELTA_CC, _S_ABS, _ASPECT_RATIO)}

# Stage A: activity coefficients, then the pH nuisances (bath id, species), spec D1 7.1.
_ACTIVITY_CLASSES = ("sulfite", "borate", "organic", "amine")
_STAGE_A_NUISANCE = (
    ("d-76", "boric-oxide"),
    ("t-max-dev", "sulfur-dioxide"),
    ("hc-110-pre2019", "potassium-hydroxide"),
)

# Stage B (per film): Urbach edge and dye bands, relative to the emulsion layer.
_INTRINSIC = (
    ("optics.intrinsic_absorption.alpha0", "log10"),
    ("optics.intrinsic_absorption.edge_energy", "linear"),
    ("optics.intrinsic_absorption.urbach_energy", "linear"),
)
_DYE_BAND = (("amplitude", "linear"), ("peak_wavenumber", "linear"), ("width", "linear"))

# Stage C (per film), in the order of the spec D1 7.1 block F row. "d_50" stands for
# grain.mean_diameter, or the aspect ratio for tabular grains; "s_abs" for the calibration-only
# scale. Tabular grains add grain.thickness at the end.
_FILM_BLOCK = (
    ("d_50", "log10"),
    ("grain.sigma_log", "linear"),
    ("latent.efficiency", "logit"),
    ("development.covering_power", "linear"),
    ("grain.silver_coverage", "log10"),
    ("development.fog_hazard_ratio", "log10"),
    ("optics.scattering_550", "linear"),
    ("s_abs", "linear"),
    ("optics.psf_weight", "logit"),
    ("optics.psf_length", "linear"),
    ("latent.tau0_median", "log10"),
    ("latent.tau0_sigma_ln", "linear"),
    ("latent.ionic_time", "log10"),
)
_TABULAR_EXTRA = (("grain.thickness", "linear"),)

# Stage D: developing-agent constants (spec D1 7.1 block C), shared chemistry, nuisances.
_AGENT_FIELDS = (
    ("metol", ("k_et", "K_ads", "lambda0", "k_des", "activation_energy")),
    ("hydroquinone", ("k_et", "K_ads", "lambda0", "k_reg", "activation_energy")),
    ("pyrazolidone", ("k_et", "K_ads", "lambda0", "k_des", "pka")),
    ("ascorbate", ("k_et", "K_ads", "lambda0", "k_reg")),
)
_AGENT_LOG_FIELDS = frozenset({"k_et", "K_ads", "lambda0", "k_des", "k_reg"})
_SHARED_CHEMISTRY = (
    ("halides.K_Br", "log10"),
    ("barriers.u_gel", "linear"),
    ("barriers.u_halide", "linear"),
    ("latent.u_li", "linear"),
    ("solvent.latent_attrition", "linear"),
    ("solvent.k_diss", "log10"),
    ("solvent.c_sat_ref", "log10"),
    ("solvent.sulfite_exponent", "linear"),
    ("solvent.k_pd", "log10"),
    ("transport.delta_max", "linear"),
    ("transport.delta_agit", "linear"),
)
_STAGE_D_NUISANCE = (
    ("t-max-dev", "hydroquinone"),
    ("t-max-dev", "dimezone-s"),
    ("t-max-dev", "sodium-bromide"),
    ("hc-110-pre2019", "hydroquinone"),
    ("hc-110-pre2019", "methylphenidone"),
    ("hc-110-pre2019", "potassium-bromide"),
)


# --------------------------------------------------------------------------------------
# Field paths
# --------------------------------------------------------------------------------------

_PART = re.compile(r"([^\[\]]+)((?:\[\d+\])*)")
_INDEX = re.compile(r"\[(\d+)\]")

Step = tuple[str, Any]


def _tokens(path: str) -> list[str | int]:
    if not isinstance(path, str) or not path:
        raise ParameterError(f"empty field path {path!r}")
    out: list[str | int] = []
    for part in path.split("."):
        match = _PART.fullmatch(part)
        if match is None:
            raise ParameterError(f"malformed field path {path!r}")
        out.append(match.group(1))
        out.extend(int(i) for i in _INDEX.findall(match.group(2)))
    return out


def _item_id(item: Any) -> str | None:
    for attr in ("id", "species"):
        value = getattr(item, attr, None)
        if isinstance(value, str):
            return value
    return None


def _greedy(tokens: list[str | int], start: int, known: Callable[[str], bool]) -> int:
    """End of the longest run ``tokens[start:end]`` of names whose dotted join is ``known``."""
    end = start
    while end < len(tokens) and isinstance(tokens[end], str):
        end += 1
    for stop in range(end, start, -1):
        if known(".".join(str(t) for t in tokens[start:stop])):
            return stop
    raise ParameterError(f"no field {tokens[start]!r}")


def _locate(obj: Any, tokens: list[str | int]) -> tuple[list[Step], Any]:
    """Follow ``tokens`` from ``obj``; return the copy steps and the object reached."""
    steps: list[Step] = []
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if isinstance(obj, BaseModel):
            if not isinstance(token, str) or token not in type(obj).model_fields:
                raise ParameterError(f"no field {token!r} in {type(obj).__name__}")
            steps.append(("attr", token))
            obj = getattr(obj, token)
            i += 1
        elif isinstance(obj, dict):
            stop = _greedy(tokens, i, lambda key, d=obj: key in d)
            key = ".".join(str(t) for t in tokens[i:stop])
            steps.append(("key", key))
            obj = obj[key]
            i = stop
        elif isinstance(obj, (list, tuple)):
            if isinstance(token, int):
                if not 0 <= token < len(obj):
                    raise ParameterError(f"index [{token}] out of range (length {len(obj)})")
                index = token
                i += 1
            else:
                ids = [_item_id(item) for item in obj]
                stop = _greedy(tokens, i, lambda key, known=ids: key in known)
                index = ids.index(".".join(str(t) for t in tokens[i:stop]))
                i = stop
            steps.append(("item", index))
            obj = obj[index]
        else:
            raise ParameterError(f"the path goes past a value at {token!r}")
    return steps, obj


def _replace(obj: Any, steps: list[Step], leaf: Any) -> Any:
    """A copy of ``obj`` with the object at ``steps`` replaced by ``leaf``; ``obj`` is unchanged."""
    if not steps:
        return leaf
    kind, key = steps[0]
    if kind == "attr":
        return obj.model_copy(update={key: _replace(getattr(obj, key), steps[1:], leaf)})
    if kind == "key":
        new = dict(obj)
        new[key] = _replace(obj[key], steps[1:], leaf)
        return new
    items = list(obj)
    items[key] = _replace(obj[key], steps[1:], leaf)
    return items if isinstance(obj, list) else tuple(items)


def _scaled(q: Quantity, value: float) -> Quantity:
    return q.model_copy(update={"value": float(value)})


def _set_value(obj: Any, tokens: list[str | int], value: float) -> Any:
    steps, q = _locate(obj, tokens)
    return _replace(obj, steps, _scaled(q, value))


def _root(library: Library) -> Library:
    """The library under every :class:`LibraryOverlay` of ``library``: the data files' values."""
    while isinstance(library, LibraryOverlay):
        library = library.base
    return library


def _raw(library: Library, kind: str, id_: str | None) -> Any:
    """The model ``kind``/``id_`` with the data overrides of every overlay of ``library`` merged
    in, before the derived parametrisations (``s_abs``, aspect ratio) are applied."""
    if isinstance(library, LibraryOverlay):
        return library._raw_model(kind, id_)
    if kind == "material":
        return library.material(id_)
    if kind == "bath":
        return library.bath(id_)
    return library.kinetics()


def _file_id(library: Library, kind: str, id_: str) -> str:
    path = library.path_of(kind, id_)
    return path.relative_to(library.root).with_suffix("").as_posix()


def _split_name(name: str) -> tuple[str, str]:
    if not isinstance(name, str):
        raise ParameterError(f"parameter names are strings, got {name!r}")
    file_id, sep, path = name.partition(":")
    if not sep or not file_id or not path:
        raise ParameterError(f"unknown parameter {name!r}: expected '<file id>:<field path>'")
    return file_id, path


def _open(library: Library, file_id: str, name: str) -> tuple[str, str | None, Any, Path]:
    """``(kind, id, model, path)`` of the data file ``file_id`` read through ``library``."""
    if file_id == KINETICS_FILE_ID:
        try:
            return "kinetics", None, _raw(library, "kinetics", None), library.root / _KINETICS_PATH
        except DataNotFound:
            raise ParameterError(f"unknown parameter {name!r}: no kinetics file") from None
    head, _, rest = file_id.partition("/")
    kind = "material" if head in _MATERIAL_DIRS else "bath" if head in _BATH_DIRS else None
    if kind is None or not rest:
        raise ParameterError(f"unknown parameter {name!r}: unknown file id {file_id!r}")
    try:
        if _file_id(library, kind, rest) != file_id:
            raise ParameterError(
                f"unknown parameter {name!r}: {file_id!r} is not the file of {rest!r}"
            )
        model = _raw(library, kind, rest)
    except DataNotFound:
        raise ParameterError(f"unknown parameter {name!r}: no file {file_id!r}") from None
    return kind, rest, model, library.path_of(kind, rest)


def _field_tokens(kind: str, path: str) -> list[str | int]:
    tokens = _tokens(path)
    return ["layers", *tokens] if kind == "material" else tokens


def _quantity(library: Library, name: str) -> tuple[Quantity, Path, str]:
    """The scalar :class:`Quantity` that data parameter ``name`` refers to, its file and path.

    Under an overlay this is the parameter's own value, before ``s_abs`` or the aspect ratio
    derive ``alpha0``, the dye amplitudes or ``mean_diameter`` from it.
    """
    file_id, path = _split_name(name)
    kind, _, model, file_path = _open(library, file_id, name)
    try:
        _, leaf = _locate(model, _field_tokens(kind, path))
    except ParameterError as exc:
        raise ParameterError(f"unknown parameter {name!r}: {exc}") from None
    if not isinstance(leaf, Quantity):
        raise ParameterError(f"unknown parameter {name!r}: {path!r} is not a quantity")
    if isinstance(leaf.value, list):
        raise ParameterError(f"parameter {name!r} is a list quantity, not a single value")
    return leaf, file_path, path


def _calibration_param(library: Library, name: str) -> tuple[_CalibrationOnly, str | None]:
    """The declaration and material id (or ``None``) of calibration-only parameter ``name``."""
    rest = name[len(CALIBRATION_PREFIX) :]
    film, _, key = rest.rpartition("/")
    decl = _CALIBRATION_ONLY.get(key)
    if decl is None or decl.per_film != bool(film):
        raise ParameterError(f"unknown parameter {name!r}")
    if not decl.per_film:
        return decl, None
    try:
        material = _raw(library, "material", film)
    except DataNotFound:
        raise ParameterError(f"unknown parameter {name!r}: no material {film!r}") from None
    if decl is _ASPECT_RATIO and _film_layer(material, library).grain.morphology != "tabular":
        raise ParameterError(f"unknown parameter {name!r}: {film!r} does not have tabular grains")
    return decl, film


def _film_layer(material: Material, library: Library) -> EmulsionLayer:
    layers = [
        layer
        for layer in material.emulsion_layers
        if layer.grain is not None and layer.optics is not None and layer.latent is not None
    ]
    if len(layers) != 1:
        path = library.path_of("material", material.id)
        raise DataError(
            f"{path}: layers: a film block needs exactly one emulsion layer with grain, optics "
            f"and latent fields (spec D1 7.1), found {len(layers)}"
        )
    return layers[0]


def _check_value(name: str, value: Any, lower: float, upper: float) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ParameterError(f"{name}: value must be a number, got {value!r}")
    v = float(value)
    if not math.isfinite(v):
        raise ParameterError(f"{name}: value must be a finite number, got {value!r}")
    if not lower <= v <= upper:
        raise ParameterError(f"{name}: value {v} is outside its bounds [{lower}, {upper}]")
    return v


# --------------------------------------------------------------------------------------
# Transforms
# --------------------------------------------------------------------------------------


def _forward(transform: str, x: float) -> float:
    if transform == "linear":
        return x
    if transform == "log10":
        return -math.inf if x == 0.0 else math.log10(x)
    if x == 0.0:
        return -math.inf
    if x == 1.0:
        return math.inf
    return math.log(x) - math.log1p(-x)


def _inverse(transform: str, t: float) -> float:
    if transform == "linear":
        return t
    if transform == "log10":
        return 10.0**t
    if t >= 0.0:
        return 1.0 / (1.0 + math.exp(-t))
    e = math.exp(t)
    return e / (1.0 + e)


def to_theta(specs: Sequence[ParameterSpec], values: Mapping[str, float]) -> np.ndarray:
    """Transformed parameter vector, in the order of ``specs``.

    ``values`` must hold exactly the names of ``specs``, each finite and inside its bounds.
    A value where the transform diverges (0 under log10 or logit, 1 under logit) raises.
    """
    names = [s.name for s in specs]
    missing = [n for n in names if n not in values]
    extra = sorted(set(values) - set(names))
    if missing or extra:
        raise ParameterError(f"values do not match the specs: missing {missing}, unknown {extra}")
    theta = np.empty(len(specs), dtype=np.float64)
    for i, spec in enumerate(specs):
        v = _check_value(spec.name, values[spec.name], spec.lower, spec.upper)
        t = _forward(spec.transform, v)
        if not math.isfinite(t):
            raise ParameterError(f"{spec.name}: {spec.transform} is not finite at {v}")
        theta[i] = t
    return theta


def from_theta(specs: Sequence[ParameterSpec], theta: np.ndarray) -> dict[str, float]:
    """Parameter values in their own units from the transformed vector ``theta``.

    A value that the inverse transform puts a few ulp outside its bounds is set to the bound.
    Values further outside are returned as they are; :class:`LibraryOverlay` rejects them.
    """
    t = np.asarray(theta, dtype=np.float64)
    if t.shape != (len(specs),):
        raise ParameterError(f"theta must have shape ({len(specs)},), got {t.shape}")
    if not np.all(np.isfinite(t)):
        raise ParameterError("theta must be finite")
    out: dict[str, float] = {}
    for spec, ti in zip(specs, t, strict=True):
        v = _inverse(spec.transform, float(ti))
        if v < spec.lower and spec.lower - v <= _BOUND_SLACK * max(1.0, abs(spec.lower)):
            v = spec.lower
        elif v > spec.upper and v - spec.upper <= _BOUND_SLACK * max(1.0, abs(spec.upper)):
            v = spec.upper
        out[spec.name] = v
    return out


def theta_bounds(specs: Sequence[ParameterSpec]) -> tuple[np.ndarray, np.ndarray]:
    """Lower and upper bounds of ``theta`` (infinite where the transform diverges)."""
    lo = np.array([_forward(s.transform, s.lower) for s in specs], dtype=np.float64)
    hi = np.array([_forward(s.transform, s.upper) for s in specs], dtype=np.float64)
    return lo, hi


# --------------------------------------------------------------------------------------
# Parameter lists
# --------------------------------------------------------------------------------------


def _data_spec(
    library: Library,
    name: str,
    transform: Transform,
    block: Block,
    stage: str,
    *,
    nuisance: bool = False,
) -> ParameterSpec:
    q, file_path, path = _quantity(library, name)
    where = f"{file_path}: {path}"
    prior = q.prior
    if prior is None:
        if not nuisance:
            raise DataError(f"{where}: a calibration parameter (spec D1 7.1) needs a prior")
        prior = "uniform"  # a nuisance amount is fitted within its bounds (spec D1 5.5)
    if q.bounds is not None:
        lower, upper = q.bounds
    elif prior == "normal" and q.uncertainty is not None:
        lower, upper = -math.inf, math.inf
    else:
        raise DataError(f"{where}: a calibration parameter (spec D1 7.1) needs bounds")
    if prior == "lognormal" and transform != "log10":
        raise DataError(f"{where}: a lognormal prior needs the log10 transform, not {transform}")
    if transform == "log10" and not lower > 0.0:
        raise DataError(f"{where}: a log10 parameter needs positive bounds, got [{lower}, {upper}]")
    if transform == "logit" and not (lower >= 0.0 and upper <= 1.0):
        raise DataError(f"{where}: a logit parameter needs bounds in [0, 1]")
    value = float(q.value)  # _quantity guarantees a single value
    mu = _forward(transform, value)
    if not math.isfinite(mu):
        raise DataError(f"{where}: value {value} is where the {transform} transform diverges")
    if prior == "normal":
        if q.uncertainty is None or not q.uncertainty > 0.0:
            raise DataError(f"{where}: a normal prior needs a positive uncertainty")
        if transform != "linear":
            raise DataError(f"{where}: a normal prior is only supported with a linear transform")
        tau = float(q.uncertainty)
    else:
        tau = math.inf  # a box prior has no width (spec D1 7.3 stage F)
    return ParameterSpec(name, q.unit, transform, lower, upper, prior, mu, tau, block, stage)


def _calibration_spec(decl: _CalibrationOnly, name: str, default: float) -> ParameterSpec:
    mu = _forward(decl.transform, default)
    return ParameterSpec(
        name,
        decl.unit,
        decl.transform,
        decl.lower,
        decl.upper,
        decl.prior,
        mu,
        math.inf,
        decl.block,
        decl.stage,
    )


def _calibration_default(library: Library, decl: _CalibrationOnly, film: str | None) -> float:
    if decl is not _ASPECT_RATIO:
        return 0.0
    material = _raw(library, "material", film)
    grain = _film_layer(material, library).grain
    ratio = float(grain.mean_diameter.value) / float(grain.thickness.value)
    if not decl.lower <= ratio <= decl.upper:
        raise DataError(
            f"{library.path_of('material', film)}: grain: mean_diameter / thickness = {ratio} is "
            f"outside the aspect-ratio bounds [{decl.lower}, {decl.upper}] (spec D1 7.1)"
        )
    return ratio


def _nuisance_specs(
    library: Library, entries: tuple[tuple[str, str], ...], block: Block, stage: str
) -> list[ParameterSpec]:
    specs = []
    for bath_id, species in entries:
        bath = library.bath(bath_id)
        if species not in bath.nuisance:
            raise DataError(
                f"{library.path_of('bath', bath_id)}: nuisance: {species!r} is a fitted nuisance "
                "(spec D1 7.1) but is not listed"
            )
        name = f"{_file_id(library, 'bath', bath_id)}:components.{species}.amount"
        specs.append(_data_spec(library, name, "linear", block, stage, nuisance=True))
    return specs


def _stage_a(library: Library) -> list[ParameterSpec]:
    specs = [
        _data_spec(library, f"{KINETICS_FILE_ID}:activity.delta_eps.{c}", "linear", "A", "A")
        for c in _ACTIVITY_CLASSES
    ]
    specs += _nuisance_specs(library, _STAGE_A_NUISANCE, "A", "A")
    name = CALIBRATION_PREFIX + _DELTA_CC.key
    specs.append(_calibration_spec(_DELTA_CC, name, 0.0))
    return specs


def _stage_b(library: Library, film: str) -> list[ParameterSpec]:
    material = library.material(film)
    layer = _film_layer(material, library)
    prefix = f"{_file_id(library, 'material', film)}:{layer.id}"
    specs = [
        _data_spec(library, f"{prefix}.{path}", transform, "B", "B")
        for path, transform in _INTRINSIC
    ]
    for j in range(len(layer.optics.dye_bands)):
        for field, transform in _DYE_BAND:
            name = f"{prefix}.optics.dye_bands[{j}].{field}"
            specs.append(_data_spec(library, name, transform, "B", "B"))
    return specs


def _stage_c(library: Library, film: str) -> list[ParameterSpec]:
    material = library.material(film)
    layer = _film_layer(material, library)
    tabular = layer.grain.morphology == "tabular"
    prefix = f"{_file_id(library, 'material', film)}:{layer.id}"
    specs = []
    for key, transform in _FILM_BLOCK + (_TABULAR_EXTRA if tabular else ()):
        if key == "s_abs" or (key == "d_50" and tabular):
            decl = _S_ABS if key == "s_abs" else _ASPECT_RATIO
            name = f"{CALIBRATION_PREFIX}{film}/{decl.key}"
            specs.append(_calibration_spec(decl, name, _calibration_default(library, decl, film)))
            continue
        path = "grain.mean_diameter" if key == "d_50" else key
        specs.append(_data_spec(library, f"{prefix}.{path}", transform, "F", "C"))
    return specs


def _stage_d(library: Library) -> list[ParameterSpec]:
    specs = []
    for agent, fields in _AGENT_FIELDS:
        for field in fields:
            transform: Transform = "log10" if field in _AGENT_LOG_FIELDS else "linear"
            name = f"{KINETICS_FILE_ID}:agents.{agent}.{field}"
            specs.append(_data_spec(library, name, transform, "C", "D"))
    for path, transform in _SHARED_CHEMISTRY:
        specs.append(_data_spec(library, f"{KINETICS_FILE_ID}:{path}", transform, "C", "D"))
    specs += _nuisance_specs(library, _STAGE_D_NUISANCE, "C", "D")
    return specs


def _stage(stage: str) -> str:
    key = _STAGE_ALIASES.get(stage, stage) if isinstance(stage, str) else None
    if key not in STAGES:
        allowed = ", ".join([*STAGES, *_STAGE_ALIASES])
        raise ValueError(f"unknown stage {stage!r}; expected one of: {allowed}")
    return key


def parameter_specs(
    library: Library, stage: str, *, films: Sequence[str] = CALIBRATION_FILMS
) -> list[ParameterSpec]:
    """The parameters that spec D1 7.1 assigns to ``stage`` (spec D1 7.3), in table order.

    ``stage`` is ``"A"`` to ``"E"`` or the matching :class:`StageRecord` name (``"A_activity"``,
    ...). B and C list each of ``films`` in turn. E is C followed by D. Every field, ``prior_mu``
    included, comes from the data files under any :class:`LibraryOverlay`: the prior stays
    anchored at the file value whatever an overlay holds (spec D1 7.1), so the result for an
    overlay equals the result for its root library.
    """
    library = _root(library)
    key = _stage(stage)
    if key == "A":
        specs = _stage_a(library)
    elif key == "B":
        specs = [s for film in films for s in _stage_b(library, film)]
    elif key == "C":
        specs = [s for film in films for s in _stage_c(library, film)]
    elif key == "D":
        specs = _stage_d(library)
    else:
        specs = parameter_specs(library, "C", films=films) + parameter_specs(library, "D")
    names = [s.name for s in specs]
    if len(set(names)) != len(names):  # pragma: no cover - the tables above are disjoint
        raise AssertionError("duplicate parameter names")
    return specs


def parameter_value(library: Library, name: str) -> float:
    """The current value of parameter ``name`` as ``library`` sees it.

    This is the parameter's own value: under ``s_abs`` it is the unscaled ``alpha0`` or dye
    amplitude, so the values read back feed a new overlay without scaling them twice. A calibration-only parameter that no overlay sets has its neutral value: 0 for
    ``delta_cc`` and ``s_abs``, and ``mean_diameter / thickness`` for the aspect ratio.
    """
    if isinstance(name, str) and name.startswith(CALIBRATION_PREFIX):
        lib = library
        if isinstance(lib, LibraryOverlay) and name in lib._calibration_chain():
            return lib._calibration_chain()[name]
        decl, film = _calibration_param(library, name)
        return _calibration_default(library, decl, film)
    q, _, _ = _quantity(library, name)
    return float(q.value)


def current_values(library: Library, specs: Sequence[ParameterSpec]) -> dict[str, float]:
    """``{name: value}`` of every spec, read through ``library``."""
    return {s.name: parameter_value(library, s.name) for s in specs}


# --------------------------------------------------------------------------------------
# Overlay
# --------------------------------------------------------------------------------------


class LibraryOverlay(Library):
    """A library that returns ``base``'s models with ``overrides`` applied in memory.

    ``overrides`` maps parameter names (see the module docstring) to values in their own
    units. Any scalar :class:`Quantity` of a material, bath or the kinetics file may be
    overridden, and so may the calibration-only parameters. An unknown name, a non-finite
    value, or a value outside the quantity's ``bounds`` (or the spec's box for a calibration-only
    parameter) raises :class:`ParameterError`. Every read goes through ``base`` (a
    ``TrainingLibrary`` in the calibration, spec D1 7.7), and nothing is written to disk.

    ``base`` may itself be an overlay. The data overrides of the whole chain are merged first
    (an outer value wins), and ``s_abs`` and the aspect ratio are applied once, on the merged
    material, with the outermost value of each.
    """

    def __init__(self, base: Library, overrides: Mapping[str, float]) -> None:
        if not isinstance(base, Library):
            raise TypeError(f"base must be a Library, got {type(base).__name__}")
        super().__init__(base.root)
        data: dict[str, dict[str, float]] = {}
        calibration: dict[str, float] = {}
        for name, value in dict(overrides).items():
            if isinstance(name, str) and name.startswith(CALIBRATION_PREFIX):
                decl, _ = _calibration_param(base, name)
                calibration[name] = _check_value(name, value, decl.lower, decl.upper)
                continue
            q, _, path = _quantity(base, name)
            lower, upper = q.bounds if q.bounds is not None else (-math.inf, math.inf)
            file_id = name.partition(":")[0]
            data.setdefault(file_id, {})[path] = _check_value(name, value, lower, upper)
        self._base = base
        self._data = {k: MappingProxyType(v) for k, v in data.items()}
        self._calibration = MappingProxyType(calibration)
        self._overrides = MappingProxyType(dict(overrides))
        chain_data = self._data_chain()
        for name in self._calibration_chain():
            film, _, key = name[len(CALIBRATION_PREFIX) :].rpartition("/")
            if key != _ASPECT_RATIO.key:
                continue
            layer = _film_layer(_raw(base, "material", film), base)
            path = f"{layer.id}.grain.mean_diameter"
            if path in chain_data.get(_file_id(base, "material", film), {}):
                raise ParameterError(
                    f"{name} and {_file_id(base, 'material', film)}:{path} both set the grain "
                    "diameter; give one of them"
                )

    def __repr__(self) -> str:
        return f"LibraryOverlay({self._base!r}, {len(self._overrides)} overrides)"

    @property
    def base(self) -> Library:
        """The library that every read is delegated to."""
        return self._base

    @property
    def overrides(self) -> Mapping[str, float]:
        """The overrides as given (read-only)."""
        return self._overrides

    @property
    def calibration_values(self) -> Mapping[str, float]:
        """This overlay's own calibration-only overrides (``calibration:...``), read-only."""
        return self._calibration

    # ---- the overlay chain ------------------------------------------------------------

    def _calibration_chain(self) -> dict[str, float]:
        """Calibration-only values of the whole chain; an outer overlay's value wins."""
        inner = self._base._calibration_chain() if isinstance(self._base, LibraryOverlay) else {}
        return {**inner, **self._calibration}

    def _data_chain(self) -> dict[str, dict[str, float]]:
        """``{file id: {path: value}}`` of the whole chain; an outer overlay's value wins."""
        inner = self._base._data_chain() if isinstance(self._base, LibraryOverlay) else {}
        merged = {k: dict(v) for k, v in inner.items()}
        for file_id, values in self._data.items():
            merged.setdefault(file_id, {}).update(values)
        return merged

    def _raw_model(self, kind: str, id_: str | None) -> Any:
        """The model with the chain's data overrides and no derived parametrisation."""
        model = _raw(self._base, kind, id_)
        file_id = KINETICS_FILE_ID if kind == "kinetics" else _file_id(self._base, kind, id_)
        for path, value in self._data.get(file_id, {}).items():
            model = _set_value(model, _field_tokens(kind, path), value)
        return model

    # ---- overridden models ------------------------------------------------------------

    def material(self, id: str) -> Material:
        material = self._raw_model("material", id)
        calibration = self._calibration_chain()
        s_abs = calibration.get(f"{CALIBRATION_PREFIX}{id}/{_S_ABS.key}")
        ratio = calibration.get(f"{CALIBRATION_PREFIX}{id}/{_ASPECT_RATIO.key}")
        if s_abs is None and ratio is None:
            return material
        layer = _film_layer(material, self._base)
        base_tokens: list[str | int] = ["layers", layer.id]
        if ratio is not None:
            thickness = float(layer.grain.thickness.value)
            tokens = [*base_tokens, "grain", "mean_diameter"]
            material = _set_value(material, tokens, ratio * thickness)
        if s_abs is not None:
            factor = 10.0**s_abs
            optics = layer.optics
            tokens = [*base_tokens, "optics", "intrinsic_absorption", "alpha0"]
            alpha0 = float(optics.intrinsic_absorption.alpha0.value)
            material = _set_value(material, tokens, alpha0 * factor)
            for j, band in enumerate(optics.dye_bands):
                tokens = [*base_tokens, "optics", "dye_bands", j, "amplitude"]
                material = _set_value(material, tokens, float(band.amplitude.value) * factor)
        return material

    def bath(self, id: str) -> Bath:
        return self._raw_model("bath", id)

    def kinetics(self) -> KineticsFile:
        return self._raw_model("kinetics", None)

    # ---- everything else is the base library's ---------------------------------------

    def ids(self, kind: Kind) -> list[str]:
        return self._base.ids(kind)

    def path_of(self, kind: Kind, id: str) -> Path:
        return self._base.path_of(kind, id)

    def process(self, id: str) -> Process:
        return self._base.process(id)

    def species(self) -> dict[str, Species]:
        return self._base.species()

    def calibration_record(self, run_id: str) -> CalibrationRecord:
        return self._base.calibration_record(run_id)

    def reference_tables(self, material_id: str) -> dict:
        return self._base.reference_tables(material_id)

    def developer_facts(self, developer: str) -> dict:
        return self._base.developer_facts(developer)

    def reference_curve(self, relpath: str | PurePosixPath) -> ReferenceCurve:
        return self._base.reference_curve(relpath)

    def resolve_curve(self, curve: SpectralCurve, grid: SpectralGrid) -> np.ndarray:
        return self._base.resolve_curve(curve, grid)

    def validate(self) -> list[str]:
        """Check the base library's files (the overrides live only in memory)."""
        return self._base.validate()
