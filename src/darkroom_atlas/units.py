"""Physical constants, unit conversion, and parsers for durations and dilutions.

Unit strings are ASCII only (``um`` for micrometre, ``C`` for degree Celsius). Conversion
factors are exact integers relative to the smallest unit of each dimension, so that
conversions such as ``mm -> um`` or ``min -> s`` are exact in floating point.
"""

import math
import re
from numbers import Real

import numpy as np

from darkroom_atlas.backend import xp_of

__all__ = [
    "C_LIGHT",
    "H_PLANCK",
    "K_BOLTZMANN",
    "K_M",
    "N_AVOGADRO",
    "R_GAS",
    "convert",
    "density_to_transmittance",
    "parse_dilution",
    "parse_duration_s",
    "photon_energy_j",
    "pixel_pitch_um",
    "transmittance_to_density",
]

H_PLANCK = 6.62607015e-34
"""Planck constant, J s (exact, SI 2019)."""

C_LIGHT = 299792458.0
"""Speed of light in vacuum, m / s (exact)."""

N_AVOGADRO = 6.02214076e23
"""Avogadro constant, 1 / mol (exact, SI 2019)."""

R_GAS = 8.314462618
"""Molar gas constant, J / (mol K)."""

K_BOLTZMANN = 1.380649e-23
"""Boltzmann constant, J / K (exact, SI 2019)."""

K_M = 683.0
"""Maximum luminous efficacy of photopic vision, lm / W."""


# --------------------------------------------------------------------------------------
# Array API helpers
# --------------------------------------------------------------------------------------


def _as_float_array(x):
    """Return ``(xp, array, was_scalar)`` with ``array`` in a real floating dtype.

    Python scalars and other non-array inputs become NumPy ``float64`` arrays; arrays that
    already belong to an Array API namespace keep it. Integer arrays are cast to float64.
    """
    try:
        xp = xp_of(x)
        was_scalar = False
    except TypeError:
        was_scalar = np.ndim(x) == 0
        x = np.asarray(x, dtype=np.float64)
        xp = xp_of(x)
    if not xp.isdtype(x.dtype, "real floating"):
        x = xp.astype(x, xp.float64)
    return xp, x, was_scalar


def _finish(result, was_scalar):
    return float(result) if was_scalar else result


def photon_energy_j(wavelength_nm):
    """Energy of one photon, in J, at ``wavelength_nm`` (``E = h c / lambda``)."""
    _, wl, was_scalar = _as_float_array(wavelength_nm)
    return _finish((H_PLANCK * C_LIGHT * 1e9) / wl, was_scalar)


def transmittance_to_density(t, *, floor=1e-12):
    """Base-10 optical density ``-log10(t)``, with ``t`` clipped below at ``floor``.

    The floor keeps the density finite for opaque samples: ``t == 0`` maps to
    ``-log10(floor)`` (12 by default).
    """
    if not floor > 0.0:
        raise ValueError(f"floor must be positive, got {floor!r}")
    xp, t, was_scalar = _as_float_array(t)
    floor_arr = xp.asarray(floor, dtype=t.dtype)
    return _finish(-xp.log10(xp.maximum(t, floor_arr)), was_scalar)


def density_to_transmittance(d):
    """Transmittance ``10 ** (-d)`` for a base-10 optical density ``d``."""
    _, d, was_scalar = _as_float_array(d)
    return _finish(10.0 ** (-d), was_scalar)


# --------------------------------------------------------------------------------------
# Unit conversion
# --------------------------------------------------------------------------------------

# Each unit maps to (dimension, factor), where ``factor`` is the size of the unit in the
# dimension's smallest unit. Integer factors keep the common conversions exact.
_LINEAR_UNITS: dict[str, tuple[str, float]] = {
    "nm": ("length", 1.0),
    "um": ("length", 1e3),
    "mm": ("length", 1e6),
    "cm": ("length", 1e7),
    "m": ("length", 1e9),
    "s": ("time", 1.0),
    "min": ("time", 60.0),
    "h": ("time", 3600.0),
    "mL": ("volume", 1.0),
    "L": ("volume", 1e3),
    "mg/L": ("mass concentration", 1.0),
    "g/L": ("mass concentration", 1e3),
    "g/mol": ("molar mass", 1.0),
    "g/m2": ("areal mass", 1.0),
    "1": ("dimensionless", 1.0),
}

# Temperatures are affine: kelvin = celsius + offset.
_TEMPERATURE_OFFSET_K: dict[str, float] = {"K": 0.0, "C": 273.15}

_KNOWN_UNITS = sorted([*_LINEAR_UNITS, *_TEMPERATURE_OFFSET_K])


def _dimension_of(unit: str) -> str:
    if unit in _TEMPERATURE_OFFSET_K:
        return "temperature"
    if unit in _LINEAR_UNITS:
        return _LINEAR_UNITS[unit][0]
    raise ValueError(f"unknown unit {unit!r}; known units: {', '.join(_KNOWN_UNITS)}")


def convert(value, from_unit: str, to_unit: str):
    """Convert ``value`` from ``from_unit`` to ``to_unit``.

    Supported: length ``nm um mm cm m``, time ``s min h``, temperature ``C K``, volume
    ``mL L``, mass concentration ``mg/L g/L``. ``value`` may be a scalar or an array; lists
    and tuples are converted to NumPy arrays. If the two units are identical, ``value`` is
    returned unchanged (without checking the unit). Unknown units and conversions between
    different dimensions raise ``ValueError``.
    """
    if from_unit == to_unit:
        return value
    from_dim = _dimension_of(from_unit)
    to_dim = _dimension_of(to_unit)
    if from_dim != to_dim:
        raise ValueError(f"cannot convert {from_unit!r} ({from_dim}) to {to_unit!r} ({to_dim})")
    if isinstance(value, (list, tuple)):
        value = np.asarray(value, dtype=np.float64)
    if from_dim == "temperature":
        kelvin = value + _TEMPERATURE_OFFSET_K[from_unit]
        return kelvin - _TEMPERATURE_OFFSET_K[to_unit]
    return value * _LINEAR_UNITS[from_unit][1] / _LINEAR_UNITS[to_unit][1]


# --------------------------------------------------------------------------------------
# Durations and dilutions
# --------------------------------------------------------------------------------------

_NUMBER = r"(?:\d+(?:\.\d*)?|\.\d+)"
_SUFFIX_SECONDS = {"": 1.0, "s": 1.0, "min": 60.0, "h": 3600.0}
_SUFFIXED_RE = re.compile(rf"^({_NUMBER})\s*(s|min|h)?$")
# h:mm:ss or m:ss; only the last field may carry a fractional part.
_CLOCK_RE = re.compile(rf"^\d+(?::\d+)?:{_NUMBER}$")


def _check_duration(seconds: float, text) -> float:
    if not math.isfinite(seconds) or seconds < 0.0:
        raise ValueError(f"duration must be a finite, non-negative number of seconds: {text!r}")
    return float(seconds)


# The parsers below raise ValueError even for inputs of the wrong type, so that pydantic
# validators calling them report every problem as an ordinary validation error.


def parse_duration_s(text: str | float | int) -> float:  # noqa: PYI041 - ints named on purpose
    """Parse a duration and return it in seconds.

    Accepted forms: ``m:ss`` (``"9:45"``), ``h:mm:ss`` (``"1:02:03"``), a number with an
    ``s``, ``min`` or ``h`` suffix (``"30s"``, ``"1.5min"``), and a bare number of seconds
    (``90`` or ``"90"``). In the colon forms every field after the first must be below 60.
    Negative, non-finite and malformed values (and booleans) raise ``ValueError``.
    """
    if isinstance(text, Real) and not isinstance(text, bool):
        return _check_duration(float(text), text)
    if not isinstance(text, str):
        raise ValueError(f"duration must be a number or a string, got {text!r}")  # noqa: TRY004

    s = text.strip()
    if not s:
        raise ValueError("empty duration")

    if ":" in s:
        if _CLOCK_RE.match(s) is None:
            raise ValueError(f"invalid duration {text!r}; expected m:ss or h:mm:ss")
        values = [float(f) for f in s.split(":")]
        if any(v >= 60.0 for v in values[1:]):
            raise ValueError(f"invalid duration {text!r}: minutes and seconds must be below 60")
        seconds = 0.0
        for v in values:
            seconds = seconds * 60.0 + v
        return _check_duration(seconds, text)

    match = _SUFFIXED_RE.match(s)
    if match is None:
        raise ValueError(
            f"invalid duration {text!r}; expected m:ss, h:mm:ss, or a number with an "
            "optional s, min or h suffix"
        )
    number, suffix = match.groups()
    return _check_duration(float(number) * _SUFFIX_SECONDS[suffix or ""], text)


_DILUTION_RE = re.compile(rf"^({_NUMBER})\s*\+\s*({_NUMBER})$")


def parse_dilution(text: str) -> float:
    """Return the fraction of stock solution in the final volume.

    ``"stock"`` is undiluted (1.0). ``"a+b"`` means ``a`` parts stock plus ``b`` parts water,
    so ``"1+1"`` is 0.5 and ``"1+50"`` is 1/51. The ratio form ``"a:b"`` is rejected because
    it is used inconsistently in the literature (``1:4`` can mean 1+4 or 1 part in 4 total).
    """
    if not isinstance(text, str):
        message = f"dilution must be a string such as 'stock' or '1+1', got {text!r}"
        raise ValueError(message)  # noqa: TRY004
    s = text.strip()
    if s.lower() == "stock":
        return 1.0
    if ":" in s:
        raise ValueError(
            f"ambiguous dilution {text!r}: write 'a+b' (a parts stock plus b parts water) "
            "or 'stock'"
        )
    match = _DILUTION_RE.match(s)
    if match is None:
        raise ValueError(f"invalid dilution {text!r}; expected 'stock' or 'a+b', e.g. '1+1'")
    stock = float(match.group(1))
    water = float(match.group(2))
    if not stock > 0.0:
        raise ValueError(f"invalid dilution {text!r}: the stock part must be positive")
    return stock / (stock + water)


# --------------------------------------------------------------------------------------
# Film-plane geometry
# --------------------------------------------------------------------------------------


def pixel_pitch_um(long_edge_mm: float, shape_hw: tuple[int, int]) -> float:
    """Pixel pitch in um when the long edge of a ``shape_hw`` image spans ``long_edge_mm``."""
    if not (math.isfinite(long_edge_mm) and long_edge_mm > 0.0):
        raise ValueError(f"long_edge_mm must be positive and finite, got {long_edge_mm!r}")
    if len(shape_hw) != 2 or any(int(n) != n or n < 1 for n in shape_hw):
        raise ValueError(f"shape_hw must be two positive integers, got {shape_hw!r}")
    return float(long_edge_mm) * 1000.0 / max(int(n) for n in shape_hw)
