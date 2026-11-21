"""SIT-type activity corrections and mixed acid-base constants (spec D1 2.3.3).

``log10 gamma_i = -z_i^2 D(I) + eps_i m_cat`` with ``D(I) = A(T) sqrt(I) / (1 + b sqrt(I))``.
pH is the hydrogen-ion activity (glass-electrode convention), so the mixed constant of the
step ``HA^z <=> A^(z-1) + H+`` relates the concentration ratio to ``a_H``::

    pK'(T, I) = pK(T) + (2 z_HA - 1) D(I) + d_eps_class m_cat

Water: ``[H+] = a_H / gamma_H`` and ``[OH-] = K_w / (a_H gamma_OH)``, with ``eps = 0`` for
both ions.
"""

import math

import numpy as np

from darkroom_atlas.chemistry.species import (
    GAS_CONSTANT,
    REFERENCE_TEMPERATURE_K,
    AcidBaseFamily,
    _check_temperature,
)
from darkroom_atlas.data.models import ActivityParams

__all__ = [
    "PKW_20C",
    "PKW_25C",
    "WATER_IONISATION_ENTHALPY_KJ_MOL",
    "conditional_pka",
    "debye_huckel_a",
    "debye_huckel_d",
    "delta_eps",
    "pkw",
    "water_constants",
]

_T20 = 293.15
_T25 = REFERENCE_TEMPERATURE_K

PKW_20C = 14.17
"""pK_w at 20 C (spec D1 2.3.3, 2.3.6; DEV 6.2, recalled -> assumed)."""

PKW_25C = 14.00
"""pK_w at 25 C (spec D1 2.3.3, 2.3.6; DEV 6.2, recalled -> assumed)."""

WATER_IONISATION_ENTHALPY_KJ_MOL = (
    (PKW_20C - PKW_25C) * GAS_CONSTANT * math.log(10.0) / (1.0 / _T20 - 1.0 / _T25) / 1e3
)
"""Delta H_w of the van 't Hoff line through the two pK_w values, about 56.9 kJ/mol
(spec D1 2.3.3: derived from the two values)."""


def _check_ionic_strength(ionic_strength: float) -> float:
    i = float(ionic_strength)
    if not math.isfinite(i) or i < 0.0:
        raise ValueError(f"ionic strength must be a finite number >= 0 mol/L, got {ionic_strength}")
    return i


def debye_huckel_a(temperature_k: float, params: ActivityParams) -> float:
    """``A(T)``: linear in T through the 20 C and 25 C values, extrapolated (spec D1 2.3.3)."""
    t = _check_temperature(temperature_k)
    a20 = float(params.debye_huckel_A_20C.to("1"))
    a25 = float(params.debye_huckel_A_25C.to("1"))
    return a20 + (a25 - a20) * (t - _T20) / (_T25 - _T20)


def debye_huckel_d(ionic_strength: float, temperature_k: float, params: ActivityParams) -> float:
    """``D(I) = A(T) sqrt(I) / (1 + b sqrt(I))``, ``b = params.sit_b`` (spec D1 2.3.3)."""
    i = _check_ionic_strength(ionic_strength)
    root = math.sqrt(i)
    return debye_huckel_a(temperature_k, params) * root / (1.0 + float(params.sit_b.to("1")) * root)


def delta_eps(activity_class: str, params: ActivityParams) -> float:
    """SIT coefficient of an activity class in kg/mol; 0 for ``none`` or a class not listed."""
    q = params.delta_eps.get(activity_class)
    return 0.0 if q is None else float(q.to("kg/mol"))


def _check_molality(cation_molality: float) -> float:
    m = float(cation_molality)
    if not math.isfinite(m) or m < 0.0:
        raise ValueError(f"cation molality must be a finite number >= 0, got {cation_molality}")
    return m


def conditional_pka(
    family: AcidBaseFamily,
    temperature_k: float,
    ionic_strength: float,
    cation_molality: float,
    params: ActivityParams,
) -> np.ndarray:
    """Mixed stepwise constants ``pK'_j``, shape ``(J,)``, read-only (spec D1 2.3.3)."""
    d = debye_huckel_d(ionic_strength, temperature_k, params)
    m = _check_molality(cation_molality)
    z_ha = np.asarray(family.charges[:-1], dtype=np.float64)
    out = (
        family.pka(temperature_k)
        + (2.0 * z_ha - 1.0) * d
        + delta_eps(family.activity_class, params) * m
    )
    out = np.array(out, dtype=np.float64)
    out.flags.writeable = False
    return out


def pkw(temperature_k: float) -> float:
    """Thermodynamic pK_w from the van 't Hoff line through 14.17 (20 C) and 14.00 (25 C)."""
    t = _check_temperature(temperature_k)
    return PKW_25C + WATER_IONISATION_ENTHALPY_KJ_MOL * 1e3 / (GAS_CONSTANT * math.log(10.0)) * (
        1.0 / t - 1.0 / _T25
    )


def water_constants(
    temperature_k: float, ionic_strength: float, params: ActivityParams
) -> tuple[float, float, float]:
    """``(log10 K_w, log10 gamma_H, log10 gamma_OH)``; both ions have ``z = 1``, ``eps = 0``."""
    d = debye_huckel_d(ionic_strength, temperature_k, params)
    return -pkw(temperature_k), -d, -d
