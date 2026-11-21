"""Acid-base families and transport data of the solution components (spec D1 2.3.1, 2.3.7).

A component's acid-base family lists its protonation states, most protonated first. State
``j`` has charge ``z_0 - j`` and ``n_j = reference_protons - j`` protons relative to the
reference state (spec D1 2.3.1). The family has ``len(pka25) + 1`` states.

Temperature laws of pK (spec D1 2.3.3):

* ``van_t_hoff``: ``pK(T) = pK(T_r) + dH / (R ln 10) * (1/T - 1/T_r)``, ``T_r = 298.15 K``;
* ``manov_boric``: ``pK(T) = 2237.94/T + 0.016883 T - 3.305`` (Manov et al. 1944), used for
  every step of the family; the stored 25 C value is not used;
* ``linear``: ``pK(T) = pK(25 C) + dpK/dT (T - 298.15)``.
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Literal

import numpy as np

from darkroom_atlas.data.models import AgentRole, Species

__all__ = [
    "GAS_CONSTANT",
    "REFERENCE_TEMPERATURE_K",
    "AcidBaseFamily",
    "ComponentInfo",
    "component_table",
    "manov_boric_pka",
]

GAS_CONSTANT = 8.314462618
"""Molar gas constant, J/(mol K) (CODATA exact value; spec D1 2.3.3 van 't Hoff law)."""

REFERENCE_TEMPERATURE_K = 298.15
"""Reference temperature of the tabulated pK values, 25 C (spec D1 2.3.3, 5.2)."""

# Manov, DeLollis & Acree 1944 boric acid law (spec D1 2.3.3; DEV 6.2).
_MANOV_A = 2237.94
_MANOV_B = 0.016883
_MANOV_C = -3.305


def manov_boric_pka(temperature_k: float) -> float:
    """Thermodynamic pK of boric acid at ``temperature_k`` (spec D1 2.3.3)."""
    return _MANOV_A / temperature_k + _MANOV_B * temperature_k + _MANOV_C


def _check_temperature(temperature_k: float) -> float:
    t = float(temperature_k)
    if not math.isfinite(t) or t <= 0.0:
        raise ValueError(f"temperature must be a positive number of kelvin, got {temperature_k}")
    return t


@dataclass(frozen=True)
class AcidBaseFamily:
    """Protonation states and pK values of one component (spec D1 2.3.1, 10 D1c)."""

    component: str
    charges: tuple[int, ...]  # per state, most protonated first
    protons_rel: tuple[int, ...]  # per state, relative to the reference state
    pka25: tuple[float, ...]
    temperature_law: Literal["van_t_hoff", "manov_boric", "linear"]
    enthalpy_kj_mol: tuple[float, ...]
    dpka_dt: tuple[float, ...]
    activity_class: str

    def __post_init__(self) -> None:
        n = len(self.pka25)
        if len(self.charges) != n + 1 or len(self.protons_rel) != n + 1:
            raise ValueError(
                f"{self.component}: charges and protons_rel need {n + 1} entries "
                f"(len(pka25) + 1), got {len(self.charges)} and {len(self.protons_rel)}"
            )
        if self.temperature_law == "van_t_hoff" and len(self.enthalpy_kj_mol) != n:
            raise ValueError(f"{self.component}: enthalpy_kj_mol needs one entry per pK")
        if self.temperature_law == "linear" and len(self.dpka_dt) != n:
            raise ValueError(f"{self.component}: dpka_dt needs one entry per pK")
        if self.temperature_law not in ("van_t_hoff", "manov_boric", "linear"):
            raise ValueError(f"{self.component}: unknown temperature law {self.temperature_law!r}")

    @property
    def n_states(self) -> int:
        """Number of protonation states."""
        return len(self.charges)

    @classmethod
    def from_species(cls, species: Species) -> "AcidBaseFamily":
        """Family of a component species; raises ``ValueError`` if it has no ``acid_base``."""
        ab = species.acid_base
        if ab is None:
            raise ValueError(f"species {species.id!r} has no acid_base")
        n = len(ab.pka)
        enthalpy = (
            tuple(float(q.to("kJ/mol")) for q in ab.ionisation_enthalpy)
            if ab.ionisation_enthalpy
            else (0.0,) * n
        )
        return cls(
            component=species.id,
            charges=tuple(ab.charge_most_protonated - j for j in range(n + 1)),
            protons_rel=tuple(ab.reference_protons - j for j in range(n + 1)),
            pka25=tuple(float(q.to("1")) for q in ab.pka),
            temperature_law=ab.temperature_law,
            enthalpy_kj_mol=enthalpy,
            dpka_dt=tuple(float(q.to("1/K")) for q in ab.dpka_dt),
            activity_class=ab.activity_class,
        )

    def pka(self, temperature_k: float) -> np.ndarray:
        """Thermodynamic stepwise pK at ``temperature_k`` (I -> 0), shape ``(J,)``, read-only."""
        t = _check_temperature(temperature_k)
        pk25 = np.asarray(self.pka25, dtype=np.float64)
        if self.temperature_law == "van_t_hoff":
            dh = np.asarray(self.enthalpy_kj_mol, dtype=np.float64) * 1e3  # J/mol
            out = pk25 + dh / (GAS_CONSTANT * math.log(10.0)) * (
                1.0 / t - 1.0 / REFERENCE_TEMPERATURE_K
            )
        elif self.temperature_law == "manov_boric":
            out = np.full(pk25.shape, manov_boric_pka(t))
        else:
            out = pk25 + np.asarray(self.dpka_dt, dtype=np.float64) * (t - REFERENCE_TEMPERATURE_K)
        out = np.array(out, dtype=np.float64)
        out.flags.writeable = False
        return out


@dataclass(frozen=True)
class ComponentInfo:
    """Transport and acid-base data of one component (spec D1 2.3.7, 10 D1c).

    ``d_water_um2_s`` is the diffusion coefficient in water at 20 C and ``gel_ratio`` is
    ``D_gel / D_w``; both are 0 for an immobile component (gelatin).
    """

    id: str
    family: AcidBaseFamily | None
    d_water_um2_s: float  # 20 C
    gel_ratio: float
    immobile: bool
    agent: AgentRole | None


def component_table(species: Mapping[str, Species]) -> dict[str, ComponentInfo]:
    """All species of kind ``component``, by id, in the order of ``species``."""
    table: dict[str, ComponentInfo] = {}
    for key, s in species.items():
        if s.kind != "component":
            continue
        if key != s.id:
            raise ValueError(f"species mapping key {key!r} differs from the species id {s.id!r}")
        d_water = (
            0.0 if s.diffusion_coefficient is None else float(s.diffusion_coefficient.to("um2/s"))
        )
        gel_ratio = 0.0 if s.gel_diffusion_ratio is None else float(s.gel_diffusion_ratio.to("1"))
        if s.immobile:
            d_water, gel_ratio = 0.0, 0.0
        table[s.id] = ComponentInfo(
            id=s.id,
            family=None if s.acid_base is None else AcidBaseFamily.from_species(s),
            d_water_um2_s=d_water,
            gel_ratio=gel_ratio,
            immobile=s.immobile,
            agent=s.developing_agent,
        )
    return table
