"""Development kinetics: resolved constants, surface state and initiation hazards
(spec D1 2.4.2, 2.4.3, 2.4.5, 2.4.6, 10 D1c).

Surface state of the competitive Langmuir-Hinshelwood rate law with ETA turnover
(spec D1 2.4.3), with ``a_i`` the active-form concentration of agent ``i``::

    Xi      = 1 + sum_E K_E a_E (1 + Phi_E) + sum_R K_R a_R + K_Br [Br-] + K_I [I-]
    theta_i = K_i a_i / Xi ,   theta_X = (K_Br [Br-] + K_I [I-]) / Xi
    Phi_E   = k_et,E / (k_des,E + sum_R k_reg,R a_R)
    g_E     = sum_R k_reg,R a_R / (k_des,E + sum_R k_reg,R a_R)
    kappa~  = sum_i k_et,i theta_i                                   (spec D1 2.4.4)
    lambda~ = sum_i lambda0_i theta_i exp(-|z_i| (u_gel sqrt(I_ref/I) + u_X theta_X))
    f_i     = k_et,i theta_i / kappa~                                (spec D1 2.4.6)

E runs over the electron-transfer agents (role ``eta``, and ``single``, which is an ETA with
no regenerator, spec D1 2.4.3), R over the regenerators. Every rate constant of
:class:`Kinetics` already contains its agent's Arrhenius factor (spec D1 2.4.5).

Sums over agents are explicit loops of elementwise operations in a fixed order, never BLAS
reductions, so a column's result does not depend on the batch it is computed in, on thread
count or on memory alignment (spec D1 3.6).
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from darkroom_atlas.chemistry.solution import WorkingSolution
from darkroom_atlas.chemistry.species import GAS_CONSTANT
from darkroom_atlas.data.models import (
    GelatinParams,
    KineticsFile,
    SolventParams,
    Species,
    TransportParams,
)

__all__ = [
    "ARRHENIUS_REFERENCE_K",
    "POSITIVE_PART_EPS",
    "ROLES",
    "Kinetics",
    "SurfaceState",
    "arrhenius_factor",
    "cluster_hazards",
    "positive_part",
    "resolve_kinetics",
    "silver_saturation",
    "surface_state",
]

ARRHENIUS_REFERENCE_K = 293.15
"""Reference temperature T_r of the Arrhenius and halide van 't Hoff laws, 20 C (spec D1
2.4.5). Not the 25 C reference of the pK values."""

POSITIVE_PART_EPS = 1e-9
"""Width eps_c of the smooth positive part, mol/L (spec D1 2.4.6, assumed)."""

ROLES = ("eta", "regenerator", "single")
"""Agent roles (spec D1 2.3.5, 5.2)."""

_CELSIUS_TO_KELVIN = 273.15
_SULFITE_REF = 0.1  # mol/L, C_sat reference sulfite (spec D1 2.4.6)
_BROMIDE_REF = 1e-3  # mol/L, C_sat reference bromide (spec D1 2.4.6)
_BROMIDE_FLOOR = 1e-5  # mol/L, C_sat bromide floor (spec D1 2.4.6)
_SULFITE_UNIT = 1.0  # mol/L, unit of [SO3 2-]_bulk in the attrition factor (spec D1 2.4.2)


def _readonly(a, dtype=np.float64) -> np.ndarray:
    a = np.array(a, dtype=dtype)
    a.flags.writeable = False
    return a


def positive_part(c: np.ndarray) -> np.ndarray:
    """Smooth positive part ``C+ = (C + sqrt(C^2 + eps_c^2)) / 2`` (spec D1 2.4.6).

    For ``C < 0`` the identical form ``eps_c^2 / (2 (sqrt(C^2 + eps_c^2) - C))`` is used, which
    has no cancellation, so ``C+ > 0`` for every finite ``C``.
    """
    c = np.asarray(c, dtype=np.float64)
    eps2 = POSITIVE_PART_EPS * POSITIVE_PART_EPS
    root = np.sqrt(c * c + eps2)
    with np.errstate(divide="ignore", invalid="ignore"):
        negative = 0.5 * eps2 / (root - c)
    return np.where(c >= 0.0, 0.5 * (c + root), negative)


def arrhenius_factor(activation_energy_kj_mol: float, temperature_k: float) -> float:
    """``A(T) = exp(-(Ea/R) (1/T - 1/T_r))``, ``T_r = 293.15 K`` (spec D1 2.4.5)."""
    t = float(temperature_k)
    if not (math.isfinite(t) and t > 0.0):
        raise ValueError(f"temperature must be a positive number of kelvin, got {temperature_k}")
    ea = float(activation_energy_kj_mol) * 1e3  # J/mol
    if not math.isfinite(ea):
        raise ValueError(f"activation energy must be finite, got {activation_energy_kj_mol}")
    return math.exp(-(ea / GAS_CONSTANT) * (1.0 / t - 1.0 / ARRHENIUS_REFERENCE_K))


@dataclass(frozen=True, eq=False)
class Kinetics:
    """Numeric constants resolved for one solution and temperature (spec D1 10 D1c).

    Per-agent arrays have shape ``(A,)`` in the order of ``agent_ids`` (component ids).
    ``charges`` is the signed charge of each agent's active form. ``k_et``, ``lambda0``,
    ``k_des`` (1/s) and ``k_reg`` (L/(mol s)) include the agent's Arrhenius factor;
    ``K_ads`` (L/mol) does not depend on temperature (spec D1 2.4.5). ``k_des`` is NaN for
    regenerators, ``k_reg`` NaN for ETAs and single agents. ``K_br`` and ``K_i`` (L/mol) are
    the halide adsorption constants at the solution temperature (``Delta H_Br``, spec D1
    2.4.5). ``u_gel`` is the gelatin barrier at I = 1 mol/L, i.e. the file value times
    ``sqrt(I_ref / (1 mol/L))``, so that the barrier is ``u_gel / sqrt(I)`` with I in mol/L
    (spec D1 2.4.2). Arrays are read-only.
    """

    agent_ids: tuple[str, ...]
    roles: tuple[str, ...]
    charges: np.ndarray
    k_et: np.ndarray
    K_ads: np.ndarray
    lambda0: np.ndarray
    k_des: np.ndarray  # NaN for regenerators
    k_reg: np.ndarray  # NaN for ETAs
    sulfonated: np.ndarray  # bool
    K_br: float  # at the solution temperature (Delta H_Br, spec D1 2.4.5)
    K_i: float  # same van 't Hoff factor
    u_gel: float
    u_halide: float
    u_li: float
    n_sat: int
    solvent: SolventParams
    transport: TransportParams
    gelatin: GelatinParams

    def __post_init__(self) -> None:
        a = len(self.agent_ids)
        object.__setattr__(self, "agent_ids", tuple(self.agent_ids))
        object.__setattr__(self, "roles", tuple(self.roles))
        if len(self.roles) != a:
            raise ValueError(f"roles needs {a} entries, got {len(self.roles)}")
        unknown = sorted(set(self.roles) - set(ROLES))
        if unknown:
            raise ValueError(f"unknown agent role(s) {unknown}; roles are {ROLES}")
        for name in ("charges", "k_et", "K_ads", "lambda0", "k_des", "k_reg"):
            arr = _readonly(getattr(self, name))
            if arr.shape != (a,):
                raise ValueError(f"{name} must have shape ({a},), got {arr.shape}")
            object.__setattr__(self, name, arr)
        sulfonated = _readonly(self.sulfonated, dtype=bool)
        if sulfonated.shape != (a,):
            raise ValueError(f"sulfonated must have shape ({a},), got {sulfonated.shape}")
        object.__setattr__(self, "sulfonated", sulfonated)
        for name in ("k_et", "K_ads", "lambda0"):
            arr = getattr(self, name)
            if not np.all(np.isfinite(arr) & (arr >= 0.0)):
                raise ValueError(f"{name} must be finite and >= 0")
        etas = self.eta_mask
        if not np.all(np.isfinite(self.k_des[etas]) & (self.k_des[etas] > 0.0)):
            raise ValueError("k_des must be finite and > 0 for every ETA and single agent")
        if not np.all(np.isfinite(self.k_reg[~etas]) & (self.k_reg[~etas] >= 0.0)):
            raise ValueError("k_reg must be finite and >= 0 for every regenerator")
        for name in ("K_br", "K_i", "u_gel", "u_halide", "u_li"):
            value = float(getattr(self, name))
            if not (math.isfinite(value) and value >= 0.0):
                raise ValueError(f"{name} must be finite and >= 0, got {value}")
            object.__setattr__(self, name, value)
        if isinstance(self.n_sat, bool) or int(self.n_sat) != self.n_sat or self.n_sat < 1:
            raise ValueError(f"n_sat must be an integer >= 1, got {self.n_sat!r}")
        object.__setattr__(self, "n_sat", int(self.n_sat))

    @property
    def eta_mask(self) -> np.ndarray:
        """True for agents that occupy sites when oxidised (``eta`` and ``single``)."""
        return np.array([r != "regenerator" for r in self.roles], dtype=bool)


def resolve_kinetics(
    kinetics: KineticsFile, solution: WorkingSolution, species: Mapping[str, Species]
) -> Kinetics:
    """Constants of ``kinetics`` for ``solution`` at its temperature (spec D1 2.4.3, 2.4.5).

    The agents are the components of ``solution.component_ids`` (in that order) that have a
    ``developing_agent`` role and a total > 0. Raises ``ValueError`` if an agent's kinetic
    class is not in ``kinetics.agents``, if an ETA or single agent has no ``k_des``, or if a
    regenerator has no ``k_reg``.
    """
    t_k = float(solution.temperature_c) + _CELSIUS_TO_KELVIN
    ids, roles, charges, sulfonated = [], [], [], []
    k_et, k_ads, lambda0, k_des, k_reg = [], [], [], [], []
    for cid, total in zip(solution.component_ids, solution.totals, strict=True):
        s = species.get(cid)
        if s is None:
            raise ValueError(f"component {cid!r} of the solution is not in species")
        agent = s.developing_agent
        if agent is None or not total > 0.0:
            continue
        entry = kinetics.agents.get(agent.kinetic_class)
        if entry is None:
            raise ValueError(
                f"{cid}: kinetic class {agent.kinetic_class!r} is not in kinetics.agents "
                f"({', '.join(sorted(kinetics.agents))})"
            )
        where = f"agents.{agent.kinetic_class}"
        factor = arrhenius_factor(float(entry.activation_energy.to("kJ/mol")), t_k)
        if agent.role == "regenerator":
            if entry.k_reg is None:
                raise ValueError(f"{cid}: regenerator needs {where}.k_reg")
            k_reg.append(float(entry.k_reg.to("L/(mol*s)")) * factor)
            k_des.append(math.nan)
        else:
            if entry.k_des is None:
                raise ValueError(f"{cid}: {agent.role} agent needs {where}.k_des")
            k_des.append(float(entry.k_des.to("1/s")) * factor)
            k_reg.append(math.nan)
        ids.append(cid)
        roles.append(agent.role)
        charges.append(float(s.acid_base.charge_most_protonated - agent.active_state))
        sulfonated.append(agent.product_sulfonated)
        k_et.append(float(entry.k_et.to("1/s")) * factor)
        k_ads.append(float(entry.K_ads.to("L/mol")))
        lambda0.append(float(entry.lambda0.to("1/s")) * factor)

    halides = kinetics.halides
    dh = float(halides.adsorption_enthalpy.to("kJ/mol"))
    # van 't Hoff for K_Br and K_I with the same factor (spec D1 2.4.5, 2.4.11).
    k_br = float(halides.K_Br.to("L/mol")) * arrhenius_factor(dh, t_k)
    k_i = k_br * float(halides.iodide_ratio.to("1"))
    i_ref = float(kinetics.barriers.ionic_strength_ref.to("mol/L"))
    if not i_ref > 0.0:
        raise ValueError(f"barriers.ionic_strength_ref must be > 0, got {i_ref}")
    n_sat = float(kinetics.latent.n_sat.to("1"))
    if n_sat != round(n_sat):
        raise ValueError(f"latent.n_sat must be an integer, got {n_sat}")
    return Kinetics(
        agent_ids=tuple(ids),
        roles=tuple(roles),
        charges=np.array(charges, dtype=np.float64),
        k_et=np.array(k_et, dtype=np.float64),
        K_ads=np.array(k_ads, dtype=np.float64),
        lambda0=np.array(lambda0, dtype=np.float64),
        k_des=np.array(k_des, dtype=np.float64),
        k_reg=np.array(k_reg, dtype=np.float64),
        sulfonated=np.array(sulfonated, dtype=bool),
        K_br=k_br,
        K_i=k_i,
        u_gel=float(kinetics.barriers.u_gel.to("1")) * math.sqrt(i_ref),
        u_halide=float(kinetics.barriers.u_halide.to("1")),
        u_li=float(kinetics.latent.u_li.to("1")),
        n_sat=round(n_sat),
        solvent=kinetics.solvent,
        transport=kinetics.transport,
        gelatin=kinetics.gelatin,
    )


@dataclass(frozen=True, eq=False)
class SurfaceState:
    """Adsorption state and rate factors of ``N`` cells (spec D1 2.4.3, 10 D1c).

    ``theta`` (A, N) and ``theta_halide`` (N,) are coverages; ``kappa_tilde`` and
    ``lambda_tilde`` (N,) are the growth and initiation rate factors in 1/s; ``partition``
    (A, N) is each agent's share ``f_i`` of the electron demand (0 where ``kappa_tilde`` is
    0); ``regenerated`` (A, N) is ``g_E``, NaN for regenerators. Arrays are read-only.
    """

    theta: np.ndarray  # (A, N)
    theta_halide: np.ndarray  # (N,)
    kappa_tilde: np.ndarray  # (N,) 1/s
    lambda_tilde: np.ndarray  # (N,) 1/s
    partition: np.ndarray  # (A, N)
    regenerated: np.ndarray  # (A, N), NaN for regenerators


def _check_concentration(name: str, value: np.ndarray) -> None:
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must be finite")
    if np.any(value < 0.0):
        raise ValueError(f"{name} must be >= 0 (apply positive_part first)")


def surface_state(
    active: np.ndarray,
    bromide: np.ndarray,
    iodide: np.ndarray,
    kin: Kinetics,
    ionic_strength: float,
) -> SurfaceState:
    """Surface state of ``N`` cells (spec D1 2.4.3).

    ``active`` (A, N) holds the active-form concentrations ``a_i`` (mol/L) in the order of
    ``kin.agent_ids``; ``bromide`` and ``iodide`` (N,) are mol/L. All are the rate-law
    concentrations: the caller applies :func:`positive_part` to the totals before the
    speciation (spec D1 2.4.6); negative or non-finite inputs raise ``ValueError``.
    ``ionic_strength`` (mol/L, > 0) screens the gelatin barrier (spec D1 2.4.2).
    """
    a_agents = len(kin.agent_ids)
    active = np.asarray(active, dtype=np.float64)
    bromide = np.asarray(bromide, dtype=np.float64)
    iodide = np.asarray(iodide, dtype=np.float64)
    if active.ndim != 2 or active.shape[0] != a_agents:
        raise ValueError(f"active must have shape ({a_agents}, N), got {active.shape}")
    n = active.shape[1]
    for name, arr in (("bromide", bromide), ("iodide", iodide)):
        if arr.shape != (n,):
            raise ValueError(f"{name} must have shape ({n},), got {arr.shape}")
    _check_concentration("active", active)
    _check_concentration("bromide", bromide)
    _check_concentration("iodide", iodide)
    i_strength = float(ionic_strength)
    if not (math.isfinite(i_strength) and i_strength > 0.0):
        raise ValueError(f"ionic_strength must be finite and > 0, got {ionic_strength}")

    etas = kin.eta_mask
    # Regeneration supply sum_R k_reg,R a_R (1/s), per cell.
    supply = np.zeros(n)
    for i in range(a_agents):
        if not etas[i]:
            supply = supply + kin.k_reg[i] * active[i]

    loaded = np.empty((a_agents, n))  # K_i a_i
    regenerated = np.full((a_agents, n), np.nan)
    xi = np.ones(n)
    for i in range(a_agents):
        loaded[i] = kin.K_ads[i] * active[i]
        if etas[i]:
            denom = kin.k_des[i] + supply
            phi = kin.k_et[i] / denom
            regenerated[i] = supply / denom
            xi = xi + loaded[i] * (1.0 + phi)
        else:
            xi = xi + loaded[i]
    halide = kin.K_br * bromide + kin.K_i * iodide
    xi = xi + halide

    theta = loaded / xi
    theta_halide = halide / xi
    barrier = kin.u_gel / math.sqrt(i_strength) + kin.u_halide * theta_halide
    kappa = np.zeros(n)
    lam = np.zeros(n)
    rate = np.empty((a_agents, n))
    for i in range(a_agents):
        rate[i] = kin.k_et[i] * theta[i]
        kappa = kappa + rate[i]
        lam = lam + kin.lambda0[i] * theta[i] * np.exp(-abs(kin.charges[i]) * barrier)
    positive = kappa > 0.0
    partition = np.where(positive, rate / np.where(positive, kappa, 1.0), 0.0)

    return SurfaceState(
        theta=_readonly(theta),
        theta_halide=_readonly(theta_halide),
        kappa_tilde=_readonly(kappa),
        lambda_tilde=_readonly(lam),
        partition=_readonly(partition),
        regenerated=_readonly(regenerated),
    )


def silver_saturation(
    free_sulfite: np.ndarray, bromide: np.ndarray, solvent: SolventParams
) -> np.ndarray:
    """Dissolved-silver saturation ``C_sat`` (mol/L) of spec D1 2.4.6::

        C_sat = C_sat,ref ([SO3 2-]+ / 0.1 M)^nu_s (1 mM / sqrt(([Br-]+)^2 + (1e-5 M)^2))^nu_b

    ``free_sulfite`` is the SO3 2- concentration and ``bromide`` the Br- concentration
    (mol/L, any matching shapes); the smooth positive part is applied here.
    """
    so3 = positive_part(free_sulfite)
    br = positive_part(bromide)
    c_ref = float(solvent.c_sat_ref.to("mol/L"))
    nu_s = float(solvent.sulfite_exponent.to("1"))
    nu_b = float(solvent.bromide_exponent.to("1"))
    floor = np.sqrt(br * br + _BROMIDE_FLOOR * _BROMIDE_FLOOR)
    return c_ref * (so3 / _SULFITE_REF) ** nu_s * (_BROMIDE_REF / floor) ** nu_b


def cluster_hazards(
    kin: Kinetics, fog_ratio: float, cluster_cap: int, free_sulfite_bulk: float
) -> np.ndarray:
    """Cluster factors ``h_n`` for ``n = 0 .. cluster_cap``, shape ``(cluster_cap + 1,)``
    aligned with the last axis of the latent table (spec D1 2.4.2)::

        h_0 = h_fog
        h_n = exp(-u_LI (n_sat - n)) exp(-zeta_s [SO3 2-]_bulk / (1 mol/L))   (n < n_sat)
        h_n = 1                                                               (n >= n_sat)

    ``free_sulfite_bulk`` is the initial free SO3 2- of the bath (mol/L). ``n = 1`` (no
    stable centre; its latent probability is 0) takes the ``n < n_sat`` formula. The top bin
    ``n = cluster_cap`` holds every larger cluster and takes the value at ``cluster_cap``.
    Read-only.
    """
    if isinstance(cluster_cap, bool) or not isinstance(cluster_cap, (int, np.integer)):
        raise ValueError(f"cluster_cap must be an integer >= 2, got {cluster_cap!r}")  # noqa: TRY004
    cap = int(cluster_cap)
    if cap < 2:
        raise ValueError(f"cluster_cap must be an integer >= 2, got {cluster_cap!r}")
    fog = float(fog_ratio)
    if not (math.isfinite(fog) and fog >= 0.0):
        raise ValueError(f"fog_ratio must be finite and >= 0, got {fog_ratio}")
    sulfite = float(free_sulfite_bulk)
    if not (math.isfinite(sulfite) and sulfite >= 0.0):
        raise ValueError(f"free_sulfite_bulk must be finite and >= 0, got {free_sulfite_bulk}")
    zeta = float(kin.solvent.latent_attrition.to("1"))
    attrition = math.exp(-zeta * sulfite / _SULFITE_UNIT)
    h = np.empty(cap + 1)
    h[0] = fog
    for size in range(1, cap + 1):
        if size >= kin.n_sat:
            h[size] = 1.0
        else:
            h[size] = math.exp(-kin.u_li * (kin.n_sat - size)) * attrition
    h.flags.writeable = False
    return h
