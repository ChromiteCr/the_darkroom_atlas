"""Working solutions: bath formula to component totals, TOTH, ionic strength and pH
(spec D1 2.3.1, 2.3.2, 5.5).

A bath lists ingredients in g/L of its base solution (stock or concentrate). With the
dilution fraction ``f = a/(a+b)`` of ``"a+b"`` (D0 ``parse_dilution``)::

    T_k   = f * sum_i (w_i / M_i) * s_i,k        [mol/L]
    TOTH  = f * sum_i (w_i / M_i) * tau_i        [mol/L]
    m_cat = f * sum_i (w_i / M_i) * (n_Na,i + n_K,i)

A component weighed directly counts as an ingredient with ``s = {self: 1}``, ``tau = 0`` and
no cations (spec D1 5.2). The ionic strength ``I = 1/2 sum_s z_s^2 [s]`` (cations, spectator
ions, every acid-base state, the charged components without a family and the water ions) is
found by fixed-point iteration from ``I0 = m_cat`` until ``|dI| < 1e-6 M`` (spec D1 2.3.2).
"""

import math
import numbers
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from darkroom_atlas.chemistry.speciation import Speciation, fractions, make_speciation, solve_ph
from darkroom_atlas.chemistry.speciation import water_ions as _water_ions
from darkroom_atlas.chemistry.species import BUILTIN_COMPONENTS, ComponentInfo, component_table
from darkroom_atlas.data.models import Bath, KineticsFile, Species
from darkroom_atlas.units import parse_dilution

__all__ = [
    "ALWAYS_COMPONENTS",
    "IONIC_STRENGTH_TOL",
    "WorkingSolution",
    "component_charges",
    "working_solution",
]

IONIC_STRENGTH_TOL = 1e-6
"""Convergence of the ionic-strength fixed point, mol/L (spec D1 2.3.2)."""

ALWAYS_COMPONENTS = ("sulfite", "bromide", "iodide", "silver")
"""Components every working solution carries, even at zero total: the products and the
solvent of development (bromide, iodide, dissolved silver, spec D1 2.4 source terms) and
sulfite, which the solvent and sulfonation laws read by name (spec D1 2.4.11)."""

_MAX_IONIC_ITER = 100
_CELSIUS_TO_KELVIN = 273.15
_CHARGE_TOL = 1e-9


def _readonly(a) -> np.ndarray:
    a = np.array(a, dtype=np.float64)
    a.flags.writeable = False
    return a


@dataclass(frozen=True, eq=False)
class WorkingSolution:
    """A developer (or other bath) at one dilution and temperature (spec D1 2.3.2, 10 D1c).

    ``component_ids`` are the components the bath's ingredients dissolve to plus
    ``ALWAYS_COMPONENTS``, in species-table order (the built-in ``toth`` and ``gelatin``
    excluded; about 8 for a developer, spec D1 3.10). Membership follows the formula, not
    the amounts, so an override of 0 g/L keeps the shapes. ``totals`` are their totals in
    mol/L. ``toth`` is the proton excess (mol/L)
    relative to the reference states, ``cation_molality`` is ``m_cat`` (Na+ plus K+, mol/L)
    and ``ionic_strength`` the converged ``I`` (mol/L). ``speciation`` holds the conditional
    constants of the families of ``component_ids`` that have one, in that order, at
    ``(T, I, m_cat)``; ``ph`` is the activity pH that solves the proton condition with them.
    """

    bath_id: str
    dilution: str
    dilution_fraction: float
    temperature_c: float
    component_ids: tuple[str, ...]
    totals: np.ndarray  # (K_c,) mol/L
    toth: float
    cation_molality: float
    ionic_strength: float
    speciation: Speciation
    ph: float

    def __post_init__(self) -> None:
        totals = _readonly(self.totals)
        if totals.shape != (len(self.component_ids),):
            raise ValueError(
                f"totals must have shape ({len(self.component_ids)},), got {totals.shape}"
            )
        object.__setattr__(self, "totals", totals)
        object.__setattr__(self, "component_ids", tuple(self.component_ids))

    @property
    def family_index(self) -> tuple[int, ...]:
        """Index into ``component_ids`` of each family of ``speciation``, in its order."""
        position = {cid: i for i, cid in enumerate(self.component_ids)}
        return tuple(position[f.component] for f in self.speciation.families)

    @property
    def family_totals(self) -> np.ndarray:
        """Totals of the speciation families, ``(F,)`` mol/L (the ``totals`` of solve_ph)."""
        return _readonly(self.totals[list(self.family_index)])

    def total(self, component: str) -> float:
        """Total of one component in mol/L; raises ``KeyError`` for an unknown id."""
        try:
            return float(self.totals[self.component_ids.index(component)])
        except ValueError:
            raise KeyError(
                f"{component!r} is not a component; components: {', '.join(self.component_ids)}"
            ) from None


def _reference_charge(info: ComponentInfo) -> int:
    family = info.family
    return family.charges[family.protons_rel.index(0)]


def component_charges(
    species: Mapping[str, Species], table: Mapping[str, ComponentInfo] | None = None
) -> dict[str, float]:
    """Charge of the reference state of every component (spec D1 2.3.1).

    Components with an acid-base family take it from the family. For the others (bromide,
    iodide, ...) it follows from the electroneutrality of the ingredients that dissolve to
    them: ``sum_k s_k z_k + tau + n_cat + sum(spectator charges) = 0``. Raises
    ``ValueError`` if two ingredients imply different charges. Components that no
    ingredient determines are left out.
    """
    table = component_table(species) if table is None else table
    charges: dict[str, float] = {
        cid: float(_reference_charge(info))
        for cid, info in table.items()
        if cid not in BUILTIN_COMPONENTS and info.family is not None
    }
    pending = [s for s in species.values() if s.kind == "ingredient" and s.dissolves_to]
    progress = True
    while pending and progress:
        progress = False
        rest = []
        for s in pending:
            unknown = [k for k in s.dissolves_to if k not in charges]
            if len(unknown) > 1:
                rest.append(s)
                continue
            known = sum(n * charges[k] for k, n in s.dissolves_to.items() if k in charges)
            balance = known + s.toth + sum(s.cations.values()) + sum(s.spectator_charge.values())
            if unknown:
                charges[unknown[0]] = -balance / s.dissolves_to[unknown[0]]
                progress = True
        pending = rest
    for s in species.values():
        if s.kind != "ingredient" or not s.dissolves_to:
            continue
        if not all(k in charges for k in s.dissolves_to):
            continue
        balance = sum(n * charges[k] for k, n in s.dissolves_to.items())
        balance += s.toth + sum(s.cations.values()) + sum(s.spectator_charge.values())
        if abs(balance) > _CHARGE_TOL:
            raise ValueError(
                f"ingredient {s.id!r} is not electroneutral with the component charges "
                f"{ {k: charges[k] for k in s.dissolves_to} } (net charge {balance:+.6g})"
            )
    return charges


_NO_OVERRIDE = object()


def _amount_mol_per_l(bath: Bath, i: int, species: Species, grams_per_l: object) -> float:
    if grams_per_l is not _NO_OVERRIDE:
        if isinstance(grams_per_l, bool) or not isinstance(grams_per_l, numbers.Real):
            raise ValueError(
                f"overrides[{species.id!r}]: amount must be a real number in g/L, "
                f"got {type(grams_per_l).__name__} {grams_per_l!r}"
            )
        w = float(grams_per_l)
        if not (math.isfinite(w) and w >= 0.0):
            raise ValueError(
                f"overrides[{species.id!r}]: amount must be a finite number >= 0 g/L, got {w}"
            )
        return w / float(species.molar_mass.to("g/mol"))
    amount = bath.components[i].amount
    if isinstance(amount.value, list):
        raise ValueError(f"bath {bath.id!r}: components[{i}].amount must be a single value")  # noqa: TRY004
    try:
        return float(amount.to("g/L")) / float(species.molar_mass.to("g/mol"))
    except ValueError:
        pass
    try:
        return float(amount.to("mol/L"))
    except ValueError:
        raise ValueError(
            f"bath {bath.id!r}: components[{i}].amount: unit {amount.unit!r} is neither a mass "
            "concentration (g/L) nor a concentration (mol/L)"
        ) from None


def working_solution(
    bath: Bath,
    species: Mapping[str, Species],
    kinetics: KineticsFile,
    *,
    dilution: str,
    temperature_c: float,
    overrides: Mapping[str, float] | None = None,
) -> WorkingSolution:
    """Working solution of ``bath`` at ``dilution`` and ``temperature_c`` (spec D1 2.3.2).

    ``dilution`` is relative to the bath's base solution (``"stock"``, ``"1+1"``,
    ``"1+31"``, ...). ``overrides`` maps ingredient ids of the bath to amounts in g/L of the
    base solution that replace the file's values (nuisance parameters, spec D1 2.3.2, 5.5).
    Conditional constants use ``kinetics.activity``; the class pK of developing agents and
    the gelatin pK come from ``kinetics`` through ``component_table``.
    """
    f = parse_dilution(dilution)
    t_c = float(temperature_c)
    if not math.isfinite(t_c) or t_c + _CELSIUS_TO_KELVIN <= 0.0:
        raise ValueError(f"temperature_c must be a finite temperature in C, got {temperature_c}")
    t_k = t_c + _CELSIUS_TO_KELVIN

    overrides = dict(overrides or {})
    names = [c.species for c in bath.components]
    unknown = sorted(set(overrides) - set(names))
    if unknown:
        raise ValueError(
            f"overrides: {', '.join(map(repr, unknown))} not an ingredient of bath "
            f"{bath.id!r}; ingredients: {', '.join(names)}"
        )

    table = component_table(species, kinetics=kinetics)
    charges = component_charges(species, table)

    ingredients = []
    members = set(ALWAYS_COMPONENTS)
    for i, name in enumerate(names):
        s = species.get(name)
        if s is None:
            raise ValueError(f"bath {bath.id!r}: components[{i}].species: unknown species {name!r}")
        if s.id in BUILTIN_COMPONENTS:
            raise ValueError(f"bath {bath.id!r}: {name!r} is a built-in pseudo-component")
        targets = (s.id,) if s.kind == "component" else tuple(s.dissolves_to)
        for target in targets:
            if target not in table or target in BUILTIN_COMPONENTS:
                raise ValueError(f"ingredient {s.id!r} dissolves to unknown component {target!r}")
        members.update(targets)
        ingredients.append(s)
    missing = sorted(members - set(table))
    if missing:
        raise ValueError(f"species table lacks the components {', '.join(map(repr, missing))}")
    component_ids = tuple(cid for cid in table if cid in members)
    index = {cid: k for k, cid in enumerate(component_ids)}

    totals = np.zeros(len(component_ids))
    toth = 0.0
    cations = 0.0
    spectator_i2 = 0.0  # sum of z^2 c over spectator ions, mol/L
    for i, (name, s) in enumerate(zip(names, ingredients, strict=True)):
        n = f * _amount_mol_per_l(bath, i, s, overrides.get(name, _NO_OVERRIDE))
        if s.kind == "component":
            if charges.get(s.id, 0.0) != 0.0:
                raise ValueError(
                    f"bath {bath.id!r}: component {s.id!r} is weighed directly but its "
                    f"reference state has charge {charges[s.id]:+g}; weigh a salt instead"
                )
            totals[index[s.id]] += n
            continue
        for target, stoich in s.dissolves_to.items():
            totals[index[target]] += n * stoich
        toth += n * s.toth
        cations += n * sum(s.cations.values())
        spectator_i2 += n * sum(q * q for q in s.spectator_charge.values())

    for cid in component_ids:
        if totals[index[cid]] > 0.0 and table[cid].family is None and cid not in charges:
            raise ValueError(
                f"component {cid!r} has no acid-base family and no ingredient fixes its charge"
            )

    family_ids = tuple(cid for cid in component_ids if table[cid].family is not None)
    families = tuple(table[cid].family for cid in family_ids)
    family_totals = totals[[index[cid] for cid in family_ids]]
    # Charged components without a family (Br-, I-, ...): z^2 c, mol/L.
    plain_i2 = sum(
        charges.get(cid, 0.0) ** 2 * totals[index[cid]]
        for cid in component_ids
        if table[cid].family is None
    )
    z2 = [np.asarray(fam.charges, dtype=np.float64) ** 2 for fam in families]
    params = kinetics.activity

    def solve(i_strength: float, ph0: float | None) -> tuple[Speciation, float]:
        spec = make_speciation(families, t_k, i_strength, cations, params)
        start = None if ph0 is None else np.asarray(ph0)
        return spec, float(solve_ph(spec, family_totals, np.asarray(toth), ph0=start))

    def ionic_strength(spec: Speciation, ph: float) -> float:
        alpha = fractions(spec, np.asarray(ph))
        h, oh = _water_ions(spec, np.asarray(ph))
        acc = cations + spectator_i2 + plain_i2 + float(h) + float(oh)
        for fi, w in enumerate(z2):
            acc += family_totals[fi] * float(np.dot(w, alpha[fi, : w.size]))
        return 0.5 * acc

    i_strength = cations
    spec, ph = solve(i_strength, None)
    for _ in range(_MAX_IONIC_ITER):
        i_new = ionic_strength(spec, ph)
        converged = abs(i_new - i_strength) < IONIC_STRENGTH_TOL
        i_strength = i_new
        spec, ph = solve(i_strength, ph)
        if converged:
            break
    else:
        raise RuntimeError(
            f"working_solution({bath.id!r}, {dilution!r}): ionic strength did not converge "
            f"to {IONIC_STRENGTH_TOL} mol/L in {_MAX_IONIC_ITER} iterations"
        )

    return WorkingSolution(
        bath_id=bath.id,
        dilution=dilution,
        dilution_fraction=f,
        temperature_c=t_c,
        component_ids=component_ids,
        totals=totals,
        toth=float(toth),
        cation_molality=float(cations),
        ionic_strength=float(i_strength),
        speciation=spec,
        ph=ph,
    )
