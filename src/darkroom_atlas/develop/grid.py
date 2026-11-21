"""Finite-volume column grid for development (spec D1 2.4.1, 3.1).

A column is, from bottom to top: ``n_z`` swollen emulsion cells (z = n_z ... 1), one swollen
overcoat cell and ``n_l`` liquid boundary-layer cells spanning ``delta_max``. The liquid is two
geometric sections, ``n_a = n_l // 2`` cells in ``[0, delta_agit]`` and ``n_l - n_a`` in
``[delta_agit, delta_max]``, so ``delta_agit`` is a cell face (the agitated Dirichlet face,
spec D1 2.4.8). Above that is the well-mixed tank (top Dirichlet value).
"""

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np

from darkroom_atlas._frozen import readonly_array

__all__ = ["PROFILE_LEVELS", "ColumnGrid", "column_grid"]

Kind = Literal["emulsion", "overcoat", "liquid"]

# Spec D1 3.1: profile -> (emulsion cells n_z, liquid cells n_l, liquid geometric ratio r_l).
PROFILE_LEVELS: dict[str, tuple[int, int, float]] = {
    "reference": (8, 16, 1.25),
    "calibration": (4, 8, 1.4),
    "image": (2, 4, 2.0),
}


@dataclass(frozen=True)
class ColumnGrid:
    """Cell layout of one development column, ordered bottom (base) to top (tank).

    ``thickness_um`` and ``height_um`` (cell centres above the film base) have shape ``(L,)``
    with ``L = n_z + 1 + n_l``. Heights in the liquid are measured from the base; the liquid
    heights of spec D1 0.3 / 2.4.8 (``delta_agit``, ``delta_max``) are measured from the
    overcoat's outer surface, i.e. ``height_um - gel_thickness_um``. ``n_a`` is the number of
    liquid cells below ``delta_agit``.
    """

    n_z: int
    n_l: int
    thickness_um: np.ndarray
    kinds: tuple[Kind, ...]
    height_um: np.ndarray
    delta_max_um: float
    delta_agit_um: float
    n_a: int

    @property
    def agit_cell(self) -> int:
        """Index (into the ``(L,)`` arrays) of the cell whose top face is at ``delta_agit``."""
        return self.n_z + self.n_a

    @property
    def gel_thickness_um(self) -> float:
        """Swollen emulsion plus overcoat thickness (the liquid starts here)."""
        return float(self.thickness_um[: self.n_z + 1].sum())


def column_grid(
    emulsion_dry_um: float,
    overcoat_dry_um: float,
    swelling_ratio: float,
    delta_agit_um: float,
    delta_max_um: float,
    profile: Literal["image", "calibration", "reference"],
) -> ColumnGrid:
    """Build the column grid of ``profile`` (spec D1 2.4.1, 3.1).

    Emulsion cells are ``swelling_ratio * emulsion_dry_um / n_z`` thick, the overcoat cell
    ``swelling_ratio * overcoat_dry_um``. The liquid section ``[a, b]`` with ``n`` cells has
    cells ``h1 * r_l**k`` with ``h1 = (b - a)(r_l - 1) / (r_l**n - 1)``; the sections are
    ``[0, delta_agit]`` (``n_l // 2`` cells) and ``[delta_agit, delta_max]``. ``TypeError`` on
    a non-number; ``ValueError`` on a non-positive or non-finite thickness or ratio,
    ``swelling_ratio < 1``, ``delta_agit_um >= delta_max_um`` or an unknown profile.
    """
    if profile not in PROFILE_LEVELS:
        raise ValueError(f"profile must be one of {sorted(PROFILE_LEVELS)}, got {profile!r}")
    _positive("emulsion_dry_um", emulsion_dry_um)
    # A zero-thickness overcoat cell would break the harmonic-mean face D (spec D1 2.4.8).
    _positive("overcoat_dry_um", overcoat_dry_um)
    _positive("delta_agit_um", delta_agit_um)
    _positive("delta_max_um", delta_max_um)
    if delta_agit_um >= delta_max_um:
        raise ValueError(
            f"delta_agit_um ({delta_agit_um!r}) must be < delta_max_um ({delta_max_um!r})"
        )
    _positive("swelling_ratio", swelling_ratio)
    if swelling_ratio < 1.0:
        raise ValueError(f"swelling_ratio must be >= 1, got {swelling_ratio!r}")

    n_z, n_l, r_l = PROFILE_LEVELS[profile]
    emulsion = np.full(n_z, swelling_ratio * emulsion_dry_um / n_z)
    overcoat = np.array([swelling_ratio * overcoat_dry_um])
    # Spec D1 2.4.1 / 3.1: two geometric liquid sections with a face exactly at delta_agit.
    n_a = n_l // 2
    lower = _geometric(float(delta_agit_um), n_a, r_l)
    upper = _geometric(float(delta_max_um) - float(delta_agit_um), n_l - n_a, r_l)
    liquid = np.concatenate([lower, upper])

    thickness = np.concatenate([emulsion, overcoat, liquid])
    tops = np.cumsum(thickness)
    height = tops - 0.5 * thickness
    kinds: tuple[Kind, ...] = ("emulsion",) * n_z + ("overcoat",) + ("liquid",) * n_l
    return ColumnGrid(
        n_z=n_z,
        n_l=n_l,
        thickness_um=readonly_array(thickness),
        kinds=kinds,
        height_um=readonly_array(height),
        delta_max_um=float(delta_max_um),
        delta_agit_um=float(delta_agit_um),
        n_a=n_a,
    )


def _geometric(length: float, n: int, ratio: float) -> np.ndarray:
    """``n`` cells growing by ``ratio`` that sum to ``length`` (last cell absorbs rounding)."""
    h1 = length * (ratio - 1.0) / (ratio**n - 1.0)
    cells = h1 * ratio ** np.arange(n, dtype=np.float64)
    cells[-1] = length - cells[:-1].sum()
    return cells


def _positive(name: str, value: float) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.floating, np.integer)):
        raise TypeError(f"{name} must be a real number, got {value!r}")
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive and finite, got {value!r}")
