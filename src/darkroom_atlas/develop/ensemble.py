"""Development engine: film column model and the ensemble integrator (spec D1 2.4, 3.1-3.5, 3.7).

A *cell* is one depth column (spec D1 2.4.1): ``n_z`` swollen emulsion levels, one overcoat
level and ``n_l`` liquid levels, ordered bottom (base) to top (tank) as in
:class:`~darkroom_atlas.develop.grid.ColumnGrid`. Its state is the totals ``C_k`` of the
transported components and ``TOTH`` on every level, and on every emulsion level the clocks
``Lambda_z`` (initiation), ``D_z`` (dissolved halide) and ``K_z`` (integral of kappa-tilde,
for :class:`CellHistories`), plus the cluster kernels ``G, H, Q`` per size class and cluster
size (spec D1 2.4.4). Members (strip steps or pixels) never enter the time integration: each
cell is driven by ``p-bar``, the area-weighted mean latent distribution of its members, and
the members' ``X, Y, P, F`` are synthesised from the kernels at the end and at snapshots.

One time step is the Strang split ``L(h/2) C(h) L(h/2)`` (spec D1 3.2): ``L`` the exact
lateral diffusion of :mod:`darkroom_atlas.develop.lateral`, ``C`` one ROS2 step of the
coupled vertical transport, reactions and kernels (spec D1 3.3) with the kernels eliminated by
a Schur complement and the column solved by the block-Thomas solver of
:mod:`darkroom_atlas.develop.rosenbrock`.

Conservation. ROS2 with an approximate Jacobian keeps a linear invariant ``w`` exactly when
``w^T J = 0``. The reaction part of ``J`` is therefore assembled as stoichiometric
coefficient times the derivative of three rates (``R_dev``, ``R_diss``, ``R_pd``), with the
same coefficients as ``f`` (spec D1 2.4.6), and the kernel rows are the exact derivatives of
which ``R_dev`` and ``R_pd`` are the ``p-bar`` weighted sums; the transport part is in flux
form. The flux booked for the tank is the ROS2-consistent flux through the Dirichlet face
(an accumulator row eliminated after the solve), so silver, halide, developer electrons and
TOTH balance to rounding.

Tank and the held upper liquid section (spec D1 2.4.8, 2.4.9). While agitated, the liquid
cells above ``delta_agit`` are held at the bulk value and belong to the tank. To make the
balance exact this module treats the held section as part of the tank volume: the
agitation-start jump is ``coef dI / (1 + coef sum_up Delta)`` and during agitation the bulk
moves by ``coef h J / (1 + coef sum_up Delta)``; ``flux_mean`` books the face flux, the
mixing impulse ``-dI/h`` and the held section's change ``sum_up Delta (c_{n+1} - c_n) / h``, so
:func:`darkroom_atlas.develop.tank.integrate_bulk` with the plain coefficient reproduces the
bulk trajectory bit for bit.

Per-block code uses no BLAS and reduces only along the last axis of contiguous arrays or by
explicit loops, so a cell's result does not depend on the batch it runs in (spec D1 3.6).
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.special import expit

from darkroom_atlas._frozen import FrozenDict, readonly_array
from darkroom_atlas.chemistry.kinetics import Kinetics, cluster_hazards, positive_part
from darkroom_atlas.chemistry.solution import WorkingSolution
from darkroom_atlas.chemistry.speciation import Speciation, fractions, solve_ph
from darkroom_atlas.chemistry.species import AcidBaseFamily, component_table
from darkroom_atlas.data.models import (
    BaseLayer,
    EmulsionLayer,
    Material,
    OvercoatLayer,
    Species,
)
from darkroom_atlas.develop import lateral as lat
from darkroom_atlas.develop import rosenbrock as ros
from darkroom_atlas.develop import tank as tk
from darkroom_atlas.develop.grid import ColumnGrid, column_grid
from darkroom_atlas.develop.schedule import Schedule
from darkroom_atlas.emulsion.classes import M_AG, SizeClasses, size_classes
from darkroom_atlas.latent.chain import ChainParams
from darkroom_atlas.latent.table import LatentTable, latent_table

__all__ = [
    "CHEM_CELL_UM",
    "PROFILE_DISCRETISATION",
    "BulkTrajectory",
    "CellHistories",
    "EnsembleResult",
    "FilmColumnModel",
    "Lateral",
    "Members",
    "TankSpec",
    "film_column_model",
    "run_ensemble",
]

CHEM_CELL_UM = 24.0
"""Chemistry cell edge of the image profile, um (spec D1 3.1, 3.8)."""

PROFILE_DISCRETISATION: dict[str, tuple[int, int, int]] = {
    "reference": (12, 21, 401),
    "calibration": (6, 9, 201),
    "image": (6, 9, 201),
}
"""Spec D1 3.1: profile -> (size classes K, Gauss-Hermite nodes Q_GH, latent table nodes)."""

_ELL_REF_UM = 1.0  # l_ref of the specific-surface growth scaling (spec D1 2.4.4, assumed)
_SOFTPLUS_WIDTH = 1e-3  # softplus width of max(C_sat - C_Ag, 0), times C_sat (spec D1 2.4.6)
_SULFITE_REF = 0.1  # mol/L, C_sat sulfite reference (spec D1 2.4.6)
_BROMIDE_REF = 1e-3  # mol/L, C_sat bromide reference (spec D1 2.4.6)
_BROMIDE_FLOOR = 1e-5  # mol/L, C_sat bromide floor (spec D1 2.4.6)
_WATER_ION_GEL_RATIO = 0.5  # r_gel of H+ and OH- (spec D1 2.3.7, salt value, assumed)
_T_REF_K = 293.15  # reference temperature of D_w (spec D1 2.3.7)
_CELSIUS = 273.15
_HALIDE_G_M2 = M_AG * 1e-15 * 1e12  # mol L^-1 um -> g/m2 of Ag, 0.10787 (spec D1 2.4.10)
_LN10 = math.log(10.0)
_ADAPTIVE_H0_S = 0.01  # first trial step of the adaptive profile after a reset point
_HISTORY_DEFAULT_SPACING_UM = CHEM_CELL_UM

_LATENT_CACHE: dict[tuple, LatentTable] = {}
_LATENT_CACHE_MAX = 16


# ---------------------------------------------------------------------------------------
# Data containers (spec D1 10, D1c develop/ensemble.py)
# ---------------------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class FilmColumnModel:
    """Everything a development column needs that does not depend on the exposure pattern.

    ``hazards`` ``(cap+1,)`` are the cluster factors ``h_n`` (entry 1 unused, spec D1 2.4.2);
    ``silver_per_cell`` ``(K,)`` is ``C_c,z`` in mol/L (equal for all z); ``gelatin_total``
    ``(L,)`` is ``T_gel`` per level in mol/L (0 in the liquid). The remaining fields extend
    spec D1 10: ``component_ids`` (the solution's transported components),
    ``diffusion_um2_s`` ``(K_c,)`` their liquid ``D(T, I)`` (spec D1 2.3.7), ``gel_ratio``
    ``(K_c,)`` ``D_gel / D_w``, ``water_ion_diffusion_um2_s`` ``(D_H, D_OH)`` in the liquid
    and ``profile``.
    """

    classes: SizeClasses
    latent: LatentTable
    hazards: np.ndarray
    silver_per_cell: np.ndarray
    gelatin_total: np.ndarray
    fog_ratio: float
    grid: ColumnGrid
    component_ids: tuple[str, ...]
    diffusion_um2_s: np.ndarray
    gel_ratio: np.ndarray
    water_ion_diffusion_um2_s: tuple[float, float]
    profile: str


@dataclass(frozen=True, eq=False)
class Members:
    """Members of a run: their exposure descriptors and how they map onto the cells.

    ``"columns"``: member ``i`` is cell ``i``, arrays ``(n_cells, n_z)``. ``"raster"``: a
    padded pixel window ``(h, w, n_z)`` over the cell grid of ``Lateral.shape_cells``, pixel
    pitch ``pixel_um``, ``origin_um`` the corner of pixel ``[0, 0]`` relative to the corner of
    cell ``[0, 0]``. Depth index ``j`` is ``z = j + 1`` (top first, spec D1 0.3).
    """

    exposure_volume: np.ndarray
    exposure_surface: np.ndarray
    layout: Literal["columns", "raster"]
    pixel_um: float | None = None
    origin_um: tuple[float, float] = (0.0, 0.0)


@dataclass(frozen=True)
class Lateral:
    """Lateral coupling of the cells (spec D1 3.5): grid shape, spacing and boundary."""

    shape_cells: tuple[int, int]
    spacing_um: float
    boundary: Literal["neumann", "periodic", "pair"]
    wavenumber_per_um: np.ndarray | None = None


@dataclass(frozen=True)
class TankSpec:
    """Tank volume (L) and the film area it serves (um^2) (spec D1 2.4.9)."""

    volume_l: float
    film_area_um2: float


@dataclass(frozen=True, eq=False)
class BulkTrajectory:
    """Bulk composition on the schedule: ``times_s`` (S+1,), ``totals`` (S+1, K_c), ``toth``
    (S+1,)."""

    times_s: np.ndarray
    totals: np.ndarray
    toth: np.ndarray


@dataclass(frozen=True, eq=False)
class CellHistories:
    """Monotone cumulative cell quantities on R recorded steps (spec D1 2.7.2).

    ``lambda_clock`` and ``kappa_integral`` ``(R, n_cells, n_z)`` (depth index ``z - 1``).
    """

    times_s: np.ndarray
    lambda_clock: np.ndarray
    kappa_integral: np.ndarray
    cell_spacing_um: float


@dataclass(frozen=True, eq=False)
class EnsembleResult:
    """Per-member development outputs (spec D1 2.4.10) and the run's tank bookkeeping."""

    silver: np.ndarray
    initiated: np.ndarray
    silver_m2: np.ndarray
    silver_physical: np.ndarray
    halide: np.ndarray
    growth_rate: np.ndarray
    snapshots: Mapping[float, "EnsembleResult"]
    flux_mean: np.ndarray
    bulk: BulkTrajectory
    min_scaled_concentration: float
    cell_history: CellHistories | None


# ---------------------------------------------------------------------------------------
# Film column model
# ---------------------------------------------------------------------------------------


def _scalar(quantity, unit: str, name: str) -> float:
    value = quantity.to(unit)
    if isinstance(value, np.ndarray):
        raise ValueError(f"{name} must be a single value, got a list")  # noqa: TRY004
    return float(value)


def _viscosity(table, temperature_k: float) -> float:
    """Water viscosity (cP): ``ln eta`` linear in ``1/T`` through the nodes, extrapolated at
    both ends (spec D1 2.3.7)."""
    t = np.asarray(table.temperature.to("C"), dtype=np.float64) + _CELSIUS
    ln_eta = np.log(np.asarray(table.viscosity.to("cP"), dtype=np.float64))
    x = 1.0 / t[::-1]  # increasing
    y = ln_eta[::-1]
    xi = 1.0 / temperature_k
    i = int(np.clip(np.searchsorted(x, xi) - 1, 0, x.size - 2))
    return math.exp(y[i] + (y[i + 1] - y[i]) * (xi - x[i]) / (x[i + 1] - x[i]))


def diffusion_factor(transport, temperature_c: float, ionic_strength: float) -> float:
    """``D(T, I) / D_w(T_r)`` = ``(T/T_r) (eta(T_r)/eta(T)) / (1 + b_eta I)`` (spec D1 2.3.7)."""
    t_k = float(temperature_c) + _CELSIUS
    b_eta = _scalar(transport.viscosity_ionic_coefficient, "L/mol", "viscosity_ionic_coefficient")
    eta_ratio = _viscosity(transport.water_viscosity, _T_REF_K) / _viscosity(
        transport.water_viscosity, t_k
    )
    return (t_k / _T_REF_K) * eta_ratio / (1.0 + b_eta * float(ionic_strength))


def _layers(material: Material) -> tuple[OvercoatLayer, EmulsionLayer, BaseLayer]:
    overcoats = [x for x in material.layers if isinstance(x, OvercoatLayer)]
    emulsions = [x for x in material.layers if isinstance(x, EmulsionLayer)]
    bases = [x for x in material.layers if isinstance(x, BaseLayer)]
    if len(overcoats) != 1 or len(emulsions) != 1 or len(bases) != 1:
        raise ValueError(
            f"material {material.id!r}: the bw-mechanistic column needs exactly one overcoat, "
            "one emulsion and one base layer (spec D1 2.0.1)"
        )
    return overcoats[0], emulsions[0], bases[0]


def _liquid_section(length: float, n: int, ratio: float) -> np.ndarray:
    h1 = length * (ratio - 1.0) / (ratio**n - 1.0)
    cells = h1 * ratio ** np.arange(n, dtype=np.float64)
    cells[-1] = length - cells[:-1].sum()
    return cells


def _custom_grid(
    emulsion_dry_um: float,
    overcoat_dry_um: float,
    swelling_ratio: float,
    delta_agit_um: float,
    delta_max_um: float,
    levels: tuple[int, int, float],
) -> ColumnGrid:
    """A :class:`ColumnGrid` with ``(n_z, n_l, r_l)`` other than a profile's (convergence
    tests, spec D1 3.9); same construction as :func:`column_grid`."""
    n_z, n_l, r_l = int(levels[0]), int(levels[1]), float(levels[2])
    if n_z < 1 or n_l < 2 or not r_l > 1.0:
        raise ValueError(f"levels must be (n_z >= 1, n_l >= 2, r_l > 1), got {levels!r}")
    base = column_grid(
        emulsion_dry_um, overcoat_dry_um, swelling_ratio, delta_agit_um, delta_max_um, "image"
    )
    emulsion = np.full(n_z, swelling_ratio * emulsion_dry_um / n_z)
    overcoat = np.array([swelling_ratio * overcoat_dry_um])
    n_a = n_l // 2
    liquid = np.concatenate(
        [
            _liquid_section(delta_agit_um, n_a, r_l),
            _liquid_section(delta_max_um - delta_agit_um, n_l - n_a, r_l),
        ]
    )
    thickness = np.concatenate([emulsion, overcoat, liquid])
    height = np.cumsum(thickness) - 0.5 * thickness
    return ColumnGrid(
        n_z=n_z,
        n_l=n_l,
        thickness_um=readonly_array(thickness),
        kinds=("emulsion",) * n_z + ("overcoat",) + ("liquid",) * n_l,
        height_um=readonly_array(height),
        delta_max_um=base.delta_max_um,
        delta_agit_um=base.delta_agit_um,
        n_a=n_a,
    )


def _sulfite_state(solution: WorkingSolution) -> tuple[int, int]:
    """(family index, state index) of free SO3 2- (the most deprotonated sulfite state)."""
    for f, fam in enumerate(solution.speciation.families):
        if fam.component == "sulfite":
            return f, fam.n_states - 1
    raise ValueError("the working solution has no sulfite family (spec D1 2.3.1)")


def _free_sulfite_bath(solution: WorkingSolution) -> float:
    f, s = _sulfite_state(solution)
    alpha = fractions(solution.speciation, np.asarray([solution.ph]))
    return float(alpha[f, s, 0] * solution.total("sulfite"))


def film_column_model(
    material: Material,
    film_format: str,
    exposure_time_s: float,
    kin: Kinetics,
    solution: WorkingSolution,
    profile: str,
    *,
    species: Mapping[str, Species] | None = None,
    n_classes: int | None = None,
    quadrature: int | None = None,
    latent_nodes: int | None = None,
    levels: tuple[int, int, float] | None = None,
) -> FilmColumnModel:
    """Build the column model of ``material`` for one solution and exposure time (spec D1
    2.0.2, 2.2.5, 2.4.1, 2.4.2, 2.4.7, 3.1).

    ``species`` supplies the components' diffusion data (default: the packaged library). The
    keyword overrides change the discretisation of ``profile`` (convergence tests, spec D1
    3.9): ``n_classes`` K, ``quadrature`` Q_GH, ``latent_nodes`` and ``levels`` ``(n_z, n_l,
    r_l)``. Latent tables are cached per material, exposure time and discretisation.
    ``ValueError`` on an unknown profile or format, a layer structure other than overcoat /
    emulsion / base, or missing emulsion parameters.
    """
    if profile not in PROFILE_DISCRETISATION:
        raise ValueError(
            f"profile must be one of {sorted(PROFILE_DISCRETISATION)}, got {profile!r}"
        )
    t_exp = float(exposure_time_s)
    if not (math.isfinite(t_exp) and t_exp > 0.0):
        raise ValueError(f"exposure_time_s must be finite and > 0, got {exposure_time_s!r}")
    overcoat, emulsion, base = _layers(material)
    if film_format not in base.formats:
        raise ValueError(
            f"material {material.id!r} has no format {film_format!r}; formats: "
            f"{', '.join(sorted(base.formats))}"
        )
    if (
        emulsion.swelling_ratio is None
        or emulsion.development is None
        or overcoat.thickness is None
    ):
        raise ValueError(
            f"material {material.id!r}: swelling_ratio, development and the overcoat thickness "
            "are required (spec D1 2.4.11)"
        )
    k_def, q_def, nodes_def = PROFILE_DISCRETISATION[profile]
    k = k_def if n_classes is None else n_classes
    q = q_def if quadrature is None else quadrature
    nodes = nodes_def if latent_nodes is None else latent_nodes
    classes = size_classes(emulsion, k)

    chain = ChainParams.from_layer(emulsion)
    key = (material.model_dump_json(), t_exp, k, q, nodes)
    table = _LATENT_CACHE.get(key)
    if table is None:
        table = latent_table(classes, chain, t_exp, nodes=nodes, quadrature=q)
        if len(_LATENT_CACHE) >= _LATENT_CACHE_MAX:
            _LATENT_CACHE.pop(next(iter(_LATENT_CACHE)))
        _LATENT_CACHE[key] = table

    swelling = _scalar(emulsion.swelling_ratio, "1", "swelling_ratio")
    overcoat_um = _scalar(overcoat.thickness, "um", "overcoat thickness")
    delta_agit = _scalar(kin.transport.delta_agit, "um", "delta_agit")
    delta_max = _scalar(kin.transport.delta_max, "um", "delta_max")
    if levels is None:
        grid = column_grid(
            classes.dry_thickness_um, overcoat_um, swelling, delta_agit, delta_max, profile
        )
    else:
        grid = _custom_grid(
            classes.dry_thickness_um, overcoat_um, swelling, delta_agit, delta_max, levels
        )
    n_z = grid.n_z
    delta_z = float(grid.thickness_um[0])

    # Spec D1 2.0.2: C_c,z = 1e15 n_c,z m_Ag,c / 107.87 with n_c,z = nu_c / (n_z Delta_z).
    silver = 1e15 * classes.number_per_um2 * classes.silver_g_per_grain / (M_AG * n_z * delta_z)

    # Spec D1 2.4.7: c_gel = 1e3 m_gel / (alpha_sw L_e) g/L (overcoat: rho_gel L_oc in place of
    # m_gel over alpha_sw L_oc), T_gel = c_gel g_OH.
    g_oh = _scalar(kin.gelatin.groups_per_gram, "mol/g", "groups_per_gram")
    rho_gel = _scalar(kin.gelatin.dry_density, "g/cm3", "dry_density")
    m_gel = _scalar(emulsion.gelatin_coverage, "g/m2", "gelatin_coverage")
    gel = np.zeros(grid.thickness_um.size)
    gel[:n_z] = 1e3 * m_gel / (swelling * classes.dry_thickness_um) * g_oh
    gel[n_z] = 1e3 * rho_gel / swelling * g_oh

    fog = _scalar(emulsion.development.fog_hazard_ratio, "1", "fog_hazard_ratio")
    cap = int(chain.cluster_cap)
    hazards = cluster_hazards(kin, fog, cap, _free_sulfite_bath(solution))

    if species is None:
        from darkroom_atlas.data.library import Library  # default data only

        species = Library().species()
    table_c = component_table(species)
    factor = diffusion_factor(kin.transport, solution.temperature_c, solution.ionic_strength)
    d_liq, ratio = [], []
    for cid in solution.component_ids:
        info = table_c.get(cid)
        if info is None:
            raise ValueError(f"component {cid!r} of the solution is not in the species table")
        d_liq.append(info.d_water_um2_s * factor)
        ratio.append(info.gel_ratio)
    d_h = _scalar(kin.transport.proton_diffusion, "um2/s", "proton_diffusion") * factor
    d_oh = _scalar(kin.transport.hydroxide_diffusion, "um2/s", "hydroxide_diffusion") * factor
    return FilmColumnModel(
        classes=classes,
        latent=table,
        hazards=hazards,
        silver_per_cell=readonly_array(silver),
        gelatin_total=readonly_array(gel),
        fog_ratio=fog,
        grid=grid,
        component_ids=tuple(solution.component_ids),
        diffusion_um2_s=readonly_array(d_liq),
        gel_ratio=readonly_array(ratio),
        water_ion_diffusion_um2_s=(float(d_h), float(d_oh)),
        profile=profile,
    )


# ---------------------------------------------------------------------------------------
# Engine set-up
# ---------------------------------------------------------------------------------------


def _lsum(a: np.ndarray) -> np.ndarray:
    """Sum over the last axis of a contiguous copy (batch-size independent, spec D1 3.6)."""
    return np.add.reduce(np.ascontiguousarray(a), axis=-1)


class _Setup:
    """Constants of one run, resolved once (no pydantic access inside ``f``)."""

    def __init__(
        self, film: FilmColumnModel, solution: WorkingSolution, kin: Kinetics, presoaked: bool
    ) -> None:
        if tuple(film.component_ids) != tuple(solution.component_ids):
            raise ValueError(
                "film column model and solution disagree on the components: "
                f"{film.component_ids} vs {solution.component_ids}"
            )
        grid = film.grid
        self.grid = grid
        self.n_z, self.n_l = grid.n_z, grid.n_l
        self.L = grid.thickness_um.size
        self.ids = tuple(solution.component_ids)
        self.kc = len(self.ids)
        self.nc = self.kc + 1  # components then TOTH
        self.it = self.kc  # TOTH index
        self.b = self.nc + 3  # dense block: components, TOTH, Lambda, D, K
        self.iL, self.iD, self.iK = self.nc, self.nc + 1, self.nc + 2
        self.layout = ros.ColumnLayout(
            n_dense=self.n_z, n_diag=self.L - self.n_z, n_comp=self.nc, n_local=3
        )
        self.delta = np.asarray(grid.thickness_um, dtype=np.float64)
        self.hf = 0.5 * (self.delta[:-1] + self.delta[1:])
        self.gel = np.array([kind != "liquid" for kind in grid.kinds])
        self.agit_cell = grid.agit_cell
        self.n_up = self.L - 1 - self.agit_cell
        self.up_thickness = self.delta[self.agit_cell + 1 :].copy()
        self.sum_up = float(self.up_thickness.sum())
        self.d_liq = np.asarray(film.diffusion_um2_s, dtype=np.float64)
        self.r_gel = np.asarray(film.gel_ratio, dtype=np.float64)
        self.d_h, self.d_oh = film.water_ion_diffusion_um2_s
        self.presoaked = bool(presoaked)
        self.half_time = _scalar(kin.transport.swelling_half_time, "s", "swelling_half_time")
        self.solution = solution
        self.kin = kin

        idx = {cid: i for i, cid in enumerate(self.ids)}
        for need in ("sulfite", "bromide", "iodide", "silver"):
            if need not in idx:
                raise ValueError(f"the solution lacks the component {need!r} (spec D1 2.3.2)")
        self.i_sulf, self.i_br, self.i_io, self.i_ag = (
            idx["sulfite"],
            idx["bromide"],
            idx["iodide"],
            idx["silver"],
        )

        # Speciation with the gelatin family appended (spec D1 2.3.4, 2.4.7).
        spec = solution.speciation
        gel_pk = _scalar(kin.gelatin.pka, "1", "gelatin.pka")
        d_ionic = -spec.log10_gamma_h  # D(I) of spec D1 2.3.3 (log10 gamma_H = -D)
        gel_family = AcidBaseFamily(
            "gelatin", (1, 0), (1, 0), (gel_pk,), "van_t_hoff", (0.0,), (), "none"
        )
        j_max = max(spec.j_max, 1)
        pka = np.full((spec.n_families + 1, j_max), np.nan)
        pka[: spec.n_families, : spec.j_max] = spec.pka
        pka[-1, 0] = gel_pk + d_ionic  # (2 z_HA - 1) D with z_HA = +1, delta eps = 0
        self.spec_bath = spec
        self.spec = Speciation(
            families=(*spec.families, gel_family),
            pka=pka,
            log10_kw=spec.log10_kw,
            log10_gamma_h=spec.log10_gamma_h,
            log10_gamma_oh=spec.log10_gamma_oh,
        )
        self.nf = spec.n_families  # solution families; gelatin is row nf
        self.fam_comp = np.array([idx[f.component] for f in spec.families], dtype=np.intp)
        self.protons = [np.asarray(f.protons_rel, dtype=np.float64) for f in self.spec.families]
        self.n_states = [f.n_states for f in self.spec.families]
        self.fam_of_comp = {int(c): f for f, c in enumerate(self.fam_comp)}
        self.mobile = [f for f in range(self.nf) if self.d_liq[self.fam_comp[f]] > 0.0]
        self.f_sulf, self.s_sulf = _sulfite_state(solution)
        self.gel_ph0 = _scalar(kin.gelatin.initial_ph, "1", "gelatin.initial_ph")

        # Developing agents and their active forms (spec D1 2.3.5, 2.4.3).
        self.n_agents = len(kin.agent_ids)
        self.agent_comp, self.agent_fam, self.agent_state = [], [], []
        for i, aid in enumerate(kin.agent_ids):
            c = idx.get(aid)
            if c is None or c not in self.fam_of_comp:
                raise ValueError(f"agent {aid!r} is not a component with an acid-base family")
            f = self.fam_of_comp[c]
            charges = spec.families[f].charges
            target = round(float(kin.charges[i]))
            if target not in charges:
                raise ValueError(f"agent {aid!r}: no state with the active charge {target}")
            self.agent_comp.append(c)
            self.agent_fam.append(f)
            self.agent_state.append(charges.index(target))
        self.eta = kin.eta_mask
        self.k_et = np.asarray(kin.k_et, dtype=np.float64)
        self.k_ads = np.asarray(kin.K_ads, dtype=np.float64)
        self.lambda0 = np.asarray(kin.lambda0, dtype=np.float64)
        self.k_des = np.asarray(kin.k_des, dtype=np.float64)
        self.k_reg = np.asarray(kin.k_reg, dtype=np.float64)
        self.abs_z = np.abs(np.asarray(kin.charges, dtype=np.float64))
        self.sulfonated = np.asarray(kin.sulfonated, dtype=bool)
        self.k_br, self.k_i = kin.K_br, kin.K_i
        self.u_gel = kin.u_gel / math.sqrt(solution.ionic_strength)  # spec D1 2.4.2
        self.u_x = kin.u_halide
        sol = kin.solvent
        self.k_diss = _scalar(sol.k_diss, "um/s", "k_diss")
        self.k_pd = _scalar(sol.k_pd, "um/s", "k_pd")
        self.c_sat_ref = _scalar(sol.c_sat_ref, "mol/L", "c_sat_ref")
        self.nu_s = _scalar(sol.sulfite_exponent, "1", "sulfite_exponent")
        self.nu_b = _scalar(sol.bromide_exponent, "1", "bromide_exponent")
        self.n_lig = _scalar(sol.ligands, "1", "ligands")
        self.c_s = _scalar(sol.sulfonation_switch, "mol/L", "sulfonation_switch")

        # Grains (spec D1 2.0.2, 2.4.4).
        classes = film.classes
        self.K = classes.count
        self.x_i = classes.iodide_fraction
        self.s_c = _ELL_REF_UM * classes.surface_um2 / classes.volume_um3
        self.c_c = np.asarray(film.silver_per_cell, dtype=np.float64)
        n_cz = classes.number_per_um2 / (self.n_z * self.delta[0])  # um^-3
        self.ns_c = n_cz * classes.surface_um2  # um^-1
        self.c_tot = float(_lsum(self.c_c))
        cap = film.latent.probabilities.shape[2] - 1
        self.n_index = np.array([0, *range(2, cap + 1)], dtype=np.intp)
        self.h_n = np.asarray(film.hazards, dtype=np.float64)[self.n_index]
        self.nn = self.n_index.size
        self.gelatin = np.asarray(film.gelatin_total, dtype=np.float64)
        self.film = film

    # -- time-dependent diffusivities ----------------------------------------------------

    def swelling(self, t: float) -> float:
        """``s(t) = 1 - exp(-t ln2 / T_half)``; 1 after a presoak (spec D1 2.4.7)."""
        if self.presoaked or self.half_time <= 0.0:
            return 1.0
        return -math.expm1(-max(t, 0.0) * math.log(2.0) / self.half_time)

    def swelling_mean(self, t0: float, t1: float) -> float:
        """Mean of ``s(t)`` over ``[t0, t1]``: the swelling factor of one column step. It is a
        second-order discretisation of the time-dependent D and tends to 1 as ``T_half -> 0``
        on every step (the presoak limit, spec D1 2.4.13)."""
        if self.presoaked or self.half_time <= 0.0:
            return 1.0
        k = math.log(2.0) / self.half_time
        return 1.0 - (math.exp(-k * max(t0, 0.0)) - math.exp(-k * max(t1, 0.0))) / (k * (t1 - t0))

    def level_diffusion(self, s: float) -> np.ndarray:
        """``(L, K_c)``: ``r_k D_k s`` in gel levels, ``D_k`` in the liquid (spec D1 2.4.7)."""
        gel_d = self.r_gel * self.d_liq * s
        return np.where(self.gel[:, None], gel_d[None, :], self.d_liq[None, :])

    def carrier_diffusion(self, s: float) -> np.ndarray:
        """``(C, L)``: D of the mobile proton carriers (families, then H+, OH-)."""
        rows = []
        for f in self.mobile:
            c = self.fam_comp[f]
            rows.append(np.where(self.gel, self.r_gel[c] * self.d_liq[c] * s, self.d_liq[c]))
        for d in (self.d_h, self.d_oh):
            rows.append(np.where(self.gel, _WATER_ION_GEL_RATIO * d * s, d))
        return np.array(rows)


def _face_mean(d: np.ndarray, delta: np.ndarray) -> np.ndarray:
    """Thickness-weighted harmonic mean of D on the faces (spec D1 2.4.8); 0 if both are 0.

    ``d`` ``(..., L)`` -> ``(..., L - 1)``.
    """
    da, db = d[..., :-1], d[..., 1:]
    la, lb = delta[:-1], delta[1:]
    den = la * db + lb * da
    with np.errstate(invalid="ignore", divide="ignore"):
        out = (la + lb) * da * db / den
    return np.where(den > 0.0, out, 0.0)


# ---------------------------------------------------------------------------------------
# Speciation and surface state with derivatives
# ---------------------------------------------------------------------------------------


@dataclass
class _Species:
    ph: np.ndarray  # (N, La)
    alpha: list  # per family (n_states, N, La)
    nbar: np.ndarray  # (F+1, N, La)
    var: np.ndarray  # (F+1, N, La)
    h: np.ndarray
    oh: np.ndarray
    beta: np.ndarray  # (N, La) total buffer capacity


def _speciate(st: _Setup, cpos: np.ndarray, toth: np.ndarray, gel: np.ndarray, ph0) -> _Species:
    """pH and moments on levels: ``cpos`` (N, La, K_c) positive totals, ``toth`` (N, La)."""
    n, la = toth.shape
    totals = np.empty((st.nf + 1, n, la))
    for f in range(st.nf):
        totals[f] = cpos[:, :, st.fam_comp[f]]
    totals[st.nf] = gel[None, :]
    ph = solve_ph(st.spec, totals, toth, ph0=ph0)
    alpha_all = fractions(st.spec, ph)  # (F+1, J+1, N, La)
    alpha, nbar, var = [], np.empty((st.nf + 1, n, la)), np.empty((st.nf + 1, n, la))
    for f in range(st.nf + 1):
        ns = st.n_states[f]
        a = alpha_all[f, :ns]
        pr = st.protons[f]
        m = pr[0] * a[0]
        for j in range(1, ns):
            m = m + pr[j] * a[j]
        v = (pr[0] - m) ** 2 * a[0]
        for j in range(1, ns):
            v = v + (pr[j] - m) ** 2 * a[j]
        alpha.append(a)
        nbar[f] = m
        var[f] = v
    h = 10.0 ** (-ph - st.spec.log10_gamma_h)
    oh = 10.0 ** (st.spec.log10_kw + ph - st.spec.log10_gamma_oh)
    beta = h + oh
    for f in range(st.nf + 1):
        beta = beta + totals[f] * var[f]
    beta = _LN10 * beta
    return _Species(ph=ph, alpha=alpha, nbar=nbar, var=var, h=h, oh=oh, beta=beta)


def _bath_carriers(st: _Setup, bulk: np.ndarray) -> np.ndarray:
    """Proton carriers ``(C,)`` of the bath (mobile ``P_f``, ``[H+]``, ``-[OH-]``)."""
    cpos = positive_part(bulk[: st.kc])
    spec = st.spec_bath
    totals = np.array([cpos[c] for c in st.fam_comp])
    ph = solve_ph(
        spec, totals[:, None], np.asarray([bulk[st.it]]), ph0=np.asarray([st.solution.ph])
    )
    alpha = fractions(spec, ph)
    out = []
    for f in st.mobile:
        pr = st.protons[f]
        m = 0.0
        for j in range(st.n_states[f]):
            m = m + pr[j] * alpha[f, j, 0]
        out.append(bulk[st.fam_comp[f]] * m)
    h = 10.0 ** (-ph[0] - spec.log10_gamma_h)
    oh = 10.0 ** (spec.log10_kw + ph[0] - spec.log10_gamma_oh)
    out.extend([h, -oh])
    return np.array(out, dtype=np.float64)


@dataclass
class _Surface:
    kappa: np.ndarray  # (N, nz)
    lam: np.ndarray
    dkappa: np.ndarray  # (A + 2, N, nz) d/d(a_0..a_{A-1}, Br, I)
    dlam: np.ndarray
    partition: np.ndarray  # (A, N, nz)
    regen: np.ndarray  # (A, N, nz) g_E (0 for regenerators)
    share: np.ndarray  # (A, N, nz) w_R = k_reg a_R / supply (0 for ETAs)


def _surface(st: _Setup, a: np.ndarray, br: np.ndarray, io: np.ndarray) -> _Surface:
    """Surface state of spec D1 2.4.3 and its derivatives, ``a`` (A, N, nz) >= 0."""
    na = st.n_agents
    shape = br.shape
    supply = np.zeros(shape)
    for i in range(na):
        if not st.eta[i]:
            supply = supply + st.k_reg[i] * a[i]
    halide = st.k_br * br + st.k_i * io
    xi = 1.0 + halide
    phi = [None] * na
    denom = [None] * na
    for i in range(na):
        if st.eta[i]:
            denom[i] = st.k_des[i] + supply
            phi[i] = st.k_et[i] / denom[i]
            xi = xi + st.k_ads[i] * a[i] * (1.0 + phi[i])
        else:
            xi = xi + st.k_ads[i] * a[i]
    # d xi / d u for u = a_0 .. a_{A-1}, Br, I
    dxi = np.zeros((na + 2, *shape))
    for j in range(na):
        if st.eta[j]:
            dxi[j] = st.k_ads[j] * (1.0 + phi[j])
        else:
            acc = np.full(shape, st.k_ads[j])
            for e in range(na):
                if st.eta[e]:
                    dphi = -st.k_et[e] * st.k_reg[j] / (denom[e] * denom[e])
                    acc = acc + st.k_ads[e] * a[e] * dphi
            dxi[j] = acc
    dxi[na] = st.k_br
    dxi[na + 1] = st.k_i
    theta = np.empty((na, *shape))
    for i in range(na):
        theta[i] = st.k_ads[i] * a[i] / xi
    thx = halide / xi
    dthx = np.empty((na + 2, *shape))
    for u in range(na + 2):
        dh = 0.0 if u < na else (st.k_br if u == na else st.k_i)
        dthx[u] = (dh - thx * dxi[u]) / xi
    barrier = st.u_gel + st.u_x * thx
    kappa = np.zeros(shape)
    lam = np.zeros(shape)
    dk = np.zeros((na + 2, *shape))
    dl = np.zeros((na + 2, *shape))
    rate = np.empty((na, *shape))
    for i in range(na):
        e_i = np.exp(-st.abs_z[i] * barrier)
        rate[i] = st.k_et[i] * theta[i]
        kappa = kappa + rate[i]
        lam = lam + st.lambda0[i] * theta[i] * e_i
        for u in range(na + 2):
            dth = -theta[i] * dxi[u] / xi
            if u == i:
                dth = dth + st.k_ads[i] / xi
            dk[u] = dk[u] + st.k_et[i] * dth
            dl[u] = dl[u] + st.lambda0[i] * e_i * (dth - theta[i] * st.abs_z[i] * st.u_x * dthx[u])
    positive = kappa > 0.0
    safe = np.where(positive, kappa, 1.0)
    partition = np.where(positive[None], rate / safe[None], 0.0)
    regen = np.zeros((na, *shape))
    share = np.zeros((na, *shape))
    has_supply = supply > 0.0
    safe_supply = np.where(has_supply, supply, 1.0)
    for i in range(na):
        if st.eta[i]:
            regen[i] = supply / denom[i]
        else:
            share[i] = np.where(has_supply, st.k_reg[i] * a[i] / safe_supply, 0.0)
    return _Surface(kappa, lam, dk, dl, partition, regen, share)


# ---------------------------------------------------------------------------------------
# Right-hand side and Jacobian of the column step
# ---------------------------------------------------------------------------------------


@dataclass
class _Kernels:
    """Kernel weights of the cells: ``omega`` (N, nz, K, Nb) couples the column to the
    kernels; ``basis`` (N, nz, K, Nb, N_n) maps the cluster terms E_n to the kernel source
    (``None``: the identity, one kernel per cluster size)."""

    omega: np.ndarray
    basis: np.ndarray | None = None


@dataclass
class _Eval:
    f_conc: np.ndarray  # (N, L, nc)
    f_loc: np.ndarray  # (N, nz, 3)
    f_g: np.ndarray  # (N, nz, K, Nb)
    f_h: np.ndarray
    f_q: np.ndarray
    phi: np.ndarray  # (N, nc) flux into the film through the Dirichlet face
    ph: np.ndarray  # (N, La)
    kappa: np.ndarray  # (N, nz)
    jac: object = None


@dataclass
class _Jac:
    column: ros.ColumnJacobian
    dphi: np.ndarray  # (N, nc) d phi / d C_top (diagonal)
    coef: np.ndarray  # (N, nz, b, 3) row coefficients of R_dev, R_diss, R_pd
    g_term: np.ndarray  # (3, N, nz, K) d R_t / d G_cb = g_t,c omega_cb
    a_g: np.ndarray  # (N, nz, K, Nb) 1/(-dfG/dG) = kappa_c, i.e. dfG/dG = -kappa_c
    kappa_c: np.ndarray  # (N, nz, K, 1)
    s_eg: np.ndarray  # (N, nz, K, Nb) s_c (E_b - G_b)
    s_gh: np.ndarray  # (N, nz, K, Nb) 2 s_c (G_b - H_b)
    de: np.ndarray  # (N, nz, K or 1, Nb) dE_b/dLambda
    grad_kappa: np.ndarray  # (N, nz, nc)
    q_g: np.ndarray  # (N, nz, K, 1) dfQ/dG
    q_ag: np.ndarray  # (N, nz, K, Nb) dfQ/dC_Ag


def _evaluate(
    st: _Setup,
    s: float,
    conc: np.ndarray,
    loc: np.ndarray,
    g: np.ndarray,
    hk: np.ndarray,
    kern: _Kernels,
    bulk: np.ndarray,
    bath_carriers: np.ndarray,
    n_act: int,
    ph0: np.ndarray,
    want_jac: bool,
) -> _Eval:
    """``f`` (and, with ``want_jac``, ``J``) of the column step with swelling factor ``s``.

    ``conc`` (N, L, nc) totals incl. TOTH; levels ``>= n_act`` are held at the bulk value
    (``f = 0`` there). ``ph0`` (N, n_act) is the previous pH.
    """
    n = conc.shape[0]
    nz, kc, nc, it = st.n_z, st.kc, st.nc, st.it
    la = n_act
    delta = st.delta[:la]
    d_lev = st.level_diffusion(s)[:la]  # (La, kc)
    d_car = st.carrier_diffusion(s)[:, :la]  # (C, La)

    c_act = conc[:, :la]
    cpos = positive_part(c_act[:, :, :kc])
    sp = _speciate(st, cpos, c_act[:, :, it], st.gelatin[:la], ph0)

    f_conc = np.zeros(conc.shape)
    # --- transport of the components (spec D1 2.4.8) -----------------------------------
    d_face = _face_mean(d_lev.T, delta).T  # (La-1, kc)
    hf = st.hf[: la - 1]
    flux = -d_face[None] * (c_act[:, 1:, :kc] - c_act[:, :-1, :kc]) / hf[None, :, None]
    d_top = d_lev[la - 1]  # liquid D at the Dirichlet face
    half = 0.5 * delta[la - 1]
    phi = np.empty((n, nc))
    phi[:, :kc] = d_top[None] * (bulk[None, :kc] - c_act[:, la - 1, :kc]) / half
    div = np.zeros((n, la, kc))
    div[:, :-1] -= flux
    div[:, 1:] += flux
    div[:, la - 1] += phi[:, :kc]
    f_conc[:, :la, :kc] = div / delta[None, :, None]

    # --- TOTH transport from the mobile carriers (spec D1 2.4.8) ------------------------
    carriers = []
    for f in st.mobile:
        carriers.append(c_act[:, :, st.fam_comp[f]] * sp.nbar[f])
    carriers.append(sp.h)
    carriers.append(-sp.oh)
    d_car_face = _face_mean(d_car, delta)  # (C, La-1)
    j_toth = np.zeros((n, la - 1))
    phi_t = np.zeros(n)
    for c, car in enumerate(carriers):
        j_toth = j_toth - d_car_face[c][None] * (car[:, 1:] - car[:, :-1]) / hf[None]
        phi_t = phi_t + d_car[c, la - 1] * (bath_carriers[c] - car[:, la - 1]) / half
    phi[:, it] = phi_t
    div_t = np.zeros((n, la))
    div_t[:, :-1] -= j_toth
    div_t[:, 1:] += j_toth
    div_t[:, la - 1] += phi_t
    f_conc[:, :la, it] = div_t / delta[None]

    # --- reactions on the emulsion levels (spec D1 2.4.2-2.4.6) --------------------------
    ce = c_act[:, :nz]
    cpe = cpos[:, :nz]
    dpos = cpe / np.sqrt(ce[:, :, :kc] ** 2 + 1e-18)  # dC+/dC = C+/sqrt(C^2 + eps^2)
    beta = sp.beta[:, :nz]
    # grad pH over (components, TOTH): dpH/dT_f = nbar_f / beta, dpH/dTOTH = -1/beta
    grad_ph = np.zeros((n, nz, nc))
    for f in range(st.nf):
        c = st.fam_comp[f]
        grad_ph[:, :, c] = sp.nbar[f, :, :nz] * dpos[:, :, c] / beta
    grad_ph[:, :, it] = -1.0 / beta

    def active(fam: int, state: int, comp: int):
        alpha = sp.alpha[fam][state][:, :nz]
        dalpha = _LN10 * alpha * (sp.nbar[fam, :, :nz] - st.protons[fam][state])
        value = alpha * cpe[:, :, comp]
        grad = (cpe[:, :, comp] * dalpha)[:, :, None] * grad_ph
        grad[:, :, comp] += alpha * dpos[:, :, comp]
        return value, grad

    na = st.n_agents
    a_val = np.empty((na, n, nz))
    a_grad = []
    for i in range(na):
        v, gr = active(st.agent_fam[i], st.agent_state[i], st.agent_comp[i])
        a_val[i] = v
        a_grad.append(gr)
    br = cpe[:, :, st.i_br]
    io = cpe[:, :, st.i_io]
    surf = _surface(st, a_val, br, io)
    grad_k = np.zeros((n, nz, nc))
    grad_l = np.zeros((n, nz, nc))
    for i in range(na):
        grad_k = grad_k + surf.dkappa[i][:, :, None] * a_grad[i]
        grad_l = grad_l + surf.dlam[i][:, :, None] * a_grad[i]
    grad_k[:, :, st.i_br] += surf.dkappa[na] * dpos[:, :, st.i_br]
    grad_k[:, :, st.i_io] += surf.dkappa[na + 1] * dpos[:, :, st.i_io]
    grad_l[:, :, st.i_br] += surf.dlam[na] * dpos[:, :, st.i_br]
    grad_l[:, :, st.i_io] += surf.dlam[na + 1] * dpos[:, :, st.i_io]
    kappa = surf.kappa

    # Initiation terms (spec D1 2.4.2): E_n = 1 - exp(-h_n Lambda), shared by all classes.
    lam_clock = loc[:, :, 0]
    expo = np.exp(-st.h_n[None, None, :] * lam_clock[:, :, None])  # (N, nz, Nn)
    e_n = -np.expm1(-st.h_n[None, None, :] * lam_clock[:, :, None])
    de_n = st.h_n[None, None, :] * expo
    if kern.basis is None:
        e_b = e_n[:, :, None, :]  # (N, nz, 1, Nb)
        de_b = de_n[:, :, None, :]
    else:
        e_b = _lsum(kern.basis * e_n[:, :, None, None, :])
        de_b = _lsum(kern.basis * de_n[:, :, None, None, :])
    omega = kern.omega
    kappa_c = st.s_c[None, None, :] * kappa[:, :, None]  # (N, nz, K)
    f_bar = _lsum(omega * e_b)  # (N, nz, K)
    x_bar = _lsum(omega * g)
    df_bar = _lsum(omega * de_b)
    r_dev = _lsum(st.c_c * kappa_c * (f_bar - x_bar))
    a_undev = _lsum(st.ns_c * (1.0 - x_bar))
    a_dev = _lsum(st.ns_c * x_bar)

    free, grad_free = active(st.f_sulf, st.s_sulf, st.i_sulf)
    br_floor = br * br + _BROMIDE_FLOOR * _BROMIDE_FLOOR
    c_sat = (
        st.c_sat_ref
        * (free / _SULFITE_REF) ** st.nu_s
        * (_BROMIDE_REF / np.sqrt(br_floor)) ** st.nu_b
    )
    grad_csat = (st.nu_s * c_sat / free)[:, :, None] * grad_free
    grad_csat[:, :, st.i_br] -= st.nu_b * c_sat * br / br_floor * dpos[:, :, st.i_br]
    ag = cpe[:, :, st.i_ag]
    width = _SOFTPLUS_WIDTH * c_sat
    x_sp = (c_sat - ag) / width
    softplus = width * np.logaddexp(0.0, x_sp)
    dsoft = expit(x_sp)
    r_diss = st.k_diss * a_undev * softplus
    r_pd = st.k_pd * a_dev * ag

    # Stoichiometric coefficients of (R_dev, R_diss, R_pd) per dense row (spec D1 2.4.6).
    coef = np.zeros((n, nz, st.b, 3))
    regen_supply = np.zeros((n, nz))
    for i in range(na):
        if st.eta[i]:
            regen_supply = regen_supply + surf.regen[i] * surf.partition[i]
    sigma = free / (free + st.c_s)
    c_sulf = np.zeros((n, nz))
    for i in range(na):
        if st.eta[i]:
            ci = 0.5 * (1.0 - surf.regen[i]) * surf.partition[i]  # E_ox per R_e
        else:
            ci = 0.5 * (surf.partition[i] + surf.share[i] * regen_supply)  # Q_R per R_e
        col = st.agent_comp[i]
        coef[:, :, col, 0] -= ci
        coef[:, :, col, 2] -= ci
        if st.sulfonated[i]:
            c_sulf = c_sulf + sigma * ci
    coef[:, :, st.i_sulf, 0] -= c_sulf
    coef[:, :, st.i_sulf, 1] -= st.n_lig
    coef[:, :, st.i_sulf, 2] += st.n_lig - c_sulf
    coef[:, :, st.i_br, 0] += 1.0 - st.x_i
    coef[:, :, st.i_br, 1] += 1.0
    coef[:, :, st.i_io, 0] += st.x_i
    coef[:, :, st.i_ag, 1] += 1.0
    coef[:, :, st.i_ag, 2] -= 1.0
    coef[:, :, it, 0] += 1.0 - c_sulf
    coef[:, :, it, 2] += 1.0 - c_sulf
    coef[:, :, st.iD, 1] += 1.0
    rates = (r_dev, r_diss, r_pd)
    react = (
        coef[..., 0] * r_dev[..., None]
        + coef[..., 1] * r_diss[..., None]
        + coef[..., 2] * r_pd[..., None]
    )  # (N, nz, b)
    f_conc[:, :nz] += react[:, :, :nc]
    f_loc = np.empty((n, nz, 3))
    f_loc[:, :, 0] = surf.lam
    f_loc[:, :, 1] = react[:, :, st.iD]
    f_loc[:, :, 2] = kappa

    # Kernels (spec D1 2.4.4).
    kc_b = kappa_c[..., None]  # (N, nz, K, 1)
    f_g = kc_b * (e_b - g)
    f_h = 2.0 * kc_b * (g - hk)
    q_g = (st.k_pd * st.ns_c / st.c_tot)[None, None, :, None] * ag[:, :, None, None]
    f_q = q_g * g

    out = _Eval(
        f_conc=f_conc, f_loc=f_loc, f_g=f_g, f_h=f_h, f_q=f_q, phi=phi, ph=sp.ph, kappa=kappa
    )
    if not want_jac:
        return out

    # --- Jacobian -------------------------------------------------------------------------
    grad_r = np.zeros((3, n, nz, st.b))
    grad_r[0, :, :, :nc] = _lsum(st.c_c * st.s_c * (f_bar - x_bar))[:, :, None] * grad_k
    grad_r[0, :, :, st.iL] = _lsum(st.c_c * kappa_c * df_bar)
    grad_r[1, :, :, :nc] = (st.k_diss * a_undev * dsoft)[:, :, None] * grad_csat
    grad_r[1, :, :, st.i_ag] -= st.k_diss * a_undev * dsoft * dpos[:, :, st.i_ag]
    grad_r[2, :, :, st.i_ag] = st.k_pd * a_dev * dpos[:, :, st.i_ag]
    dense = (
        coef[..., :, 0:1] * grad_r[0][:, :, None, :]
        + coef[..., :, 1:2] * grad_r[1][:, :, None, :]
        + coef[..., :, 2:3] * grad_r[2][:, :, None, :]
    )
    dense[:, :, st.iL, :nc] += grad_l
    dense[:, :, st.iK, :nc] += grad_k
    g_term = np.empty((3, n, nz, st.K))
    g_term[0] = -st.c_c * kappa_c
    g_term[1] = -(st.k_diss * softplus)[:, :, None] * st.ns_c
    g_term[2] = (st.k_pd * ag)[:, :, None] * st.ns_c

    # Transport entries.
    diag = np.zeros((n, st.L - nz, nc))
    lower = np.zeros((n, st.L - 1, nc))
    upper = np.zeros((n, st.L - 1, nc))
    self_c = np.zeros((la, kc))
    lo_c = d_face / (hf[:, None] * delta[1:, None])
    up_c = d_face / (hf[:, None] * delta[:-1, None])
    self_c[:-1] -= d_face / (hf[:, None] * delta[:-1, None])
    self_c[1:] -= d_face / (hf[:, None] * delta[1:, None])
    self_c[la - 1] -= d_top / (half * delta[la - 1])
    lower[:, : la - 1, :kc] = lo_c[None]
    upper[:, : la - 1, :kc] = up_c[None]
    beta_car = np.empty((len(carriers), n, la))
    k = 0
    for f in st.mobile:
        beta_car[k] = _LN10 * c_act[:, :, st.fam_comp[f]].clip(min=0.0) * sp.var[f]
        k += 1
    beta_car[k] = _LN10 * sp.h
    beta_car[k + 1] = _LN10 * sp.oh
    t_lo, t_up, t_self = ros.toth_transport_jacobian(
        beta_car, sp.beta, d_car_face[:, None, :], hf, delta
    )
    g_top = np.zeros(n)
    for c in range(len(carriers)):
        g_top = g_top + d_car[c, la - 1] * beta_car[c, :, la - 1]
    g_top = g_top / (sp.beta[:, la - 1] * half)
    t_self[:, la - 1] -= g_top / delta[la - 1]
    lower[:, : la - 1, it] = t_lo
    upper[:, : la - 1, it] = t_up
    idx = np.arange(nc)
    for lev in range(nz):
        dense[:, lev, idx[:kc], idx[:kc]] += self_c[lev][None]
        dense[:, lev, it, it] += t_self[:, lev]
    diag[:, : la - nz, :kc] = self_c[nz:][None]
    diag[:, : la - nz, it] = t_self[:, nz:]
    dphi = np.empty((n, nc))
    dphi[:, :kc] = -d_top[None] / half
    dphi[:, it] = -g_top
    column = ros.ColumnJacobian(st.layout, dense, diag, lower, upper)
    out.jac = _Jac(
        column=column,
        dphi=dphi,
        coef=coef,
        g_term=g_term,
        a_g=kc_b,
        kappa_c=kc_b,
        s_eg=st.s_c[None, None, :, None] * (e_b - g),
        s_gh=2.0 * st.s_c[None, None, :, None] * (g - hk),
        de=de_b,
        grad_kappa=grad_k,
        q_g=q_g,
        q_ag=(st.k_pd * st.ns_c / st.c_tot)[None, None, :, None]
        * g
        * dpos[:, :, st.i_ag][:, :, None, None],
    )
    out.rates = rates  # type: ignore[attr-defined]
    out.grad_rates = grad_r  # type: ignore[attr-defined]
    return out


# ---------------------------------------------------------------------------------------
# State container and the ROS2 column step with Schur elimination of the kernels
# ---------------------------------------------------------------------------------------


@dataclass
class _State:
    conc: np.ndarray  # (N, L, nc)
    loc: np.ndarray  # (N, nz, 3): Lambda, D, K
    g: np.ndarray  # (N, nz, K, Nb)
    h: np.ndarray
    q: np.ndarray
    ph: np.ndarray  # (N, L) last pH

    def copy(self) -> "_State":
        return _State(
            self.conc.copy(),
            self.loc.copy(),
            self.g.copy(),
            self.h.copy(),
            self.q.copy(),
            self.ph.copy(),
        )


def _pack(st: _Setup, conc: np.ndarray, loc: np.ndarray) -> np.ndarray:
    dense = np.concatenate([conc[:, : st.n_z], loc], axis=2)
    return st.layout.join(dense, conc[:, st.n_z :])


def _unpack(st: _Setup, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    dense, diag = st.layout.split(y)
    conc = np.concatenate([dense[:, :, : st.nc], diag], axis=1)
    return conc, np.ascontiguousarray(dense[:, :, st.nc :])


@dataclass
class _StepOut:
    state: _State
    flux: np.ndarray  # (N, nc) step-mean flux into the film
    err: tuple | None  # (err_y, y_old, y_new, err_kernels) for the adaptive controller


def _make_solver(st: _Setup, jac: _Jac, kern: _Kernels, gh: float):
    """Solver of ``(I - gh J) k = r`` for the column and its kernels (spec D1 3.3).

    The kernels' block ``J_uu`` is lower triangular (G, H, Q) per (c, z, b) and is inverted in
    closed form; the column's dense blocks get ``-(gh)^2 J_cu inv(I - gh J_uu) J_uc`` and the
    right-hand side ``gh J_cu inv(I - gh J_uu) r_u``. Returns ``solve(r_conc, r_loc, r_g,
    r_h, r_q) -> (k_conc, k_loc, k_g, k_h, k_q)``.
    """
    n, nz = jac.coef.shape[0], st.n_z
    kc_b = jac.kappa_c
    a_g = 1.0 + gh * kc_b  # (N, nz, K, 1)
    a_h = 1.0 + 2.0 * gh * kc_b
    omega = kern.omega
    # Schur update of the dense blocks: -(gh)^2 sum_t coef_t (A_t grad_kappa + B_t e_Lambda).
    w_t = [jac.g_term[t_][..., None] * omega / a_g for t_ in range(3)]  # (N, nz, K, Nb)
    flat = (n, nz, -1)
    a_t = [_lsum((w * jac.s_eg).reshape(flat)) for w in w_t]
    b_t = [_lsum((w * kc_b * jac.de).reshape(flat)) for w in w_t]
    upd = np.zeros((n, nz, st.b, st.b))
    for t_ in range(3):
        row = np.zeros((n, nz, st.b))
        row[:, :, : st.nc] = a_t[t_][:, :, None] * jac.grad_kappa
        row[:, :, st.iL] += b_t[t_]
        upd = upd + jac.coef[..., :, t_ : t_ + 1] * row[:, :, None, :]
    factor = ros.factor_columns(jac.column, gh, dense_update=-(gh * gh) * upd)

    def solve(r_conc, r_loc, r_g, r_h, r_q):
        rhs = _pack(st, r_conc, r_loc)
        dense_rhs, diag_rhs = st.layout.split(rhs)
        extra = np.zeros((n, nz, st.b))
        for t_ in range(3):
            s_t = _lsum((w_t[t_] * r_g).reshape(flat))
            extra = extra + jac.coef[..., t_] * s_t[:, :, None]
        dense_rhs = dense_rhs + gh * extra
        k_col = ros.solve_columns(factor, st.layout.join(dense_rhs, diag_rhs))
        k_conc, k_loc = _unpack(st, k_col)
        dot_k = _lsum(jac.grad_kappa * k_conc[:, :nz])  # (N, nz)
        k_g = (
            r_g
            + gh * (jac.s_eg * dot_k[:, :, None, None] + kc_b * jac.de * k_loc[:, :, 0, None, None])
        ) / a_g
        k_h = (r_h + gh * (jac.s_gh * dot_k[:, :, None, None] + 2.0 * kc_b * k_g)) / a_h
        k_q = r_q + gh * (jac.q_ag * k_conc[:, :nz, st.i_ag][:, :, None, None] + jac.q_g * k_g)
        return k_conc, k_loc, k_g, k_h, k_q

    return solve


def _c_step(
    st: _Setup,
    state: _State,
    kern: _Kernels,
    t: float,
    h: float,
    bulk0: np.ndarray,
    bulk1: np.ndarray,
    n_act: int,
    want_err: bool,
) -> _StepOut:
    """One ROS2 column step (spec D1 3.3) with the kernels eliminated by a Schur complement."""
    gh = ros.GAMMA * h
    car0 = _bath_carriers(st, bulk0)
    swell = st.swelling_mean(t, t + h)
    ev0 = _evaluate(
        st,
        swell,
        state.conc,
        state.loc,
        state.g,
        state.h,
        kern,
        bulk0,
        car0,
        n_act,
        state.ph[:, :n_act],
        True,
    )
    jac: _Jac = ev0.jac
    solve = _make_solver(st, jac, kern, gh)
    if callable(bulk1):  # tank coupling: the bulk at t + h predicted from the face flux at t
        bulk1 = bulk1(ev0.phi, h)
    car1 = car0 if bulk1 is bulk0 else _bath_carriers(st, bulk1)

    top = n_act - 1
    k1 = solve(ev0.f_conc, ev0.f_loc, ev0.f_g, ev0.f_h, ev0.f_q)
    kphi1 = ev0.phi + gh * jac.dphi * k1[0][:, top]
    s1 = _State(
        state.conc + h * k1[0],
        state.loc + h * k1[1],
        state.g + h * k1[2],
        state.h + h * k1[3],
        state.q + h * k1[4],
        state.ph,
    )
    ev1 = _evaluate(st, swell, s1.conc, s1.loc, s1.g, s1.h, kern, bulk1, car1, n_act, ev0.ph, False)
    k2 = solve(
        ev1.f_conc - 2.0 * k1[0],
        ev1.f_loc - 2.0 * k1[1],
        ev1.f_g - 2.0 * k1[2],
        ev1.f_h - 2.0 * k1[3],
        ev1.f_q - 2.0 * k1[4],
    )
    kphi2 = ev1.phi - 2.0 * kphi1 + gh * jac.dphi * k2[0][:, top]
    parts = [state.conc, state.loc, state.g, state.h, state.q]
    new = [p + (1.5 * h) * a + (0.5 * h) * b for p, a, b in zip(parts, k1, k2, strict=True)]
    ph = state.ph.copy()
    ph[:, :n_act] = ev1.ph
    new_state = _State(new[0], new[1], new[2], new[3], new[4], ph)
    flux = 1.5 * kphi1 + 0.5 * kphi2
    err = None
    if want_err:
        err = tuple((0.5 * h) * (a + b) for a, b in zip(k1, k2, strict=True))
    return _StepOut(new_state, flux, err)


# ---------------------------------------------------------------------------------------
# Members: latent distributions, p-bar aggregation and output synthesis
# ---------------------------------------------------------------------------------------


def _check_exposure(members: Members, n_z: int) -> tuple[np.ndarray, np.ndarray]:
    if not isinstance(members, Members):
        raise TypeError("members must be a Members")
    ev = np.asarray(members.exposure_volume, dtype=np.float64)
    es = np.asarray(members.exposure_surface, dtype=np.float64)
    if ev.shape != es.shape or ev.ndim < 2 or ev.shape[-1] != n_z:
        raise ValueError(
            f"exposure_volume and exposure_surface must share a shape (..., n_z = {n_z}), got "
            f"{ev.shape} and {es.shape}"
        )
    if not (np.all(np.isfinite(ev)) and np.all(np.isfinite(es))):
        raise ValueError("exposures must be finite")
    if np.any(ev < 0.0) or np.any(es < 0.0):
        raise ValueError("exposures must be >= 0")
    return ev, es


def _latent(st: _Setup, ev: np.ndarray, es: np.ndarray) -> np.ndarray:
    """``p_c,z(n)`` on the kernel cluster sizes, ``(..., nz, K, N_n)`` in *level* order
    (bottom first); ``ev``, ``es`` ``(..., n_z)`` in depth order (z = 1 first)."""
    classes = st.film.classes
    table = st.film.latent
    ev = ev[..., ::-1]
    es = es[..., ::-1]
    out = np.empty((*ev.shape, st.K, st.nn))
    for c in range(st.K):
        m = classes.volume_um3[c] * ev + classes.surface_um2[c] * es  # spec D1 2.2.1
        p = table.lookup(m, c)
        out[..., c, :] = p[..., st.n_index]
    return out


def _zero_latent(st: _Setup) -> np.ndarray:
    """Latent distribution of unexposed film ``(K, N_n)``."""
    zero = np.zeros((1, st.n_z))
    return _latent(st, zero, zero)[0, 0]


def _synthesise(st: _Setup, p, g, hk, q, lam, dz, kappa):
    """Depth-merged member outputs (spec D1 2.4.10) from member ``p`` and cell kernels.

    ``p``, ``g``, ``hk``, ``q`` ``(M, nz, K, N_n)``; ``lam``, ``dz``, ``kappa`` ``(M, nz)``.
    """
    classes = st.film.classes
    nz = st.n_z
    x = _lsum(p * g)  # (M, nz, K)
    y = _lsum(p * hk)
    pq = _lsum(_lsum(p * q))  # (M, nz)
    e_n = -np.expm1(-st.h_n[None, None, :] * lam[:, :, None])
    f = _lsum(p * e_n[:, :, None, :])
    xz = _lsum(np.moveaxis(x, 1, 2))  # (M, K) sum over z
    silver = 1e12 * classes.silver_g_per_grain * classes.number_per_um2 / nz * xz
    initiated = _lsum(np.moveaxis(f, 1, 2)) / nz
    silver_m2 = _lsum(np.moveaxis(y, 1, 2)) / nz
    m_coat = classes.silver_coverage_g_m2
    physical = m_coat / nz * _lsum(pq)
    dissolved = _HALIDE_G_M2 * _lsum(dz * st.delta[:nz])
    halide = m_coat - _lsum(silver) - dissolved
    growth = _lsum(kappa) / nz
    return silver, initiated, silver_m2, physical, halide, growth


@dataclass
class _Raster:
    """Pixel window over the cell grid (layout ``"raster"``)."""

    ny: int
    nx: int
    spacing: float
    pixel: float
    origin: tuple[float, float]
    h: int
    w: int
    inner_cells: tuple[slice, slice]
    rows: tuple[int, int]  # inner pixel rows [r0, r1)
    cols: tuple[int, int]


_ROW_BATCH = 64


def _overlap(n_pix: int, pixel: float, origin: float, n_cells: int, spacing: float) -> np.ndarray:
    """1-D overlap lengths ``(n_pix, n_cells)`` of pixels with cells, um."""
    lo = origin + pixel * np.arange(n_pix)
    hi = lo + pixel
    c_lo = spacing * np.arange(n_cells)
    c_hi = c_lo + spacing
    return np.clip(
        np.minimum(hi[:, None], c_hi[None]) - np.maximum(lo[:, None], c_lo[None]), 0.0, None
    )


def _raster_pbar(st: _Setup, ras: _Raster, ev: np.ndarray, es: np.ndarray) -> np.ndarray:
    """``p-bar`` of every cell ``(ny * nx, nz, K, N_n)`` by area overlap (spec D1 2.4.4, 3.6):
    cell area not covered by pixels counts as unexposed film."""
    oy = _overlap(ras.h, ras.pixel, ras.origin[0], ras.ny, ras.spacing)
    ox = _overlap(ras.w, ras.pixel, ras.origin[1], ras.nx, ras.spacing)
    acc = np.zeros((ras.ny, ras.nx, st.n_z, st.K, st.nn))
    for r0 in range(0, ras.h, _ROW_BATCH):
        r1 = min(r0 + _ROW_BATCH, ras.h)
        p = _latent(st, ev[r0:r1], es[r0:r1])  # (rows, w, nz, K, Nn)
        part = np.einsum("rwzkn,wb->rbzkn", p, ox, optimize=False)
        acc += np.einsum("rbzkn,ra->abzkn", part, oy[r0:r1], optimize=False)
    area = ras.spacing * ras.spacing
    covered = np.einsum("ra,wb->ab", oy, ox, optimize=False)
    rest = np.clip(area - covered, 0.0, None)
    acc += rest[:, :, None, None, None] * _zero_latent(st)[None, None, None]
    return (acc / area).reshape(ras.ny * ras.nx, st.n_z, st.K, st.nn)


def _bilinear(n_cells: int, spacing: float, centres: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Lower cell index and weight of the upper cell for bilinear interpolation between cell
    centres (clamped at the grid edges)."""
    u = centres / spacing - 0.5
    if n_cells == 1:
        return np.zeros(centres.size, dtype=np.intp), np.zeros(centres.size)
    i0 = np.clip(np.floor(u).astype(np.intp), 0, n_cells - 2)
    w = np.clip(u - i0, 0.0, 1.0)
    return i0, w


def _raster_outputs(st: _Setup, ras: _Raster, ev, es, state: _State, kappa: np.ndarray):
    """Member outputs of the inner pixels (row-major) by bilinear interpolation of the cell
    kernels at the pixel centres (spec D1 3.8)."""
    r0, r1 = ras.rows
    c0, c1 = ras.cols
    ys = ras.origin[0] + ras.pixel * (np.arange(r0, r1) + 0.5)
    xs = ras.origin[1] + ras.pixel * (np.arange(c0, c1) + 0.5)
    iy, wy = _bilinear(ras.ny, ras.spacing, ys)
    ix, wx = _bilinear(ras.nx, ras.spacing, xs)
    nyx = (ras.ny, ras.nx)
    fields = {
        "g": state.g.reshape(*nyx, *state.g.shape[1:]),
        "h": state.h.reshape(*nyx, *state.h.shape[1:]),
        "q": state.q.reshape(*nyx, *state.q.shape[1:]),
        "lam": state.loc[:, :, 0].reshape(*nyx, st.n_z),
        "dz": state.loc[:, :, 1].reshape(*nyx, st.n_z),
        "kappa": kappa.reshape(*nyx, st.n_z),
    }
    outs = []
    for b0 in range(0, r1 - r0, _ROW_BATCH):
        b1 = min(b0 + _ROW_BATCH, r1 - r0)
        rows = slice(b0, b1)
        interp = {}
        for name, a in fields.items():
            a_lo = a[iy[rows]]
            a_hi = a[iy[rows] + min(1, ras.ny - 1)]
            wy_b = wy[rows].reshape((-1,) + (1,) * (a.ndim - 1))
            col = a_lo + wy_b * (a_hi - a_lo)  # (rows, nx, ...)
            c_lo = col[:, ix]
            c_hi = col[:, ix + min(1, ras.nx - 1)]
            wx_b = wx.reshape((1, -1) + (1,) * (a.ndim - 2))
            interp[name] = (c_lo + wx_b * (c_hi - c_lo)).reshape((-1, *a.shape[2:]))
        p = _latent(st, ev[r0 + b0 : r0 + b1, c0:c1], es[r0 + b0 : r0 + b1, c0:c1])
        p = p.reshape(-1, *p.shape[2:])
        outs.append(
            _synthesise(
                st,
                p,
                interp["g"],
                interp["h"],
                interp["q"],
                interp["lam"],
                interp["dz"],
                interp["kappa"],
            )
        )
    return tuple(np.concatenate([o[i] for o in outs]) for i in range(6))


# ---------------------------------------------------------------------------------------
# Lateral half step (spec D1 3.2, 3.5)
# ---------------------------------------------------------------------------------------


def _lateral_half(
    st: _Setup,
    state: _State,
    lateral: Lateral,
    inner: tuple[slice, ...],
    tau: float,
    t_mid: float,
    n_act: int,
) -> None:
    la = n_act
    conc = state.conc[:, :la]
    cpos = positive_part(conc[:, :, : st.kc])
    sp = _speciate(st, cpos, conc[:, :, st.it], st.gelatin[:la], state.ph[:, :la])
    state.ph[:, :la] = sp.ph
    s = st.swelling(t_mid)
    d = np.zeros((la, st.nc))
    d[:, : st.kc] = st.level_diffusion(s)[:la]
    d_car = st.carrier_diffusion(s)[:, :la]  # (C, La)
    beta = []
    for f in st.mobile:
        beta.append(_LN10 * cpos[:, :, st.fam_comp[f]] * sp.var[f])
    beta.append(_LN10 * sp.h)
    beta.append(_LN10 * sp.oh)
    beta = np.array(beta)  # (C, N, La)
    gel_bound = st.gelatin[None, :la] * sp.nbar[st.nf]  # (N, La)
    n = conc.shape[0]
    if lateral.boundary == "pair":
        cells_shape: tuple[int, ...] = (n,)
    else:
        cells_shape = tuple(lateral.shape_cells)
    beta_cells = beta.reshape(beta.shape[0], *cells_shape, la)
    d_toth = lat.toth_lateral_diffusion(beta_cells, d_car, inner=inner)
    toth = lat.TothLateral(
        index=st.it, gel_bound=gel_bound.reshape(*cells_shape, la), diffusion_um2_s=d_toth
    )
    c = conc.reshape(*cells_shape, la, st.nc)
    out = lat.lateral_diffuse(
        c,
        d,
        tau,
        boundary=lateral.boundary,
        spacing_um=None if lateral.boundary == "pair" else lateral.spacing_um,
        wavenumber_per_um=lateral.wavenumber_per_um,
        toth=toth,
    )
    state.conc[:, :la] = out.reshape(n, la, st.nc)


# ---------------------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------------------


@dataclass
class _Run:
    """Everything a run produced (the public result plus the final cell state, for tests)."""

    result: EnsembleResult
    state: _State
    setup: _Setup
    bulk: np.ndarray  # (S+1, nc)
    held_final: bool


def _initial_state(st: _Setup, n: int, bulk0: np.ndarray, nb: int) -> _State:
    conc = np.zeros((n, st.L, st.nc))
    conc[:, ~st.gel] = bulk0[None, None, :]
    # Gel levels hold only gelatin at pH_gel,0 (spec D1 2.4.8): TOTH = T_gel nbar_gel + H - OH.
    ph0 = st.gel_ph0
    h = 10.0 ** (-ph0 - st.spec.log10_gamma_h)
    oh = 10.0 ** (st.spec.log10_kw + ph0 - st.spec.log10_gamma_oh)
    gel_pk = st.spec.pka[st.nf, 0]
    nbar_gel = 1.0 / (1.0 + 10.0 ** (ph0 - gel_pk))
    gel_toth = st.gelatin * nbar_gel + h - oh
    conc[:, st.gel, st.it] = gel_toth[st.gel][None, :]
    ph = np.where(st.gel, ph0, st.solution.ph)
    ph = np.broadcast_to(ph, (n, st.L)).copy()
    zeros = np.zeros((n, st.n_z, st.K, nb))
    return _State(conc, np.zeros((n, st.n_z, 3)), zeros, zeros.copy(), zeros.copy(), ph)


def _inner_slices(lateral: Lateral | None, weights: np.ndarray) -> tuple[slice, ...]:
    if lateral is None or lateral.boundary == "pair":
        return (slice(None),)
    w = weights.reshape(lateral.shape_cells) > 0.0
    rows = np.nonzero(w.any(axis=1))[0]
    cols = np.nonzero(w.any(axis=0))[0]
    return (slice(int(rows[0]), int(rows[-1]) + 1), slice(int(cols[0]), int(cols[-1]) + 1))


def _history_steps(times: np.ndarray, count: int) -> np.ndarray:
    targets = np.linspace(times[0], times[-1], count)
    idx = np.abs(times[None, :] - targets[:, None]).argmin(axis=1)
    return np.unique(idx)


def _error_norm(st: _Setup, err, old: _State, new: _State, rtol: float, scale_c: np.ndarray):
    loc_scale = np.array([1.0, max(st.c_sat_ref, 1e-6), 1.0])
    parts = [
        (err[0], old.conc, new.conc, scale_c),
        (err[1], old.loc, new.loc, loc_scale),
        (err[2], old.g, new.g, np.ones(1)),
        (err[3], old.h, new.h, np.ones(1)),
        (err[4], old.q, new.q, np.full(1, 1e-3)),
    ]
    return max(ros.error_norm(e, a, b, rtol=rtol, scale=s) for e, a, b, s in parts)


def _integrate(
    film: FilmColumnModel,
    solution: WorkingSolution,
    kin: Kinetics,
    schedule: Schedule,
    members: Members,
    *,
    n_cells: int,
    cell_weights: np.ndarray,
    lateral: Lateral | None = None,
    tank: TankSpec | None = None,
    bulk: BulkTrajectory | None = None,
    presoaked: bool = False,
    method: str = "fixed",
    rtol: float = 1e-4,
    record_cells: int = 0,
    initial: tuple[np.ndarray, np.ndarray] | None = None,
    kernels: _Kernels | None = None,
) -> _Run:
    """The integrator behind :func:`run_ensemble`; ``initial`` ``(conc (N, L, K_c), toth (N,
    L))`` replaces the dry / presoaked initial state and ``kernels`` the member coupling (both
    for tests)."""
    if (tank is None) == (bulk is None):
        raise ValueError("give exactly one of tank (exact coupling) and bulk (given trajectory)")
    if method not in ("fixed", "adaptive"):
        raise ValueError(f"method must be 'fixed' or 'adaptive', got {method!r}")
    rtol = float(rtol)
    if not (math.isfinite(rtol) and rtol > 0.0):
        raise ValueError(f"rtol must be finite and > 0, got {rtol!r}")
    if (
        isinstance(record_cells, bool)
        or not isinstance(record_cells, (int, np.integer))
        or (record_cells < 0)
    ):
        raise ValueError(f"record_cells must be an int >= 0, got {record_cells!r}")
    if isinstance(n_cells, bool) or not isinstance(n_cells, (int, np.integer)) or n_cells < 1:
        raise ValueError(f"n_cells must be an int >= 1, got {n_cells!r}")
    n = int(n_cells)
    weights = np.asarray(cell_weights, dtype=np.float64)
    if weights.shape != (n,) or not np.all(np.isfinite(weights)) or np.any(weights < 0.0):
        raise ValueError(f"cell_weights must be finite and >= 0 with shape ({n},)")
    if not np.sum(weights) > 0.0:
        raise ValueError("cell_weights must have a positive sum")
    if not isinstance(schedule, Schedule):
        raise TypeError("schedule must be a Schedule")
    times = np.asarray(schedule.times_s, dtype=np.float64)
    agitated = np.asarray(schedule.agitated, dtype=bool)
    s_steps = agitated.size

    st = _Setup(film, solution, kin, presoaked)
    ev, es = _check_exposure(members, st.n_z)
    raster = None
    if members.layout == "columns":
        if ev.shape != (n, st.n_z):
            raise ValueError(
                f"columns layout needs exposures of shape ({n}, {st.n_z}), got {ev.shape}"
            )
        p_member = _latent(st, ev, es)
        p_bar = p_member
    elif members.layout == "raster":
        if lateral is None or lateral.boundary == "pair":
            raise ValueError("the raster layout needs a grid Lateral (neumann or periodic)")
        if ev.ndim != 3 or members.pixel_um is None or not members.pixel_um > 0.0:
            raise ValueError("raster members need exposures (h, w, n_z) and pixel_um > 0")
        ny, nx = (int(v) for v in lateral.shape_cells)
        inner = _inner_slices(lateral, weights)
        sp_um = float(lateral.spacing_um)
        pix = float(members.pixel_um)
        oy, ox = (float(v) for v in members.origin_um)

        def owned(sl: slice, origin: float, size: int) -> tuple[int, int]:
            a = math.ceil((sl.start * sp_um - origin) / pix - 0.5)
            b = math.ceil((sl.stop * sp_um - origin) / pix - 0.5)
            return max(a, 0), min(b, size)

        raster = _Raster(
            ny,
            nx,
            sp_um,
            pix,
            (oy, ox),
            ev.shape[0],
            ev.shape[1],
            inner,
            owned(inner[0], oy, ev.shape[0]),
            owned(inner[1], ox, ev.shape[1]),
        )
        p_member = None
        p_bar = _raster_pbar(st, raster, ev, es)
    else:
        raise ValueError(f"members.layout must be 'columns' or 'raster', got {members.layout!r}")

    if lateral is not None:
        if lateral.boundary not in lat.BOUNDARIES:
            raise ValueError(f"lateral.boundary must be one of {lat.BOUNDARIES}")
        if lateral.boundary == "pair":
            k = np.asarray(lateral.wavenumber_per_um, dtype=np.float64)
            if k.ndim != 1 or 2 * k.size != n:
                raise ValueError(f"pair lateral needs n_cells = 2 * n_pairs, got {n} and {k.shape}")
        elif int(np.prod(lateral.shape_cells)) != n:
            raise ValueError(f"lateral.shape_cells {lateral.shape_cells} does not hold {n} cells")
    inner_sl = _inner_slices(lateral, weights) if lateral is not None else (slice(None),)

    if kernels is None:
        kernels = _Kernels(omega=np.ascontiguousarray(p_bar))
    nb = kernels.omega.shape[-1]

    if tank is not None:
        coef = tk.tank_coefficient(film_area_um2=tank.film_area_um2, volume_l=tank.volume_l)
        bulk0 = np.concatenate([np.asarray(solution.totals, dtype=np.float64), [solution.toth]])
        traj = None
    else:
        coef = 0.0
        if not isinstance(bulk, BulkTrajectory):
            raise TypeError("bulk must be a BulkTrajectory")
        bt = np.asarray(bulk.times_s, dtype=np.float64)
        if bt.shape != times.shape or not np.array_equal(bt, times):
            raise ValueError("bulk.times_s must equal schedule.times_s")
        traj = np.concatenate(
            [
                np.asarray(bulk.totals, dtype=np.float64),
                np.asarray(bulk.toth, dtype=np.float64)[:, None],
            ],
            axis=1,
        )
        if traj.shape != (times.size, st.nc) or not np.all(np.isfinite(traj)):
            raise ValueError(
                f"bulk trajectory must be finite with shape ({times.size}, {st.kc}) + toth"
            )
        bulk0 = traj[0].copy()

    state = _initial_state(st, n, bulk0, nb)
    if initial is not None:
        c0, t0 = (np.asarray(a, dtype=np.float64) for a in initial)
        state.conc[:, :, : st.kc] = c0
        state.conc[:, :, st.it] = t0
    scale_c = ros.concentration_scale(bulk0[: st.kc])
    scale_all = np.concatenate([scale_c, [max(abs(float(bulk0[st.it])), ros.MIN_SCALE_MOL_L)]])

    flux_mean = np.zeros((s_steps, st.nc))
    bulk_hist = np.empty((times.size, st.nc))
    bulk_hist[0] = bulk0
    bulk_n = bulk0.copy()
    min_scaled = ros.min_scaled_concentration(state.conc[:, :, : st.kc], scale_c)
    record_at = {int(i): float(t) for t, i in schedule.record_steps.items()}
    snapshots: dict[float, EnsembleResult] = {}
    hist_steps = _history_steps(times, int(record_cells)) if record_cells else np.zeros(0, int)
    hist_lam, hist_k = [], []
    held = False
    h_try = _ADAPTIVE_H0_S

    def snapshot(step: int) -> EnsembleResult:
        return _result(
            st,
            state,
            schedule,
            step,
            raster,
            ev,
            es,
            p_member,
            flux_mean,
            bulk_hist,
            min_scaled,
            held,
            None,
            bulk_n,
        )

    if 0 in record_at:
        snapshots[record_at[0]] = snapshot(0)
    if hist_steps.size and hist_steps[0] == 0:
        hist_lam.append(state.loc[:, :, 0].copy())
        hist_k.append(state.loc[:, :, 2].copy())

    for i in range(s_steps):
        t0, t1 = float(times[i]), float(times[i + 1])
        h = t1 - t0
        agit = bool(agitated[i])
        start = agit and (i == 0 or not agitated[i - 1])
        n_act = st.agit_cell + 1 if agit else st.L
        d_inv = np.zeros(st.nc)
        if start and st.n_up:
            inv = tk.mixing_inventory(state.conc[:, st.agit_cell + 1 :], st.up_thickness, bulk_n)
            d_inv = tk.area_mean(inv, weights)
        if traj is None:
            c_plus = bulk_n + coef * d_inv / (1.0 + coef * st.sum_up) if start else bulk_n
            stage0 = c_plus
            damp = 1.0 + coef * st.sum_up if agit else 1.0

            def stage1(phi, hs, c=c_plus, damp=damp):
                # Predicted bulk at the end of the step for the second ROS2 stage; the tank
                # itself is advanced once per step with the booked flux (spec D1 3.2).
                if coef == 0.0:
                    return c
                return c - coef * hs * tk.area_mean(phi, weights) / damp
        else:
            stage0, stage1 = traj[i], traj[i + 1]
        if agit:
            state.conc[:, st.agit_cell + 1 :] = stage0[None, None, :]
        if method == "fixed":
            if lateral is not None:
                _lateral_half(st, state, lateral, inner_sl, 0.5 * h, t0 + 0.25 * h, n_act)
            out = _c_step(st, state, kernels, t0, h, stage0, stage1, n_act, False)
            state = out.state
            if lateral is not None:
                _lateral_half(st, state, lateral, inner_sl, 0.5 * h, t0 + 0.75 * h, n_act)
            col_flux = out.flux
        else:
            if start or i == 0:
                h_try = min(h_try, _ADAPTIVE_H0_S)
            t = t0
            acc = np.zeros((n, st.nc))
            while t < t1:
                hs = min(h_try, t1 - t)
                if t1 - (t + hs) < 1e-9 * h:
                    hs = t1 - t
                if traj is None:
                    b0, b1 = stage0, stage1
                else:
                    b0 = tk.bulk_at(times, traj, t)
                    b1 = tk.bulk_at(times, traj, min(t + hs, t1))
                trial = state.copy()
                if lateral is not None:
                    _lateral_half(st, trial, lateral, inner_sl, 0.5 * hs, t + 0.25 * hs, n_act)
                out = _c_step(st, trial, kernels, t, hs, b0, b1, n_act, True)
                new = out.state
                if lateral is not None:
                    _lateral_half(st, new, lateral, inner_sl, 0.5 * hs, t + 0.75 * hs, n_act)
                e_norm = _error_norm(st, out.err, trial, new, rtol, scale_all)
                m_s = ros.min_scaled_concentration(new.conc[:, :n_act, : st.kc], scale_c)
                ok = ros.adaptive_accept(e_norm, m_s)
                if ok:
                    state = new
                    acc = acc + hs * out.flux
                    t = t1 if hs == t1 - t else t + hs
                h_try = ros.next_step_size(hs, e_norm, ok)
            col_flux = acc / h
        j_mean = tk.area_mean(col_flux, weights)
        if traj is None:
            if agit:
                bulk_next = c_plus - coef * h * j_mean / (1.0 + coef * st.sum_up)
            else:
                bulk_next = bulk_n - coef * h * j_mean
        else:
            bulk_next = traj[i + 1]
        booked = j_mean - d_inv / h
        if agit:
            booked = booked + st.sum_up * (bulk_next - bulk_n) / h
        if traj is None:
            bulk_next = tk.advance_bulk(bulk_n, booked, h, coef)
        flux_mean[i] = booked
        bulk_n = np.array(bulk_next, dtype=np.float64)
        bulk_hist[i + 1] = bulk_n
        held = agit
        if agit:
            state.conc[:, st.agit_cell + 1 :] = bulk_n[None, None, :]
        min_scaled = min(
            min_scaled, ros.min_scaled_concentration(state.conc[:, :, : st.kc], scale_c)
        )
        if np.any(hist_steps == i + 1):
            hist_lam.append(state.loc[:, :, 0].copy())
            hist_k.append(state.loc[:, :, 2].copy())
        if (i + 1) in record_at:
            snapshots[record_at[i + 1]] = snapshot(i + 1)

    history = None
    if record_cells:
        spacing = _HISTORY_DEFAULT_SPACING_UM if lateral is None else float(lateral.spacing_um)
        history = CellHistories(
            times_s=readonly_array(times[hist_steps]),
            lambda_clock=readonly_array(np.array(hist_lam)[:, :, ::-1]),
            kappa_integral=readonly_array(np.array(hist_k)[:, :, ::-1]),
            cell_spacing_um=spacing,
        )
    result = _result(
        st,
        state,
        schedule,
        s_steps,
        raster,
        ev,
        es,
        p_member,
        flux_mean,
        bulk_hist,
        min_scaled,
        held,
        history,
        bulk_n,
        snapshots,
    )
    return _Run(result=result, state=state, setup=st, bulk=bulk_hist, held_final=held)


def _result(
    st,
    state,
    schedule,
    step,
    raster,
    ev,
    es,
    p_member,
    flux_mean,
    bulk_hist,
    min_scaled,
    held,
    history,
    bulk_n,
    snapshots=None,
) -> EnsembleResult:
    t = float(schedule.times_s[step])
    n_act = st.agit_cell + 1 if held else st.L
    car = _bath_carriers(st, bulk_n)
    ev_end = _evaluate(
        st,
        st.swelling(t),
        state.conc,
        state.loc,
        state.g,
        state.h,
        _Kernels(omega=np.zeros_like(state.g)),
        bulk_n,
        car,
        n_act,
        state.ph[:, :n_act],
        False,
    )
    kappa = ev_end.kappa
    if raster is None:
        outs = _synthesise(
            st, p_member, state.g, state.h, state.q, state.loc[:, :, 0], state.loc[:, :, 1], kappa
        )
    else:
        outs = _raster_outputs(st, raster, ev, es, state, kappa)
    traj = bulk_hist[: step + 1]
    bulk = BulkTrajectory(
        times_s=readonly_array(schedule.times_s[: step + 1]),
        totals=readonly_array(traj[:, : st.kc]),
        toth=readonly_array(traj[:, st.it]),
    )
    return EnsembleResult(
        silver=readonly_array(outs[0]),
        initiated=readonly_array(outs[1]),
        silver_m2=readonly_array(outs[2]),
        silver_physical=readonly_array(outs[3]),
        halide=readonly_array(outs[4]),
        growth_rate=readonly_array(outs[5]),
        snapshots=FrozenDict(snapshots or {}),
        flux_mean=readonly_array(flux_mean[:step]),
        bulk=bulk,
        min_scaled_concentration=float(min_scaled),
        cell_history=history,
    )


def run_ensemble(
    film: FilmColumnModel,
    solution: WorkingSolution,
    kin: Kinetics,
    schedule: Schedule,
    members: Members,
    *,
    n_cells: int,
    cell_weights: np.ndarray,
    lateral: Lateral | None = None,
    tank: TankSpec | None = None,
    bulk: BulkTrajectory | None = None,
    presoaked: bool = False,
    method: Literal["fixed", "adaptive"] = "fixed",
    rtol: float = 1e-4,
    record_cells: int = 0,
) -> EnsembleResult:
    """Develop ``n_cells`` cells on ``schedule`` (spec D1 2.4, 3.2-3.5, 3.7).

    Give exactly one of ``tank`` (exact coupling: strips, the prepass) and ``bulk`` (a given
    trajectory: image blocks). ``cell_weights`` ``(n_cells,)`` weigh the cells in the
    area-mean flux (0 for halo cells). ``lateral`` adds the lateral diffusion split; without it
    the cells are independent. ``method="fixed"`` steps on ``schedule.times_s``;
    ``"adaptive"`` (reference profile) subdivides every schedule step under ROS2 error
    control with ``rtol``. ``record_cells > 0`` records :class:`CellHistories` on the steps
    nearest ``record_cells`` equally spaced times. Outputs are per member: ``(n_cells, ...)``
    for ``"columns"``, the inner pixels row-major for ``"raster"``.
    """
    return _integrate(
        film,
        solution,
        kin,
        schedule,
        members,
        n_cells=n_cells,
        cell_weights=cell_weights,
        lateral=lateral,
        tank=tank,
        bulk=bulk,
        presoaked=presoaked,
        method=method,
        rtol=rtol,
        record_cells=record_cells,
    ).result
