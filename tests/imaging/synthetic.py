"""Synthetic developed ``bw-mechanistic`` states for the D1d tests (fields of spec 2.0.3).

The development step does not exist yet, so the tests build states with the fields it
writes directly. Loaded by path from tests/grain and tests/develop as well.
"""

import copy
from pathlib import Path

import numpy as np

from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import Material
from darkroom_atlas.emulsion.classes import size_classes
from darkroom_atlas.spectral import DEFAULT_GRID
from darkroom_atlas.state import EmulsionState, FrameGeometry

REPO = Path(__file__).resolve().parents[2]
FIXTURES = REPO / "tests" / "fixtures" / "library"
DATA = REPO / "data"
FILM = "test/bw-mechanistic"


def _strip_bounds(node: dict) -> None:
    node.pop("bounds", None)
    node.pop("prior", None)


def material(grain: dict | None = None, development: dict | None = None) -> Material:
    """The fixture film with grain / development values replaced (values in the file's units).

    ``grain`` keys: ``morphology``, ``mean_diameter``, ``sigma_log``, ``silver_coverage``,
    ``thickness``; ``development`` keys: any :class:`EmulsionDevelopment` field. Replaced
    values lose their bounds and prior.
    """
    data = copy.deepcopy(Library(FIXTURES).material(FILM).model_dump(mode="json"))
    emulsion = data["layers"][1]
    for key, value in (grain or {}).items():
        if key == "morphology":
            emulsion["grain"]["morphology"] = value
        else:
            node = emulsion["grain"][key]
            node["value"] = value
            _strip_bounds(node)
    for key, value in (development or {}).items():
        node = emulsion["development"][key]
        node["value"] = value
        _strip_bounds(node)
    return Material.model_validate(data)


def class_capacity(mat: Material, k: int) -> np.ndarray:
    """``1e12 m_Ag,c nu_c`` (g/m2) of the ``k`` classes of ``mat``."""
    classes = size_classes(mat.emulsion_layers[0], k)
    return 1e12 * classes.silver_g_per_grain * classes.number_per_um2


def state(
    mat: Material,
    shape: tuple[int, int],
    pitch_um: float,
    *,
    k: int = 6,
    initiated,
    developed,
    developed_m2,
    physical=0.0,
    halide=0.0,
    growth_rate=None,
    residual_dye=1.0,
    seed: int = 0,
    film_format: str = "135",
    origin_px: tuple[int, int] = (0, 0),
    frame_shape_hw: tuple[int, int] | None = None,
) -> EmulsionState:
    """A developed state. ``initiated`` (F), ``developed`` (X) and ``developed_m2`` (Y)
    broadcast to ``(h, w, k)``; ``physical``, ``halide``, ``growth_rate`` and
    ``residual_dye`` to ``(h, w)``. ``residual_dye`` given as a float becomes a broadcast
    view, as the exposure model writes it."""
    h, w = shape
    cap = class_capacity(mat, k)
    f = np.broadcast_to(np.asarray(initiated, dtype=np.float64), (h, w, k)).copy()
    x = np.broadcast_to(np.asarray(developed, dtype=np.float64), (h, w, k)).copy()
    y = np.broadcast_to(np.asarray(developed_m2, dtype=np.float64), (h, w, k)).copy()
    fields = {
        "emulsion/silver": x * cap,
        "emulsion/initiated": f,
        "emulsion/silver_m2": y,
        "emulsion/silver_physical": np.broadcast_to(
            np.asarray(physical, dtype=np.float64), (h, w)
        ).copy(),
        "emulsion/halide": np.broadcast_to(np.asarray(halide, dtype=np.float64), (h, w)).copy(),
    }
    if isinstance(residual_dye, float):
        one = np.full((1, 1), residual_dye)
        one.flags.writeable = False
        fields["emulsion/residual_dye"] = np.broadcast_to(one, (h, w))
    else:
        fields["emulsion/residual_dye"] = np.asarray(residual_dye, dtype=np.float64)
    if growth_rate is not None:
        fields["emulsion/growth_rate"] = np.broadcast_to(
            np.asarray(growth_rate, dtype=np.float64), (h, w)
        ).copy()
    geometry = FrameGeometry((h, w), pitch_um, origin_px, frame_shape_hw)
    return EmulsionState(
        mat,
        geometry,
        DEFAULT_GRID,
        seed,
        fields,
        meta={"exposure": {"film_format": film_format}},
    )


def initiated_for_density(mat: Material, k: int, density: float, x_bar: float) -> float:
    """Uniform F giving ``D_Ag = density`` with ``X = x_bar F`` in every class, no physical silver."""
    from darkroom_atlas.imaging.silver import covering_powers

    classes = size_classes(mat.emulsion_layers[0], k)
    cp = covering_powers(classes, mat.emulsion_layers[0].development)
    per_unit = float(np.sum(cp * class_capacity(mat, k)))
    return density / (x_bar * per_unit)
