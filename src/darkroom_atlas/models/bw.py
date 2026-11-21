"""The ``bw-mechanistic`` model: the D1 black-and-white chain (spec D1 section 2).

D1a provides its exposure model (spec §2.1.9). From :func:`~darkroom_atlas.pipeline.
unexposed_state` it writes

* ``emulsion/exposure_volume`` ``(H, W, n_z)``: ``e_V,z``, photons absorbed per um^3 of AgX;
* ``emulsion/exposure_surface`` ``(H, W, n_z)``: ``e_S,z``, photons absorbed by the dye per
  um^2 of grain surface;
* ``emulsion/residual_dye``: 1.0 everywhere, as a zero-byte broadcast view;

and ``meta["exposure"]`` (spec §2.1.9). ``n_z`` and the number of size classes follow
``RunContext.profile`` (spec §3.1); the film format is ``RunContext.film_format`` or the
material's ``default_format``. Large fields are memory-mapped under
``RunContext.scratch_dir`` (spec §4.3).
"""

import uuid

import numpy as np

from darkroom_atlas.data.library import DataError, DataNotFound
from darkroom_atlas.data.models import BW_MECHANISTIC, Material, bw_mechanistic_problems
from darkroom_atlas.develop.grid import PROFILE_LEVELS
from darkroom_atlas.emulsion.classes import size_classes
from darkroom_atlas.fieldstore import FieldStore
from darkroom_atlas.light import LightField
from darkroom_atlas.optics.expose import BLOCK_PX, exposure_fields
from darkroom_atlas.optics.kubelka import depth_response
from darkroom_atlas.optics.stack import optical_stack
from darkroom_atlas.pipeline import RunContext, unexposed_state
from darkroom_atlas.registry import register_exposure
from darkroom_atlas.state import EmulsionState

__all__ = [
    "EXPOSURE_TEMPERATURE_C",
    "MODEL",
    "PROFILE_CLASSES",
    "check_material",
    "expose_bw",
    "film_format_of",
]

MODEL = BW_MECHANISTIC
"""Name of the model (spec §2.0.4)."""

EXPOSURE_TEMPERATURE_C = 23.0
"""Exposure environment of ISO 6 §5.3.1, deg C (spec §2.1.9; CAL §1.1)."""

PROFILE_CLASSES: dict[str, int] = {"reference": 12, "calibration": 6, "image": 6}
"""Number of grain size classes ``K`` per discretisation profile (spec §3.1)."""


def check_material(material: Material, ctx: RunContext) -> None:
    """Raise :class:`DataError` if ``material`` lacks a field the model needs (spec §5).

    The message starts with the material's file (when the run's library has it) and lists
    one ``<field path>: <problem>`` line per problem.
    """
    problems = bw_mechanistic_problems(material)
    if not problems:
        return
    try:
        where = str(ctx.library.path_of("material", material.id))
    except DataNotFound:
        where = f"material {material.id!r}"
    raise DataError("\n".join(f"{where}: {line}" for line in problems))


def film_format_of(material: Material, ctx: RunContext) -> str:
    """The run's film format: ``ctx.film_format``, else the material's default format."""
    film_format = ctx.film_format if ctx.film_format is not None else material.default_format
    formats = material.formats
    if film_format is None or film_format not in formats:
        raise ValueError(
            f"material {material.id!r}: unknown film format {film_format!r}; available: "
            f"{', '.join(formats) or '(none)'}"
        )
    return film_format


@register_exposure(MODEL)
def expose_bw(
    material: Material, light: LightField, time_s: float, ctx: RunContext
) -> EmulsionState:
    """Expose a ``bw-mechanistic`` film: light transport in the dry emulsion (spec §2.1)."""
    check_material(material, ctx)
    film_format = film_format_of(material, ctx)
    state = unexposed_state(material, light, ctx)
    (layer,) = material.emulsion_layers
    n_z = PROFILE_LEVELS[ctx.profile][0]
    classes = size_classes(layer, PROFILE_CLASSES[ctx.profile])
    stack = optical_stack(material, film_format, classes, ctx.grid)
    response = depth_response(stack, n_z)
    store = (
        None
        if ctx.scratch_dir is None
        else FieldStore(ctx.scratch_dir / f"expose-{uuid.uuid4().hex[:12]}")
    )
    # The store is not closed: the state refers to its files, which the runner deletes
    # with the scratch directory.
    volume, surface, optics_meta = exposure_fields(
        light,
        time_s,
        stack,
        response,
        state.geometry,
        block_px=BLOCK_PX,
        workers=ctx.workers,
        store=store,
    )
    one = np.ones((1, 1))
    one.flags.writeable = False
    residual_dye = np.broadcast_to(one, state.geometry.shape_hw)
    meta = {
        "time_s": float(time_s),
        "temperature_c": EXPOSURE_TEMPERATURE_C,
        "n_z": n_z,
        "profile": ctx.profile,
        "film_format": film_format,
        **optics_meta,
    }
    return state.with_fields(
        {
            "emulsion/exposure_volume": volume,
            "emulsion/exposure_surface": surface,
            "emulsion/residual_dye": residual_dye,
        }
    ).with_meta(exposure=meta)
