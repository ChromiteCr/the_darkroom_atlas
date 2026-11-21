"""Stop, fix, wash and dry for the ``bw-mechanistic`` model (stage 5, spec D1 2.6).

All four steps are pointwise first-order updates (``interaction_range_um = 0``):

* ``stop``: the developer is neutralised with time constant ``tau_s``, so the growth rate
  decays as ``kappa(t) = kappa_c,end exp(-t / tau_s)``; with
  ``K(t) = kappa_c,end tau_s (1 - exp(-t / tau_s))`` and ``F`` fixed (no new initiation),
  ``X = F - (F - X0) e^-K`` and ``Y = F - 2 (F - X0) e^-K + (F - 2 X0 + Y0) e^-2K``.
  ``kappa_c,end = (l_ref S_c / V_c) growth_rate``. Updates ``emulsion/silver`` and
  ``emulsion/silver_m2`` and drops ``emulsion/growth_rate``.
* ``fix``: ``halide <- halide exp(-t / (f_fix tau_fix))``, ``residual_dye <- residual_dye
  exp(-t / tau_dye)``, ``frame/fixer_residual <- 1``.
* ``wash``: ``fixer_residual <- fixer_residual exp(-t / tau_w)`` (when present),
  ``residual_dye <- residual_dye exp(-t / tau_dye,w)``.
* ``dry``: no change.

Each step appends one entry (its index, kind, bath, time and constants) to
``meta["post"]``, a tuple in step order. Uniform (broadcast) fields stay broadcast views.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from darkroom_atlas._frozen import freeze_mapping
from darkroom_atlas.data.models import BW_MECHANISTIC, StepKind
from darkroom_atlas.imaging.silver import imaging_constants
from darkroom_atlas.registry import register_step

if TYPE_CHECKING:
    from darkroom_atlas.data.models import ProcessStep
    from darkroom_atlas.pipeline import StepContext
    from darkroom_atlas.state import EmulsionState

__all__ = ["L_REF_UM", "MODEL", "dry_bw", "fix_bw", "stop_bw", "stop_update", "wash_bw"]

MODEL = BW_MECHANISTIC
"""Material model the steps are registered for (spec 2.0.4)."""

L_REF_UM = 1.0
"""Reference length of the specific-surface growth-rate scaling, um (spec 2.4.4)."""

_LAYER = "emulsion"


def _finite_array(name: str, value) -> np.ndarray:
    array = np.asarray(value, dtype=np.float64)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _positive(name: str, value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number, got {value!r}")
    if not value > 0.0:
        raise ValueError(f"{name} must be positive, got {value!r}")
    return float(value)


def stop_update(
    x: np.ndarray,
    y: np.ndarray,
    f: np.ndarray,
    kappa_end: np.ndarray,
    tau_s: float,
    time_s: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Developed fraction ``X`` and mean squared fraction ``Y`` after a stop of ``time_s``.

    Closed form of ``dX/dK = F - X``, ``dY/dK = 2 (X - Y)`` at
    ``K = kappa_end tau_s (1 - exp(-time_s / tau_s))`` (spec 2.6); ``F`` does not change.
    Arguments broadcast; ``kappa_end`` (s^-1) must be non-negative, ``tau_s`` positive and
    ``time_s`` non-negative. Returns new arrays ``(X, Y)``.
    """
    x = _finite_array("x", x)
    y = _finite_array("y", y)
    f = _finite_array("f", f)
    kappa = _finite_array("kappa_end", kappa_end)
    if np.any(kappa < 0.0):
        raise ValueError("kappa_end must not be negative")
    tau_s = _positive("tau_s", tau_s)
    if isinstance(time_s, bool) or not isinstance(time_s, (int, float)):
        raise ValueError(f"time_s must be a number, got {time_s!r}")  # noqa: TRY004
    if not math.isfinite(time_s) or time_s < 0.0:
        raise ValueError(f"time_s must be finite and non-negative, got {time_s!r}")
    k = kappa * (tau_s * -math.expm1(-float(time_s) / tau_s))
    e = np.exp(-k)
    undeveloped = f - x
    x_new = f - undeveloped * e
    y_new = f - 2.0 * undeveloped * e + (f - 2.0 * x + y) * (e * e)
    return x_new, y_new


def _bath_constant(ctx: StepContext, step: ProcessStep, name: str) -> float:
    bath = ctx.bath
    if bath is None:
        raise ValueError(
            f"a {step.kind.value!r} step of the bw-mechanistic model needs a bath with "
            f"constants.{name}; the step names no bath"
        )
    quantity = None if bath.constants is None else getattr(bath.constants, name)
    if quantity is None:
        raise ValueError(
            f"bath {bath.id!r} has no constants.{name}, which the bw-mechanistic "
            f"{step.kind.value!r} step needs"
        )
    value = quantity.to("s")
    if isinstance(value, np.ndarray):
        raise ValueError(f"bath {bath.id!r}: constants.{name} must be a single value")  # noqa: TRY004
    return _positive(f"bath {bath.id!r} constants.{name}", value)


def _time(step: ProcessStep) -> float:
    t = step.time_s
    if t is None:
        raise ValueError(f"a {step.kind.value!r} step needs a time")
    return float(t)


def _read_only(array: np.ndarray) -> np.ndarray:
    array = np.ascontiguousarray(array, dtype=np.float64)
    array.flags.writeable = False
    return array


def _uniform(value: float, shape_hw: tuple[int, int]) -> np.ndarray:
    one = np.full((1, 1), value, dtype=np.float64)
    one.flags.writeable = False
    return np.broadcast_to(one, shape_hw)


def _scaled(field: np.ndarray, factor: float) -> np.ndarray:
    """``field * factor``; a broadcast (uniform) field stays a broadcast view."""
    if field.ndim >= 2 and field.strides[:2] == (0, 0):
        one = _finite_array("field", field[:1, :1]) * factor
        one.flags.writeable = False
        return np.broadcast_to(one, field.shape)
    return _read_only(_finite_array("field", field) * factor)


def _with_post_meta(
    state: EmulsionState, ctx: StepContext, step: ProcessStep, constants: dict
) -> EmulsionState:
    entry = {
        "index": ctx.step_index,
        "kind": step.kind.value,
        "bath": step.bath,
        "time_s": step.time_s,
        **constants,
    }
    previous = tuple(state.meta.get("post", ()))
    return state.with_meta(post=(*previous, freeze_mapping(entry)))


@register_step(MODEL, StepKind.STOP, interaction_range_um=0.0)
def stop_bw(state: EmulsionState, step: ProcessStep, ctx: StepContext) -> EmulsionState:
    """Stop bath: development continues while the developer is neutralised (spec 2.6)."""
    tau_s = _bath_constant(ctx, step, "neutralisation_time")
    time_s = _time(step)
    const = imaging_constants(state)
    classes = const.classes
    silver = _finite_array("emulsion/silver", state.get(_LAYER, "silver"))
    f = _finite_array("emulsion/initiated", state.get(_LAYER, "initiated"))
    y0 = _finite_array("emulsion/silver_m2", state.get(_LAYER, "silver_m2"))
    growth = _finite_array("emulsion/growth_rate", state.get(_LAYER, "growth_rate"))
    if f.shape != silver.shape or y0.shape != silver.shape:
        raise ValueError(
            f"emulsion/initiated {f.shape} and emulsion/silver_m2 {y0.shape} must have the "
            f"shape of emulsion/silver {silver.shape}"
        )
    if growth.shape != silver.shape[:2]:
        raise ValueError(f"emulsion/growth_rate must have shape (H, W), got {growth.shape}")
    if np.any(growth < 0.0):
        raise ValueError("emulsion/growth_rate must not be negative")
    capacity = const.class_silver_g_m2
    specific = L_REF_UM * np.asarray(classes.surface_um2) / np.asarray(classes.volume_um3)
    kappa_end = growth[..., None] * specific
    x0 = silver / capacity
    x, y = stop_update(x0, y0, f, kappa_end, tau_s, time_s)
    out = state.with_fields(
        {"emulsion/silver": _read_only(x * capacity), "emulsion/silver_m2": _read_only(y)}
    ).drop("emulsion/growth_rate")
    return _with_post_meta(out, ctx, step, {"neutralisation_time_s": tau_s})


@register_step(MODEL, StepKind.FIX, interaction_range_um=0.0)
def fix_bw(state: EmulsionState, step: ProcessStep, ctx: StepContext) -> EmulsionState:
    """Fixer: remaining halide and residual dye decay; the film takes up fixer (spec 2.6)."""
    tau_fix = _bath_constant(ctx, step, "clearing_time_constant")
    tau_dye = _bath_constant(ctx, step, "dye_removal_time")
    time_s = _time(step)
    layer = next(layer for layer in state.material.emulsion_layers if layer.id == _LAYER)
    if layer.development is None:
        raise ValueError(f"material {state.material.id!r}: emulsion layer needs development")
    f_fix = layer.development.fixing_factor.to("1")
    if isinstance(f_fix, np.ndarray):
        raise ValueError("development.fixing_factor must be a single value")  # noqa: TRY004
    f_fix = _positive("development.fixing_factor", f_fix)
    halide = _finite_array("emulsion/halide", state.get(_LAYER, "halide"))
    residual = state.get(_LAYER, "residual_dye")
    out = state.with_fields(
        {
            "emulsion/halide": _read_only(halide * math.exp(-time_s / (f_fix * tau_fix))),
            "emulsion/residual_dye": _scaled(residual, math.exp(-time_s / tau_dye)),
            "frame/fixer_residual": _uniform(1.0, state.geometry.shape_hw),
        }
    )
    constants = {
        "clearing_time_constant_s": tau_fix,
        "dye_removal_time_s": tau_dye,
        "fixing_factor": f_fix,
    }
    return _with_post_meta(out, ctx, step, constants)


@register_step(MODEL, StepKind.WASH, interaction_range_um=0.0)
def wash_bw(state: EmulsionState, step: ProcessStep, ctx: StepContext) -> EmulsionState:
    """Wash: residual fixer and residual dye decay (spec 2.6)."""
    tau_w = _bath_constant(ctx, step, "wash_time_constant")
    tau_dye_w = _bath_constant(ctx, step, "wash_dye_time_constant")
    time_s = _time(step)
    updates = {
        "emulsion/residual_dye": _scaled(
            state.get(_LAYER, "residual_dye"), math.exp(-time_s / tau_dye_w)
        )
    }
    if state.has("frame", "fixer_residual"):
        updates["frame/fixer_residual"] = _scaled(
            state.get("frame", "fixer_residual"), math.exp(-time_s / tau_w)
        )
    constants = {"wash_time_constant_s": tau_w, "wash_dye_time_constant_s": tau_dye_w}
    return _with_post_meta(state.with_fields(updates), ctx, step, constants)


@register_step(MODEL, StepKind.DRY, interaction_range_um=0.0)
def dry_bw(state: EmulsionState, step: ProcessStep, ctx: StepContext) -> EmulsionState:
    """Dry: no change to the fields; the step is recorded in ``meta["post"]`` (spec 2.6)."""
    return _with_post_meta(state, ctx, step, {})
