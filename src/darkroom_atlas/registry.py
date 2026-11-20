"""Registries of exposure models, process-step handlers and readout devices.

Physics plugs into the pipeline here. A material's ``model`` field selects the exposure model
and, together with each step's ``kind``, the step handler; readout devices are registered by
name. Registration is done with decorators when a model module is imported; importing
:mod:`darkroom_atlas.models` registers the built-in models, and the getters below make sure
that has happened.

Calling conventions (see :mod:`darkroom_atlas.pipeline` for the context objects):

* exposure model: ``fn(material, light, time_s, ctx: RunContext) -> EmulsionState``
* step handler: ``fn(state, step: ProcessStep, ctx: StepContext) -> EmulsionState``
* readout device: ``fn(state, ctx: RunContext, **params) -> Readout``

A step handler declares its lateral reach on the film, ``interaction_range_um``:

* ``0.0``: pointwise; each pixel depends only on itself.
* a positive number: pixels interact within this distance (um). The pipeline may run the
  step tile by tile with an overlap of at least this range, so the result is seamless.
* ``None``: global coupling (for example a finite volume of developer shared by the whole
  frame). The pipeline never tiles such a step; the handler receives the whole frame and
  tiles whatever local parts it has itself, with :func:`darkroom_atlas.tiling.map_tiles`.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from numbers import Real
from typing import TYPE_CHECKING, Any, Protocol

from darkroom_atlas.data.models import StepKind

if TYPE_CHECKING:
    from darkroom_atlas.data.models import Material, ProcessStep
    from darkroom_atlas.light import Light
    from darkroom_atlas.pipeline import Readout, RunContext, StepContext
    from darkroom_atlas.state import EmulsionState

__all__ = [
    "ExposureModel",
    "ReadoutDevice",
    "RegisteredStep",
    "RegistryError",
    "StepHandler",
    "get_exposure",
    "get_readout",
    "get_step",
    "register_exposure",
    "register_readout",
    "register_step",
    "temporary_registrations",
]


class RegistryError(LookupError):
    """Nothing is registered under the requested name. The message lists what is."""


class ExposureModel(Protocol):
    """Turns light on the film into the latent state of a material."""

    def __call__(
        self, material: Material, light: Light, time_s: float, ctx: RunContext
    ) -> EmulsionState: ...


class StepHandler(Protocol):
    """Applies one process step to a state."""

    def __call__(
        self, state: EmulsionState, step: ProcessStep, ctx: StepContext
    ) -> EmulsionState: ...


class ReadoutDevice(Protocol):
    """Reads a processed state out as an image, spectra or densities."""

    def __call__(self, state: EmulsionState, ctx: RunContext, **params: Any) -> Readout: ...


@dataclass(frozen=True)
class RegisteredStep:
    """A step handler and its lateral interaction range in um (``None``: global)."""

    fn: StepHandler
    interaction_range_um: float | None


_EXPOSURES: dict[str, ExposureModel] = {}
_STEPS: dict[tuple[str, StepKind], RegisteredStep] = {}
_READOUTS: dict[str, ReadoutDevice] = {}


# --------------------------------------------------------------------------------------
# Registration
# --------------------------------------------------------------------------------------


def register_exposure(model: str) -> Callable[[ExposureModel], ExposureModel]:
    """Decorator registering an exposure model under ``model``; the function is returned as is.

    Registering a model name twice raises ``ValueError``.
    """
    model = _name("model", model)

    def decorator(fn: ExposureModel) -> ExposureModel:
        _add(_EXPOSURES, model, fn, f"exposure model {model!r}")
        return fn

    return decorator


def register_step(
    model: str, kind: StepKind | str, *, interaction_range_um: float | None = 0.0
) -> Callable[[StepHandler], StepHandler]:
    """Decorator registering the handler of ``kind`` steps for materials of ``model``.

    ``interaction_range_um`` is ``0.0`` (pointwise, the default), a positive lateral range in
    um, or ``None`` for a step with global coupling; see the module docs. Registering the
    same ``(model, kind)`` twice raises ``ValueError``.
    """
    model = _name("model", model)
    kind = _kind(kind)
    if interaction_range_um is not None:
        r = interaction_range_um
        if isinstance(r, bool) or not isinstance(r, Real) or not math.isfinite(r) or r < 0.0:
            raise ValueError(
                f"interaction_range_um must be None or a finite number >= 0, got {r!r}"
            )
        interaction_range_um = float(r)

    def decorator(fn: StepHandler) -> StepHandler:
        entry = RegisteredStep(fn=fn, interaction_range_um=interaction_range_um)
        _add(_STEPS, (model, kind), entry, f"{kind.value!r} step for model {model!r}")
        return fn

    return decorator


def register_readout(name: str) -> Callable[[ReadoutDevice], ReadoutDevice]:
    """Decorator registering a readout device under ``name``; the function is returned as is.

    Registering a name twice raises ``ValueError``.
    """
    name = _name("readout device name", name)

    def decorator(fn: ReadoutDevice) -> ReadoutDevice:
        _add(_READOUTS, name, fn, f"readout device {name!r}")
        return fn

    return decorator


@contextmanager
def temporary_registrations() -> Iterator[None]:
    """Context manager that undoes, on exit, every registration made inside it.

    Meant for tests that register throw-away models. The built-in models are imported
    first, so that registering them for the first time inside the block cannot undo them.
    """
    _load_builtin_models()
    saved = (dict(_EXPOSURES), dict(_STEPS), dict(_READOUTS))
    try:
        yield
    finally:
        for table, contents in zip((_EXPOSURES, _STEPS, _READOUTS), saved, strict=True):
            table.clear()
            table.update(contents)


# --------------------------------------------------------------------------------------
# Lookup
# --------------------------------------------------------------------------------------


def get_exposure(model: str) -> ExposureModel:
    """The exposure model registered as ``model``; :class:`RegistryError` if there is none."""
    _load_builtin_models()
    try:
        return _EXPOSURES[model]
    except (KeyError, TypeError):
        raise RegistryError(
            f"no exposure model {model!r} is registered; registered: {_listing(_EXPOSURES)}"
        ) from None


def get_step(model: str, kind: StepKind | str) -> RegisteredStep:
    """The handler of ``kind`` steps for ``model``; :class:`RegistryError` if there is none."""
    _load_builtin_models()
    try:
        return _STEPS[(model, StepKind(kind))]
    except (KeyError, TypeError, ValueError):
        for_model = sorted(k.value for m, k in _STEPS if m == model)
        registered = ", ".join(f"{m}:{k.value}" for m, k in sorted(_STEPS)) or "(none)"
        raise RegistryError(
            f"no {str(kind)!r} step is registered for model {model!r} (its steps: "
            f"{', '.join(for_model) or 'none'}); registered model:step pairs: {registered}"
        ) from None


def get_readout(name: str) -> ReadoutDevice:
    """The readout device registered as ``name``; :class:`RegistryError` if there is none."""
    _load_builtin_models()
    try:
        return _READOUTS[name]
    except (KeyError, TypeError):
        raise RegistryError(
            f"no readout device {name!r} is registered; registered: {_listing(_READOUTS)}"
        ) from None


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------


def _load_builtin_models() -> None:
    # Imported here, not at the top: the model modules import this one.
    import darkroom_atlas.models  # noqa: F401


def _add(table: dict, key: Any, value: Any, what: str) -> None:
    if key in table:
        raise ValueError(f"{what} is already registered")
    table[key] = value


def _name(what: str, value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string, got {value!r}")
    return value


def _kind(kind: object) -> StepKind:
    try:
        return StepKind(kind)
    except ValueError:
        allowed = ", ".join(k.value for k in StepKind)
        raise ValueError(f"unknown step kind {kind!r}; expected one of: {allowed}") from None


def _listing(table: dict[str, Any]) -> str:
    return ", ".join(sorted(table)) or "(none)"
