import pytest

import darkroom_atlas.pipeline  # noqa: F401  (importing the pipeline registers the built-ins)
from darkroom_atlas import registry
from darkroom_atlas.data.models import StepKind
from darkroom_atlas.registry import (
    RegisteredStep,
    RegistryError,
    get_exposure,
    get_readout,
    get_step,
    register_exposure,
    register_readout,
    register_step,
)


@pytest.fixture
def scratch():
    """Registrations made in a test are undone afterwards, even if the test fails."""
    with registry.temporary_registrations():
        yield


def _noop(*args, **kwargs):
    raise AssertionError("not meant to be called")


def test_ideal_registered():
    assert callable(get_exposure("ideal"))
    step = get_step("ideal", StepKind.identity)
    assert isinstance(step, RegisteredStep)
    assert callable(step.fn)
    assert step.interaction_range_um == 0.0
    assert get_step("ideal", "identity") is step  # the kind may be given as its string value
    assert callable(get_readout("ideal_colorimetric"))


def test_unknown_lists_registered():
    with pytest.raises(RegistryError, match="ideal") as info:
        get_exposure("nope")
    assert "nope" in str(info.value)
    assert isinstance(info.value, LookupError)
    with pytest.raises(RegistryError, match="identity"):
        get_step("ideal", StepKind.develop)
    with pytest.raises(RegistryError, match="ideal"):
        get_step("nope", StepKind.identity)
    with pytest.raises(RegistryError, match="ideal_colorimetric"):
        get_readout("nope")


def test_duplicate_rejected(scratch):
    register_exposure("duptest")(_noop)
    with pytest.raises(ValueError, match="duptest"):
        register_exposure("duptest")(_noop)
    register_step("duptest", StepKind.identity)(_noop)
    with pytest.raises(ValueError, match="duptest"):
        register_step("duptest", "identity")(_noop)
    register_readout("duptest")(_noop)
    with pytest.raises(ValueError, match="duptest"):
        register_readout("duptest")(_noop)
    with pytest.raises(ValueError):
        register_exposure("ideal")(_noop)  # the built-ins are protected too


def test_step_range_stored(scratch):
    fn = register_step("rangetest", StepKind.develop, interaction_range_um=12.0)(_noop)
    assert fn is _noop  # the decorator returns the function unchanged
    step = get_step("rangetest", StepKind.develop)
    assert step.interaction_range_um == 12.0
    assert step.fn is _noop
    register_step("rangetest", StepKind.wash, interaction_range_um=None)(_noop)
    assert get_step("rangetest", StepKind.wash).interaction_range_um is None


def test_register_rejects_bad_arguments(scratch):
    with pytest.raises(ValueError):
        register_step("badtest", "nope")
    with pytest.raises(ValueError):
        register_step("badtest", StepKind.identity, interaction_range_um=-1.0)
    with pytest.raises(ValueError):
        register_step("badtest", StepKind.identity, interaction_range_um=float("nan"))
    with pytest.raises(ValueError):
        register_exposure("")
    with pytest.raises(ValueError):
        register_readout("")


def test_temporary_registrations_cleaned_up():
    with registry.temporary_registrations():
        register_exposure("cleanuptest")(_noop)
        register_step("cleanuptest", StepKind.identity)(_noop)
        register_readout("cleanuptest")(_noop)
        assert get_exposure("cleanuptest") is _noop
    with pytest.raises(RegistryError):
        get_exposure("cleanuptest")
    with pytest.raises(RegistryError):
        get_step("cleanuptest", StepKind.identity)
    with pytest.raises(RegistryError):
        get_readout("cleanuptest")
    assert callable(get_exposure("ideal"))
    assert callable(get_readout("ideal_colorimetric"))


def test_temporary_registrations_cleaned_up_on_error():
    with pytest.raises(RuntimeError), registry.temporary_registrations():
        register_exposure("errortest")(_noop)
        raise RuntimeError
    with pytest.raises(RegistryError):
        get_exposure("errortest")
