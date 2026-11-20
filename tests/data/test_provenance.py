import numpy as np
import pytest
from pydantic import ValidationError

from darkroom_atlas.data.provenance import (
    CurveOrigin,
    FormulaConfidence,
    Quantity,
    Source,
    VerificationStatus,
)


def test_enum_values():
    assert {s.value for s in Source} == {"datasheet", "literature", "fitted", "assumed"}
    assert {c.value for c in FormulaConfidence} == {"published", "patent", "msds-estimated"}
    assert {v.value for v in VerificationStatus} == {"validated", "predicted", "failed"}
    assert {o.value for o in CurveOrigin} == {"datasheet", "measured"}


def test_enum_member_names_both_cases():
    assert CurveOrigin.DATASHEET is CurveOrigin.datasheet is CurveOrigin("datasheet")
    assert FormulaConfidence.MSDS_ESTIMATED is FormulaConfidence.msds_estimated
    assert FormulaConfidence.msds_estimated == "msds-estimated"
    with pytest.raises(AttributeError):
        _ = Source.NOPE


def test_quantity_ok():
    q = Quantity(value=1.5, unit="um", source="assumed")
    assert q.source is Source.ASSUMED
    assert q.ref is None and q.uncertainty is None and q.note is None
    converted = q.to("nm")
    assert converted == 1500.0
    assert isinstance(converted, float)
    assert q.to("um") == 1.5


def test_quantity_list_to_array():
    q = Quantity(value=[1, 2], unit="mm", source="assumed")
    out = q.to("um")
    assert isinstance(out, np.ndarray)
    assert out.dtype == np.float64
    np.testing.assert_array_equal(out, [1000.0, 2000.0])
    same = q.to("mm")
    assert isinstance(same, np.ndarray)
    np.testing.assert_array_equal(same, [1.0, 2.0])


def test_quantity_to_rejects_wrong_dimension():
    with pytest.raises(ValueError):
        Quantity(value=1.0, unit="um", source="assumed").to("s")


def test_ref_required():
    with pytest.raises(ValidationError, match="ref"):
        Quantity(value=1.0, unit="1", source="datasheet")
    with pytest.raises(ValidationError, match="ref"):
        Quantity(value=1.0, unit="1", source="literature", ref="")
    with pytest.raises(ValidationError, match="ref"):
        Quantity(value=1.0, unit="1", source="literature", ref="   ")
    ok = Quantity(value=1.0, unit="1", source="datasheet", ref="test reference")
    assert ok.ref == "test reference"
    # fitted and assumed values need no reference
    Quantity(value=1.0, unit="1", source="fitted")
    Quantity(value=1.0, unit="1", source="assumed")


def test_negative_uncertainty():
    with pytest.raises(ValidationError):
        Quantity(value=1.0, unit="1", source="assumed", uncertainty=-1)
    assert Quantity(value=1.0, unit="1", source="assumed", uncertainty=0).uncertainty == 0.0


def test_extra_field_forbidden():
    with pytest.raises(ValidationError, match="foo"):
        Quantity(value=1.0, unit="um", source="assumed", foo=1)


def test_frozen():
    q = Quantity(value=1.0, unit="um", source="assumed")
    with pytest.raises(ValidationError):
        q.value = 2
    assert q.value == 1.0


@pytest.mark.parametrize(
    "value",
    [float("nan"), float("inf"), [1.0, float("nan")], [], True, "abc"],
)
def test_bad_values_rejected(value):
    with pytest.raises(ValidationError):
        Quantity(value=value, unit="1", source="assumed")


@pytest.mark.parametrize("unit", ["", "µm", "u m"])
def test_unit_must_be_plain_ascii(unit):
    with pytest.raises(ValidationError):
        Quantity(value=1.0, unit=unit, source="assumed")


def test_unknown_source_rejected():
    with pytest.raises(ValidationError):
        Quantity(value=1.0, unit="1", source="guessed")
