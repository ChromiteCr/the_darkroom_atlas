"""Conditions and the ``# conditions:`` header of the reference CSVs (spec D1 6.1, D1e)."""

import pytest

from darkroom_atlas.data.models import Agitation, Tank
from darkroom_atlas.validate.conditions import (
    Condition,
    fahrenheit_to_c,
    parse_conditions_header,
    parse_temperature_c,
)


def test_conditions_header_parsing(reference_kodak_dir):
    files = sorted(reference_kodak_dir.glob("*/*.csv"))
    assert len(files) >= 50
    for path in files:
        text = path.read_text(encoding="utf-8")
        fields = parse_conditions_header(text)
        assert fields, path.name
        assert all(value for value in fields.values()), path.name
        assert "densitometry" in fields, path.name
        if path.name.startswith(("char-", "ci-", "mtf-")):
            assert fields["developer"], path.name
            assert "agitation" in fields, path.name
            temperature = parse_temperature_c(fields["temperature"])
            assert temperature is None or temperature == pytest.approx(20.0), path.name
        if path.name.startswith("char-"):
            assert fields["time"].endswith("min"), path.name
            assert "format" in fields and "dilution" in fields, path.name
        # Values that contain ";" or "=" inside remarks stay in one field.
        if "legend_check" in fields:
            assert "plot_floor" in fields
            assert "plot_floor" not in fields["legend_check"]


def test_conditions_header_exact_fields():
    line = (
        "# conditions: developer=d-76; dilution=stock (the figure names no dilution); "
        "time=8min; temperature=20 C; agitation=large tank, agitation at 1 minute intervals; "
        "note=a; b=c remark"
    )
    assert parse_conditions_header(line) == {
        "developer": "d-76",
        "dilution": "stock (the figure names no dilution)",
        "time": "8min",
        "temperature": "20 C",
        "agitation": "large tank, agitation at 1 minute intervals",
        "note": "a",
        "b": "c remark",
    }


def test_conditions_header_errors():
    with pytest.raises(ValueError, match="no '# conditions:'"):
        parse_conditions_header("# material: kodak/t-max-100\n1,2\n")
    with pytest.raises(ValueError, match="empty"):
        parse_conditions_header("# conditions:   ")
    with pytest.raises(ValueError, match="twice"):
        parse_conditions_header("# conditions: a=1; a=2")
    with pytest.raises(ValueError, match="key=value"):
        parse_conditions_header("# conditions: Developer D-76")


@pytest.mark.parametrize(
    ("f", "c"),
    [(65, 18.33), (68, 20.00), (70, 21.11), (72, 22.22), (75, 23.89), (80, 26.67)],
)
def test_fahrenheit_conversion_is_exact(f, c):
    # Spec 6.1: Kodak's °F converted exactly, not the rounded °C labels.
    assert fahrenheit_to_c(f) == (f - 32) * 5 / 9
    assert fahrenheit_to_c(f) == pytest.approx(c, abs=0.005)
    assert parse_temperature_c(f"{f} F") == (f - 32) * 5 / 9
    assert parse_temperature_c(f"{f} F ({c:.0f} C)") == (f - 32) * 5 / 9


def test_parse_temperature_c():
    assert parse_temperature_c("20 C") == 20.0
    assert parse_temperature_c("not stated") is None
    with pytest.raises(ValueError):
        parse_temperature_c("room temperature")


def _kwargs(**over):
    base = {
        "material_id": "kodak/t-max-100",
        "film_format": "135",
        "developer": "d-76",
        "dilution": "stock",
        "temperature_c": 20.0,
        "agitation": Agitation(method="inversion", interval_s=30.0, duration_s=5.0),
        "times_s": (360.0, 450.0, 600.0),
        "tank": Tank(volume_ml=473.0),
    }
    return base | over


def test_condition_defaults_and_validation():
    cond = Condition(**_kwargs(times_s=[360, 450]))
    assert cond.times_s == (360.0, 450.0)
    assert cond.exposure_time_s == 0.02
    assert cond.illuminant == "iso_daylight"
    assert hash(cond) == hash(Condition(**_kwargs(times_s=(360.0, 450.0))))
    with pytest.raises(ValueError, match="increasing"):
        Condition(**_kwargs(times_s=(450.0, 360.0)))
    with pytest.raises(ValueError, match="positive"):
        Condition(**_kwargs(times_s=(0.0, 360.0)))
    with pytest.raises(ValueError, match="at least one"):
        Condition(**_kwargs(times_s=()))
    with pytest.raises(ValueError, match="exposure_time_s"):
        Condition(**_kwargs(exposure_time_s=0.0))
    with pytest.raises(ValueError, match="illuminant"):
        Condition(**_kwargs(illuminant="D65"))
    with pytest.raises(ValueError, match="temperature_c"):
        Condition(**_kwargs(temperature_c=float("nan")))


@pytest.mark.xfail(
    strict=True,
    reason="spec 10 D1e: header -> Condition mapping (agitation, tank) is conditions_for, "
    "a later task",
)
def test_conditions_header_parsing_to_condition(reference_kodak_dir):
    # Spec 10 D1e: every staged char-/ci- CSV header becomes a Condition, with Kodak's
    # Fahrenheit temperatures converted exactly (spec 6.1).
    from darkroom_atlas.data.library import Library
    from darkroom_atlas.validate.conditions import conditions_for

    library = Library()
    for film_dir in sorted(p for p in reference_kodak_dir.iterdir() if p.is_dir()):
        material = f"kodak/{film_dir.name}"
        conditions = [
            c
            for split in ("train", "within", "held_out")
            for c in conditions_for(material, library, split=split)
        ]
        sources = " ".join(c.source for c in conditions)
        for path in sorted(film_dir.glob("*.csv")):
            if path.name.startswith(("char-", "ci-")):
                assert path.name in sources, path.name
        for c in conditions:
            assert isinstance(c, Condition)
            assert (
                c.temperature_c
                == pytest.approx(round(c.temperature_c * 1.8 + 32.0) / 1.8 - 32.0 / 1.8, abs=1e-9)
                or float(c.temperature_c).is_integer()
            )
