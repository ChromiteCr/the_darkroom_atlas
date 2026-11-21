"""Recipes: one YAML file per figure, saying exactly how its curves were read.

A recipe names the publication (number, edition, URL, sha256 of the PDF), the page and the
figure region; how each axis is calibrated (which labelled ticks, which values, linear or
log10, least squares or piecewise); how each curve's paths are selected (stroke width and
dash array, checked against the legend sample next to the curve's legend text); and the
output file names. :func:`run_recipe` applies a recipe to a page and returns the curves;
it never writes files (see :mod:`tools.digitize.runner`).
"""

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated, Literal, Self

import numpy as np
import yaml
from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator

from tools.digitize.calibrate import Axis, Reference
from tools.digitize.figure import (
    AssembleStats,
    DigitizeError,
    assemble,
    curve_candidates,
    find_phrase,
    has_bar,
    label_words,
    legend_sample,
    match_style,
    tick_candidates,
    to_data,
)
from tools.digitize.page import BBox, Page
from tools.digitize.raster import mask_boxes, trace_curve

__all__ = [
    "CurveResult",
    "FigureResult",
    "Recipe",
    "csv_text",
    "load_recipe",
    "run_recipe",
]

FILE_NAME = r"^[A-Za-z0-9][A-Za-z0-9._+-]*$"


class _Strict(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")


Box = Annotated[list[float], Field(min_length=4, max_length=4)]


class Publication(_Strict):
    """The source document. ``file`` is the PDF's name under ``.cache/datasheets/pdf/``."""

    publisher: str
    number: str
    title: str
    edition: str
    url: str
    file: Annotated[str, Field(pattern=FILE_NAME)]
    sha256: Annotated[str, Field(pattern=r"^[0-9a-f]{64}$")]


class TickRef(_Strict):
    """One calibration point.

    Normally ``label`` (the tick label as the text layer spells it) and ``value`` (what it
    means). With ``bar: true`` the label is in bar notation and must carry the bar. Only
    when there is no usable label, ``coord`` gives the page coordinate directly and
    ``note`` says where it was read.
    """

    value: float
    label: str | None = None
    bar: bool = False
    coord: float | None = None
    note: str | None = None

    @model_validator(mode="after")
    def _check(self) -> Self:
        if (self.label is None) == (self.coord is None):
            raise ValueError("give exactly one of label or coord")
        if self.coord is not None and not self.note:
            raise ValueError("a tick given by coord needs a note saying where it was read")
        return self


class AxisSpec(_Strict):
    """How one axis is calibrated. ``title`` and ``unit`` are as printed on the figure."""

    title: str
    unit: str
    side: Literal["bottom", "top", "left", "right"]
    scale: Literal["linear", "log10"] = "linear"
    fit: Literal["lsq", "piecewise"] = "lsq"
    relative: bool = False
    bar_notation: bool = False
    label_band_pt: float = 24.0
    snap_pt: float = 4.0
    ticks: list[TickRef] = Field(min_length=2)
    note: str | None = None

    @model_validator(mode="after")
    def _check(self) -> Self:
        if self.relative and "relative" not in self.unit:
            raise ValueError("a relative axis must say 'relative' in its unit")
        return self


class Style(_Strict):
    """Stroke style of a curve: width and dash array in page points, optional colour."""

    width: float
    dash: list[float] = Field(default_factory=list)
    stroke_color: list[float] | None = None


class Expect(_Strict):
    """Checks on the result. Spans are ``[first x, last x]`` in data units."""

    paths: int = 1
    x_span: Annotated[list[float], Field(min_length=2, max_length=2)] | None = None
    y_span: Annotated[list[float], Field(min_length=2, max_length=2)] | None = None
    tolerance: Annotated[list[float], Field(min_length=2, max_length=2)] = Field(
        default_factory=lambda: [0.05, 0.05]
    )


class RasterSpec(_Strict):
    """Raster tracing of one curve. Points and boxes are in page points."""

    seed: Annotated[list[float], Field(min_length=2, max_length=2)]
    threshold: int = 128
    max_jump_pt: float = 1.0
    max_gap_pt: float = 4.0
    max_height_pt: float = 6.0
    mask: list[Box] = Field(default_factory=list)


class CurveSpec(_Strict):
    """One output curve."""

    output: Annotated[str, Field(pattern=FILE_NAME)]
    quantity: str
    legend: str | None = None
    style: Style | None = None
    raster: RasterSpec | None = None
    monotonic: Literal["strict", "drop-backtracking"] = "strict"
    monotonic_reason: str | None = None
    expect: Expect = Field(default_factory=Expect)
    conditions: dict[str, str]
    note: str | None = None

    @model_validator(mode="after")
    def _check(self) -> Self:
        if (self.style is None) == (self.raster is None):
            raise ValueError("give exactly one of style (vector) or raster")
        if self.monotonic == "drop-backtracking" and not self.monotonic_reason:
            raise ValueError("monotonic: drop-backtracking needs a monotonic_reason")
        for key, value in self.conditions.items():
            if ";" in value or "=" in key or "\n" in value:
                raise ValueError(f"conditions[{key!r}] must not contain ';', '=' or newlines")
        return self


class Recipe(_Strict):
    """A whole recipe file."""

    schema_version: Literal[1]
    publication: Publication
    page: int = Field(ge=1)
    printed_page: str
    figure: str
    caption: str
    material: str
    method: Literal["vector", "raster"] = "vector"
    raster_dpi: float = 600.0
    region: Box
    frame: Box
    exclude: list[Box] = Field(default_factory=list)
    flatten_tolerance_pt: float = Field(default=0.01, gt=0)
    decimals: Annotated[list[int], Field(min_length=2, max_length=2)] = Field(
        default_factory=lambda: [4, 4]
    )
    x: AxisSpec
    y: AxisSpec
    curves: list[CurveSpec] = Field(min_length=1)
    notes: list[str] = Field(default_factory=list)
    # A cross-check recipe reads the same figure in another edition. It stages nothing:
    # its curves are compared with the staged curves of the same output names, and the
    # differences must stay within max_difference = [x, y] (data units).
    cross_check_of: str | None = None
    max_difference: Annotated[list[float], Field(min_length=2, max_length=2)] | None = None

    @model_validator(mode="after")
    def _check(self) -> Self:
        outputs = [c.output for c in self.curves]
        if len(set(outputs)) != len(outputs):
            raise ValueError("curve outputs must be unique")
        if (self.cross_check_of is None) != (self.max_difference is None):
            raise ValueError("cross_check_of and max_difference go together")
        for curve in self.curves:
            if (self.method == "raster") != (curve.raster is not None):
                raise ValueError(f"curve {curve.output!r}: its kind does not match method")
        return self


def load_recipe(path: Path) -> Recipe:
    """Read and validate a recipe file."""
    try:
        data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        return Recipe.model_validate(data)
    except (yaml.YAMLError, ValidationError) as exc:
        raise DigitizeError(f"{path}: {exc}") from None


# --------------------------------------------------------------------------------------
# Running
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class CurveResult:
    """One digitised curve."""

    spec: CurveSpec
    x: np.ndarray
    y: np.ndarray
    error_x: float  # axis-variable units (decades on a log10 axis)
    error_y: float
    method: str  # "vector" or "raster"
    paths: tuple[int, ...]
    stats: AssembleStats
    stroke_pt: float


@dataclass(frozen=True)
class FigureResult:
    """Everything a recipe produced."""

    recipe: Recipe
    x_axis: Axis
    y_axis: Axis
    curves: tuple[CurveResult, ...]
    legend_paths: dict[str, int] = field(default_factory=dict)


def _calibrate(page: Page, recipe: Recipe, spec: AxisSpec) -> Axis:
    frame: BBox = tuple(recipe.frame)
    refs = []
    words = None
    candidates = None
    horizontal = spec.side in ("bottom", "top")
    for tick in spec.ticks:
        if tick.coord is not None:
            refs.append(Reference(tick.coord, tick.value, f"explicit coordinate: {tick.note}"))
            continue
        if words is None:
            words = label_words(page, frame, spec.side, spec.label_band_pt, spec.label_band_pt)
            candidates = tick_candidates(page, frame, spec.side)
        found = [w for w in words if w.text == tick.label]
        if spec.bar_notation:
            found = [w for w in found if has_bar(w, page) == tick.bar]
        elif tick.bar:
            raise DigitizeError(
                f"{spec.title}: tick {tick.label!r} has bar: true but the axis "
                "does not set bar_notation"
            )
        if len(found) != 1:
            bar = " with a bar" if tick.bar else (" without a bar" if spec.bar_notation else "")
            raise DigitizeError(
                f"{spec.title}: expected one label {tick.label!r}{bar} beside the "
                f"{spec.side} edge, found {len(found)}"
            )
        cx, cy = found[0].centre
        at = cx if horizontal else cy
        near = sorted(candidates, key=lambda c: abs(c[0] - at))
        if not near or abs(near[0][0] - at) > spec.snap_pt:
            raise DigitizeError(
                f"{spec.title}: no tick within {spec.snap_pt} pt of label {tick.label!r} "
                f"(label centre {at:.2f})"
            )
        if len(near) > 1 and abs(near[1][0] - at) < 1.5 * abs(near[0][0] - at) + 0.3:
            raise DigitizeError(
                f"{spec.title}: label {tick.label!r} at {at:.2f} is about as near the ticks at "
                f"{near[0][0]:.2f} and {near[1][0]:.2f}; give its coord explicitly"
            )
        coord, source = near[0]
        refs.append(Reference(coord, tick.value, f"label {tick.label!r} -> {source}"))
    try:
        return Axis(tuple(refs), scale=spec.scale, fit=spec.fit)
    except ValueError as exc:
        raise DigitizeError(f"{spec.title}: {exc}") from None


def _check_expect(curve: CurveSpec, xs: np.ndarray, ys: np.ndarray) -> None:
    tol_x, tol_y = curve.expect.tolerance
    for name, span, values, tol in (
        ("x", curve.expect.x_span, xs, tol_x),
        ("y", curve.expect.y_span, ys, tol_y),
    ):
        if span is None:
            continue
        got = (float(values[0]), float(values[-1]))
        if abs(got[0] - span[0]) > tol or abs(got[1] - span[1]) > tol:
            raise DigitizeError(
                f"curve {curve.output!r}: {name} runs from {got[0]:g} to {got[1]:g}, expected "
                f"{span[0]:g} to {span[1]:g} (tolerance {tol:g})"
            )


def _vector_curve(page, recipe, curve, x_axis, y_axis, skip, legend_paths) -> CurveResult:
    style = curve.style
    frame: BBox = tuple(recipe.frame)
    if curve.legend is not None:
        label = find_phrase(page, curve.legend, frame)
        sample = legend_sample(page, label)
        legend_paths[curve.output] = sample.index
        if not match_style(sample, style.width, style.dash, style.stroke_color):
            raise DigitizeError(
                f"curve {curve.output!r}: the legend sample for {curve.legend!r} "
                f"(path #{sample.index}, {sample.style()}) does not have the recipe's style "
                f"width={style.width} dash={style.dash}"
            )
    candidates = curve_candidates(page, frame, [tuple(b) for b in recipe.exclude], skip)
    chosen = [p for p in candidates if match_style(p, style.width, style.dash, style.stroke_color)]
    if len(chosen) != curve.expect.paths:
        listing = ", ".join(f"#{p.index}" for p in chosen) or "none"
        raise DigitizeError(
            f"curve {curve.output!r}: expected {curve.expect.paths} path(s) with width="
            f"{style.width} dash={style.dash}, found {len(chosen)} ({listing})"
        )
    tol = recipe.flatten_tolerance_pt
    polylines = [poly for p in chosen for poly in p.polylines(tol)]
    data = to_data(polylines, x_axis, y_axis)
    xs, ys, stats = assemble(data, tuple(recipe.decimals), curve.monotonic)
    stroke = max(p.line_width for p in chosen)
    error_x = x_axis.width_to_axis(stroke + tol) + x_axis.residual_max
    error_y = y_axis.width_to_axis(stroke + tol) + y_axis.residual_max
    return CurveResult(
        curve, xs, ys, error_x, error_y, "vector", tuple(p.index for p in chosen), stats, stroke
    )


def _raster_curve(recipe, curve, gray, page_height, x_axis, y_axis) -> CurveResult:
    spec = curve.raster
    scale = recipe.raster_dpi / 72.0
    ink = gray < spec.threshold

    def px(x, y):
        return round(x * scale - 0.5), round((page_height - y) * scale - 0.5)

    # Only the inside of the frame counts; masked boxes are cleared.
    fx0, fy0, fx1, fy1 = recipe.frame
    c0, r1 = px(fx0, fy0)
    c1, r0 = px(fx1, fy1)
    inside = np.zeros_like(ink)
    inside[r0 + 2 : r1 - 1, c0 + 2 : c1 - 1] = True
    boxes = []
    for x0, y0, x1, y1 in [*spec.mask, *recipe.exclude]:
        a0, b1 = px(x0, y0)
        a1, b0 = px(x1, y1)
        boxes.append((a0, b0, a1 + 1, b1 + 1))
    ink = mask_boxes(ink & inside, boxes)
    seed = px(*spec.seed)
    try:
        trace = trace_curve(
            ink,
            seed,
            max_jump=spec.max_jump_pt * scale,
            max_gap=math.ceil(spec.max_gap_pt * scale),
            max_height=math.ceil(spec.max_height_pt * scale),
        )
    except ValueError as exc:
        raise DigitizeError(f"curve {curve.output!r}: {exc}") from None
    xs_pt = (trace.cols + 0.5) / scale
    ys_pt = page_height - (trace.rows + 0.5) / scale
    data = to_data([np.column_stack([xs_pt, ys_pt])], x_axis, y_axis)
    xs, ys, stats = assemble(data, tuple(recipe.decimals), curve.monotonic)
    stroke = float(np.median(trace.heights)) / scale
    pixel = 1.0 / scale
    error_x = x_axis.width_to_axis(stroke + pixel) + x_axis.residual_max
    error_y = y_axis.width_to_axis(stroke + pixel) + y_axis.residual_max
    return CurveResult(curve, xs, ys, error_x, error_y, "raster", (), stats, stroke)


def run_recipe(recipe: Recipe, page: Page, gray: np.ndarray | None = None) -> FigureResult:
    """Apply ``recipe`` to ``page``; ``gray`` is the page rendering for raster recipes."""
    if page.number != recipe.page:
        raise DigitizeError(f"the recipe is for page {recipe.page}, got page {page.number}")
    x_axis = _calibrate(page, recipe, recipe.x)
    y_axis = _calibrate(page, recipe, recipe.y)
    region: BBox = tuple(recipe.region)
    if recipe.method == "vector" and page.images_in(region):
        raise DigitizeError(
            "the figure region contains an embedded image; the curves may be raster, "
            "use method: raster"
        )
    curves = []
    legend_paths: dict[str, int] = {}
    if recipe.method == "vector":
        # Legend samples are never curves; find them all first.
        for curve in recipe.curves:
            if curve.legend is not None:
                label = find_phrase(page, curve.legend, tuple(recipe.frame))
                legend_paths[curve.output] = legend_sample(page, label).index
        skip = set(legend_paths.values())
        for curve in recipe.curves:
            result = _vector_curve(page, recipe, curve, x_axis, y_axis, skip, legend_paths)
            _check_expect(curve, result.x, result.y)
            curves.append(result)
        used = [i for c in curves for i in c.paths]
        if len(used) != len(set(used)):
            raise DigitizeError("two curves selected the same path; make the styles distinct")
    else:
        if gray is None:
            raise DigitizeError("a raster recipe needs the rendered page")
        for curve in recipe.curves:
            result = _raster_curve(recipe, curve, gray, page.height, x_axis, y_axis)
            _check_expect(curve, result.x, result.y)
            curves.append(result)
    return FigureResult(recipe, x_axis, y_axis, tuple(curves), legend_paths)


# --------------------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------------------


def _fmt(value: float, decimals: int) -> str:
    text = f"{value:.{decimals}f}"
    if text.startswith("-") and float(text) == 0.0:
        text = text[1:]
    return text


def _error_text(axis: Axis, error: float, unit: str) -> str:
    text = axis.describe_error(error)
    if axis.scale == "log10" or unit == "1":
        return text
    return f"{text} {unit}"


def citation(recipe: Recipe, curve: CurveSpec) -> str:
    """The publication, edition, page and figure a curve comes from."""
    pub = recipe.publication
    which = f" (curve: {curve.legend})" if curve.legend else ""
    return (
        f"{pub.publisher}, {pub.title}, Publication {pub.number}, {pub.edition}, "
        f'p. {recipe.printed_page}, figure "{recipe.figure}"{which}'
    )


def csv_text(result: FigureResult, curve: CurveResult, recipe_relpath: str) -> str:
    """The reference-CSV text for one curve (Library format, LF line endings)."""
    recipe = result.recipe
    pub = recipe.publication
    if curve.method == "vector":
        how = "digitised from the PDF's vector paths"
    else:
        how = f"traced from a {recipe.raster_dpi:g} dpi raster rendering (embedded image)"
    ref = (
        f"{citation(recipe, curve.spec)}; {pub.url} (sha256 {pub.sha256}); {how} by "
        f"{recipe_relpath}; estimated error x +/-"
        f"{_error_text(result.x_axis, curve.error_x, recipe.x.unit)}, y +/-"
        f"{_error_text(result.y_axis, curve.error_y, recipe.y.unit)}"
    )
    conditions = "; ".join(f"{k}={v}" for k, v in curve.spec.conditions.items())
    header = {
        "material": recipe.material,
        "quantity": curve.spec.quantity,
        "x_unit": recipe.x.unit,
        "y_unit": recipe.y.unit,
        "origin": "datasheet",
        "ref": ref,
        "conditions": conditions,
    }
    lines = []
    for key, value in header.items():
        value = " ".join(str(value).split())
        lines.append(f"# {key}: {value}")
    dx, dy = recipe.decimals
    for x, y in zip(curve.x, curve.y):
        lines.append(f"{_fmt(float(x), dx)},{_fmt(float(y), dy)}")
    return "\n".join(lines) + "\n"
