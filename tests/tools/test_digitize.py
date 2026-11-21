"""Tests for the datasheet digitising tool (tools/digitize).

Nothing here needs a downloaded PDF. Most tests build a page (paths, glyphs) by hand; the
last group writes a tiny PDF in memory and runs only when pdfminer.six is installed
(``uv run --group digitize pytest tests/tools``).
"""

import math

import numpy as np
import pytest
from PIL import Image, ImageDraw

from darkroom_atlas.data.library import Library
from tools.digitize.calibrate import Axis, Reference, axis_from_pairs
from tools.digitize.figure import DigitizeError, assemble, find_phrase, has_bar
from tools.digitize.geometry import Op, flatten_cubic, flatten_subpath, polyline_deviation
from tools.digitize.page import Char, Page, VectorPath, group_words
from tools.digitize.raster import trace_curve
from tools.digitize.recipe import Recipe, csv_text, run_recipe
from tools.digitize.runner import url_allowed

# --------------------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------------------


def _bezier(p0, p1, p2, p3, t):
    t = np.asarray(t)[:, None]
    p0, p1, p2, p3 = (np.asarray(p, float) for p in (p0, p1, p2, p3))
    return (1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1 + 3 * (1 - t) * t**2 * p2 + t**3 * p3


@pytest.mark.parametrize("tolerance", [0.5, 0.05, 0.005])
def test_flatten_cubic_stays_within_tolerance(tolerance):
    ctrl = [(0, 0), (30, 90), (80, -40), (100, 20)]
    points = np.array([ctrl[0], *flatten_cubic(*ctrl, tolerance)])
    assert tuple(points[-1]) == (100.0, 20.0)
    exact = _bezier(*ctrl, np.linspace(0, 1, 2001))
    # Every point of the true curve is within the tolerance of the polyline.
    assert polyline_deviation(points, exact) <= tolerance
    # Tighter tolerance, more points.
    assert len(points) > 3


def test_flatten_cubic_is_deterministic_and_needs_positive_tolerance():
    ctrl = [(0, 0), (1, 2), (3, 2), (4, 0)]
    assert flatten_cubic(*ctrl, 0.01) == flatten_cubic(*ctrl, 0.01)
    with pytest.raises(ValueError):
        flatten_cubic(*ctrl, 0.0)


def test_flatten_subpath_handles_v_y_and_close():
    ops = [
        Op("m", ((0.0, 0.0),)),
        Op("l", ((10.0, 0.0),)),
        Op("v", ((20.0, 10.0), (20.0, 20.0))),
        Op("y", ((10.0, 30.0), (0.0, 20.0))),
        Op("h", ()),
    ]
    poly = flatten_subpath(ops, 0.01)
    assert tuple(poly[0]) == (0.0, 0.0)
    assert tuple(poly[1]) == (10.0, 0.0)
    assert (20.0, 20.0) in [tuple(p) for p in poly]
    assert tuple(poly[-1]) == (0.0, 0.0)  # closed back to the start
    v_exact = _bezier((10, 0), (10, 0), (20, 10), (20, 20), np.linspace(0, 1, 500))
    assert polyline_deviation(poly, v_exact) <= 0.01


def test_op_checks_its_point_count():
    with pytest.raises(ValueError):
        Op("c", ((0.0, 0.0),))
    with pytest.raises(ValueError):
        Op("q", ())


# --------------------------------------------------------------------------------------
# Calibration
# --------------------------------------------------------------------------------------


def test_linear_axis_maps_and_inverts():
    axis = axis_from_pairs([(100.0, -4.0), (136.886, -3.0), (284.43, 1.0)])
    assert axis.to_data(100.0) == pytest.approx(-4.0)
    assert axis.to_data(192.215) == pytest.approx(-1.5, abs=1e-3)
    values = np.array([-3.7, 0.0, 0.35])
    assert axis.to_data(axis.to_page(values)) == pytest.approx(values)
    assert axis.residual_max < 1e-3
    assert axis.max_slope == pytest.approx(5.0 / 184.43, rel=1e-3)


def test_log10_axis():
    axis = axis_from_pairs([(0.0, 1.0), (50.0, 10.0), (100.0, 100.0)], scale="log10")
    assert axis.to_data(25.0) == pytest.approx(math.sqrt(10.0))
    assert axis.to_page(1000.0) == pytest.approx(150.0)
    assert axis.residual_max == pytest.approx(0.0, abs=1e-12)
    assert "decades" in axis.describe_error(0.01)


def test_lsq_residual_reports_a_misplaced_tick():
    # The 5 tick is drawn 1 pt too far right: with 10 pt per unit that is 0.1 units.
    pairs = [(0.0, 0.0), (51.0, 5.0), (100.0, 10.0)]
    axis = axis_from_pairs(pairs)
    coords, values = np.array(pairs).T
    expected = values - np.polyval(np.polyfit(coords, values, 1), coords)
    assert axis.residuals == pytest.approx(expected)
    # Least squares spreads the 0.1-unit error: the largest residual is about 2/3 of it.
    assert axis.residual_max == pytest.approx(0.0667, abs=1e-4)
    assert axis.residual_rms == pytest.approx(np.sqrt(np.mean(expected**2)))
    # A straight line through two of the three ticks would hide the error entirely.
    assert axis_from_pairs(pairs[:2]).residual_max == pytest.approx(0.0, abs=1e-12)


def test_piecewise_axis_passes_through_every_reference():
    pairs = [(0.0, 1.0), (71.7, 10.0), (144.9, 100.0), (166.3, 200.0)]
    axis = axis_from_pairs(pairs, scale="log10", fit="piecewise")
    for coord, value in pairs:
        assert axis.to_data(coord) == pytest.approx(value)
        assert axis.to_page(value) == pytest.approx(coord)
    # Beyond the last reference the last segment is extended.
    slope = (math.log10(200) - 2.0) / (166.3 - 144.9)
    assert axis.to_data(176.3) == pytest.approx(10 ** (math.log10(200) + 10 * slope))
    # The least-squares residual still measures the uneven grid.
    assert axis.residual_max > 0.001


@pytest.mark.parametrize(
    ("pairs", "scale"),
    [
        ([(0.0, 0.0)], "linear"),  # one point
        ([(0.0, 0.0), (0.0, 1.0)], "linear"),  # same coordinate
        ([(0.0, 1.0), (1.0, 1.0)], "linear"),  # same value
        ([(0.0, 0.0), (1.0, 10.0)], "log10"),  # non-positive value on a log axis
        ([(0.0, 0.0), (1.0, 2.0), (2.0, 1.0)], "linear"),  # not monotonic
    ],
)
def test_axis_rejects_unusable_references(pairs, scale):
    with pytest.raises(ValueError):
        axis_from_pairs(pairs, scale=scale)


def test_reversed_axis_works():
    # A y axis drawn downwards (value grows as the coordinate falls).
    axis = Axis((Reference(300.0, 0.0), Reference(100.0, 4.0)))
    assert axis.to_data(200.0) == pytest.approx(2.0)
    assert axis.to_page(1.0) == pytest.approx(250.0)


# --------------------------------------------------------------------------------------
# Assembling points
# --------------------------------------------------------------------------------------


def test_assemble_orders_pieces_and_drops_duplicates():
    a = np.array([[2.0, 2.0], [1.0, 1.0]])  # drawn right to left
    b = np.array([[2.0, 2.0], [3.0, 3.0]])  # starts where a ends
    xs, ys, stats = assemble([b, a], (4, 4), "strict")
    assert xs.tolist() == [1.0, 2.0, 3.0]
    assert ys.tolist() == [1.0, 2.0, 3.0]
    assert stats.reversed_pieces == 1
    assert stats.duplicates == 1


def test_assemble_backtracking_is_an_error_unless_allowed():
    poly = np.array([[0.0, 0.0], [1.0, 1.0], [0.9, 0.5], [2.0, 2.0]])
    with pytest.raises(DigitizeError, match="backwards"):
        assemble([poly], (4, 4), "strict")
    xs, _, stats = assemble([poly], (4, 4), "drop-backtracking")
    assert xs.tolist() == [0.0, 1.0, 2.0]
    assert stats.backtracking == 1


def test_assemble_drops_points_that_round_to_the_same_x():
    poly = np.array([[0.0, 0.0], [1.00001, 1.0], [1.00002, 1.5], [2.0, 2.0]])
    xs, ys, stats = assemble([poly], (4, 4), "strict")
    assert xs.tolist() == [0.0, 1.0, 2.0]
    assert ys.tolist() == [0.0, 1.0, 2.0]
    assert stats.same_x == 1


# --------------------------------------------------------------------------------------
# A synthetic chart
# --------------------------------------------------------------------------------------
#
# Frame (100, 100)-(300, 300) pt. x: log exposure -2..+2, 50 pt per unit, labels in bar
# notation for negative values. y: density 0..4, 50 pt per unit. Two curves, told apart by
# dash pattern, each with a legend sample.

FRAME = (100.0, 100.0, 300.0, 300.0)


def _px(v):
    return 200.0 + 50.0 * v


def _py(d):
    return 100.0 + 50.0 * d


class _Paths:
    def __init__(self):
        self.items: list[VectorPath] = []

    def add(self, subpaths, width=0.432, dash=(), stroke=True, fill=False):
        self.items.append(
            VectorPath(
                index=len(self.items), stroke=stroke, fill=fill, line_width=width,
                dash=tuple(dash), dash_phase=0.0, stroke_color=(0.0,), fill_color=(0.0,),
                subpaths=tuple(tuple(sub) for sub in subpaths),
            )
        )  # fmt: skip

    def line(self, p, q, **kw):
        self.add([[Op("m", (p,)), Op("l", (q,))]], **kw)

    def rect(self, x0, y0, x1, y1, **kw):
        self.add(
            [[Op("m", ((x0, y0),)), Op("l", ((x1, y0),)), Op("l", ((x1, y1),)),
              Op("l", ((x0, y1),)), Op("h", ())]],
            **kw,
        )  # fmt: skip


def _text(chars, text, x0, y_centre, size=7.0):
    width = 0.55 * size
    for i, ch in enumerate(text):
        left = x0 + i * width
        chars.append(Char(ch, left, y_centre - size / 2, left + width, y_centre + size / 2, size))


def _text_centred(chars, text, cx, cy, size=7.0):
    _text(chars, text, cx - 0.55 * size * len(text) / 2, cy, size)


SOLID_DATA = [(-1.5, 0.5), (-0.5, 1.0), (0.5, 2.0), (1.5, 3.5)]
BEZIER = [(-1.5, 0.3), (-0.5, 0.3), (0.5, 1.8), (1.5, 2.6)]  # data-space control points


def _chart(bar_on_minus_one=True) -> Page:
    paths = _Paths()
    chars: list[Char] = []
    paths.rect(*FRAME)
    # x ticks at -1, 0, +1 (frame edges serve -2 and +2) with labels below.
    for v in (-1, 0, 1):
        paths.line((_px(v), 100.0), (_px(v), 104.0))
    for v, text in ((-2, "2.0"), (-1, "1.0"), (0, "0.0"), (1, "1.0"), (2, "2.0")):
        _text_centred(chars, text, _px(v), 90.0)
        if v < 0 and (v != -1 or bar_on_minus_one):
            glyph_x0 = _px(v) - 0.55 * 7 * 3 / 2
            paths.line((glyph_x0 + 0.3, 94.3), (glyph_x0 + 3.5, 94.3))  # the bar
    # y ticks and labels.
    for d in (1, 2, 3):
        paths.line((100.0, _py(d)), (104.0, _py(d)))
    for d in range(5):
        _text_centred(chars, f"{d}.0", 88.0, _py(d))
    # Legend: samples left of their labels, inside the frame.
    paths.line((200.0, 132.0), (225.0, 132.0), width=0.72)
    _text(chars, "solid line", 228.0, 132.0)
    paths.line((200.0, 122.0), (225.0, 122.0), width=0.72, dash=(3.0, 1.5))
    _text(chars, "dashed line", 228.0, 122.0)
    # Curves.
    pts = [(_px(x), _py(y)) for x, y in SOLID_DATA]
    paths.add([[Op("m", (pts[0],)), *[Op("l", (p,)) for p in pts[1:]]]], width=0.72)
    ctrl = [(_px(x), _py(y)) for x, y in BEZIER]
    paths.add([[Op("m", (ctrl[0],)), Op("c", tuple(ctrl[1:]))]], width=0.72, dash=(3.0, 1.5))
    return Page(number=1, width=400.0, height=400.0, paths=tuple(paths.items), chars=tuple(chars))


def _recipe(**overrides) -> Recipe:
    data = {
        "schema_version": 1,
        "publication": {
            "publisher": "Test",
            "number": "T-1",
            "title": "Synthetic",
            "edition": "test",
            "url": "https://www.kodak.com/test.pdf",
            "file": "test.pdf",
            "sha256": "0" * 64,
        },
        "page": 1,
        "printed_page": "1",
        "figure": "Synthetic chart",
        "caption": "none",
        "material": "test/film",
        "region": [50.0, 50.0, 350.0, 350.0],
        "frame": list(FRAME),
        "x": {
            "title": "LOG EXPOSURE",
            "unit": "log10(lux*s)",
            "side": "bottom",
            "bar_notation": True,
            "ticks": [
                {"label": "2.0", "bar": True, "value": -2.0},
                {"label": "1.0", "bar": True, "value": -1.0},
                {"label": "0.0", "value": 0.0},
                {"label": "1.0", "value": 1.0},
                {"label": "2.0", "value": 2.0},
            ],
        },
        "y": {
            "title": "DENSITY",
            "unit": "1",
            "side": "left",
            "ticks": [{"label": f"{d}.0", "value": float(d)} for d in range(5)],
        },
        "curves": [
            {
                "output": "solid",
                "quantity": "characteristic_curve",
                "legend": "solid line",
                "style": {"width": 0.72, "dash": []},
                "expect": {"paths": 1, "x_span": [-1.5, 1.5], "tolerance": [1e-6, 1e-6]},
                "conditions": {"developer": "test", "time": "6min"},
            },
            {
                "output": "dashed",
                "quantity": "characteristic_curve",
                "legend": "dashed line",
                "style": {"width": 0.72, "dash": [3.0, 1.5]},
                "expect": {"paths": 1, "x_span": [-1.5, 1.5], "tolerance": [1e-6, 1e-6]},
                "conditions": {"developer": "test", "time": "9min"},
            },
        ],
    }
    data.update(overrides)
    return Recipe.model_validate(data)


def test_group_words_splits_on_spaces_and_gaps():
    chars: list[Char] = []
    _text(chars, "11 minutes", 10.0, 50.0)
    _text(chars, "3", 100.0, 50.0)
    words = group_words(chars)
    assert [w.text for w in words] == ["11", "minutes", "3"]


def test_bar_notation_is_read_from_paths():
    page = _chart()
    barred = {w.centre[0]: has_bar(w, page) for w in page.words if w.centre[1] < 95}
    assert barred[_px(-2)] and barred[_px(-1)]
    assert not barred[_px(0)] and not barred[_px(1)] and not barred[_px(2)]


def test_legend_phrase_must_match_a_whole_text_run():
    page = _chart()
    find_phrase(page, "solid line", FRAME)
    with pytest.raises(DigitizeError):
        find_phrase(page, "line", FRAME)  # appears in two runs


def test_recipe_digitises_the_synthetic_chart_exactly():
    page = _chart()
    result = run_recipe(_recipe(), page)
    assert [r.coord for r in result.x_axis.references] == pytest.approx(
        [100.0, 150.0, 200.0, 250.0, 300.0]
    )
    solid, dashed = result.curves
    assert solid.paths != dashed.paths
    np.testing.assert_allclose(np.column_stack([solid.x, solid.y]), SOLID_DATA, atol=1e-4)
    # The dashed Bezier, flattened, lies on the true Bezier (in data space) within the
    # flattening tolerance converted to data units, plus rounding.
    exact = _bezier(*BEZIER, np.linspace(0, 1, 2001))
    poly = np.column_stack([dashed.x, dashed.y])
    assert polyline_deviation(poly, exact) <= 0.01 / 50.0 + 1e-4
    assert dashed.x[0] == -1.5 and dashed.x[-1] == 1.5
    # Error estimate: line width plus flattening tolerance, in data units (residual is 0).
    assert solid.error_x == pytest.approx((0.72 + 0.01) / 50.0)
    assert solid.error_y == pytest.approx((0.72 + 0.01) / 50.0)


def test_missing_bar_is_caught():
    page = _chart(bar_on_minus_one=False)
    with pytest.raises(DigitizeError, match="with a bar"):
        run_recipe(_recipe(), page)


def test_wrong_style_is_caught_by_the_legend_check():
    recipe = _recipe()
    curves = [c.model_dump() for c in recipe.curves]
    curves[0]["style"] = {"width": 0.72, "dash": [3.0, 1.5]}
    with pytest.raises(DigitizeError, match="legend sample"):
        run_recipe(_recipe(curves=curves), _chart())


def test_unexpected_span_is_caught():
    recipe = _recipe()
    curves = [c.model_dump() for c in recipe.curves]
    curves[0]["expect"]["x_span"] = [-1.0, 1.5]
    with pytest.raises(DigitizeError, match="expected"):
        run_recipe(_recipe(curves=curves), _chart())


def test_csv_is_read_back_by_the_library(tmp_path):
    result = run_recipe(_recipe(), _chart())
    curve = result.curves[0]
    text = csv_text(result, curve, "tools/digitize/recipes/test/chart.yaml")
    assert text.endswith("\n") and "\r" not in text
    assert "digitised from the PDF's vector paths" in text
    target = tmp_path / "reference" / "test" / "film" / "solid.csv"
    target.parent.mkdir(parents=True)
    target.write_text(text, encoding="utf-8")
    loaded = Library(tmp_path).reference_curve("test/film/solid.csv")
    assert loaded.material == "test/film"
    assert loaded.x_unit == "log10(lux*s)"
    assert loaded.conditions == "developer=test; time=6min"
    np.testing.assert_allclose(np.column_stack([loaded.x, loaded.y]), SOLID_DATA, atol=1e-4)
    # Same input, same bytes.
    again = run_recipe(_recipe(), _chart())
    assert csv_text(again, again.curves[0], "tools/digitize/recipes/test/chart.yaml") == text


# --------------------------------------------------------------------------------------
# Raster fallback
# --------------------------------------------------------------------------------------


def test_trace_curve_follows_a_dashed_line():
    image = Image.new("L", (400, 200), 255)
    draw = ImageDraw.Draw(image)
    xs = np.arange(20, 380)
    ys = 100 + 40 * np.sin(xs / 50.0)
    for x0 in range(20, 380, 12):  # 8 px dashes, 4 px gaps
        seg = (xs >= x0) & (xs < x0 + 8)
        draw.line(list(zip(xs[seg], ys[seg])), fill=0, width=3)
    ink = np.asarray(image) < 128
    trace = trace_curve(ink, (100, round(100 + 40 * math.sin(2.0))), max_gap=6)
    # The first dash starts at column 20, the last ends at column 375.
    assert trace.cols[0] <= 20 and trace.cols[-1] >= 375
    truth = 100 + 40 * np.sin(trace.cols / 50.0)
    error = np.abs(trace.rows - truth)
    assert np.median(error) <= 0.75
    # Dash ends give slivers; the worst point stays inside the error the recipe states for
    # raster curves (median run height plus one pixel).
    assert np.max(error) <= np.median(trace.heights) + 1.0


def test_raster_recipe_states_a_larger_error():
    page = _chart()
    scale = 600 / 72
    image = Image.new("L", (int(400 * scale), int(400 * scale)), 255)
    draw = ImageDraw.Draw(image)
    pts = [(_px(x) * scale, (400 - _py(y)) * scale) for x, y in SOLID_DATA]
    draw.line(pts, fill=0, width=int(0.72 * scale))
    gray = np.asarray(image)
    recipe = _recipe(
        method="raster",
        curves=[
            {
                "output": "solid", "quantity": "characteristic_curve",
                "raster": {"seed": [_px(0.0), _py(1.5)], "mask": [[195.0, 115.0, 290.0, 140.0]]},
                "conditions": {"developer": "test"},
            }
        ],
    )  # fmt: skip
    result = run_recipe(recipe, page, gray)
    curve = result.curves[0]
    assert curve.method == "raster"
    truth = np.interp(curve.x, *zip(*SOLID_DATA))
    assert np.max(np.abs(curve.y - truth)) <= curve.error_y
    assert curve.error_y > (0.72 + 0.01) / 50.0  # larger than the vector estimate
    assert "raster rendering" in csv_text(result, curve, "r.yaml")


# --------------------------------------------------------------------------------------
# Fetch policy
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("url", "ok"),
    [
        ("https://kodakprofessional.com/sites/default/files/x.pdf", True),
        ("https://imaging.kodakalaris.com/x.pdf", True),
        ("https://web.archive.org/web/2009id_/http://www.kodak.com:80/x.pdf", True),
        ("https://125px.com/docs/film/kodak/f4017-2016.pdf", False),
        ("https://web.archive.org/web/2009id_/http://125px.com/x.pdf", False),
        ("https://kodak.com.evil.example/x.pdf", False),
        ("ftp://www.kodak.com/x.pdf", False),
    ],
)
def test_url_allowed(url, ok):
    assert url_allowed(url) is ok


# --------------------------------------------------------------------------------------
# A tiny real PDF (needs pdfminer.six from the digitize group)
# --------------------------------------------------------------------------------------


def _pdf(content: bytes) -> bytes:
    objects = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        (
            b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 400 400] /Contents 4 0 R "
            b"/Resources << /Font << /F1 5 0 R >> >> >>"
        ),
        b"<< /Length %d >>\nstream\n" % len(content) + content + b"\nendstream",
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
    ]
    out = bytearray(b"%PDF-1.4\n")
    offsets = []
    for i, body in enumerate(objects, start=1):
        offsets.append(len(out))
        out += b"%d 0 obj\n" % i + body + b"\nendobj\n"
    xref = len(out)
    out += b"xref\n0 %d\n0000000000 65535 f \n" % (len(objects) + 1)
    for offset in offsets:
        out += b"%010d 00000 n \n" % offset
    out += b"trailer\n<< /Size %d /Root 1 0 R >>\nstartxref\n%d\n%%%%EOF\n" % (
        len(objects) + 1, xref,
    )  # fmt: skip
    return bytes(out)


def _label(text, cx, cy, size=7):
    # Helvetica digits and '.' are 0.556 and 0.278 em wide.
    width = sum(0.278 if ch == "." else 0.556 for ch in text) * size
    return b"BT /F1 %d Tf %.3f %.3f Td (%s) Tj ET\n" % (
        size, cx - width / 2, cy - 0.35 * size, text.encode(),
    )  # fmt: skip


def test_reads_paths_and_text_from_a_pdf(tmp_path):
    pytest.importorskip("pdfminer")
    from tools.digitize.pdfread import load_page

    ctrl = [(_px(x), _py(y)) for x, y in BEZIER]
    content = bytearray()
    content += b"q 1 0 0 1 0 0 cm 0.432 w 100 100 200 200 re S Q\n"
    for v in (-1, 0, 1):
        content += b"q 0.432 w %g 100 m %g 104 l S Q\n" % (_px(v), _px(v))
    for d in (1, 2, 3):
        content += b"q 0.432 w 100 %g m 104 %g l S Q\n" % (_py(d), _py(d))
    for v in (-2, -1, 0, 1, 2):
        content += _label(f"{abs(v)}.0", _px(v), 90.0)
        if v < 0:
            # A bar over the first digit: "2.0" is 1.39 em wide, the digit 0.556 em.
            x0 = _px(v) - 1.39 * 7 / 2 + 0.4
            content += b"q 0.432 w %.3f 94.6 m %.3f 94.6 l S Q\n" % (x0, x0 + 3.0)
    for d in range(5):
        content += _label(f"{d}.0", 88.0, _py(d))
    # The curve is drawn at half scale so the tool must apply the matrix to width and dash.
    content += b"q 0.5 0 0 0.5 0 0 cm 1.44 w [6 3] 0 d %g %g m %g %g %g %g %g %g c S Q\n" % (
        2 * ctrl[0][0], 2 * ctrl[0][1], 2 * ctrl[1][0], 2 * ctrl[1][1],
        2 * ctrl[2][0], 2 * ctrl[2][1], 2 * ctrl[3][0], 2 * ctrl[3][1],
    )  # fmt: skip
    pdf = tmp_path / "synthetic.pdf"
    pdf.write_bytes(_pdf(bytes(content)))

    page = load_page(pdf, 1)
    assert page.width == 400 and page.height == 400
    curve = [p for p in page.paths if p.shape == "mc"]
    assert len(curve) == 1
    assert curve[0].line_width == pytest.approx(0.72)
    assert curve[0].dash == pytest.approx((3.0, 1.5))
    assert any(p.is_rect for p in page.paths)
    labels = sorted(w.text for w in page.words)
    assert labels.count("1.0") == 3 and labels.count("2.0") == 3  # two on x, one on y

    recipe = _recipe(
        curves=[
            {
                "output": "dashed", "quantity": "characteristic_curve",
                "style": {"width": 0.72, "dash": [3.0, 1.5]},
                "expect": {"paths": 1, "x_span": [-1.5, 1.5], "tolerance": [1e-3, 1e-3]},
                "conditions": {"developer": "test"},
            }
        ]
    )  # fmt: skip
    result = run_recipe(recipe, page)
    dashed = result.curves[0]
    exact = _bezier(*BEZIER, np.linspace(0, 1, 2001))
    assert polyline_deviation(np.column_stack([dashed.x, dashed.y]), exact) <= 1e-3
