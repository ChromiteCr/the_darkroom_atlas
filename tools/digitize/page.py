"""What the tool knows about one PDF page: vector paths, characters, words and images.

Plain dataclasses with no PDF library behind them, so tests can build pages by hand.
:mod:`tools.digitize.pdfread` fills them from a real PDF. Coordinates are PDF page
coordinates in points, origin at the bottom-left corner, y up.
"""

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

import numpy as np

from tools.digitize.geometry import Op, flatten_subpath, straight_segments

__all__ = [
    "BBox",
    "Char",
    "Page",
    "VectorPath",
    "Word",
    "bbox_contains",
    "bbox_overlaps",
    "dash_equal",
    "fmt_color",
    "group_words",
]

BBox = tuple[float, float, float, float]
Color = tuple[float, ...] | str | None


def bbox_contains(outer: BBox, inner: BBox, tolerance: float = 0.0) -> bool:
    """True if ``inner`` lies inside ``outer`` grown by ``tolerance`` on every side."""
    return (
        inner[0] >= outer[0] - tolerance
        and inner[1] >= outer[1] - tolerance
        and inner[2] <= outer[2] + tolerance
        and inner[3] <= outer[3] + tolerance
    )


def bbox_overlaps(a: BBox, b: BBox) -> bool:
    """True if the two boxes share any area or edge."""
    return a[0] <= b[2] and b[0] <= a[2] and a[1] <= b[3] and b[1] <= a[3]


def fmt_color(color: Color) -> str:
    if color is None:
        return "-"
    if isinstance(color, str):
        return color
    return "(" + ",".join(f"{c:g}" for c in color) + ")"


@dataclass(frozen=True)
class VectorPath:
    """One painted path, in drawing order.

    ``line_width`` and ``dash`` are in page points: the tool multiplies the values in the
    graphics state by the scale of the current transformation matrix. ``dash`` is the dash
    array (empty for a solid line) and ``dash_phase`` its phase.
    """

    index: int
    stroke: bool
    fill: bool
    line_width: float
    dash: tuple[float, ...]
    dash_phase: float
    stroke_color: Color
    fill_color: Color
    subpaths: tuple[tuple[Op, ...], ...]

    @property
    def points(self) -> np.ndarray:
        """Every anchor and control point, shape ``(N, 2)``."""
        pts = [p for sub in self.subpaths for op in sub for p in op.points]
        return np.asarray(pts, dtype=np.float64).reshape(-1, 2)

    @property
    def bbox(self) -> BBox:
        """Box around all anchor and control points (control points may lie outside)."""
        pts = self.points
        return (
            float(pts[:, 0].min()),
            float(pts[:, 1].min()),
            float(pts[:, 0].max()),
            float(pts[:, 1].max()),
        )

    @property
    def shape(self) -> str:
        """Operator letters, e.g. ``mlllh`` or ``mlccc``."""
        return "".join(op.kind for sub in self.subpaths for op in sub)

    @property
    def n_segments(self) -> int:
        """Number of drawing operators (everything except ``m``)."""
        return sum(op.kind != "m" for sub in self.subpaths for op in sub)

    def polylines(self, tolerance: float) -> list[np.ndarray]:
        """One flattened polyline per subpath."""
        return [flatten_subpath(sub, tolerance) for sub in self.subpaths]

    def segments(self) -> list[tuple[tuple[float, float], tuple[float, float]]]:
        """All straight segments of all subpaths."""
        return [seg for sub in self.subpaths for seg in straight_segments(sub)]

    @property
    def is_rect(self) -> bool:
        """True for one closed axis-aligned rectangle (``re`` or ``m l l l h``)."""
        if len(self.subpaths) != 1:
            return False
        pts = [op.points[0] for op in self.subpaths[0] if op.kind in ("m", "l")]
        if any(op.kind in ("c", "v", "y") for op in self.subpaths[0]):
            return False
        if len(pts) == 5 and pts[4] == pts[0]:
            pts = pts[:4]
        if len(pts) != 4:
            return False
        xs = sorted({round(p[0], 3) for p in pts})
        ys = sorted({round(p[1], 3) for p in pts})
        return len(xs) == 2 and len(ys) == 2

    def straight_line(self) -> tuple[tuple[float, float], tuple[float, float]] | None:
        """The segment if the path is exactly one straight ``m l`` line, else ``None``."""
        if len(self.subpaths) != 1:
            return None
        ops = self.subpaths[0]
        if len(ops) != 2 or ops[0].kind != "m" or ops[1].kind != "l":
            return None
        return ops[0].points[0], ops[1].points[0]

    def style(self) -> str:
        """Width, dash and colour, for listings."""
        dash = "solid" if not self.dash else "[" + ", ".join(f"{d:g}" for d in self.dash) + "]"
        return f"width={self.line_width:.4g} dash={dash} stroke={fmt_color(self.stroke_color)}"


@dataclass(frozen=True)
class Char:
    """One glyph and its box."""

    text: str
    x0: float
    y0: float
    x1: float
    y1: float
    size: float

    @property
    def bbox(self) -> BBox:
        return (self.x0, self.y0, self.x1, self.y1)


@dataclass(frozen=True)
class Word:
    """Consecutive glyphs on one line without a space or a wide gap between them."""

    chars: tuple[Char, ...]

    @property
    def text(self) -> str:
        return "".join(c.text for c in self.chars)

    @property
    def bbox(self) -> BBox:
        return (
            min(c.x0 for c in self.chars),
            min(c.y0 for c in self.chars),
            max(c.x1 for c in self.chars),
            max(c.y1 for c in self.chars),
        )

    @property
    def centre(self) -> tuple[float, float]:
        x0, y0, x1, y1 = self.bbox
        return ((x0 + x1) / 2, (y0 + y1) / 2)

    @property
    def size(self) -> float:
        return max(c.size for c in self.chars)


def group_words(chars: Iterable[Char]) -> list[Word]:
    """Split glyphs into words: same baseline, no whitespace and no gap wider than 0.25 em.

    Words come out sorted top to bottom, then left to right.
    """
    remaining = sorted(chars, key=lambda c: (-c.y0, c.x0))
    lines: list[list[Char]] = []
    for ch in remaining:
        for line in lines:
            ref = line[0]
            if abs(ch.y0 - ref.y0) <= 0.3 * max(ref.size, ch.size, 1e-9):
                line.append(ch)
                break
        else:
            lines.append([ch])
    words: list[Word] = []
    for line in sorted(lines, key=lambda ln: -ln[0].y0):
        line.sort(key=lambda c: c.x0)
        current: list[Char] = []
        for ch in line:
            if ch.text.isspace():
                if current:
                    words.append(Word(tuple(current)))
                current = []
                continue
            if current and ch.x0 - current[-1].x1 > 0.25 * max(ch.size, 1e-9):
                words.append(Word(tuple(current)))
                current = []
            current.append(ch)
        if current:
            words.append(Word(tuple(current)))
    return words


@dataclass(frozen=True)
class Page:
    """One page: size, paths in drawing order, glyphs and embedded-image boxes."""

    number: int
    width: float
    height: float
    paths: tuple[VectorPath, ...]
    chars: tuple[Char, ...] = ()
    images: tuple[BBox, ...] = ()
    words: tuple[Word, ...] = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "words", tuple(group_words(self.chars)))

    def paths_in(self, region: BBox, tolerance: float = 0.0) -> list[VectorPath]:
        """Paths lying entirely inside ``region``."""
        return [p for p in self.paths if bbox_contains(region, p.bbox, tolerance)]

    def words_in(self, region: BBox) -> list[Word]:
        """Words whose centre lies inside ``region``."""
        out = []
        for w in self.words:
            cx, cy = w.centre
            if region[0] <= cx <= region[2] and region[1] <= cy <= region[3]:
                out.append(w)
        return out

    def images_in(self, region: BBox) -> list[BBox]:
        """Embedded images overlapping ``region``."""
        return [b for b in self.images if bbox_overlaps(region, b)]


def dash_equal(a: Sequence[float], b: Sequence[float], tolerance: float) -> bool:
    """True if two dash arrays have the same length and agree element by element."""
    return len(a) == len(b) and all(math.isclose(x, y, abs_tol=tolerance) for x, y in zip(a, b))
