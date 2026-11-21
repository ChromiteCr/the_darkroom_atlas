"""Read one PDF page into a :class:`~tools.digitize.page.Page` with pdfminer.six.

pdfminer's own layout objects drop the transformation matrix, so line widths and dash
arrays would come out in whatever local units the drawing used. This module records every
painted path itself, together with the matrix in force, and scales width and dash to page
points. Text comes from pdfminer's glyph boxes.

Only parsing happens here: a PDF is untrusted input. Run the tool with ``python -I``.
"""

import math
from pathlib import Path

from tools.digitize.geometry import Op
from tools.digitize.page import Char, Page, VectorPath

__all__ = ["load_page", "page_count"]


def _imports():
    try:
        from pdfminer.converter import PDFPageAggregator
        from pdfminer.layout import LTChar, LTContainer, LTImage
        from pdfminer.pdfdocument import PDFDocument
        from pdfminer.pdfinterp import PDFPageInterpreter, PDFResourceManager
        from pdfminer.pdfpage import PDFPage
        from pdfminer.pdfparser import PDFParser
        from pdfminer.utils import apply_matrix_pt
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise SystemExit(
            "pdfminer.six is missing: run the tool with `uv run --group digitize ...`"
        ) from exc
    return (
        PDFPageAggregator, LTChar, LTContainer, LTImage, PDFDocument, PDFPageInterpreter,
        PDFResourceManager, PDFPage, PDFParser, apply_matrix_pt,
    )  # fmt: skip


def _number(value) -> float:
    return float(value)


def _color(value):
    if value is None:
        return None
    if isinstance(value, int | float):
        return (round(float(value), 4),)
    if isinstance(value, list | tuple):
        try:
            return tuple(round(float(v), 4) for v in value)
        except (TypeError, ValueError):
            return str(value)
    return str(value)


def _ctm_scale(ctm) -> float:
    a, b, c, d, _, _ = ctm
    return math.sqrt(abs(a * d - b * c))


def page_count(pdf_path: Path | str) -> int:
    """Number of pages in the PDF."""
    *_, PDFDocument, _, _, PDFPage, PDFParser, _ = _imports()
    with open(pdf_path, "rb") as handle:
        document = PDFDocument(PDFParser(handle))
        return sum(1 for _ in PDFPage.create_pages(document))


def load_page(pdf_path: Path | str, number: int) -> Page:
    """Paths, glyphs and image boxes of page ``number`` (1-based) of ``pdf_path``."""
    (
        PDFPageAggregator, LTChar, LTContainer, LTImage, PDFDocument, PDFPageInterpreter,
        PDFResourceManager, PDFPage, PDFParser, apply_matrix_pt,
    ) = _imports()  # fmt: skip

    class _Recorder(PDFPageAggregator):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.recorded: list[VectorPath] = []
            self._inside = False

        def paint_path(self, gstate, stroke, fill, evenodd, path):
            # pdfminer splits multi-subpath paths by calling paint_path again for each
            # subpath; record the path once, as painted.
            if self._inside:
                super().paint_path(gstate, stroke, fill, evenodd, path)
                return
            self._record(gstate, stroke, fill, path)
            self._inside = True
            try:
                super().paint_path(gstate, stroke, fill, evenodd, path)
            finally:
                self._inside = False

        def _record(self, gstate, stroke, fill, path):
            ctm = self.ctm
            scale = _ctm_scale(ctm)
            subpaths: list[list[Op]] = []
            for operation in path:
                kind = str(operation[0])
                nums = [_number(v) for v in operation[1:]]
                pts = tuple(
                    tuple(float(v) for v in apply_matrix_pt(ctm, (nums[i], nums[i + 1])))
                    for i in range(0, len(nums) - 1, 2)
                )
                if kind == "m" or not subpaths:
                    if kind != "m":
                        break  # malformed path without a starting point
                    subpaths.append([])
                subpaths[-1].append(Op(kind, pts))
            dash_array, dash_phase = gstate.dash or ((), 0)
            self.recorded.append(
                VectorPath(
                    index=len(self.recorded),
                    stroke=bool(stroke),
                    fill=bool(fill),
                    line_width=float(gstate.linewidth) * scale,
                    dash=tuple(_number(v) * scale for v in (dash_array or ())),
                    dash_phase=_number(dash_phase or 0) * scale,
                    stroke_color=_color(gstate.scolor),
                    fill_color=_color(gstate.ncolor),
                    # A lone 'm' paints nothing; drop such subpaths.
                    subpaths=tuple(tuple(sub) for sub in subpaths if len(sub) > 1),
                )
            )

    class _Interpreter(PDFPageInterpreter):
        # pdfminer multiplies the line width by the matrix in force when `w` runs; the PDF
        # rule is the matrix in force when the path is stroked. Keep the raw width here and
        # scale it (and the dash array) in _record.
        def do_w(self, linewidth):
            try:
                self.graphicstate.linewidth = float(linewidth)
            except (TypeError, ValueError):
                pass

    with open(pdf_path, "rb") as handle:
        document = PDFDocument(PDFParser(handle))
        manager = PDFResourceManager()
        device = _Recorder(manager, laparams=None)
        interpreter = _Interpreter(manager, device)
        for i, pdf_page in enumerate(PDFPage.create_pages(document), start=1):
            if i == number:
                interpreter.process_page(pdf_page)
                layout = device.get_result()
                break
        else:
            raise ValueError(f"{pdf_path}: has no page {number}")

    chars: list[Char] = []
    images: list[tuple[float, float, float, float]] = []

    def walk(item):
        if isinstance(item, LTChar):
            chars.append(Char(item.get_text(), item.x0, item.y0, item.x1, item.y1, item.size))
        elif isinstance(item, LTImage):
            images.append((item.x0, item.y0, item.x1, item.y1))
        if isinstance(item, LTContainer):
            for child in item:
                walk(child)

    walk(layout)
    paths = tuple(p for p in device.recorded if p.subpaths)
    return Page(
        number=number,
        width=float(layout.width),
        height=float(layout.height),
        paths=paths,
        chars=tuple(chars),
        images=tuple(images),
    )
