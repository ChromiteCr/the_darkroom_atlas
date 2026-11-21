"""Checking images: the rendered figure region with the digitised points drawn on top."""

from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

from tools.digitize.recipe import FigureResult

__all__ = ["draw_overlay", "render_rgb"]

# Saturated colours that stand out against black-on-white charts; cycled per curve.
PALETTE = [
    (230, 0, 120),
    (0, 160, 255),
    (255, 120, 0),
    (0, 190, 90),
    (150, 60, 255),
    (220, 30, 30),
    (0, 200, 200),
    (200, 170, 0),
]
TICK_COLOUR = (0, 120, 255)


def render_rgb(pdf_path: Path | str, page_number: int, scale: float) -> Image.Image:
    """RGB rendering of a page at ``scale`` pixels per point."""
    try:
        import pypdfium2 as pdfium
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise SystemExit("pypdfium2 is missing: run with `uv run --group digitize ...`") from exc
    document = pdfium.PdfDocument(str(pdf_path))
    try:
        return document[page_number - 1].render(scale=scale).to_pil().convert("RGB")
    finally:
        document.close()


def draw_overlay(
    page_image: Image.Image, page_height: float, scale: float, result: FigureResult
) -> Image.Image:
    """Crop the recipe region from ``page_image`` and mark points, ticks and the frame."""
    recipe = result.recipe
    x0, y0, x1, y1 = recipe.region

    def to_px(x, y):
        return (x - x0) * scale, (y1 - y) * scale

    left, top = int(np.floor(x0 * scale)), int(np.floor((page_height - y1) * scale))
    right, bottom = int(np.ceil(x1 * scale)), int(np.ceil((page_height - y0) * scale))
    image = page_image.crop((left, top, right, bottom)).convert("RGB")
    # Fade the page so the marks stand out.
    image = Image.blend(image, Image.new("RGB", image.size, (255, 255, 255)), 0.35)
    draw = ImageDraw.Draw(image)
    fx0, fy0, fx1, fy1 = recipe.frame
    draw.rectangle([to_px(fx0, fy1), to_px(fx1, fy0)], outline=(0, 200, 0), width=1)
    for ref in result.x_axis.references:
        px, py = to_px(ref.coord, fy0)
        draw.line([(px, py - 8), (px, py + 8)], fill=TICK_COLOUR, width=2)
    for ref in result.y_axis.references:
        px, py = to_px(fx0, ref.coord)
        draw.line([(px - 8, py), (px + 8, py)], fill=TICK_COLOUR, width=2)
    radius = max(1.5, scale * 0.45)
    for i, curve in enumerate(result.curves):
        colour = PALETTE[i % len(PALETTE)]
        xs = result.x_axis.to_page(curve.x)
        ys = result.y_axis.to_page(curve.y)
        pts = [to_px(float(a), float(b)) for a, b in zip(xs, ys)]
        if len(pts) > 1:
            draw.line(pts, fill=colour, width=1)
        for px, py in pts:
            draw.ellipse([px - radius, py - radius, px + radius, py + radius], outline=colour)
        lx, ly = pts[-1]
        draw.text((lx + 4, ly - 6), curve.spec.output, fill=colour)
    return image
