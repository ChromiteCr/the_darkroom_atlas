"""Command line: ``uv run --group digitize python -I -m tools.digitize <command> ...``.

Commands:

* ``paths PDF PAGE [--region x0,y0,x1,y1] [--all]``: list vector paths (index, stroke
  width, dash array, colours, operators, bounding box). Without ``--all`` the frame-free
  filled boxes and zero-size paths are skipped.
* ``text PDF PAGE [--region ...]``: list words with their boxes, marking bar-notation labels.
* ``images PDF PAGE``: list embedded images (a figure that is an image needs raster mode).
* ``render PDF PAGE OUT.png [--scale 4] [--region ...]``: render a page or region.
* ``run [RECIPE ...]``: run recipes (all when none given), write the reference CSVs and the
  overlays, print a summary.
* ``check [RECIPE ...]``: run recipes in memory and compare with the reference CSVs byte for
  byte; exit status 1 on any difference.
* ``fetch [RECIPE ...]``: download the PDFs recipes name (Kodak hosts only, sha256 checked).
"""

import argparse
import sys
from pathlib import Path

from tools.digitize.figure import DigitizeError, has_bar

__all__ = ["main"]


def _region(text: str | None):
    if text is None:
        return None
    parts = [float(v) for v in text.split(",")]
    if len(parts) != 4:
        raise argparse.ArgumentTypeError("a region is x0,y0,x1,y1 in points")
    return tuple(parts)


def _recipes(names: list[str]) -> list[Path]:
    from tools.digitize.runner import all_recipes

    return [Path(n) for n in names] if names else all_recipes()


def cmd_paths(args) -> int:
    from tools.digitize.page import bbox_contains, fmt_color
    from tools.digitize.pdfread import load_page

    page = load_page(args.pdf, args.page)
    print(f"page {page.number}: {page.width:g} x {page.height:g} pt, {len(page.paths)} paths")
    for path in page.paths:
        x0, y0, x1, y1 = path.bbox
        if args.region and not bbox_contains(args.region, path.bbox, 0.5):
            continue
        if not args.all and not path.stroke and (x1 - x0 < 0.5 or y1 - y0 < 0.5):
            continue
        mode = ("S" if path.stroke else "-") + ("F" if path.fill else "-")
        shape = path.shape if len(path.shape) <= 14 else f"{path.shape[:10]}...[{path.n_segments}]"
        print(
            f"#{path.index:<4d} {mode} {path.style()} fill={fmt_color(path.fill_color)} {shape} "
            f"bbox=({x0:.2f},{y0:.2f},{x1:.2f},{y1:.2f})"
        )
    return 0


def cmd_text(args) -> int:
    from tools.digitize.pdfread import load_page

    page = load_page(args.pdf, args.page)
    words = page.words_in(args.region) if args.region else page.words
    for word in words:
        x0, y0, x1, y1 = word.bbox
        bar = " [bar]" if has_bar(word, page) else ""
        print(
            f"{word.text!r:28s} centre=({word.centre[0]:.2f},{word.centre[1]:.2f}) "
            f"bbox=({x0:.2f},{y0:.2f},{x1:.2f},{y1:.2f}) size={word.size:.2f}{bar}"
        )
    return 0


def cmd_images(args) -> int:
    from tools.digitize.pdfread import load_page

    page = load_page(args.pdf, args.page)
    if not page.images:
        print("no embedded images")
    for box in page.images:
        print("image bbox=({:.2f},{:.2f},{:.2f},{:.2f})".format(*box))
    return 0


def cmd_render(args) -> int:
    from tools.digitize.overlay import render_rgb
    from tools.digitize.pdfread import load_page

    image = render_rgb(args.pdf, args.page, args.scale)
    if args.region:
        height = load_page(args.pdf, args.page).height
        x0, y0, x1, y1 = args.region
        s = args.scale
        image = image.crop(
            (int(x0 * s), int((height - y1) * s), int(x1 * s), int((height - y0) * s))
        )
    image.save(args.out)
    print(f"wrote {args.out} ({image.size[0]} x {image.size[1]})")
    return 0


def _summary(run) -> None:
    result = run.result
    x, y = result.x_axis, result.y_axis
    print(f"{run.path}")
    print(
        f"  x axis: {x.scale}/{x.fit}, {len(x.references)} refs, lsq residual max "
        f"{x.residual_max:.4g}; y axis: {y.scale}/{y.fit}, {len(y.references)} refs, "
        f"lsq residual max {y.residual_max:.4g}"
    )
    for comparison in run.comparisons:
        verdict = "within" if comparison.within else "OUTSIDE"
        print(
            f"  cross-check {comparison.output}: dx {comparison.dx:.4g}, dy {comparison.dy:.4g} "
            f"({verdict} max_difference {list(result.recipe.max_difference)}) vs "
            f"{comparison.staged.name}"
        )
    for curve in result.curves:
        s = curve.stats
        dropped = f"dup {s.duplicates}, same-x {s.same_x}, backtrack {s.backtracking}"
        if curve.edge_dropped:
            dropped += f", on frame edge {curve.edge_dropped}"
        print(
            f"  {curve.spec.output}: {len(curve.x)} pts ({curve.method}, paths "
            f"{list(curve.paths)}; {dropped}); x {curve.x[0]:g}..{curve.x[-1]:g}, "
            f"y {curve.y.min():g}..{curve.y.max():g}; err x "
            f"{x.describe_error(curve.error_x)}, y {y.describe_error(curve.error_y)}"
        )


def cmd_run(args) -> int:
    from tools.digitize.runner import run_file, write_overlay

    for path in _recipes(args.recipes):
        run = run_file(path)
        for target, text in run.csvs.items():
            target.parent.mkdir(parents=True, exist_ok=True)
            with open(target, "w", encoding="utf-8", newline="\n") as handle:
                handle.write(text)
        _summary(run)
        if not args.no_overlay:
            print(f"  overlay: {write_overlay(run, args.scale)}")
    return 0


def cmd_check(args) -> int:
    from tools.digitize.runner import run_file

    failures = 0
    for path in _recipes(args.recipes):
        run = run_file(path)
        for target, text in run.csvs.items():
            current = target.read_text(encoding="utf-8") if target.is_file() else None
            if current != text:
                failures += 1
                state = "missing" if current is None else "differs"
                print(f"{target}: {state} (from {path})")
            else:
                print(f"{target}: ok")
        for comparison in run.comparisons:
            if comparison.within:
                print(f"{path}: cross-check {comparison.output} ok")
            else:
                failures += 1
                print(
                    f"{path}: cross-check {comparison.output} differs by dx {comparison.dx:.4g}, "
                    f"dy {comparison.dy:.4g}"
                )
    if failures:
        print(f"{failures} reference file(s) do not match their recipes")
        return 1
    return 0


def cmd_fetch(args) -> int:
    from tools.digitize.recipe import load_recipe
    from tools.digitize.runner import fetch_pdf

    done = set()
    for path in _recipes(args.recipes):
        recipe = load_recipe(path)
        if recipe.publication.file in done:
            continue
        done.add(recipe.publication.file)
        print(f"{fetch_pdf(recipe)} ok ({recipe.publication.url})")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -I -m tools.digitize",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("paths", help="list vector paths on a page")
    p.add_argument("pdf", type=Path)
    p.add_argument("page", type=int)
    p.add_argument("--region", type=_region)
    p.add_argument("--all", action="store_true")
    p.set_defaults(func=cmd_paths)

    p = sub.add_parser("text", help="list words on a page")
    p.add_argument("pdf", type=Path)
    p.add_argument("page", type=int)
    p.add_argument("--region", type=_region)
    p.set_defaults(func=cmd_text)

    p = sub.add_parser("images", help="list embedded images on a page")
    p.add_argument("pdf", type=Path)
    p.add_argument("page", type=int)
    p.set_defaults(func=cmd_images)

    p = sub.add_parser("render", help="render a page or a region to PNG")
    p.add_argument("pdf", type=Path)
    p.add_argument("page", type=int)
    p.add_argument("out", type=Path)
    p.add_argument("--scale", type=float, default=4.0)
    p.add_argument("--region", type=_region)
    p.set_defaults(func=cmd_render)

    p = sub.add_parser("run", help="run recipes, write reference CSVs and overlays")
    p.add_argument("recipes", nargs="*")
    p.add_argument("--no-overlay", action="store_true")
    p.add_argument("--scale", type=float, default=4.0)
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("check", help="verify reference CSVs reproduce exactly")
    p.add_argument("recipes", nargs="*")
    p.set_defaults(func=cmd_check)

    p = sub.add_parser("fetch", help="download the PDFs recipes name")
    p.add_argument("recipes", nargs="*")
    p.set_defaults(func=cmd_fetch)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except DigitizeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
