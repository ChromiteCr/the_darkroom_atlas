"""File handling: where PDFs, recipes, staged CSVs and overlays live, and running recipes.

Layout (relative to the repository root)::

    .cache/datasheets/pdf/<file>                    publication PDFs (git-ignored)
    .cache/datasheets/overlays/<material>/<n>.png   checking images (git-ignored)
    tools/digitize/recipes/<pub>-<edition>/<n>.yaml recipes
    data/reference/<material>/<output>.csv          digitised curves (<material> = kodak/<film>)
"""

import hashlib
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from tools.digitize.figure import DigitizeError
from tools.digitize.recipe import FigureResult, Recipe, csv_text, load_recipe, run_recipe

__all__ = [
    "ALLOWED_HOSTS",
    "OVERLAY_DIR",
    "PDF_DIR",
    "RECIPE_DIR",
    "REPO",
    "STAGED_DIR",
    "Comparison",
    "Run",
    "fetch_pdf",
    "run_file",
    "sha256_of",
]

REPO = Path(__file__).resolve().parents[2]
PDF_DIR = REPO / ".cache" / "datasheets" / "pdf"
OVERLAY_DIR = REPO / ".cache" / "datasheets" / "overlays"
RECIPE_DIR = REPO / "tools" / "digitize" / "recipes"
STAGED_DIR = REPO / "data" / "reference"

ALLOWED_HOSTS = (
    "kodak.com",
    "kodakalaris.com",
    "kodakprofessional.com",
    "business.kodakmoments.com",
)
"""Hosts (and their subdomains) publications may be fetched from, directly or via the
Wayback Machine (web.archive.org) copies of their URLs."""


def sha256_of(path: Path) -> str:
    """Hex SHA-256 of a file."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _host_allowed(host: str) -> bool:
    host = host.lower().rstrip(".")
    return any(host == h or host.endswith("." + h) for h in ALLOWED_HOSTS)


def url_allowed(url: str) -> bool:
    """True for https/http URLs on a Kodak host, or Wayback copies of such URLs."""
    parts = urllib.parse.urlsplit(url)
    if parts.scheme not in ("http", "https") or not parts.hostname:
        return False
    if parts.hostname == "web.archive.org":
        # /web/<timestamp>[id_]/<original url>
        segments = parts.path.split("/", 3)
        if len(segments) < 4 or segments[1] != "web":
            return False
        original = segments[3]
        if "://" not in original:
            original = "http://" + original
        inner = urllib.parse.urlsplit(original)
        return bool(inner.hostname) and _host_allowed(inner.hostname)
    return _host_allowed(parts.hostname)


def pdf_path(recipe: Recipe, check: bool = True) -> Path:
    """The recipe's PDF; with ``check``, it must exist and match the recorded sha256."""
    path = PDF_DIR / recipe.publication.file
    if not check:
        return path
    if not path.is_file():
        raise DigitizeError(
            f"{path} is missing; fetch it with `python -I -m tools.digitize fetch <recipe>`"
        )
    digest = sha256_of(path)
    if digest != recipe.publication.sha256:
        raise DigitizeError(
            f"{path}: sha256 {digest} does not match the recipe's {recipe.publication.sha256}"
        )
    return path


def fetch_pdf(recipe: Recipe) -> Path:
    """Download the recipe's PDF into the cache unless it is already there and intact."""
    target = PDF_DIR / recipe.publication.file
    if target.is_file() and sha256_of(target) == recipe.publication.sha256:
        return target
    url = recipe.publication.url
    if not url_allowed(url):
        raise DigitizeError(f"refusing to fetch {url}: not a Kodak host or a Wayback copy")
    PDF_DIR.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        final = response.geturl()
        if not url_allowed(final):
            raise DigitizeError(f"{url} redirected to {final}, which is not an allowed host")
        partial.write_bytes(response.read())
    digest = sha256_of(partial)
    if digest != recipe.publication.sha256:
        raise DigitizeError(
            f"{url}: downloaded file has sha256 {digest}, the recipe expects "
            f"{recipe.publication.sha256}; kept as {partial}"
        )
    partial.replace(target)
    return target


@dataclass(frozen=True)
class Comparison:
    """A cross-check curve against the staged curve of the same name.

    ``dx`` is the larger difference of the first and last x; ``dy`` the largest difference
    in y where both curves have data (cross-check interpolated linearly at the staged x).
    Both are in axis-variable units: data units, or decades on a log10 axis.
    """

    output: str
    staged: Path
    dx: float
    dy: float
    within: bool


@dataclass(frozen=True)
class Run:
    """A recipe that has been run: its result and the CSV text of every curve.

    A cross-check recipe has no CSVs; it has comparisons instead.
    """

    path: Path
    result: FigureResult
    csvs: dict[Path, str]
    comparisons: tuple[Comparison, ...] = ()


def recipe_relpath(path: Path) -> str:
    """Path of a recipe relative to the repository root, with forward slashes."""
    try:
        return Path(path).resolve().relative_to(REPO).as_posix()
    except ValueError:
        return Path(path).as_posix()


def staged_csv(material: str, output: str) -> Path:
    """Where the CSV for ``output`` of ``material`` is staged."""
    return STAGED_DIR / material / f"{output}.csv"


def run_file(path: Path) -> Run:
    """Load, check and run one recipe file; nothing is written."""
    from tools.digitize.pdfread import load_page

    recipe = load_recipe(path)
    pdf = pdf_path(recipe)
    page = load_page(pdf, recipe.page)
    gray = None
    if recipe.method == "raster":
        from tools.digitize.raster import render_gray

        gray = render_gray(pdf, recipe.page, recipe.raster_dpi / 72.0)
    result = run_recipe(recipe, page, gray)
    if recipe.cross_check_of is not None:
        return Run(Path(path), result, {}, _compare(result))
    rel = recipe_relpath(path)
    csvs = {
        staged_csv(recipe.material, curve.spec.output): csv_text(result, curve, rel)
        for curve in result.curves
    }
    return Run(Path(path), result, csvs)


def read_staged(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """The x and y columns of a staged CSV."""
    xs, ys = [], []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip() and not line.startswith("#"):
            x, y = line.split(",")
            xs.append(float(x))
            ys.append(float(y))
    return np.asarray(xs), np.asarray(ys)


def _compare(result: FigureResult) -> tuple[Comparison, ...]:
    recipe = result.recipe
    tx = np.log10 if recipe.x.scale == "log10" else np.asarray
    ty = np.log10 if recipe.y.scale == "log10" else np.asarray
    out = []
    for curve in result.curves:
        target = staged_csv(recipe.material, curve.spec.output)
        if not target.is_file():
            raise DigitizeError(f"cross-check: {target} is not staged yet")
        sx, sy = read_staged(target)
        cx, cy = tx(curve.x), ty(curve.y)
        sx, sy = tx(sx), ty(sy)
        dx = max(abs(cx[0] - sx[0]), abs(cx[-1] - sx[-1]))
        inside = (sx >= cx[0]) & (sx <= cx[-1])
        dy = float(np.max(np.abs(np.interp(sx[inside], cx, cy) - sy[inside])))
        limit_x, limit_y = recipe.max_difference
        out.append(
            Comparison(curve.spec.output, target, float(dx), dy, dx <= limit_x and dy <= limit_y)
        )
    return tuple(out)


def write_overlay(run: Run, scale: float = 4.0) -> Path:
    """Write the checking image for a run and return its path."""
    from tools.digitize.overlay import draw_overlay, render_rgb
    from tools.digitize.pdfread import load_page

    recipe = run.result.recipe
    pdf = pdf_path(recipe, check=False)
    page = load_page(pdf, recipe.page)
    image = draw_overlay(render_rgb(pdf, recipe.page, scale), page.height, scale, run.result)
    target = OVERLAY_DIR / recipe.material / f"{run.path.parent.name}-{run.path.stem}.png"
    target.parent.mkdir(parents=True, exist_ok=True)
    image.save(target)
    return target


def all_recipes() -> list[Path]:
    """Every recipe file: recipes that stage curves first, then cross-checks; sorted."""
    paths = sorted(RECIPE_DIR.rglob("*.yaml"))
    return sorted(paths, key=lambda p: load_recipe(p).cross_check_of is not None)
