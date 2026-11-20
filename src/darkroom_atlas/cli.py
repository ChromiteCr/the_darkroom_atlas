"""The ``atlas`` command line.

* ``atlas run INPUT -m MATERIAL -p PROCESS -o DIR``: run an image through a material and a
  process (:func:`darkroom_atlas.runner.run_image`) and print where the image and its run
  manifest were written, then any warnings.
* ``atlas check-data``: check every file of the data library and print each problem.

Errors the user can fix (an unknown material, process or readout device, a missing input,
a bad ``--data`` directory, an output path that is a file, an invalid setting) are printed
to stderr as ``error: <message>`` without a traceback, and the command exits with
:data:`EXIT_USAGE` (2). ``check-data`` exits with :data:`EXIT_PROBLEMS` (1) when the library
has problems.
"""

from pathlib import Path
from typing import Annotated, NoReturn

import typer

from darkroom_atlas.data.library import DataError, Library
from darkroom_atlas.registry import RegistryError
from darkroom_atlas.runner import DEFAULT_ISO, manifest_path, run_image

__all__ = ["EXIT_PROBLEMS", "EXIT_USAGE", "app"]

EXIT_USAGE = 2
"""Exit code for errors in what the user asked for (names, paths, settings)."""

EXIT_PROBLEMS = 1
"""Exit code of ``check-data`` when the library has problems."""

app = typer.Typer(
    name="atlas",
    help="The Darkroom Atlas: physically based film exposure, development and readout.",
    no_args_is_help=True,
    add_completion=False,
    pretty_exceptions_show_locals=False,
)

DataOption = Annotated[
    Path | None,
    typer.Option(
        "--data",
        help="Data library directory (default: $DARKROOM_ATLAS_DATA, else the repository's data/).",
        show_default=False,
    ),
]


@app.command("run")
def run_command(
    input_path: Annotated[
        Path, typer.Argument(metavar="INPUT", help="Input image (PNG or TIFF).", show_default=False)
    ],
    material: Annotated[
        str, typer.Option("--material", "-m", help="Material id, e.g. ideal.", show_default=False)
    ],
    process: Annotated[
        str, typer.Option("--process", "-p", help="Process id, e.g. identity.", show_default=False)
    ],
    output: Annotated[
        Path,
        typer.Option(
            "--output", "-o", help="Output directory (created if missing).", show_default=False
        ),
    ],
    readout: Annotated[
        str | None,
        typer.Option(
            "--readout", "-r", help="Readout device (default: the material's).", show_default=False
        ),
    ] = None,
    illuminant: Annotated[
        str, typer.Option("--illuminant", help="Scene illuminant: D65, D50, A or E.")
    ] = "D65",
    iso: Annotated[
        float | None,
        typer.Option(
            "--iso",
            help=f"Metering ISO (default: the material's, else {DEFAULT_ISO:g}).",
            show_default=False,
        ),
    ] = None,
    ev: Annotated[float, typer.Option("--ev", help="Exposure compensation in stops.")] = 0.0,
    time_s: Annotated[float, typer.Option("--time", help="Exposure time in seconds.")] = 1 / 125,
    seed: Annotated[int, typer.Option("--seed", help="Seed of all random numbers.")] = 0,
    tile: Annotated[
        int | None,
        typer.Option(
            "--tile",
            help="Tile size in pixels for processing (default: no tiling).",
            show_default=False,
        ),
    ] = None,
    long_edge_mm: Annotated[
        float, typer.Option("--long-edge-mm", help="Long edge of the frame on the film, mm.")
    ] = 36.0,
    data: DataOption = None,
) -> None:
    """Run an image through a material and a process; write the image and its run manifest."""
    try:
        out_path, manifest = run_image(
            input_path,
            output,
            material_id=material,
            process_id=process,
            readout_device=readout,
            illuminant=illuminant,
            iso=iso,
            ev=ev,
            time_s=time_s,
            seed=seed,
            tile_px=tile,
            long_edge_mm=long_edge_mm,
            library=Library(data),
        )
    except (DataError, RegistryError, OSError, ValueError) as exc:
        _fail(exc)
    typer.echo(f"image:    {out_path}")
    typer.echo(f"manifest: {manifest_path(out_path)}")
    typer.echo(f"recipe:   {manifest.recipe_hash()}")
    for warning in manifest.warnings:
        typer.echo(f"warning: {warning}", err=True)


@app.command("check-data")
def check_data_command(data: DataOption = None) -> None:
    """Check every file of the data library; print each problem, or OK and what was found."""
    try:
        library = Library(data)
    except DataError as exc:
        _fail(exc)
    problems = library.validate()
    if problems:
        for problem in problems:
            typer.echo(problem)
        typer.echo(f"{_count(len(problems), 'problem')} in {library.root}", err=True)
        raise typer.Exit(EXIT_PROBLEMS)
    counts = ", ".join(
        [
            _count(len(library.ids("material")), "material"),
            _count(len(library.ids("bath")), "bath"),
            _count(len(library.ids("process")), "process", "processes"),
            _count(len(library.species()), "species", "species"),
        ]
    )
    typer.echo(f"OK: {library.root}: {counts}")


def _fail(exc: Exception) -> NoReturn:
    message = str(exc).strip() or type(exc).__name__
    for line in message.splitlines():
        typer.echo(f"error: {line}", err=True)
    raise typer.Exit(EXIT_USAGE)


def _count(n: int, singular: str, plural: str | None = None) -> str:
    return f"{n} {singular if n == 1 else (plural or singular + 's')}"
