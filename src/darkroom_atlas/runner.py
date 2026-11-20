"""One whole run: read an image, expose, process and read it out, write the image and manifest.

:func:`run_image` is what ``atlas run`` calls. Everything that can be looked up is looked up
before any work is done or any file is written: the material, the process and its baths,
the curve files the material names (hashed for the manifest), the exposure model, step
handlers and readout device, the input file and the output directory. A run that cannot complete therefore fails at once and leaves nothing behind.

Output files are named ``<input stem>__<material>__<process>.png`` (``/`` in ids becomes
``-``); the manifest has the same name with ``.run.json`` instead of ``.png``.
"""

from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from darkroom_atlas import ATLAS_VERSION
from darkroom_atlas.data.library import _SPECIES_FILE, DataError, Kind, Library
from darkroom_atlas.data.models import Material, Process
from darkroom_atlas.manifest import (
    DataRef,
    FileRef,
    ReadoutRef,
    RunManifest,
    current_environment,
    current_git_commit,
)
from darkroom_atlas.pipeline import RunContext, expose, process, readout
from darkroom_atlas.registry import get_exposure, get_readout, get_step
from darkroom_atlas.scene.camera import scene_light
from darkroom_atlas.scene.io import read_image, write_png16

__all__ = ["DEFAULT_ISO", "MANIFEST_SUFFIX", "manifest_path", "output_name", "run_image"]

DEFAULT_ISO = 100.0
"""ISO speed used for metering when neither the caller nor the material gives one."""

MANIFEST_SUFFIX = ".run.json"
"""Suffix of the manifest file, which replaces the output image's ``.png``."""


def output_name(input_path: Path | str, material_id: str, process_id: str) -> str:
    """File name of the output image: ``<input stem>__<material>__<process>.png``."""
    material = material_id.replace("/", "-")
    process_ = process_id.replace("/", "-")
    return f"{Path(input_path).stem}__{material}__{process_}.png"


def manifest_path(image_path: Path | str) -> Path:
    """Path of the manifest written next to the output image ``image_path``."""
    return Path(image_path).with_suffix(MANIFEST_SUFFIX)


def run_image(
    input_path: Path,
    output_dir: Path,
    *,
    material_id: str,
    process_id: str,
    readout_device: str | None = None,
    illuminant: str = "D65",
    iso: float | None = None,
    ev: float = 0.0,
    time_s: float = 1 / 125,
    seed: int = 0,
    tile_px: int | None = None,
    long_edge_mm: float = 36.0,
    library: Library | None = None,
) -> tuple[Path, RunManifest]:
    """Run ``input_path`` through a material and a process; return the image path and manifest.

    The scene is lit by ``illuminant`` and metered for ``iso`` (default: the material's
    ``iso_speed``, else :data:`DEFAULT_ISO`) at ``ev``, exposed for ``time_s`` seconds,
    processed (tiled with ``tile_px`` when given) and read out with ``readout_device``
    (default: the material's ``default_readout``, else ``ValueError``), which must give an
    ``image_linear_srgb`` readout. The image (:func:`output_name`) and its manifest
    (:func:`manifest_path`) are written to ``output_dir``, which is created if needed;
    existing files of the same name are replaced.

    Raises :class:`~darkroom_atlas.data.library.DataNotFound` for unknown ids,
    :class:`~darkroom_atlas.data.library.DataError` for a curve file the material names but
    the library lacks, :class:`~darkroom_atlas.registry.RegistryError` when a model, step handler or readout
    device is not registered, ``FileNotFoundError`` for a missing input,
    ``NotADirectoryError`` when ``output_dir`` cannot be a directory, and ``ValueError`` for
    invalid settings or input images. All of these happen before anything is written.
    """
    input_path = Path(input_path)
    output_dir = Path(output_dir)
    library = Library() if library is None else library
    ctx = RunContext(library, seed=seed, tile_px=tile_px, long_edge_mm=long_edge_mm)

    material = library.material(material_id)
    proc = library.process(process_id)
    device = _readout_device(material, readout_device)
    bath_ids = _check_runnable(material, proc, device, library)
    iso_value = _iso(material, iso)
    curve_refs = _curve_refs(library, material)
    _check_input(input_path)
    _check_output_dir(output_dir)

    image = read_image(input_path)
    light = scene_light(
        image.linear_rgb,
        grid=ctx.grid,
        illuminant=illuminant,
        iso=iso_value,
        ev=ev,
        time_s=time_s,
    )
    state = process(expose(material, light, time_s, ctx), proc, ctx)
    result = readout(state, device, ctx)
    if result.kind != "image_linear_srgb":
        raise ValueError(
            f"readout device {device!r} gives a {result.kind!r} readout; run_image writes "
            "images and needs an 'image_linear_srgb' readout"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / output_name(input_path, material_id, process_id)
    write_png16(out_path, result.data)

    species_path = library.root / _SPECIES_FILE
    geometry = state.geometry
    manifest = RunManifest(
        atlas_version=ATLAS_VERSION,
        created_utc=datetime.now(UTC).replace(microsecond=0),
        input=_file_ref(input_path, str(input_path)),
        material=_data_ref(library, "material", material_id),
        process=_data_ref(library, "process", process_id),
        baths=[_data_ref(library, "bath", bath_id) for bath_id in bath_ids],
        species_file=(
            _file_ref(species_path, _SPECIES_FILE.as_posix()) if species_path.is_file() else None
        ),
        curve_files=curve_refs,
        readout=ReadoutRef(device=device, params={}),
        exposure=dict(light.meta),
        geometry={
            "shape_hw": list(geometry.shape_hw),
            "long_edge_mm": ctx.long_edge_mm,
            "pixel_pitch_um": geometry.pixel_pitch_um,
        },
        spectral_grid={
            "start_nm": ctx.grid.start_nm,
            "stop_nm": ctx.grid.stop_nm,
            "step_nm": ctx.grid.step_nm,
        },
        seed=ctx.seed,
        tile_px=ctx.tile_px,
        environment=current_environment(),
        git_commit=current_git_commit(),
        output=_file_ref(out_path, str(out_path)),
        warnings=list(image.warnings),
    )
    manifest.write(manifest_path(out_path))
    return out_path, manifest


# --------------------------------------------------------------------------------------
# Checks done before any work
# --------------------------------------------------------------------------------------


def _readout_device(material: Material, readout_device: str | None) -> str:
    if readout_device is not None:
        return readout_device
    if material.default_readout:
        return material.default_readout
    raise ValueError(
        f"material {material.id!r} has no default_readout; name a readout device "
        "(readout_device, or -r/--readout on the command line)"
    )


def _check_runnable(material: Material, proc: Process, device: str, library: Library) -> list[str]:
    """Look up every registered function and bath the run needs; return the bath ids."""
    get_exposure(material.model)
    for step in proc.steps:
        get_step(material.model, step.kind)
    get_readout(device)
    bath_ids = list(dict.fromkeys(step.bath for step in proc.steps if step.bath is not None))
    for bath_id in bath_ids:
        library.bath(bath_id)
    return bath_ids


def _curve_refs(library: Library, material: Material) -> list[FileRef]:
    """Hash every curve file ``material`` names; :class:`DataError` if one is missing."""
    refs = []
    for rel in library.curve_files(material):
        path = library.root / rel
        if not path.is_file():
            raise DataError(
                f"{library.path_of('material', material.id)}: curve file {rel!r} not found "
                f"under {library.root}"
            )
        refs.append(_file_ref(path, rel))
    return refs


def _iso(material: Material, iso: float | None) -> float:
    if iso is not None:
        return iso
    if material.iso_speed is None:
        return DEFAULT_ISO
    value = material.iso_speed.to("1")
    if np.ndim(value) != 0:
        raise ValueError(f"material {material.id!r}: iso_speed must be a single number")
    return float(value)


def _check_input(input_path: Path) -> None:
    if input_path.is_file():
        return
    if input_path.is_dir():
        raise IsADirectoryError(f"{input_path}: the input is a directory, not an image file")
    raise FileNotFoundError(f"{input_path}: input image not found")


def _check_output_dir(output_dir: Path) -> None:
    """Fail unless ``output_dir`` is a directory or could be created as one."""
    for path in (output_dir, *output_dir.parents):
        if not path.exists():
            continue
        if path.is_dir():
            return
        if path == output_dir:
            raise NotADirectoryError(f"{output_dir}: the output path exists and is not a directory")
        raise NotADirectoryError(
            f"{output_dir}: cannot create the output directory, {path} is not a directory"
        )


# --------------------------------------------------------------------------------------
# Manifest parts
# --------------------------------------------------------------------------------------


def _file_ref(path: Path, recorded_as: str) -> FileRef:
    return FileRef(path=recorded_as, sha256=Library.sha256(path))


def _data_ref(library: Library, kind: Kind, id_: str) -> DataRef:
    path = library.path_of(kind, id_)
    return DataRef(id=id_, file=_file_ref(path, path.relative_to(library.root).as_posix()))
