"""One whole run: read an image, expose, process and read it out, write the image and manifest.

:func:`run_image` is what ``atlas run`` calls. Everything that can be looked up is looked up
before any work is done or any file is written: the material, the process and its baths,
the curve files the material names (hashed for the manifest), the exposure model, step
handlers and readout device, the input file and the output directory. A run that cannot
complete therefore fails at once and leaves nothing behind.

D1 (spec D1 4.5, 4.7): the scene light of a model other than ``ideal`` is low-rank (its
coefficients are the image itself); the run gets one scratch directory, a fresh subdirectory
of ``scratch_dir`` (default: the system temporary directory) holding its memory-mapped
fields, which is deleted at the end unless ``keep_scratch``; the number of worker processes
follows the memory budget (:func:`~darkroom_atlas.blocks.workers_for_budget`).

Output files are named ``<input stem>__<material>__<process>.png`` (``/`` in ids becomes
``-``); the manifest has the same name with ``.run.json`` instead of ``.png``.
"""

import shutil
import tempfile
from dataclasses import dataclass, replace
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

import numpy as np

from darkroom_atlas import ATLAS_VERSION
from darkroom_atlas.blocks import workers_for_budget
from darkroom_atlas.data.library import _KINETICS_FILE, _SPECIES_FILE, DataError, Kind, Library
from darkroom_atlas.data.models import Material, Process
from darkroom_atlas.fieldstore import FieldStore
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
from darkroom_atlas.scene.io import ImageInput, read_image, write_png16

__all__ = [
    "DEFAULT_ISO",
    "MANIFEST_SUFFIX",
    "SCRATCH_PREFIX",
    "manifest_path",
    "output_name",
    "run_image",
]

DEFAULT_ISO = 100.0
"""ISO speed used for metering when neither the caller nor the material gives one."""

MANIFEST_SUFFIX = ".run.json"
"""Suffix of the manifest file, which replaces the output image's ``.png``."""

SCRATCH_PREFIX = "atlas-run-"
"""Name prefix of a run's scratch directory."""

_IDEAL_MODEL = "ideal"


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
    film_format: str | None = None,
    profile: Literal["image", "calibration", "reference"] = "image",
    workers: int | None = None,
    memory_budget_gb: float = 6.0,
    scratch_dir: Path | None = None,
    keep_scratch: bool = False,
) -> tuple[Path, RunManifest]:
    """Run ``input_path`` through a material and a process; return the image path and manifest.

    The scene is lit by ``illuminant`` and metered for ``iso`` (default: the material's
    ``iso_speed``, else :data:`DEFAULT_ISO`) at ``ev``, exposed for ``time_s`` seconds,
    processed (tiled with ``tile_px`` when given) and read out with ``readout_device``
    (default: the material's ``default_readout``, else ``ValueError``), which must give an
    ``image_linear_srgb`` readout. The image (:func:`output_name`) and its manifest
    (:func:`manifest_path`) are written to ``output_dir``, which is created if needed;
    existing files of the same name are replaced.

    ``film_format`` (default: the material's) and ``profile`` (spec D1 3.1) configure the
    ``bw-mechanistic`` model; ``workers`` (default: as many as ``memory_budget_gb`` allows)
    is the number of worker processes. Memory-mapped fields live in a new subdirectory of
    ``scratch_dir`` (default: the system temporary directory), removed at the end unless
    ``keep_scratch``.

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
    if not isinstance(keep_scratch, bool):
        raise TypeError(f"keep_scratch must be True or False, got {keep_scratch!r}")
    n_workers = workers_for_budget(memory_budget_gb, workers)
    ctx = RunContext(
        library,
        seed=seed,
        tile_px=tile_px,
        long_edge_mm=long_edge_mm,
        film_format=film_format,
        profile=profile,
        workers=n_workers,
        memory_budget_gb=memory_budget_gb,
    )

    material = library.material(material_id)
    proc = library.process(process_id)
    device = _readout_device(material, readout_device)
    bath_ids = _check_runnable(material, proc, device, library)
    _check_format(material, film_format)
    iso_value = _iso(material, iso)
    curve_refs = _curve_refs(library, material)
    _check_input(input_path)
    _check_output_dir(output_dir)
    image = read_image(input_path)

    if scratch_dir is not None:
        Path(scratch_dir).mkdir(parents=True, exist_ok=True)
    run_dir = Path(tempfile.mkdtemp(prefix=SCRATCH_PREFIX, dir=scratch_dir))
    store = FieldStore(run_dir / "fields")
    try:
        ctx = replace(ctx, scratch_dir=run_dir, store=store)
        out_path, manifest = _run(
            _Plan(input_path, output_dir, image, material, proc, device, bath_ids, curve_refs),
            illuminant=illuminant,
            iso=iso_value,
            ev=ev,
            time_s=time_s,
            ctx=ctx,
        )
    finally:
        store.close()
        if not keep_scratch:
            shutil.rmtree(run_dir, ignore_errors=True)
    return out_path, manifest


@dataclass(frozen=True)
class _Plan:
    """What the checks before any work established."""

    input_path: Path
    output_dir: Path
    image: ImageInput
    material: Material
    proc: Process
    device: str
    bath_ids: list[str]
    curve_refs: list[FileRef]


def _run(
    plan: _Plan, *, illuminant: str, iso: float, ev: float, time_s: float, ctx: RunContext
) -> tuple[Path, RunManifest]:
    """Expose, process and read out; write the image and the manifest."""
    input_path, output_dir, image = plan.input_path, plan.output_dir, plan.image
    material, proc, device = plan.material, plan.proc, plan.device
    library = ctx.library
    light = scene_light(
        image.linear_rgb,
        grid=ctx.grid,
        illuminant=illuminant,
        iso=iso,
        ev=ev,
        time_s=time_s,
        low_rank=material.model != _IDEAL_MODEL,
    )
    exposure_meta = dict(light.meta)
    exposed = expose(material, light, time_s, ctx)
    del light  # spec D1 4.5: the light is dropped after exposure
    state = process(exposed, proc, ctx)
    del exposed
    result = readout(state, device, ctx)
    if result.kind != "image_linear_srgb":
        raise ValueError(
            f"readout device {device!r} gives a {result.kind!r} readout; run_image writes "
            "images and needs an 'image_linear_srgb' readout"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / output_name(input_path, material.id, proc.id)
    write_png16(out_path, result.data)

    species_path = library.root / _SPECIES_FILE
    kinetics_path = library.root / _KINETICS_FILE
    geometry = state.geometry
    mechanistic = material.model != _IDEAL_MODEL
    exposure_state = state.meta.get("exposure") or {}
    manifest = RunManifest(
        atlas_version=ATLAS_VERSION,
        created_utc=datetime.now(UTC).replace(microsecond=0),
        input=_file_ref(input_path, str(input_path)),
        material=_data_ref(library, "material", material.id),
        process=_data_ref(library, "process", proc.id),
        baths=[_data_ref(library, "bath", bath_id) for bath_id in plan.bath_ids],
        species_file=(
            _file_ref(species_path, _SPECIES_FILE.as_posix()) if species_path.is_file() else None
        ),
        curve_files=plan.curve_refs,
        readout=ReadoutRef(device=device, params={}),
        exposure=exposure_meta,
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
        film_format=exposure_state.get("film_format", ctx.film_format),
        model_settings=_model_settings(state, ctx) if mechanistic else {},
        kinetics_file=(
            _file_ref(kinetics_path, _KINETICS_FILE.as_posix()) if kinetics_path.is_file() else None
        ),
        calibration_runs=_calibration_runs(library) if mechanistic else [],
        workers=ctx.workers,
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


def _check_format(material: Material, film_format: str | None) -> None:
    if film_format is None or film_format in material.formats:
        return
    raise ValueError(
        f"material {material.id!r}: unknown film format {film_format!r}; available: "
        f"{', '.join(material.formats) or '(none)'}"
    )


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


def _model_settings(state, ctx: RunContext) -> dict:
    """Discretisation of a ``bw-mechanistic`` run (spec D1 3.1, 3.6)."""
    from darkroom_atlas.develop import step as dstep
    from darkroom_atlas.develop.grid import PROFILE_LEVELS
    from darkroom_atlas.models.bw import PROFILE_CLASSES

    n_z, n_l, r_l = PROFILE_LEVELS[ctx.profile]
    h0, ratio, h_max = dstep.PROFILE_SCHEDULES[ctx.profile]
    develop = state.meta.get("develop") or {}
    return {
        "profile": ctx.profile,
        "n_z": n_z,
        "n_l": n_l,
        "r_l": r_l,
        "size_classes": PROFILE_CLASSES[ctx.profile],
        "time_steps": {"h0_s": h0, "ratio": ratio, "h_max_s": h_max},
        "develop_block_cells": dstep.BLOCK_CELLS,
        "develop_halo_um": develop.get("halo_um"),
    }


def _calibration_runs(library: Library) -> list[str]:
    """Run ids of the calibrations that wrote fitted values into the kinetics file."""
    if not (library.root / _KINETICS_FILE).is_file():
        return []
    calibration = library.kinetics().calibration
    return [] if calibration is None else [calibration.run_id]


def _file_ref(path: Path, recorded_as: str) -> FileRef:
    return FileRef(path=recorded_as, sha256=Library.sha256(path))


def _data_ref(library: Library, kind: Kind, id_: str) -> DataRef:
    path = library.path_of(kind, id_)
    return DataRef(id=id_, file=_file_ref(path, path.relative_to(library.root).as_posix()))
