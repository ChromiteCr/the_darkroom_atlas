"""The run manifest: everything needed to reproduce one output image, written next to it.

A :class:`RunManifest` records the input and output files (with their SHA-256), the data
files that defined the run (material, process, baths, species), the readout device, the
exposure (the light's metadata), the frame geometry, the spectral grid, the seed and the
tile size, plus the software environment and the git commit for reference.

:meth:`RunManifest.recipe_hash` condenses what determines the result into one SHA-256. It
leaves out what may differ between equivalent runs: the creation time, the environment, the
git commit, the tile size (tiling never changes the result) and where the input and output
files are (their hashes are kept). Two runs with the same recipe hash used the same inputs,
data and settings and produced the same output bytes.

Data-file paths are relative to the library root, so the same data in two checkouts gives
the same hash. The free-form ``dict`` fields hold JSON values only: tuples become lists and
NaN or infinity is rejected, so a manifest read back from disk equals the one written.
"""

import hashlib
import json
import math
import platform
import re
import subprocess
from collections.abc import Mapping
from importlib import metadata
from pathlib import Path
from typing import Annotated, Any, Literal, Self

import numpy as np
from pydantic import AfterValidator, AwareDatetime, BaseModel, ConfigDict, Field, ValidationError

__all__ = [
    "DataRef",
    "FileRef",
    "ReadoutRef",
    "RunManifest",
    "current_environment",
    "current_git_commit",
]

_SHA256 = re.compile(r"[0-9a-f]{64}")
_GIT_COMMIT = re.compile(r"[0-9a-f]{40}([0-9a-f]{24})?")
_VOLATILE_FIELDS = ("created_utc", "environment", "git_commit", "tile_px")
_PATH_FIELDS = ("input", "output")


def _plain(value: Any) -> Any:
    """``value`` as plain JSON data (dicts, lists, str, int, float, bool, None)."""
    if isinstance(value, Mapping):
        return {_key(k): _plain(v) for k, v in value.items()}
    if isinstance(value, list | tuple):
        return [_plain(v) for v in value]
    if isinstance(value, np.ndarray):
        return _plain(value.tolist())
    if isinstance(value, np.generic):
        return _plain(value.item())
    if isinstance(value, bool | int | str) or value is None:
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"values must be finite, got {value!r}")
        return float(value)
    raise ValueError(f"values must be JSON data, got a {type(value).__name__}: {value!r}")


def _key(key: Any) -> str:
    # ValueError, not TypeError, so that pydantic reports it as a validation error.
    if not isinstance(key, str):
        raise ValueError(f"keys must be strings, got {key!r}")  # noqa: TRY004
    return key


def _json_dict(value: dict[str, Any]) -> dict[str, Any]:
    return _plain(value)


def _sha256(value: str) -> str:
    if _SHA256.fullmatch(value) is None:
        raise ValueError(f"expected 64 lowercase hexadecimal digits, got {value!r}")
    return value


JsonDict = Annotated[dict[str, Any], AfterValidator(_json_dict)]
Sha256 = Annotated[str, AfterValidator(_sha256)]


class _Model(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")


class FileRef(_Model):
    """A file and the SHA-256 of its bytes."""

    path: str
    sha256: Sha256


class DataRef(_Model):
    """A data-library entry: its id and the file that defined it (path relative to the root)."""

    id: str
    file: FileRef


class ReadoutRef(_Model):
    """The readout device and the parameters it was called with."""

    device: str
    params: JsonDict = Field(default_factory=dict)


class RunManifest(_Model):
    """Record of one run of :func:`darkroom_atlas.runner.run_image`; see the module docs.

    ``exposure`` is the light's metadata (illuminant, ISO, EV, exposure time, metering and
    energy scale), ``geometry`` holds ``shape_hw long_edge_mm pixel_pitch_um``,
    ``spectral_grid`` holds ``start_nm stop_nm step_nm`` and ``environment`` the versions of
    ``python numpy scipy colour`` and the ``platform``.
    """

    schema_version: Literal[1] = 1
    atlas_version: str = Field(min_length=1)
    created_utc: AwareDatetime
    input: FileRef
    material: DataRef
    process: DataRef
    baths: list[DataRef] = Field(default_factory=list)
    species_file: FileRef | None = None
    readout: ReadoutRef
    exposure: JsonDict
    geometry: JsonDict
    spectral_grid: JsonDict
    seed: int = Field(ge=0)
    tile_px: int | None = Field(default=None, ge=1)
    environment: JsonDict
    git_commit: str | None = None
    output: FileRef
    warnings: list[str] = Field(default_factory=list)

    def recipe_hash(self) -> str:
        """SHA-256 (hex) of the canonical JSON of everything that determines the result.

        Excluded: ``created_utc``, ``environment``, ``git_commit``, ``tile_px`` and the paths
        (not the hashes) of ``input`` and ``output``.
        """
        data = self.model_dump(mode="json")
        for name in _VOLATILE_FIELDS:
            del data[name]
        for name in _PATH_FIELDS:
            del data[name]["path"]
        canonical = json.dumps(
            data, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
        )
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    def write(self, path: Path | str) -> None:
        """Write the manifest to ``path`` as indented UTF-8 JSON."""
        Path(path).write_text(self.model_dump_json(indent=2) + "\n", encoding="utf-8")

    @classmethod
    def read(cls, path: Path | str) -> Self:
        """Read a manifest written by :meth:`write`; ``ValueError`` naming ``path`` if invalid."""
        text = Path(path).read_text(encoding="utf-8")
        try:
            return cls.model_validate_json(text)
        except ValidationError as exc:
            problems = "; ".join(
                f"{'.'.join(map(str, error['loc'])) or '(top level)'}: {error['msg']}"
                for error in exc.errors(include_url=False)
            )
            raise ValueError(f"{path}: not a valid run manifest: {problems}") from None


# --------------------------------------------------------------------------------------
# Environment
# --------------------------------------------------------------------------------------


def current_environment() -> dict[str, str]:
    """Versions of Python, NumPy, SciPy and colour-science, and the platform string."""
    return {
        "python": platform.python_version(),
        "numpy": _version("numpy"),
        "scipy": _version("scipy"),
        "colour": _version("colour-science"),
        "platform": platform.platform(),
    }


def current_git_commit() -> str | None:
    """``git rev-parse HEAD`` of the checkout holding this package; ``None`` if unavailable."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(__file__).resolve().parent,
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    commit = result.stdout.strip()
    if result.returncode != 0 or _GIT_COMMIT.fullmatch(commit) is None:
        return None
    return commit


def _version(distribution: str) -> str:
    try:
        return metadata.version(distribution)
    except metadata.PackageNotFoundError:
        return "unknown"
