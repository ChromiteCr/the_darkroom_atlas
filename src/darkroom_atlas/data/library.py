"""Loading and checking the data library.

Layout of a library root (``data/`` in the repository by default)::

    chemistry/species.yaml          every chemical species
    chemistry/kinetics.yaml         shared kinetics constants (D1)
    developers/<id>.yaml            developers           \\  both hold baths;
    baths/<id>.yaml                 other baths          /   developers/ is searched first
    films/<id>.yaml                 films                \\  both hold materials;
    papers/<id>.yaml                papers               /   films/ is searched first
    processes/<id>.yaml             processes
    reference/<material id>/*.csv   reference curves
    reference/<material id>/tables.yaml   datasheet tables and facts (D1)
    reference/developers/<name>.yaml      developer-publication facts (D1)
    calibration/<run id>.yaml       calibration records (D1)

Ids may contain ``/`` (``kodak/tri-x-400`` lives in ``films/kodak/tri-x-400.yaml``), and
the ``id`` written in a file must equal the id implied by its path.

Reference curve CSV files start with ``# key: value`` header lines (keys ``material quantity
x_unit y_unit origin ref conditions``), followed by ``x,y`` data rows without a header row.
Files named by :attr:`SpectralCurve.file` use the same format; their header is optional and
their ``x`` values are wavelengths in nm unless ``x_unit`` names another length unit.

Every problem is reported as a :class:`DataError` whose message lines read
``<file path>: <field path>: <problem>``.
"""

import hashlib
import os
from collections.abc import Hashable, Iterator
from pathlib import Path, PurePosixPath
from typing import Any, Literal

import numpy as np
import yaml
from pydantic import BaseModel, ValidationError

from darkroom_atlas.data.models import (
    BW_MECHANISTIC,
    DATA_ID_PATTERN,
    NAME_PATTERN,
    Bath,
    CalibrationRecord,
    KineticsFile,
    Material,
    Process,
    ReferenceCurve,
    Species,
    SpeciesFile,
    SpectralCurve,
    bw_mechanistic_problems,
    iter_quantities,
)
from darkroom_atlas.data.provenance import Source
from darkroom_atlas.spectral import DEFAULT_GRID, SpectralGrid, cmf_xyz
from darkroom_atlas.units import convert

__all__ = [
    "DEVELOPER_FACTS_FILES",
    "ENV_DATA_ROOT",
    "DataError",
    "DataNotFound",
    "Library",
    "default_data_root",
]

ENV_DATA_ROOT = "DARKROOM_ATLAS_DATA"
"""Environment variable that overrides the default data root."""

Kind = Literal["material", "bath", "process"]

_KIND_DIRS: dict[str, tuple[str, ...]] = {
    "material": ("films", "papers"),
    "bath": ("developers", "baths"),
    "process": ("processes",),
}
_KIND_MODELS: dict[str, type[BaseModel]] = {
    "material": Material,
    "bath": Bath,
    "process": Process,
}
_SPECIES_FILE = PurePosixPath("chemistry/species.yaml")
_KINETICS_FILE = PurePosixPath("chemistry/kinetics.yaml")
_CALIBRATION_DIR = "calibration"
_REFERENCE_DIR = "reference"
_TABLES_FILE = "tables.yaml"
_DEVELOPER_FACTS_DIR = "developers"
# Built-in pseudo-components added by chemistry.species.component_table; they are not
# entries of species.yaml (spec D1 2.3.1, 5.2).
_BUILTIN_COMPONENTS = frozenset({"toth", "gelatin"})

DEVELOPER_FACTS_FILES: dict[str, str] = {
    "d-76": "d-76",
    "xtol": "xtol",
    "hc-110-pre2019": "hc-110",
    "t-max-dev": "t-max-dev",
}
"""Bath id -> file name (without ``.yaml``) under ``reference/developers/`` (spec D1 5.9).

The 2017 J-24 (HC-110) predates the 2019 reformulation, so its facts belong to
``hc-110-pre2019``; ``t-max-rs`` is not used in D1.
"""
_CSV_HEADER_KEYS = ("material", "quantity", "x_unit", "y_unit", "origin", "ref", "conditions")
_BUILTIN_CMF_COLUMN = {"cie1931_x": 0, "cie1931_y": 1, "cie1931_z": 2}


class DataError(Exception):
    """A data file is missing, malformed or inconsistent.

    The message starts with the path of the offending file; each problem is on its own line.
    """


class DataNotFound(DataError):
    """A requested id or file does not exist. The message lists what is available."""


def default_data_root() -> Path:
    """The data root: ``$DARKROOM_ATLAS_DATA`` if set, else the repository's ``data/``."""
    override = os.environ.get(ENV_DATA_ROOT)
    if override:
        return Path(override).expanduser()
    return Path(__file__).resolve().parents[3] / "data"


# --------------------------------------------------------------------------------------
# Reading files
# --------------------------------------------------------------------------------------


class _StrictLoader(yaml.SafeLoader):
    """Safe YAML loader that rejects duplicate keys instead of keeping the last one."""


def _construct_mapping(loader: _StrictLoader, node: yaml.MappingNode, deep: bool = False):
    mapping: dict[Any, Any] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if not isinstance(key, Hashable):
            raise yaml.constructor.ConstructorError(
                "while constructing a mapping", node.start_mark, "found unhashable key",
                key_node.start_mark,
            )  # fmt: skip
        if key in mapping:
            raise yaml.constructor.ConstructorError(
                "while constructing a mapping", node.start_mark, f"duplicate key {key!r}",
                key_node.start_mark,
            )  # fmt: skip
        mapping[key] = loader.construct_object(value_node, deep=deep)
    return mapping


_StrictLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_mapping)


def _read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8-sig")
    except UnicodeDecodeError:
        raise DataError(f"{path}: file is not valid UTF-8 text") from None
    except OSError as exc:
        raise DataError(f"{path}: cannot read file: {exc.strerror or exc}") from None


def _read_yaml(path: Path) -> dict[str, Any]:
    text = _read_text(path)
    try:
        data = yaml.load(text, Loader=_StrictLoader)
    except yaml.MarkedYAMLError as exc:
        mark = exc.problem_mark or exc.context_mark
        where = f" at line {mark.line + 1}, column {mark.column + 1}" if mark else ""
        what = "; ".join(part for part in (exc.context, exc.problem) if part)
        raise DataError(f"{path}: YAML error{where}: {what}") from None
    except yaml.YAMLError as exc:
        message = " ".join(str(exc).split())
        raise DataError(f"{path}: YAML error: {message}") from None
    if data is None:
        raise DataError(f"{path}: file is empty")
    if not isinstance(data, dict):
        raise DataError(
            f"{path}: the top level must be a mapping of field names to values, "
            f"got a {type(data).__name__}"
        )
    return data


def _format_loc(loc: tuple[int | str, ...]) -> str:
    text = ""
    for part in loc:
        if isinstance(part, int):
            text += f"[{part}]"
        else:
            text += f".{part}" if text else str(part)
    return text


def _format_validation_error(path: Path, exc: ValidationError) -> str:
    lines = []
    for error in exc.errors(include_url=False):
        message = error["msg"]
        if error["type"] == "value_error":
            message = message.removeprefix("Value error, ")
        elif error["type"] == "extra_forbidden":
            message = "unknown field"
        loc = _format_loc(error["loc"])
        lines.append(f"{path}: {loc}: {message}" if loc else f"{path}: {message}")
    return "\n".join(lines)


def _validate_model[M: BaseModel](model: type[M], data: Any, path: Path) -> M:
    try:
        return model.model_validate(data)
    except ValidationError as exc:
        raise DataError(_format_validation_error(path, exc)) from None


def _read_curve_csv(path: Path) -> tuple[dict[str, str], list[float], list[float]]:
    """Read a curve CSV: ``# key: value`` header lines, then ``x,y`` rows."""
    header: dict[str, str] = {}
    xs: list[float] = []
    ys: list[float] = []
    for number, raw in enumerate(_read_text(path).splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            if xs:
                raise DataError(f"{path}: line {number}: header lines must precede the data")
            key, sep, value = line[1:].partition(":")
            key, value = key.strip(), value.strip()
            if not sep or not key:
                raise DataError(f"{path}: line {number}: expected a header line '# key: value'")
            if key not in _CSV_HEADER_KEYS:
                raise DataError(
                    f"{path}: line {number}: unknown header key {key!r}; "
                    f"allowed: {', '.join(_CSV_HEADER_KEYS)}"
                )
            if key in header:
                raise DataError(f"{path}: line {number}: duplicate header key {key!r}")
            if value:
                header[key] = value
            continue
        fields = line.split(",")
        try:
            if len(fields) != 2:
                raise ValueError
            x, y = float(fields[0]), float(fields[1])
        except ValueError:
            raise DataError(
                f"{path}: line {number}: expected two comma-separated numbers 'x,y', got {raw!r}"
            ) from None
        if not (np.isfinite(x) and np.isfinite(y)):
            raise DataError(f"{path}: line {number}: values must be finite, got {raw!r}")
        xs.append(x)
        ys.append(y)
    if not xs:
        raise DataError(f"{path}: no data rows")
    return header, xs, ys


def _available(names: list[str]) -> str:
    return ", ".join(names) if names else "(none)"


def _material_curves(material: Material) -> Iterator[tuple[str, SpectralCurve]]:
    """Yield ``(field path, curve)`` for every spectral curve in ``material``."""
    for i, layer in enumerate(material.layers):
        for name in ("spectral_sensitivity", "spectral_density"):
            curve = getattr(layer, name, None)
            if curve is not None:
                yield f"layers[{i}].{name}", curve
        for j, coupler in enumerate(getattr(layer, "couplers", ())):
            if coupler.dye_spectral_density is not None:
                yield (
                    f"layers[{i}].couplers[{j}].dye_spectral_density",
                    coupler.dye_spectral_density,
                )


def _check_reference_yaml(path: Path, data: dict[str, Any], key: str, expected: str) -> None:
    version = data.get("schema_version", 1)
    if version != 1:
        raise DataError(f"{path}: schema_version: expected 1, got {version!r}")
    declared = data.get(key, expected)
    if declared != expected:
        raise DataError(
            f"{path}: {key}: the file declares {declared!r} but its path implies {expected!r}"
        )


# --------------------------------------------------------------------------------------
# Library
# --------------------------------------------------------------------------------------


class Library:
    """A data library rooted at ``root`` (default: :func:`default_data_root`).

    Files are read on every call, so a library always reflects what is on disk. A ``root``
    that is not a directory raises :class:`DataNotFound`.
    """

    def __init__(self, root: Path | str | None = None) -> None:
        self._root = Path(root) if root is not None else default_data_root()
        if not self._root.is_dir():
            raise DataNotFound(f"{self._root}: data directory does not exist")

    @property
    def root(self) -> Path:
        """The library root directory."""
        return self._root

    def __repr__(self) -> str:
        return f"Library({str(self._root)!r})"

    # ---- listing ---------------------------------------------------------------------

    @staticmethod
    def _dirs(kind: str) -> tuple[str, ...]:
        try:
            return _KIND_DIRS[kind]
        except (KeyError, TypeError):
            raise ValueError(
                f"unknown kind {kind!r}; expected one of: {', '.join(_KIND_DIRS)}"
            ) from None

    def _files(self, kind: str, suffix: str = ".yaml") -> list[tuple[str, Path]]:
        """``(id, path)`` for every file of ``kind``, directory by directory, ids sorted."""
        found = []
        for directory in self._dirs(kind):
            base = self._root / directory
            if not base.is_dir():
                continue
            entries = []
            for path in base.rglob(f"*{suffix}"):
                rel = path.relative_to(base)
                if any(part.startswith(".") for part in rel.parts) or not path.is_file():
                    continue
                entries.append((rel.with_suffix("").as_posix(), path))
            found.extend(sorted(entries))
        return found

    def ids(self, kind: Kind) -> list[str]:
        """Sorted ids of every ``kind`` (``"material"``, ``"bath"`` or ``"process"``) file."""
        return sorted({id_ for id_, _ in self._files(kind)})

    def path_of(self, kind: Kind, id: str) -> Path:
        """Path of the file defining ``id``; raises :class:`DataNotFound` if there is none."""
        files = self._files(kind)
        for file_id, path in files:
            if file_id == id:
                return path
        looked_in = " and ".join(f"{d}/" for d in self._dirs(kind))
        available = sorted({file_id for file_id, _ in files})
        raise DataNotFound(
            f"{self._root}: no {kind} {id!r} in {looked_in}; available: {_available(available)}"
        )

    # ---- loading ---------------------------------------------------------------------

    def _load_file(self, kind: str, path: Path, expected_id: str) -> Any:
        obj = _validate_model(_KIND_MODELS[kind], _read_yaml(path), path)
        if obj.id != expected_id:
            raise DataError(
                f"{path}: id: the file declares id {obj.id!r} but its path implies {expected_id!r}"
            )
        return obj

    def material(self, id: str) -> Material:
        """Load ``films/<id>.yaml``, or ``papers/<id>.yaml`` if there is no such film."""
        return self._load_file("material", self.path_of("material", id), id)

    def bath(self, id: str) -> Bath:
        """Load ``developers/<id>.yaml``, or ``baths/<id>.yaml`` if there is no such developer."""
        return self._load_file("bath", self.path_of("bath", id), id)

    def process(self, id: str) -> Process:
        """Load ``processes/<id>.yaml``."""
        return self._load_file("process", self.path_of("process", id), id)

    def species(self) -> dict[str, Species]:
        """All species by id; empty if ``chemistry/species.yaml`` does not exist."""
        path = self._root / _SPECIES_FILE
        if not path.is_file():
            return {}
        species_file = _validate_model(SpeciesFile, _read_yaml(path), path)
        return {s.id: s for s in species_file.species}

    def kinetics(self) -> KineticsFile:
        """Load ``chemistry/kinetics.yaml``; raises :class:`DataNotFound` if it does not exist."""
        path = self._root / _KINETICS_FILE
        if not path.is_file():
            raise DataNotFound(f"{path}: kinetics file does not exist")
        return _validate_model(KineticsFile, _read_yaml(path), path)

    def _calibration_files(self) -> list[tuple[str, Path]]:
        base = self._root / _CALIBRATION_DIR
        if not base.is_dir():
            return []
        return sorted(
            (path.stem, path)
            for path in base.glob("*.yaml")
            if path.is_file() and not path.name.startswith(".")
        )

    def calibration_record(self, run_id: str) -> CalibrationRecord:
        """Load ``calibration/<run_id>.yaml``; its ``run_id`` must equal the file name."""
        path = self._root / _CALIBRATION_DIR / f"{run_id}.yaml"
        if not (isinstance(run_id, str) and NAME_PATTERN.fullmatch(run_id) and path.is_file()):
            available = [name for name, _ in self._calibration_files()]
            raise DataNotFound(
                f"{self._root / _CALIBRATION_DIR}: no calibration record {run_id!r}; "
                f"available: {_available(available)}"
            )
        record = _validate_model(CalibrationRecord, _read_yaml(path), path)
        if record.run_id != run_id:
            raise DataError(
                f"{path}: run_id: the file declares {record.run_id!r} but its name implies "
                f"{run_id!r}"
            )
        return record

    def _tables_files(self) -> list[tuple[str, Path]]:
        """``(material id, path)`` of every ``reference/<material id>/tables.yaml``."""
        base = self._root / _REFERENCE_DIR
        if not base.is_dir():
            return []
        entries = []
        for path in base.rglob(_TABLES_FILE):
            rel = path.relative_to(base)
            if any(part.startswith(".") for part in rel.parts) or not path.is_file():
                continue
            if rel.parts[0] == _DEVELOPER_FACTS_DIR:
                continue
            entries.append((rel.parent.as_posix(), path))
        return sorted(entries)

    def reference_tables(self, material_id: str) -> dict:
        """Read ``reference/<material_id>/tables.yaml`` as plain data (spec D1 5.9).

        The structure is the digitising team's (``tools/digitize/staged/kodak/<film>/``);
        only ``schema_version`` (1 if given) and ``material`` (equal to ``material_id`` if
        given) are checked.
        """
        path = self._root / _REFERENCE_DIR / material_id / _TABLES_FILE
        safe = isinstance(material_id, str) and DATA_ID_PATTERN.fullmatch(material_id)
        if not (safe and path.is_file()):
            available = [id_ for id_, _ in self._tables_files()]
            raise DataNotFound(
                f"{self._root / _REFERENCE_DIR}: no reference tables for {material_id!r}; "
                f"available: {_available(available)}"
            )
        data = _read_yaml(path)
        _check_reference_yaml(path, data, "material", material_id)
        return data

    def developer_facts(self, developer: str) -> dict:
        """Read the developer-publication facts of bath ``developer`` as plain data.

        ``developer`` is a bath id of :data:`DEVELOPER_FACTS_FILES` (``hc-110-pre2019`` reads
        ``reference/developers/hc-110.yaml``). The structure is the one described in
        ``tools/digitize/staged/developers/README.md``; only ``schema_version`` (1 if given)
        and ``developer`` (equal to the file name if given) are checked.
        """
        name = DEVELOPER_FACTS_FILES.get(developer) if isinstance(developer, str) else None
        base = self._root / _REFERENCE_DIR / _DEVELOPER_FACTS_DIR
        path = base / f"{name}.yaml"
        if name is None or not path.is_file():
            available = [
                b for b, f in DEVELOPER_FACTS_FILES.items() if (base / f"{f}.yaml").is_file()
            ]
            raise DataNotFound(
                f"{base}: no developer facts for {developer!r}; available: {_available(available)}"
            )
        data = _read_yaml(path)
        _check_reference_yaml(path, data, "developer", name)
        return data

    def _reference_files(self) -> list[tuple[str, Path]]:
        base = self._root / _REFERENCE_DIR
        if not base.is_dir():
            return []
        entries = []
        for path in base.rglob("*.csv"):
            rel = path.relative_to(base)
            if any(part.startswith(".") for part in rel.parts) or not path.is_file():
                continue
            entries.append((rel.as_posix(), path))
        return sorted(entries)

    def reference_curve(self, relpath: str | PurePosixPath) -> ReferenceCurve:
        """Load ``reference/<relpath>``; the header's ``material`` must match its directory."""
        rel = PurePosixPath(str(relpath))
        path = self._root / _REFERENCE_DIR / rel
        safe = str(relpath).strip() and not rel.is_absolute() and ".." not in rel.parts
        if not (safe and path.is_file()):
            available = [name for name, _ in self._reference_files()]
            raise DataNotFound(
                f"{self._root / _REFERENCE_DIR}: no reference curve {str(relpath)!r}; "
                f"available: {_available(available)}"
            )
        header, xs, ys = _read_curve_csv(path)
        curve = _validate_model(ReferenceCurve, {**header, "x": xs, "y": ys}, path)
        implied = rel.parent.as_posix()
        if implied == ".":
            raise DataError(
                f"{path}: reference curves must be stored under {_REFERENCE_DIR}/<material id>/"
            )
        if curve.material != implied:
            raise DataError(
                f"{path}: material: the header declares {curve.material!r} but the file's "
                f"directory implies {implied!r}"
            )
        return curve

    # ---- curves ----------------------------------------------------------------------

    def resolve_curve(self, curve: SpectralCurve, grid: SpectralGrid) -> np.ndarray:
        """Sample ``curve`` at the band centres of ``grid``, shape ``(B,)``.

        Builtin curves are the CIE 1931 colour-matching functions. Inline and file curves are
        interpolated linearly and are 0 outside the wavelengths they cover.
        """
        if curve.builtin is not None:
            return cmf_xyz(grid)[:, _BUILTIN_CMF_COLUMN[curve.builtin]].copy()
        if curve.inline is not None:
            return grid.resample(curve.inline.wavelengths_nm, curve.inline.values, outside="zero")

        path = self._root / curve.file
        if not path.is_file():
            raise DataError(f"{path}: spectral curve file not found")
        header, xs, ys = _read_curve_csv(path)
        try:
            wavelengths = convert(np.asarray(xs), header.get("x_unit", "nm"), "nm")
            return grid.resample(wavelengths, ys, outside="zero")
        except ValueError as exc:
            raise DataError(f"{path}: {exc}") from None

    @staticmethod
    def curve_files(material: Material) -> list[str]:
        """Sorted paths, relative to the root, of every curve file that ``material`` names.

        These are the :attr:`SpectralCurve.file` values of all its layers and couplers
        (``./`` and repeated ``/`` removed, duplicates listed once), whether or not the
        material's model reads them.
        """
        return sorted(
            {
                PurePosixPath(curve.file).as_posix()
                for _, curve in _material_curves(material)
                if curve.file is not None
            }
        )

    # ---- hashing ---------------------------------------------------------------------

    @staticmethod
    def sha256(path: Path | str) -> str:
        """Hex SHA-256 digest of the file at ``path``."""
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        return digest.hexdigest()

    # ---- whole-library check ---------------------------------------------------------

    def validate(self) -> list[str]:
        """Check every file in the library and return one line per problem (empty if none).

        Reports files that fail to load, ids defined twice, ``.yml`` files (which are not
        read), reference curves for unknown materials, steps that name unknown baths, bath
        components that name unknown species, and curve files that are missing or malformed.
        Never stops at the first problem.

        D1 checks (spec D1 section 5): ingredients whose ``dissolves_to`` names something
        other than a component; developing agents whose ``kinetic_class`` is not an agent of
        ``chemistry/kinetics.yaml`` (which must then load); ``bw-mechanistic`` materials
        missing a field the model needs; ``reference/<material>/tables.yaml`` and
        ``reference/developers/*.yaml`` that fail to load or name another material or
        developer; calibration records that fail to load or whose ``pair_statuses`` name an
        unknown material or bath; fitted values whose ``fit.run_id`` names no calibration
        record; and species files that define the built-in pseudo-components ``toth`` or
        ``gelatin``.
        """
        problems: list[str] = []
        loaded: dict[str, list[tuple[Path, Any]]] = {kind: [] for kind in _KIND_DIRS}
        known: dict[str, list[str]] = {}

        for kind in _KIND_DIRS:
            first_path: dict[str, Path] = {}
            for id_, path in self._files(kind):
                if id_ in first_path:
                    problems.append(
                        f"{path}: id: {kind} {id_!r} is also defined by {first_path[id_]}"
                    )
                    continue
                first_path[id_] = path
                try:
                    loaded[kind].append((path, self._load_file(kind, path, id_)))
                except DataError as exc:
                    problems.extend(str(exc).splitlines())
            known[kind] = sorted(first_path)
            for _, path in self._files(kind, suffix=".yml"):
                problems.append(f"{path}: not read: data files must use the .yaml extension")

        species_ids: list[str] | None
        species: dict[str, Species] = {}
        try:
            species = self.species()
            species_ids = sorted(species)
        except DataError as exc:
            problems.extend(str(exc).splitlines())
            species_ids = None

        kinetics: KineticsFile | None = None
        kinetics_path = self._root / _KINETICS_FILE
        if kinetics_path.is_file():
            try:
                kinetics = self.kinetics()
            except DataError as exc:
                problems.extend(str(exc).splitlines())
        problems.extend(self._species_problems(species, kinetics, kinetics_path))

        calibration_ids: set[str] = set()
        records: list[tuple[Path, CalibrationRecord]] = []
        for run_id, path in self._calibration_files():
            try:
                records.append((path, self.calibration_record(run_id)))
            except DataError as exc:
                problems.extend(str(exc).splitlines())
                continue
            calibration_ids.add(run_id)

        for material_id, path in self._tables_files():
            try:
                self.reference_tables(material_id)
            except DataError as exc:
                problems.extend(str(exc).splitlines())
        facts_dir = self._root / _REFERENCE_DIR / _DEVELOPER_FACTS_DIR
        if facts_dir.is_dir():
            # Files no bath id maps to (t-max-rs) are checked too, though D1 does not read them.
            for path in sorted(facts_dir.glob("*.yaml")):
                if path.name.startswith(".") or not path.is_file():
                    continue
                try:
                    _check_reference_yaml(path, _read_yaml(path), "developer", path.stem)
                except DataError as exc:
                    problems.extend(str(exc).splitlines())

        fitted_sources: list[tuple[Path, Any]] = [
            (path, obj) for kind in loaded for path, obj in loaded[kind]
        ]
        if species:
            fitted_sources.append((self._root / _SPECIES_FILE, {"species": list(species.values())}))
        if kinetics is not None:
            fitted_sources.append((kinetics_path, kinetics))
        for path, obj in fitted_sources:
            for loc, quantity in iter_quantities(obj):
                if quantity.source is Source.FITTED and quantity.fit.run_id not in calibration_ids:
                    problems.append(
                        f"{path}: {loc}.fit.run_id: no calibration record "
                        f"{quantity.fit.run_id!r} in {self._root / _CALIBRATION_DIR}"
                    )

        for relpath, path in self._reference_files():
            try:
                curve = self.reference_curve(relpath)
            except DataError as exc:
                problems.extend(str(exc).splitlines())
                continue
            if curve.material not in known["material"]:
                problems.append(
                    f"{path}: material: unknown material {curve.material!r}; "
                    f"available: {_available(known['material'])}"
                )

        for path, record in records:
            for i, pair in enumerate(record.pair_statuses):
                for field, kind, value in (
                    ("material_id", "material", pair.material_id),
                    ("developer_id", "bath", pair.developer_id),
                ):
                    if value not in known[kind]:
                        problems.append(
                            f"{path}: pair_statuses[{i}].{field}: unknown {kind} {value!r}; "
                            f"available: {_available(known[kind])}"
                        )

        for path, proc in loaded["process"]:
            for i, step in enumerate(proc.steps):
                if step.bath is not None and step.bath not in known["bath"]:
                    problems.append(
                        f"{path}: steps[{i}].bath: unknown bath {step.bath!r}; "
                        f"available: {_available(known['bath'])}"
                    )

        if species_ids is not None:
            for path, bath in loaded["bath"]:
                for i, component in enumerate(bath.components):
                    if component.species not in species_ids:
                        problems.append(
                            f"{path}: components[{i}].species: unknown species "
                            f"{component.species!r}; available: {_available(species_ids)}"
                        )

        for path, material in loaded["material"]:
            if material.model == BW_MECHANISTIC:
                problems.extend(f"{path}: {line}" for line in bw_mechanistic_problems(material))
            for loc, curve in _material_curves(material):
                if curve.file is None:
                    continue
                if not (self._root / curve.file).is_file():
                    problems.append(
                        f"{path}: {loc}.file: curve file {curve.file!r} not found under "
                        f"{self._root}"
                    )
                    continue
                try:
                    self.resolve_curve(curve, DEFAULT_GRID)
                except DataError as exc:
                    problems.append(f"{path}: {loc}.file: {exc}")

        return problems

    def _species_problems(
        self, species: dict[str, Species], kinetics: KineticsFile | None, kinetics_path: Path
    ) -> list[str]:
        """D1 cross-reference checks of the species file (spec D1 5.2)."""
        path = self._root / _SPECIES_FILE
        problems = []
        index = {s.id: i for i, s in enumerate(species.values())}
        components = sorted(s.id for s in species.values() if s.kind == "component")
        for s in species.values():
            i = index[s.id]
            if s.id in _BUILTIN_COMPONENTS:
                problems.append(
                    f"{path}: species[{i}].id: {s.id!r} is a built-in pseudo-component and "
                    "must not be defined in the species file (spec D1 2.3.1, 5.2)"
                )
            for target in s.dissolves_to:
                if target not in components:
                    problems.append(
                        f"{path}: species[{i}].dissolves_to.{target}: {target!r} is not a "
                        f"component species; components: {_available(components)}"
                    )
            if s.developing_agent is None:
                continue
            loc = f"{path}: species[{i}].developing_agent.kinetic_class"
            if not kinetics_path.is_file():
                problems.append(f"{loc}: needs {kinetics_path}, which does not exist")
            elif kinetics is not None and s.developing_agent.kinetic_class not in kinetics.agents:
                problems.append(
                    f"{loc}: unknown agent class {s.developing_agent.kinetic_class!r}; "
                    f"available in {kinetics_path}: {_available(sorted(kinetics.agents))}"
                )
        return problems
