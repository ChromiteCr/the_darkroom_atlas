"""Shared fixtures for the validation tests: the digitised Kodak reference curves."""

from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest

# The one place that says where the reference CSVs live (spec 0.1, 5.9).
_REFERENCE_KODAK = Path(__file__).resolve().parents[2] / "data" / "reference" / "kodak"


@pytest.fixture(scope="session")
def reference_kodak_dir() -> Path:
    """Directory holding ``<film>/*.csv`` of the digitised Kodak data."""
    if not _REFERENCE_KODAK.is_dir():
        pytest.skip(f"reference data not found at {_REFERENCE_KODAK}")
    return _REFERENCE_KODAK


def _read(path: Path) -> tuple[str, np.ndarray, np.ndarray]:
    text = path.read_text(encoding="utf-8")
    rows = [line for line in text.splitlines() if line.strip() and not line.startswith("#")]
    data = np.array([[float(v) for v in row.split(",")] for row in rows], dtype=np.float64)
    return text, data[:, 0].copy(), data[:, 1].copy()


@pytest.fixture(scope="session")
def read_reference(
    reference_kodak_dir: Path,
) -> Callable[[str], tuple[str, np.ndarray, np.ndarray]]:
    """``read_reference("t-max-100/char-d76-06min.csv") -> (text, x, y)``."""

    def read(relpath: str) -> tuple[str, np.ndarray, np.ndarray]:
        return _read(reference_kodak_dir / relpath)

    return read
