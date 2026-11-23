"""Weighted SVD of the calibration Jacobian (spec D1 7.5)."""

import numpy as np
import pytest

from darkroom_atlas.calibrate.identifiability import weighted_svd


def _weighted(j: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    return j / sigma[:, None]


def test_weighted_svd_matches_numpy_on_random_matrix() -> None:
    rng = np.random.default_rng(0)
    j = rng.normal(size=(40, 7))
    sigma = rng.uniform(0.01, 0.5, size=40)
    s, v = weighted_svd(j, sigma)
    a = _weighted(j, sigma)
    s_np = np.linalg.svd(a, compute_uv=False)
    np.testing.assert_allclose(s, s_np, rtol=1e-13, atol=0)
    assert v.shape == (7, 7)
    np.testing.assert_allclose(v.T @ v, np.eye(7), atol=1e-13)
    # Columns are the right singular vectors: ||A v_k|| = s_k and A^T A v_k = s_k^2 v_k.
    np.testing.assert_allclose(np.linalg.norm(a @ v, axis=0), s, rtol=1e-12)
    np.testing.assert_allclose(a.T @ a @ v, v * s**2, rtol=0, atol=1e-10 * s[0] ** 2)
    # Same subspaces as numpy's vectors, up to sign.
    _, _, vt_np = np.linalg.svd(a, full_matrices=False)
    np.testing.assert_allclose(np.abs(np.sum(v * vt_np.T, axis=0)), 1.0, atol=1e-12)


def test_weighted_svd_rank_deficient() -> None:
    rng = np.random.default_rng(1)
    base = rng.normal(size=(30, 3))
    # Column 3 = column 0 + 2 column 1; column 4 is zero (a parameter only the prior constrains).
    j = np.column_stack([base, base[:, 0] + 2 * base[:, 1], np.zeros(30)])
    sigma = rng.uniform(0.05, 0.2, size=30)
    s, v = weighted_svd(j, sigma)
    a = _weighted(j, sigma)
    np.testing.assert_allclose(s, np.linalg.svd(a, compute_uv=False), rtol=1e-12, atol=1e-12 * s[0])
    assert np.all(s[:3] > 1e-3 * s[0])
    assert np.all(s[3:] < 1e-12 * s[0])
    np.testing.assert_allclose(v.T @ v, np.eye(5), atol=1e-12)
    np.testing.assert_allclose(a @ v[:, 3:], 0.0, atol=1e-12 * s[0])
    # The null space holds the zero column and the combination (1, 2, 0, -1, 0)/sqrt(6).
    null = v[:, 3:]
    for direction in (np.array([0, 0, 0, 0, 1.0]), np.array([1, 2, 0, -1, 0]) / np.sqrt(6)):
        np.testing.assert_allclose(np.linalg.norm(null.T @ direction), 1.0, atol=1e-12)


def test_weighted_svd_weights() -> None:
    j = np.array([[1.0, 0.0], [0.0, 1.0]])
    s, _ = weighted_svd(j, np.array([0.5, 2.0]))
    np.testing.assert_allclose(s, [2.0, 0.5], rtol=1e-15)


def test_weighted_svd_wide_matrix() -> None:
    rng = np.random.default_rng(2)
    j = rng.normal(size=(3, 6))
    s, v = weighted_svd(j, np.ones(3))
    assert s.shape == (3,)
    assert v.shape == (6, 3)


@pytest.mark.parametrize(
    ("jacobian", "sigma"),
    [
        (np.ones(3), np.ones(3)),
        (np.ones((2, 2, 2)), np.ones(2)),
        (np.ones((0, 2)), np.ones(0)),
        (np.ones((3, 2)), np.ones(2)),
        (np.ones((3, 2)), np.ones((3, 1))),
        (np.ones((3, 2)), np.array([1.0, 0.0, 1.0])),
        (np.ones((3, 2)), np.array([1.0, -1.0, 1.0])),
        (np.ones((3, 2)), np.array([1.0, np.nan, 1.0])),
        (np.ones((3, 2)), np.array([1.0, np.inf, 1.0])),
        (np.array([[1.0, np.nan], [0.0, 1.0]]), np.ones(2)),
    ],
)
def test_weighted_svd_bad_input(jacobian: np.ndarray, sigma: np.ndarray) -> None:
    with pytest.raises(ValueError):
        weighted_svd(jacobian, sigma)
