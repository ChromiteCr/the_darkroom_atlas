"""Identifiability of the calibration parameters (spec D1 7.5).

:func:`weighted_svd` gives the singular values and right singular vectors of the weighted
Jacobian ``W^{1/2} J`` with ``W = diag(1 / sigma^2)``. Spec D1 7.5 runs it on the data residuals
only, in transform space, block by block.
"""

import numpy as np

__all__ = ["weighted_svd"]


def weighted_svd(jacobian: np.ndarray, sigma: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Singular values and right singular vectors of ``jacobian / sigma[:, None]``.

    ``jacobian`` has shape ``(n_residuals, n_parameters)``, holding the derivatives of the
    *unweighted* residuals with respect to the transformed parameters. ``sigma`` has shape
    ``(n_residuals,)``, holding each residual's standard deviation (finite and positive).
    The result is ``(s, V)``. ``s`` has ``n_parameters`` values in decreasing order. The columns
    of ``V`` (shape ``(n_parameters, n_parameters)``, orthonormal) are the matching right
    singular vectors. With fewer residuals than parameters, ``s`` is padded with zeros and ``V``
    holds the null directions too, so that every unconstrained direction of a block appears with
    singular value 0 for the spec D1 7.5 rule. The sign of each vector is arbitrary, as in
    :func:`numpy.linalg.svd`.
    """
    j = np.asarray(jacobian, dtype=np.float64)
    s = np.asarray(sigma, dtype=np.float64)
    if j.ndim != 2 or j.shape[0] == 0 or j.shape[1] == 0:
        raise ValueError(f"jacobian must be a non-empty 2-D array, got shape {j.shape}")
    if s.shape != (j.shape[0],):
        raise ValueError(f"sigma must have shape ({j.shape[0]},), got {s.shape}")
    if not np.all(np.isfinite(j)):
        raise ValueError("jacobian must be finite")
    if not (np.all(np.isfinite(s)) and np.all(s > 0.0)):
        raise ValueError("sigma must be finite and positive")
    n_residuals, n_parameters = j.shape
    # A full U is only needed (and only cheap) when it is n_residuals x n_residuals with
    # n_residuals < n_parameters; then full_matrices=True gives the null rows of V^T.
    wide = n_residuals < n_parameters
    _, values, vt = np.linalg.svd(j / s[:, None], full_matrices=wide)
    if wide:
        values = np.concatenate([values, np.zeros(n_parameters - n_residuals)])
    return values, vt.T.copy()
