"""Exposure fields: light on the film to the absorption descriptors of each depth cell (§2.1.4).

For each depth cell ``z`` the two descriptors are

    e_X,z = P * [sum_b W+_X,z,b N_b] + H * [sum_b W-_X,z,b N_b],     X in {V, S}
    W+_X,z,b = w^X_b g+_b,z,   W-_X,z,b = w^X_b g-_b,z tau_b,   w^V = alpha, w^S = s

with ``N_b = Phi_b t_exp``, ``P`` the scatter PSF and ``H`` the halation kernel (both
spectrally neutral, so the spectral projection comes before the convolution). For a
:class:`~darkroom_atlas.light.LowRankLight` with ``R < 2 n_z`` the ``R`` coefficient channels
are convolved instead and combined afterwards with ``t_exp Q W`` (``2R`` convolutions); for
any other light the ``2 n_z`` projections of each path are convolved (``4 n_z``
convolutions).

The convolutions run by FFT on the canonical block grid of spec §3.6: 512-pixel inner
blocks anchored at the frame origin, a halo of ``ceil(max(R_H, R_P)/p)`` pixels, zero fill
outside the frame and one FFT size for every block. Each block's result is a pure function
of its padded input, and the per-block code uses only FFTs and element-wise arithmetic (no
BLAS reductions), so the fields are bit-identical for any worker count.
"""

import math
import shutil
import tempfile
import uuid
from functools import lru_cache, partial
from numbers import Integral, Real
from pathlib import Path
from typing import Any

import numpy as np
from scipy import fft

from darkroom_atlas.blocks import Block, block_grid, padded_window, run_blocks
from darkroom_atlas.fieldstore import FieldStore
from darkroom_atlas.light import LightField, LowRankLight
from darkroom_atlas.optics.halation import halation_fraction, halation_kernel
from darkroom_atlas.optics.kubelka import DepthResponse
from darkroom_atlas.optics.psf import scatter_kernel, tail_radius_ratio
from darkroom_atlas.optics.stack import OpticalStack
from darkroom_atlas.state import FieldSpec, FrameAssembler, FrameGeometry

__all__ = [
    "BLOCK_PX",
    "PSF_TAIL",
    "depth_weights",
    "exposure_fields",
    "uniform_exposure",
]

BLOCK_PX = 512
"""Inner size of the canonical optical blocks, pixels (spec §2.1.7, §3.6)."""

PSF_TAIL = 1e-4
"""Tail mass beyond the scatter PSF's truncation radius ``R_P`` (spec §2.1.6, §12.7)."""

# Blocks handed to one worker pool at a time, per worker (bounds the results held at once).
_BLOCKS_PER_WORKER_BATCH = 8

_VOLUME = "emulsion/exposure_volume"
_SURFACE = "emulsion/exposure_surface"


def depth_weights(stack: OpticalStack, response: DepthResponse) -> tuple[np.ndarray, np.ndarray]:
    """Spectral weights ``(W+, W-)`` of spec §2.1.4, each ``(B, 2 n_z)``.

    Column ``z`` (0-based) is ``X = V`` in cell ``z``, column ``n_z + z`` is ``X = S``.
    """
    _check_stack_response(stack, response)
    n_z = response.n_z
    w_x = np.stack([stack.alpha_per_um, stack.dye_cross_section])  # (2, B)
    plus = w_x[:, None, :] * response.plus[None, :, :]  # (2, n_z, B)
    minus = w_x[:, None, :] * (response.minus * response.into_base)[None, :, :]
    return (
        np.ascontiguousarray(plus.reshape(2 * n_z, -1).T),
        np.ascontiguousarray(minus.reshape(2 * n_z, -1).T),
    )


_halation_fraction = lru_cache(maxsize=64)(halation_fraction)


def uniform_exposure(
    photons_per_band: np.ndarray, stack: OpticalStack, response: DepthResponse
) -> tuple[np.ndarray, np.ndarray]:
    """Descriptors ``(e_V, e_S)`` of a laterally uniform exposure, each ``(n_z,)``.

    ``photons_per_band`` is ``N_b`` (photons um^-2 per band), shape ``(B,)``; leading
    dimensions ``(..., B)`` give ``(..., n_z)``. ``P`` sums to 1 and the halation returns the
    untruncated fraction ``f_H`` of the light entering the base (spec §6.1), so
    ``e_X,z = sum_b N_b (W+_X,z,b + f_H W-_X,z,b)``. This is the one-column path of the
    sensitometer; the interior of a uniform frame in :func:`exposure_fields` differs only by
    the truncated halation tail (``meta["halation_truncated"]``).
    """
    _check_stack_response(stack, response)
    return _uniform(photons_per_band, stack, response, _halation_fraction(stack.base))


def _uniform(
    photons_per_band: np.ndarray,
    stack: OpticalStack,
    response: DepthResponse,
    halation_weight: float,
) -> tuple[np.ndarray, np.ndarray]:
    n_bands = stack.grid.n_bands
    try:
        photons = np.asarray(photons_per_band, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"photons_per_band must be a real numeric array: {exc}") from exc
    if photons.ndim < 1 or photons.shape[-1] != n_bands:
        raise ValueError(f"photons_per_band must have shape (..., {n_bands}), got {photons.shape}")
    if not np.all(np.isfinite(photons)) or np.any(photons < 0.0):
        raise ValueError("photons_per_band must be finite and non-negative")
    plus, minus = depth_weights(stack, response)
    weights = plus + halation_weight * minus
    out = np.einsum("...b,bm->...m", photons, weights, optimize=False)
    n_z = response.n_z
    return out[..., :n_z], out[..., n_z:]


def exposure_fields(
    light: LightField,
    time_s: float,
    stack: OpticalStack,
    response: DepthResponse,
    geometry: FrameGeometry,
    *,
    block_px: int = BLOCK_PX,
    workers: int = 1,
    store: FieldStore | None = None,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Exposure descriptors of a whole frame (spec §2.1.4-§2.1.7).

    Returns ``(e_V, e_S, meta)``: ``e_V`` (um^-3) and ``e_S`` (um^-2) have shape
    ``(H, W, n_z)`` and are read-only (memory-mapped when ``store`` maps them); ``meta``
    holds ``halation_fraction`` (``f_H``), ``halation_radius_um`` (``R_H``),
    ``halation_truncated``, ``psf_radius_um`` (``R_P``), ``block_px`` and ``halo_px``.
    ``geometry`` must describe the whole frame of ``light``. ``workers > 1`` runs the blocks
    in worker processes (the convolution source is passed through a scratch file).
    """
    if not isinstance(light, LightField):
        raise TypeError(f"light must be a LightField, got {type(light).__name__}")
    _check_stack_response(stack, response)
    if not isinstance(geometry, FrameGeometry):
        raise TypeError(f"geometry must be a FrameGeometry, got {type(geometry).__name__}")
    if light.grid != stack.grid:
        raise ValueError(f"light is on grid {light.grid}, the optical stack on {stack.grid}")
    if geometry.shape_hw != tuple(light.shape_hw):
        raise ValueError(
            f"geometry has shape {geometry.shape_hw}, the light {tuple(light.shape_hw)}"
        )
    if geometry.origin_px != (0, 0) or geometry.frame_shape_hw != geometry.shape_hw:
        raise ValueError("exposure_fields needs the geometry of a whole frame, not a window")
    if (
        isinstance(time_s, bool)
        or not isinstance(time_s, Real)
        or not math.isfinite(time_s)
        or time_s <= 0.0
    ):
        raise ValueError(f"time_s must be a positive, finite number, got {time_s!r}")
    block_px = _positive_int("block_px", block_px)
    workers = _positive_int("workers", workers)
    if store is not None and not isinstance(store, FieldStore):
        raise TypeError(f"store must be a FieldStore or None, got {type(store).__name__}")

    n_z = response.n_z
    pitch = geometry.pixel_pitch_um
    psf = scatter_kernel(pitch, stack.psf_weight, stack.psf_length_um, tail=PSF_TAIL)
    halation = halation_kernel(pitch, stack.base)
    psf_radius = tail_radius_ratio(PSF_TAIL) * stack.psf_length_um
    # Halo of spec §2.1.7; never smaller than either kernel's half width.
    halo = max(
        math.ceil(max(halation.radius_um, psf_radius) / pitch),
        psf.shape[0] // 2,
        halation.kernel.shape[0] // 2,
    )
    side = fft.next_fast_len(block_px + 2 * halo, real=True)

    plus, minus = depth_weights(stack, response)
    time_s = float(time_s)
    if isinstance(light, LowRankLight) and light.rank < 2 * n_z:
        source = light.coefficients
        mix: tuple[np.ndarray, np.ndarray] | None = (
            time_s * (light.basis @ plus),
            time_s * (light.basis @ minus),
        )
    else:
        source = light.project(time_s * np.concatenate([plus, minus], axis=1))
        mix = None

    blocks = block_grid(geometry.shape_hw, block_px, halo)
    assembler = FrameAssembler(
        geometry, [FieldSpec(_VOLUME, (n_z,)), FieldSpec(_SURFACE, (n_z,))], store
    )
    job = partial(
        _block_job,
        psf=psf,
        halation=halation.kernel,
        side=side,
        n_z=n_z,
        mix=mix,
    )
    if workers == 1 or len(blocks) == 1:
        for block in blocks:
            _write(assembler, block, job(block, source=source))
    else:
        scratch, owned = _scratch_dir(store)
        try:
            path = scratch / f"exposure-source-{uuid.uuid4().hex[:12]}.npy"
            np.save(path, source)
            remote = partial(job, source=str(path))
            batch = _BLOCKS_PER_WORKER_BATCH * workers
            for start in range(0, len(blocks), batch):
                chunk = blocks[start : start + batch]
                for block, result in zip(
                    chunk, run_blocks(remote, chunk, workers=workers), strict=True
                ):
                    _write(assembler, block, result)
            path.unlink(missing_ok=True)
        finally:
            if owned:
                shutil.rmtree(scratch, ignore_errors=True)
    fields = assembler.finish()
    meta = {
        "halation_fraction": float(halation.fraction),
        "halation_radius_um": float(halation.radius_um),
        "halation_truncated": float(halation.truncated),
        "psf_radius_um": float(psf_radius),
        "block_px": block_px,
        "halo_px": int(halo),
    }
    return fields[_VOLUME], fields[_SURFACE], meta


def _write(assembler: FrameAssembler, block: Block, result: tuple[np.ndarray, np.ndarray]) -> None:
    volume, surface = result
    assembler.write(_VOLUME, block.inner, volume)
    assembler.write(_SURFACE, block.inner, surface)


def _block_job(
    block: Block,
    *,
    source: np.ndarray | str,
    psf: np.ndarray,
    halation: np.ndarray,
    side: int,
    n_z: int,
    mix: tuple[np.ndarray, np.ndarray] | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Convolve one padded block and return its inner ``(e_V, e_S)``.

    Only FFTs and element-wise arithmetic: the result does not depend on the process, the
    thread count or the memory alignment of the inputs (spec §3.6).
    """
    array = np.load(source, mmap_mode="r") if isinstance(source, str) else source
    padded = padded_window(array, block)
    shape = (side, side)
    spectrum = fft.rfft2(padded, s=shape, axes=(0, 1))
    rows, cols = block.inner_in_padded

    def convolve(kernel: np.ndarray, channels: slice) -> np.ndarray:
        product = spectrum[:, :, channels] * _kernel_spectrum(kernel, side)[:, :, None]
        return fft.irfft2(product, s=shape, axes=(0, 1))[rows, cols]

    if mix is None:
        total = convolve(psf, slice(0, 2 * n_z)) + convolve(halation, slice(2 * n_z, 4 * n_z))
    else:
        mix_plus, mix_minus = mix
        scattered = convolve(psf, slice(None))
        returned = convolve(halation, slice(None))
        total = scattered[:, :, 0, None] * mix_plus[0]
        total += returned[:, :, 0, None] * mix_minus[0]
        for r in range(1, mix_plus.shape[0]):
            total += scattered[:, :, r, None] * mix_plus[r]
            total += returned[:, :, r, None] * mix_minus[r]
    return np.ascontiguousarray(total[:, :, :n_z]), np.ascontiguousarray(total[:, :, n_z:])


def _kernel_spectrum(kernel: np.ndarray, side: int) -> np.ndarray:
    """``rfft2`` of an odd square kernel placed with its centre at ``[0, 0]`` (circular)."""
    half = kernel.shape[0] // 2
    if kernel.shape[0] > side:
        raise ValueError(f"kernel of size {kernel.shape[0]} exceeds the FFT size {side}")
    placed = np.zeros((side, side))
    placed[: kernel.shape[0], : kernel.shape[1]] = kernel
    placed = np.roll(placed, (-half, -half), axis=(0, 1))
    return fft.rfft2(placed)


def _scratch_dir(store: FieldStore | None) -> tuple[Path, bool]:
    """Directory for the worker source file: the store's root, else a new temporary one."""
    if store is not None and store.root is not None:
        store.root.mkdir(parents=True, exist_ok=True)
        return store.root, False
    return Path(tempfile.mkdtemp(prefix="atlas-expose-")), True


def _check_stack_response(stack: OpticalStack, response: DepthResponse) -> None:
    if not isinstance(stack, OpticalStack):
        raise TypeError(f"stack must be an OpticalStack, got {type(stack).__name__}")
    if not isinstance(response, DepthResponse):
        raise TypeError(f"response must be a DepthResponse, got {type(response).__name__}")
    if response.plus.shape[1] != stack.grid.n_bands:
        raise ValueError(
            f"response has {response.plus.shape[1]} bands, the stack {stack.grid.n_bands}"
        )


def _positive_int(name: str, value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
        raise ValueError(f"{name} must be a positive integer, got {value!r}")
    return int(value)
