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
of its padded input. Following the coding rule of §3.6, the per-block code uses only FFTs
and element-wise arithmetic, and the spectral projections and mixing weights use
``einsum(optimize=False)``: no BLAS call whose result could depend on the thread count or
on memory alignment. The fields are therefore bit-identical for any tiling and worker count.

With several workers the blocks never come back through the main process (spec §4.4): each
worker opens the convolution source (``"r"``) and the output fields (``"r+"``) from
:class:`~darkroom_atlas.fieldstore.FieldDescriptor` objects, writes its inner block and
returns only the block's region, which the main process registers with
:meth:`~darkroom_atlas.state.FrameAssembler.mark`. Outputs that the store keeps in RAM are
written by the workers into a temporary mapped file and copied into RAM at the end.
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
from darkroom_atlas.fieldstore import FieldDescriptor, FieldStore, descriptor_of, open_field
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
    "mixing_weights",
    "uniform_exposure",
]

BLOCK_PX = 512
"""Inner size of the canonical optical blocks, pixels (spec §2.1.7, §3.6)."""

PSF_TAIL = 1e-4
"""Tail mass beyond the scatter PSF's truncation radius ``R_P`` (spec §2.1.6, §12.7)."""

# Bytes of full-spectrum light read per batch of rows when projecting onto the weights.
_PROJECT_BATCH_BYTES = 64 * 2**20

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
    # einsum(optimize=False), not "@": no BLAS reduction (spec §3.6).
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
    The fields' store names carry a fresh ``"expose-<hex>:"`` prefix; a later step releases
    them with ``store.release(store.name_of(field))`` (spec §4.3, §4.5).
    ``geometry`` must describe the whole frame of ``light``. ``workers > 1`` runs the blocks
    in worker processes, which write their blocks straight into memory-mapped outputs
    (spec §4.4); the convolution source reaches them as a mapped file under the store's
    ``root`` (or a temporary directory), removed afterwards.
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
        mix: tuple[np.ndarray, np.ndarray] | None = mixing_weights(light.basis, plus, minus, time_s)
    else:
        source = _project(light, time_s * np.concatenate([plus, minus], axis=1))
        mix = None

    blocks = block_grid(geometry.shape_hw, block_px, halo)
    specs = [FieldSpec(_VOLUME, (n_z,)), FieldSpec(_SURFACE, (n_z,))]
    # Exposure is not a process step: a fresh prefix keeps store names unique when one run
    # exposes several times (spec §4.3 names steps' fields only).
    prefix = f"expose-{uuid.uuid4().hex[:8]}:"
    assembler = FrameAssembler(geometry, specs, store, name_prefix=prefix)
    kernels = {"psf": psf, "halation": halation.kernel, "side": side, "n_z": n_z, "mix": mix}
    if workers == 1 or len(blocks) == 1:
        spectra = (_kernel_spectrum(psf, side), _kernel_spectrum(halation.kernel, side))
        for block in blocks:
            volume, surface = _block_job(block, source=source, spectra=spectra, **kernels)
            assembler.write(_VOLUME, block.inner, volume)
            assembler.write(_SURFACE, block.inner, surface)
        fields = assembler.finish()
    else:
        fields = _run_workers(
            blocks, source, kernels, workers, geometry, specs, store, assembler, prefix
        )
    meta = {
        "halation_fraction": float(halation.fraction),
        "halation_radius_um": float(halation.radius_um),
        "halation_truncated": float(halation.truncated),
        "psf_radius_um": float(psf_radius),
        "block_px": block_px,
        "halo_px": int(halo),
    }
    return fields[_VOLUME], fields[_SURFACE], meta


def mixing_weights(
    basis: np.ndarray, plus: np.ndarray, minus: np.ndarray, time_s: float
) -> tuple[np.ndarray, np.ndarray]:
    """``(t_exp Q W+, t_exp Q W-)``, each ``(R, 2 n_z)``, for low-rank light (spec §2.1.4).

    ``einsum(optimize=False)`` rather than ``@`` (spec §3.6).
    """
    return (
        time_s * np.einsum("rb,bm->rm", basis, plus, optimize=False),
        time_s * np.einsum("rb,bm->rm", basis, minus, optimize=False),
    )


def _project(light: LightField, weights: np.ndarray) -> np.ndarray:
    """``sum_b Phi_b(x) weights_b,m`` over the frame, ``(H, W, M)``, in batches of rows.

    Low-rank light is projected through its basis (``(Q weights)`` first); any other light
    is read with ``window`` a batch of rows at a time. ``einsum(optimize=False)`` only.
    """
    h, w = light.shape_hw
    if isinstance(light, LowRankLight):
        reduced = np.einsum("rb,bm->rm", light.basis, weights, optimize=False)
        return np.einsum("hwr,rm->hwm", light.coefficients, reduced, optimize=False)
    out = np.empty((h, w, weights.shape[1]))
    rows = max(1, _PROJECT_BATCH_BYTES // (w * light.grid.n_bands * 8))
    for r0 in range(0, h, rows):
        band = slice(r0, min(r0 + rows, h))
        flux = light.window(band, slice(0, w)).photon_flux
        out[band] = np.einsum("hwb,bm->hwm", flux, weights, optimize=False)
    return out


def _run_workers(
    blocks: tuple[Block, ...],
    source: np.ndarray,
    kernels: dict[str, Any],
    workers: int,
    geometry: FrameGeometry,
    specs: list[FieldSpec],
    store: FieldStore | None,
    assembler: FrameAssembler,
    prefix: str,
) -> dict[str, np.ndarray]:
    """Run the blocks in worker processes that write their inner blocks themselves (§4.4)."""
    h, w = geometry.shape_hw
    # The assembler (that is, the store) decides which outputs are mapped; fields it keeps
    # in RAM are written by the workers into transient maps and copied in once at the end.
    outputs = assembler.descriptors()
    mapped = set(outputs) == {_VOLUME, _SURFACE}
    scratch, owned = _scratch_dir(store)
    # Transient files: the source (unless it already is a mapped file) and, when the store
    # keeps the outputs in RAM, the outputs the workers write.
    transient = FieldStore(scratch, ram_threshold_bytes=0)
    try:
        try:
            source_desc = descriptor_of(source)
        except ValueError:
            copy = transient.allocate(f"{prefix}source", source.shape)
            copy[...] = source
            del copy
            transient.freeze(f"{prefix}source")
            source_desc = transient.descriptor(f"{prefix}source")
        target = assembler if mapped else FrameAssembler(geometry, specs, transient)
        if not mapped:
            outputs = target.descriptors()
        job = partial(_worker_job, source=source_desc, outputs=outputs, **kernels)
        for inner in run_blocks(job, blocks, workers=workers):
            target.mark(inner)
        fields = target.finish()
        if not mapped:
            for key, value in fields.items():
                # One copy, from the transient map into the RAM output.
                assembler.write(key, (slice(0, h), slice(0, w)), value)
            fields = assembler.finish()
        return fields
    finally:
        transient.close()
        if owned:
            shutil.rmtree(scratch, ignore_errors=True)


def _worker_job(
    block: Block,
    *,
    source: FieldDescriptor,
    outputs: dict[str, FieldDescriptor],
    **kernels: Any,
) -> tuple[slice, slice]:
    """Worker side: convolve one block, write it into the mapped outputs, return its region."""
    volume, surface = _block_job(block, source=open_field(source, "r"), **kernels)
    for key, values in ((_VOLUME, volume), (_SURFACE, surface)):
        out = open_field(outputs[key], "r+")
        out[block.inner] = values
        out.flush()
        del out
    return block.inner


def _block_job(
    block: Block,
    *,
    source: np.ndarray,
    psf: np.ndarray,
    halation: np.ndarray,
    side: int,
    n_z: int,
    mix: tuple[np.ndarray, np.ndarray] | None,
    spectra: tuple[np.ndarray, np.ndarray] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Convolve one padded block and return its inner ``(e_V, e_S)``.

    Only FFTs and element-wise arithmetic: the result does not depend on the process, the
    thread count or the memory alignment of the inputs (spec §3.6). ``spectra`` are the
    kernels' ``rfft2`` (computed here when None).
    """
    padded = padded_window(source, block)
    shape = (side, side)
    spectrum = fft.rfft2(padded, s=shape, axes=(0, 1))
    rows, cols = block.inner_in_padded
    if spectra is None:
        spectra = (_kernel_spectrum(psf, side), _kernel_spectrum(halation, side))
    psf_spectrum, halation_spectrum = spectra

    def convolve(kernel_spectrum: np.ndarray, channels: slice) -> np.ndarray:
        product = spectrum[:, :, channels] * kernel_spectrum[:, :, None]
        return fft.irfft2(product, s=shape, axes=(0, 1))[rows, cols]

    if mix is None:
        total = convolve(psf_spectrum, slice(0, 2 * n_z)) + convolve(
            halation_spectrum, slice(2 * n_z, 4 * n_z)
        )
    else:
        mix_plus, mix_minus = mix
        scattered = convolve(psf_spectrum, slice(None))
        returned = convolve(halation_spectrum, slice(None))
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
    """Directory for transient worker files: the store's root, else a new temporary one."""
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
