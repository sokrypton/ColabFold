"""Route the upstream sampler's denoised predictions through RGI's minimizer.

alphafold3-colabfold 3.1.18 and later expose this point as
``diffusion_head.DENOISED_HOOK``. For 3.1.11, the pinned release before it,
checked hook points are added to the sampler's source in memory instead. Either
way the original sampler owns all schedules, random keys, Euler/Heun steps and
alignment. No installed file is modified, and vanilla keeps the original function.
"""

from __future__ import annotations

import hashlib
import inspect
from contextlib import contextmanager
from functools import lru_cache


def has_denoised_hook(diffusion_head):
    """True when the installed release exposes the guidance point (3.1.18+)."""
    return hasattr(diffusion_head, "DENOISED_HOOK")


@contextmanager
def guided(diffusion_head, minimize):
    """Within the block, ``minimize(coords, sigma, step)`` sees every denoised prediction.

    Traced code reads it, so call (and jit) the model inside the block.
    """
    if has_denoised_hook(diffusion_head):
        name, value = "DENOISED_HOOK", minimize
    else:
        name, value = "sample", make_guided_sampler(diffusion_head.sample, minimize)
    original = getattr(diffusion_head, name)
    setattr(diffusion_head, name, value)
    try:
        yield
    finally:
        setattr(diffusion_head, name, original)


@lru_cache(maxsize=4)
def _instrument(original):
    source = inspect.getsource(original)
    digest = hashlib.sha256(source.encode()).hexdigest()
    if digest != EXPECTED_SAMPLE_SHA256:
        raise RuntimeError(
            "The installed diffusion sampler differs from the tested release. "
            "Re-run the install cell with alphafold3-colabfold==3.1.11. "
            f"Sampler SHA256: {digest}"
        )
    changes = (
        (
            "noise_level, noise_level_prev = step_noise",
            "noise_level, noise_level_prev, rgi_step = step_noise",
        ),
        (
            "positions_denoised = denoising_step(positions_noisy, t_hat)",
            (
                "positions_denoised = denoising_step(positions_noisy, t_hat)\n"
                "    positions_denoised = _rgi_minimize("
                "positions_denoised, noise_level_prev, rgi_step)"
            ),
        ),
        (
            "denoising_step(positions_out, noise_level)",
            (
                "_rgi_minimize(denoising_step(positions_out, noise_level), "
                "noise_level_prev, rgi_step)"
            ),
        ),
        (
            "(noise_levels[1:], noise_levels[:-1]),",
            (
                "(noise_levels[1:], noise_levels[:-1], "
                "jnp.arange(noise_levels.shape[0] - 1)),"
            ),
        ),
    )
    for old, new in changes:
        if source.count(old) != 1:
            raise RuntimeError(f"Cannot locate exactly one RGI sampler hook: {old}")
        source = source.replace(old, new)
    return compile(source, inspect.getsourcefile(original), "exec")


def make_guided_sampler(original, minimize):
    namespace = dict(original.__globals__, _rgi_minimize=minimize)
    exec(_instrument(original), namespace)  # noqa: S102 - SHA256-checked upstream source.
    return namespace[original.__name__]


def restraint_layout(example, model):
    """Match the restraint atom mask to the predictor's sampling coordinates."""
    import numpy as np

    bioemu = model == "bioemu2"
    if bioemu and any(key.startswith("bioemu2steer/preset_") for key in example):
        raise ValueError("BioEmu2 RGI requires BIOEMU2_SAMPLER=edm.")
    key = "bioemu2/af3_index" if bioemu else "structbook/residue_atom_gather"
    gather = example.get(key)
    if gather is None:
        if bioemu:
            raise ValueError(
                "BioEmu2 RGI needs its atom map from alphafold3-colabfold 3.1.18 or later."
            )
        return example, None
    gather = np.asarray(gather)
    mask = np.asarray(example["ref_mask"], bool)
    if gather.shape != mask.shape or not np.issubdtype(gather.dtype, np.integer):
        raise ValueError(
            f"{model} atom map must be an integer array matching ref_mask."
        )
    valid = gather[gather >= 0]
    if len(np.unique(valid)) != len(valid):
        raise ValueError(f"{model} atom map contains duplicate atoms.")
    if bioemu:
        atom_mask = np.asarray(example["bioemu2/atom_mask"], bool)
        if np.any(valid >= len(atom_mask)) or not np.all(atom_mask[valid]):
            raise ValueError("BioEmu2 atom map points outside its predicted atoms.")
        # BioEmu2 omits terminal OXT/OP3 and leaving atoms retained in ref_mask.
        example = dict(example, ref_mask=mask & (gather >= 0))
    elif np.any(mask & (gather < 0)):
        raise ValueError("OpenDDE structural atom map does not cover the RGI atoms.")
    return example, gather


def coordinate_minimizer(minimize, gather=None):
    """Bridge AF3, OpenDDE, or BioEmu2 sampling coordinates to the engine."""
    import jax.numpy as jnp

    def apply(coords, sigma, step):
        flat = coords.reshape(-1, 3)
        if gather is None:
            return minimize(flat, sigma, step).reshape(coords.shape)
        # RGI retains residue numbering across reordered sampling atom layouts.
        # Invalid padded entries scatter out of bounds with mode='drop'.
        indices = gather.reshape(-1)
        valid = indices >= 0
        residue = jnp.where(valid[:, None], flat[jnp.maximum(indices, 0)], 0)
        corrected = minimize(residue, sigma, step)
        destinations = jnp.where(valid, indices, flat.shape[0])
        return flat.at[destinations].set(corrected, mode="drop").reshape(coords.shape)

    return apply


EXPECTED_SAMPLE_SHA256 = (
    "d238e26b0b7f4f86efaa121c1ea2b2279b673c97dc9b7ed332f4ceffbe43ef1c"
)
