"""Add checked hook points to the pinned upstream sampler in memory.

The original sampler owns all schedules, random keys, Euler/Heun steps and
alignment. No installed file is modified, and vanilla keeps the original function.
"""

from __future__ import annotations

import hashlib
import inspect
from functools import lru_cache


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


def coordinate_minimizer(minimize, gather=None):
    """Bridge dense AF3 or OpenDDE structural-token coordinates to the engine."""
    import jax.numpy as jnp

    def apply(coords, sigma, step):
        flat = coords.reshape(-1, 3)
        if gather is None:
            return minimize(flat, sigma, step).reshape(coords.shape)
        # OpenDDE diffuses structural subtokens; RGI retains per-residue numbering.
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
