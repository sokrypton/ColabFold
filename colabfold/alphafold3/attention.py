"""
Run alphafold3's attention on ColabFold's fused kernels instead of tokamax.

alphafold3 calls one function for every attention
(``alphafold3.model.components.attention.dot_product_attention``) with q/k/v as
(..., seq, heads, dim); ColabFold's kernels are heads-major [batch, heads, seq, dim].
:func:`install` swaps the two, so the AF3 graph runs on tri_flash on Ampere and on the
prebuilt sm_70/sm_75 kernels below it.
"""
import logging

import numpy as np

logger = logging.getLogger(__name__)

NEG = -1e9


def _to_heads_major(x):
    """(..., seq, heads, dim) -> [batch, heads, seq, dim], with the batch flattened."""
    import jax.numpy as jnp

    lead = tuple(x.shape[:-3])
    seq, heads, dim = x.shape[-3:]
    x = jnp.reshape(x, (int(np.prod(lead)) if lead else 1, seq, heads, dim))
    return jnp.swapaxes(x, 1, 2), lead


def _mask_bias(mask, batch: int, keys: int, dtype):
    """alphafold3's boolean mask -> ColabFold's additive [batch, 1, 1, S_kv]."""
    import jax.numpy as jnp

    if mask is None:
        return jnp.zeros((batch, 1, 1, keys), dtype)
    mask = jnp.asarray(mask)
    flat = jnp.reshape(mask, (-1, mask.shape[-1]))
    if flat.shape[0] not in (1, batch):
        flat = flat[:1]
    flat = jnp.broadcast_to(flat, (batch, mask.shape[-1]))
    if flat.shape[-1] != keys:
        flat = jnp.broadcast_to(flat[..., :1], (batch, keys))
    return jnp.where(flat, 0.0, NEG).astype(dtype)[:, None, None, :]


def _nonbatched_bias(bias, heads: int, queries: int, keys: int, dtype):
    """alphafold3's (..., heads, q, k) bias -> ColabFold's shared [heads, q, k]."""
    import jax.numpy as jnp

    if bias is None:
        return None
    bias = jnp.asarray(bias)
    if bias.ndim > 3:
        lead = int(np.prod(bias.shape[:-3]))
        if lead != 1:
            return None  # per-row bias: the kernel only takes a shared one
        bias = jnp.reshape(bias, bias.shape[-3:])
    bias = jnp.broadcast_to(bias, (heads, queries, keys))
    return bias.astype(dtype)


def colabfold_attention(q, k, v, *, mask=None, bias=None, implementation=None, scale=None,
                        kernel=None):
    """alphafold3's attention signature, computed by ColabFold's kernel."""
    import jax.numpy as jnp

    if kernel is None:
        from alphafold.model.tri_flash import pallas_attention as kernel

    qh, lead = _to_heads_major(q)
    kh, _ = _to_heads_major(k)
    vh, _ = _to_heads_major(v)
    batch, heads, queries, dim = qh.shape
    keys = kh.shape[2]

    nonbatched = _nonbatched_bias(bias, heads, queries, keys, qh.dtype)
    if nonbatched is None and bias is not None:
        return None  # caller falls back
    if scale is None:
        scale = dim ** -0.5

    out = kernel(qh, kh, vh, _mask_bias(mask, batch, keys, qh.dtype), nonbatched, scale)
    out = jnp.swapaxes(out, 1, 2)
    return jnp.reshape(out, lead + (queries, heads, dim))


def install(fallback=True) -> bool:
    """Point alphafold3's attention at ColabFold's kernels. True if it took."""
    try:
        from alphafold3.model.components import attention as af3_attention
    except ModuleNotFoundError:
        return False

    original = getattr(af3_attention, "_colabfold_original", af3_attention.dot_product_attention)

    def dot_product_attention(q, k, v, *, mask=None, bias=None, implementation=None, scale=None):
        try:
            out = colabfold_attention(q, k, v, mask=mask, bias=bias, scale=scale)
            if out is not None:
                return out
        except Exception as e:
            if not fallback:
                raise
            logger.warning(f"colabfold kernels declined this attention, using alphafold3's: {e}")
        return original(q, k, v, mask=mask, bias=bias, implementation=implementation, scale=scale)

    af3_attention._colabfold_original = original
    af3_attention.dot_product_attention = dot_product_attention
    logger.info("alphafold3 attention is running on ColabFold's fused kernels")
    return True


def uninstall() -> None:
    from alphafold3.model.components import attention as af3_attention

    original = getattr(af3_attention, "_colabfold_original", None)
    if original is not None:
        af3_attention.dot_product_attention = original
        del af3_attention._colabfold_original
