"""
Run alphafold3's attention on ColabFold's fused kernels instead of tokamax.

alphafold3 calls one function for every attention
(``alphafold3.model.components.attention.dot_product_attention``) with q/k/v as
(..., seq, heads, dim); ColabFold's kernels read either that layout or heads-major [batch, heads, seq, dim].
:func:`colabfold_attention` adapts the masks and the bias; the routing itself lives in
:mod:`colabfold.alphafold3.tokamax_shim`.
"""
import logging

import numpy as np

logger = logging.getLogger(__name__)

NEG = -1e9


def _flatten_batch(x):
    """(..., seq, heads, dim) -> [batch, seq, heads, dim], with the batch flattened."""
    import jax.numpy as jnp

    lead = tuple(x.shape[:-3])
    return jnp.reshape(x, (int(np.prod(lead)) if lead else 1,) + x.shape[-3:]), lead


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
                        kernel=None, seq_major=False):
    """alphafold3's attention signature, computed by ColabFold's kernel.

    seq_major says the kernel reads (batch, seq, heads, dim) as alphafold3 hands it
    over, so q/k/v and the output need no transpose.
    """
    import jax.numpy as jnp

    if kernel is None:
        from colabfold_kernels.tri_flash import pallas_attention as kernel

        seq_major = False

    qh, lead = _flatten_batch(q)
    kh, _ = _flatten_batch(k)
    vh, _ = _flatten_batch(v)
    batch, queries, heads, dim = qh.shape
    keys = kh.shape[1]
    if not seq_major:
        qh, kh, vh = (jnp.swapaxes(x, 1, 2) for x in (qh, kh, vh))

    nonbatched = _nonbatched_bias(bias, heads, queries, keys, qh.dtype)
    if nonbatched is None and bias is not None:
        return None  # caller falls back
    if scale is None:
        scale = dim ** -0.5

    out = kernel(qh, kh, vh, _mask_bias(mask, batch, keys, qh.dtype), nonbatched, scale)
    if not seq_major:
        out = jnp.swapaxes(out, 1, 2)
    return jnp.reshape(out, lead + (queries, heads, dim))
