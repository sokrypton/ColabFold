"""
Stand in for tokamax so alphafold3 runs on ColabFold's kernels.

Routing only: each name alphafold3 imports is mapped onto a ColabFold kernel
through ``colabfold_kernels.fused_ops``, then onto tokamax where it runs, then
onto jax's own, for the shapes each one declines. No maths lives here.
"""
import functools
import importlib.util
import logging
import sys

logger = logging.getLogger(__name__)

DotProductAttentionImplementation = str  # alphafold3 only uses it as an annotation

_REAL = None  # the real tokamax, kept because importing it would find us instead

# the modules that call tokamax, and so bind it at import; model_config only annotates with it
_CALLERS = ("alphafold3.model.network.modules", "alphafold3.model.components.attention")


@functools.lru_cache(maxsize=None)
def _tokamax_runs_here() -> bool:
    """tokamax serves sm_80+ NVIDIA only: ROCm and sm_70/sm_75 have no implementation."""
    import jax

    try:
        return float(jax.devices()[0].compute_capability) >= 8.0
    except Exception:
        return False


def _next(name: str, why: str):
    """Hand on what ColabFold's kernels decline: tokamax if it runs here, else jax."""
    logger.debug(f"{name}: {why}, routing on")
    if _REAL is not None and _tokamax_runs_here():
        return getattr(_REAL, name)
    return _XLA[name]


def dot_product_attention(q, k, v, *, mask=None, bias=None, implementation=None,
                          scale=None, **kwargs):
    if kwargs:
        raise TypeError(f"the tokamax shim does not implement {sorted(kwargs)}")
    from colabfold_kernels import dispatch as _dispatch, fused_ops

    from colabfold.alphafold3.attention import colabfold_attention

    # alphafold3 hands attention over seq-major, so prefer the kernel that reads it
    # that way and leave q/k/v where they are
    dispatch = _dispatch()
    seq_major = True
    kernel = fused_ops.attention(dispatch, q.dtype, q.shape[-1], v.shape[-1], layout="seq")
    if kernel is None:
        seq_major = False
        kernel = fused_ops.attention(dispatch, q.dtype, q.shape[-1], v.shape[-1])
    if kernel is not None:
        out = colabfold_attention(q, k, v, mask=mask, bias=bias, scale=scale,
                                  kernel=kernel, seq_major=seq_major)
        if out is not None:
            return out
        why = "the kernel takes a shared bias only"
    else:
        why = f"no kernel for {q.dtype} at head dim {q.shape[-1]}"
    return _next("dot_product_attention", why)(
        q, k, v, mask=mask, bias=bias, implementation=implementation, scale=scale)


def gated_linear_unit(x, weights, activation=None, precision=None, **kwargs):
    """``activation(x @ weights[:, 0]) * (x @ weights[:, 1])``, as tokamax defines it."""
    if kwargs:
        raise TypeError(f"the tokamax shim does not implement {sorted(kwargs)}")
    import jax.numpy as jnp

    from colabfold_kernels import dispatch as _dispatch, fused_ops

    channels, two, out_dim = weights.shape
    kernel = fused_ops.gated_dual_proj(_dispatch(), x.dtype, activation)
    if kernel is None or two != 2 or channels % 32 or out_dim % 32:
        return _next("gated_linear_unit",
                     f"no kernel for {x.dtype} at [{channels}, {two}, {out_dim}]")(
            x=x, weights=weights, activation=activation, precision=precision)

    gate, projection = weights[:, 0, :], weights[:, 1, :]
    zero = jnp.zeros((out_dim,), jnp.float32)
    lead = x.shape[:-1]
    flat = jnp.reshape(x, (-1, channels))
    keep = jnp.ones((flat.shape[0],), x.dtype)
    out = kernel(flat, projection, zero, gate, zero, keep, split=False)
    return jnp.reshape(out, lead + (out_dim,))


def _xla_attention(q, k, v, *, mask=None, bias=None, implementation=None, scale=None):
    """jax's own attention, which already takes alphafold3's (..., seq, heads, dim)."""
    import jax

    return jax.nn.dot_product_attention(q, k, v, bias=bias, mask=mask, scale=scale)


def _xla_gated_linear_unit(x, weights, activation=None, precision=None):
    import jax.numpy as jnp

    gate = jnp.einsum("...k,kp->...p", x, weights[:, 0, :])
    projection = jnp.einsum("...k,kp->...p", x, weights[:, 1, :])
    return (activation(gate) if activation is not None else gate) * projection


_XLA = {"dot_product_attention": _xla_attention,
        "gated_linear_unit": _xla_gated_linear_unit}


def install(force: bool = False) -> str:
    """Stand in unless the real tokamax is there. Returns what will run."""
    global _REAL
    if any(name in sys.modules for name in _CALLERS):
        logger.warning("alphafold3 already imported tokamax; the shim will not take effect")
    here = sys.modules[__name__]
    try:
        import tokamax

        if tokamax is not here:
            _REAL = tokamax
    except Exception:
        # a tokamax that cannot import is no different from one that is not there
        _REAL = None
    if importlib.util.find_spec("colabfold_kernels") is None:
        logger.info("colabfold-kernels is not installed, leaving tokamax in place")
        return "tokamax"
    if not force and _REAL is not None and _tokamax_runs_here():
        return "tokamax"
    sys.modules["tokamax"] = here
    return "colabfold"
