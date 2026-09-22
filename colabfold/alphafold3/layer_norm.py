"""
Run alphafold3's LayerNorm on ColabFold's fused kernel.

alphafold3 upcasts 16-bit activations to float32, normalises, then casts back.
ColabFold's kernel does the same arithmetic in one pass without materialising the
float32 tensor. :func:`interceptor` swaps the two where the kernel fits; it is a
haiku method interceptor because the parameters have to be fetched inside the
module's own name scope.
"""
import logging

logger = logging.getLogger(__name__)


def _fused(module, x):
    """The kernel's answer, or None when it cannot serve this layer norm."""
    import haiku as hk
    import jax.numpy as jnp

    from colabfold_kernels import dispatch as _dispatch, fused_ops

    half = x.dtype in (jnp.bfloat16, jnp.float16)
    if half and not getattr(module, "upcast", False):
        return None  # the kernel always normalises in float32
    if module.param_axis not in (None, (-1,), [-1]):
        return None
    kernel = fused_ops.layer_norm(_dispatch(), x.dtype)
    if kernel is None:
        return None

    # alphafold3 creates these after its upcast, so they are float32 in the checkpoint
    channels = x.shape[-1]
    scale = (hk.get_parameter("scale", (channels,), jnp.float32, init=module.scale_init)
             if module._temp_create_scale else jnp.ones((channels,), jnp.float32))
    offset = (hk.get_parameter("offset", (channels,), jnp.float32, init=module.offset_init)
              if module._temp_create_offset else jnp.zeros((channels,), jnp.float32))
    return kernel(x, scale, offset, eps=module.eps)


def interceptor(next_fun, args, kwargs, context):
    """Hand alphafold3's LayerNorm calls to the kernel, or pass them through."""
    try:
        from alphafold3.model.components import haiku_modules as hm
    except ModuleNotFoundError:
        return next_fun(*args, **kwargs)

    if (context.method_name == "__call__" and isinstance(context.module, hm.LayerNorm)
            and len(args) == 1 and not kwargs):
        out = _fused(context.module, args[0])
        if out is not None:
            return out
    return next_fun(*args, **kwargs)
