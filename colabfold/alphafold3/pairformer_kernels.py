"""Give alphafold3's pairformer the fused colabfold_kernels paths: packed triangle attention between LayerNorm-fused
projections, the LayerNorm-fused triangle multiplication, and the fused transition."""

import math

# the triangle contraction of incoming edges, as alphafold3 writes it channel-major
INCOMING = "ckj,cki->cij"


def install(enable: bool = True, fp8: bool = False) -> bool:
    """Set alphafold3's pairformer hooks, or clear them; False when nothing was installed.

    fp8 is AF2's global_config.fp8: e4m3 k and v where Triton runs the attention, and the triangle contraction through
    colabfold-cuda-kernels' cuBLASLt GEMM, on compute capability 8.9 on. Lossy, so opt-in.
    """
    try:
        from alphafold3.model.network import modules
    except ImportError:
        return False
    modules.PACKED_ATTENTION = modules.PAIR_PROJECTIONS = modules.TRIMUL_KERNELS = modules.GLU_TRANSITION = None
    modules.PAD_KEEP = False
    if not enable:
        return False
    try:
        from colabfold_kernels import dispatch, fused_ops
    except ImportError:
        return False
    # trim_padding: the packed attention skips the key blocks, query blocks and rows the bucket's padding adds; the
    # pairformer masks padding everywhere downstream, so the zeros it leaves in padded rows are never read
    config = dict(dispatch(), fp8=fp8, trim_padding=True)

    # colabfold-kernels before the LayerNorm-fused kernels: the packed attention alone, the stock triangle multiplication
    fused_ln = hasattr(fused_ops, "ln_projection")

    def pair_projections(dtype, channels, num_head, head_dim):
        if not fused_ln:
            return None
        project = fused_ops.ln_projection(config, dtype, channels, num_head * head_dim)
        output_projection = fused_ops.output_projection(config, dtype, num_head * head_dim, channels)
        if project is None or output_projection is None:
            return None
        slots = getattr(fused_ops, "fp8_attention_slots", None)

        def ln_projection(x, ln_scale, ln_offset, w, wb, *keep):
            # the e4m3 slots the attention kernel that runs reads (colabfold-kernels 0.7); the pair is N x N rows
            if slots is not None:
                f8 = slots(config, dtype, head_dim, math.isqrt(x.shape[0]))
            else:
                f8 = ((1, 2), None) if fused_ops.fp8_attention(config, head_dim) else None
            f8_slots, f8_scales = f8 or ((), None)
            p, bias, *kv = project(x, ln_scale, ln_offset, w, wb, *keep, f8_slots=f8_slots,
                                   **({} if f8_scales is None else {"f8_scales": f8_scales}))
            return p, bias, kv[0] if kv else None

        return ln_projection, output_projection

    def trimul_kernels(dtype, channels, hidden):
        if not fused_ln:
            return None
        project = fused_ops.ln_gated_dual_proj(config, dtype, channels, hidden)
        epilogue = fused_ops.trimul_epilogue_ln(config, dtype)
        if project is None or epilogue is None:
            return None

        def contract(x, ln_scale, ln_offset, wp, bp, wg, bg, mask, *, equation):
            import jax.numpy as jnp
            # not the cheap 64-channel template stacks: chai1's has near-constant operands that e4m3 cannot resolve
            fp8 = fused_ops.fp8_contraction(config, dtype, channels, hidden, x.shape[1]) if channels >= 128 else None
            if fp8 is not None:
                # cuBLASLt contracts e4m3 operands along their last axis: incoming edges are projected transposed
                proj, gemm = fp8
                incoming = equation == INCOMING
                a, b = proj(x, ln_scale, ln_offset, wp, bp, wg, bg, mask, transpose=incoming)
                return gemm(*((b, a) if incoming else (a, b))).astype(x.dtype)
            a, b = project(x.reshape(-1, channels), ln_scale, ln_offset, wp, bp, wg, bg, mask.reshape(-1))
            return jnp.einsum(equation, a.reshape(hidden, *x.shape[:2]), b.reshape(hidden, *x.shape[:2]))

        return contract, epilogue

    modules.PACKED_ATTENTION = lambda dtype, head_dim: fused_ops.packed_attention(config, dtype, head_dim, head_dim)
    modules.PAIR_PROJECTIONS = pair_projections
    modules.TRIMUL_KERNELS = trimul_kernels
    modules.GLU_TRANSITION = lambda dtype, channels, hidden: fused_ops.glu_transition(config, dtype, channels, hidden)
    # the projections, the epilogue and the transition skip padded row blocks given the pair mask
    modules.PAD_KEEP = True
    return True
