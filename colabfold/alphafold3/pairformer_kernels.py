"""Give alphafold3's pairformer the fused colabfold_kernels paths: packed triangle attention, fused triangle multiplication
and transition."""

# the epilogue holds whole [hidden, channels] weight tiles: 128 fits 99 KiB of shared memory, 256 does not
MAX_TRIMUL_WIDTH = 128


def install(enable: bool = True) -> bool:
    """Set alphafold3's pairformer hooks, or clear them; False when nothing was installed."""
    try:
        from alphafold3.model.network import modules
    except ImportError:
        return False
    modules.PACKED_ATTENTION = modules.TRIMUL_KERNELS = modules.GLU_TRANSITION = None
    if not enable:
        return False
    try:
        from colabfold_kernels import dispatch, fused_ops
    except ImportError:
        return False
    config = dispatch()

    def trimul_kernels(dtype, channels, hidden):
        if max(channels, hidden) > MAX_TRIMUL_WIDTH:
            return None
        gated_dual_proj = fused_ops.gated_dual_proj(config, dtype)
        epilogue = fused_ops.trimul_epilogue(config, dtype)
        return None if gated_dual_proj is None or epilogue is None else (gated_dual_proj, epilogue)

    modules.PACKED_ATTENTION = lambda dtype, head_dim: fused_ops.packed_attention(config, dtype, head_dim, head_dim)
    modules.TRIMUL_KERNELS = trimul_kernels
    modules.GLU_TRANSITION = lambda dtype, channels, hidden: fused_ops.glu_transition(config, dtype, channels, hidden)
    return True
