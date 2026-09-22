"""
Make tokamax's kernels fit a GPU with Ampere-sized shared memory.

It only shrinks its staging on sm80, and a workstation Blackwell has the same
shared memory with a capability number none of its branches match.
"""
import logging

logger = logging.getLogger(__name__)

# 99 KiB on Ampere and workstation Blackwell, 228 on the datacenter parts
SMALL_SHARED_MEMORY = 128 * 1024
MAX_STAGES = 2


def install() -> bool:
    """Cap how deeply tokamax stages, wherever it picks a config."""
    try:
        from tokamax._src.ops import op as tokamax_op
    except ImportError:
        return False
    if getattr(tokamax_op, "_colabfold_capped_stages", False):
        return True

    from colabfold_kernels import shared_memory_limit

    limit = shared_memory_limit()
    if limit is None or limit >= SMALL_SHARED_MEMORY:
        return False

    import dataclasses

    original = tokamax_op.BoundArguments.get_config

    def capped(self, *args, **kwargs):
        config = original(self, *args, **kwargs)
        if getattr(config, "num_stages", 0) > MAX_STAGES:
            config = dataclasses.replace(config, num_stages=MAX_STAGES)
        return config

    tokamax_op.BoundArguments.get_config = capped
    tokamax_op._colabfold_capped_stages = True
    logger.info(f"capping tokamax at {MAX_STAGES} stages, for {limit >> 10} KiB of shared memory")
    return True
