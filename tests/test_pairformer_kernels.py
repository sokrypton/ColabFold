"""The pairformer hooks, run on the kernels' float32 references, compute what alphafold3's stock graph does, per lineage."""
import functools

import numpy as np
import pytest

hk = pytest.importorskip("haiku")
jax = pytest.importorskip("jax")
jnp = jax.numpy
modules = pytest.importorskip("alphafold3.model.network.modules")
model_config = pytest.importorskip("alphafold3.model.model_config")
fused_ops = pytest.importorskip("colabfold_kernels.fused_ops")
from colabfold_kernels import contract, packed_attn, tri_mul  # noqa: E402

from colabfold.alphafold3 import pairformer_kernels  # noqa: E402

# each lineage's own pair-bias and bias conventions, and chai1's parallel block
LINEAGES = ["alphafold3", "openfold3", "openbind0", "boltz2", "protenix2", "rosettafold3", "chai1"]


def counted(calls, name, f):
    """f, noting each call and whether it added the residual: add_input, or output_projection's fourth argument."""
    def wrapped(*args, **kwargs):
        calls.append((name, kwargs.get("add_input", len(args) > 3)))
        return f(*args, **kwargs)
    return wrapped


def install_references(monkeypatch, mode, calls):
    """Every selector hands back its float32 reference, so the hooks run on the CPU.

    mode: fused, fp8, or packed (no LayerNorm-fused projections: the attention's earlier packed path).
    """
    fp8 = mode == "fp8"
    def ref(name, f):
        monkeypatch.setattr(fused_ops, name, lambda *a, **k: counted(calls, name, f))

    ref("packed_attention", packed_attn.ref_packed_attention)
    ref("ln_projection", packed_attn.ref_ln_projection)
    ref("ln_gated_dual_proj", tri_mul.ref_ln_gdp)
    ref("trimul_epilogue_ln", tri_mul.ref_trimul_epilogue_ln)
    ref("output_projection", packed_attn.ref_output_projection)
    if mode == "packed":
        monkeypatch.setattr(fused_ops, "ln_projection", lambda *a: None)
    monkeypatch.setattr(fused_ops, "glu_transition", lambda *a: None)
    monkeypatch.setattr(fused_ops, "fp8_attention", lambda config, head_dim: fp8)
    proj = functools.partial(tri_mul.ref_ln_gdp, pad=16)
    monkeypatch.setattr(fused_ops, "fp8_contraction",
                        lambda *a: (counted(calls, "fp8_projection", proj), contract.ref_fp8_contract) if fp8 else None)
    assert pairformer_kernels.install(enable=True, fp8=fp8)


def pairformer(lineage, n=40, c=128, seed=0):
    # tokamax's kernels do not run on the CPU; unchunked attention, as ColabFold sets it
    gc = model_config.GlobalConfig(model=lineage, final_init="linear", flash_attention_implementation="xla",
                                   pair_attention_chunk_size=((2048, None), (None, 256)))
    no_glu = lambda equation: modules.TriangleMultiplication.Config(equation=equation, use_glu_kernel=False)
    cfg = modules.PairFormerIteration.Config(
        num_layer=1, pair_attention=modules.GridSelfAttention.Config(num_head=2),
        pair_transition=modules.TransitionBlock.Config(use_glu_kernel=False),
        triangle_multiplication_incoming=no_glu("kjc,kic->ijc"), triangle_multiplication_outgoing=no_glu("ikc,jkc->ijc"))
    f = hk.transform(lambda x, m: modules.PairFormerIteration(cfg, gc, name="pf")(x, m))
    keys = jax.random.split(jax.random.PRNGKey(seed), 4)
    x = jax.random.normal(keys[0], (n, n, c))
    seq = (jnp.arange(n) < n - 5).astype(jnp.float32)
    m = seq[:, None] * seq[None, :]
    params = f.init(keys[1], x, m)
    # random LayerNorm scales and biases too, so a misplaced one shows
    leaves, tree = jax.tree.flatten(params)
    noise = [jax.random.normal(k, p.shape) * 0.3 + (p if p.ndim == 1 else 0)
             for k, p in zip(jax.random.split(keys[2], len(leaves)), leaves)]
    return f, jax.tree.unflatten(tree, noise), x, m


@pytest.mark.parametrize("mode", ["fused", "fp8", "packed"])
@pytest.mark.parametrize("lineage", LINEAGES)
def test_hooks_match_stock(monkeypatch, lineage, mode):
    f, params, x, m = pairformer(lineage)
    with jax.default_matmul_precision("highest"):
        pairformer_kernels.install(enable=False)
        stock = np.asarray(f.apply(params, jax.random.PRNGKey(5), x, m))
        calls = []
        install_references(monkeypatch, mode, calls)
        try:
            fused = np.asarray(f.apply(params, jax.random.PRNGKey(5), x, m))
        finally:
            pairformer_kernels.install(enable=False)
    names = [name for name, _ in calls]
    parallel = lineage == "chai1"
    assert names.count("trimul_epilogue_ln") == 2
    assert names.count("fp8_projection") == (2 if mode == "fp8" else 0)
    assert names.count("packed_attention") == 2
    assert names.count("ln_projection") == (0 if mode == "packed" else 2)
    # the residual goes into the kernels of a sequential block, chai1's parallel one sums the updates itself
    assert all(add_input != parallel for name, add_input in calls if name in ("trimul_epilogue_ln", "output_projection"))
    rel = np.abs(fused - stock).max() / np.abs(stock).max()
    assert rel < 1e-4, rel
