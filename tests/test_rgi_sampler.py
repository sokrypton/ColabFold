"""Exercise the actual pinned sampler under JIT, vmap, and scan."""

from types import SimpleNamespace

import numpy as np
import pytest

from colabfold.rgi.sampler import coordinate_minimizer, make_guided_sampler


def test_structural_token_coordinates_preserve_residue_mapping():
    jax = pytest.importorskip("jax")
    jnp = pytest.importorskip("jax.numpy")
    coords = jnp.arange(18, dtype=jnp.float32).reshape(3, 2, 3)
    gather = jnp.array([[4, 1], [3, -1]])
    apply = coordinate_minimizer(lambda x, sigma, step: x + 10, gather)
    actual = jax.jit(apply)(coords, 1.0, 0).reshape(-1, 3)
    expected = np.asarray(coords).reshape(-1, 3).copy()
    expected[[4, 1, 3]] += 10
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("model", ["openbind0", "chai1", "esmfold2"])
def test_identity_hook_preserves_upstream_sampling_paths(model):
    dh = pytest.importorskip("alphafold3.model.network.diffusion_head")
    import haiku as hk
    import jax
    import jax.numpy as jnp

    config = dh.SampleConfig(
        steps=6, num_samples=2, max_sigma=256 if model.startswith("esmfold2") else 0
    )
    batch = SimpleNamespace(
        predicted_structure_info=SimpleNamespace(atom_mask=jnp.ones((2, 2)))
    )
    gc = SimpleNamespace(model=model)

    def make(sample):
        @hk.transform
        def forward():
            return sample(lambda x, sigma: x / 2, batch, hk.next_rng_key(), config, gc)[
                "atom_positions"
            ]

        return jax.jit(forward.apply)

    original = dh.sample
    guided = make_guided_sampler(original, lambda x, sigma, step: x)
    key = jax.random.PRNGKey(42)
    np.testing.assert_array_equal(make(original)({}, key), make(guided)({}, key))
    assert dh.sample is original


def test_schedule_passes_pre_churn_sigma_and_zero_based_steps():
    dh = pytest.importorskip("alphafold3.model.network.diffusion_head")
    import haiku as hk
    import jax
    import jax.numpy as jnp

    seen = []

    def minimize(x, sigma, step):
        jax.debug.callback(lambda s, i: seen.append((float(s), int(i))), sigma, step)
        return x

    guided = make_guided_sampler(dh.sample, minimize)
    batch = SimpleNamespace(
        predicted_structure_info=SimpleNamespace(atom_mask=jnp.ones((2, 2)))
    )
    config = dh.SampleConfig(steps=6, num_samples=1)

    @hk.transform
    def forward():
        return guided(
            lambda x, sigma: x / 2,
            batch,
            hk.next_rng_key(),
            config,
            SimpleNamespace(model="openbind0"),
        )

    jax.block_until_ready(jax.jit(forward.apply)({}, jax.random.PRNGKey(42)))
    expected = dh.noise_schedule(np.linspace(0, 1, 7, dtype=np.float32))[:-1]
    assert sorted(i for s, i in seen) == list(range(6))
    np.testing.assert_allclose(
        [s for s, i in sorted(seen, key=lambda x: x[1])], expected
    )
