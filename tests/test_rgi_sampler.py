"""Exercise the actual upstream sampler under JIT, vmap, and scan.

`guided` takes alphafold3-colabfold's DENOISED_HOOK when the installed release
has it (3.1.18+) and instruments the pinned 3.1.11 sampler otherwise, so these
hold the same expectations against whichever is installed."""

from types import SimpleNamespace

import numpy as np
import pytest

from colabfold.rgi.sampler import coordinate_minimizer, guided, restraint_layout


def test_bioemu_atom_layout_excludes_missing_terminal_atoms_without_mutation():
    example = {
        "ref_mask": np.array([[1, 1, 1], [1, 0, 0]]),
        "bioemu2/af3_index": np.array([[2, 0, -1], [1, -1, -1]]),
        "bioemu2/atom_mask": np.array([1, 1, 1, 0]),
    }
    adapted, gather = restraint_layout(example, "bioemu2")
    np.testing.assert_array_equal(adapted["ref_mask"], [[1, 1, 0], [1, 0, 0]])
    assert example["ref_mask"][0, 2] == 1
    np.testing.assert_array_equal(gather, example["bioemu2/af3_index"])


@pytest.mark.parametrize("indices", [[[0, 0]], [[0, 3]], [[0, 2]], [[0.0, 1.0]]])
def test_bioemu_rejects_invalid_atom_maps(indices):
    example = {
        "ref_mask": np.ones((1, 2)),
        "bioemu2/af3_index": np.asarray(indices),
        "bioemu2/atom_mask": np.array([1, 1, 0]),
    }
    with pytest.raises(ValueError, match="atom map"):
        restraint_layout(example, "bioemu2")


def test_bioemu_rejects_samplers_that_bypass_the_guidance_hook():
    with pytest.raises(ValueError, match="BIOEMU2_SAMPLER=edm"):
        restraint_layout({"bioemu2steer/preset_edm_gd": np.zeros(0)}, "bioemu2")
    with pytest.raises(ValueError, match="3.1.18"):
        restraint_layout({}, "bioemu2")


def test_bioemu_coordinate_mapping_is_preserved_under_jit_and_vmap():
    jax = pytest.importorskip("jax")
    jnp = pytest.importorskip("jax.numpy")
    coords = jnp.arange(24, dtype=jnp.float32).reshape(2, 4, 1, 3)
    gather = jnp.array([[2, 0, -1], [1, -1, -1]])
    offsets = jnp.arange(6, dtype=jnp.float32)[:, None]
    apply = coordinate_minimizer(lambda x, sigma, step: x + offsets, gather)
    actual = jax.jit(jax.vmap(lambda x: apply(x, 1.0, 0)))(coords)
    expected = np.asarray(coords).copy()
    expected[:, 0, 0] += 1
    expected[:, 1, 0] += 3
    np.testing.assert_array_equal(actual, expected)


def test_structural_token_coordinates_preserve_residue_mapping():
    jax = pytest.importorskip("jax")
    jnp = pytest.importorskip("jax.numpy")
    coords = jnp.arange(18, dtype=jnp.float32).reshape(3, 2, 3)
    example = {
        "ref_mask": np.array([[1, 1], [1, 0]]),
        "structbook/residue_atom_gather": np.array([[4, 1], [3, -1]]),
    }
    adapted, gather = restraint_layout(example, "opendde")
    assert adapted is example
    gather = jnp.asarray(gather)
    apply = coordinate_minimizer(lambda x, sigma, step: x + 10, gather)
    actual = jax.jit(apply)(coords, 1.0, 0).reshape(-1, 3)
    expected = np.asarray(coords).reshape(-1, 3).copy()
    expected[[4, 1, 3]] += 10
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize(
    "model", ["openbind0", "chai1", "esmfold2", "bioemu2", "protenix1"]
)
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

    @hk.transform
    def forward():
        return dh.sample(lambda x, sigma: x / 2, batch, hk.next_rng_key(), config, gc)[
            "atom_positions"
        ]

    original = dh.sample
    key = jax.random.PRNGKey(42)
    plain = jax.jit(forward.apply)({}, key)
    with guided(dh, lambda x, sigma, step: x):
        hooked = jax.jit(forward.apply)({}, key)
    np.testing.assert_array_equal(plain, hooked)
    assert dh.sample is original
    assert getattr(dh, "DENOISED_HOOK", None) is None


def test_schedule_passes_pre_churn_sigma_and_zero_based_steps():
    dh = pytest.importorskip("alphafold3.model.network.diffusion_head")
    import haiku as hk
    import jax
    import jax.numpy as jnp

    seen = []

    def minimize(x, sigma, step):
        jax.debug.callback(lambda s, i: seen.append((float(s), int(i))), sigma, step)
        return x

    batch = SimpleNamespace(
        predicted_structure_info=SimpleNamespace(atom_mask=jnp.ones((2, 2)))
    )
    config = dh.SampleConfig(steps=6, num_samples=1)

    @hk.transform
    def forward():
        return dh.sample(
            lambda x, sigma: x / 2,
            batch,
            hk.next_rng_key(),
            config,
            SimpleNamespace(model="openbind0"),
        )

    with guided(dh, minimize):
        jax.block_until_ready(jax.jit(forward.apply)({}, jax.random.PRNGKey(42)))
    expected = dh.noise_schedule(np.linspace(0, 1, 7, dtype=np.float32))[:-1]
    assert sorted(i for s, i in seen) == list(range(6))
    np.testing.assert_allclose(
        [s for s, i in sorted(seen, key=lambda x: x[1])], expected
    )
