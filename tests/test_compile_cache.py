import os

import jax
import jax.numpy as jnp
import pytest

from colabfold.alphafold import compile_cache
from colabfold.alphafold.compile_cache import CachedApply, clean, install_hash

needs_permissions = pytest.mark.skipif(not hasattr(os, "geteuid") or os.geteuid() == 0,
                                       reason="root ignores file permissions")


def apply(x):
    return jnp.sin(x) * 2


def entries(cache):
    return sorted(cache.glob("*/af2-*.exec"))


def test_stores_then_loads(tmp_path, caplog):
    x = jnp.ones((4, 8))
    out = CachedApply(jax.jit(apply), tmp_path, "env")(x)
    assert len(entries(tmp_path)) == 1 and entries(tmp_path)[0].parent.name == install_hash()
    with caplog.at_level("INFO"):
        again = CachedApply(jax.jit(apply), tmp_path, "env")(x)
    assert "Loaded the compiled model" in caplog.text and (again == out).all()
    CachedApply(jax.jit(apply), tmp_path, "env")(jnp.ones((5, 8)))
    CachedApply(jax.jit(apply), tmp_path, "other config")(x)
    assert len(entries(tmp_path)) == 3


@needs_permissions
@pytest.mark.parametrize("mode", [0o555, 0o000])
def test_read_only_or_unreadable_cache_runs(tmp_path, mode, caplog):
    CachedApply(jax.jit(apply), tmp_path, "env")(jnp.ones((4, 8)))
    dirs = [tmp_path, *(p for p in tmp_path.rglob("*") if p.is_dir())]
    for d in reversed(dirs):
        d.chmod(mode)
    try:
        cached = CachedApply(jax.jit(apply), tmp_path, "env")
        with caplog.at_level("INFO"):
            cached(jnp.ones((4, 8)))
            cached(jnp.ones((5, 8)))
            cached(jnp.ones((6, 8)))
        assert ("Loaded the compiled model" in caplog.text) == (mode == 0o555)
        assert caplog.text.count("Could not store") == 1
    finally:
        for d in dirs:
            d.chmod(0o755)


def test_uncreatable_cache_runs(tmp_path):
    (tmp_path / "file").write_text("")
    assert CachedApply(jax.jit(apply), tmp_path / "file" / "cache", "env")(jnp.ones(3)).shape == (3,)


def test_trusts_only_what_the_cache_owner_alone_could_write(tmp_path, monkeypatch, caplog):
    CachedApply(jax.jit(apply), tmp_path, "env")(jnp.ones(3))
    monkeypatch.setattr(compile_cache, "_owned", lambda path: False)  # now another user's entries
    for path in [tmp_path, *tmp_path.rglob("*")]:
        path.chmod(0o755 if path.is_dir() else 0o644)
    with caplog.at_level("INFO"):
        CachedApply(jax.jit(apply), tmp_path, "env")(jnp.ones(3))
    assert "Loaded the compiled model" in caplog.text
    caplog.clear()
    tmp_path.chmod(0o777)
    with caplog.at_level("INFO"):
        CachedApply(jax.jit(apply), tmp_path, "env")(jnp.ones(3))
    assert "Loaded the compiled model" not in caplog.text and "not loaded" in caplog.text
    tmp_path.chmod(0o755)


def test_clean_keeps_this_installation_and_clean_all_nothing(tmp_path):
    CachedApply(jax.jit(apply), tmp_path, "env")(jnp.ones(3))
    old = tmp_path / "0123456789abcdef"
    old.mkdir()
    (old / "af2-0.exec").write_bytes(b"")
    assert clean(tmp_path, keep=install_hash()) == 1 and not old.exists() and len(entries(tmp_path)) == 1
    assert clean(tmp_path) == 1 and not entries(tmp_path) and not any(tmp_path.iterdir())
