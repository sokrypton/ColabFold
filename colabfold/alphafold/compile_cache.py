"""AF2's compiled forward pass on disk, so a shape compiles once per machine.

JAX's persistent cache keys on the lowered module, so even its hits trace and lower the model (~2-3 s a shape). This
one loads the serialised executable (~0.5 s). Entries live in a directory per installation (GPU, driver, jax/jaxlib,
XLA_FLAGS, the sources of the model, kernel and colabfold packages), keyed there by the model config and input shapes.
Nothing is evicted: an upgrade starts a new directory, and colabfold_compile --clean removes the old ones.
"""
import functools
import hashlib
import json
import logging
import os
import pickle
import tempfile
import time
from pathlib import Path

logger = logging.getLogger(__name__)

# read by the host loop around the model, not traced: a run with more recycles reuses the executable
_HOST_ONLY = (("model", "num_recycle"), ("model", "recycle_early_stop_tolerance"), ("model", "stop_at_score"),
              ("data", "common", "num_recycle"))

_warned = set()


def _warn_once(key, message):
    if key not in _warned:
        _warned.add(key)
        logger.warning(message)


def _drop_crashed_writes(directory):
    """Our .tmp files a day old: writes that crashed before their rename."""
    for path in directory.glob("af2-*.tmp"):
        try:
            if time.time() - path.stat().st_mtime > 24 * 3600 and _owned(path):
                path.unlink()
        except OSError:  # removed by a concurrent run
            pass


def _owned(path) -> bool:
    return not hasattr(os, "getuid") or path.stat().st_uid == os.getuid()


def _trusted(path, cache_dir) -> bool:
    """Unpickling runs code: load what this user wrote, or what only the cache directory's owner could have written
    (a shared cache an administrator filled): same owner, nothing on the way writable by group or others."""
    if _owned(path):
        return True
    owner = Path(cache_dir).stat().st_uid
    return all(p.stat().st_uid == owner and not p.stat().st_mode & 0o022 for p in (path, path.parent, Path(cache_dir)))


class Precompiled(Exception):
    """Raised in place of running the model when only compiling (colabfold_compile)."""


@functools.lru_cache(maxsize=None)
def install_hash() -> str:
    """A hash of what every executable of this installation depends on: its entries' directory name."""
    import importlib.util
    import jax
    import jaxlib
    h = hashlib.sha256()
    for package in ("alphafold", "colabfold", "colabfold_kernels", "colabfold_legacy_kernels", "colabfold_cuda_kernels"):
        spec = importlib.util.find_spec(package)
        if spec is None or spec.origin is None:
            continue
        root = Path(spec.origin).parent
        for path in sorted(root.rglob("*")):
            if path.suffix == ".py":
                h.update(path.read_bytes())
            elif ".so" in path.suffixes or path.suffix in (".hsaco", ".cubin"):
                # kernel libraries by size and modification time: reading them would cost every start
                stat = path.stat()
                h.update(f"{path.relative_to(root)} {stat.st_size} {stat.st_mtime_ns}".encode())
    device = jax.devices()[0]
    h.update(json.dumps([jax.__version__, jaxlib.__version__, device.device_kind, device.client.platform_version,
                         os.environ.get("XLA_FLAGS", ""), str(jax.config.jax_default_matmul_precision)]).encode())
    return h.hexdigest()[:16]


def environment(config, jit_kwargs) -> str:
    """A hash of the model config and compiler options an executable was built with."""
    config = json.loads(config.to_json_best_effort())
    for *parents, key in _HOST_ONLY:
        node = config
        for parent in parents:
            node = node.get(parent, {})
        node.pop(key, None)
    return hashlib.sha256(json.dumps([config, repr(sorted((jit_kwargs or {}).items()))], sort_keys=True).encode()
                          ).hexdigest()


def clean(cache_dir, keep=None) -> int:
    """Remove this user's entries from cache_dir, all but those in the installation directory keep. -> how many."""
    removed = 0
    for directory in Path(cache_dir).glob("*/"):
        if directory.name == keep:
            continue
        for path in directory.glob("af2-*"):
            try:
                if _owned(path):
                    path.unlink()
                    removed += path.suffix == ".exec"
            except OSError:  # removed by a concurrent run, or not ours to remove
                pass
        try:
            directory.rmdir()  # once empty: another user's entries keep theirs
        except OSError:
            pass
    return removed


class CachedApply:
    """RunModel.apply's replacement: the executable for these input shapes from memory, else disk, else compiled
    (and stored). compile_only raises Precompiled once the executable exists instead of running it. A cache it cannot
    read or write costs only the compile."""

    def __init__(self, jitted, cache_dir, env, compile_only=False):
        self.jitted, self.root, self.env, self.compile_only = jitted, Path(cache_dir), env, compile_only
        self.dir = self.root / install_hash()
        self.loaded = {}

    def _load(self, path):
        from jax.experimental import serialize_executable
        try:
            if not path.is_file():
                return None
            if not _trusted(path, self.root):
                _warn_once(("untrusted", self.root), f"Compiled models in {self.root} that another user could have "
                                                      "written are not loaded; compiling them again")
                return None
            compiled = serialize_executable.deserialize_and_load(*pickle.loads(path.read_bytes()))
        except Exception as e:  # unreadable, torn or stale: compile again
            logger.warning(f"Compiled model {path} did not load ({type(e).__name__}: {e}); compiling it again")
            return None
        logger.info(f"Loaded the compiled model from {path}")
        return compiled

    def _store(self, compiled, path):
        from jax.experimental import serialize_executable
        try:
            self.dir.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(dir=self.dir, prefix="af2-", suffix=".tmp", delete=False) as handle:
                pickle.dump(serialize_executable.serialize(compiled), handle)
            os.replace(handle.name, path)  # atomic: concurrent writers of one shape leave one whole entry
            _drop_crashed_writes(self.dir)
        except Exception as e:  # read-only or shared cache: runs on, compiling what it lacks
            _warn_once(("store", self.dir), f"Could not store compiled models in {self.dir} "
                                            f"({type(e).__name__}: {e}); they compile each run")

    def _executable(self, args):
        import jax
        leaves, treedef = jax.tree_util.tree_flatten(args)
        signature = repr((self.env, str(treedef), [(x.shape, str(x.dtype), getattr(x, "weak_type", False))
                                                   for x in leaves]))
        if signature not in self.loaded:
            path = self.dir / f"af2-{hashlib.sha256(signature.encode()).hexdigest()[:32]}.exec"
            compiled = self._load(path)
            if compiled is None:
                compiled = self.jitted.lower(*jax.tree.map(jax.ShapeDtypeStruct.like, args)).compile()
                self._store(compiled, path)
            self.loaded[signature] = compiled
        return self.loaded[signature]

    def __call__(self, *args):
        compiled = self._executable(args)
        if self.compile_only:
            raise Precompiled()
        return compiled(*args)


def install(model_runner, cache_dir, jit_kwargs=None, compile_only=False):
    """Route model_runner.apply through the cache. A loaded executable calls kernels untraced, so first load what
    tracing would: Mosaic GPU's runtime, and the FFI kernels the config can call."""
    try:
        import jax.experimental.mosaic.gpu  # noqa: F401
    except ImportError:  # no Mosaic GPU (ROCm): nothing of it to load
        pass
    from colabfold_kernels import fused_ops
    fused_ops.register_ffi_targets(model_runner.config.model.global_config)
    env = environment(model_runner.config, jit_kwargs)
    model_runner.apply = CachedApply(model_runner.apply, cache_dir, env, compile_only)
    if model_runner.multimer_mode and model_runner.config.model.get("split_msa_sampling", False):
        # the multimer's MSA sampler, compiled per depth bucket: kept beside the model, run even when only compiling
        from alphafold.model import msa_sampling
        sampler = msa_sampling.make_sampler(model_runner.config.model.embeddings_and_evoformer)
        model_runner._sampler = CachedApply(sampler, cache_dir, env + "-msa-sampler")
