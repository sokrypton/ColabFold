#!/bin/sh
# Preserve the caller's working directory so relative inputs and outputs work.
set -eu
# Resolve the cache before Python imports ColabFold's default_data_dir.
# Honor explicit XDG_CACHE_HOME, including pre-populated read-only caches.
if [ -z "${XDG_CACHE_HOME:-}" ]; then
    caat_cache_writable() {
        mkdir -p "$1/colabfold" 2>/dev/null || return 1
        CAAT_CACHE_PROBE=$(mktemp "$1/colabfold/.write-check.XXXXXX" 2>/dev/null) || return 1
        rm -f "$CAAT_CACHE_PROBE"
    }
    if [ -n "${HOME:-}" ] && caat_cache_writable "$HOME/.cache"; then
        XDG_CACHE_HOME="$HOME/.cache"
    elif caat_cache_writable "$PWD/.cache"; then
        XDG_CACHE_HOME="$PWD/.cache"
    else
        XDG_CACHE_HOME=$(mktemp -d "${TMPDIR:-/tmp}/caat-cache.XXXXXX")
        echo "CAAT: using temporary cache $XDG_CACHE_HOME; set --data-dir for persistent weights." >&2
    fi
    export XDG_CACHE_HOME
fi
if [ "$#" -eq 0 ]; then
    set -- --help
fi
case "$1" in
    scripts/run_e2e_pipeline.py|run_e2e_pipeline.py|scripts/run_attention_heads.py|run_attention_heads.py|scripts/run_analysis_pipeline.py|run_analysis_pipeline.py)
        CAAT_SCRIPT=${1##*/}
        shift
        exec python "/opt/CAAT/scripts/$CAAT_SCRIPT" "$@"
        ;;
    python|python3)
        shift
        exec python "$@"
        ;;
    colabfold_batch|colabfold_search|colabfold_split_msas|colabfold_relax)
        exec "$@"
        ;;
    -*)
        exec python /opt/CAAT/scripts/run_e2e_pipeline.py "$@"
        ;;
    *)
        exec python "$@"
        ;;
esac
