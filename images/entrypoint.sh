#!/bin/sh
# Preserve the caller's working directory so relative inputs and outputs work.
set -eu
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
