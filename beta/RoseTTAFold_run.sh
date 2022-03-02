#!/usr/bin/env bash
set -e
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
INPUT="$(readlink -f "$1")"
shift
OUTPUT="$(readlink -f "$2")"
shift
USERNAME="$3"
shift
PASSWORD="$4"
shift
unset CONDA_SHLVL
eval "$(conda shell.bash hook)"
conda activate rosettafold
mkdir -p ${OUTPUT}
HEADER="$(cat "$INPUT" | awk 'NR == 1 { print substr($1, 2); exit; }')"
SEQUENCE="$(cat "$INPUT" | awk 'NR == 2 { print $1; exit; }')" 
if [ "${INPUT##*.}" = "a3m" ]; then
    A3M_FILE="$INPUT"
    MSA_METHOD=custom_a3m
else
    A3M_FILE=""
    MSA_METHOD=mmseqs2
fi
PYROSETTA=True
/usr/bin/time -v papermill "${SCRIPT_DIR}/RoseTTAFold.ipynb" "${OUTPUT}.ipynb" --inject-input-path -p use_pyrosetta "${PYROSETTA}" -p username "${USERNAME}" -p password "${PASSWORD}" -p jobname "$HEADER" -p sequence "$SEQUENCE" -p a3m_file "${A3M_FILE}" -p msa_method "${MSA_METHOD}" -p host_url "https://a3m.mmseqs.com" "$@" 2>&1 | tee "${OUTPUT}.log"
