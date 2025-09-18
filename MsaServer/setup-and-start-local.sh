#!/usr/bin/env bash
set -e

# run this script from the MsaServer folder
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "${SCRIPT_DIR}"

# enable this if you want to start a GPU MSA server 
#export GPU=1

# choose which pdb rsync server to use
#PDB_SERVER=rsync.wwpdb.org::ftp                                   # RCSB PDB server name
#PDB_PORT=33444                                                    # port RCSB PDB server is using
#
#PDB_SERVER=rsync.ebi.ac.uk::pub/databases/rcsb/pdb-remediated     # PDBe server name
#PDB_PORT=873                                                      # port PDBe server is using
#
#PDB_SERVER=pdb.protein.osaka-u.ac.jp::ftp                         # PDBj server name
#PDB_PORT=873                                                      # port PDBj server is using
if [ -z "${PDB_SERVER}" ] || [ -z "${PDB_PORT}" ]; then
    echo "PDB rsync server was not chosen, please edit this script to choose which PDB download server you want to use"
    exit 1
fi

# set which commits to use
MMSEQS_COMMIT=${1:-05ae20cbc628d9911ad0aa421fba029cc457b76e}
BACKEND_COMMIT=${2:-01365aa4735539ba95b417f73fb5326c77410394}

# check if all dependencies are there
for i in curl aria2c rsync aws; do
  if ! command -v "${i}" > /dev/null 2>&1; then
    echo "${i} is not installed, please install it first"
    exit 1
  fi 
done

ARCH=avx2
if [ -n "${GPU}" ]; then
  ARCH=gpu
fi
ARCH_SERVER=x86_64;
case "$(uname -m)" in
  aarch64|arm64)
    ARCH=arm64
    if [ -n "${GPU}" ]; then
      ARCH=gpu-arm64
    fi
    ARCH_SERVER=arm64
  ;;
esac

OS=linux
case "$(uname -s)" in
  Darwin)
    OS=osx
    ARCH=universal
    ARCH_SERVER=universal
    unset GPU
  ;;
esac

# check that the correct mmseqs commit is there
if [ -x ./mmseqs/bin/mmseqs ]; then
  # download it again if its a different commit
  if [ $(./mmseqs/bin/mmseqs version) != ${MMSEQS_COMMIT} ]; then 
    rm -rf -- mmseqs
    curl -s -o- https://mmseqs.com/archive/${MMSEQS_COMMIT}/mmseqs-${OS}-${ARCH}.tar.gz | tar -xzf - mmseqs/bin/mmseqs
  fi
else
  curl -s -o- https://mmseqs.com/archive/${MMSEQS_COMMIT}/mmseqs-${OS}-${ARCH}.tar.gz | tar -xzf - mmseqs/bin/mmseqs
fi

# check that the correct mmseqs commit is there
if [ -x ./mmseqs-server/bin/mmseqs-server ]; then
  # download it again if its a different commit
  if [ $(./mmseqs-server/bin/mmseqs-server -version) != ${BACKEND_COMMIT} ]; then 
    rm -rf -- mmseqs-server
    curl -s -o- https://mmseqs.com/archive/${BACKEND_COMMIT}/mmseqs-server-${OS}-${ARCH_SERVER}.tar.gz | tar -xzf - mmseqs-server/bin/mmseqs-server
  fi
else
  curl -s -o- https://mmseqs.com/archive/${BACKEND_COMMIT}/mmseqs-server-${OS}-${ARCH_SERVER}.tar.gz | tar -xzf - mmseqs-server/bin/mmseqs-server
fi


# mmseqs needs to be in PATH for the setup_databases script to work
PATH="${SCRIPT_DIR}/mmseqs/bin:$PATH"

# This downloads a tiny swissprot based database for testing only
#export DEBUG_MINI_DB=1
# don't re-download databases if they already exist as they are large
if [ ! -d databases ]; then
  mkdir -p databases
  if [ -n "${DEBUG_MINI_DB}" ]; then
    touch databases/SKIP_TEMPLATES
  fi
  ../setup_databases.sh databases "${PDB_SERVER}" "${PDB_PORT}"
fi

# Extra GPU parameters
GPU_PARAMS=""
if [ -n "${GPU}" ]; then
  GPU_PARAMS="--paths.colabfold.gpu.gpu 1 --paths.colabfold.gpu.server 1"
fi

# start the server in local mode
# meaning both workers and server run from the same process
if [ -n "${DEBUG_MINI_DB}" ]; then
  ./mmseqs-server/bin/mmseqs-server -local -config config.json -paths.colabfold.pdb70 '' ${GPU_PARAMS}
else
  ./mmseqs-server/bin/mmseqs-server -local -config config.json ${GPU_PARAMS}
fi
