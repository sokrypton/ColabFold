#!/usr/bin/env bash
set -e

# run this script from the MsaServer folder
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "${SCRIPT_DIR}"

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
MMSEQS_COMMIT=${1:-4589151554eb83a70ff0c4d04d21b83cabc203e4}
BACKEND_COMMIT=${2:-14e087560f309f989a5e1feb54fd1f9c988076d5}

# check if all dependencies are there
for i in go curl git aria2c rsync; do
  if ! command -v "${i}" > /dev/null 2>&1; then
    echo "${i} is not installed, please install it first"
    exit 1
  fi 
done

# check that the correct mmseqs commit is there
if [ -x ./mmseqs/bin/mmseqs ]; then
  # download it again if its a different commit
  if [ $(./mmseqs/bin/mmseqs version) != ${MMSEQS_COMMIT} ]; then 
    rm -rf -- mmseqs
    curl -s -o- https://mmseqs.com/archive/${MMSEQS_COMMIT}/mmseqs-linux-avx2.tar.gz | tar -xzf - mmseqs/bin/mmseqs
  fi
else
  curl -s -o- https://mmseqs.com/archive/${MMSEQS_COMMIT}/mmseqs-linux-avx2.tar.gz | tar -xzf - mmseqs/bin/mmseqs
fi

# mmseqs needs to be in PATH for the setup_databases script to work
PATH="${SCRIPT_DIR}/mmseqs/bin:$PATH"

# don't re-download databases if they already exist as they are quite large
if [ ! -d databases ]; then
  mkdir -p databases
  ../setup_databases.sh databases "${PDB_SERVER}" "${PDB_PORT}"
fi

# make sure the API server is checked out
if [ -d mmseqs-server ]; then
  pushd mmseqs-server
  git pull
  git checkout ${BACKEND_COMMIT}
  popd
else
  git clone https://github.com/soedinglab/MMseqs2-App.git mmseqs-server
  pushd mmseqs-server
  git checkout ${BACKEND_COMMIT}
  popd
fi

# compile the api server
pushd mmseqs-server/backend
go build -o ../../msa-server
popd

# start the server in local mode
# meaning both workers and server run from the same process
./msa-server -local -config config.json
