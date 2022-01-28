#!/usr/bin/env bash
set -e

# run this script from the MsaServer folder
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "${SCRIPT_DIR}"

# set which commits to use
MMSEQS_COMMIT=${1:-edb8223d1ea07385ffe63d4f103af0eb12b2058e}
BACKEND_COMMIT=${2:-1d84d23ec1199a9e46df9a5cb3301f6d73f3530d}

# check if all dependencies are there
for i in go curl git aria2c; do
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
    curl -s -o- https://mmseqs.com/archive/${COMMIT}/mmseqs-linux-avx2.tar.gz | tar -xf - mmseqs/bin/mmseqs
  fi
else
  curl -s -o- https://mmseqs.com/archive/${COMMIT}/mmseqs-linux-avx2.tar.gz | tar -xf - mmseqs/bin/mmseqs
fi

# mmseqs needs to be in PATH for the setup_databases script to work
PATH="${SCRIPT_DIR}/mmseqs/bin:$PATH"

# don't re-download databases if they already exist as they are quite large
if [ ! -d databases ]; then
  mkdir -p databases
  ../setup_databases.sh databases
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
