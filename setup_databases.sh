#!/bin/bash -e
# Setup everything for using mmseqs locally
# Set MMSEQS_NO_INDEX to skip the index creation step (not useful for colabfold_search in most cases)
ARIA_NUM_CONN=8
WORKDIR="${1:-$(pwd)}"

PDB_SERVER="${2:-"rsync.wwpdb.org::ftp"}"
PDB_PORT="${3:-"33444"}"
# do initial download of the PDB through aws?
# still syncs latest structures through rsync
PDB_AWS_DOWNLOAD="${4:-}"
PDB_AWS_SNAPSHOT="20240101"

UNIREF30DB="uniref30_2302"
MMSEQS_NO_INDEX=${MMSEQS_NO_INDEX:-}
DOWNLOADS_ONLY=${DOWNLOADS_ONLY:-}
# download prebuilt MMseqs2 databases that support both CPU and GPU
# these can result in different results on MMseqs2-CPU due to database ordering
FAST_PREBUILT_DATABASES=${FAST_PREBUILT_DATABASES:-"1"}
GPU=${GPU:-}
mkdir -p -- "${WORKDIR}"
cd "${WORKDIR}"

hasCommand () {
    command -v "$1" >/dev/null 2>&1
}

fail() {
    echo "Error: $1"
    exit 1
}

if ! hasCommand mmseqs; then
  fail "mmseqs command not found in PATH. Please install MMseqs2."
fi

STRATEGY=""
if hasCommand aria2c; then STRATEGY="$STRATEGY ARIA"; fi
if hasCommand curl;   then STRATEGY="$STRATEGY CURL"; fi
if hasCommand wget;   then STRATEGY="$STRATEGY WGET"; fi
if [ "$STRATEGY" = "" ]; then
	    fail "No download tool found in PATH. Please install aria2c, curl or wget."
fi

if [ -n "${PDB_AWS_DOWNLOAD}" ]; then
    if ! hasCommand aws; then
        fail "aws command not found in PATH. Please install the aws command line tool."
    fi
fi

downloadFile() {
    URL="$1"
    OUTPUT="$2"
    set +e
    for i in $STRATEGY; do
        case "$i" in
        ARIA)
            FILENAME=$(basename "${OUTPUT}")
            DIR=$(dirname "${OUTPUT}")
            aria2c --max-connection-per-server="$ARIA_NUM_CONN" --allow-overwrite=true -o "$FILENAME" -d "$DIR" "$URL" && set -e && return 0
            ;;
        CURL)
            curl -L -o "$OUTPUT" "$URL" && set -e && return 0
            ;;
        WGET)
            wget -O "$OUTPUT" "$URL" && set -e && return 0
            ;;
        esac
    done
    set -e
    fail "Could not download $URL to $OUTPUT"
}

if [ ! -f DOWNLOADS_READY ]; then
  if [ "${FAST_PREBUILT_DATABASES}" = "1" ]; then
    # new prebuilt GPU+CPU databases, that don't require calling tsv2exprofiledb
    downloadFile "https://opendata.mmseqs.org/colabfold/${UNIREF30DB}.db.tar.gz" "${UNIREF30DB}.tar.gz"
    downloadFile "https://opendata.mmseqs.org/colabfold/colabfold_envdb_202108.db.tar.gz" "colabfold_envdb_202108.tar.gz"
  else
    # old .tsv + tsv2exprofiledb databases
    downloadFile "https://opendata.mmseqs.org/colabfold/${UNIREF30DB}.tar.gz" "${UNIREF30DB}.tar.gz"
    downloadFile "https://opendata.mmseqs.org/colabfold/colabfold_envdb_202108.tar.gz" "colabfold_envdb_202108.tar.gz"
  fi
  if [ "${UNIREF30DB}" = "uniref30_2302" ]; then
    downloadFile "https://opendata.mmseqs.org/colabfold/uniref30_2302_newtaxonomy.tar.gz" "uniref30_2302_newtaxonomy.tar.gz"
  fi
  downloadFile "https://opendata.mmseqs.org/colabfold/pdb100_230517.fasta.gz" "pdb100_230517.fasta.gz"
  downloadFile "https://opendata.mmseqs.org/colabfold/pdb100_foldseek_230517.tar.gz" "pdb100_foldseek_230517.tar.gz"
  touch DOWNLOADS_READY
fi

if [ ! -f PDB_MMCIF_READY ]; then
  mkdir -p pdb/divided
  mkdir -p pdb/obsolete
  if [ -n "${PDB_AWS_DOWNLOAD}" ]; then
    aws s3 cp --no-sign-request --recursive s3://pdbsnapshots/${PDB_AWS_SNAPSHOT}/pub/pdb/data/structures/divided/mmCIF/ pdb/divided/
    aws s3 cp --no-sign-request --recursive s3://pdbsnapshots/${PDB_AWS_SNAPSHOT}/pub/pdb/data/structures/obsolete/mmCIF/ pdb/obsolete/
  fi
  rsync -rlpt -v -z --delete --port=${PDB_PORT} ${PDB_SERVER}/data/structures/divided/mmCIF/ pdb/divided
  rsync -rlpt -v -z --delete --port=${PDB_PORT} ${PDB_SERVER}/data/structures/obsolete/mmCIF/ pdb/obsolete
  touch PDB_MMCIF_READY
fi

if [ -n "$DOWNLOADS_ONLY" ]; then
  exit 0
fi


# Make MMseqs2 merge the databases to avoid spamming the folder with files
export MMSEQS_FORCE_MERGE=1

GPU_PAR=""
GPU_INDEX_PAR=""
if [ -n "${GPU}" ]; then
  GPU_PAR="--gpu 1"
  GPU_INDEX_PAR=" --split 1 --index-subset 2"

  if ! mmseqs --help | grep -q 'gpuserver'; then
    echo "The installed MMseqs2 has no GPU support, update to at least release 16"
    exit 1
  fi
fi

if [ ! -f UNIREF30_READY ]; then
  tar -xzvf "${UNIREF30DB}.tar.gz"
  if [ "${FAST_PREBUILT_DATABASES}" != "1" ]; then
    mmseqs tsv2exprofiledb "${UNIREF30DB}" "${UNIREF30DB}_db" ${GPU_PAR}
  fi
  if [ -z "$MMSEQS_NO_INDEX" ]; then
    mmseqs createindex "${UNIREF30DB}_db" tmp1 --remove-tmp-files 1 ${GPU_INDEX_PAR}
  fi

  # replace mapping and taxonomy with rebuilt versions, see:
  # https://github.com/sokrypton/ColabFold/wiki/MSA-Server-Database-History#2025-08-04-updated-uniref100_2302-taxonomypairing-files
  if [ -e "uniref30_2302_newtaxonomy.tar.gz" ]; then
    tar -xzvf "uniref30_2302_newtaxonomy.tar.gz"
  fi

  if [ -e ${UNIREF30DB}_db_mapping ]; then
    # create binary, mmap-able taxonomy mapping, saves a few seconds of load time during pairing
    TAXHEADER=$(od -An -N4 -t x4 "${UNIREF30DB}_db_mapping" | tr -d ' ')
    # check if the file is already binary, it has a binary-encoded header that spells TAXM if so
    if [ "${TAXHEADER}" != "0c170013" ]; then
      mmseqs createbintaxmapping "${UNIREF30DB}_db_mapping" "${UNIREF30DB}_db_mapping.bin"
      mv -f -- "${UNIREF30DB}_db_mapping.bin" "${UNIREF30DB}_db_mapping"
    fi
    ln -sf ${UNIREF30DB}_db_mapping ${UNIREF30DB}_db.idx_mapping
  fi
  if [ -e ${UNIREF30DB}_db_taxonomy ]; then
    ln -sf ${UNIREF30DB}_db_taxonomy ${UNIREF30DB}_db.idx_taxonomy
  fi
  touch UNIREF30_READY
fi

if [ ! -f COLABDB_READY ]; then
  tar -xzvf "colabfold_envdb_202108.tar.gz"
  if [ "${FAST_PREBUILT_DATABASES}" != "1" ]; then
    mmseqs tsv2exprofiledb "colabfold_envdb_202108" "colabfold_envdb_202108_db" ${GPU_PAR}
  fi
  # TODO: split memory value for createindex?
  if [ -z "$MMSEQS_NO_INDEX" ]; then
    mmseqs createindex "colabfold_envdb_202108_db" tmp2 --remove-tmp-files 1 ${GPU_INDEX_PAR}
  fi
  touch COLABDB_READY
fi

if [ ! -f PDB_READY ]; then
  # for consistency with the other prebuilt databases
  # make pdb also compatible with both gpu and cpu
  if [ -n "${GPU}" ] || [ "${FAST_PREBUILT_DATABASES}" = "1" ]; then
    mmseqs createdb pdb100_230517.fasta.gz pdb100_230517_tmp
    mmseqs makepaddedseqdb pdb100_230517_tmp pdb100_230517
    mmseqs rmdb pdb100_230517_tmp
  else
    mmseqs createdb pdb100_230517.fasta.gz pdb100_230517
  fi
  if [ -z "$MMSEQS_NO_INDEX" ]; then
    mmseqs createindex pdb100_230517 tmp3 --remove-tmp-files 1 ${GPU_INDEX_PAR}
  fi
  touch PDB_READY
fi

if [ ! -f PDB100_READY ]; then
  tar -xzvf pdb100_foldseek_230517.tar.gz pdb100_a3m.ffdata pdb100_a3m.ffindex
  touch PDB100_READY
fi
