#!/bin/bash -ex
# Setup everything for using mmseqs locally
ARIA_NUM_CONN=8
WORKDIR="${1:-$(pwd)}"

PDB_SERVER="${2:-"rsync.wwpdb.org::ftp"}"
PDB_PORT="${3:-"33444"}"

cd "${WORKDIR}"

hasCommand () {
    command -v "$1" >/dev/null 2>&1
}

STRATEGY=""
if hasCommand aria2c; then STRATEGY="$STRATEGY ARIA"; fi
if hasCommand curl;   then STRATEGY="$STRATEGY CURL"; fi
if hasCommand wget;   then STRATEGY="$STRATEGY WGET"; fi
if [ "$STRATEGY" = "" ]; then
	    fail "No download tool found in PATH. Please install aria2c, curl or wget."
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

if [ ! -f UNIREF30_READY ]; then
  downloadFile "https://wwwuser.gwdg.de/~compbiol/colabfold/uniref30_2202.tar.gz" "uniref30_2202.tar.gz"
  tar xzvf "uniref30_2202.tar.gz"
  mmseqs tsv2exprofiledb "uniref30_2202" "uniref30_2202_db"
  mmseqs createindex "uniref30_2202_db" tmp1 --remove-tmp-files 1
  if [ -e uniref30_2202_db_mapping ]; then
    ln -sf uniref30_2202_db_mapping uniref30_2202_db.idx_mapping
  fi
  if [ -e uniref30_2202_db_taxonomy ]; then
    ln -sf uniref30_2202_db_taxonomy uniref30_2202_db.idx_taxonomy
  fi
  touch UNIREF30_READY
fi

if [ ! -f COLABDB_READY ]; then
  downloadFile "https://wwwuser.gwdg.de/~compbiol/colabfold/colabfold_envdb_202108.tar.gz" "colabfold_envdb_202108.tar.gz"
  tar xzvf "colabfold_envdb_202108.tar.gz"
  mmseqs tsv2exprofiledb "colabfold_envdb_202108" "colabfold_envdb_202108_db"
  # TODO: split memory value for createindex?
  mmseqs createindex "colabfold_envdb_202108_db" tmp2 --remove-tmp-files 1
  touch COLABDB_READY
fi

if [ ! -f PDB_READY ]; then
  downloadFile "https://wwwuser.gwdg.de/~compbiol/colabfold/pdb70_220313.fasta.gz" "pdb70_220313.fasta.gz"
  mmseqs createdb pdb70_220313.fasta.gz pdb70_220313
  mmseqs createindex pdb70_220313 tmp3 --remove-tmp-files 1
  touch PDB_READY
fi

if [ ! -f PDB70_READY ]; then
  downloadFile "https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb70_from_mmcif_220313.tar.gz" "pdb70_from_mmcif_220313.tar.gz"
  tar xzvf pdb70_from_mmcif_220313.tar.gz pdb70_a3m.ffdata pdb70_a3m.ffindex
  touch PDB70_READY
fi
if [ ! -f PDB_MMCIF_READY ]; then
  mkdir -p pdb/divided
  mkdir -p pdb/obsolete
  rsync -rlpt -v -z --delete --port=${PDB_PORT} ${PDB_SERVER}/data/structures/divided/mmCIF/ pdb/divided
  rsync -rlpt -v -z --delete --port=${PDB_PORT} ${PDB_SERVER}/data/structures/obsolete/mmCIF/ pdb/obsolete
  touch PDB_MMCIF_READY
fi