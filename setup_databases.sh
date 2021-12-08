#!/bin/bash -ex
# Setup everything for using mmseqs locally
ARIA_NUM_CONN=8
WORKDIR="${1:-$(pwd)}"

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
            aria2c --max-connection-per-server="$ARIA_NUM_CONN" --allow-overwrite=true -o "$FILENAME" -d "$DIR" "$URL" && return 0
            ;;
        CURL)
            curl -L -o "$OUTPUT" "$URL" && return 0
            ;;
        WGET)
            wget -O "$OUTPUT" "$URL" && return 0
            ;;
        esac
    done
    set -e
    fail "Could not download $URL to $OUTPUT"
}

if [ ! -f UNIREF30_READY ]; then
  downloadFile "http://wwwuser.gwdg.de/~compbiol/colabfold/uniref30_2103.tar.gz" "uniref30_2103.tar.gz"
  tar xzvf "uniref30_2103.tar.gz"
  mmseqs tsv2exprofiledb "uniref30_2103" "uniref30_2103_db"
  mmseqs createindex "uniref30_2103_db" tmp1 --remove-tmp-files 1
  touch UNIREF30_READY
fi

if [ ! -f COLABDB_READY ]; then
  downloadFile "http://wwwuser.gwdg.de/~compbiol/colabfold/colabfold_envdb_202108.tar.gz" "colabfold_envdb_202108.tar.gz"
  tar xzvf "colabfold_envdb_202108.tar.gz"
  mmseqs tsv2exprofiledb "colabfold_envdb_202108" "colabfold_envdb_202108_db"
  # TODO: split memory value for createindex?
  mmseqs createindex "colabfold_envdb_202108_db" tmp2 --remove-tmp-files 1
  touch COLABDB_READY
fi
