#!/bin/bash
# Setup everything for using mmseqs locally

set -ex

cd "$1"

if ! [ -x "$(command -v aria2c)" ]; then
  alias download="wget"
else
  alias download="aria2c -x 16" # TODO: Is 16 the right default?
fi

if [ ! -f UNIREF30_READY ]; then
  download http://wwwuser.gwdg.de/~compbiol/colabfold/uniref30_2103.tar.gz
  tar xzvf uniref30_2103.tar.gz
  mmseqs tsv2exprofiledb uniref30_2103 uniref30_2103_db
  mmseqs createindex uniref30_2103_db tmp1 --remove-tmp-files 1
  touch UNIREF30_READY
fi

if [ ! -f COLABDB_READY ]; then
  download http://wwwuser.gwdg.de/~compbiol/colabfold/colabfold_envdb_202108.tar.gz
  tar xzvf colabfold_envdb_202108.tar.gz
  mmseqs tsv2exprofiledb colabfold_envdb_202108 colabfold_envdb_202108_db
  # TODO: split memory value for createindex?
  mmseqs createindex colabfold_envdb_202108_db tmp2 --remove-tmp-files 1
  touch COLABDB_READY
fi
