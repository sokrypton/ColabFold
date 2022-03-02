#!/bin/sh -e
sudo systemctl stop mmseqs-server
rm -rf -- jobs
mkdir -p jobs
sudo systemctl start mmseqs-server
