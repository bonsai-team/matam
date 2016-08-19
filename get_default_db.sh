#! /usr/bin/env bash

MATAMOG_DIR=$(pwd)

git lfs pull

echo "-- Extracting default ref db --"
cd $MATAMOG_DIR/db
tar jxvf SILVA_123_SSURef_rdNs_NR95.tar.bz2

echo "-- Indexing default ref db --"
cd $MATAMOG_DIR
$MATAMOG_DIR/scripts/index_ref_db.py -v
