#! /usr/bin/env bash

MATAMOG_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $MATAMOG_DIR

git lfs pull

echo "\n-- Extracting default ref db --"
cd $MATAMOG_DIR/db
tar jxvf SILVA_123_SSURef_rdNs_NR95.tar.bz2

echo "\n-- Indexing default ref db --"
cd $MATAMOG_DIR
$MATAMOG_DIR/scripts/index_ref_db.py -v
