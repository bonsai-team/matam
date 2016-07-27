#! /usr/bin/env bash

MATAMOG_DIR=$(pwd)

echo "-- Updating git submodules --"
git submodule update --init --recursive

echo "-- Compiling ComponentSearch --"
cd $MATAMOG_DIR/componentsearch
./compile.sh

echo "-- Compiling ovgraphbuild --"
mkdir $MATAMOG_DIR/ovgraphbuild/build
cd $MATAMOG_DIR/ovgraphbuild/build
cmake -G 'CodeBlocks - Unix Makefiles' ..
cmake ..
make

echo "-- Building SortMeRNA --"
cd $MATAMOG_DIR/sortmerna
./build.sh

echo "-- Compiling Sumaclust --"
cd $MATAMOG_DIR/sumaclust
make

echo "-- Compiling Bamtools lib (for SGA) --"
cd $MATAMOG_DIR/lib/bamtools
mkdir build
cd build/
cmake ..
make

echo "-- Compiling SGA --"
cd $MATAMOG_DIR/sga/src
./autogen.sh
./configure --with-bamtools=$MATAMOG_DIR/lib/bamtools
make

mkdir $MATAMOG_DIR/bin
cd $MATAMOG_DIR/bin

echo "-- Creating links into bin dir --"
ln -sf $MATAMOG_DIR/scripts/matam_*.py $MATAMOG_DIR/bin/.

echo "-- Extracting default ref db --"
cd $MATAMOG_DIR/db
tar jxvf SILVA_123_SSURef_rdNs_NR95.tar.bz2

echo "-- MATAM building complete --"
