#! /usr/bin/env bash

MATAMOG_DIR=$(pwd)

git submodule update --init --recursive

cd $MATAMOG_DIR/contigsearch
./compile.sh

mkdir $MATAMOG_DIR/ovgraphbuild/build
cd $MATAMOG_DIR/ovgraphbuild/build
cmake ..
make

cd $MATAMOG_DIR/sortmerna
./build.sh

cd $MATAMOG_DIR/sumaclust
make

mkdir $MATAMOG_DIR/bin
cd $MATAMOG_DIR/bin

ln -sf $MATAMOG_DIR/scripts/* $MATAMOG_DIR/bin/.
ln -sf $MATAMOG_DIR/sumaclust/sumaclust $MATAMOG_DIR/bin/.
ln -sf $MATAMOG_DIR/ovgraphbuild/bin/ovgraphbuild $MATAMOG_DIR/bin/.
ln -sf $MATAMOG_DIR/sortmerna/indexdb_rna $MATAMOG_DIR/bin/.
ln -sf $MATAMOG_DIR/sortmerna/sortmerna $MATAMOG_DIR/bin/.
