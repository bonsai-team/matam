#! /usr/bin/env bash

git submodule update --init --recursive

cd contigsearch
./compile.sh
cd ..

mkdir ovgraphbuild/build
cd ovgraphbuild/build
cmake ..
make
cd ../../

cd sortmerna
./build.sh
cd ..

cd sumaclust
make
cd ../

mkdir bin

ln -sf scripts/* bin/.
ln -sf sumaclust/sumaclust bin/.
ln -sf ovgraphbuild/bin/ovgraphbuild bin/.
ln -sf sortmerna/indexdb_rna bin/.
ln -sf sortmerna/sortmerna bin/.
