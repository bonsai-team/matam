# MATAM

## Install

`git submodule update --init --recursive`

`./build.sh`

## Dependencies

### ovgraphbuild

* gcc v5.1 or superior (full C++14 support)
* libraries: rt, pthread, zlib
* Seqan v2.0.1 or superior (provided as a submodule of ovgraphbuild)

### SGA

* automake
* google sparse hash library (sparsehash paquet on debian)
* Bamtools library (http://github.com/pezmaster31/bamtools, provided as a submodule of matam)
* zlib
* (optional but suggested) the jemalloc memory allocator

### bamtools

* zlib
