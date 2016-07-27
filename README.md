# MATAM

*Mapping-Assisted Targeted-Assembly for Metagenomics* 

## Getting Started

Compiling MATAM and dependencies

`./build.sh`

Indexing default reference database

`./scripts/index_ref_db.py -v`

Assembling

`./bin/matam_assembly.py -i reads.fastq --cpu 4 -v`

## Dependencies

### ovgraphbuild

* gcc v4.9.0 or superior (full C++11 support, \<regex\> included, and partial C++14 support)
* libraries: rt, pthread, zlib
* Seqan v2.0.1 or superior (v2.2.0 is provided as a submodule of ovgraphbuild)

### SGA

* automake
* google sparse hash library (sparsehash paquet on debian)
* Bamtools library (http://github.com/pezmaster31/bamtools, provided as a submodule of matam)
* zlib
* (optional but suggested) the jemalloc memory allocator

### bamtools

* zlib
