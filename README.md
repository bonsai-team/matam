# MATAM

*Mapping-Assisted Targeted-Assembly for Metagenomics* 

## Getting Started

`./build.sh`

`matam.py -s 0 1 2 3 4 5 6 7 8 -i reads.fastq -d SILVA_123_SSURef.fasta`

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
