# MATAM

*Mapping-Assisted Targeted-Assembly for Metagenomics* 

## Getting Started

Compiling MATAM and dependencies

`./build.sh`

Getting and indexing default reference database

`./get_default_db.sh`

Assembling

`./bin/matam_assembly.py -i reads.fastq --cpu 4 --max_memory 10000 -v`

## Hardware requirements

Recommanded free RAM is 10 Go, but MATAM should also work with less RAM if --max\_memory is set to a lower value (eg. --max\_memory 4000 for 4Go)

Some steps of MATAM are highly paralelized. You can get a significant speed increase during these steps by setting the --cpu option to a higher value

## Dependencies

### Quick install

For Debian-like distributions:

`sudo apt-get update && sudo apt-get install curl git gcc g++ default-jdk automake make cmake libsparsehash-dev zlib1g-dev && curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash && sudo apt-get update && sudo apt-get install git-lfs`

### General dependencies

* gcc v4.9.0 or superior (full C++11 support, \<regex\> included, and partial C++14 support)
* c++ libraries: rt, pthread, zlib
* automake, make, cmake
* google sparse hash library (sparsehash paquet on debian)
* git large file system (git-lfs) to get the supplied reference database
