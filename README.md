# MATAM

*Mapping-Assisted Targeted-Assembly for Metagenomics* 

## Getting Started

1. Cloning MATAM repository

  `git clone https://github.com/ppericard/matam.git`

2. Compiling MATAM and dependencies

  `./build.py`

3. (Optional) Getting and indexing default SSU rRNA reference database

  `./index_default_ssu_rrna_db.py --max_memory 10000`

4. Assembling

  `$MATAMDIR/bin/matam_assembly.py -i reads.fastq --cpu 4 --max_memory 10000 -v`

## Hardware requirements

Recommanded free RAM when using provided default reference database is 10 Go, but MATAM should also work with less RAM if --max\_memory is set to a lower value (eg. --max\_memory 4000 for 4Go)

Some steps of MATAM are highly paralelized. You can get a significant speed increase during these steps by setting the --cpu option to a higher value

## Dependencies

### Quick install

To install all the needed depencies, you can run the following command-line in Debian-like distributions :

  `sudo apt-get update && sudo apt-get install curl git gcc g++ default-jdk automake make cmake libsparsehash-dev zlib1g-dev && curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash && sudo apt-get update && sudo apt-get install git-lfs`

### Full dependencies list

* **gcc v4.9.0 or superior** (full C++11 support, \<regex\> included, and partial C++14 support)
* C++ libraries: rt, pthread, zlib
* Java SE 7 JDK. OpenJDK is ok (openjdk-7-jdk paquet on debian)
* automake, make, cmake
* google sparse hash library (libsparsehash-dev paquet on debian)
* git large file system (git-lfs), if you want to get the supplied reference database

## Running example datasets

The following example datasets are provided:

### 16 bacterial species simulated dataset

* Running assembly

  `$MATAMDIR/bin/matam_assembly.py -i $MATAMDIR/db/16sp.art_HS25_pe_100bp_50x.fq --cpu 4 --max_memory 10000 -v`
  
* Running validation mode (Exonerate must be available in $PATH)

  `$MATAMDIR/bin/matam_assembly.py -i $MATAMDIR/db/16sp.art_HS25_pe_100bp_50x.fq  --true_references $MATAMDIR/db/16sp.fasta --true_ref_taxo $MATAMDIR/db/16sp.taxo.tab --cpu 4 --max_memory 10000 -v`
