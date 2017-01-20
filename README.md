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

Most MATAM steps will run with less than 10 Go of RAM and can run with less RAM if --max\_memory is set to a lower value (eg. --max\_memory 4000 for 4Go). However if the analysed sample is complex, the resulting graph can be huge and need to be stored in RAM by the graph_compaction step. In those cases you may need from 10s Go to 100s Go RAM.

Some steps of MATAM are highly paralelized. You can get a significant speed increase during these steps by setting the --cpu option to a higher value

## Dependencies

### Quick install

To install all of the needed depencies except samtools, you can run the following command-line in Debian-like distributions :

  `sudo apt-get update && sudo apt-get install curl git gcc g++ python3 default-jdk automake make cmake libsparsehash-dev zlib1g-dev bzip2`
  
Since the samtools package in current Ubuntu-like distributions is usualy a deprecated version (v0.1.19), you probably have to get a more recent version. We recommand getting samtools through bioconda (https://bioconda.github.io/)

### Full dependencies list

* **gcc v4.9.0 or superior**, (full C++11 support, \<regex\> included, and partial C++14 support)
* C++ libraries: rt, pthread, zlib
* Samtools v1.x or superior
* Python 3
* Java SE 7 JDK. OpenJDK is ok (openjdk-7-jdk paquet on debian)
* automake, make, cmake
* bzip2
* google sparse hash library (libsparsehash-dev paquet on debian)

## MATAM in Docker

To run MATAM using docker, just run:

`docker run ppericard/matam matam_assembly.py`

## Running example datasets

The following example datasets are provided:

### 16 bacterial species simulated dataset

* Running de-novo assembly

  `$MATAMDIR/bin/matam_assembly.py -i $MATAMDIR/examples/16sp_simulated_dataset/16sp.art_HS25_pe_100bp_50x.fq --cpu 4 --max_memory 10000 -v`
  
* Running assembly in validation mode (For developpers. Exonerate must be available in $PATH)

  `$MATAMDIR/bin/matam_assembly.py -i $MATAMDIR/examples/16sp_simulated_dataset/16sp.art_HS25_pe_100bp_50x.fq  --true_references $MATAMDIR/examples/16sp_simulated_dataset/16sp.fasta --true_ref_taxo $MATAMDIR/examples/16sp_simulated_dataset/16sp.taxo.tab --cpu 4 --max_memory 10000 -v`
