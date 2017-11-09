# MATAM

*Mapping-Assisted Targeted-Assembly for Metagenomics*

## Getting Started

The recommended way of getting MATAM is through conda (see below). For __SSU rRNA__ assembly, run:

1. Getting and indexing default SSU rRNA reference database

 `index_default_ssu_rrna_db.py -d $DBDIR --max_memory 10000`

2. Running MATAM

    * SSU rRNA recovery only

      `matam_assembly.py -d $DBDIR/SILVA_128_SSURef_NR95 -i reads.fastq --cpu 4 --max_memory 10000 -v`

    * SSU rRNA recovery and taxonomic assignment

      `matam_assembly.py -d $DBDIR/SILVA_128_SSURef_NR95 -i reads.fastq --cpu 4 --max_memory 10000 -v --perform_taxonomic_assignment`

       The taxonomic assignment is done with [RDP classifier](https://rdp.cme.msu.edu/) and the training model used by default is "16srrna"



## Compiling MATAM from source code

1. Cloning MATAM repository

  `git clone https://github.com/bonsai-team/matam.git && cd matam`

2. Compiling MATAM and dependencies

  `./build.py`

3. (Optional) Getting and indexing default SSU rRNA reference database

  `./index_default_ssu_rrna_db.py --max_memory 10000`

4. Assembling

  `$MATAMDIR/bin/matam_assembly.py -i reads.fastq --cpu 4 --max_memory 10000 -v`

## Hardware requirements

We recommand running MATAM with at least 10Go of free RAM. You can try running MATAM with less RAM if --max\_memory is set to a lower value (eg. --max\_memory 4000 for 4Go).

Some steps of MATAM are highly paralelized. You can get a significant speed increase during these steps by setting the --cpu option to a higher value

## Dependencies

### Quick install

To install all of the needed dependencies except samtools, you can run the following command-lines in Debian-like distributions :
```bash
sudo apt-get update && sudo apt-get install curl git gcc g++ python3 python3-pip default-jdk automake make cmake ant libsparsehash-dev zlib1g-dev bzip2
sudo pip install numpy
```

Since the samtools package in current Ubuntu-like distributions is usualy a deprecated version (v0.1.19), you probably have to get a more recent version. We recommand getting samtools through bioconda (https://bioconda.github.io/)

### Full dependencies list

* **gcc v4.9.0 or superior**, (full C++11 support, \<regex\> included, and partial C++14 support)
* C++ libraries: rt, pthread, zlib
* Samtools v1.x or superior
* automake, make, **cmake v3.1 or superior**
* Python 3
* pip
* numpy
* Apache Ant
* Java SE 7 JDK. OpenJDK is ok (openjdk-7-jdk paquet on debian)
* bzip2
* google sparse hash library (libsparsehash-dev paquet on debian)

## MATAM in Docker
[![Docker Build Status](https://img.shields.io/docker/build/bonsaiteam/matam.svg)](https://hub.docker.com/r/bonsaiteam/matam/)

To run MATAM using docker, just run:

`docker run bonsaiteam/matam matam_assembly.py`

## MATAM with conda

[![Anaconda-Server Badge](https://anaconda.org/bonsai-team/matam/badges/installer/conda.svg)](https://conda.anaconda.org/bonsai-team)
[![Anaconda-Server Badge](https://anaconda.org/bonsai-team/matam/badges/version.svg)](https://anaconda.org/bonsai-team/matam)

A conda package is available here: https://anaconda.org/bonsai-team/matam

Before you begin, you should have installed Miniconda or Anaconda. See https://conda.io/docs/installation.html for more details.  
Then you will need to add the followings channels:
```
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda config --add channels bonsai-team
conda config --add channels salford_systems
```
Finally, matam can be installed with:
`conda install matam`

## Indexing a custom reference database

To run MATAM on a custom reference database, run:

`matam_db_preprocessing.py -i ref_db.fasta -d my_ref_db --cpu 4 --max_memory 10000 -v`

## Running example datasets

The following example datasets are provided:

### 16 bacterial species simulated dataset

* Running de-novo assembly

  `$MATAMDIR/bin/matam_assembly.py -i $MATAMDIR/examples/16sp_simulated_dataset/16sp.art_HS25_pe_100bp_50x.fq --cpu 4 --max_memory 10000 -v`

* Running assembly in validation mode (For developpers. Exonerate must be available in $PATH)

  `$MATAMDIR/bin/matam_assembly.py -i $MATAMDIR/examples/16sp_simulated_dataset/16sp.art_HS25_pe_100bp_50x.fq  --true_references $MATAMDIR/examples/16sp_simulated_dataset/16sp.fasta --true_ref_taxo $MATAMDIR/examples/16sp_simulated_dataset/16sp.taxo.tab --cpu 4 --max_memory 10000 --debug`

## Release versioning

MATAM releases will be following the Semantic Versioning 2.0.0 rules described here: http://semver.org/spec/v2.0.0.html
