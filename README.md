MATAM [![Build Status](https://travis-ci.org/bonsai-team/matam.svg?branch=master)](https://travis-ci.org/bonsai-team/matam)
=========
*Mapping-Assisted Targeted-Assembly for Metagenomics*

MATAM is a software dedicated to the fast and accurate targeted assembly of short reads sequenced from a genomic marker of interest.
It has been applied to the assembly of 16S rRNA markers and is validated on simulated, synthetic and genuine metagenomes.

The related article of this method is available [here](https://doi.org/10.1093/bioinformatics/btx644).


# Table of contents
* [Hardware requirements](#hardware-requirements)
* [Installation](#installation)
  * [MATAM with conda](#matam-with-conda) (recommanded)
  * [MATAM in Docker](#matam-in-docker)
  * [MATAM from source code](#matam-from-source-code)
    * [Full dependencies list](#full-dependencies-list)
    * [Install dependencies](#install-dependencies)
    * [Compile MATAM](#compile-matam)
* [Run MATAM](#run-matam)
  * [Database preparation](#database-preparation)
    * [Provided database](#provided-database)
    * [Custom database](#custom-database)
  * [De-novo assembly](#de-novo-assembly)
  * [Example with default database and provided dataset](#example-with-default-database-and-provided-dataset)
* [Samples Comparaison](#samples-comparaison)
* [Release versioning](#release-versioning)


# <a id="hardware-requirements"></a>Hardware requirements

We recommand running MATAM with at least 10Go of free RAM. You can try running MATAM with less RAM if --max\_memory is set to a lower value (eg. --max\_memory 4000 for 4Go).

Some steps of MATAM are highly paralelized. You can get a significant speed increase during these steps by setting the --cpu option to a higher value


# <a id="installation"></a>Installation

## <a id="matam-with-conda"></a>MATAM with conda

[![Anaconda-Server Badge](https://anaconda.org/bonsai-team/matam/badges/installer/conda.svg)](https://conda.anaconda.org/bonsai-team)
[![Anaconda-Server Badge](https://anaconda.org/bonsai-team/matam/badges/version.svg)](https://anaconda.org/bonsai-team/matam)

Before you begin, you should have installed Miniconda or Anaconda. See https://conda.io/docs/installation.html for more details.  
Then you will need to add the followings channels:
```
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda config --add channels bonsai-team
```
Finally, matam can be installed with:

`conda install matam`

All the commands used in this README will be available in your PATH.


## <a id="matam-in-docker"></a>MATAM in Docker

[![Docker Build Status](https://img.shields.io/docker/build/bonsaiteam/matam.svg)](https://hub.docker.com/r/bonsaiteam/matam/)

To retrieve the docker image, run the following command:

`docker pull bonsaiteam/matam`

Then all the commands used in this README will be available as:

`docker run -v host_directory:/workdir bonsaiteam/matam CMD`

Noticed that you have to specify a docker volume to share data between the host and the container and use this workdir for your analysis. Otherwise your data will be lost when exiting the container.

Finally, if you prefer an interactive session with the container, run:

`docker run -it bonsaiteam/matam`


## <a id="matam-from-source-code"></a>MATAM from source code

### <a id="full-dependencies-list"></a>Full dependencies list

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


### <a id="install-dependencies"></a>Install dependencies

To install all of the needed dependencies except samtools, you can run the following command-lines in Debian-like distributions :
```bash
sudo apt-get update && sudo apt-get install curl git gcc g++ python3 python3-pip default-jdk automake make cmake ant libsparsehash-dev zlib1g-dev bzip2
sudo pip install numpy
```

Since the samtools package in current Ubuntu-like distributions is usualy a deprecated version (v0.1.19), you probably have to get a more recent version. We recommand getting samtools through bioconda (https://bioconda.github.io/)


### <a id="compile-matam"></a>Compile MATAM

1. Cloning MATAM repository

  `git clone https://github.com/bonsai-team/matam.git && cd matam`

2. Compile MATAM and dependencies

  `./build.py`

3.  Update your PATH to make MATAM's commands available:

  ```bash
  echo 'export PATH=$MATAMDIR/bin:$PATH' >> ~/.profile
  source ~/.profile
  ```


# <a id="run-matam"></a>Run MATAM

## <a id="database-preparation"></a>1. Database preparation (clusterization & indexation)

### <a id="provided-database"></a>Provided database

By default, MATAM provides a SSU rRNA reference database where the clusterisation step has already been done (i.e. the sequences sharing 95% of identity have been clusterised with [Sumaclust](https://git.metabarcoding.org/obitools/sumaclust/wikis/home)).  
The  [FASTA](https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz) file used for this database comes from [SILVA 128 release](https://www.arb-silva.de/documentation/release-128/).

To use the default SSU rRNA reference database, run the following command:

`index_default_ssu_rrna_db.py -d $DBDIR --max_memory 10000`

### <a id="custom-database"></a>Custom database

If the provided database does not fulfill your needs, you can prepare a custom database of your own by running the following command:

`matam_db_preprocessing.py -i ref_db.fasta -d $DBDIR --cpu 4 --max_memory 10000 -v`

## <a id="de-novo-assembly"></a>2. De-novo assembly

When your database is ready, then you will be able to reconstruct your markers:
* SSU rRNA recovery only

  `matam_assembly.py -d $DBDIR/prefix -i reads.fastq --cpu 4 --max_memory 10000 -v`

* SSU rRNA recovery and taxonomic assignment

  `matam_assembly.py -d $DBDIR/prefix -i reads.fastq --cpu 4 --max_memory 10000 -v --perform_taxonomic_assignment`  
    The taxonomic assignment is done with [RDP classifier](https://rdp.cme.msu.edu/) and the training model used by default is "16srrna"


## <a id="example-with-default-database-and-provided-dataset"></a>Example with default database and provided dataset

1. Getting and indexing default SSU rRNA reference database

  `index_default_ssu_rrna_db.py -d $DBDIR --max_memory 10000`

2. Retrieve the example dataset: [16 bacterial species simulated dataset](examples/16sp_simulated_dataset/16sp.art_HS25_pe_100bp_50x.fq)

  `wget https://raw.githubusercontent.com/bonsai-team/matam/master/examples/16sp_simulated_dataset/16sp.art_HS25_pe_100bp_50x.fq`

3. De-novo assembly

  `matam_assembly.py -d $DBDIR/SILVA_128_SSURef_NR95 -i 16sp.art_HS25_pe_100bp_50x.fq --cpu 4 --max_memory 10000 -v --perform_taxonomic_assignment`

For developpers, the commande to run assembly in validation mode (Exonerate must be available in $PATH) is:

`$MATAMDIR/bin/matam_assembly.py -i 16sp.art_HS25_pe_100bp_50x.fq  --true_references $MATAMDIR/examples/16sp_simulated_dataset/16sp.fasta --true_ref_taxo $MATAMDIR/examples/16sp_simulated_dataset/16sp.taxo.tab --cpu 4 --max_memory 10000 --debug`

# <a id="samples-comparaison"></a>Samples Comparaison

An utilatary script is provided to compare the abundance of different samples. (Available only if the `--perform_taxonomic_assignment` was used when running MATAM).

`matam_compare_samples.py -s samples_to_compare.tsv -t contingency_table.tsv -c comparaison_table.tsv`

The `samples_to_compare.tsv` file is a tabulated file listing the FASTA & RDP files of each sample to compare (see example below).  
The `contingency_table.tsv` file will report the abundance for each sequence.  
The `comparaison_table.tsv` file will report a comparaison by "taxonomic path" of the abundance for the samples.


Example
```
sample1 <tab> $WORKDIR/matam_sample1/final_assembly.fa <tab> $WORKDIR/matam_sample1/rdp.tab
sample2 <tab> $WORKDIR/matam_sample2/final_assembly.fa <tab> $WORKDIR/matam_sample2/rdp.tab
```

The first column is the ID of the sample and it must be unique among the file.  



# <a id="release-versioning"></a>Release versioning

MATAM releases will be following the Semantic Versioning 2.0.0 rules described here: http://semver.org/spec/v2.0.0.html
