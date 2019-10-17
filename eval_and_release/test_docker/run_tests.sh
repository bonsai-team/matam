#!/usr/bin/env bash

# THIS SCRIPT IS FOR DEBIAN OR REDHAT SYSTEM ONLY
#######################################################################

set -e; #Exit immediately if a command exits with a non-zero status
#set -o xtrace

if [ $# -eq 1 ] && [ "$1" != "advanced" ]
then
        echo "Usage : $0 [advanced]"
        exit 1
fi

# update package information & install packages
if [ -n "$(command -v apt)" ]
then
    apt update && apt install -y wget
elif [ -n "$(command -v yum)" ]
then
    yum update -y && yum install -y wget bzip2 ca-certificates
else
    echo 'Not a debian or redhat distribution ??'
    exit 1
fi

# install anaconda
wget https://repo.anaconda.com/archive/Anaconda2-2019.07-Linux-x86_64.sh -O ~/anaconda.sh && \
bash ~/anaconda.sh -b -p $HOME/anaconda && \
export PATH="$HOME/anaconda/bin:$PATH"

# conda env
conda create -n matam_env -y && source activate matam_env

# conda channels (the order can inpact the outcome)
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# install matam
# conda install /tmp/conda_package/matam-*tar.bz2
# conda install --use-local matam # install dependencies
conda install -c bonsai-team/label/dev matam -y

## run matam
#matam_assembly.py -i /tmp/16sp_simulated_dataset/16sp.art_HS25_pe_100bp_50x.fq -d /tmp/db/SILVA_128_SSURef_NR95 --cpu 4 --max_memory 10000 -v -o /tmp/matam_out --perform_taxonomic_assignment


## SIMPLE TESTS ##
#sortmerna -h  rc not properly set
echo 'Testing SGA'
sga --help
echo 'Testing VSEARCH'
vsearch --help
echo 'Testing SAMTOOLS'
samtools --help
echo 'Testing MATAM'
matam_assembly.py --help
index_default_ssu_rrna_db.py --help
matam_db_preprocessing.py --help
matam_compare_samples.py --help


## ADVANCED TESTS ##
if [ "$1" == "advanced" ]
then

# (avoid metaquast UnicodeDecodeError: 'ascii' codec can't decode byte ...)
if [ -n "$(command -v locale)" ]
then
    apt-get install -y locales && sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    locale-gen && dpkg-reconfigure --frontend=noninteractive locales && \
    update-locale LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

elif [ -n "$(command -v localedef)" ]
then
    localedef -v -c -i en_US -f UTF-8 en_US.UTF-8
else
    echo 'Not a debian or redhat distribution ??'
    exit 1
fi


# prepare test env
conda create -n test_env pytest statsmodels quast==5.0.2 -y
export PATH="$HOME/anaconda/envs/test_env/bin:$PATH"

# the functional test is expecting to find the database at a given location. Makes it accessible to the right location
ln -s /tmp/db/ $(realpath ~/anaconda/envs/matam_env/opt/matam-*)/db

# run tests
cd $(realpath ~/anaconda/envs/matam_env/opt/matam-*)/tests

echo 'Testing MATAM with pytest'
python -m pytest -rsx # this command must be executed in the tests dir
fi
