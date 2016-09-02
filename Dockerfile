############################################################
# Dockerfile to build MATAM container images
# Based on Ubuntu
############################################################

# Set the base image to Ubuntu
FROM ubuntu

# File Author / Maintainer
MAINTAINER Pierre Pericard

# Update the repository sources list
RUN apt-get update

################## BEGIN INSTALLATION ######################
# Install MATAM following the instructions from Github repository
# Ref: https://github.com/ppericard/matam

# Install dependencies
RUN apt-get install -y curl
RUN apt-get install -y git

#RUN apt-get install -y gcc-4.9 g++-4.9
#RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 60 --slave /usr/bin/g++ g++ /usr/bin/g++-4.9
RUN apt-get install -y gcc g++

RUN apt-get install -y default-jdk
RUN apt-get install -y automake make
RUN apt-get install -y cmake
RUN apt-get install -y libsparsehash-dev
RUN apt-get install -y zlib1g-dev

RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt-get install git-lfs

RUN apt-get install bzip2

# Tests
RUN gcc --version
RUN g++ --version

# Get MATAM
RUN git clone https://github.com/ppericard/matam.git

# Build MATAM
WORKDIR /matam
RUN ./build.py
RUN ./index_default_ssu_rrna_db.py --max_memory 4000

##################### INSTALLATION END #####################
