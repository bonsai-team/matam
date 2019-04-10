############################################################
# Dockerfile to build MATAM container images
# Based on Debian
############################################################

# Set the base image to Debian
FROM debian

# File Author / Maintainer
MAINTAINER Bonsai Team

################## BEGIN INSTALLATION ######################
# Install MATAM following the instructions from Github repository
# Ref: https://github.com/bonsai-team/matam

# Install dependencies
RUN apt-get update && apt-get install --no-install-recommends -y \
    git \
    gcc \
    g++ \
    wget \
    bash \
    default-jdk \
    automake \
    make \
    cmake \
    libsparsehash-dev \
    zlib1g-dev \
    bzip2 \
    ant
    # Python 3 now comes from conda

# Install Conda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b
ENV PATH /root/miniconda3/bin:$PATH

# Install Bioconda and samtools

# Keep channels in this order. See https://github.com/bioconda/bioconda-recipes/issues/12100
# for more details.
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda install --update-deps -y samtools
RUN conda install --update-deps -y numpy

# Clean apt cache
RUN rm -rf /var/lib/apt/lists/*

# Cloning MATAM
RUN git clone https://github.com/bonsai-team/matam.git

# Build MATAM
WORKDIR /matam
RUN ./build.py
#RUN ./index_default_ssu_rrna_db.py --max_memory 4000

# Add index_default_ssu_rrna_db.py in the PATH
# ln -s /matam/index_default_ssu_rrna_db.py /matam/bin/index_default_ssu_rrna_db.py

# Set PATH
ENV PATH /matam/bin:$PATH

##################### INSTALLATION END #####################
