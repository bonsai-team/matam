############################################################
# Dockerfile to build MATAM container images
# Based on Debian
############################################################

# Set the base image to Debian
FROM debian

# File Author / Maintainer
MAINTAINER Pierre Pericard

################## BEGIN INSTALLATION ######################
# Install MATAM following the instructions from Github repository
# Ref: https://github.com/ppericard/matam

# Install dependencies
RUN apt-get update && apt-get install --no-install-recommends -y \
    curl \
    git \
    gcc \
    g++ \
    python3 \
    default-jdk \
    automake \
    make \
    cmake \
    libsparsehash-dev \
    zlib1g-dev \
    bzip2

# Install git lfs repository and paquet
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash \
&& apt-get install -y git-lfs

# Clean apt cache
RUN rm -rf /var/lib/apt/lists/*

# Build MATAM
RUN ./build.py
#RUN ./index_default_ssu_rrna_db.py --max_memory 4000

# Set PATH
ENV PATH /bin:$PATH

##################### INSTALLATION END #####################
