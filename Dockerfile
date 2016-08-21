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
RUN apt-get install gcc-4.9
RUN apt-get install libsparsehash-dev

# Build MATAM
RUN ./build.sh

##################### INSTALLATION END #####################
