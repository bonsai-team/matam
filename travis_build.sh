#!/bin/bash

# An activated env is needed for $CONDA_PREFIX
#CC=${CONDA_PREFIX}/bin/gcc
#CXX=${CONDA_PREFIX}/bin/g++

# Build componentsearch
cd componentsearch && make && cd -

# Build ovgraphbuild
build_dir=ovgraphbuild/build
mkdir -p  $build_dir && cd $build_dir && cmake .. -G"CodeBlocks - Unix Makefiles" && make && cd -