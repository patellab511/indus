#!/bin/bash -eE

# Sample install script
# - Run from the top-level directory in the indus repo

# Use MPI compilers
export CC=mpicc
export CXX=mpic++

# xdrfile installation
export XDRFILE_DIR="${HOME}/programs/xdrfile/1.1.4"

# GPTL installation
export GPTL_DIR="${HOME}/programs/gptl/git/mpi"


############################################


### Configure ###

# Create/reset build directory
build_dir="$PWD/build"
if [[ -d $build_dir ]]; then
	rm -r $build_dir
fi
mkdir -p $build_dir

# Configure from build directory
cd $build_dir
cmake .. \
	-DCMAKE_INSTALL_PREFIX="${HOME}/programs/indus" \
	-DXDRFILE_DIR="$XDRFILE_DIR" \
	-DMPI_ENABLED=ON \
	-DOPENMP_ENABLED=ON \
	-DGPTL_ENABLED=ON \
	-DGPTL_DIR="$GPTL_DIR"


### Build ###

make -j 8
