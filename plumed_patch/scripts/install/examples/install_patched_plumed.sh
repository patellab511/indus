#!/bin/bash -eE

# Run this script **after** running the INDUS patching script

# Number of Makefile jobs to run simultaneously
num_processors=8


####################
### Set up paths ###
####################

# Where to install PLUMED
# - Once PLUMED is installed, if your system has environment modules support, you
#   can put the necessary variables in your paths by running the following command:
#      module load $plumed_install_dir/lib/plumed_mpi/modulefile
plumed_install_dir="${HOME}/programs/plumed/2.4.0/indus/mpi"

# Root directory of PLUMED repo
plumed_dir="plumed2"

# xdrfile installation
have_xdrfile=1
XDRFILE_DIR="${HOME}/programs/xdrfile/1.1.4"


#################
### Compilers ###
#################

module li

# Use MPI compiler wrappers
export CC=mpicc
export CXX=mpic++
have_mpi=1

# General flags and libraries
FLAGS="-fPIC -g"                              # basic flags
FLAGS="$FLAGS -O3 -ffast-math -march=native"  # optimizations
if [[ $have_mpi -eq 1 ]]; then
	FLAGS="$FLAGS -DMPI_ENABLED"  # MPI flag for Indus
fi
export CFLAGS="${FLAGS}"
export CXXFLAGS="${FLAGS} -std=c++11"
export LIBS="-lstdc++"

# OpenMP
have_openmp=1
export OPENMP_FLAG=""
if [[ $have_openmp -eq 1 ]]; then
	export OPENMP_FLAG="-fopenmp"  # for Intel compilers, this is usually '-qopenmp' instead
	export CFLAGS="$OPENMP_FLAG ${CFLAGS}"
	export CXXFLAGS="$OPENMP_FLAG ${CXXFLAGS}"
	export LDFLAGS="$OPENMP_FLAG $LDFLAGS"
fi

# xdrfile
if [[ $have_xdrfile -eq 1 ]]; then
	XDRFILE_FLAGS="-I${XDRFILE_DIR}/include/xdrfile"
	export CFLAGS="${CFLAGS} ${XDRFILE_FLAGS}"
	export CXXFLAGS="${CXXFLAGS} ${XDRFILE_FLAGS}"
	export LDFLAGS="${LDFLAGS} -L${XDRFILE_DIR}/lib"
	export LIBS="-lxdrfile $LIBS"
fi


#################
### Configure ###
#################

echo "DATE: $( date )"

# Checks
if [[ ! -d $plumed_dir ]]; then
	echo "Unable to find PLUMED repo at location: $plumed_dir"
	echo "Exiting."
	exit 1
fi
patch_files_dir="$plumed_dir/src/orderparameters"
if [[ ! -d $patch_files_dir ]]; then
	echo "Unable to find INDUS patch files. Make sure you patch PLUMED before using this script."
	echo "Exiting."
	exit 1
fi

cd $plumed_dir

# Configure from main PLUMED directory
echo "PHASE: CONFIGURING"
module list
./configure --prefix=${plumed_install_dir} \
            --enable-modules=crystallization \
            --program-suffix=_mpi

# Modify PLUMED Makefiles as needed
if [[ $have_openmp -eq 1 ]]; then
	# Need to append the OpenMP flag to Makefile recipes for some of the
	# executables that PLUMED makes, or you get a link error
	# - This is necessary for the version of PLUMED that I use: it may not apply to other versions
	targets=( '$(PLUMED_MAIN_SHARED):'  '$(PLUMED_MAIN_RUNTIME):' )
	makefile_to_edit="src/lib/Makefile"
	for target in "${targets[@]}"; do
		echo "${target}"
		sed -i "/${target}/{n;s/$/ $OPENMP_FLAG/}" $makefile_to_edit
	done
fi


###########################
### Compile and install ###
###########################

echo "PHASE: COMPILING"
make -j $num_processors

echo "PHASE: INSTALLING"
make install

echo "DONE"
echo "DATE: $( date )"
