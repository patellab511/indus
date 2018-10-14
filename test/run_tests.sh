#!/bin/bash

# Runs *all* of the regtests

# Driver script used to perform each test
test_driver=$( realpath "run_calc.sh" )


#############
### Input ###
#############

# Default: MPI-enabled executable
if [ "$#" -lt 1 ]; then
	echo "Usage:  ./run_tests <executable>  [<echo_failed_diffs{0,1}> <is_mpi_enabled{0,1}> <mpi_cmd>]"
fi

# Project executable
indus=$(realpath $1)

# Whether to print the file diff when two files are different
echo_failed_diffs=$2
if [[ -z $echo_failed_diffs ]]; then
	echo_failed_diffs=0
fi

# Whether to test for MPI (0 or 1)
is_mpi_enabled=$3
if [[ -z $is_mpi_enabled ]]; then
	is_mpi_enabled=0
fi

# Command to use for mpirun (includes -np)
mpirun_cmd=$4
if [[ -z $is_mpi_enabled ]]; then
	mpirun_cmd="mpirun -np 2"
fi


#####################
### List of tests ###
#####################

# Test registry
declare -a regtests=(
	"indus/sphere_near_box_edge/bulk_water" \
	"indus/cylinder/bulk_water" \
	"indus/box/bulk_water" \
)

# MPI test registry
declare -a mpi_regtests=(\
	"indus_mpi/sphere_near_box_edge/bulk_water" \
)


#################
### Run tests ###
#################

start_dir=$PWD

# Register whether the test requires MPI
declare -a is_mpi_test
for (( i=0; $i<${#regtests[@]}; ++i )); do
	is_mpi_test+=(0)
done

# Include MPI regtests if MPI is enabled
if [[ $is_mpi_enabled -eq 1 ]]; then
	for mpi_regtest in ${mpi_regtests[@]}; do
		regtests+=($mpi_regtest)
		is_mpi_test+=(1)
	done
fi

# Perform each test
for (( i=0; $i<${#regtests[@]}; ++i )); do
	# Go to test directory
	regtest=${regtests[$i]}
	cd $regtest

	echo "Running test $regtest ..."
	result=$( $test_driver $indus $echo_failed_diffs ${is_mpi_test[$i]} "$mpirun_cmd" )
	if [[ -z "$result" ]]; then 
		echo "  PASSED" 
	else
		# Failed the test
		if [[ $echo_failed_diffs -eq 1 ]]; then
			# Print diffs
			echo "  $result"
		else
			echo "  FAILED"
		fi
	fi
	echo ""

	# !!! DANGER !!!
	# Uncomment one of these lines to update the reference output files
	### find . -maxdepth 1 -type f -exec cp -t "ref" {} +
	### cp * ref

	cd $start_dir
done
