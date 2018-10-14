#!/bin/bash

# Use this script for inidividual tests.

# Prints "FAILED" if one of the output files doesn't match the reference
# file located in './ref', else doesn't print anything

#############
### Input ###
#############

# Project executable
indus=$1
if [[ -z $indus ]]; then
	# Try default path
	indus="${HOME}/source/indus/build/bin/indus"
	if [[ ! $( command -v $indus ) ]]; then
		echo "FAILED (couldn't find indus executable)"
		exit
	fi
fi

# Whether to print full diffs when a comparison fails
echo_failed_diffs=$2
if [[ -z $echo_failed_diffs ]]; then
	echo_failed_diffs=1
fi

# Whether this test requires MPI
is_mpi_test=$3
if [[ -z $is_mpi_test ]]; then
	# Default: non-MPI test
	is_mpi_test=0
fi

# mpirun command (includes -np <#>) for MPI tests
mpirun_cmd=$4
if [[ -z "$mpirun_cmd" ]]; then
	# Default MPI commmand: automatically select number of ranks
	mpirun_cmd="mpirun -np 0"
fi


################
### Run test ###
################

# Trap runtime errors
# - e.g. segfaults, uncaught exceptions, etc.
# - Note: set -eE  -->  set -o errexit -o errtrace
set -eE 
trap "echo \"FAILED (Program exited with error)\"" ERR

# Run program
if [[ $is_mpi_test -eq 1 ]]; then
	$mpirun_cmd $indus indus.input &> stdout.log
else
	$indus indus.input &> stdout.log
fi

# Disable error trapping
trap '' ERR
set +eE

# Output files to check (read into array)
read -d '' -ra output_files < "output_files_to_check.input"

# Check output
for output_file in ${output_files[@]}; do
	diff_out=$( diff $output_file ref/$output_file )
	if [[ ! -z $diff_out ]]; then
		echo "FAILED"

		if [[ $echo_failed_diffs -eq 1 ]]; then
			echo "  diff $output_file ref/$output_file"
			diff $output_file ref/$output_file
			echo ""
		fi
	fi
done
