#!/bin/bash

#############
### Input ###
#############

# Project executable
md_analysis=$1
if [[ -z $md_analysis ]]; then
	# Try default path
	md_analysis="${HOME}/source/MDAnalysis++/build/bin/MDAnalysis++"
	if [[ ! $( command -v $md_analysis ) ]]; then
		echo "FAILED (couldn't find MDAnalysis++ executable)"
		exit
	fi
fi


################
### Run test ###
################

# Trap runtime errors
# - e.g. segfaults, uncaught exceptions, etc.
# - Note: set -eE  -->  set -o errexit -o errtrace
set -eE 
trap "echo \"FAILED (Program exited with error)\"" ERR

export OMP_NUM_THREADS=1
$md_analysis test &> stdout.log

# Disable error trapping
trap '' ERR
set +eE
