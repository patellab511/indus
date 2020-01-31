#!/bin/bash

# Run tests on plumed patched with MDAnalysis++
# - Input
#    (1) plumed executable


#############
### Input ###
#############

min_num_args=1
if [ "$#" -lt $min_num_args ]; then
	echo "Error in test script: Must supply the PLUMED executable for testing"
	exit 1
fi
plumed_exe=$1

# Make sure the executable is valid
if [[ -z $( command -v $plumed_exe ) ]]; then
	echo "Error in test script: could not find plumed executable \"$plumed_exe\""
	exit 1
fi


#################
### Functions ###
#################



#############
### Tests ###
#############

# Location of this script
main_test_dir=$( realpath $( dirname $0 ) )


### Bias in a cylinder in bulk water ###

test_subdir="bias_ntilde_v/sphere/RESTRAINT"
test_dir="${main_test_dir}/${test_subdir}"
echo "Running INDUS test: ($test_subdir) ..."
cd $test_dir
./run_driver.sh $plumed_exe &> stdout.log
test_1=$( diff -q plumed.out ref/plumed.out )
test_2=$( diff -q forces.out ref/forces.out )
if [[ -z $test_1 && -z $test_2 ]]; then
	echo "PASSED"
else
	echo "FAILED"
fi
cd $main_test_dir
