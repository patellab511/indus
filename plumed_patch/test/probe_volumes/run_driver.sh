#!/bin/bash
# Run 'plumed driver' for one probe-volume test. Called from inside the test directory.
#   usage: ../run_driver.sh <plumed executable>
# The test directory provides plumed.dat, indus.input and a file 'trajectory' holding the
# xtc path (relative to the test directory). Writes plumed.out and forces.out.
plumed_exe=$1
if [[ $# -lt 1 || -z $( command -v "$plumed_exe" ) ]]; then
	echo "usage: $0 <plumed executable>"; exit 1
fi
xtc=$( cat trajectory )
if [[ ! -f $xtc ]]; then
	echo "trajectory not found: $xtc"; exit 1
fi
rm -f plumed.out forces.out bck.*
"$plumed_exe" driver --plumed plumed.dat --ixtc "$xtc" --timestep 0.002 --trajectory-stride 500 \
	--dump-forces forces.out --dump-forces-fmt %.5f
status=$?
rm -f bck.*
exit $status
