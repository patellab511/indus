#!/bin/bash

plumed_exe=$1
plumed_input="plumed.dat"
xtc_file="../../../../test/sample_traj/bulk_water_278K/traj_water_278K_b1000_e1010.xtc"

if [[ $# -lt 1 ]]; then
	echo "Missing PLUMED driver"
	exit 1
elif [[ ! $( command -v $plumed_exe ) ]]; then
	echo "Invalid PLUMED driver: $plumed_exe"
	exit 1
fi

$plumed_exe driver \
	--plumed ${plumed_input} \
	--timestep 0.002 --trajectory-stride 500 \
	--ixtc ${xtc_file} \
	--dump-forces "forces.out" --dump-full-virial

if [[ $( ls bck.* ) ]]; then
	rm bck.*
fi
