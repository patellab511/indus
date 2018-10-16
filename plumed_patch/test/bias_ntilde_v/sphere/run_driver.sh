#!/bin/bash

plumed_exe=$1
plumed_input="plumed.dat"
xtc_file="../../../../test/sample_traj/bulk_water_278K/traj_water_278K_b1000_e1010.xtc"

$plumed_exe driver \
	--plumed ${plumed_input} \
	--timestep 0.002 --trajectory-stride 500 \
	--ixtc ${xtc_file} \
	--dump-forces "forces.out" --dump-full-virial

if [ -f bck.* ]; then
	rm bck.*
fi
