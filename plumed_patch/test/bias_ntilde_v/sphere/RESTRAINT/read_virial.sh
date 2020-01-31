#!/bin/bash

# Extract estimates of the virial from INDUS and PLUMED output files

# PLUMED stores 2x the virial in its output (a quirk)
# - Match lines beginning with the num_atoms token (16500) and print the following line,
#   which contains the components of the virial
plumed_file="forces.out"
sed -n -e '/^16500/{n;p}' $plumed_file > plumed_2x_virial.out

# Reference from INDUS code
# - Search for the token 'virial=' and print the rest of the line after it,
#   which contains the components of the virial
indus_file="${HOME}/source/indus/public/indus/test/indus/sphere_near_box_edge/bulk_water/forces_debug.out"
sed -n -e 's/virial= //p' $indus_file > indus_virial.out
