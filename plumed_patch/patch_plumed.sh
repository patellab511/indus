#!/bin/bash

#############
### Input ###
#############

min_num_args=1
if [ "$#" -lt $min_num_args ]; then
	echo "Error - Must supply the location of the PLUMED build to patch."
	exit 1
fi

# Main PLUMED directory (contains src directory)
plumed_root_dir=$1


#############################
### Files and directories ###
#############################

# Location of MDAnalysis++ source code
patch_dir=$( realpath $( dirname $0 ) )
mda_dir=$( dirname $patch_dir )
src_dir="${mda_dir}/src"

# INDUS ++ objects needed for patch
obj_list_file="${patch_dir}/external_objects.list"

# Interface object
interface_obj="IndusInterface"
interface_dir="${patch_dir}/IndusInterface"


##################
### Copy files ###
##################

# Colvar directory in plumed source
colvar_dir="$plumed_root_dir/src/colvar"
if [ ! -d $colvar_dir ]; then
	echo "Couldnt' find plumed colvar directory \"$colvar_dir\" - exiting."
	exit
fi

# Copy interface object
if [ ! -d $interface_dir ]; then
	echo "Couldnt' find interface directory \"$interface_dir\" - exiting."
	exit
fi
cp -v ${interface_dir}/${interface_obj}* $colvar_dir

# Read list of Indus objects
if [ ! -f $obj_list_file ]; then
	echo "Couldnt' find object list \"$obj_list_file\" - exiting."
	exit
fi
#read -r -a objects <<< $( cat $obj_list_file )
objects=($( cat "$obj_list_file" ))


# Copy over Indus objects
echo "Copying the following ${#objects[@]} objects in addition to $interface_obj:"
for obj in ${objects[@]}; do
	echo "  ${obj}"
done

for obj in ${objects[@]}; do
	# Loop over extensions for header and source
	for ext in h cpp; do
		# Some "objects" are header-only
		file_to_copy="${src_dir}/${obj}.${ext}"
		if [ -f $file_to_copy ]; then
			cp -v $file_to_copy $colvar_dir
		fi
	done
done

echo "Done"
