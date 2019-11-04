#!/bin/bash -eE

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

# Template module directory
op_module="${patch_dir}/src/orderparameters"
module_name=$( basename $op_module )

# A list of file names (without extensions) to ignore
# - The vast majority are included, so only ignore particular ones
exclude_files_list="${patch_dir}/exclude_files.list"


##################
### Copy files ###
##################

# Colvar directory in plumed source
plumed_src="$plumed_root_dir/src"
if [[ ! -d $plumed_src ]]; then
	echo "Couldn't find plumed source directory, \"$plumed_src\" - exiting."
	exit 1
fi

# Copy module files
plumed_module="${plumed_src}/${module_name}"
if [[ ! -d $op_module ]]; then
	echo "Couldn't find directory with module files, \"$op_module\" - exiting."
	exit 1
fi
if [[ -d $plumed_module ]]; then
	echo "Removing old module directory from PLUMED repo"
	rm -rf $plumed_module
fi
cp -rv $op_module $plumed_module

# Read list of OrderParameters excluded_files
if [[ ! -f $exclude_files_list ]]; then
	echo "Couldn't find list of files to exclude, \"$exclude_files_list\" - exiting."
	exit 1
fi
excluded_files=($( cat "$exclude_files_list" ))

# Get a list of all src files
IFS=$'\n' src_files=($( find ${src_dir} -type f ))

# Copy non-excluded files to module directory in PLUMED source
for src_file in ${src_files[@]}; do
	src_file_name=$( basename $src_file )  # trim off path
	src_name=${src_file_name%.*}           # trim off file extension

	if [[ ! " ${excluded_files[@]} " =~ " ${src_name} " ]]; then
		cp -v $src_file $plumed_module
	else
		echo "Ignoring file $src_file_name"
	fi
done

echo "Done"
