#!/bin/bash

# - Solution for parsing arguments adapted from user 'John' on StackOverflow
#   - Source: https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash  


#############
### Input ###
#############

start_dir=$PWD

### Defaults ###

program=""             # Project executable
test_dir="test"        # Directory in repo with test files
echo_failed_diffs=1    # whether to echo failed 'diffs'

# Parallelization
num_ranks=1
num_threads=1


### Parse Input ###

# Save positional args
POSITIONAL=()

while [[ $# -gt 0 ]]; do
	key="$1"
	case $key in
		# Options with a value
		-e|--exe|-program|--program)
			program=$( realpath $2 )
			shift; shift
			;;
		-n|-np|--np)
			num_ranks="$2"
			shift; shift
			;;
		-t|-nt|--nt)
			num_threads="$2"
			shift; shift
			;;
		-d|--dir|--test_dir)
			test_dir="$2"
			shift; shift
			;;
		# Flags
		--diff)
			echo_failed_diffs=1
			shift;
			;;
		#--default)
		#	DEFAULT=YES
		#	shift
		#	;;
		*)
			# Save positional args
			POSITIONAL+=("$1")
			shift
			;;
	esac
done

set -- "${POSITIONAL[@]}"  # restore positional parameters


### Check Input ###

if [[ -z $program ]]; then
	echo "FAILED: no executable provided"
	exit 1
elif [[ ! $( command -v $program ) ]]; then
	echo "FAILED: couldn't find executable"
	echo "(input: $program)"
	exit 1
fi

if [[ ! -d $test_dir ]]; then
	echo "FAILED: couldn't find test directory"
	echo "(input: relpath $test_dir, fullpath $start_dir/$test_dir)"
	exit 1
fi

is_parallel=0
if [[ $num_threads -gt 1 || $num_threads -gt 1 ]]; then
	is_parallel=1
	export MPI_NUM_RANKS=$num_ranks
	export OMP_NUM_THREADS=$num_threads
	program="mpirun --use-hwthread-cpus -np $MPI_NUM_RANKS $program"
fi


###############
### Logging ###
###############

# DEBUG: print info
echo "DEBUG: Inside test script: $0"
echo "  PWD:                $PWD"
echo "  program:            $program"
echo "  test_dir:           $test_dir"
echo "  echo_failed_diffs:  $echo_failed_diffs"
echo "  is_parallel:        $is_parallel"
if [[ $is_parallel -eq 1 ]]; then
	echo "  MPI_NUM_RANKS:      $MPI_NUM_RANKS"
	echo "  OMP_NUM_THREADS:    $OMP_NUM_THREADS"
fi
echo ""


################
### Run test ###
################

# Go to testing directory
cd $test_dir

# !!! DANGER !!!
# TODO: How to update this procedure?
# Uncomment one of these lines to update the reference output files
### find . -maxdepth 1 -type f -exec cp -t "ref" {} +
### cp * ref

# Run test
echo "Running test  ..."
passed_test=1
$program "indus.input"
if [[ $? -ne 0 ]]; then
	echo "  FAILED: Program exited with error"
	passed_test=0
fi

if [[ $passed_test -eq 1 ]]; then
	# Provided that the program exited successfully, check output
	read -d '' -ra output_files < "output_files_to_check.input"
	for output_file in ${output_files[@]}; do
		# Check that the files exist
		if [[ ! -f $output_file ]]; then
			echo "  FAILED: output file \"$output_file\" not found"
			passed_test=0
			continue
		elif [[ ! -f ref/$output_file ]]; then
			echo "  FAILED: reference output file \"ref/$output_file\" not found"
			passed_test=0
			continue
		fi

		# Compare files
		diff_out=$( diff $output_file ref/$output_file )
		if [[ ! -z $diff_out ]]; then
			passed_test=0
			echo "  FAILED"

			if [[ $echo_failed_diffs -eq 1 ]]; then
				echo "  diff $output_file ref/$output_file"
				diff $output_file ref/$output_file
				echo ""
			fi
		fi
	done

	# If all is well by this point, the test has been passed
	if [[ $passed_test -eq 1 ]]; then
		echo "  PASSED"
	else
		echo "  FAILED"
	fi
fi

# Convert from 'passed_test' flag to return code
#   ret=0: passed
#   ret=1: failed
ret=0
if [[ $passed_test -ne 1 ]]; then
	ret=1
fi

exit $ret
