#!/bin/bash
# Run the INDUS-in-PLUMED regression tests with 'plumed driver'.
#
# usage: run_tests.sh <plumed executable>
#
# Each test directory holds plumed.dat, indus.input and run_driver.sh, plus a ref/
# directory with the expected plumed.out and forces.out. Old outputs are removed
# before each run so a crashed driver cannot pass by leaving stale files behind.
# Exit status: 0 if every test passed, 1 otherwise.

if [[ $# -lt 1 ]]; then
	echo "Error: must supply the PLUMED executable for testing"
	exit 1
fi
plumed_exe=$1
if [[ -z $( command -v "$plumed_exe" ) ]]; then
	echo "Error: could not find plumed executable \"$plumed_exe\""
	exit 1
fi

main_test_dir=$( realpath "$( dirname "$0" )" )

tests=(
	"bias_ntilde_v/sphere/RESTRAINT"
)

num_failed=0
for test_subdir in "${tests[@]}"; do
	test_dir="${main_test_dir}/${test_subdir}"
	echo "Running INDUS test: ($test_subdir) ..."
	cd "$test_dir" || { echo "  FAILED (missing directory)"; num_failed=$((num_failed+1)); continue; }

	rm -f plumed.out forces.out
	if ! ./run_driver.sh "$plumed_exe" &> stdout.log; then
		echo "  FAILED (plumed driver exited with an error; see $test_dir/stdout.log)"
		num_failed=$((num_failed+1)); continue
	fi

	passed=1
	for f in plumed.out forces.out; do
		if [[ ! -f $f ]]; then
			echo "  FAILED (no $f produced)"; passed=0
		elif ! diff -q "$f" "ref/$f" > /dev/null; then
			echo "  FAILED ($f differs from ref/$f)"; passed=0
		fi
	done
	if [[ $passed -eq 1 ]]; then echo "  PASSED"; else num_failed=$((num_failed+1)); fi
	cd "$main_test_dir"
done

if [[ $num_failed -eq 0 ]]; then
	echo "All PLUMED tests passed"
	exit 0
else
	echo "$num_failed PLUMED test(s) FAILED"
	exit 1
fi
