#!/bin/bash
#SBATCH --job-name=indus-stack
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=4                 # the parallel tests use TEST_RANKS=4 MPI ranks
#SBATCH --cpus-per-task=8          # ntasks x cpus-per-task = cores used for the build
#SBATCH --time=01:55:00
#SBATCH --output=indus-stack.%j.out
# Build and test the INDUS + PLUMED + GROMACS stack inside a SLURM job (CPU only here).
# Before submitting, on the login node (compute nodes usually have no internet):
#   scripts/install/install_indus_stack.sh --no-gpu --root $STACK_ROOT --only fetch
# Adjust the module names, STACK_ROOT and GMX_SIMD for your site.
STACK_ROOT=/scratch/$USER/indus-stack

source /etc/profile.d/modules.sh
module purge
module load gcc-9.2.0/openmpi/4.1.0          # a C++17 compiler and the site's MPI
export MPICC=$(command -v mpicc) MPICXX=$(command -v mpicxx) MPIRUN=$(command -v mpirun)
export CC_SERIAL=$(command -v gcc) CXX_SERIAL=$(command -v g++)
export JOBS=$SLURM_CPUS_ON_NODE TEST_RANKS=4 TEST_THREADS=2
export GROMACS_CMAKE_ARGS="-DGMX_SIMD=AVX2_256"   # the SIMD level of the nodes you will run on

cd "$(dirname "$0")/../../.."                     # repository root
scripts/install/install_indus_stack.sh --no-gpu --root "$STACK_ROOT" --print-config
scripts/install/install_indus_stack.sh --no-gpu --root "$STACK_ROOT"
