#!/bin/bash
# install_indus_stack.sh - build INDUS + PLUMED + GROMACS from pristine sources
#
# What it builds, in order (each step is a "phase" that can be rerun on its own):
#   fetch          download PLUMED and GROMACS tarballs into $STACK/src
#   patch-plumed   copy this repo's INDUS sources into PLUMED (plumed_patch/patch_plumed.sh)
#   build-plumed   configure + make + install PLUMED            -> $STACK/plumed
#   patch-gromacs  'plumed patch' the GROMACS sources (runtime mode)
#   build-gromacs  cmake + make + install GROMACS (gmx_mpi)     -> $STACK/gromacs
#   build-indus    cmake + make the standalone INDUS driver     -> $STACK/indus
#   test           INDUS ctest, PLUMED-driver tests, short MD with INDUS active
#   env            write $STACK/env.sh and $STACK/manifest.txt
#
# Usage:
#   scripts/install/install_indus_stack.sh [options]
#     --root DIR       parent directory for stacks        (default: $HOME/programs/indus-stack)
#     --name NAME      stack name, a subdirectory of root (default: gromacs-<ver>_plumed-<ver>)
#     --jobs N         parallel build jobs                (default: all cores)
#     --no-gpu         build GROMACS without CUDA
#     --no-mpi         build everything without MPI (GROMACS then uses thread-MPI)
#     --from PHASE     start at PHASE and run the rest
#     --only PHASE     run just PHASE
#     --list           list phases and exit
#     --print-config   show the resolved settings and exit (nothing is built)
#     --config FILE    read settings (NAME=value lines, shell syntax) from FILE first
#   Every setting in the CONFIG block below can also be overridden with an
#   environment variable of the same name, e.g.  GROMACS_VERSION=2023.5 PLUMED_VERSION=2.9.1 ...
#   Precedence: command-line flags > environment > --config file > defaults.
#   An annotated example is scripts/install/stack.conf.example.
#   Worked examples (GPU workstation, CPU laptop, cluster with modules): scripts/install/README.md
#
# Layout of a finished stack:
#   $STACK/{src,build,logs,stamps}   sources, build trees, one log per phase, done-markers
#   $STACK/{plumed,gromacs,indus}    installed components
#   $STACK/env.sh                    source this to use the stack
#   $STACK/manifest.txt              versions, flags, INDUS commit, test results
set -euo pipefail

# --config FILE is applied before the defaults below, so the file can set any CONFIG variable.
# Variables already set in the environment keep their value (environment beats file).
for (( i=1; i<=$#; i++ )); do
	if [[ "${!i}" == --config ]]; then
		j=$((i+1)); config_file="${!j}"
		[[ -f "$config_file" ]] || { echo "ERROR: config file not found: $config_file" >&2; exit 1; }
		set -a; eval "$(sed -E 's/^([A-Za-z_][A-Za-z0-9_]*)=/\1=${\1:-}; [[ -n "${\1}" ]] || \1=/' "$config_file")"; set +a
	fi
done

############################################################
### CONFIG  (override any of these from the environment) ###
############################################################

PLUMED_VERSION="${PLUMED_VERSION:-2.9.4}"
GROMACS_VERSION="${GROMACS_VERSION:-2024.3}"
# Name of the PLUMED patch for this GROMACS; 'plumed patch -l' lists them
PLUMED_PATCH_ENGINE="${PLUMED_PATCH_ENGINE:-gromacs-${GROMACS_VERSION}}"

PLUMED_URL="${PLUMED_URL:-https://github.com/plumed/plumed2/releases/download/v${PLUMED_VERSION}/plumed-${PLUMED_VERSION}.tgz}"
GROMACS_URL="${GROMACS_URL:-https://ftp.gromacs.org/gromacs/gromacs-${GROMACS_VERSION}.tar.gz}"
# SHA-256 of the tarballs. Known for the default versions (computed from the official
# downloads, 2026-10-05); set them yourself for other versions, or leave empty to skip the check.
PLUMED_SHA256="${PLUMED_SHA256:-$( [[ $PLUMED_VERSION == 2.9.4 ]] && echo 032c99bda66b20f0710e38872b7f3fa069ad647b5f73e83a26228f60052bb54d || true )}"
GROMACS_SHA256="${GROMACS_SHA256:-$( [[ $GROMACS_VERSION == 2024.3 ]] && echo bbda056ee59390be7d58d84c13a9ec0d4e3635617adf2eb747034922cba1f029 || true )}"

STACK_ROOT="${STACK_ROOT:-$HOME/programs/indus-stack}"
STACK_NAME="${STACK_NAME:-gromacs-${GROMACS_VERSION}_plumed-${PLUMED_VERSION}}"
JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"

USE_MPI="${USE_MPI:-yes}"            # yes|no
USE_GPU="${USE_GPU:-auto}"           # yes|no|auto (auto = yes if nvcc is found)
USE_OPENMP="${USE_OPENMP:-yes}"      # yes|no

CC_SERIAL="${CC_SERIAL:-gcc}"        # compilers used when USE_MPI=no
CXX_SERIAL="${CXX_SERIAL:-g++}"
MPICC="${MPICC:-mpicc}"              # MPI wrappers used when USE_MPI=yes
MPICXX="${MPICXX:-mpicxx}"
CUDA_HOME="${CUDA_HOME:-/usr/local/cuda}"
GMX_CUDA_TARGET_SM="${GMX_CUDA_TARGET_SM:-}"   # e.g. "86" for an RTX 3080; empty = GROMACS default set

OPT_FLAGS="${OPT_FLAGS:--O3 -g -fPIC}"          # shared by PLUMED and the INDUS driver

# FFTW for GROMACS. Default: GROMACS builds its own from this tarball, which the fetch phase
# downloads (compute nodes are often offline). Set USE_OWN_FFTW=no to use a site FFTW instead
# and give its location through GROMACS_CMAKE_ARGS (-DGMX_FFT_LIBRARY=fftw3 ...).
USE_OWN_FFTW="${USE_OWN_FFTW:-yes}"
FFTW_VERSION="${FFTW_VERSION:-3.3.8}"            # the version GROMACS 2024 expects
FFTW_URL="${FFTW_URL:-http://www.fftw.org/fftw-${FFTW_VERSION}.tar.gz}"
FFTW_MD5="${FFTW_MD5:-$( [[ $FFTW_VERSION == 3.3.8 ]] && echo 8aac833c943d8e90d51b697b27d4384d || true )}"
FFTW_SHA256="${FFTW_SHA256:-$( [[ $FFTW_VERSION == 3.3.8 ]] && echo 6113262f6e92c5bd474f2875fa1b01054c4ad5040f6b0da7c03c98821d9ae303 || true )}"

# Extra, site-specific arguments appended verbatim (e.g. a cluster's FFTW or SIMD settings)
GROMACS_CMAKE_ARGS="${GROMACS_CMAKE_ARGS:-}"        # e.g. "-DGMX_SIMD=AVX2_256 -DGMX_FFT_LIBRARY=fftw3"
PLUMED_CONFIGURE_ARGS="${PLUMED_CONFIGURE_ARGS:-}"  # e.g. "--enable-modules=all"

# Test phase: MPI ranks x OpenMP threads for the parallel tests
TEST_RANKS="${TEST_RANKS:-4}"
TEST_THREADS="${TEST_THREADS:-2}"

############################################################
### Derived paths and helpers                            ###
############################################################

INDUS_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PHASES=(fetch patch-plumed build-plumed patch-gromacs build-gromacs build-indus test env)

from_phase=""; only_phase=""; print_config=0
while [[ $# -gt 0 ]]; do
	case "$1" in
		--root)   STACK_ROOT="$2"; shift 2 ;;
		--name)   STACK_NAME="$2"; shift 2 ;;
		--jobs)   JOBS="$2"; shift 2 ;;
		--no-gpu) USE_GPU=no; shift ;;
		--no-mpi) USE_MPI=no; shift ;;
		--from)   from_phase="$2"; shift 2 ;;
		--only)   only_phase="$2"; shift 2 ;;
		--list)   printf '%s\n' "${PHASES[@]}"; exit 0 ;;
		--print-config) print_config=1; shift ;;
		--config) shift 2 ;;   # already applied above
		-h|--help) sed -n '2,36p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
		*) echo "unknown option: $1" >&2; exit 1 ;;
	esac
done

STACK="$STACK_ROOT/$STACK_NAME"
SRC="$STACK/src";   BUILD="$STACK/build";  LOGS="$STACK/logs";  STAMPS="$STACK/stamps"
PLUMED_PREFIX="$STACK/plumed"; GROMACS_PREFIX="$STACK/gromacs"; INDUS_PREFIX="$STACK/indus"
PLUMED_SRC="$SRC/plumed-$PLUMED_VERSION"
GROMACS_SRC="$SRC/gromacs-$GROMACS_VERSION"

if [[ "$USE_GPU" == auto ]]; then
	if [[ -x "$CUDA_HOME/bin/nvcc" ]]; then USE_GPU=yes; else USE_GPU=no; fi
fi
if [[ "$USE_MPI" == yes ]]; then CC="$MPICC"; CXX="$MPICXX"; else CC="$CC_SERIAL"; CXX="$CXX_SERIAL"; fi
# mpirun from the same MPI installation as the wrappers (overridable)
MPIRUN="${MPIRUN:-$(dirname "$(command -v "$MPICC" || echo /usr/bin/mpicc)")/mpirun}"

say()  { printf '\n==> %s\n' "$*"; }
die()  { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "required tool not found: $1"; }

# The GROMACS executable of this stack (gmx_mpi for MPI builds, gmx otherwise)
gmx_bin() {
	if [[ -x "$GROMACS_PREFIX/bin/gmx_mpi" ]]; then echo "$GROMACS_PREFIX/bin/gmx_mpi"; else echo "$GROMACS_PREFIX/bin/gmx"; fi
}

# Environment needed to run the freshly built PLUMED / GROMACS
stack_env() {
	export PATH="$GROMACS_PREFIX/bin:$PLUMED_PREFIX/bin:$INDUS_PREFIX/bin:$PATH"
	# the stack's own mpirun must win over any other MPI on the PATH (tests call bare 'mpirun')
	[[ "$USE_MPI" == yes ]] && export PATH="$(dirname "$MPIRUN"):$PATH"
	export LD_LIBRARY_PATH="$PLUMED_PREFIX/lib:$GROMACS_PREFIX/lib:${LD_LIBRARY_PATH:-}"
	export PLUMED_KERNEL="$PLUMED_PREFIX/lib/libplumedKernel.so"
}

show_config() {
	cat <<EOF
stack:            $STACK
PLUMED:           $PLUMED_VERSION  $PLUMED_URL
GROMACS:          $GROMACS_VERSION  $GROMACS_URL  (patch engine: $PLUMED_PATCH_ENGINE)
checksums:        plumed=${PLUMED_SHA256:-unchecked}  gromacs=${GROMACS_SHA256:-unchecked}
FFTW:             own build = $USE_OWN_FFTW  (fftw $FFTW_VERSION, $FFTW_URL, sha256=${FFTW_SHA256:-unchecked})
MPI:              $USE_MPI  (CC=$CC CXX=$CXX mpirun=$MPIRUN)
GPU:              $USE_GPU  (CUDA_HOME=$CUDA_HOME, target SM: ${GMX_CUDA_TARGET_SM:-GROMACS default})
OpenMP:           $USE_OPENMP
flags:            OPT_FLAGS='$OPT_FLAGS'  GROMACS_CMAKE_ARGS='$GROMACS_CMAKE_ARGS'  PLUMED_CONFIGURE_ARGS='$PLUMED_CONFIGURE_ARGS'
jobs:             $JOBS      tests: $TEST_RANKS ranks x $TEST_THREADS threads
INDUS source:     $INDUS_ROOT
EOF
}

# Fail early, with the setting to change, rather than deep inside a build
check_prerequisites() {
	need tar; need make; need cmake
	command -v curl >/dev/null 2>&1 || command -v wget >/dev/null 2>&1 || die "need curl or wget to download sources"
	command -v "$CC"  >/dev/null 2>&1 || die "C compiler '$CC' not found (set MPICC, or CC_SERIAL with --no-mpi)"
	command -v "$CXX" >/dev/null 2>&1 || die "C++ compiler '$CXX' not found (set MPICXX, or CXX_SERIAL with --no-mpi)"
	if [[ "$USE_MPI" == yes ]]; then
		[[ -x "$MPIRUN" ]] || die "mpirun not found at '$MPIRUN' (set MPIRUN, or use --no-mpi)"
	fi
	if [[ "$USE_GPU" == yes ]]; then
		[[ -x "$CUDA_HOME/bin/nvcc" ]] || die "nvcc not found at '$CUDA_HOME/bin/nvcc' (set CUDA_HOME, or use --no-gpu)"
	fi
	local cmake_ver; cmake_ver=$(cmake --version | awk 'NR==1{print $3}')
	echo "    tools: cmake $cmake_ver, $("$CXX" --version 2>/dev/null | head -1)"
}

sha256_of() {
	if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{print $1}'
	else shasum -a 256 "$1" | awk '{print $1}'; fi
}

fetch_file() {   # fetch_file DEST URL [SHA256]
	local dest="$1" url="$2" sum="${3:-}"
	if [[ ! -f "$dest" ]]; then
		if command -v curl >/dev/null 2>&1; then curl -L --fail -o "$dest" "$url" || true; else wget -O "$dest" "$url" || true; fi
		[[ -s "$dest" ]] || { rm -f "$dest"; die "download failed: $url
    No network here? On clusters the compute nodes often cannot reach the internet: run
    '$0 --only fetch' on a login node first (or copy the tarballs into $SRC), then rerun."; }
	fi
	if [[ -n "$sum" ]]; then
		local got; got=$(sha256_of "$dest")
		[[ "$got" == "$sum" ]] || die "checksum mismatch for $dest: expected $sum, got $got (delete the file to re-download)"
		echo "checksum ok: $(basename "$dest")"
	else
		echo "checksum not verified (no SHA256 given): $(basename "$dest")"
	fi
}

# run_phase NAME: skip if already done, else run phase_NAME with output in logs/NAME.log
run_phase() {
	local name="$1" fn="phase_${1//-/_}"
	if [[ -f "$STAMPS/$name.done" && -z "$only_phase" ]]; then
		say "$name: already done (remove $STAMPS/$name.done to redo)"; return
	fi
	say "$name  (log: $LOGS/$name.log)"
	local t0=$SECONDS
	if "$fn" > "$LOGS/$name.log" 2>&1; then
		touch "$STAMPS/$name.done"
		echo "    done in $(( SECONDS - t0 )) s"
	else
		echo "    FAILED after $(( SECONDS - t0 )) s - see $LOGS/$name.log"; tail -n 30 "$LOGS/$name.log"; exit 1
	fi
}

############################################################
### Phases                                               ###
############################################################

phase_fetch() {
	cd "$SRC"
	fetch_file "plumed-$PLUMED_VERSION.tgz"      "$PLUMED_URL"  "$PLUMED_SHA256"
	fetch_file "gromacs-$GROMACS_VERSION.tar.gz" "$GROMACS_URL" "$GROMACS_SHA256"
	[[ "$USE_OWN_FFTW" == yes ]] && fetch_file "fftw-$FFTW_VERSION.tar.gz" "$FFTW_URL" "$FFTW_SHA256"
	[[ -d "$PLUMED_SRC" ]]  || tar xzf "plumed-$PLUMED_VERSION.tgz"
	[[ -d "$GROMACS_SRC" ]] || tar xzf "gromacs-$GROMACS_VERSION.tar.gz"
	ls -la "$SRC"
}

phase_patch_plumed() {
	# Copies $INDUS_ROOT/src/* into $PLUMED_SRC/src/orderparameters (a PLUMED module)
	"$INDUS_ROOT/plumed_patch/patch_plumed.sh" "$PLUMED_SRC"
	# Record which INDUS source went in; "-dirty" means uncommitted changes were included
	{ git -C "$INDUS_ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown
	  git -C "$INDUS_ROOT" diff --quiet HEAD 2>/dev/null || echo "-dirty"; } | tr -d '\n' > "$STACK/indus_commit.txt"
	echo >> "$STACK/indus_commit.txt"
}

phase_build_plumed() {
	cd "$PLUMED_SRC"
	local flags="$OPT_FLAGS -std=c++11" conf=(--prefix="$PLUMED_PREFIX")
	if [[ "$USE_MPI" == yes ]]; then flags="$flags -DMPI_ENABLED"; conf+=(--enable-mpi); else conf+=(--disable-mpi); fi
	if [[ "$USE_OPENMP" == yes ]]; then flags="$flags -fopenmp"; conf+=(--enable-openmp); else conf+=(--disable-openmp); fi
	# -DMPI_ENABLED / -fopenmp are what the INDUS module keys its MPI and OpenMP code on
	local ldflags=""; [[ "$USE_OPENMP" == yes ]] && ldflags="-fopenmp"
	# shellcheck disable=SC2086  (PLUMED_CONFIGURE_ARGS is meant to word-split)
	CC="$CC" CXX="$CXX" CFLAGS="$flags" CXXFLAGS="$flags" LDFLAGS="$ldflags" LIBS="-lstdc++" \
		./configure "${conf[@]}" $PLUMED_CONFIGURE_ARGS
	make -j "$JOBS"
	make install
	# Warnings from the INDUS module only, for the record
	grep -n 'orderparameters/.*warning' "$LOGS/build-plumed.log" > "$LOGS/indus-module-warnings.log" || true
	stack_env; plumed info --version
}

phase_patch_gromacs() {
	stack_env
	if ! plumed patch -l 2>/dev/null | grep -qw "$PLUMED_PATCH_ENGINE"; then
		echo "PLUMED $PLUMED_VERSION has no patch named '$PLUMED_PATCH_ENGINE'. Available:"
		plumed patch -l 2>/dev/null | grep -E '^\s*[a-z]' | tr -s ' \n' ' '; echo
		die "set GROMACS_VERSION to a supported release, or PLUMED_PATCH_ENGINE to one of the names above"
	fi
	cd "$GROMACS_SRC"
	[[ -f .indus-plumed-patched ]] && plumed patch -r --engine "$PLUMED_PATCH_ENGINE" || true
	# runtime mode: GROMACS loads the kernel named by $PLUMED_KERNEL, so kernels can be swapped later
	plumed patch -p --engine "$PLUMED_PATCH_ENGINE" --runtime
	touch .indus-plumed-patched
}

phase_build_gromacs() {
	local bdir="$BUILD/gromacs-$GROMACS_VERSION"
	mkdir -p "$bdir"; cd "$bdir"
	local args=(-DCMAKE_INSTALL_PREFIX="$GROMACS_PREFIX" -DCMAKE_BUILD_TYPE=Release -DREGRESSIONTEST_DOWNLOAD=OFF)
	if [[ "$USE_OWN_FFTW" == yes ]]; then
		# built from the tarball the fetch phase downloaded, so this works offline
		args+=(-DGMX_BUILD_OWN_FFTW=ON -DGMX_BUILD_OWN_FFTW_URL="$SRC/fftw-$FFTW_VERSION.tar.gz")
		[[ -n "$FFTW_MD5" ]] && args+=(-DGMX_BUILD_OWN_FFTW_MD5="$FFTW_MD5")
	else
		args+=(-DGMX_BUILD_OWN_FFTW=OFF)
	fi
	if [[ "$USE_MPI" == yes ]]; then
		args+=(-DGMX_MPI=ON -DMPI_C_COMPILER="$MPICC" -DMPI_CXX_COMPILER="$MPICXX")
	else
		args+=(-DGMX_MPI=OFF -DGMX_THREAD_MPI=ON)
	fi
	[[ "$USE_OPENMP" == yes ]] && args+=(-DGMX_OPENMP=ON) || args+=(-DGMX_OPENMP=OFF)
	if [[ "$USE_GPU" == yes ]]; then
		args+=(-DGMX_GPU=CUDA -DCMAKE_CUDA_COMPILER="$CUDA_HOME/bin/nvcc" -DCMAKE_CUDA_HOST_COMPILER="$CXX_SERIAL")
		[[ -n "$GMX_CUDA_TARGET_SM" ]] && args+=(-DGMX_CUDA_TARGET_SM="$GMX_CUDA_TARGET_SM")
		export PATH="$CUDA_HOME/bin:$PATH"
	else
		args+=(-DGMX_GPU=OFF)
	fi
	# Plain system compilers here; MPI comes in through the MPI_*_COMPILER hints (nvcc dislikes wrappers)
	# shellcheck disable=SC2086  (GROMACS_CMAKE_ARGS is meant to word-split)
	CC="$CC_SERIAL" CXX="$CXX_SERIAL" cmake "$GROMACS_SRC" "${args[@]}" $GROMACS_CMAKE_ARGS
	make -j "$JOBS"
	make install
	stack_env; "$(gmx_bin)" --version | grep -E 'GROMACS version|MPI library|GPU support|OpenMP'
}

phase_build_indus() {
	local bdir="$BUILD/indus"
	mkdir -p "$bdir" "$INDUS_PREFIX/bin"; cd "$bdir"
	local mpi=OFF omp=OFF
	[[ "$USE_MPI" == yes ]] && mpi=ON; [[ "$USE_OPENMP" == yes ]] && omp=ON
	CC="$CC" CXX="$CXX" CXXFLAGS="$OPT_FLAGS" cmake "$INDUS_ROOT" -DCMAKE_BUILD_TYPE=Release \
		-DMPI_ENABLED=$mpi -DOPENMP_ENABLED=$omp -DGPTL_ENABLED=OFF
	make -j "$JOBS"
	cp bin/indus "$INDUS_PREFIX/bin/indus"
}

phase_test() {
	stack_env
	local results="$STACK/test-results.txt"; : > "$results"
	# 1) standalone driver regression tests (test/CMakeLists.txt)
	( cd "$BUILD/indus" && ctest --output-on-failure ) && echo "indus ctest: PASS" >> "$results" || echo "indus ctest: FAIL" >> "$results"
	# 2) INDUS inside PLUMED, via 'plumed driver' on stored trajectories
	( "$INDUS_ROOT/plumed_patch/test/run_tests.sh" plumed | tee "$LOGS/test-plumed-driver.log" | grep -q FAILED ) \
		&& echo "plumed driver tests: FAIL" >> "$results" || echo "plumed driver tests: PASS" >> "$results"
	# 3) a short MD run with INDUS biasing a sphere in a small water box
	if md_smoke_test; then echo "gromacs+plumed MD smoke: PASS" >> "$results"; else echo "gromacs+plumed MD smoke: FAIL" >> "$results"; fi
	cat "$results"
	! grep -q FAIL "$results"
}

# Build a 3 nm SPC/E water box with the fresh GROMACS, run 500 steps with an INDUS
# RESTRAINT on Ntilde in a sphere, and check that PLUMED's step-0 Ntilde equals the
# standalone driver's value for the same configuration.
md_smoke_test() {
	local d="$BUILD/md-smoke" gmx; gmx="$(gmx_bin)"
	rm -rf "$d"; mkdir -p "$d"; cd "$d"
	"$gmx" solvate -cs spc216.gro -box 3 3 3 -o water.gro > solvate.log 2>&1 || { cat solvate.log; return 1; }
	local nsol; nsol=$(grep -c OW water.gro)
	printf '#include "oplsaa.ff/forcefield.itp"\n#include "oplsaa.ff/spce.itp"\n[ system ]\nwater\n[ molecules ]\nSOL %d\n' "$nsol" > topol.top
	# short minimization first: a freshly solvated box has close contacts that break SETTLE
	printf 'integrator=steep\nnsteps=200\nemtol=1000\nnstlist=10\nrcoulomb=1.0\nrvdw=1.0\ncoulombtype=PME\n' > em.mdp
	"$gmx" grompp -f em.mdp -c water.gro -p topol.top -o em.tpr -maxwarn 2 > grompp_em.log 2>&1 || { cat grompp_em.log; return 1; }
	local ntmpi=(); [[ "$USE_MPI" == yes ]] || ntmpi=(-ntmpi 1)   # thread-MPI GROMACS wants -ntmpi with -ntomp
	"$gmx" mdrun -deffnm em "${ntmpi[@]}" -ntomp "$TEST_THREADS" > mdrun_em.log 2>&1 || { tail -20 mdrun_em.log; return 1; }
	# every step goes to md.xtc at high precision so the standalone driver can recompute Ntilde per frame
	printf 'integrator=md\nnsteps=500\ndt=0.002\nnstlist=10\nrcoulomb=1.0\nrvdw=1.0\ncoulombtype=PME\ntcoupl=v-rescale\ntc-grps=System\ntau_t=0.5\nref_t=300\nconstraints=h-bonds\ngen_vel=yes\ngen_temp=300\nnstcalcenergy=1\nnstxout-compressed=1\ncompressed-x-precision=100000\n' > md.mdp
	printf 'Target = [ atom_index 1-%d:3 ]\nProbeVolume = {\n  type = sphere\n  center = [ 1.5 1.5 1.5 ]\n  r_max = 0.6\n  sigma = 0.01\n  alpha_c = 0.02\n}\nBias = {\n  order_parameter = ntilde\n  x_star = 0.0\n  kappa = 0.0\n}\n' "$((nsol*3))" > indus.input
	printf 'indus: INDUS INPUTFILE=indus.input\nr: RESTRAINT ARG=indus.ntilde AT=20.0 KAPPA=0.5\nPRINT ARG=indus.n,indus.ntilde,r.bias STRIDE=1 FILE=plumed.out\n' > plumed.dat
	"$gmx" grompp -f md.mdp -c em.gro -p topol.top -o md.tpr -maxwarn 2 > grompp.log 2>&1 || { cat grompp.log; return 1; }
	if [[ "$USE_MPI" == yes ]]; then
		"$MPIRUN" -np "$TEST_RANKS" "$gmx" mdrun -deffnm md -plumed plumed.dat -ntomp "$TEST_THREADS" > mdrun.log 2>&1 || { tail -30 mdrun.log; return 1; }
	else
		# A PLUMED built without MPI cannot take the communicator GROMACS passes when it runs
		# more than one thread-MPI rank, so a --no-mpi stack runs GROMACS with one rank
		"$gmx" mdrun -deffnm md -plumed plumed.dat -ntmpi 1 -ntomp $((TEST_RANKS * TEST_THREADS)) > mdrun.log 2>&1 || { tail -30 mdrun.log; return 1; }
	fi
	# Standalone driver on the MD trajectory: Ntilde must match what PLUMED printed at every step
	# (xtc positions carry 1e-5 nm rounding, so allow 1e-2 on Ntilde)
	printf 'GroFile = em.gro\nXtcFile = md.xtc\n' > indus_standalone.input; cat indus.input >> indus_standalone.input
	indus indus_standalone.input > indus_standalone.log 2>&1 || { tail -20 indus_standalone.log; return 1; }
	paste <(awk '!/^#/ {print $3}' plumed.out) <(awk '!/^#/ {print $3}' time_samples_indus.out) \
		| awk 'NF==2 { n++; d = $1 - $2; if (d < 0) d = -d; if (d > max) max = d }
		       END { printf "Ntilde PLUMED vs standalone over %d frames: max |diff| = %.2e\n", n, max; exit !(n > 100 && max < 1e-2) }'
}

phase_env() {
	cat > "$STACK/env.sh" <<EOF
# INDUS stack: GROMACS $GROMACS_VERSION + PLUMED $PLUMED_VERSION + INDUS ($(cat "$STACK/indus_commit.txt"))
export INDUS_STACK="$STACK"
export PATH="$GROMACS_PREFIX/bin:$PLUMED_PREFIX/bin:$INDUS_PREFIX/bin:\$PATH"
export LD_LIBRARY_PATH="$PLUMED_PREFIX/lib:$GROMACS_PREFIX/lib:\${LD_LIBRARY_PATH:-}"
export PLUMED_KERNEL="$PLUMED_PREFIX/lib/libplumedKernel.so"
$([[ "$USE_MPI" == yes ]] && echo "export PATH=\"$(dirname "$MPIRUN"):\$PATH\"   # the MPI this stack was built with")
EOF
	{
		echo "built on $(date) by $(whoami)@$(hostname)"
		echo "INDUS commit:    $(cat "$STACK/indus_commit.txt")"
		echo "PLUMED:          $PLUMED_VERSION  ($PLUMED_URL)"
		echo "GROMACS:         $GROMACS_VERSION  ($GROMACS_URL), patch engine $PLUMED_PATCH_ENGINE, runtime mode"
		echo "MPI=$USE_MPI  GPU=$USE_GPU  OpenMP=$USE_OPENMP  CC=$CC  CXX=$CXX  OPT_FLAGS='$OPT_FLAGS'"
		echo "tests:"; sed 's/^/  /' "$STACK/test-results.txt" 2>/dev/null || echo "  (test phase not run)"
	} > "$STACK/manifest.txt"
	cat "$STACK/manifest.txt"
	echo; echo "To use this stack:  source $STACK/env.sh"
}

############################################################
### Main                                                 ###
############################################################

if [[ $print_config -eq 1 ]]; then show_config; exit 0; fi
say "stack: $STACK"
echo "    PLUMED $PLUMED_VERSION, GROMACS $GROMACS_VERSION (patch $PLUMED_PATCH_ENGINE), INDUS from $INDUS_ROOT"
echo "    MPI=$USE_MPI GPU=$USE_GPU OpenMP=$USE_OPENMP  CC=$CC CXX=$CXX  jobs=$JOBS  (--print-config shows everything)"
check_prerequisites
mkdir -p "$SRC" "$BUILD" "$LOGS" "$STAMPS"
show_config > "$STACK/config.txt"

started=0
for p in "${PHASES[@]}"; do
	if [[ -n "$only_phase" ]]; then [[ "$p" == "$only_phase" ]] && run_phase "$p"; continue; fi
	[[ -n "$from_phase" && "$p" == "$from_phase" ]] && started=1
	[[ -n "$from_phase" && $started -eq 0 ]] && continue
	run_phase "$p"
done
