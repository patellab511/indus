# Installing INDUS with PLUMED and GROMACS

`install_indus_stack.sh` builds a complete, tested stack from pristine sources into one
directory: PLUMED with INDUS patched in, GROMACS patched with that PLUMED, and the standalone
INDUS driver, all from the same checkout of this repository. It then runs the regression
tests against what it built and writes a manifest.

```
scripts/install/install_indus_stack.sh --print-config    # what would be built, nothing runs
scripts/install/install_indus_stack.sh                   # build everything
source ~/programs/indus-stack/gromacs-2024.3_plumed-2.9.4/env.sh
```

Nothing outside the stack directory is touched; several stacks can coexist.

## What you need

| Requirement | Notes |
| --- | --- |
| C and C++ compiler with C++17 (GROMACS 2024) | gcc 9 or newer, clang 10 or newer |
| CMake 3.18 or newer | GROMACS 2024 needs 3.18; the INDUS driver needs 3.15 |
| make, tar, and curl or wget | |
| MPI (optional, default on) | OpenMPI or MPICH; the wrappers `mpicc`/`mpicxx` and the matching `mpirun` |
| CUDA toolkit (optional) | 11.0 or newer for GROMACS 2024; used when `nvcc` is found or `USE_GPU=yes` |
| Internet access during the `fetch` phase | or pre-downloaded tarballs placed in `<stack>/src` |
| Disk | about 6 GB per stack including sources and build trees |

## Worked examples

**Workstation with a GPU and a system OpenMPI**

```
MPICC=/usr/local/bin/mpicc MPICXX=/usr/local/bin/mpicxx \
CUDA_HOME=/usr/local/cuda-11.7 GMX_CUDA_TARGET_SM=86 \
    scripts/install/install_indus_stack.sh --root $HOME/programs/indus-stack
```

Give the MPI wrappers explicitly whenever more than one MPI is installed (Anaconda ships
MPICH, which may shadow a system OpenMPI on the PATH). `GMX_CUDA_TARGET_SM` is the compute
capability of your GPU (86 for an RTX 3080); leaving it empty builds GROMACS's default set,
which takes longer.

**Laptop or container, CPU only, no MPI**

```
scripts/install/install_indus_stack.sh --no-mpi --no-gpu --jobs 4
```

GROMACS is built with thread-MPI, PLUMED and the driver without MPI. The MPI regression tests
are skipped; everything else runs. Note that a PLUMED without MPI only works with a single
GROMACS rank: run `gmx mdrun -ntmpi 1 -ntomp <threads> -plumed ...` and use OpenMP threads
for parallelism (with INDUS this is the fastest layout anyway, see the README). For
multi-rank GROMACS, build with MPI.

**Cluster with environment modules**

```
module load gcc/12 openmpi/4.1 cuda/12 cmake
MPICC=mpicc MPICXX=mpicxx CUDA_HOME=$CUDA_HOME \
GROMACS_CMAKE_ARGS="-DGMX_SIMD=AVX2_256 -DGMX_FFT_LIBRARY=fftw3 -DFFTWF_LIBRARY=$FFTW_ROOT/lib/libfftw3f.so -DFFTWF_INCLUDE_DIR=$FFTW_ROOT/include" \
TEST_RANKS=4 TEST_THREADS=2 \
    scripts/install/install_indus_stack.sh --root $HOME/software/indus-stack --name gpu-node
```

Build on a node of the same architecture you will run on, or set `GMX_SIMD` for the target.
If the site forbids `mpirun` on login nodes, run with `--from fetch --only build-indus` style
phases on the login node and the `test` phase inside a job. Load the same modules before
`source <stack>/env.sh`.

**Other PLUMED or GROMACS versions**

```
PLUMED_VERSION=2.9.1 GROMACS_VERSION=2023.5 scripts/install/install_indus_stack.sh --name legacy
```

The PLUMED patch name defaults to `gromacs-<GROMACS_VERSION>`; the script checks that the
chosen PLUMED ships it and lists the available names if not. Checksums are known for the
default versions only; for others, pass `PLUMED_SHA256` and `GROMACS_SHA256` from the
projects' download pages, or leave them unset to skip verification (the script says so).

## Settings

All settings are environment variables with defaults; `--print-config` shows the resolved
values. The most useful ones:

| Setting | Default | Meaning |
| --- | --- | --- |
| `PLUMED_VERSION`, `GROMACS_VERSION` | 2.9.4, 2024.3 | release to download |
| `PLUMED_URL`, `GROMACS_URL` | official URLs | alternative download locations |
| `PLUMED_SHA256`, `GROMACS_SHA256` | known for the defaults | tarball checksums |
| `PLUMED_PATCH_ENGINE` | `gromacs-<version>` | name from `plumed patch -l` |
| `STACK_ROOT`, `STACK_NAME` | `~/programs/indus-stack`, `gromacs-<v>_plumed-<v>` | where the stack goes |
| `USE_MPI`, `USE_GPU`, `USE_OPENMP` | yes, auto, yes | `--no-mpi` / `--no-gpu` set the first two |
| `MPICC`, `MPICXX`, `MPIRUN` | `mpicc`, `mpicxx`, next to `MPICC` | the MPI to build and test with |
| `CC_SERIAL`, `CXX_SERIAL` | `gcc`, `g++` | compilers for GROMACS's CMake and non-MPI builds |
| `CUDA_HOME`, `GMX_CUDA_TARGET_SM` | `/usr/local/cuda`, GROMACS default | CUDA toolkit and target architecture |
| `OPT_FLAGS` | `-O3 -g -fPIC` | PLUMED and driver optimization flags |
| `GROMACS_CMAKE_ARGS`, `PLUMED_CONFIGURE_ARGS` | empty | appended verbatim to CMake / configure |
| `JOBS` | all cores | parallel build jobs |
| `TEST_RANKS`, `TEST_THREADS` | 4, 2 | layout for the parallel tests |

## Phases

`fetch`, `patch-plumed`, `build-plumed`, `patch-gromacs`, `build-gromacs`, `build-indus`,
`test`, `env`. Each writes `<stack>/logs/<phase>.log` and a done-marker in `<stack>/stamps/`.
A finished phase is skipped on rerun; `--from <phase>` resumes, `--only <phase>` reruns one.
After changing INDUS sources: `--only patch-plumed`, `--only build-plumed`, `--only test`
(GROMACS loads the PLUMED kernel at runtime and needs no rebuild).

## What the test phase checks

1. The standalone driver's `ctest` suite: every probe volume against stored references,
   serially and (with MPI) under `TEST_RANKS` x `TEST_THREADS`.
2. INDUS inside PLUMED via `plumed driver`, one test per probe volume, against references
   that were themselves validated against the standalone ones.
3. A 500-step MD of a small water box under the built GROMACS with INDUS restraining Ntilde:
   the standalone driver recomputes Ntilde from every frame and must match PLUMED's output.

Results go to `<stack>/test-results.txt` and the manifest.

## Layout of a finished stack

```
<stack>/
  env.sh          source this: PATH, LD_LIBRARY_PATH, PLUMED_KERNEL, the matching mpirun
  manifest.txt    versions, URLs, flags, INDUS commit, test results
  config.txt      the resolved settings this stack was built with
  plumed/ gromacs/ indus/      installed components
  src/ build/ logs/ stamps/    sources, build trees, one log per phase, done-markers
```

## Troubleshooting

- **`checksum mismatch`**: the download is corrupt or the mirror differs. Delete the file in
  `<stack>/src` to re-download, or set the `*_SHA256` variable if you know the new value.
- **`no patch named gromacs-X`**: that PLUMED release does not support that GROMACS. Pick a
  pair from the list the script prints (PLUMED 2.9.4: GROMACS 2020.7, 2021.7, 2022.5, 2023.5,
  2024.3).
- **MPI tests fail with four identical single-rank outputs**: `mpirun` does not match the MPI
  the binaries were built with. Set `MPIRUN` (and `MPICC`/`MPICXX`) to one installation.
- **`No CMAKE_CUDA_COMPILER could be found`**: set `CUDA_HOME` to the toolkit directory that
  contains `bin/nvcc`, or build with `--no-gpu`.
- **PLUMED aborts at startup with "component ... has not been added to the manual"**: an
  INDUS source tree older than this repository's `plumed_patch`; rerun `patch-plumed` from
  the current checkout.
- **A phase failed**: read `<stack>/logs/<phase>.log`, fix the cause, rerun with
  `--from <phase>`.
