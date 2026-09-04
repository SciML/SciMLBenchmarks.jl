#!/bin/bash
# Per-benchmark setup hook. Run by .github/scripts/build_benchmark.sh BEFORE the
# project is instantiated. Its job: make MPI.jl, PETSc_jll, and (later) HYPRE_jll
# agree on a single MPI backend so multi-rank runs don't crash on ABI mismatch.
#
# Any variable written to $BENCHMARK_ENV_FILE is sourced back into the build env.
set -euo pipefail

BENCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "--- LinearSolveDistributed setup: pinning MPI backend"

# Load the project containing the committed MPICH preference before downstream setup.
# PETSc_jll's default artifact is built against MPICH, so this keeps the stack on one ABI.
julia --project="${BENCH_DIR}" -e '
    using MPIPreferences
'

# Expose mpiexec to any downstream step that wants it (the .jmd itself gets it
# from MPI.jl via `using MPI; mpiexec()`, so this is belt-and-suspenders).
# BENCHMARK_ENV_FILE is set by build_benchmark.sh in CI; guard so the script is
# also runnable by hand under `set -u`.
if [[ -n "${BENCHMARK_ENV_FILE:-}" ]]; then
    echo "JULIA_MPI_BINARY=MPICH_jll" >> "$BENCHMARK_ENV_FILE"
fi

echo "--- MPI backend pinned"
