#!/bin/bash
set -euo pipefail

TARGET="${1}"

# Determine the benchmark directory
if [[ -d "${TARGET}" ]]; then
    BENCH_DIR="${TARGET}"
else
    BENCH_DIR="$(dirname "${TARGET}")"
fi

# Run per-benchmark setup script if one exists
# Setup scripts can install system dependencies, configure registries, etc.
BENCHMARK_ENV_FILE="$(mktemp)"
export BENCHMARK_ENV_FILE
if [[ -x "${BENCH_DIR}/setup.sh" ]]; then
    echo "--- Running setup script for ${BENCH_DIR}"
    "${BENCH_DIR}/setup.sh"
    # Source any environment variables exported by the setup script
    if [[ -s "${BENCHMARK_ENV_FILE}" ]]; then
        source "${BENCHMARK_ENV_FILE}"
    fi
fi
rm -f "${BENCHMARK_ENV_FILE}"

# Instantiate, to install the overall project dependencies, and build() for conda
echo "--- Instantiate"
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.build()'

# Run benchmark
echo "+++ Run benchmark for ${TARGET}"
julia --threads=auto --project=. benchmark.jl "${TARGET}"
