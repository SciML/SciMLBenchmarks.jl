#!/bin/bash
set -euo pipefail

BENCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHON="" R_HOME="*"
export CONDA_JL_HOME="${CONDA_JL_HOME:-${HOME}/.julia/conda/SciMLBenchmarks/MultiLanguage}"
export JULIA_PKG_PRECOMPILE_AUTO=0
# RCall needs Conda's C++ runtime before locating R.
R_HOME="_" julia --project="${BENCH_DIR}" -e '
        using Pkg
        Pkg.instantiate()
        Pkg.build("Conda")
        using Conda
        Conda.add(["scipy", "r-desolve"]; channel = "conda-forge")
    '

# PyCall and RCall use Conda; R also needs Conda's newer libstdc++ ahead of Julia's.
LD_LIBRARY_PATH="${CONDA_JL_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}" \
    julia --project="${BENCH_DIR}" -e 'using Pkg; Pkg.build(["PyCall", "RCall"])'

if [[ -n "${BENCHMARK_ENV_FILE:-}" ]]; then
    echo 'export PYTHON=""' >> "${BENCHMARK_ENV_FILE}"
    echo 'export R_HOME="*"' >> "${BENCHMARK_ENV_FILE}"
    echo 'export CONDA_JL_HOME="${CONDA_JL_HOME:-${HOME}/.julia/conda/SciMLBenchmarks/MultiLanguage}"' >> "${BENCHMARK_ENV_FILE}"
    echo 'export LD_LIBRARY_PATH="${CONDA_JL_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"' >> "${BENCHMARK_ENV_FILE}"
fi
