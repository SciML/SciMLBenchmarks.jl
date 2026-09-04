#!/bin/bash
set -euo pipefail

# PyCall and RCall use Conda; R also needs Conda's newer libstdc++ ahead of Julia's.
if [[ -n "${BENCHMARK_ENV_FILE:-}" ]]; then
    echo 'export PYTHON=""' >> "${BENCHMARK_ENV_FILE}"
    echo 'export R_HOME="*"' >> "${BENCHMARK_ENV_FILE}"
    echo 'export CONDA_JL_HOME="${CONDA_JL_HOME:-${HOME}/.julia/conda/3/x86_64}"' >> "${BENCHMARK_ENV_FILE}"
    echo 'export LD_LIBRARY_PATH="${CONDA_JL_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"' >> "${BENCHMARK_ENV_FILE}"
fi
