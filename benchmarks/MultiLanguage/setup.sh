#!/bin/bash
set -euo pipefail

# An empty PYTHON makes PyCall use its Conda environment, where SciPyDiffEq can install SciPy.
if [[ -n "${BENCHMARK_ENV_FILE:-}" ]]; then
    echo 'export PYTHON=""' >> "${BENCHMARK_ENV_FILE}"
fi
