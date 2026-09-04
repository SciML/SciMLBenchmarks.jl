#!/bin/bash
set -euo pipefail

# These sentinel values make PyCall and RCall use their managed Conda environments.
if [[ -n "${BENCHMARK_ENV_FILE:-}" ]]; then
    echo 'export PYTHON=""' >> "${BENCHMARK_ENV_FILE}"
    echo 'export R_HOME="*"' >> "${BENCHMARK_ENV_FILE}"
fi
