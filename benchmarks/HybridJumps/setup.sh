#!/bin/bash
set -euo pipefail

BENCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

JULIA_PKG_PRECOMPILE_AUTO=0 julia --startup-file=no --project="${BENCH_DIR}" -e '
    using Pkg
    Pkg.instantiate()
    # Rebuild Conda before loading it: its module body errors when the depot
    # conda env dir is missing (persistent CI depots can lose it).
    Pkg.build("Conda")
    using Conda
    Conda.add(["python=3.8", "conda<24", "numpy"], Conda.ROOTENV)
    Conda.pip_interop(true, Conda.ROOTENV)
    Conda.pip("install", "tick==0.7.0.1", Conda.ROOTENV)
    ENV["PYTHON"] = ""
    Pkg.build("PyCall")
'

if [[ -n "${BENCHMARK_ENV_FILE:-}" ]]; then
    echo 'export PYTHON=""' >> "${BENCHMARK_ENV_FILE}"
fi
