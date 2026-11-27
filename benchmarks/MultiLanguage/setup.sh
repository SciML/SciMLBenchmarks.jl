#!/bin/bash
set -euo pipefail

bench_dir="$(dirname "${BASH_SOURCE[0]}")"
export CONDA_JL_HOME="${CONDA_JL_HOME:-${HOME}/.julia/conda/SciMLBenchmarks/MultiLanguage}"
export PYTHON=""
export JULIA_PKG_PRECOMPILE_AUTO=0

R_HOME="_" julia --project="${bench_dir}" -e '
    using Pkg
    Pkg.instantiate()
    Pkg.build("Conda")
    using Conda
    Conda.add(["scipy", "r-desolve"]; channel = "conda-forge", satisfied_skip_solve = true)
'

export R_HOME="${CONDA_JL_HOME}/lib/R"
# R needs Conda's C++ runtime. Julia's Downloads needs Julia's own libcurl.
if [[ -f "${CONDA_JL_HOME}/lib/libstdc++.so.6" ]]; then
    export LD_PRELOAD="${CONDA_JL_HOME}/lib/libstdc++.so.6${LD_PRELOAD:+:${LD_PRELOAD}}"
fi
julia --project="${bench_dir}" -e 'using Pkg; Pkg.build(["PyCall", "RCall"])'

if [[ -n "${BENCHMARK_ENV_FILE:-}" ]]; then
    printf 'export CONDA_JL_HOME=%q\nexport PYTHON=%q\nexport R_HOME=%q\n' \
        "${CONDA_JL_HOME}" "${PYTHON}" "${R_HOME}" >> "${BENCHMARK_ENV_FILE}"
    if [[ -n "${LD_PRELOAD:-}" ]]; then
        printf 'export LD_PRELOAD=%q\n' "${LD_PRELOAD}" >> "${BENCHMARK_ENV_FILE}"
    fi
fi
