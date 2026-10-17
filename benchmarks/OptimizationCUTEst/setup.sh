#!/bin/bash
set -euo pipefail

# CUTEst.jl requires a Fortran compiler: __init__ errors when `gfortran` is not
# on PATH and decoded SIF problems are compiled with `gfortran -O3 -fPIC` then
# linked with `gfortran -shared`. The self-hosted runners have no system
# gfortran, so install a user-local conda-forge toolchain and export it for the
# benchmark Julia process.

CACHE_DIR="${HOME}/.cache/sciml-benchmarks"
PREFIX="${CACHE_DIR}/cutest-gfortran"
MICROMAMBA="${CACHE_DIR}/micromamba"

micromamba_platform() {
    case "$(uname -s)-$(uname -m)" in
        Linux-x86_64) echo linux-64 ;;
        Linux-aarch64) echo linux-aarch64 ;;
        *)
            echo "unsupported $(uname -s) $(uname -m) for conda-forge gfortran" >&2
            return 1
            ;;
    esac
}

if [[ ! -x "${MICROMAMBA}" ]]; then
    mkdir -p "${CACHE_DIR}"
    plat="$(micromamba_platform)"
    curl -fsSL "https://micro.mamba.pm/api/micromamba/${plat}/latest" |
        tar -xj -C "${CACHE_DIR}" bin/micromamba
    mv "${CACHE_DIR}/bin/micromamba" "${MICROMAMBA}"
    rmdir "${CACHE_DIR}/bin" 2>/dev/null || true
    chmod +x "${MICROMAMBA}"
fi

if [[ ! -x "${PREFIX}/bin/gfortran" ]]; then
    echo "--- Installing gfortran into ${PREFIX}"
    export MAMBA_ROOT_PREFIX="${CACHE_DIR}/mamba-root"
    mkdir -p "${MAMBA_ROOT_PREFIX}"
    "${MICROMAMBA}" create -y -p "${PREFIX}" -c conda-forge --override-channels \
        gfortran
fi

"${PREFIX}/bin/gfortran" --version | head -1
echo "gfortran target: $("${PREFIX}/bin/gfortran" -dumpmachine)"

if [[ -n "${BENCHMARK_ENV_FILE:-}" ]]; then
    {
        echo "export PATH=\"${PREFIX}/bin:\${PATH}\""
        # Libraries produced by `gfortran -shared` record NEEDED libgfortran
        # without an rpath; make the conda's runtime lib dir visible to dlopen.
        echo "export LD_LIBRARY_PATH=\"${PREFIX}/lib:\${LD_LIBRARY_PATH:-}\""
    } >> "${BENCHMARK_ENV_FILE}"
fi
