#!/bin/bash
set -euo pipefail

# OptimizationSciPy / OptimizationPyCMA need numpy+scipy+cma. CondaPkg's pixi
# backend on this host produced a 0-byte libffi (Python ctypes ImportError).
# Install a user-local conda-forge prefix and tell PythonCall to use it.

CACHE_DIR="${HOME}/.cache/sciml-benchmarks"
PREFIX="${CACHE_DIR}/globalopt-python"
MICROMAMBA="${CACHE_DIR}/micromamba"
BENCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

micromamba_platform() {
    case "$(uname -s)-$(uname -m)" in
        Linux-x86_64) echo linux-64 ;;
        Linux-aarch64) echo linux-aarch64 ;;
        *)
            echo "unsupported $(uname -s) $(uname -m) for conda-forge python" >&2
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

if [[ ! -x "${PREFIX}/bin/python" ]]; then
    echo "--- Installing python/numpy/scipy/matplotlib/cma into ${PREFIX}"
    export MAMBA_ROOT_PREFIX="${CACHE_DIR}/mamba-root"
    mkdir -p "${MAMBA_ROOT_PREFIX}"
    "${MICROMAMBA}" create -y -p "${PREFIX}" -c conda-forge --override-channels \
        python=3.12 numpy scipy matplotlib cma
fi

"${PREFIX}/bin/python" -c "import ctypes, numpy, scipy, cma"

# Drop a half-built pixi env so PythonCall cannot pick it up.
rm -rf "${BENCH_DIR}/.CondaPkg"

if [[ -n "${BENCHMARK_ENV_FILE:-}" ]]; then
    {
        echo 'export JULIA_LOAD_PATH="@:@stdlib"'
        echo "export JULIA_PYTHONCALL_EXE=\"${PREFIX}/bin/python\""
        echo 'export JULIA_CONDAPKG_BACKEND="Null"'
    } >> "${BENCHMARK_ENV_FILE}"
fi
