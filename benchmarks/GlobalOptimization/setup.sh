#!/bin/bash
set -euo pipefail

# OptimizationSciPy / OptimizationPyCMA need numpy+scipy+cma. CondaPkg's pixi
# backend on this host produced a 0-byte libffi (Python ctypes ImportError).
# Install a user-local conda-forge prefix and tell PythonCall to use it.

CACHE_DIR="${HOME}/.cache/sciml-benchmarks"
PREFIX="${CACHE_DIR}/globalopt-python"
BENCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

is_windows() {
    case "$(uname -s)" in
        MINGW* | MSYS* | CYGWIN*) return 0 ;;
        *) return 1 ;;
    esac
}

# micromamba.exe and python.exe are native Windows binaries: values handed to
# them (env vars, -p prefix) must be Windows paths even though this script runs
# under Git Bash, and JULIA_PYTHONCALL_EXE must be one for Julia to load it.
native_path() {
    if is_windows; then
        if command -v cygpath >/dev/null 2>&1; then
            cygpath -w "$1"
        else
            echo "cygpath not found; cannot convert $1 to a Windows path" >&2
            return 1
        fi
    else
        printf '%s\n' "$1"
    fi
}

micromamba_platform() {
    case "$(uname -s)-$(uname -m)" in
        Linux-x86_64) echo linux-64 ;;
        Linux-aarch64) echo linux-aarch64 ;;
        MINGW*-x86_64 | MSYS*-x86_64 | CYGWIN*-x86_64) echo win-64 ;;
        MINGW*-aarch64 | MSYS*-aarch64 | CYGWIN*-aarch64) echo win-arm64 ;;
        *)
            echo "unsupported $(uname -s) $(uname -m) for conda-forge python" >&2
            return 1
            ;;
    esac
}

if is_windows; then
    # win-64 tarballs ship micromamba under Library/bin and conda-forge puts
    # python.exe at the env root, not bin/python.
    MICROMAMBA="${CACHE_DIR}/micromamba.exe"
    TAR_MEMBER="Library/bin/micromamba.exe"
    PYTHON="${PREFIX}/python.exe"
else
    MICROMAMBA="${CACHE_DIR}/micromamba"
    TAR_MEMBER="bin/micromamba"
    PYTHON="${PREFIX}/bin/python"
fi

if [[ ! -x "${MICROMAMBA}" ]]; then
    mkdir -p "${CACHE_DIR}"
    plat="$(micromamba_platform)"
    curl -fsSL "https://micro.mamba.pm/api/micromamba/${plat}/latest" |
        tar -xj -C "${CACHE_DIR}" "${TAR_MEMBER}"
    mv "${CACHE_DIR}/${TAR_MEMBER}" "${MICROMAMBA}"
    if is_windows; then
        rmdir "${CACHE_DIR}/Library/bin" "${CACHE_DIR}/Library" 2>/dev/null || true
    else
        rmdir "${CACHE_DIR}/bin" 2>/dev/null || true
    fi
    chmod +x "${MICROMAMBA}"
fi

if [[ ! -x "${PYTHON}" ]]; then
    echo "--- Installing python/numpy/scipy/matplotlib/cma into ${PREFIX}"
    # Standalone assignments so set -e aborts if a conversion fails.
    NATIVE_MAMBA_ROOT="$(native_path "${CACHE_DIR}/mamba-root")"
    NATIVE_PREFIX="$(native_path "${PREFIX}")"
    export MAMBA_ROOT_PREFIX="${NATIVE_MAMBA_ROOT}"
    mkdir -p "${CACHE_DIR}/mamba-root"
    "${MICROMAMBA}" create -y -p "${NATIVE_PREFIX}" \
        -c conda-forge --override-channels \
        python=3.12 numpy scipy matplotlib cma
fi

"${PYTHON}" -c "import ctypes, numpy, scipy, cma"

# Drop a half-built pixi env so PythonCall cannot pick it up.
rm -rf "${BENCH_DIR}/.CondaPkg"

if [[ -n "${BENCHMARK_ENV_FILE:-}" ]]; then
    NATIVE_PYTHON="$(native_path "${PYTHON}")"
    {
        echo 'export JULIA_LOAD_PATH="@:@stdlib"'
        echo "export JULIA_PYTHONCALL_EXE=\"${NATIVE_PYTHON}\""
        echo 'export JULIA_CONDAPKG_BACKEND="Null"'
        if is_windows; then
            # python.exe resolves its DLLs from the env root and Library/bin.
            # Entries are POSIX-style here; MSYS rewrites PATH into native
            # form when it spawns julia.exe.
            echo "export PATH=\"${PREFIX}:${PREFIX}/Library/bin:${PREFIX}/Scripts:\${PATH}\""
        fi
    } >> "${BENCHMARK_ENV_FILE}"
fi
