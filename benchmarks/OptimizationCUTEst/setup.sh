#!/bin/bash
set -euo pipefail

BENCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX="${BENCH_DIR}/.gfortran"
MICROMAMBA="${BENCH_DIR}/.micromamba"

if ! command -v gfortran >/dev/null 2>&1; then
    if [[ -x "${PREFIX}/bin/gfortran" ]]; then
        export PATH="${PREFIX}/bin:${PATH}"
    else
        case "$(uname -s)-$(uname -m)" in
            Linux-x86_64) platform=linux-64 ;;
            Linux-aarch64) platform=linux-aarch64 ;;
            *)
                echo "ERROR: gfortran is required for OptimizationCUTEst; no conda-forge installer is configured for $(uname -s) $(uname -m)." >&2
                exit 1
                ;;
        esac

        mkdir -p "${BENCH_DIR}/.mamba-root"
        if [[ ! -x "${MICROMAMBA}" ]]; then
            echo "--- Bootstrapping micromamba for the CUTEst Fortran compiler"
            curl -fsSL "https://micro.mamba.pm/api/micromamba/${platform}/latest" |
                tar -xj -C "${BENCH_DIR}" bin/micromamba
            mv "${BENCH_DIR}/bin/micromamba" "${MICROMAMBA}"
            rmdir "${BENCH_DIR}/bin"
            chmod +x "${MICROMAMBA}"
        fi

        echo "--- Installing gfortran from conda-forge into ${PREFIX}"
        export MAMBA_ROOT_PREFIX="${BENCH_DIR}/.mamba-root"
        "${MICROMAMBA}" create -y -p "${PREFIX}" -c conda-forge --override-channels gfortran
        export PATH="${PREFIX}/bin:${PATH}"
    fi
fi

if ! command -v gfortran >/dev/null 2>&1; then
    echo "ERROR: gfortran is not on PATH after OptimizationCUTEst setup" >&2
    exit 1
fi
echo "--- gfortran: $(command -v gfortran) ($(gfortran -dumpmachine))"

if [[ -n "${BENCHMARK_ENV_FILE:-}" && -x "${PREFIX}/bin/gfortran" ]]; then
    echo "export PATH=\"${PREFIX}/bin:\${PATH}\"" >> "${BENCHMARK_ENV_FILE}"
fi
