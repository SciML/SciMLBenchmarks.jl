#!/bin/bash
set -euo pipefail

# Put omc on PATH for OMJulia.OMCSession() in RCCircuit/ThermalFluid.
# Official install is apt `omc` from the OpenModelica repo (needs passwordless
# sudo, as on the self-hosted runners). This machine and any runner without
# sudo fall back to a user-local conda-forge prefix under $HOME.

CACHE_DIR="${HOME}/.cache/sciml-benchmarks"
CONDA_PREFIX="${CACHE_DIR}/openmodelica"

export_omc_env() {
    local omc_bin omc_dir prefix
    omc_bin="$(command -v omc)"
    omc_dir="$(cd "$(dirname "${omc_bin}")" && pwd)"
    prefix="$(cd "${omc_dir}/.." && pwd)"
    echo "--- omc: ${omc_bin} ($(omc --version 2>/dev/null | head -1 || true))"
    if [[ -z "${BENCHMARK_ENV_FILE:-}" ]]; then
        return 0
    fi
    {
        echo "export PATH=\"${omc_dir}:\${PATH}\""
        if [[ -d "${prefix}/conda-meta" ]]; then
            # omc's RPATH already includes $ORIGIN/../lib; do not export
            # LD_LIBRARY_PATH (it breaks Julia/CairoMakie against conda libs).
            echo "export OPENMODELICAHOME=\"${prefix}\""
        fi
    } >> "${BENCHMARK_ENV_FILE}"
}

install_omc_apt() {
    sudo -n true 2>/dev/null || return 1
    echo "--- Installing omc via the OpenModelica apt repository"
    export DEBIAN_FRONTEND=noninteractive
    sudo -n apt-get update -qq
    sudo -n apt-get install -y -qq ca-certificates curl gnupg
    if [[ ! -f /usr/share/keyrings/openmodelica-keyring.gpg ]]; then
        curl -fsSL https://build.openmodelica.org/apt/openmodelica.asc |
            sudo -n gpg --dearmor -o /usr/share/keyrings/openmodelica-keyring.gpg
    fi
    local codename
    # shellcheck disable=SC1091
    codename="$(. /etc/os-release && printf '%s\n' "${UBUNTU_CODENAME:-${VERSION_CODENAME:-}}")"
    if [[ -z "${codename}" ]]; then
        echo "--- could not detect Debian/Ubuntu codename" >&2
        return 1
    fi
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/openmodelica-keyring.gpg] https://build.openmodelica.org/apt ${codename} release" |
        sudo -n tee /etc/apt/sources.list.d/openmodelica.list >/dev/null
    sudo -n apt-get update -qq
    sudo -n apt-get install -y --no-install-recommends omc
    hash -r
    command -v omc >/dev/null 2>&1
}

micromamba_platform() {
    case "$(uname -s)-$(uname -m)" in
        Linux-x86_64) echo linux-64 ;;
        Linux-aarch64) echo linux-aarch64 ;;
        Darwin-x86_64) echo osx-64 ;;
        Darwin-arm64) echo osx-arm64 ;;
        *)
            echo "unsupported $(uname -s) $(uname -m) for conda-forge openmodelica" >&2
            return 1
            ;;
    esac
}

install_omc_conda() {
    mkdir -p "${CACHE_DIR}"
    if [[ -x "${CONDA_PREFIX}/bin/omc" ]]; then
        echo "--- reusing conda-forge omc at ${CONDA_PREFIX}"
        export PATH="${CONDA_PREFIX}/bin:${PATH}"
        export OPENMODELICAHOME="${CONDA_PREFIX}"
        return 0
    fi
    local plat micromamba
    plat="$(micromamba_platform)"
    micromamba="${CACHE_DIR}/micromamba"
    if [[ ! -x "${micromamba}" ]]; then
        echo "--- bootstrapping micromamba (${plat})"
        curl -fsSL "https://micro.mamba.pm/api/micromamba/${plat}/latest" |
            tar -xj -C "${CACHE_DIR}" bin/micromamba
        mv "${CACHE_DIR}/bin/micromamba" "${micromamba}"
        rmdir "${CACHE_DIR}/bin" 2>/dev/null || true
        chmod +x "${micromamba}"
    fi
    echo "--- Installing openmodelica from conda-forge into ${CONDA_PREFIX}"
    export MAMBA_ROOT_PREFIX="${CACHE_DIR}/mamba-root"
    mkdir -p "${MAMBA_ROOT_PREFIX}"
    "${micromamba}" create -y -p "${CONDA_PREFIX}" -c conda-forge --override-channels openmodelica
    export PATH="${CONDA_PREFIX}/bin:${PATH}"
    export OPENMODELICAHOME="${CONDA_PREFIX}"
    hash -r
}

if command -v omc >/dev/null 2>&1; then
    export_omc_env
    exit 0
fi

if ! install_omc_apt; then
    echo "--- apt omc install unavailable; falling back to conda-forge"
    install_omc_conda
fi

if ! command -v omc >/dev/null 2>&1; then
    echo "ERROR: omc is not on PATH after ModelingToolkit setup" >&2
    exit 1
fi
export_omc_env
