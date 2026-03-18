#!/bin/bash
set -euo pipefail

# Add JuliaHubRegistry for ModelingToolkit benchmarks
echo "--- Adding JuliaHubRegistry"

export JULIA_PKG_SERVER="juliahub.com"
mkdir -p "${HOME}/.julia/servers/${JULIA_PKG_SERVER}"
if [[ -n "${JULIAHUB_TOKEN:-}" ]]; then
    echo "${JULIAHUB_TOKEN}" > "${HOME}/.julia/servers/${JULIA_PKG_SERVER}/auth.toml"
fi
julia -e 'using Pkg; Pkg.Registry.add(); Pkg.Registry.status()'
