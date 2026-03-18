#!/bin/bash
set -euo pipefail

JULIAHUBREGISTRY_BENCHMARK_TARGETS=(benchmarks/ModelingToolkit/)
OPENMODELICA_BENCHMARK_TARGETS=(benchmarks/ModelingToolkit/)

if [[ "${JULIAHUBREGISTRY_BENCHMARK_TARGETS[*]}" =~ "${1}" ]]; then
	echo "--- Adding JuliaHubRegistry"

	export JULIA_PKG_SERVER="juliahub.com"
	mkdir -p "${HOME}/.julia/servers/${JULIA_PKG_SERVER}"
	if [[ -n "${JULIAHUB_TOKEN:-}" ]]; then
		echo "${JULIAHUB_TOKEN}" > "${HOME}/.julia/servers/${JULIA_PKG_SERVER}/auth.toml"
	fi
	julia -e 'using Pkg; Pkg.Registry.add(); Pkg.Registry.status()'
fi

# Instantiate, to install the overall project dependencies, and build() for conda
echo "--- Instantiate"
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.build()'

if [[ "${1}" == *BayesianInference* ]]; then
	export CMDSTAN_HOME="$(pwd)/cmdstan-2.38.0/"
	curl -LO https://github.com/stan-dev/cmdstan/releases/download/v2.38.0/cmdstan-2.38.0.tar.gz
	tar -xzpf cmdstan-2.38.0.tar.gz --no-same-owner
	echo "STAN_THREADS=true" > ./cmdstan-2.38.0/make/local
	echo "g++ version"
	g++ --version
	make -C cmdstan-2.38.0 build
fi

# Run benchmark
echo "+++ Run benchmark for ${1}"
julia --threads=auto --project=. benchmark.jl "${1}"
