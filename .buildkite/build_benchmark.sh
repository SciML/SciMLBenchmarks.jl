#!/bin/bash

set -euo pipefail

JULIAHUBREGISTRY_BENCHMARK_TARGETS=(benchmarks/ModelingToolkit/)
OPENMODELICA_BENCHMARK_TARGETS=(benchmarks/ModelingToolkit/)

if [[ "${JULIAHUBREGISTRY_BENCHMARK_TARGETS[*]}" =~ "${1}" ]]; then
	echo "--- :julia: Adding JuliaHubRegistry"

	export JULIA_PKG_SERVER="juliahub.com"
	mkdir -p "${JULIA_DEPOT_PATH}/servers/${JULIA_PKG_SERVER}"
	cp .buildkite/secrets/token.toml "${JULIA_DEPOT_PATH}/servers/${JULIA_PKG_SERVER}/auth.toml"
	julia -e 'using Pkg; Pkg.Registry.add(); Pkg.Registry.status()'
fi

# Instantiate, to install the overall project dependencies, and `build()` for conda
echo "--- :julia: Instantiate"
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.build()'

if [[ "${1}" == *BayesianInference* ]]; then
	export CMDSTAN_HOME="$(pwd)/cmdstan-2.29.2/"
	curl -LO https://github.com/stan-dev/cmdstan/releases/download/v2.29.2/cmdstan-2.29.2.tar.gz
	tar -xzpf cmdstan-2.29.2.tar.gz --no-same-owner
	echo "STAN_THREADS=true" > ./cmdstan-2.29.2/make/local
	echo "g++ version"
	echo $(g++ -v)
	make -C cmdstan-2.29.2 build
fi

# Run benchmark
echo "+++ :julia: Run benchmark for ${1}"
julia --threads=auto --project=. benchmark.jl "${1}"

