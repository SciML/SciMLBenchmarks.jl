#!/bin/bash

set -euo pipefail

# Instantiate, to install the overall project dependencies, and `build()` for conda
echo "--- Instantiate"
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
echo "+++ Run benchmark for ${1}"
julia --threads=auto --project=. benchmark.jl "${1}"

