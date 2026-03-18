#!/bin/bash
set -euo pipefail

# Install CmdStan for BayesianInference benchmarks
CMDSTAN_VERSION="2.38.0"

export CMDSTAN_HOME="$(pwd)/cmdstan-${CMDSTAN_VERSION}/"
echo "--- Installing CmdStan ${CMDSTAN_VERSION}"

if [[ ! -d "cmdstan-${CMDSTAN_VERSION}" ]]; then
    curl -LO "https://github.com/stan-dev/cmdstan/releases/download/v${CMDSTAN_VERSION}/cmdstan-${CMDSTAN_VERSION}.tar.gz"
    tar -xzpf "cmdstan-${CMDSTAN_VERSION}.tar.gz" --no-same-owner
    rm -f "cmdstan-${CMDSTAN_VERSION}.tar.gz"
fi

echo "STAN_THREADS=true" > "./cmdstan-${CMDSTAN_VERSION}/make/local"
echo "g++ version:"
g++ --version
make -C "cmdstan-${CMDSTAN_VERSION}" build

# Export CMDSTAN_HOME for the benchmark Julia process
echo "export CMDSTAN_HOME=\"${CMDSTAN_HOME}\"" >> "${BENCHMARK_ENV_FILE:-/dev/null}"
