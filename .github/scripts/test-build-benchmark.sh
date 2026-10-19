#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_SCRIPT="${BUILD_SCRIPT:-${SCRIPT_DIR}/build_benchmark.sh}"
FIXTURE="$(mktemp -d)"
trap 'rm -rf "${FIXTURE}"' EXIT

mkdir -p "${FIXTURE}/benchmarks/Example" "${FIXTURE}/benchmarks/Default" "${FIXTURE}/bin"
mkdir -p "${FIXTURE}/.github/scripts"
cp "${SCRIPT_DIR}/read-benchmark-config.sh" "${FIXTURE}/.github/scripts/"
printf 'runner = ["self-hosted", "benchmark"]\ntimeout = 12000\nthreads = 32\n' \
    > "${FIXTURE}/.github/benchmark_defaults.toml"
printf 'threads = 128\n' > "${FIXTURE}/benchmarks/Example/benchmark_config.toml"
(
    cd "${FIXTURE}"
    source .github/scripts/read-benchmark-config.sh
    if [[ "$(get_runner_config benchmarks/Example)" != '["self-hosted", "benchmark"]|12000' ]]; then
        echo "get_runner_config should return runner and timeout only" >&2
        exit 1
    fi
)
printf 'julia_version = "1.11.9"\n[[deps.StyledStrings]]\n' > "${FIXTURE}/Manifest.toml"
printf '%s\n' \
    '#!/bin/bash' \
    'set -euo pipefail' \
    "printf \"%s\\n\" \"\$*\" >> \"\${FAKE_JULIA_LOG}\"" \
    "if [[ \"\$*\" == *\"Pkg.instantiate\"* && -f Manifest.toml ]]; then" \
    '    echo "ERROR: Could not locate the source code for the StyledStrings package." >&2' \
    '    exit 1' \
    'fi' \
    "if [[ \"\$*\" == *\"Pkg.instantiate\"* ]]; then" \
    "    printf 'julia_version = \"1.10.11\"\\n' > Manifest.toml" \
    'fi' > "${FIXTURE}/bin/julia"
chmod +x "${FIXTURE}/bin/julia"

(
    cd "${FIXTURE}"
    PATH="${FIXTURE}/bin:${PATH}" FAKE_JULIA_LOG="${FIXTURE}/julia.log" \
        "${BUILD_SCRIPT}" benchmarks/Default
    PATH="${FIXTURE}/bin:${PATH}" FAKE_JULIA_LOG="${FIXTURE}/julia.log" \
        "${BUILD_SCRIPT}" benchmarks/Example
)

grep -Fq 'julia_version = "1.10.11"' "${FIXTURE}/Manifest.toml"
if grep -Fq 'StyledStrings' "${FIXTURE}/Manifest.toml"; then
    exit 1
fi
test "$(wc -l < "${FIXTURE}/julia.log")" -eq 4
grep -Fq -- '--project=. -e using Pkg; Pkg.instantiate(); Pkg.build()' \
    "${FIXTURE}/julia.log"
grep -Fq -- '--threads=32 --project=. benchmark.jl benchmarks/Default' \
    "${FIXTURE}/julia.log"
grep -Fq -- '--threads=128 --project=. benchmark.jl benchmarks/Example' \
    "${FIXTURE}/julia.log"

(
    cd "${SCRIPT_DIR}/../.."
    source .github/scripts/read-benchmark-config.sh
    failed=0
    for page in benchmarks/AstroChem/{astrochem,nelson}.jmd \
        benchmarks/Bio/Bidkhori2012.jmd \
        benchmarks/BayesianInference/DiffEqBayesLorenz.jmd \
        benchmarks/ComplicatedPDE/Filament.jmd \
        benchmarks/MethodOfLinesPDE/{MOL_fdm,MOLxPDESystemLibrary}.jmd \
        benchmarks/NonStiffBVP/{linear_wpd,nonlinear_wpd}.jmd \
        benchmarks/ParameterEstimation/BioPreDynB{1..6}ParameterEstimation.jmd; do
        test -f "${page}"
        for target in "${page%/*}" "${page}"; do
            config=$(get_benchmark_config "${target}")
            if [[ "${config##*|}" != 128 ]]; then
                echo "${target}: expected 128 threads, got ${config##*|}" >&2
                failed=1
            fi
        done
    done
    test "${failed}" -eq 0
)

echo "build_benchmark.sh tests passed"
