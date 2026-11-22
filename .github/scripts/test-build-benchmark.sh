#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/read-benchmark-config.sh"
BUILD_SCRIPT="${BUILD_SCRIPT:-${SCRIPT_DIR}/build_benchmark.sh}"
FIXTURE="$(mktemp -d)"
trap 'rm -rf "${FIXTURE}"' EXIT

mkdir -p "${FIXTURE}/benchmarks/Example" "${FIXTURE}/benchmarks/Default" "${FIXTURE}/bin"
mkdir -p "${FIXTURE}/.github"
printf 'runner = ["self-hosted", "benchmark"]\ntimeout = 12000\nthreads = 32\n' \
    > "${FIXTURE}/.github/benchmark_defaults.toml"
printf 'threads = 128  # page and folder jobs\n' > "${FIXTURE}/benchmarks/Example/benchmark_config.toml"
(
    cd "${FIXTURE}"
    if [[ "$(get_runner_config benchmarks/Example)" != '["self-hosted", "benchmark"]|12000' ]]; then
        echo "get_runner_config should return runner and timeout only" >&2
        exit 1
    fi
)
printf 'julia_version = "1.11.9"\n[[deps.StyledStrings]]\n' > "${FIXTURE}/Manifest.toml"
cat > "${FIXTURE}/bin/julia" <<'JULIA'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${FAKE_JULIA_LOG}"
printf '%s|%s\n' "$JULIA_NUM_THREADS" "$SCIMLBENCHMARKS_BUILD_SENTINEL" >> "$FAKE_JULIA_ENV_LOG"
if [[ "$*" == *"Pkg.instantiate"* && -f Manifest.toml ]]; then
    echo "ERROR: Could not locate the source code for the StyledStrings package." >&2
    exit 1
fi
if [[ "$*" == *"Pkg.instantiate"* ]]; then
    printf 'julia_version = "1.10.11"\n' > Manifest.toml
fi
JULIA
chmod +x "${FIXTURE}/bin/julia"

(
    cd "${FIXTURE}"
    export JULIA_NUM_THREADS=3 SCIMLBENCHMARKS_BUILD_SENTINEL=kept
    export FAKE_JULIA_ENV_LOG="${FIXTURE}/environment.log"
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
grep -F 'benchmark.jl' "${FIXTURE}/julia.log" > "${FIXTURE}/benchmark-commands"
printf '%s\n' '--threads=32 --project=. benchmark.jl benchmarks/Default' \
    '--threads=128 --project=. benchmark.jl benchmarks/Example' > "${FIXTURE}/expected-commands"
diff -u "${FIXTURE}/expected-commands" "${FIXTURE}/benchmark-commands"
printf '%s\n' '32|kept' '32|kept' '128|kept' '128|kept' > "${FIXTURE}/expected-environment.log"
diff -u "${FIXTURE}/expected-environment.log" "${FIXTURE}/environment.log"

(
    cd "${SCRIPT_DIR}/../.."
    shopt -s nullglob dotglob
    for project in benchmarks/*/Project.toml; do
        folder=${project%/*}
        case "${folder##*/}" in
            AdaptiveSDE|AstroChem|BayesianInference|Bio|ComplicatedPDE|GlobalOptimization|\
            Jumps|LinearSolve|MethodOfLinesPDE|NonStiffBVP|NonStiffODE|NonStiffSDE|\
            NonlinearProblem|Optimization|ParameterEstimation|StiffODE|StiffSDE|Testing)
                expected=128 ;;
            DiffEqGPU|LinearSolveGPU)
                expected=58 ;;
            *) expected=32 ;;
        esac
        for target in "${folder}" "${folder}"/*.jmd; do
            config=$(get_benchmark_config "${target}")
            if [[ "${config##*|}" != "${expected}" ]]; then
                echo "${target}: expected ${expected} threads, got ${config##*|}" >&2
                exit 1
            fi
        done
    done
)

(
    cd "${FIXTURE}"
    printf 'threads = 7\n' > .github/benchmark_defaults.toml
    config=$(get_benchmark_config benchmarks/Default)
    [[ "${config##*|}" == 7 ]]
    for invalid in 0 -1 auto 2,2; do
        printf 'threads = %s\n' "${invalid}" > benchmarks/Default/benchmark_config.toml
        if get_benchmark_config benchmarks/Default; then
            echo "Accepted invalid thread count: ${invalid}" >&2
            exit 1
        fi
    done
)

echo "build_benchmark.sh tests passed"

printf '[deps]\n' > "${FIXTURE}/benchmarks/Example/Project.toml"
printf 'julia_version = "1.11.9"\n' > "${FIXTURE}/benchmarks/Example/Manifest.toml"
cp "${FIXTURE}/benchmarks/Example/Project.toml" "${FIXTURE}/expected-project"
cp "${FIXTURE}/benchmarks/Example/Manifest.toml" "${FIXTURE}/expected-manifest"
cat > "${FIXTURE}/benchmarks/Example/setup.sh" <<'SETUP'
#!/bin/bash
set -euo pipefail
[[ "${JULIA_NUM_THREADS}" == 128 ]]
printf 'export SHARED_SETUP=ready\n' >> "${BENCHMARK_ENV_FILE}"
printf 'setup\n' >> "${FAKE_SETUP_LOG}"
SETUP
chmod +x "${FIXTURE}/benchmarks/Example/setup.sh"
cat > "${FIXTURE}/bin/julia" <<'JULIA'
#!/bin/bash
set -euo pipefail
[[ "${SHARED_SETUP}" == ready ]]
printf '%s\n' "$*" >> "${FAKE_JULIA_LOG}"
JULIA
chmod +x "${FIXTURE}/bin/julia"

for target in benchmarks/Example benchmarks/Example/page.jmd; do
    : > "${FIXTURE}/julia.log"
    : > "${FIXTURE}/setup.log"
    (
        cd "${FIXTURE}"
        PATH="${FIXTURE}/bin:${PATH}" FAKE_JULIA_LOG="${FIXTURE}/julia.log" \
            FAKE_SETUP_LOG="${FIXTURE}/setup.log" "${BUILD_SCRIPT}" "${target}"
    )
    test "$(wc -l < "${FIXTURE}/setup.log")" -eq 1
    test "$(wc -l < "${FIXTURE}/julia.log")" -eq 2
    grep -Fxq -- '--project=. -e using Pkg; Pkg.instantiate(); Pkg.build()' "${FIXTURE}/julia.log"
    grep -Fxq -- "--threads=128 --project=. benchmark.jl ${target}" "${FIXTURE}/julia.log"
    cmp "${FIXTURE}/expected-project" "${FIXTURE}/benchmarks/Example/Project.toml"
    cmp "${FIXTURE}/expected-manifest" "${FIXTURE}/benchmarks/Example/Manifest.toml"
done

echo 'Folder and file builds both run shared setup, propagate its environment, and preserve Project/Manifest.'
