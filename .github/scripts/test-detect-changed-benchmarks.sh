#!/bin/bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TEST_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/detect-changed-benchmarks.XXXXXX")
trap 'rm -rf "${TEST_ROOT}"' EXIT

FIXTURE="${TEST_ROOT}/fixture"
mkdir -p "${FIXTURE}/.github/scripts" "${FIXTURE}/benchmarks/Alpha" \
    "${FIXTURE}/benchmarks/Beta" "${TEST_ROOT}/bin"
cp "${REPO_ROOT}/.github/scripts/detect-changed-benchmarks.sh" \
    "${REPO_ROOT}/.github/scripts/read-benchmark-config.sh" "${FIXTURE}/.github/scripts/"
cp "${REPO_ROOT}/.github/benchmark_defaults.toml" "${FIXTURE}/.github/"

printf '[deps]\n' > "${FIXTURE}/benchmarks/Alpha/Project.toml"
printf 'julia_version = "1.11.9"\n' > "${FIXTURE}/benchmarks/Alpha/Manifest.toml"
printf '[deps]\n' > "${FIXTURE}/benchmarks/Beta/Project.toml"
printf 'julia_version = "1.11.9"\n' > "${FIXTURE}/benchmarks/Beta/Manifest.toml"

git -C "${FIXTURE}" init -q -b master
git -C "${FIXTURE}" config user.name "Detector Test"
git -C "${FIXTURE}" config user.email "detector-test@example.com"
git -C "${FIXTURE}" add .
git -C "${FIXTURE}" commit -qm "Initial fixture"
INITIAL_SHA=$(git -C "${FIXTURE}" rev-parse HEAD)

printf 'version = "0.1.0"\n' >> "${FIXTURE}/benchmarks/Alpha/Project.toml"
git -C "${FIXTURE}" add benchmarks/Alpha/Project.toml
git -C "${FIXTURE}" commit -qm "Update Alpha"
PUSH_BASE_SHA=$(git -C "${FIXTURE}" rev-parse HEAD)

printf 'version = "0.1.0"\n' >> "${FIXTURE}/benchmarks/Beta/Project.toml"
git -C "${FIXTURE}" add benchmarks/Beta/Project.toml
git -C "${FIXTURE}" commit -qm "Update Beta"

printf '#!/bin/sh\nprintf '\''[{"commit":{"message":"Published by build of: SciML/SciMLBenchmarks.jl@%%s"}}]\\n'\'' "%s"\n' \
    "${INITIAL_SHA}" > "${TEST_ROOT}/bin/curl"
chmod +x "${TEST_ROOT}/bin/curl"

OUTPUT_FILE="${TEST_ROOT}/github-output"
(
    cd "${FIXTURE}"
    GITHUB_EVENT_NAME=push \
        PUSH_BASE_SHA="${PUSH_BASE_SHA}" \
        GITHUB_OUTPUT="${OUTPUT_FILE}" \
        PATH="${TEST_ROOT}/bin:${PATH}" \
        .github/scripts/detect-changed-benchmarks.sh
)

TARGETS=$(sed -n 's/^matrix=//p' "${OUTPUT_FILE}" \
    | grep -o '"target":"[^"]*"' \
    | cut -d'"' -f4 \
    | sort)

if [[ "${TARGETS}" != "benchmarks/Beta" ]]; then
    printf 'Expected only benchmarks/Beta for the second push, got:\n%s\n' "${TARGETS}" >&2
    exit 1
fi

SUPPORT_BASE_SHA=$(git -C "${FIXTURE}" rev-parse HEAD)
printf '#!/bin/sh\n' > "${FIXTURE}/benchmarks/Alpha/setup.sh"
git -C "${FIXTURE}" add benchmarks/Alpha/setup.sh
git -C "${FIXTURE}" commit -qm "Update Alpha setup"

SUPPORT_OUTPUT_FILE="${TEST_ROOT}/support-github-output"
(
    cd "${FIXTURE}"
    GITHUB_EVENT_NAME=push \
        PUSH_BASE_SHA="${SUPPORT_BASE_SHA}" \
        GITHUB_OUTPUT="${SUPPORT_OUTPUT_FILE}" \
        PATH="${TEST_ROOT}/bin:${PATH}" \
        .github/scripts/detect-changed-benchmarks.sh
)

SUPPORT_TARGETS=$(sed -n 's/^matrix=//p' "${SUPPORT_OUTPUT_FILE}" \
    | grep -o '"target":"[^"]*"' \
    | cut -d'"' -f4 \
    | sort || true)

if [[ "${SUPPORT_TARGETS}" != "benchmarks/Alpha" ]]; then
    printf 'Expected benchmarks/Alpha for a support-file change, got:\n%s\n' \
        "${SUPPORT_TARGETS}" >&2
    exit 1
fi

printf '# Deleted benchmark\n' > "${FIXTURE}/benchmarks/Alpha/deleted.jmd"
git -C "${FIXTURE}" add benchmarks/Alpha/deleted.jmd
git -C "${FIXTURE}" commit -qm "Add benchmark to delete"
DELETION_BASE_SHA=$(git -C "${FIXTURE}" rev-parse HEAD)
git -C "${FIXTURE}" rm -q benchmarks/Alpha/deleted.jmd
git -C "${FIXTURE}" commit -qm "Delete benchmark"

DELETION_OUTPUT_FILE="${TEST_ROOT}/deletion-github-output"
(
    cd "${FIXTURE}"
    GITHUB_EVENT_NAME=push \
        PUSH_BASE_SHA="${DELETION_BASE_SHA}" \
        GITHUB_OUTPUT="${DELETION_OUTPUT_FILE}" \
        PATH="${TEST_ROOT}/bin:${PATH}" \
        .github/scripts/detect-changed-benchmarks.sh
)

if ! grep -Fxq 'has_changes=false' "${DELETION_OUTPUT_FILE}"; then
    printf 'Expected a deleted .jmd file to schedule no benchmark, got:\n' >&2
    sed -n 's/^matrix=//p' "${DELETION_OUTPUT_FILE}" >&2
    exit 1
fi

echo "Push detection selected only changes since the event's before SHA."

assert_detection() {
    local expected="$1"
    shift
    local base_sha output_file targets
    base_sha=$(git -C "${FIXTURE}" rev-parse HEAD)
    git -C "${FIXTURE}" add .
    git -C "${FIXTURE}" commit -qm "Update detection fixture" --allow-empty
    output_file="${TEST_ROOT}/case-github-output"
    : > "${output_file}"
    (
        cd "${FIXTURE}"
        GITHUB_EVENT_NAME=push PUSH_BASE_SHA="${base_sha}" \
            GITHUB_OUTPUT="${output_file}" \
            .github/scripts/detect-changed-benchmarks.sh "$@"
    )
    sed -n 's/^matrix=/matrix=/p' "${output_file}"
    targets=$(sed -n 's/^matrix=//p' "${output_file}" \
        | grep -o '"target":"[^"]*"' | cut -d'"' -f4 | sort || true)
    if [[ "${targets}" != "${expected}" ]]; then
        printf 'Expected targets:\n%s\nGot:\n%s\n' "${expected}" "${targets}" >&2
        exit 1
    fi
}

printf '# First page\n' > "${FIXTURE}/benchmarks/Alpha/first.jmd"
printf '# Second page\n' > "${FIXTURE}/benchmarks/Alpha/second.jmd"
printf '# Normal page\n' > "${FIXTURE}/benchmarks/Beta/page.jmd"
printf '\n' >> "${FIXTURE}/benchmarks/Beta/Project.toml"
printf 'split = true\nrunner = ["split-test"]\ntimeout = 90\n' \
    > "${FIXTURE}/benchmarks/Alpha/benchmark_config.toml"
SPLIT_TARGETS=$'benchmarks/Alpha/first.jmd\nbenchmarks/Alpha/second.jmd'
assert_detection "${SPLIT_TARGETS}"$'\nbenchmarks/Beta'
if [[ $(grep -o '"runner":\["split-test"\],"timeout":90,"julia_version":"1.11"' \
    "${TEST_ROOT}/case-github-output" | wc -l) -ne 2 ]]; then
    echo 'Split targets must inherit their folder runner, timeout, and Julia version.' >&2
    exit 1
fi

# A shared environment change coalesces file changes into one job per page.
printf '\n' >> "${FIXTURE}/benchmarks/Alpha/Manifest.toml"
printf '\n' >> "${FIXTURE}/benchmarks/Alpha/first.jmd"
printf '\n' >> "${FIXTURE}/benchmarks/Beta/Project.toml"
printf '\n' >> "${FIXTURE}/benchmarks/Beta/page.jmd"
assert_detection "${SPLIT_TARGETS}"$'\nbenchmarks/Beta'

printf '\n' >> "${FIXTURE}/benchmarks/Alpha/first.jmd"
assert_detection 'benchmarks/Alpha/first.jmd'

assert_detection "${SPLIT_TARGETS}" benchmarks/Alpha
assert_detection 'benchmarks/Alpha/first.jmd' benchmarks/Alpha/first.jmd
assert_detection 'benchmarks/Beta' benchmarks/Beta

printf 'split = false\n' > "${FIXTURE}/benchmarks/Alpha/benchmark_config.toml"
assert_detection 'benchmarks/Alpha'
assert_detection 'benchmarks/Alpha' benchmarks/Alpha

printf 'split = true   # One job per page\n' > "${FIXTURE}/benchmarks/Alpha/benchmark_config.toml"
printf '# Hidden page\n' > "${FIXTURE}/benchmarks/Alpha/.hidden.jmd"
assert_detection $'benchmarks/Alpha/.hidden.jmd\n'"${SPLIT_TARGETS}"

rm "${FIXTURE}/benchmarks/Alpha/"*.jmd "${FIXTURE}/benchmarks/Alpha/.hidden.jmd"
assert_detection '' benchmarks/Alpha
if ! grep -Fxq 'has_changes=false' "${TEST_ROOT}/case-github-output"; then
    echo 'An empty split folder must not schedule a job.' >&2
    exit 1
fi

if (cd "${FIXTURE}" && .github/scripts/detect-changed-benchmarks.sh benchmarks/Missing); then
    echo 'An invalid manual target must fail.' >&2
    exit 1
fi

echo 'Split and normal folders, file coalescing, config inheritance, and manual targets passed.'
