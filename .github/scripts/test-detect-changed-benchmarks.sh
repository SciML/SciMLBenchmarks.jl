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

echo "Push detection selected only changes since the event's before SHA."
