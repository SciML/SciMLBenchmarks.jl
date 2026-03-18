#!/bin/bash
set -eo pipefail

# Detect changed benchmark files and produce a JSON matrix for GHA.
# Replicates the logic from .buildkite/path_processors/project-coalescing:
#   - A changed .jmd file triggers a rebuild of just that file
#   - A changed .toml file triggers a rebuild of its entire benchmark directory
#   - If a directory is already being rebuilt, individual .jmd files in it are suppressed

if [[ "${GITHUB_EVENT_NAME}" == "pull_request" ]]; then
    BASE_SHA="${GITHUB_BASE_REF}"
    git fetch origin "${BASE_SHA}" --depth=1 2>/dev/null || true
    CHANGED_FILES=$(git diff --name-only "origin/${BASE_SHA}...HEAD" -- 'benchmarks/')
else
    # Push event: compare with parent commit
    CHANGED_FILES=$(git diff --name-only HEAD~1 -- 'benchmarks/' 2>/dev/null || true)
fi

# Also check if src/ or root toml changed — if so, we should note it but
# the forerunner in buildkite only watched benchmarks/** for launching benchmark jobs.
# src/ changes triggered test_sciml.yml instead. We replicate that behavior.

declare -A FILES     # .jmd file -> its project directory
declare -A PROJECTS  # project directories that need full rebuild

# Find the Project.toml-containing directory for a file
find_project() {
    local d="$1"
    while [[ "${d}" =~ benchmarks/.* ]]; do
        if [[ -f "${d}/Project.toml" ]]; then
            echo "${d}"
            return
        fi
        d="$(dirname "${d}")"
    done
}

while IFS= read -r f; do
    [[ -z "${f}" ]] && continue
    proj=$(find_project "$(dirname "${f}")")
    if [[ -z "${proj}" ]]; then
        echo "::warning::Unable to find project for ${f}"
        continue
    fi

    if [[ "${f}" == *.jmd ]]; then
        FILES["${f}"]="${proj}"
    elif [[ "${f}" == *.toml ]]; then
        PROJECTS["${proj}"]=1
    fi
done <<< "${CHANGED_FILES}"

# Build target list: project directories first, then individual .jmd files
# whose projects are not already being fully rebuilt
BUILD_TARGETS=()

if [[ ${#PROJECTS[@]} -gt 0 ]]; then
    for proj in "${!PROJECTS[@]}"; do
        BUILD_TARGETS+=("${proj}")
    done
fi

if [[ ${#FILES[@]} -gt 0 ]]; then
    for f in "${!FILES[@]}"; do
        proj="${FILES[$f]}"
        if [[ -z "${PROJECTS[$proj]+x}" ]]; then
            BUILD_TARGETS+=("${f}")
        fi
    done
fi

# Output as JSON array for GHA matrix
if [[ ${#BUILD_TARGETS[@]} -eq 0 ]]; then
    echo "has_changes=false" >> "$GITHUB_OUTPUT"
    echo "matrix=[]" >> "$GITHUB_OUTPUT"
    echo "No benchmark changes detected."
else
    echo "has_changes=true" >> "$GITHUB_OUTPUT"
    # Build JSON array
    JSON="["
    for i in "${!BUILD_TARGETS[@]}"; do
        if [[ $i -gt 0 ]]; then
            JSON+=","
        fi
        JSON+="\"${BUILD_TARGETS[$i]}\""
    done
    JSON+="]"
    echo "matrix=${JSON}" >> "$GITHUB_OUTPUT"
    echo "Detected benchmark targets: ${BUILD_TARGETS[*]}"
fi
