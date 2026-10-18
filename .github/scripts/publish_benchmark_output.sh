#!/bin/bash
set -euo pipefail

# Set up SSH for pushing to SciMLBenchmarksOutput
mkdir -p ~/.ssh
ssh-keyscan github.com >> ~/.ssh/known_hosts
eval "$(ssh-agent -s)"

# Write deploy key from secret
DEPLOY_KEY_FILE=$(mktemp)
echo "${DEPLOY_KEY}" > "${DEPLOY_KEY_FILE}"
chmod 600 "${DEPLOY_KEY_FILE}"
ssh-add "${DEPLOY_KEY_FILE}"

git config --global user.email "github-actions@sciml.ai"
git config --global user.name "SciML Benchmarks CI"

# Nothing to do if no benchmark job uploaded artifacts (download-artifact is
# continue-on-error, so this script still runs).
have_output=false
for d in docs html notebook pdf script markdown; do
    if [[ -d "${d}" ]]; then
        have_output=true
    fi
done
if [[ "${have_output}" != true ]]; then
    echo "No benchmark output directories found; nothing to publish."
    rm -f "${DEPLOY_KEY_FILE}"
    exit 0
fi

# Publish jobs from different master runs are not serialized by a concurrency
# group (that would let a newly queued publish cancel a pending one and lose
# its results). Instead, tolerate concurrent publishes: if the push is
# rejected because another run pushed first, re-clone and try again. Each
# attempt copies the same files on top of the newest SciMLBenchmarksOutput, so
# retrying is safe; there is nothing to merge.
MAX_ATTEMPTS=${PUBLISH_MAX_ATTEMPTS:-6}
OUTPUT_REPO=${PUBLISH_OUTPUT_REPO:-git@github.com:SciML/SciMLBenchmarksOutput}
# Optional hook run after the clone, before the push (used by the local test to
# inject a concurrent push). Unset in CI.
PUBLISH_PRE_PUSH_HOOK=${PUBLISH_PRE_PUSH_HOOK:-}
attempt=1
while true; do
    temp_dir=$(mktemp -d)
    git -C "${temp_dir}" clone --depth 1 "${OUTPUT_REPO}" .

    # Copy output artifacts into it
    for d in docs html notebook pdf script markdown; do
        if [[ -d "${d}" ]]; then
            cp -vRa "${d}/" "${temp_dir}"
        fi
    done

    git -C "${temp_dir}" add .
    if git -C "${temp_dir}" diff --cached --quiet; then
        echo "No changes to publish."
        rm -rf "${temp_dir}"
        break
    fi

    git -C "${temp_dir}" commit -q -m "Automatic build
Published by build of: ${GITHUB_REPOSITORY}@${GITHUB_SHA}"

    if [[ -n "${PUBLISH_PRE_PUSH_HOOK}" ]]; then
        "${PUBLISH_PRE_PUSH_HOOK}" "${attempt}"
    fi

    if git -C "${temp_dir}" push; then
        echo "Published on attempt ${attempt}."
        rm -rf "${temp_dir}"
        break
    fi

    rm -rf "${temp_dir}"
    if [[ "${attempt}" -ge "${MAX_ATTEMPTS}" ]]; then
        echo "::error::Failed to push to SciMLBenchmarksOutput after ${attempt} attempts."
        rm -f "${DEPLOY_KEY_FILE}"
        exit 1
    fi
    # Jitter so two publishers that collided don't retry in lock-step.
    sleep $((10 + RANDOM % 20))
    attempt=$((attempt + 1))
    echo "Push rejected (another publish likely landed first); retrying (attempt ${attempt}/${MAX_ATTEMPTS})..."
done

# Clean up
rm -f "${DEPLOY_KEY_FILE}"
