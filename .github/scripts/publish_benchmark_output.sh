#!/bin/bash
set -euo pipefail

have_output=false
for d in docs html notebook pdf script markdown; do
    if [[ -d "${d}" ]]; then
        have_output=true
    fi
done
if [[ "${have_output}" != true ]]; then
    echo "No benchmark output directories found; nothing to publish."
    exit 0
fi

publish_dir=$(mktemp -d)
trap 'rm -rf "$publish_dir"' EXIT
(umask 077; printf '%s\n' "$DEPLOY_KEY" > "$publish_dir/deploy_key")
unset DEPLOY_KEY
script_dir=$(cd "$(dirname "$0")" && pwd)
cat > "$publish_dir/ssh_config" <<EOF
BatchMode yes
IdentitiesOnly yes
IdentityAgent none
IdentityFile "$publish_dir/deploy_key"
StrictHostKeyChecking yes
GlobalKnownHostsFile /dev/null
UserKnownHostsFile "$script_dir/github_known_hosts"
EOF

# Use Git Bash's OpenSSH on Windows so Git and SSH resolve the same paths.
printf -v GIT_SSH_COMMAND '/usr/bin/ssh -F %q' "$publish_dir/ssh_config"
export GIT_SSH_COMMAND

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
    temp_dir="$publish_dir/output"
    mkdir "$temp_dir"
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
        break
    fi

    git -C "${temp_dir}" -c user.email=github-actions@sciml.ai -c user.name="SciML Benchmarks CI" commit -q -m "Automatic build
Published by build of: ${GITHUB_REPOSITORY}@${GITHUB_SHA}"

    if [[ -n "${PUBLISH_PRE_PUSH_HOOK}" ]]; then
        "${PUBLISH_PRE_PUSH_HOOK}" "${attempt}"
    fi

    if git -C "${temp_dir}" push; then
        echo "Published on attempt ${attempt}."
        break
    fi

    rm -rf "${temp_dir}"
    if [[ "${attempt}" -ge "${MAX_ATTEMPTS}" ]]; then
        echo "::error::Failed to push to SciMLBenchmarksOutput after ${attempt} attempts."
        exit 1
    fi
    # Jitter so two publishers that collided don't retry in lock-step.
    sleep $((10 + RANDOM % 20))
    attempt=$((attempt + 1))
    echo "Push rejected (another publish likely landed first); retrying (attempt ${attempt}/${MAX_ATTEMPTS})..."
done
