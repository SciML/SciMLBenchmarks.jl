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

# Clone SciMLBenchmarksOutput to temporary directory
temp_dir=$(mktemp -d)
git -C "${temp_dir}" clone git@github.com:SciML/SciMLBenchmarksOutput .

# Copy output artifacts into it
for d in docs html notebook pdf script markdown; do
    if [[ -d "${d}" ]]; then
        cp -vRa "${d}/" "${temp_dir}"
    fi
done

# Commit and push
set -e
git -C "${temp_dir}" add .
if git -C "${temp_dir}" diff --cached --quiet; then
    echo "No changes to publish."
else
    git -C "${temp_dir}" commit -m "Automatic build
Published by build of: ${GITHUB_REPOSITORY}@${GITHUB_SHA}"
    git -C "${temp_dir}" push
fi

# Clean up
rm -f "${DEPLOY_KEY_FILE}"
rm -rf "${temp_dir}"
