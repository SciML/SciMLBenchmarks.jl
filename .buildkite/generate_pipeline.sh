#!/usr/bin/env bash

## This script will trigger multiple jobs depending on what has been changed in the PR currently under test:
##  * If model defintion source files have been changed, the model will be re-trained.
##  * If SciMLBenchmarks code itself has been changed, the tests themselves will be run.

if [[ "${BUILDKITE_PULL_REQUEST}" == "false" ]]; then
    # If this is not a pull request, just return
    # TODO: We need to actually register models upon commit.
    echo "Early-exiting as this is not a pull request"
    exit 0
fi

# Compare current tip-of-branch with origin
git fetch origin "${BUILDKITE_PULL_REQUEST_BASE_BRANCH}"
TARGET_BRANCH="remotes/origin/${BUILDKITE_PULL_REQUEST_BASE_BRANCH}"
COMPARE_AGAINST=$(git merge-base --fork-point ${TARGET_BRANCH} HEAD)
git fetch origin "refs/pull/${BUILDKITE_PULL_REQUEST}/head:refs/remotes/origin/pr/${BUILDKITE_PULL_REQUEST}"

# Find the directories of everything that has been modified; if it
# contains something frmo `benchmarks/` we're probably gonna retrain
# some models, other stuff usually means we're going to run the tests
# of this package itself.
CURR_HEAD=$(git rev-parse HEAD)
MODIFIED_FILES=( $(git diff-tree --no-commit-id --name-only -r HEAD "${COMPARE_AGAINST}") )

echo "Found ${#MODIFIED_FILES[@]} modified files between ${CURR_HEAD} and ${COMPARE_AGAINST}"

## Search for modified model directories
# Coalesce these directories by finding the overall model names of all those modified directories
declare -A MODELS
for FILE in ${MODIFIED_FILES[@]}; do
    DIR="$(dirname "${FILE}")"
    # If we modify, e.g. `benchmarks/AdaptiveSDE/AdaptiveEfficiencyTests.jmd`, keep the `folder/jmdfile`
    # name until we hit `benchmark` or find a `Project.toml` file.
    while [[ "${DIR}" == benchmarks/* ]]; do
        if [[ -f "${DIR}/Project.toml" ]]; then
            # Strip `benchmarks` from the beginning to get the model name
            MODEL="${FILE#benchmarks/}"
            MODELS["${MODEL}"]=1

            echo "${FILE} triggered model '${MODEL}' rebuild"
            break
        fi
        DIR=$(dirname ${DIR})
    done
done

## Search for files that matter to the package tests
for FILE in ${MODIFIED_FILES[@]}; do
    # Anything in the `src/` or `test/` directories (or any of the `.toml` files
    # in the top-level package directory) should trigger the tests
    if [[ "${FILE}" == src/* ]] ||
       [[ "${FILE}" == test/* ]] ||
       [[ "${FILE}" =~ "[^/]+\.toml" ]]; then
        RUN_SCIML_TESTS="true"
        echo "File ${FILE} triggered package tests"
        break
    fi
done

# If we should run our own tests, send the configuration to buildkite
if [[ "${RUN_SCIML_TESTS}" == "true" ]]; then
    cat ./.buildkite/test_sciml.yml | buildkite-agent pipeline upload
fi

# Now that we've got a unique list of models, we'll template that into
# `rebuild_model.yml` and send it to the buildkite agent, one pipeline
# run per model that we've identified needs to be built.
for MODEL in "${!MODELS[@]}"; do
    echo "Rebuilding '${MODEL}'"
    sed -e "s&{MODEL}&${MODEL}&g" < ./.buildkite/rebuild_model.yml | buildkite-agent pipeline upload
done