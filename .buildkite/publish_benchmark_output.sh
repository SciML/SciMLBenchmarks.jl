#!/bin/bash

# Clone SciMLBenchmarksOutput to temporary directory
temp_dir=$(mktemp -d)
git -C "${temp_dir}" clone git@github.com:SciML/SciMLBenchmarksOutput .

# Copy our output artifacts into it:
for d in html markdown notebook pdf script; do
    cp -vRa "${d}/" "${temp_dir}"
done

# Commit the result up to output
git -C "${temp_dir}" add .
git -C "${temp_dir}" commit -m "Automatic build\nPublished by build of: ${BUILDKITE_REPO%.git}/commit/${BUILDKITE_COMMIT}"
git -C "${temp_dir}" push

rm -rf "${temp_dir}"
