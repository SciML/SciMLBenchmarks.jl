#!/bin/bash
# Shared helper: read runner configuration for a benchmark target.
# Usage: source this file, then call get_runner_config <target>
# Returns: "runner_json|timeout" e.g. '["self-hosted", "benchmark"]|12000'

# Read a key from a TOML file (simple parser for flat keys and arrays)
read_toml_key() {
    local file="$1" key="$2"
    if [[ ! -f "${file}" ]]; then
        return
    fi
    local line
    line=$(grep -E "^${key}\s*=" "${file}" 2>/dev/null | head -1) || true
    if [[ -n "${line}" ]]; then
        echo "${line#*=}" | sed 's/^[[:space:]]*//'
    fi
}

# Get runner config for a benchmark target (file or directory)
get_runner_config() {
    local target="$1"
    local bench_dir

    if [[ -d "${target}" ]]; then
        bench_dir="${target}"
    else
        bench_dir="$(dirname "${target}")"
    fi

    local runner timeout
    local config_file="${bench_dir}/benchmark_config.toml"
    local defaults_file=".github/benchmark_defaults.toml"

    runner=$(read_toml_key "${config_file}" "runner")
    if [[ -z "${runner}" ]]; then
        runner=$(read_toml_key "${defaults_file}" "runner")
    fi
    if [[ -z "${runner}" ]]; then
        runner='["self-hosted", "benchmark"]'
    fi

    timeout=$(read_toml_key "${config_file}" "timeout")
    if [[ -z "${timeout}" ]]; then
        timeout=$(read_toml_key "${defaults_file}" "timeout")
    fi
    if [[ -z "${timeout}" ]]; then
        timeout="12000"
    fi

    echo "${runner}|${timeout}"
}
