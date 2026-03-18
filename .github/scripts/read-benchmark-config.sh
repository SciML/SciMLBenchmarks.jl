#!/bin/bash
# Shared helper: read runner configuration for a benchmark target.
# Usage: source this file, then call get_benchmark_config <target>
# Returns: "runner_json|timeout|julia_version"
#   e.g. '["self-hosted", "benchmark"]|12000|1.10'

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

# Extract Julia minor version (e.g. "1.10") from a Manifest.toml's julia_version field
get_julia_version() {
    local bench_dir="$1"
    local manifest="${bench_dir}/Manifest.toml"
    local raw
    raw=$(read_toml_key "${manifest}" "julia_version")
    if [[ -n "${raw}" ]]; then
        # Strip quotes and extract major.minor (e.g. "1.10.10" -> "1.10")
        local version
        version=$(echo "${raw}" | tr -d '"' | cut -d. -f1,2)
        echo "${version}"
    fi
}

# Get full config for a benchmark target (file or directory)
get_benchmark_config() {
    local target="$1"
    local bench_dir

    if [[ -d "${target}" ]]; then
        bench_dir="${target}"
    else
        bench_dir="$(dirname "${target}")"
    fi

    local runner timeout julia_version
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

    julia_version=$(get_julia_version "${bench_dir}")
    if [[ -z "${julia_version}" ]]; then
        julia_version="1.10"
    fi

    echo "${runner}|${timeout}|${julia_version}"
}

# Backwards-compatible alias
get_runner_config() {
    local config
    config=$(get_benchmark_config "$1")
    # Return just runner|timeout for existing callers
    echo "${config%|*}"
}
