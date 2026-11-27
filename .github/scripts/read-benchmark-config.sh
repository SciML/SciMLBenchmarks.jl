#!/bin/bash
# Shared helper: read runner configuration for a benchmark target.
# Usage: source this file, then call get_benchmark_config <target>
# Returns: "runner_json|timeout|julia_version|threads"
#   e.g. '["self-hosted", "benchmark"]|12000|1.10|32'

# One Julia process resolves every needed top-level key of a folder's
# benchmark_config.toml (falling back to .github/benchmark_defaults.toml
# and then to built-in values) plus the Manifest's julia_version. Results
# are kept in a scratch directory keyed by the input files' contents, so
# repeat lookups for a folder's pages, matrix entry and summary log never
# relaunch Julia — including across command-substitution subshells, which
# inherit this variable after the file is sourced.
BENCHMARK_CONFIG_CACHE_DIR=$(mktemp -d "${TMPDIR:-/tmp}/sciml-benchmark-config.XXXXXX")
trap 'rm -rf "${BENCHMARK_CONFIG_CACHE_DIR}"' EXIT

# Emit a folder's resolved configuration, one record per line: timeout,
# threads, Julia version (or empty), the split flag, then every runner
# label raw. A malformed file or out-of-schema value fails loudly.
resolve_benchmark_config() {
    local bench_dir="$1" target="$2"
    julia --startup-file=no -e '
        using TOML
        config_file, defaults_file, manifest_file, target = ARGS
        load(path) = isfile(path) ? TOML.parsefile(path) : Dict{String, Any}()
        config = load(config_file)
        defaults = load(defaults_file)
        pick(key, fallback) = get(config, key, get(defaults, key, fallback))

        runner = pick("runner", ["self-hosted", "benchmark"])
        runner isa AbstractVector || (runner = Any[runner])
        all(label -> label isa AbstractString, runner) ||
            error("runner for $(target) must be a string or an array of strings")
        any(label -> occursin("\n", label) || occursin("\r", label), runner) &&
            error("runner labels for $(target) must not contain line breaks")
        (isempty(runner) || any(isempty, runner)) &&
            error("runner for $(target) must not be empty or contain empty labels")

        timeout = pick("timeout", 12000)
        timeout isa Integer && !(timeout isa Bool) && timeout >= 1 ||
            error("Expected a positive integer timeout for $(target), got $(repr(timeout))")

        threads = pick("threads", 32)
        threads isa Integer && !(threads isa Bool) && threads >= 1 ||
            error("Expected a positive integer thread count for $(target), got $(repr(threads))")

        split = pick("split", false)
        split isa Bool || error("split for $(target) must be a boolean")

        julia_version = ""
        if isfile(manifest_file)
            raw = get(load(manifest_file), "julia_version", nothing)
            if raw isa AbstractString
                version = tryparse(VersionNumber, raw)
                if !isnothing(version)
                    minor = VersionNumber(version.major, version.minor)
                    julia_version = minor < v"1.10" ? "1.10" : "$(version.major).$(version.minor)"
                end
            end
        end

        println(timeout)
        println(threads)
        println(julia_version)
        println(split)
        foreach(println, runner)
    ' -- "${bench_dir}/benchmark_config.toml" ".github/benchmark_defaults.toml" \
        "${bench_dir}/Manifest.toml" "${target}"
}

# Print the resolved records for a target's folder, launching Julia at most
# once per distinct configuration in this process tree.
benchmark_config_records() {
    local target="$1" bench_dir
    if [[ -d "${target}" ]]; then
        bench_dir="${target}"
    else
        bench_dir="$(dirname "${target}")"
    fi
    local file key
    key=$(
        for input in "${bench_dir}/benchmark_config.toml" \
            ".github/benchmark_defaults.toml" "${bench_dir}/Manifest.toml"; do
            [[ -f "${input}" ]] && md5sum "${input}"
        done | md5sum | cut -d' ' -f1
    )
    file="${BENCHMARK_CONFIG_CACHE_DIR}/${key}"
    if [[ ! -f "${file}" ]]; then
        local raw
        raw=$(resolve_benchmark_config "${bench_dir}" "${target}") || return 1
        printf '%s\n' "${raw}" > "${file}"
    fi
    cat "${file}"
}

# Split resolved records into named fields: CFG_TIMEOUT, CFG_THREADS,
# CFG_JULIA_VERSION, CFG_SPLIT and the CFG_LABELS array.
parse_benchmark_config() {
    local -a fields
    mapfile -t fields <<<"$1"
    ((${#fields[@]} >= 4)) || return 1
    CFG_TIMEOUT=${fields[0]}
    CFG_THREADS=${fields[1]}
    CFG_JULIA_VERSION=${fields[2]:-1.10}
    CFG_SPLIT=${fields[3]}
    CFG_LABELS=("${fields[@]:4}")
}

# Serialize raw runner labels to a JSON array in the same ["a", "b"]
# spelling the workflow matrix embeds. jq does the escaping, so labels may
# contain $, quotes, backslashes and other JSON-special characters. Pipes
# are escaped so that callers can split the |-delimited config string.
runner_labels_json() {
    if (($# == 0)); then
        printf '[]\n'
    else
        printf '%s\n' "$@" | jq -Rrn \
            '[inputs] | "[" + (map(@json | split("|") | join("\\u007c")) | join(", ")) + "]"'
    fi
}

# Get full config for a benchmark target (file or directory)
get_benchmark_config() {
    local target="$1" records runner_json
    records=$(benchmark_config_records "${target}") || return 1
    parse_benchmark_config "${records}" || return 1
    runner_json=$(runner_labels_json ${CFG_LABELS[@]+"${CFG_LABELS[@]}"}) || return 1
    echo "${runner_json}|${CFG_TIMEOUT}|${CFG_JULIA_VERSION}|${CFG_THREADS}"
}

# Read a top-level key from a TOML file — the test suites' independent
# oracle for the resolved values. A missing file or key prints nothing; a
# malformed file fails with the parser's error.
read_toml_key() {
    local file="$1" key="$2"
    if [[ ! -f "${file}" ]]; then
        return
    fi
    julia --startup-file=no -e '
        using TOML
        value = get(TOML.parsefile(ARGS[1]), ARGS[2], nothing)
        if value isa AbstractVector
            foreach(println, value)
        elseif !isnothing(value)
            println(value)
        end' -- "${file}" "${key}"
}

# Backwards-compatible alias
get_runner_config() {
    local config
    config=$(get_benchmark_config "$1") || return 1
    # Return just runner|timeout for existing callers
    echo "${config}" | cut -d'|' -f1-2
}
