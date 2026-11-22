#!/bin/bash
# Shared helper: read runner configuration for a benchmark target.
# Usage: source this file, then call get_benchmark_config <target>
# Returns: "runner_json|timeout|julia_version|threads"
#   e.g. '["self-hosted", "benchmark"]|12000|1.10|32'

# Read a top-level key from a TOML file with Julia's TOML parser, which keeps
# table and array-of-table entries table-local whatever their indentation.
# Values print in a canonical spelling — strings keep their TOML quoting,
# integers come out in decimal (1_6, 0x10 and friends normalize), booleans as
# true/false and arrays as a JSON-style ["a", "b"] list. A missing file or key
# prints nothing; a malformed file fails with the parser's error.
read_toml_key() {
    local file="$1" key="$2"
    if [[ ! -f "${file}" ]]; then
        return
    fi
    julia --startup-file=no -e '
        using TOML
        value = get(TOML.parsefile(ARGS[1]), ARGS[2], nothing)
        if value isa AbstractVector
            print("[", join((sprint(show, element) for element in value), ", "), "]")
        elseif value isa AbstractString
            print(sprint(show, value))
        elseif !isnothing(value)
            print(value)
        end' -- "${file}" "${key}"
}

# Normalize a TOML integer literal to a decimal string. Accepts the literals
# Julia's TOML parser accepts: decimal, 0x hex, 0o octal and 0b binary forms,
# underscores between digits, and a leading +. Returns 1 otherwise.
parse_toml_int() {
    local raw="${1#+}" base digits pattern
    case "${raw}" in
        0x*) base=16; digits="${raw:2}" ;;
        0o*) base=8; digits="${raw:2}" ;;
        0b*) base=2; digits="${raw:2}" ;;
        *) base=10; digits="${raw}" ;;
    esac
    case "${base}" in
        16) pattern='^[0-9a-fA-F]+(_[0-9a-fA-F]+)*$' ;;
        8) pattern='^[0-7]+(_[0-7]+)*$' ;;
        2) pattern='^[01]+(_[01]+)*$' ;;
        *) pattern='^(0|[1-9][0-9]*(_[0-9]+)*)$' ;;
    esac
    if [[ "${digits}" =~ ${pattern} ]]; then
        echo "$(( ${base}#${digits//_/} ))"
    else
        return 1
    fi
}

# Minimum Julia version — clamp to LTS so old Manifests don't pull ancient versions
JULIA_LTS="1.10"

# Compare two major.minor version strings. Returns 0 if $1 < $2.
version_lt() {
    local major1 minor1 major2 minor2
    major1=$(echo "$1" | cut -d. -f1)
    minor1=$(echo "$1" | cut -d. -f2)
    major2=$(echo "$2" | cut -d. -f1)
    minor2=$(echo "$2" | cut -d. -f2)
    if [[ "${major1}" -lt "${major2}" ]]; then
        return 0
    elif [[ "${major1}" -eq "${major2}" && "${minor1}" -lt "${minor2}" ]]; then
        return 0
    fi
    return 1
}

# Extract Julia minor version (e.g. "1.10") from a Manifest.toml's julia_version field.
# Clamps to LTS as a minimum — old Manifests (1.8, 1.9) will use LTS instead.
get_julia_version() {
    local bench_dir="$1"
    local manifest="${bench_dir}/Manifest.toml"
    local raw
    raw=$(read_toml_key "${manifest}" "julia_version") || return 1
    if [[ -n "${raw}" ]]; then
        local version
        version=$(echo "${raw}" | tr -d '"' | cut -d. -f1,2)
        # Clamp to LTS minimum
        if version_lt "${version}" "${JULIA_LTS}"; then
            echo "${JULIA_LTS}"
        else
            echo "${version}"
        fi
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

    local runner timeout julia_version threads
    local config_file="${bench_dir}/benchmark_config.toml"
    local defaults_file=".github/benchmark_defaults.toml"

    runner=$(read_toml_key "${config_file}" "runner") || return 1
    if [[ -z "${runner}" ]]; then
        runner=$(read_toml_key "${defaults_file}" "runner") || return 1
    fi
    if [[ -z "${runner}" ]]; then
        runner='["self-hosted", "benchmark"]'
    fi

    timeout=$(read_toml_key "${config_file}" "timeout") || return 1
    if [[ -z "${timeout}" ]]; then
        timeout=$(read_toml_key "${defaults_file}" "timeout") || return 1
    fi
    if [[ -z "${timeout}" ]]; then
        timeout="12000"
    fi

    threads=$(read_toml_key "${config_file}" "threads") || return 1
    if [[ -z "${threads}" ]]; then
        threads=$(read_toml_key "${defaults_file}" "threads") || return 1
    fi
    if [[ -z "${threads}" ]]; then
        threads="32"
    fi

    local parsed_threads
    if parsed_threads=$(parse_toml_int "${threads}") \
        && [[ "${parsed_threads}" -ge 1 ]]; then
        threads=${parsed_threads}
    else
        echo "Expected a positive integer thread count for ${target}, got '${threads}'" >&2
        return 1
    fi

    julia_version=$(get_julia_version "${bench_dir}") || return 1
    if [[ -z "${julia_version}" ]]; then
        julia_version="1.10"
    fi

    echo "${runner}|${timeout}|${julia_version}|${threads}"
}

# Backwards-compatible alias
get_runner_config() {
    local config
    config=$(get_benchmark_config "$1") || return 1
    # Return just runner|timeout for existing callers
    echo "${config}" | cut -d'|' -f1-2
}
