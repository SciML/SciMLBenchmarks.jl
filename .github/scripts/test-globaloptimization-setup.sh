#!/bin/bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SETUP_SCRIPT="${SETUP_SCRIPT:-${REPO_ROOT}/benchmarks/GlobalOptimization/setup.sh}"
TEST_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/globalopt-setup.XXXXXX")
trap 'rm -rf "${TEST_ROOT}"' EXIT

FAILURES=0
ok() { printf 'ok - %s\n' "$1"; }
bad() {
    printf 'not ok - %s\n' "$1"
    FAILURES=$((FAILURES + 1))
}
cond() {
    local desc="$1"
    shift
    if "$@"; then ok "${desc}"; else bad "${desc}"; fi
}
cond_not() {
    local desc="$1"
    shift
    if "$@"; then bad "${desc}"; else ok "${desc}"; fi
}

new_case() {
    CASE_DIR="${TEST_ROOT}/$1"
    CASE_HOME="${CASE_DIR}/${2:-home}"
    mkdir -p "${CASE_DIR}/bin" "${CASE_HOME}"
    : > "${CASE_DIR}/env"
}

stub_uname_mingw() {
    cat > "${CASE_DIR}/bin/uname" <<'EOF'
#!/bin/bash
case "${1:--s}" in
    -m) echo x86_64 ;;
    *) echo MINGW64_NT-10.0-26100 ;;
esac
EOF
    chmod +x "${CASE_DIR}/bin/uname"
}

stub_cygpath_ok() {
    cat > "${CASE_DIR}/bin/cygpath" <<'EOF'
#!/bin/bash
if [[ "${1:-}" != "-w" || $# -ne 2 ]]; then
    exit 9
fi
printf 'Z:%s\n' "${2//\//\\}"
EOF
    chmod +x "${CASE_DIR}/bin/cygpath"
}

stub_cygpath_fails() {
    printf '#!/bin/bash\nexit 7\n' > "${CASE_DIR}/bin/cygpath"
    chmod +x "${CASE_DIR}/bin/cygpath"
}

stub_network() {
    printf '#!/bin/bash\nprintf "curl %%s\\n" "$*" >> "${STUB_LOG}"\n' \
        > "${CASE_DIR}/bin/curl"
    printf '#!/bin/bash\nprintf "tar %%s\\n" "$*" >> "${STUB_LOG}"\n' \
        > "${CASE_DIR}/bin/tar"
    chmod +x "${CASE_DIR}/bin/curl" "${CASE_DIR}/bin/tar"
}

install_micromamba_stub() {
    local mamba_dir="${CASE_HOME}/.cache/sciml-benchmarks"
    mkdir -p "${mamba_dir}"
    cat > "${mamba_dir}/$1" <<'EOF'
#!/bin/bash
{
    printf 'MAMBA_ROOT_PREFIX=%s\n' "${MAMBA_ROOT_PREFIX-<unset>}"
    printf 'arg=%s\n' "$@"
} >> "${STUB_LOG}"
if [[ "${1:-}" == create && -n "${FAKE_PYTHON_EXE:-}" ]]; then
    mkdir -p "$(dirname "${FAKE_PYTHON_EXE}")"
    printf '#!/bin/bash\nexit 0\n' > "${FAKE_PYTHON_EXE}"
    chmod +x "${FAKE_PYTHON_EXE}"
fi
EOF
    chmod +x "${mamba_dir}/$1"
}

install_python_stub() {
    local py="${CASE_HOME}/.cache/sciml-benchmarks/$1"
    mkdir -p "$(dirname "${py}")"
    printf '#!/bin/bash\nexit 0\n' > "${py}"
    chmod +x "${py}"
}

run_setup() {
    rc=0
    (
        export HOME="${CASE_HOME}"
        export STUB_LOG="${CASE_DIR}/stub.log"
        export FAKE_PYTHON_EXE="${CASE_HOME}/.cache/sciml-benchmarks/$1"
        export BENCHMARK_ENV_FILE="${CASE_DIR}/env"
        PATH="${CASE_DIR}/bin:${PATH}" bash "${SETUP_SCRIPT}"
    ) > "${CASE_DIR}/out.log" 2>&1 || rc=$?
}

stub_log_grep() {
    [[ -f "${CASE_DIR}/stub.log" ]] && grep "$@" "${CASE_DIR}/stub.log"
}

check_windows_conversion_failure() {
    local name="$1" cygpath_stub="$2" variant
    for variant in no-cached-python cached-python; do
        new_case "${name}-${variant}"
        stub_uname_mingw
        stub_network
        if [[ "${cygpath_stub}" != none ]]; then
            "${cygpath_stub}"
        fi
        install_micromamba_stub micromamba.exe
        if [[ "${variant}" == cached-python ]]; then
            install_python_stub globalopt-python/python.exe
        fi
        run_setup globalopt-python/python.exe
        cond "${name} ${variant}: setup exits nonzero" \
            test "${rc}" -ne 0
        cond_not "${name} ${variant}: micromamba invoked" \
            stub_log_grep -q '^arg='
        cond_not "${name} ${variant}: micromamba got empty root or arg" \
            stub_log_grep -Eqx 'MAMBA_ROOT_PREFIX=|arg='
        cond_not "${name} ${variant}: empty JULIA_PYTHONCALL_EXE appended" \
            grep -Fq 'JULIA_PYTHONCALL_EXE=""' "${CASE_DIR}/env"
    done
}

# Case 1: a missing cygpath on the Windows path must fail loudly.
check_windows_conversion_failure case1-missing-cygpath none

# Case 2: a failing cygpath on the Windows path must fail loudly.
check_windows_conversion_failure case2-failing-cygpath stub_cygpath_fails

# Case 3: successful conversion with spaces in HOME.
new_case case3 "fake home"
stub_uname_mingw
stub_network
stub_cygpath_ok
install_micromamba_stub micromamba.exe
run_setup globalopt-python/python.exe

cond "case3: setup exits zero" test "${rc}" -eq 0

PREFIX_POSIX="${CASE_HOME}/.cache/sciml-benchmarks/globalopt-python"
EXPECTED_PREFIX="$("${CASE_DIR}/bin/cygpath" -w "${PREFIX_POSIX}")"
EXPECTED_ROOT="$("${CASE_DIR}/bin/cygpath" -w \
    "${CASE_HOME}/.cache/sciml-benchmarks/mamba-root")"
EXPECTED_PY="${EXPECTED_PREFIX}\\python.exe"

{
    printf 'MAMBA_ROOT_PREFIX=%s\n' "${EXPECTED_ROOT}"
    printf 'arg=%s\n' create -y -p "${EXPECTED_PREFIX}" -c conda-forge \
        --override-channels python=3.12 numpy scipy matplotlib cma
} > "${CASE_DIR}/expected.log"
cond "case3: single native -p prefix argument preserved" \
    cmp "${CASE_DIR}/expected.log" "${CASE_DIR}/stub.log"
cond "case3: env JULIA_LOAD_PATH" \
    grep -Fxq 'export JULIA_LOAD_PATH="@:@stdlib"' "${CASE_DIR}/env"
cond "case3: env JULIA_CONDAPKG_BACKEND" \
    grep -Fxq 'export JULIA_CONDAPKG_BACKEND="Null"' "${CASE_DIR}/env"
cond "case3: root-level python.exe exported" \
    grep -Fxq "export JULIA_PYTHONCALL_EXE=\"${EXPECTED_PY}\"" "${CASE_DIR}/env"
cond_not "case3: no bin/python export" \
    grep -Fq 'bin/python' "${CASE_DIR}/env"

check_path_order() {
    (
        source "${CASE_DIR}/env"
        IFS=':' read -ra parts <<< "${PATH}"
        [[ "${parts[0]}" == "$1" ]] &&
            [[ "${parts[1]}" == "$1/Library/bin" ]] &&
            [[ "${parts[2]}" == "$1/Scripts" ]]
    )
}
cond "case3: sourced PATH is prefix, Library/bin, Scripts" \
    check_path_order "${PREFIX_POSIX}"
cond "case3: env file has four exports" \
    test "$(wc -l < "${CASE_DIR}/env")" -eq 4

# Case 4: Linux control keeps the three exports and no PATH line.
new_case case4
stub_network
install_micromamba_stub micromamba
run_setup globalopt-python/bin/python

cond "case4: setup exits zero" test "${rc}" -eq 0
cond "case4: env JULIA_LOAD_PATH" \
    grep -Fxq 'export JULIA_LOAD_PATH="@:@stdlib"' "${CASE_DIR}/env"
cond "case4: env JULIA_CONDAPKG_BACKEND" \
    grep -Fxq 'export JULIA_CONDAPKG_BACKEND="Null"' "${CASE_DIR}/env"
cond "case4: env JULIA_PYTHONCALL_EXE is POSIX python" \
    grep -Fxq "export JULIA_PYTHONCALL_EXE=\"${CASE_HOME}/.cache/sciml-benchmarks/globalopt-python/bin/python\"" \
    "${CASE_DIR}/env"
cond "case4: env file has three exports" \
    test "$(wc -l < "${CASE_DIR}/env")" -eq 3
cond_not "case4: no PATH export" \
    grep -Fq 'export PATH' "${CASE_DIR}/env"

if [[ "${FAILURES}" -gt 0 ]]; then
    printf '%d check(s) failed\n' "${FAILURES}" >&2
    exit 1
fi
echo "GlobalOptimization setup.sh tests passed"
