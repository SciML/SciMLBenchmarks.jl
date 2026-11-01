#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_SCRIPT="${BUILD_SCRIPT:-${SCRIPT_DIR}/build_benchmark.sh}"
FIXTURE="$(mktemp -d)"
trap 'rm -rf "${FIXTURE}"' EXIT

mkdir -p "${FIXTURE}/benchmarks/Example" "${FIXTURE}/bin"
printf 'julia_version = "1.11.9"\n[[deps.StyledStrings]]\n' > "${FIXTURE}/Manifest.toml"
printf '%s\n' \
    '#!/bin/bash' \
    'set -euo pipefail' \
    "printf \"%s\\n\" \"\$*\" >> \"\${FAKE_JULIA_LOG}\"" \
    "if [[ \"\$*\" == *\"Pkg.instantiate\"* && -f Manifest.toml ]]; then" \
    '    echo "ERROR: Could not locate the source code for the StyledStrings package." >&2' \
    '    exit 1' \
    'fi' \
    "if [[ \"\$*\" == *\"Pkg.instantiate\"* ]]; then" \
    "    printf 'julia_version = \"1.10.11\"\\n' > Manifest.toml" \
    'fi' > "${FIXTURE}/bin/julia"
chmod +x "${FIXTURE}/bin/julia"

(
    cd "${FIXTURE}"
    PATH="${FIXTURE}/bin:${PATH}" FAKE_JULIA_LOG="${FIXTURE}/julia.log" \
        "${BUILD_SCRIPT}" benchmarks/Example
)

grep -Fq 'julia_version = "1.10.11"' "${FIXTURE}/Manifest.toml"
if grep -Fq 'StyledStrings' "${FIXTURE}/Manifest.toml"; then
    exit 1
fi
test "$(wc -l < "${FIXTURE}/julia.log")" -eq 2
grep -Fq -- '--project=. -e using Pkg; Pkg.instantiate(); Pkg.build()' \
    "${FIXTURE}/julia.log"
grep -Fq -- '--threads=auto --project=. benchmark.jl benchmarks/Example' \
    "${FIXTURE}/julia.log"

echo "build_benchmark.sh tests passed"

printf '[deps]\n' > "${FIXTURE}/benchmarks/Example/Project.toml"
printf 'julia_version = "1.11.9"\n' > "${FIXTURE}/benchmarks/Example/Manifest.toml"
cp "${FIXTURE}/benchmarks/Example/Project.toml" "${FIXTURE}/expected-project"
cp "${FIXTURE}/benchmarks/Example/Manifest.toml" "${FIXTURE}/expected-manifest"
cat > "${FIXTURE}/benchmarks/Example/setup.sh" <<'SETUP'
#!/bin/bash
set -euo pipefail
printf 'export SHARED_SETUP=ready\n' >> "${BENCHMARK_ENV_FILE}"
printf 'setup\n' >> "${FAKE_SETUP_LOG}"
SETUP
chmod +x "${FIXTURE}/benchmarks/Example/setup.sh"
cat > "${FIXTURE}/bin/julia" <<'JULIA'
#!/bin/bash
set -euo pipefail
[[ "${SHARED_SETUP}" == ready ]]
printf '%s\n' "$*" >> "${FAKE_JULIA_LOG}"
JULIA
chmod +x "${FIXTURE}/bin/julia"

for target in benchmarks/Example benchmarks/Example/page.jmd; do
    : > "${FIXTURE}/julia.log"
    : > "${FIXTURE}/setup.log"
    (
        cd "${FIXTURE}"
        PATH="${FIXTURE}/bin:${PATH}" FAKE_JULIA_LOG="${FIXTURE}/julia.log" \
            FAKE_SETUP_LOG="${FIXTURE}/setup.log" "${BUILD_SCRIPT}" "${target}"
    )
    test "$(wc -l < "${FIXTURE}/setup.log")" -eq 1
    test "$(wc -l < "${FIXTURE}/julia.log")" -eq 2
    grep -Fxq -- '--project=. -e using Pkg; Pkg.instantiate(); Pkg.build()' "${FIXTURE}/julia.log"
    grep -Fxq -- "--threads=auto --project=. benchmark.jl ${target}" "${FIXTURE}/julia.log"
    cmp "${FIXTURE}/expected-project" "${FIXTURE}/benchmarks/Example/Project.toml"
    cmp "${FIXTURE}/expected-manifest" "${FIXTURE}/benchmarks/Example/Manifest.toml"
done

echo 'Folder and file builds both run shared setup, propagate its environment, and preserve Project/Manifest.'
