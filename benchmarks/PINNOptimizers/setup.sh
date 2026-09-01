#!/bin/bash
set -euo pipefail

benchmark_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project="${benchmark_dir}/Project.toml"
if ! grep -q '^CUDA_Runtime_jll = ' "${project}"; then
    if grep -q '^\[extras\]' "${project}"; then
        sed -i '/^\[extras\]/a CUDA_Runtime_jll = "76a88914-d11a-5bdc-97e0-2f5a05c973a2"' "${project}"
    else
        printf '\n[extras]\nCUDA_Runtime_jll = "76a88914-d11a-5bdc-97e0-2f5a05c973a2"\n' >> "${project}"
    fi
fi
printf '[CUDA_Runtime_jll]\n__clear__ = ["local"]\nversion = "12.9"\n' > "${benchmark_dir}/LocalPreferences.toml"
