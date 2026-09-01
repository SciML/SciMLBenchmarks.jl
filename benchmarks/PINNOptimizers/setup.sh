#!/bin/bash
set -euo pipefail

benchmark_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
printf '[CUDA_Runtime_jll]\n__clear__ = ["local"]\nversion = "12.9"\n' > "${benchmark_dir}/LocalPreferences.toml"
