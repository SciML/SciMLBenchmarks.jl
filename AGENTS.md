# Agent instructions

## Validating benchmark changes

For changes confined to benchmark sources or their environments under `benchmarks/`, the authoritative test is that every affected benchmark Weaves without error. Run the same entrypoint used by CI:

```sh
julia --threads=auto --project=. benchmark.jl benchmarks/<folder>/<benchmark>.jmd
```

Run the folder instead when a dependency or environment change can affect multiple pages:

```sh
julia --threads=auto --project=. benchmark.jl benchmarks/<folder>
```

If the benchmark directory has a `setup.sh`, reproduce the complete CI sequence with `.github/scripts/build_benchmark.sh <target>`.

Do not add tests under `test/` that merely inspect benchmark source text or assert dependency strings in a benchmark `Project.toml` or `Manifest.toml`. Those checks do not establish that the benchmark executes successfully. Add unit tests when changing the shared benchmark harness or another independently testable code path; otherwise, report the exact Weave command and result in the pull request.
