# Agent instructions

## Validating benchmark changes

For changes confined to benchmark sources or their environments under `benchmarks/`, the authoritative test is that every affected benchmark Weaves without error. Run the CI build script to apply the configured Julia thread count and folder setup:

```sh
.github/scripts/build_benchmark.sh benchmarks/<folder>/<benchmark>.jmd
```

Run the folder instead when a dependency or environment change can affect multiple pages:

```sh
.github/scripts/build_benchmark.sh benchmarks/<folder>
```

Do not add tests under `test/` that merely inspect benchmark source text or assert dependency strings in a benchmark `Project.toml` or `Manifest.toml`. Those checks do not establish that the benchmark executes successfully. Add unit tests when changing the shared benchmark harness or another independently testable code path; otherwise, report the exact Weave command and result in the pull request.
