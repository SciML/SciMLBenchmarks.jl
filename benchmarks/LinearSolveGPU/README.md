# LinearSolveGPU

GPU benchmarks for LinearSolve.jl, running on the GPU runner (see
`benchmark_config.toml`).

- `DenseGPUOffload.jmd`: CPU LU vs GPU LU/QR/32-mixed offload across sizes, with
  the PCIe round-trip measured separately and an explicit CPU/GPU crossover N.
- `SparseGPU.jmd`: CUDSS sparse direct (via `LUFactorization` on
  `CuSparseMatrixCSR`) vs UMFPACK/KLU, factor vs cached re-solve split, with the
  CSR host-to-device transfer reported explicitly.

History note: this folder initially shipped without a Manifest because no
version set satisfied LinearSolve 5 + CUDSS 0.7 + the pre-6.0 CUDA stack
(CUDACore version conflict). The CUDA 6 modular-stack releases resolved it;
the committed Manifest pins the working combination.
