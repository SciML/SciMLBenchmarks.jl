# LinearSolveGPU — status

Both documents are written and validated end-to-end on a Tesla T4 (Modal, reduced
sizes, correctness gates green). The folder is **not yet CI-runnable**: no
Manifest can currently be resolved with `LinearSolve = "5"`, because of version
skew in the CUDA modular-stack migration:

```
LinearSolve 5  →  CUDSS 0.7  →  CUDACore ≥ 0.7
CUDA (4–5)     →  CUDACore ≤ 0.6.7        ← conflict
```

Unblock paths, in preference order:

1. A CUDA.jl release with CUDACore 0.7 compat (ecosystem timing — likely soon,
   since LinearSolve's own CUDA extension already targets the modular
   `cuSOLVER`/`CUDACore` stack).
2. Migrate these documents off `using CUDA` onto the modular packages the
   LinearSolve extension itself uses (`cuSOLVER`, `CUDACore`, `cuSPARSE`), and
   depend on those instead. Requires re-validating on a GPU.

Until then this folder intentionally ships without a Manifest so nobody burns a
demeter3 run discovering the conflict. Related: `CUDAOffload32MixedLUFactorization`
is unconditionally broken in the same migration (all nine of its GPU-array
references missed the rename) — fix validated on a T4, upstream report pending.
