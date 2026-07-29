# Findings — Distributed LinearSolve Benchmark Bring-up (Phase 0)

**Date:** 2026-07-17
**Author:** Jash (GSoC 2026)
**Context:** Standing up the `benchmarks/LinearSolveDistributed` suite to produce
performance evidence for the distributed PETSc/MPI linear-solver work in
LinearSolve.jl (branch `benchmark-1`, v5.0.0).

Phase 0's goal was a working, reproducible harness. It succeeded — and along the way
the benchmark surfaced two real bugs/limitations in the library and three methodology
traps that anyone benchmarking iterative distributed solvers will hit. This document
records them so the analysis isn't lost and so the eventual SciMLBenchmarks PR is
defensible.

---

## Environment

- LinearSolve.jl `benchmark-1` @ v5.0.0; PETSc.jl 0.4.10; MPI.jl 0.20.26
- `PETSc_jll` 3.22.1 (MPICH ABI, `ch4:ofi`); `MPICH_jll` 5.0.1; Julia 1.12.3
- Laptop: AMD Ryzen 7 4800H, 8 cores / 16 threads

---

## Correctness: confirmed working

The distributed solve path is **correct**. The replicated-`SparseMatrixCSC`
`PETScAlgorithm(comm = MPI.COMM_WORLD)` path solves and gathers the full solution on
every rank:

- 2-D FD Laplacian, GMRES+Jacobi: residual ~1e-11, `Success`, on 1 and 2 ranks.
- GAMG-CG: converges in **10 iterations** (textbook multigrid), residual drops 6 orders
  of magnitude — verified via PETSc's `ksp_monitor_true_residual`.
- Direct LU (`preonly`+`lu`) on 10k: solves in <1s, residual 3e-13.

So none of the issues below are correctness problems in the solve.

---

## Finding 1 — `solve!` discards iteration count on safety-check failure *(library bug)*

**Where:** `ext/LinearSolvePETScExt.jl`, `SciMLBase.solve!`, ~line 822–837.

**What:** The convergence metadata (`iters`, `resid`) is computed correctly at line
822–824. But when the a-posteriori residual safety check fails (line 832–834), the
function early-returns the safety-check's own solution object, which reports `iters = 0`
— discarding the real iteration count.

**Reproduction:** GAMG-CG on the 40k Laplacian with `reltol = 1e-8`. PETSc's log says
`CONVERGED_RTOL iterations 10`, but `sol.iters == 0` and
`sol.retcode == APosterioriSafetyFailure`.

**Why it happens:** preconditioned CG converges in the *preconditioned* residual norm,
so it declares convergence when that norm hits `rtol` — but the *true* residual is
looser (~1e-6 here). LinearSolve's safety check uses the true residual, correctly
notices it misses `rtol = 1e-8`, and rejects. The rejection path just doesn't carry the
iteration count forward.

**Suggested fix:** preserve `iters`/`resid` on the safety-failure return so users retain
convergence diagnostics. Small, self-contained PR against LinearSolve.jl.

**Benchmark workaround:** set `ksp_norm_type = "unpreconditioned"` so PETSc's stopping
test uses the true residual, aligning it with the safety check (applied in
`run_solve.jl`).

---

## Finding 2 — Replicated-matrix path does not scale in memory *(design limitation)*

**What:** The Weeks 1–2 deliverable assembles a replicated `SparseMatrixCSC` — **every
rank holds the full matrix** before row-owning its slice into PETSc. Memory therefore
grows with rank count.

**Symptom:** 8-rank GAMG solve of the 40k system was OOM-killed (signal 9) on a 16 GB
laptop — 8 full-matrix copies plus 8 GAMG hierarchies.

**Implication for benchmarking:** on a single machine there is **no problem size that
gives a clean 1→8 strong-scaling curve** for this path — small `N` is overhead-bound,
large `N` runs out of RAM. Honest strong/weak scaling for the *replicated* path needs
real multi-node hardware (amdci cluster) where each node has its own RAM. For
memory-scalable single-machine scaling, the `PSparseMatrix`/`PVector` path (each rank
owns only its slice) is the right target — a candidate for a future benchmark variant.

**This is not a bug** — it's inherent to the replicated design, which trades memory for
API simplicity (users hand in a plain Julia sparse matrix). Worth documenting as the
explicit trade-off it is.

---

## Finding 3 — Methodology traps (for the record)

Three ways the benchmark *looked* like it was working while measuring the wrong thing.
Each is now guarded or documented in the suite.

1. **Rank-dependent preconditioner fakes superlinear speedup.** Block-Jacobi's block
   structure follows the row partition, so iteration count changes with rank count —
   turning "same problem, more cores" into "different algorithm per rank." First run
   showed 24× speedup / 300% efficiency on 8 ranks. Fix: use a preconditioner whose
   strength is (near) partition-independent, and **always report iteration counts** so
   the invariance is auditable.

2. **Slow serial baseline inflates all speedups.** With unpreconditioned CG, the
   1-rank solve took 59s for 416 iterations (~143 ms/iter for a 40k sparse matvec —
   ~1000× slower than optimized PETSc). The "speedup" was largely parallelism chipping
   away at an overhead-bound baseline, not real scaling. Root cause: the JLL PETSc build
   is not performance-optimized, and unpreconditioned CG needs `O(√N)` iterations. Fix:
   use a scalable preconditioner (GAMG → 10 iters) so per-iteration overhead stops
   dominating.

3. **Preconditioned vs. true residual mismatch.** See Finding 1 — the stopping norm and
   the correctness check must agree, or converged solves get flagged as failures.

**Guard added:** `StrongScaling.jmd` now prints iteration counts every run and emits
`@warn` when (a) GAMG iteration spread across ranks exceeds 25%, or (b) efficiency
exceeds 110%. The document refuses to silently present a physically impossible curve.

---

## Toolchain notes (non-obvious, will recur on every runner)

- **`-launcher fork` is mandatory.** MPICH_jll 5.x's hydra fails to bootstrap PMI on a
  single node (`PMI_Init returned -1`); `-launcher fork` spawns ranks with `fork()`
  instead. System MPI here is OpenMPI 4.1.6, unusable because PETSc_jll needs the MPICH
  ABI — so `MPIPreferences.use_system_binary()` is not an option.
- **Never call `MPI.Finalize()` in a worker.** PETSc's GC finalizers touch MPI during
  teardown; MPI.jl's atexit hook finalizes MPI *after* them. An explicit `Finalize()`
  finalizes too early → "MPI routine after finalizing MPICH" → nonzero exit on an
  otherwise-successful solve.
- **`OMP_NUM_THREADS=1` + `BLAS.set_num_threads(1)`** so parallelism comes only from
  rank count, not per-rank thread oversubscription.

---

## What's next

- **Scale on real hardware.** Run `StrongScaling.jmd` at large `N` and higher `RANKS` on
  amdci (canonical, production PETSc, per-node RAM) via SciMLBenchmarks CI, and/or a
  large Modal container for fast iteration. That is where the honest "marketing" curve
  comes from — the laptop's job (harness + correctness + finding bugs) is done.
- **Consider the `PSparseMatrix` variant** for memory-scalable single-machine scaling.
- **File Finding 1** as a LinearSolve.jl issue / small PR.
