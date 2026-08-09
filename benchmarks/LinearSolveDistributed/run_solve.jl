# run_solve.jl — MPI worker for the distributed LinearSolve benchmarks.
#
# Launched as:  mpiexec -n <P> julia --project=. run_solve.jl <N> [solver] [pc]
#
#   <N>      target number of unknowns (a 2-D FD Laplacian of side ~sqrt(N))
#   solver   PETSc KSP type: gmres (default), cg, bcgs, ...
#   pc       preconditioner: none (default), jacobi, gamg, ilu, ...
#
# The matrix is a replicated SparseMatrixCSC assembled identically on every rank;
# PETScAlgorithm(comm = MPI.COMM_WORLD) then row-owns it across the communicator,
# solves, and gathers the full solution back into sol.u on every rank (this is the
# replicated-SparseMatrixCSC path). Rank 0 prints one CSV line:
#
#   ranks,N,nnz,solver,pc,time_s,residual,retcode
#
# so the orchestrating .jmd can parse it directly. All correctness checking
# (residual < tol) happens here, on the worker, and a failure exits nonzero.
#
#   ranks,N,nnz,solver,pc,time_s,residual,iters,retcode

using LinearAlgebra
using SparseArrays
using LinearSolve
using MPI
using PETSc
using SparseMatricesCSR   # loaded so the PETSc MPI extension is active
using BenchmarkTools
import SciMLBase

# Pin every rank to a single BLAS/LAPACK thread. Parallelism in this benchmark must
# come ONLY from the number of MPI ranks — not from each rank silently spawning BLAS
# threads across all cores. Without this the -n 1 baseline oversubscribes all cores
# on a bandwidth-bound sparse solve (cache/NUMA thrash), runs pathologically slow,
# and manufactures fake "superlinear" speedup vs the multi-rank runs. One thread per
# rank makes the serial baseline honest, so speedup reflects real parallel work.
# (OMP_NUM_THREADS is also exported by the caller as a belt-and-suspenders guard.)
BLAS.set_num_threads(1)

MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NPROC = MPI.Comm_size(COMM)

const PETScExt = Base.get_extension(LinearSolve, :LinearSolvePETScExt)

# 2-D 5-point Laplacian on an m×m grid → SPD, size m^2. Strictly the negative
# discrete Laplacian shifted to be diagonally dominant so simple PCs converge.
function laplacian_2d(m::Int)
    n = m * m
    I = Int[]; J = Int[]; V = Float64[]
    lin(i, j) = (j - 1) * m + i
    @inbounds for j in 1:m, i in 1:m
        k = lin(i, j)
        push!(I, k); push!(J, k); push!(V, 4.0)
        if i > 1; push!(I, k); push!(J, lin(i - 1, j)); push!(V, -1.0); end
        if i < m; push!(I, k); push!(J, lin(i + 1, j)); push!(V, -1.0); end
        if j > 1; push!(I, k); push!(J, lin(i, j - 1)); push!(V, -1.0); end
        if j < m; push!(I, k); push!(J, lin(i, j + 1)); push!(V, -1.0); end
    end
    return sparse(I, J, V, n, n)
end

function main()
    N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10_000
    solver = length(ARGS) >= 2 ? Symbol(ARGS[2]) : :gmres
    pc = length(ARGS) >= 3 ? Symbol(ARGS[3]) : :none

    m = max(2, round(Int, sqrt(N)))
    A = laplacian_2d(m)
    n = size(A, 1)
    b = ones(n)

    # ksp_norm_type=unpreconditioned makes PETSc's convergence test use the TRUE
    # residual, not the preconditioned one. Without this, preconditioned CG (e.g.
    # under GAMG) declares CONVERGED_RTOL when the *preconditioned* norm hits
    # reltol, but the true residual is looser — and LinearSolve's a-posteriori
    # safety check (which uses the true residual) then rejects the solve. Forcing
    # the unpreconditioned norm aligns PETSc's stopping test with the safety check.
    alg = PETScAlgorithm(solver; comm = COMM, pc_type = pc,
                         ksp_options = (ksp_norm_type = "unpreconditioned",))

    # Tolerances are env-overridable. Default reltol is deliberately loose (1e-8)
    # because preconditioned Krylov methods (e.g. GAMG-CG) converge in the
    # *preconditioned* residual norm, so the *true* relative residual they reach
    # is looser — asking for 1e-10 true residual makes LinearSolve's a-posteriori
    # safety check reject an otherwise-converged solve. 1e-8 is comfortably
    # achievable and still a tight benchmark target.
    rtol = parse(Float64, get(ENV, "BENCH_RELTOL", "1e-8"))
    atol = parse(Float64, get(ENV, "BENCH_ABSTOL", "1e-12"))
    # Correctness gate for the warm-up: pass if the true residual is within ~100×
    # of the requested reltol (accounts for the preconditioned-vs-true gap above).
    gate = 100 * rtol

    # Warm-up solve (excluded from timing): triggers compilation + PETSc setup.
    let cache = SciMLBase.init(LinearProblem(A, b), alg; abstol = atol, reltol = rtol)
        sol = solve!(cache)
        PETScExt.cleanup_petsc_cache!(cache)
        # Correctness gate — a fast wrong answer must not be reported.
        res = norm(A * sol.u - b) / norm(b)
        if !(res < gate)
            RANK == 0 && println(stderr, "RESIDUAL FAILURE: $res (gate $gate)")
            MPI.Barrier(COMM)
            exit(2)
        end
    end

    # Timed solve: minimum wall-clock over samples, setup rebuilds each sample so
    # we measure the full assemble→solve→gather cost (the honest end-to-end number).
    # Budget is env-overridable: keep it small for fast local iteration; the
    # dedicated benchmark runner can raise BENCH_SECONDS/BENCH_SAMPLES for tighter
    # statistics on the published numbers. Strong-scaling *ratios* (T₁/T_P) are
    # stable even with a short budget, so the default is deliberately cheap.
    bench_seconds = parse(Float64, get(ENV, "BENCH_SECONDS", "5"))
    bench_samples = parse(Int, get(ENV, "BENCH_SAMPLES", "3"))

    t = @belapsed begin
        c = SciMLBase.init(LinearProblem($A, $b), $alg; abstol = $atol, reltol = $rtol)
        s = solve!(c)
        $PETScExt.cleanup_petsc_cache!(c)
        s
    end samples=bench_samples seconds=bench_seconds evals=1

    # Recompute a residual + retcode + iteration count on a fresh solve for the
    # report line. sol.iters matters a lot for iterative solvers: the iteration
    # count is the honest measure of solver work, and it must stay ~constant
    # across ranks for a strong-scaling ratio to be meaningful — so we surface it,
    # never hide it.
    cache = SciMLBase.init(LinearProblem(A, b), alg; abstol = atol, reltol = rtol)
    sol = solve!(cache)
    res = norm(A * sol.u - b) / norm(b)
    rc = sol.retcode
    iters = sol.iters
    PETScExt.cleanup_petsc_cache!(cache)

    if RANK == 0
        println("$NPROC,$n,$(nnz(A)),$solver,$pc,$t,$res,$iters,$rc")
    end

    # NB: do NOT call MPI.Finalize() here. PETSc registers GC finalizers that
    # touch MPI during teardown; MPI.jl's own atexit hook finalizes MPI after
    # those run. Calling Finalize() ourselves finalizes MPI too early, so the
    # trailing GC hits "MPI routine after finalizing MPICH" and the rank exits
    # nonzero even though the solve succeeded. Let atexit handle it.
    MPI.Barrier(COMM)
end

main()
