# probe_phases.jl — diagnostic: split the end-to-end time into phases to find
# where the superlinear scaling comes from. NOT a benchmark; a one-off probe.
#
#   mpiexec -launcher fork -n <P> julia --project=. probe_phases.jl <N>
#
# Prints, on rank 0:
#   ranks,N,build_s,init_s,solve_s,solve2_s,iters
# where
#   build_s  = time to build the Julia SparseMatrixCSC
#   init_s   = time for SciMLBase.init  (PETSc assembly: MatSetValues + KSP setup)
#   solve_s  = time for the FIRST solve!  (includes any lazy setup)
#   solve2_s = time for a SECOND solve! on the same cache (pure solve, warm)
#
# If init_s is the term that blows up on 1 rank, the cost is matrix assembly
# (likely PETSc preallocation), not the linear solve.

using LinearAlgebra
using SparseArrays
using LinearSolve
using MPI
using PETSc
using SparseMatricesCSR
import SciMLBase

BLAS.set_num_threads(1)
MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NPROC = MPI.Comm_size(COMM)
const PETScExt = Base.get_extension(LinearSolve, :LinearSolvePETScExt)

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
    N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 40_000
    m = max(2, round(Int, sqrt(N)))
    alg = PETScAlgorithm(:cg; comm = COMM, pc_type = :none)

    # --- warm up compilation on a tiny system (excluded from timing) ---
    let As = laplacian_2d(8), bs = ones(64)
        c = SciMLBase.init(LinearProblem(As, bs), alg; abstol = 1e-10, reltol = 1e-10)
        solve!(c)
        PETScExt.cleanup_petsc_cache!(c)
    end
    MPI.Barrier(COMM)

    # --- phase 1: build Julia matrix ---
    local A, b
    build_s = @elapsed begin
        A = laplacian_2d(m)
        b = ones(size(A, 1))
    end

    # --- phase 2: init (PETSc assembly + KSP setup) ---
    local cache
    init_s = @elapsed begin
        cache = SciMLBase.init(LinearProblem(A, b), alg; abstol = 1e-10, reltol = 1e-10)
    end

    # --- phase 3: first solve ---
    local sol
    solve_s = @elapsed begin
        sol = solve!(cache)
    end
    iters = sol.iters

    # --- phase 4: second solve on same cache (warm, pure solve) ---
    solve2_s = @elapsed begin
        solve!(cache)
    end

    PETScExt.cleanup_petsc_cache!(cache)

    if RANK == 0
        println("$NPROC,$(size(A,1)),$build_s,$init_s,$solve_s,$solve2_s,$iters")
    end
    MPI.Barrier(COMM)
end

main()
