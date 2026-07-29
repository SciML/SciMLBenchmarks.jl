# probe_gamg.jl — why does GAMG report 0 iters / APosterioriSafetyFailure?
# Uses the SAME 2-D 5-point Laplacian as run_solve.jl and turns on PETSc's own
# per-iteration monitor + converged-reason so we can see what the solver does.
#
#   mpiexec -launcher fork -n 1 julia --project=. probe_gamg.jl

using LinearAlgebra, SparseArrays, LinearSolve, MPI, PETSc, SparseMatricesCSR
import SciMLBase

BLAS.set_num_threads(1)
MPI.Init()
const COMM = MPI.COMM_WORLD
const PETScExt = Base.get_extension(LinearSolve, :LinearSolvePETScExt)

function laplacian_2d(m::Int)
    n = m * m
    I = Int[]; J = Int[]; V = Float64[]
    lin(i, j) = (j - 1) * m + i
    @inbounds for j in 1:m, i in 1:m
        k = lin(i, j)
        push!(I, k); push!(J, k); push!(V, 4.0)
        i > 1 && (push!(I, k); push!(J, lin(i-1, j)); push!(V, -1.0))
        i < m && (push!(I, k); push!(J, lin(i+1, j)); push!(V, -1.0))
        j > 1 && (push!(I, k); push!(J, lin(i, j-1)); push!(V, -1.0))
        j < m && (push!(I, k); push!(J, lin(i, j+1)); push!(V, -1.0))
    end
    return sparse(I, J, V, n, n)
end

A = laplacian_2d(200)          # 40k unknowns, same as the failing case
b = ones(size(A, 1))

for (label, pc, opts) in [
        ("gamg + monitor", :gamg, (ksp_monitor_true_residual="", ksp_converged_reason="")),
        ("hypre? (skip if missing)", :gamg, (ksp_converged_reason="", pc_gamg_type="agg", pc_gamg_agg_nsmooths=1)),
    ]
    println("\n=== $label ===")
    alg = PETScAlgorithm(:cg; comm=COMM, pc_type=pc, ksp_options=opts)
    c = SciMLBase.init(LinearProblem(A, b), alg; abstol=1e-12, reltol=1e-8)
    s = solve!(c)
    res = norm(A*s.u - b)/norm(b)
    println(">> iters=", s.iters, "  retcode=", s.retcode, "  rel_res=", res)
    PETScExt.cleanup_petsc_cache!(c)
end
