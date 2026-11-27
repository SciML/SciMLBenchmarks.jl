using Random
using BlackBoxOptimizationBenchmarking: bbob_suite, BenchmarkSetup, solve_problem
using OptimizationEvolutionary: DE
using OptimizationOptimJL: NelderMead

mutable struct BestObserved{F}
    f::F
    best::Float64
    count::Int
end

function (tracked::BestObserved)(u)
    value = tracked.f(u)
    tracked.best = min(tracked.best, value)
    tracked.count += 1
    return value
end

function run_stage(setup, f, budget; u0)
    tracked = BestObserved(f, Inf, 0)
    sol = solve_problem(setup, tracked, 3, budget; u0)
    fresh = f(sol.u)
    println("stage = ", setup.method)
    println("budget = ", budget)
    println("reported objective = ", sol.objective)
    println("fresh objective at returned point = ", fresh)
    println("best observed objective = ", tracked.best)
    println("objective evaluations = ", tracked.count)
    return sol, tracked.count, fresh, tracked.best
end

Random.seed!(4471)
f = bbob_suite(Val(3))[3]
u0 = [1.0, -1.0, 0.5]
de_setup = BenchmarkSetup(DE(), isboxed = true)
de_sol, de_count, de_fresh, de_best = run_stage(de_setup, f, 9; u0)
de_matches_point = de_sol.objective ≈ de_fresh
de_matches_best = de_sol.objective ≈ de_best
println("DE objective matches returned point = ", de_matches_point)
println("DE objective matches best observation = ", de_matches_best)

nm_setup = BenchmarkSetup(NelderMead(), isboxed = false)
nm_sol, nm_count, nm_fresh, nm_best = run_stage(nm_setup, f, 1; u0 = de_sol.u)
@assert nm_sol.objective ≈ nm_fresh
@assert nm_sol.objective ≈ nm_best
println("final objective = ", nm_sol.objective)
println("total objective evaluations = ", de_count + nm_count)

@assert de_matches_point
@assert de_matches_best
@assert nm_fresh <= de_fresh "Nelder–Mead worsened the fresh DE objective"
