using Random
using BlackBoxOptimizationBenchmarking: bbob_suite, BenchmarkSetup,
    FunctionCallsCounter, solve_problem
using OptimizationEvolutionary: DE
using OptimizationOptimJL: NelderMead
const ForwardDiff = Base.require(
    Base.PkgId(Base.UUID("f6369f11-7733-5829-9624-2563aa707210"), "ForwardDiff")
)

mutable struct BestObserved{F}
    f::F
    best::Base.RefValue{Float64}
end

function (tracked::BestObserved)(u)
    value = tracked.f(u)
    tracked.best[] = min(tracked.best[], value)
    return value
end

function run_stage(setup, f, budget; u0)
    best = Ref(Inf)
    tracked = FunctionCallsCounter(BestObserved(f, best))
    sol = solve_problem(setup, tracked, 3, budget; u0)
    fresh = f(sol.u)
    println("stage = ", setup.method)
    println("budget = ", budget)
    println("reported objective = ", sol.objective)
    println("fresh objective at returned point = ", fresh)
    println("best observed objective = ", best[])
    println("objective evaluations = ", tracked.count)
    return sol, tracked, fresh, best[]
end

Random.seed!(4471)
f = bbob_suite(Val(3))[3]
u0 = [1.0, -1.0, 0.5]
de_setup = BenchmarkSetup(DE(), isboxed = true)
de_sol, de_counter, de_fresh, de_best = run_stage(de_setup, f, 9; u0)
de_matches_point = de_sol.objective ≈ de_fresh
de_matches_best = de_sol.objective ≈ de_best
println("DE objective matches returned point = ", de_matches_point)
println("DE objective matches best observation = ", de_matches_best)

nm_setup = BenchmarkSetup(NelderMead(), isboxed = false)
nm_sol, nm_counter, nm_fresh, nm_best = run_stage(nm_setup, f, 1; u0 = de_sol.u)
@assert nm_sol.objective ≈ nm_fresh
@assert nm_sol.objective ≈ nm_best
println("final objective = ", nm_sol.objective)
println("total objective evaluations = ", de_counter.count + nm_counter.count)

@assert de_matches_point
@assert de_matches_best
@assert nm_sol.objective ≈ 174.910943039 atol = 1.0e-9
