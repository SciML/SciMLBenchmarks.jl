using Printf
using Logging
using DataFrames
using Plots
using StatsPlots
using StatsBase: countmap
using Statistics

using Optimization
using NLPModels
using OptimizationNLPModels
using OptimizationOptimJL
using OptimizationOptimJL: LBFGS, ConjugateGradient, NelderMead, SimulatedAnnealing,
    ParticleSwarm, BFGS, Newton, NewtonTrustRegion
using OptimizationOptimisers: Optimisers
using OptimizationMOI
using OptimizationMOI: MOI
using Ipopt
using OptimizationNLopt
using OptimizationNLopt: NLopt
using OptimizationMadNLP
using MadNLP

const MAX_PROBLEMS_PER_CATEGORY = 50
const MAX_NVAR = 1_000
const MAX_NCON = 1_000
const SOLVE_MAXITERS = 1_000
const SOLVE_TIMEOUT_SECONDS = 90.0

# A category is considered degenerate (and the page fails) when fewer than this
# fraction of its rows complete at the harness level (`status == "OK"`), regardless
# of the solver-reported return code.
const MIN_COMPLETED_FRACTION = 0.8

# Tiny problems used to trigger JIT compilation before the first measured solve.
const WARMUP_UNCONSTRAINED_PROBLEM = "ROSENBR"
const WARMUP_CONSTRAINED_PROBLEM = "HS35"

const SUCCESS_RETCODES = Set(["Success", "Terminated", "FirstOrderOptimal"])

# Solution-quality tolerances used by the "verified" success criterion. A run is a
# verified success when its return code is in SUCCESS_RETCODES, the returned point
# violates no constraint or bound by more than FEAS_TOL, and either the projected
# gradient norm is below OPT_TOL or the objective is within GAP_TOL (relative) of the
# best feasible objective any solver found for that problem.
const FEAS_TOL = 1.0e-6
const OPT_TOL = 1.0e-5
const GAP_TOL = 1.0e-6

const KNOWN_BAD_PROBLEMS = Set(
    lowercase.(
        String[
            "BLOWEYA", "CHARDIS1", "CLEUVEN4", "CMPC3", "CMPC10", "CVXQP2",
            "DITTERT", "HIER13", "LUKVLE8", "LUKVLI7", "MPC2", "PATTERNNE",
            "READING2", "READING6", "NINENEW", "MSS1",
        ],
    ),
)

const UNCONSTRAINED_SOLVERS = [
    "LBFGS",
    "BFGS",
    "NLopt_LD_LBFGS",
    "ConjugateGradient",
    "Newton",
    "NewtonTrustRegion",
    "NelderMead",
    "NLopt_LN_BOBYQA",
    "SimulatedAnnealing",
    "ParticleSwarm",
    "Adam",
]
const CONSTRAINED_SOLVERS = ["Ipopt", "MadNLP", "NLopt-SLSQP", "NLopt-AUGLAG-LBFGS"]
const NLOPT_SOLVERS = Set(["NLopt-SLSQP", "NLopt-AUGLAG-LBFGS"])
const NLOPT_RELTOL = 1.0e-8

# Global heuristics run until `maxiters` and Optimisers.jl rules run for exactly
# `maxiters` gradient steps; neither reports convergence through the return code, so
# their rows are summarized separately from the local convergence-based solvers.
const SOLVER_CLASSES = Dict(
    "LBFGS" => "quasi-Newton",
    "BFGS" => "quasi-Newton",
    "NLopt_LD_LBFGS" => "quasi-Newton",
    "ConjugateGradient" => "conjugate gradient",
    "Newton" => "Newton",
    "NewtonTrustRegion" => "Newton",
    "NelderMead" => "derivative-free",
    "NLopt_LN_BOBYQA" => "derivative-free",
    "SimulatedAnnealing" => "global heuristic",
    "ParticleSwarm" => "global heuristic",
    "Adam" => "first-order (fixed budget)",
    "Ipopt" => "interior point",
    "MadNLP" => "interior point",
    "NLopt-SLSQP" => "sequential quadratic",
    "NLopt-AUGLAG-LBFGS" => "augmented Lagrangian",
)
const BUDGET_LIMITED_CLASSES = Set(["global heuristic", "first-order (fixed budget)"])

solver_class(name) = get(SOLVER_CLASSES, name) do
    error("Unknown optimizer class: $name")
end

convergence_based(class) = !(class in BUDGET_LIMITED_CLASSES)

function optimizer_from_name(name)
    if name == "LBFGS"
        return LBFGS()
    elseif name == "BFGS"
        return BFGS()
    elseif name == "NLopt_LD_LBFGS"
        return NLopt.LD_LBFGS()
    elseif name == "ConjugateGradient"
        return ConjugateGradient()
    elseif name == "Newton"
        return Newton()
    elseif name == "NewtonTrustRegion"
        return NewtonTrustRegion()
    elseif name == "NelderMead"
        return NelderMead()
    elseif name == "NLopt_LN_BOBYQA"
        return NLopt.LN_BOBYQA()
    elseif name == "SimulatedAnnealing"
        return SimulatedAnnealing()
    elseif name == "ParticleSwarm"
        return ParticleSwarm()
    elseif name == "Adam"
        return Optimisers.Adam()
    elseif name == "Ipopt"
        return MOI.OptimizerWithAttributes(
            Ipopt.Optimizer,
            "max_iter" => SOLVE_MAXITERS,
            "max_wall_time" => SOLVE_TIMEOUT_SECONDS,
            "hessian_approximation" => "limited-memory",
            "tol" => 1.0e-6,
            "print_level" => 0,
        )
    elseif name == "MadNLP"
        # CompactLBFGS matches Ipopt's limited-memory Hessian path and avoids
        # ExactHessian, which requires a Lagrangian Hessian OptimizationNLPModels
        # does not currently expose.
        return MadNLPOptimizer(
            hessian_approximation = MadNLP.CompactLBFGS,
            acceptable_tol = 1.0e-6,
            additional_options = Dict{Symbol, Any}(
                :print_level => MadNLP.ERROR,
            ),
        )
    elseif name == "NLopt-SLSQP"
        return NLopt.LD_SLSQP()
    elseif name == "NLopt-AUGLAG-LBFGS"
        return NLopt.AUGLAG()
    else
        error("Unknown optimizer: $name")
    end
end

# NLopt has no default stopping tolerance, so without `reltol` every NLopt run ends in
# `MaxIters`; AUGLAG additionally requires a local optimizer.
function solve_kwargs_from_name(name)
    if name == "NLopt-SLSQP"
        return (; reltol = NLOPT_RELTOL)
    elseif name == "NLopt-AUGLAG-LBFGS"
        return (;
            reltol = NLOPT_RELTOL, local_method = NLopt.LD_LBFGS(),
            local_maxiters = SOLVE_MAXITERS,
        )
    else
        return (;)
    end
end

# OptimizationNLopt hands constraints to NLopt as `cons(x) <= 0` (inequalities) and
# `cons(x) == 0` (equalities) and ignores `lcons`/`ucons`, so NLopt only solves the
# intended problem when every CUTEst constraint is already in that form.
function nlopt_constraint_form_ok(meta)
    for (l, u) in zip(meta.lcon, meta.ucon)
        l == u == 0 && continue
        l == -Inf && u == 0 && continue
        return false
    end
    return true
end

# Returns `nothing` when `solver_name` can be applied to a problem with metadata `meta`,
# otherwise a short reason recorded alongside the `SKIPPED` status.
function solver_skip_reason(solver_name, meta)
    if solver_name in NLOPT_SOLVERS && !nlopt_constraint_form_ok(meta)
        return "constraint bounds are not of the cons(x) <= 0 / cons(x) == 0 form " *
            "that OptimizationNLopt passes to NLopt"
    end
    # SLSQP's least-squares subproblem needs at most nvar equality constraints; NLopt
    # otherwise returns INVALID_ARGS before iterating.
    if solver_name == "NLopt-SLSQP" && count(meta.lcon .== meta.ucon) > meta.nvar
        return "SLSQP supports at most nvar equality constraints"
    end
    return nothing
end

function problem_metadata(name)
    nlp = nothing
    try
        nlp = CUTEstModel(name)
        return (; ok = true, nvar = nlp.meta.nvar, ncon = nlp.meta.ncon)
    catch err
        @warn "Unable to load CUTEst problem metadata" problem = name exception = (
            err, catch_backtrace(),
        )
        return (; ok = false, nvar = -1, ncon = -1)
    finally
        if nlp !== nothing
            try
                finalize(nlp)
            catch err
                @warn "Unable to finalize CUTEst problem metadata" problem = name exception = (
                    err, catch_backtrace(),
                )
            end
        end
    end
end

function select_safe_problems(
        candidates;
        max_problems = MAX_PROBLEMS_PER_CATEGORY,
        max_var = MAX_NVAR,
        max_con = MAX_NCON,
    )
    selected = String[]

    for name in candidates
        lowercase(name) in KNOWN_BAD_PROBLEMS && continue

        meta = problem_metadata(name)

        meta.ok || continue
        meta.nvar <= max_var || continue
        meta.ncon <= max_con || continue

        push!(selected, name)

        length(selected) >= max_problems && break
    end

    return selected
end

# Solver-reported time (`sol.stats.time`). Its definition differs between backends
# (some include setup, some only the iteration loop, some never set it), so it is
# recorded as `solver_reported_secs` for reference only and is NOT the primary metric.
function solve_seconds(sol, fallback)
    try
        if hasfield(typeof(sol), :stats) && hasfield(typeof(sol.stats), :time)
            secs = Float64(sol.stats.time)
            isfinite(secs) && secs >= 0 && return secs
        end
    catch
    end

    return fallback
end

function retcode_name(retcode)
    name = string(retcode)
    return startswith(name, "ReturnCode.") ? last(split(name, '.')) : name
end

function print_exception(prefix, err, bt)
    return println(prefix, ": ", sprint(showerror, err, bt))
end

elapsed_seconds(started_ns) = (time_ns() - started_ns) / 1.0e9

# `Base.cumulative_compile_time_ns()` returns `(compile_ns, recompile_ns)` on Julia >= 1.9
# (a bare integer on older versions); recompilation is a subset of compilation, so the
# first element is the total compile time (this is what `@time` reports).
compile_ns(t) = t isa Tuple ? first(t) : t
compile_seconds_since(before) = (compile_ns(Base.cumulative_compile_time_ns()) - compile_ns(before)) / 1.0e9

function solve_once(prob, solver_name; maxtime = SOLVE_TIMEOUT_SECONDS)
    return solve(
        prob, optimizer_from_name(solver_name);
        maxiters = SOLVE_MAXITERS,
        maxtime = maxtime,
        solve_kwargs_from_name(solver_name)...,
    )
end

# The second (clean) solve is skipped when the first one already hit the time cap or
# errored, so a timed-out pair costs one cap, not two.
rerun_worthwhile(first_retcode, first_solve_secs; maxtime = SOLVE_TIMEOUT_SECONDS) =
    first_retcode != "MaxTime" && first_solve_secs < maxtime

function build_problem(nlp, solver_name)
    if solver_class(solver_name) != "first-order (fixed budget)"
        return OptimizationNLPModels.OptimizationProblem(nlp)
    end

    # OptimizationOptimisers evaluates the objective through `fg` (or a one-argument
    # `f`), neither of which the NLPModels wrapper provides, so supply `objgrad!` here.
    fg(G, x, p) = first(NLPModels.objgrad!(nlp, x, G))
    return OptimizationProblem(
        OptimizationNLPModels.OptimizationFunction(nlp; fg),
        nlp.meta.x0,
    )
end

const NO_QUALITY = (; objective = NaN, grad_norm = NaN, cons_viol = NaN)

# Maximum violation of `lo <= v <= hi`, elementwise; 0.0 when everything is within bounds.
function bound_violation(v, lo, hi)
    viol = 0.0
    for i in eachindex(v)
        viol = max(viol, lo[i] - v[i], v[i] - hi[i])
    end
    return viol
end

# Evaluate the objective, the projected gradient norm and the maximum constraint/bound
# violation of the NLPModels object at `x`. For unconstrained problems the projected
# gradient reduces to the infinity norm of the gradient; for problems with general
# constraints it is only indicative and `cons_viol` is the primary check. Every metric is
# guarded independently so an evaluation failure yields `NaN` for that metric only.
function solution_quality(nlp, x)
    x = collect(Float64, x)
    meta = nlp.meta

    objective = try
        Float64(NLPModels.obj(nlp, x))
    catch
        NaN
    end

    grad_norm = try
        g = NLPModels.grad(nlp, x)
        proj = clamp.(x .- g, meta.lvar, meta.uvar)
        maximum(abs, x .- proj; init = 0.0)
    catch
        NaN
    end

    cons_viol = try
        viol = bound_violation(x, meta.lvar, meta.uvar)
        if meta.ncon > 0
            c = NLPModels.cons(nlp, x)
            viol = max(viol, bound_violation(c, meta.lcon, meta.ucon))
        end
        viol
    catch
        NaN
    end

    return (; objective, grad_norm, cons_viol)
end

# Each (problem, solver) pair is solved twice. Timing columns:
#   first_solve_secs     wall-clock time of the first `solve(...)`, includes any compilation
#   compile_secs         Julia compilation time measured during the first solve
#                        (`Base.cumulative_compile_timing`). Solver-level compilation is
#                        absorbed by `warmup_solvers`, so this column shows per-problem-shape
#                        recompilation (bounds / constraints present or not) if any
#   secs                 wall-clock time of the second `solve(...)` on a fresh
#                        `OptimizationProblem` from the same `nlp` (primary metric); equals
#                        `first_solve_secs` when the rerun is skipped
#   decode_secs          wall-clock time of the CUTEst SIF decode / `CUTEstModel(name)`
#   solver_reported_secs `sol.stats.time` when the backend provides it, otherwise NaN
# `retcode` comes from the second run, `first_retcode` from the first.
# Quality columns (`objective`, `grad_norm`, `cons_viol`) are evaluated at the final
# `sol.u` after the clean (or only) solve; NaN on load/solve failure.
function run_single_solve(problem_name, solver_name; maxtime = SOLVE_TIMEOUT_SECONDS)
    nlp = nothing
    decode_started = time_ns()
    class = solver_class(solver_name)

    try
        nlp = CUTEstModel(problem_name)
    catch err
        bt = catch_backtrace()
        print_exception("CUTEst problem load failed for $problem_name", err, bt)

        return (;
            problem = problem_name,
            solver = solver_name,
            solver_class = class,
            n_vars = -1,
            secs = 0.0,
            first_solve_secs = 0.0,
            compile_secs = 0.0,
            decode_secs = elapsed_seconds(decode_started),
            solver_reported_secs = NaN,
            retcode = "LOAD_FAILED",
            first_retcode = "LOAD_FAILED",
            status = "LOAD_FAILED",
            NO_QUALITY...,
        )
    end
    decode_secs = elapsed_seconds(decode_started)

    first_solve_secs = 0.0
    compile_secs = 0.0
    first_retcode = "FAILED"
    solve_started = time_ns()
    try
        skip_reason = solver_skip_reason(solver_name, nlp.meta)
        if skip_reason !== nothing
            println("CUTEst solve skipped for $problem_name with $solver_name: $skip_reason")

            return (;
                problem = problem_name,
                solver = solver_name,
                solver_class = class,
                n_vars = nlp.meta.nvar,
                secs = 0.0,
                first_solve_secs = 0.0,
                compile_secs = 0.0,
                decode_secs = decode_secs,
                solver_reported_secs = NaN,
                retcode = "SKIPPED",
                first_retcode = "SKIPPED",
                status = "SKIPPED",
                NO_QUALITY...,
            )
        end

        # OptimizationBase warns on every Newton solve that no `SecondOrder` ADtype was
        # given, even though the NLPModels Hessian is supplied directly.
        logger = class == "Newton" ? ConsoleLogger(stderr, Logging.Error) : current_logger()
        prob = build_problem(nlp, solver_name)
        Base.cumulative_compile_timing(true)
        compile_before = Base.cumulative_compile_time_ns()
        solve_started = time_ns()
        sol = try
            with_logger(logger) do
                solve_once(prob, solver_name; maxtime)
            end
        finally
            first_solve_secs = elapsed_seconds(solve_started)
            compile_secs = compile_seconds_since(compile_before)
            Base.cumulative_compile_timing(false)
        end
        first_retcode = retcode_name(sol.retcode)
        secs = first_solve_secs

        if rerun_worthwhile(first_retcode, first_solve_secs; maxtime)
            prob = build_problem(nlp, solver_name)
            solve_started = time_ns()
            sol = with_logger(logger) do
                solve_once(prob, solver_name; maxtime)
            end
            secs = elapsed_seconds(solve_started)
        else
            print(" [rerun skipped: first run ended with $first_retcode]")
        end

        quality = try
            solution_quality(nlp, sol.u)
        catch
            NO_QUALITY
        end

        return (;
            problem = problem_name,
            solver = solver_name,
            solver_class = class,
            n_vars = nlp.meta.nvar,
            secs = secs,
            first_solve_secs = first_solve_secs,
            compile_secs = compile_secs,
            decode_secs = decode_secs,
            solver_reported_secs = solve_seconds(sol, NaN),
            retcode = retcode_name(sol.retcode),
            first_retcode = first_retcode,
            status = "OK",
            quality...,
        )
    catch err
        bt = catch_backtrace()
        print_exception("CUTEst solve failed for $problem_name with $solver_name", err, bt)

        return (;
            problem = problem_name,
            solver = solver_name,
            solver_class = class,
            n_vars = nlp.meta.nvar,
            secs = elapsed_seconds(solve_started),
            first_solve_secs = first_solve_secs,
            compile_secs = compile_secs,
            decode_secs = decode_secs,
            solver_reported_secs = NaN,
            retcode = "FAILED",
            first_retcode = first_retcode,
            status = "FAILED",
            NO_QUALITY...,
        )
    finally
        if nlp !== nothing
            try
                finalize(nlp)
            catch err
                @warn "Unable to finalize CUTEst problem" problem = problem_name exception = (
                    err, catch_backtrace(),
                )
            end
        end
    end
end

warmup_problem_for(solver_name) =
    solver_name in CONSTRAINED_SOLVERS ? WARMUP_CONSTRAINED_PROBLEM :
    WARMUP_UNCONSTRAINED_PROBLEM

# Run one throwaway solve per solver on a tiny problem so that JIT compilation of the
# solver / NLPModels / MOI code paths is not attributed to whichever problem comes first
# in the benchmark loop. The per-solver compile time measured here is printed and
# returned; the per-pair `compile_secs` column then only shows recompilation caused by a
# change of problem shape. `problem` may be a problem name applied to every solver, or
# `nothing` to pick a tiny unconstrained/constrained problem per solver via
# `warmup_problem_for`.
function warmup_solvers(solvers; problem = nothing, maxtime = SOLVE_TIMEOUT_SECONDS)
    rows = NamedTuple[]

    println()
    println("Warming up solvers (JIT), results are not benchmarked:")
    for solver_name in solvers
        problem_name = problem === nothing ? warmup_problem_for(solver_name) : problem
        @printf("  %-18s %-24s", solver_name, problem_name)
        started = time_ns()
        row = run_single_solve(problem_name, solver_name; maxtime)
        push!(rows, row)
        @printf(
            " %s %s first %.3fs compile %.3fs clean %.3fs (total %.3fs)\n", row.status,
            row.retcode, row.first_solve_secs, row.compile_secs, row.secs,
            elapsed_seconds(started)
        )
    end
    return DataFrame(rows)
end

function run_benchmarks(
        category, problems, solvers;
        warmup = true, maxtime = SOLVE_TIMEOUT_SECONDS,
    )
    rows = NamedTuple[]

    warmup && warmup_solvers(solvers; maxtime)

    println()
    println("Running $category benchmarks")
    println("Problems: ", length(problems))
    println("Solvers: ", join(solvers, ", "))
    println("Per-solve time cap: ", maxtime, " s")

    for problem_name in problems
        for solver_name in solvers
            @printf("  %-18s %-24s", solver_name, problem_name)
            row = run_single_solve(problem_name, solver_name; maxtime)
            push!(rows, merge((category = category,), row))
            @printf(
                " %s %s %.3fs (first %.3fs, compile %.3fs, decode %.3fs)\n", row.status,
                row.retcode, row.secs, row.first_solve_secs, row.compile_secs, row.decode_secs
            )
        end
    end

    results = isempty(rows) ? DataFrame(
            category = String[], problem = String[], solver = String[],
            solver_class = String[], n_vars = Int[], secs = Float64[],
            first_solve_secs = Float64[], compile_secs = Float64[], decode_secs = Float64[],
            solver_reported_secs = Float64[], retcode = String[], first_retcode = String[],
            status = String[],
            objective = Float64[], grad_norm = Float64[], cons_viol = Float64[]
        ) : DataFrame(rows)

    add_quality_columns!(results)
    assert_has_measurements(results, category)
    return results
end

function count_distribution(values)
    counts = sort(collect(countmap(values)); by = x -> x[2], rev = true)
    return join(["$code=$n" for (code, n) in counts], ", ")
end

# Smallest objective among runs whose point is feasible to FEAS_TOL; NaN if none is.
function best_feasible_objective(objective, cons_viol)
    best = Inf
    for (f, v) in zip(objective, cons_viol)
        if isfinite(f) && isfinite(v) && v <= FEAS_TOL && f < best
            best = f
        end
    end
    return isfinite(best) ? best : NaN
end

function is_verified_success(retcode, cons_viol, grad_norm, obj_gap)
    retcode in SUCCESS_RETCODES || return false
    isfinite(cons_viol) && cons_viol <= FEAS_TOL || return false
    return (isfinite(grad_norm) && grad_norm <= OPT_TOL) ||
        (isfinite(obj_gap) && obj_gap <= GAP_TOL)
end

# Add the cross-solver quality columns: `best_objective` (per category/problem),
# the relative gap `obj_gap = (objective - best_objective) / max(1, |best_objective|)`
# and the boolean `verified_success` criterion. Idempotent.
function add_quality_columns!(results)
    if nrow(results) == 0
        results.best_objective = Float64[]
        results.obj_gap = Float64[]
        results.verified_success = Bool[]
        return results
    end

    transform!(
        groupby(results, [:category, :problem]),
        [:objective, :cons_viol] => best_feasible_objective => :best_objective,
    )
    results.obj_gap = (results.objective .- results.best_objective) ./
        max.(1.0, abs.(results.best_objective))
    results.verified_success = is_verified_success.(
        results.retcode, results.cons_viol, results.grad_norm, results.obj_gap
    )

    return results
end

# Fails the page when a category is degenerate:
#   * fewer than `MIN_COMPLETED_FRACTION` of the *attempted* (non-SKIPPED) rows completed
#     at the harness level (`status == "OK"`), or
#   * any attempted solver has zero completed rows (every non-SKIPPED call errored /
#     was killed -> broken backend).
# `SKIPPED` rows are intentional inapplicability and are excluded from both checks.
# Zero *success* (retcode-based) for a solver is only warned about, since global
# heuristics such as SimulatedAnnealing legitimately report 0% success on some categories.
function assert_has_measurements(results, category)
    nrow(results) == 0 && error("CUTEst benchmark for $category produced no rows at all.")

    attempted = filter(:status => !=("SKIPPED"), results)
    nrow(attempted) == 0 && error(
        "CUTEst benchmark for $category produced only SKIPPED rows (no solver was applicable)."
    )

    completed = count(==("OK"), attempted.status)
    completed_fraction = completed / nrow(attempted)
    status_distribution = count_distribution(results.status)
    retcode_distribution = count_distribution(results.retcode)

    failures = String[]
    if completed_fraction < MIN_COMPLETED_FRACTION
        msg = @sprintf(
            "only %d/%d attempted rows (%.1f%%) completed, below MIN_COMPLETED_FRACTION = %.1f%%",
            completed, nrow(attempted), 100 * completed_fraction, 100 * MIN_COMPLETED_FRACTION
        )
        push!(failures, msg)
    end

    solvers = unique(attempted.solver)
    broken_solvers = filter(solvers) do solver
        count(row -> row.solver == solver && row.status == "OK", eachrow(attempted)) == 0
    end
    if !isempty(broken_solvers)
        push!(failures, "solver(s) with zero completed rows: " * join(broken_solvers, ", "))
    end

    if !isempty(failures)
        error(
            "CUTEst benchmark for $category is degenerate: " * join(failures, "; ") *
                ". Status distribution: $status_distribution. " *
                "Return code distribution: $retcode_distribution",
        )
    end

    no_success_solvers = filter(solvers) do solver
        count(row -> row.solver == solver && row.retcode in SUCCESS_RETCODES, eachrow(attempted)) == 0
    end
    if !isempty(no_success_solvers)
        @warn "CUTEst benchmark for $category has solver(s) with 0% success (no retcode in SUCCESS_RETCODES); allowed, but worth a look" category solvers = join(no_success_solvers, ", ") retcode_distribution
    end

    skipped = count(==("SKIPPED"), results.status)
    completed_pct = round(100 * completed_fraction; digits = 1)
    println(
        "  $category: $completed/$(nrow(attempted)) attempted rows completed ($completed_pct%)",
        skipped > 0 ? "; $skipped skipped" : "",
        "; status: $status_distribution",
    )
    return nothing
end

# Skipped pairings never ran a solver, so they carry no timing information.
function median_attempted_secs(status, secs)
    attempted = secs[status .!= "SKIPPED"]
    return isempty(attempted) ? NaN : median(attempted)
end

function summarize_results(results)
    hasproperty(results, :verified_success) || add_quality_columns!(results)

    println()
    println("Return code distribution:")
    for (code, n) in sort(collect(countmap(results.retcode)); by = x -> x[2], rev = true)
        println("  $code: $n")
    end

    summary = combine(
        groupby(results, [:category, :solver_class, :solver]),
        :status => (x -> count(==("OK"), x)) => :completed_runs,
        :retcode => (x -> count(in(SUCCESS_RETCODES), x)) => :successful_runs,
        :verified_success => count => :verified_runs,
        :status => (x -> count(==("SKIPPED"), x)) => :skipped_runs,
        :retcode => length => :total_runs,
        [:status, :secs] => median_attempted_secs => :median_secs,
        :compile_secs => median => :median_compile_secs,
        :compile_secs => sum => :total_compile_secs,
    )

    summary.completion_rate = round.(summary.completed_runs ./ summary.total_runs .* 100; digits = 1)
    summary.success_rate = round.(summary.successful_runs ./ summary.total_runs .* 100; digits = 1)
    summary.verified_success_rate = round.(
        summary.verified_runs ./ summary.total_runs .* 100; digits = 1
    )
    summary.convergence_based = convergence_based.(summary.solver_class)
    sort!(summary, [:category, order(:convergence_based; rev = true), :solver_class, :solver])

    local_summary = filter(:convergence_based => identity, summary)
    budget_summary = filter(:convergence_based => !, summary)

    println()
    println(
        "Summary (local convergence-based solvers; success = return code in ",
        join(sort(collect(SUCCESS_RETCODES)), "/"), "):",
    )
    display(select(local_summary, Not(:convergence_based)))

    if nrow(budget_summary) > 0
        println()
        println(
            "Summary (budget-limited solvers; these run until maxiters and do not report ",
            "convergence through the return code, so only completion and time are shown):",
        )
        display(
            select(
                budget_summary, Not(
                    [
                        :successful_runs, :success_rate, :verified_runs, :verified_success_rate,
                        :convergence_based,
                    ]
                )
            )
        )
    end

    return summary
end

function plot_solve_times(results, title)
    completed = filter(:status => ==("OK"), results)
    if nrow(completed) == 0
        return nothing
    end

    solve_time_plot = @df completed scatter(
        :n_vars,
        :secs,
        group = :solver,
        xlabel = "Number of variables",
        ylabel = "Seconds",
        title = title,
        yscale = :log10,
        legend = :topleft,
        size = (900, 600),
    )

    return display(solve_time_plot)
end

# Compile time measured during the first solve of each pair. Rows with exactly zero
# compile time (nothing was compiled) cannot be drawn on a log axis and are dropped.
function plot_compile_times(results, title)
    compiled = filter([:status, :compile_secs] => (s, c) -> s == "OK" && c > 0, results)
    if nrow(compiled) == 0
        println("No rows with nonzero compile time; skipping compile time plot.")
        return nothing
    end

    compile_time_plot = @df compiled scatter(
        :n_vars,
        :compile_secs,
        group = :solver,
        xlabel = "Number of variables",
        ylabel = "Compile seconds (first solve)",
        title = title,
        yscale = :log10,
        legend = :topleft,
        size = (900, 600),
    )

    return display(compile_time_plot)
end

function plot_success_rates(summary, title)
    local_summary = filter(:convergence_based => identity, summary)
    if nrow(local_summary) == 0
        return nothing
    end

    # Alphabetical group order keeps solvers of the same class adjacent in the legend.
    local_summary.label = local_summary.solver_class .* ": " .* local_summary.solver

    success_rate_plot = @df local_summary groupedbar(
        :category,
        :success_rate,
        group = :label,
        xlabel = "Problem category",
        ylabel = "Success rate (%)",
        title = title,
        xrotation = 30,
        legend = :topright,
        size = (900, 600),
    )

    return display(success_rate_plot)
end

# Dolan–Moré performance ratios for one category: `ratios[solver][i]` is the solve time of
# `solver` on the i-th problem divided by the fastest verified solve of that problem.
# Runs that are not verified successes (and problems nobody solved) get ratio Inf.
function performance_ratios(results)
    problems = unique(results.problem)
    solvers = unique(results.solver)
    ratios = Dict(s => fill(Inf, length(problems)) for s in solvers)

    for (i, problem) in enumerate(problems)
        rows = filter(r -> r.problem == problem && r.verified_success, results)
        nrow(rows) == 0 && continue
        best = max(minimum(rows.secs), 1.0e-9)
        for r in eachrow(rows)
            ratios[r.solver][i] = max(r.secs, 1.0e-9) / best
        end
    end

    return problems, solvers, ratios
end

function performance_profile_subplot(results, title)
    problems, solvers, ratios = performance_ratios(results)
    nprob = length(problems)

    finite_log_ratios = [log2(r) for v in values(ratios) for r in v if isfinite(r)]
    tau_max = isempty(finite_log_ratios) ? 1.0 : max(1.0, maximum(finite_log_ratios) * 1.05)
    taus = sort!(unique!(vcat(0.0, finite_log_ratios, tau_max)))

    subplot = plot(
        xlabel = "log2(time ratio to best solver)",
        ylabel = "Fraction of problems solved",
        title = title,
        ylims = (0, 1.02),
        legend = :bottomright,
    )
    for solver in solvers
        log_ratios = log2.(ratios[solver])
        rho = [count(<=(tau), log_ratios) / nprob for tau in taus]
        plot!(subplot, taus, rho; seriestype = :steppost, label = solver, linewidth = 2)
    end

    return subplot
end

# Performance profile of solve time (Dolan & Moré, 2002), one panel per category.
function plot_performance_profile(results, title)
    hasproperty(results, :verified_success) || add_quality_columns!(results)
    nrow(results) == 0 && return nothing

    categories = unique(results.category)
    subplots = [
        performance_profile_subplot(
                filter(:category => ==(category), results),
                length(categories) == 1 ? title : "$title: $category",
            )
            for category in categories
    ]

    profile_plot = plot(
        subplots...; layout = (length(subplots), 1),
        size = (900, 500 * length(subplots))
    )

    return display(profile_plot)
end

# Work-precision scatter: solve time against the attained precision, where precision is
# `metric` (`:obj_gap` by default, `:grad_norm` is the natural choice for unconstrained
# problems). Only completed runs with a finite metric are shown; for `:obj_gap` the point
# must also be feasible to FEAS_TOL so that infeasible points cannot show a spurious gap.
# Values at or below `floor` are drawn at `floor` so they fit on the log axis.
function plot_work_precision(results, title; metric = :obj_gap, floor = 1.0e-16)
    hasproperty(results, :verified_success) || add_quality_columns!(results)

    completed = filter(results) do r
        r.status == "OK" && isfinite(r[metric]) && r.secs > 0 &&
            (metric != :obj_gap || (isfinite(r.cons_viol) && r.cons_viol <= FEAS_TOL))
    end
    nrow(completed) == 0 && return nothing

    precision = max.(completed[!, metric], floor)

    work_precision_plot = scatter(
        precision,
        completed.secs,
        group = completed.solver,
        xlabel = "$(metric) (floored at $(floor))",
        ylabel = "Seconds",
        title = title,
        xscale = :log10,
        yscale = :log10,
        legend = :topright,
        size = (900, 600),
    )

    return display(work_precision_plot)
end
