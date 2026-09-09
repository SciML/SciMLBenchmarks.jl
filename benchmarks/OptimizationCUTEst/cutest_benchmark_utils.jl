using Printf
using DataFrames
using Plots
using StatsPlots
using StatsBase: countmap
using Statistics

using Optimization
using OptimizationNLPModels
using OptimizationOptimJL
using OptimizationOptimJL: LBFGS, ConjugateGradient, NelderMead, SimulatedAnnealing,
                           ParticleSwarm
using OptimizationMOI
using OptimizationMOI: MOI
using Ipopt

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
    "ConjugateGradient",
    "NelderMead",
    "SimulatedAnnealing",
    "ParticleSwarm",
]
const CONSTRAINED_SOLVERS = ["Ipopt"]

function optimizer_from_name(name)
    if name == "LBFGS"
        return LBFGS()
    elseif name == "ConjugateGradient"
        return ConjugateGradient()
    elseif name == "NelderMead"
        return NelderMead()
    elseif name == "SimulatedAnnealing"
        return SimulatedAnnealing()
    elseif name == "ParticleSwarm"
        return ParticleSwarm()
    elseif name == "Ipopt"
        return MOI.OptimizerWithAttributes(
            Ipopt.Optimizer,
            "max_iter" => SOLVE_MAXITERS,
            "max_wall_time" => SOLVE_TIMEOUT_SECONDS,
            "hessian_approximation" => "limited-memory",
            "tol" => 1.0e-6,
            "print_level" => 0,
        )
    else
        error("Unknown optimizer: $name")
    end
end

function problem_metadata(name)
    nlp = nothing
    try
        nlp = CUTEstModel(name)
        return (; ok = true, nvar = nlp.meta.nvar, ncon = nlp.meta.ncon)
    catch err
        @warn "Unable to load CUTEst problem metadata" problem = name exception = (
            err, catch_backtrace())
        return (; ok = false, nvar = -1, ncon = -1)
    finally
        if nlp !== nothing
            try
                finalize(nlp)
            catch err
                @warn "Unable to finalize CUTEst problem metadata" problem = name exception = (
                    err, catch_backtrace())
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
    println(prefix, ": ", sprint(showerror, err, bt))
end

elapsed_seconds(started_ns) = (time_ns() - started_ns) / 1.0e9

# `Base.cumulative_compile_time_ns()` returns `(compile_ns, recompile_ns)` on Julia >= 1.9
# (a bare integer on older versions); recompilation is a subset of compilation, so the
# first element is the total compile time (this is what `@time` reports).
compile_ns(t) = t isa Tuple ? first(t) : t
compile_seconds_since(before) = (compile_ns(Base.cumulative_compile_time_ns()) - compile_ns(before)) / 1.0e9

solve_once(prob, solver_name) = solve(
    prob, optimizer_from_name(solver_name);
    maxiters = SOLVE_MAXITERS,
    maxtime = SOLVE_TIMEOUT_SECONDS
)

# The second (clean) solve is skipped when the first one already hit the time cap or
# errored, so a timed-out pair costs one cap, not two.
rerun_worthwhile(first_retcode, first_solve_secs) =
    first_retcode != "MaxTime" && first_solve_secs < SOLVE_TIMEOUT_SECONDS

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
function run_single_solve(problem_name, solver_name)
    nlp = nothing
    decode_started = time_ns()

    try
        nlp = CUTEstModel(problem_name)
    catch err
        bt = catch_backtrace()
        print_exception("CUTEst problem load failed for $problem_name", err, bt)

        return (;
            problem = problem_name,
            solver = solver_name,
            n_vars = -1,
            secs = 0.0,
            first_solve_secs = 0.0,
            compile_secs = 0.0,
            decode_secs = elapsed_seconds(decode_started),
            solver_reported_secs = NaN,
            retcode = "LOAD_FAILED",
            first_retcode = "LOAD_FAILED",
            status = "LOAD_FAILED",
        )
    end
    decode_secs = elapsed_seconds(decode_started)

    first_solve_secs = 0.0
    compile_secs = 0.0
    first_retcode = "FAILED"
    solve_started = time_ns()
    try
        prob = OptimizationNLPModels.OptimizationProblem(nlp)
        Base.cumulative_compile_timing(true)
        compile_before = Base.cumulative_compile_time_ns()
        solve_started = time_ns()
        sol = try
            solve_once(prob, solver_name)
        finally
            first_solve_secs = elapsed_seconds(solve_started)
            compile_secs = compile_seconds_since(compile_before)
            Base.cumulative_compile_timing(false)
        end
        first_retcode = retcode_name(sol.retcode)
        secs = first_solve_secs

        if rerun_worthwhile(first_retcode, first_solve_secs)
            prob = OptimizationNLPModels.OptimizationProblem(nlp)
            solve_started = time_ns()
            sol = solve_once(prob, solver_name)
            secs = elapsed_seconds(solve_started)
        else
            print(" [rerun skipped: first run ended with $first_retcode]")
        end

        return (;
            problem = problem_name,
            solver = solver_name,
            n_vars = nlp.meta.nvar,
            secs = secs,
            first_solve_secs = first_solve_secs,
            compile_secs = compile_secs,
            decode_secs = decode_secs,
            solver_reported_secs = solve_seconds(sol, NaN),
            retcode = retcode_name(sol.retcode),
            first_retcode = first_retcode,
            status = "OK",
        )
    catch err
        bt = catch_backtrace()
        print_exception("CUTEst solve failed for $problem_name with $solver_name", err, bt)

        return (;
            problem = problem_name,
            solver = solver_name,
            n_vars = nlp.meta.nvar,
            secs = elapsed_seconds(solve_started),
            first_solve_secs = first_solve_secs,
            compile_secs = compile_secs,
            decode_secs = decode_secs,
            solver_reported_secs = NaN,
            retcode = "FAILED",
            first_retcode = first_retcode,
            status = "FAILED",
        )
    finally
        if nlp !== nothing
            try
                finalize(nlp)
            catch err
                @warn "Unable to finalize CUTEst problem" problem = problem_name exception = (
                    err, catch_backtrace())
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
function warmup_solvers(solvers; problem = nothing)
    rows = NamedTuple[]

    println()
    println("Warming up solvers (JIT), results are not benchmarked:")
    for solver_name in solvers
        problem_name = problem === nothing ? warmup_problem_for(solver_name) : problem
        @printf("  %-18s %-24s", solver_name, problem_name)
        started = time_ns()
        row = run_single_solve(problem_name, solver_name)
        push!(rows, row)
        @printf(
            " %s %s first %.3fs compile %.3fs clean %.3fs (total %.3fs)\n", row.status,
            row.retcode, row.first_solve_secs, row.compile_secs, row.secs,
            elapsed_seconds(started)
        )
    end
    return DataFrame(rows)
end

function run_benchmarks(category, problems, solvers; warmup = true)
    rows = NamedTuple[]

    warmup && warmup_solvers(solvers)

    println()
    println("Running $category benchmarks")
    println("Problems: ", length(problems))
    println("Solvers: ", join(solvers, ", "))

    for problem_name in problems
        for solver_name in solvers
            @printf("  %-18s %-24s", solver_name, problem_name)
            row = run_single_solve(problem_name, solver_name)
            push!(rows, merge((category = category,), row))
            @printf(
                " %s %s %.3fs (first %.3fs, compile %.3fs, decode %.3fs)\n", row.status,
                row.retcode, row.secs, row.first_solve_secs, row.compile_secs, row.decode_secs
            )
        end
    end

    results = isempty(rows) ? DataFrame(
            category = String[], problem = String[], solver = String[],
            n_vars = Int[], secs = Float64[], first_solve_secs = Float64[],
            compile_secs = Float64[], decode_secs = Float64[],
            solver_reported_secs = Float64[], retcode = String[], first_retcode = String[],
            status = String[]
        ) : DataFrame(rows)

    assert_has_measurements(results, category)
    return results
end

function count_distribution(values)
    counts = sort(collect(countmap(values)); by = x -> x[2], rev = true)
    return join(["$code=$n" for (code, n) in counts], ", ")
end

# Fails the page when a category is degenerate:
#   * fewer than `MIN_COMPLETED_FRACTION` of the rows completed at the harness level
#     (`status == "OK"`), or
#   * any solver has zero completed rows (every call errored / was killed -> broken backend).
# Zero *success* (retcode-based) for a solver is only warned about, since global
# heuristics such as SimulatedAnnealing legitimately report 0% success on some categories.
function assert_has_measurements(results, category)
    nrow(results) == 0 && error("CUTEst benchmark for $category produced no rows at all.")

    completed = count(==("OK"), results.status)
    completed_fraction = completed / nrow(results)
    status_distribution = count_distribution(results.status)
    retcode_distribution = count_distribution(results.retcode)

    failures = String[]
    if completed_fraction < MIN_COMPLETED_FRACTION
        msg = @sprintf(
            "only %d/%d rows (%.1f%%) completed, below MIN_COMPLETED_FRACTION = %.1f%%",
            completed, nrow(results), 100 * completed_fraction, 100 * MIN_COMPLETED_FRACTION
        )
        push!(failures, msg)
    end

    solvers = unique(results.solver)
    broken_solvers = filter(solvers) do solver
        count(row -> row.solver == solver && row.status == "OK", eachrow(results)) == 0
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
        count(row -> row.solver == solver && row.retcode in SUCCESS_RETCODES, eachrow(results)) == 0
    end
    if !isempty(no_success_solvers)
        @warn "CUTEst benchmark for $category has solver(s) with 0% success (no retcode in SUCCESS_RETCODES); allowed, but worth a look" category solvers = join(no_success_solvers, ", ") retcode_distribution
    end

    completed_pct = round(100 * completed_fraction; digits = 1)
    println("  $category: $completed/$(nrow(results)) rows completed ($completed_pct%); status: $status_distribution")
    return nothing
end

function summarize_results(results)
    println()
    println("Return code distribution:")
    for (code, n) in sort(collect(countmap(results.retcode)); by = x -> x[2], rev = true)
        println("  $code: $n")
    end

    summary = combine(
        groupby(results, [:category, :solver]),
        :status => (x -> count(==("OK"), x)) => :completed_runs,
        :retcode => (x -> count(in(SUCCESS_RETCODES), x)) => :successful_runs,
        :retcode => length => :total_runs,
        :secs => median => :median_secs,
        :compile_secs => median => :median_compile_secs,
        :compile_secs => sum => :total_compile_secs,
    )

    summary.completion_rate = round.(summary.completed_runs ./ summary.total_runs .* 100; digits = 1)
    summary.success_rate = round.(summary.successful_runs ./ summary.total_runs .* 100; digits = 1)

    println()
    println("Summary:")
    display(summary)

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
    if nrow(summary) == 0
        return nothing
    end

    success_rate_plot = @df summary groupedbar(
        :category,
        :success_rate,
        group = :solver,
        xlabel = "Problem category",
        ylabel = "Success rate (%)",
        title = title,
        xrotation = 30,
        legend = :topright,
        size = (900, 600),
    )

    return display(success_rate_plot)
end
