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

function run_single_solve(problem_name, solver_name)
    nlp = nothing
    started = time()

    try
        nlp = CUTEstModel(problem_name)
    catch err
        bt = catch_backtrace()
        print_exception("CUTEst problem load failed for $problem_name", err, bt)

        return (;
            problem = problem_name,
            solver = solver_name,
            n_vars = -1,
            secs = time() - started,
            retcode = "LOAD_FAILED",
            status = "LOAD_FAILED",
        )
    end

    try
        prob = OptimizationNLPModels.OptimizationProblem(nlp)
        sol = solve(prob, optimizer_from_name(solver_name);
            maxiters = SOLVE_MAXITERS,
            maxtime = SOLVE_TIMEOUT_SECONDS)
        secs = solve_seconds(sol, time() - started)

        return (;
            problem = problem_name,
            solver = solver_name,
            n_vars = nlp.meta.nvar,
            secs = secs,
            retcode = retcode_name(sol.retcode),
            status = "OK",
        )
    catch err
        bt = catch_backtrace()
        print_exception("CUTEst solve failed for $problem_name with $solver_name", err, bt)

        return (;
            problem = problem_name,
            solver = solver_name,
            n_vars = nlp.meta.nvar,
            secs = time() - started,
            retcode = "FAILED",
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

function run_benchmarks(category, problems, solvers)
    rows = NamedTuple[]

    println()
    println("Running $category benchmarks")
    println("Problems: ", length(problems))
    println("Solvers: ", join(solvers, ", "))

    for problem_name in problems
        for solver_name in solvers
            @printf("  %-18s %-24s", solver_name, problem_name)
            row = run_single_solve(problem_name, solver_name)
            push!(rows, merge((category = category,), row))
            @printf(" %s %s %.3fs\n", row.status, row.retcode, row.secs)
        end
    end

    results = isempty(rows) ? DataFrame(
        category = String[], problem = String[], solver = String[],
        n_vars = Int[], secs = Float64[], retcode = String[], status = String[]
    ) : DataFrame(rows)

    assert_has_measurements(results, category)
    return results
end

function assert_has_measurements(results, category)
    completed = nrow(results) == 0 ? 0 : count(==("OK"), results.status)
    completed > 0 && return nothing

    distribution = nrow(results) == 0 ? "no rows" :
                   join(["$code=$n" for (code, n) in sort(
                       collect(countmap(results.retcode)); by = x -> x[2], rev = true)], ", ")

    error(
        "CUTEst benchmark for $category produced no completed solver measurements. " *
        "Return code distribution: $distribution",
    )
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
