using Printf
using DataFrames
using StatsBase: countmap
using Statistics

const MAX_PROBLEMS_PER_CATEGORY = 20
const MAX_NVAR = 1_000
const MAX_NCON = 1_000
const SOLVE_TIMEOUT_SECONDS = 90

const SUCCESS_RETCODES = Set(["Success", "FirstOrderOptimal", "MaxIters", "MaxTime"])

const KNOWN_BAD_PROBLEMS = Set(
    lowercase.(
        String[
            "BLOWEYA", "CHARDIS1", "CLEUVEN4", "CMPC3", "CMPC10", "CVXQP2",
            "DITTERT", "HIER13", "LUKVLE8", "LUKVLI7", "MPC2", "PATTERNNE",
            "READING2", "READING6", "NINENEW", "MSS1",
        ]
    )
)

const UNCONSTRAINED_SOLVERS = ["LBFGS", "ConjugateGradient", "NelderMead"]
const CONSTRAINED_SOLVERS = ["Ipopt"]

function problem_metadata(name)
    nlp = nothing
    try
        nlp = CUTEstModel(name)
        return (; ok = true, nvar = nlp.meta.nvar, ncon = nlp.meta.ncon)
    catch
        return (; ok = false, nvar = -1, ncon = -1)
    finally
        if nlp !== nothing
            try
                finalize(nlp)
            catch
            end
        end
    end
end

function select_safe_problems(candidates; max_problems = MAX_PROBLEMS_PER_CATEGORY)
    selected = String[]

    for name in candidates

        lowercase(name) in KNOWN_BAD_PROBLEMS && continue

        meta = problem_metadata(name)


        meta.ok || continue
        meta.nvar <= MAX_NVAR || continue
        meta.ncon <= MAX_NCON || continue

        push!(selected, name)

        length(selected) >= max_problems && break
    end

    return selected
end

function benchmark_dir()
    return isdefined(Main, :WEAVE_ARGS) ? WEAVE_ARGS[:folder] : @__DIR__
end

function parse_child_result(output, problem_name, solver_name, secs, status)
    result_lines = filter(line -> startswith(line, "CUTEST_RESULT\t"), split(output, '\n'))

    if isempty(result_lines)
        return (;
            problem = problem_name,
            solver = solver_name,
            n_vars = problem_metadata(problem_name).nvar,
            secs = secs,
            retcode = status,
            status = status,
        )
    end

    fields = split(last(result_lines), '\t')
    return (;
        problem = String(fields[2]),
        solver = String(fields[3]),
        n_vars = parse(Int, fields[4]),
        secs = parse(Float64, fields[5]),
        retcode = String(fields[6]),
        status = String(fields[7]),
    )
end

function run_child_solve(problem_name, solver_name)
    script = joinpath(benchmark_dir(), "solve_cutest_case.jl")
    project_dir = dirname(Base.active_project())
    cmd = `$(Base.julia_cmd()) --project=$project_dir $script $problem_name $solver_name`

    out = Pipe()
    err = Pipe()
    started = time()
    proc = run(pipeline(cmd; stdout = out, stderr = err); wait = false)

    close(out.in)
    close(err.in)

    out_task = @async read(out, String)
    err_task = @async read(err, String)

    while !process_exited(proc)
        if time() - started > SOLVE_TIMEOUT_SECONDS
            kill(proc)
            wait(proc)
            output = fetch(out_task)
            fetch(err_task)

            return parse_child_result(
                output,
                problem_name,
                solver_name,
                SOLVE_TIMEOUT_SECONDS,
                "TIMEOUT",
            )
        end
        sleep(0.25)
    end

    wait(proc)
    output = fetch(out_task)
    fetch(err_task)

    return parse_child_result(output, problem_name, solver_name, time() - started, "CRASHED")
end

function run_benchmarks(category, problems, solvers)
    rows = NamedTuple[]

    println()
    println("Running $category benchmarks")
    println("Problems: ", length(problems))
    println("Solvers: ", join(solvers, ", "))

    for problem_name in problems
        for solver_name in solvers
            @printf("  %-12s %-24s", solver_name, problem_name)
            row = run_child_solve(problem_name, solver_name)
            push!(rows, merge((category = category,), row))
            @printf(" %s %s %.3fs\n", row.status, row.retcode, row.secs)
        end
    end

    if isempty(rows)
        return DataFrame(
            category = String[], problem = String[], solver = String[],
            n_vars = Int[], secs = Float64[], retcode = String[], status = String[]
        )
    end

    return DataFrame(rows)
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
        :secs => median => :median_secs
    )

    summary.completion_rate = round.(summary.completed_runs ./ summary.total_runs .* 100; digits = 1)
    summary.success_rate = round.(summary.successful_runs ./ summary.total_runs .* 100; digits = 1)

    println()
    println("Summary:")
    display(summary)

    return summary
end

function plot_solve_times(results, title)
    if nrow(results) == 0
        return nothing
    end

    solve_time_plot = @df results scatter(
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
