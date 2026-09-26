using Distributed
using LinearAlgebra
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
# Extra wall-clock budget beyond the per-pass `maxtime` caps. The hard deadline
# is `2 * maxtime + HARD_TIMEOUT_GRACE_SECONDS` so a legitimate two-pass solve
# (first + clean) that stays within each cooperative `maxtime` is not killed,
# while SIF decode / wrap-up / cleanup still have a fixed allowance.
const HARD_TIMEOUT_GRACE_SECONDS = 60.0

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

function problem_metadata_inprocess(name)
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

# Decode CUTEst metadata. On the page master with `hard_timeout = true` (default),
# the decode runs on the long-lived benchmark worker so a hung SIF decode cannot
# stall selection. Workers always decode in-process to avoid recursive workers.
function problem_metadata(
        name; hard_timeout = true, maxtime = SOLVE_TIMEOUT_SECONDS
    )
    if myid() != 1 || !hard_timeout
        return problem_metadata_inprocess(name)
    end

    worker = ensure_benchmark_worker()
    # Metadata is a single decode (no two-pass solve); one `maxtime` + grace.
    deadline = maxtime + HARD_TIMEOUT_GRACE_SECONDS
    waiter = @async remotecall_fetch(
        Core.eval, worker, Main, :(problem_metadata_inprocess($name))
    )

    if timedwait(() -> istaskdone(waiter), deadline; pollint = 0.5) === :timed_out
        println(
            " HARD TIMEOUT: metadata decode for $name exceeded $(deadline)s, " *
                "killing worker $worker"
        )
        kill_benchmark_worker(worker)
        return (; ok = false, nvar = -1, ncon = -1)
    end

    try
        return fetch(waiter)
    catch err
        bt = catch_backtrace()
        print_exception("Benchmark worker $worker failed loading metadata for $name", err, bt)
        kill_benchmark_worker(worker)
        return (; ok = false, nvar = -1, ncon = -1)
    end
end

function select_safe_problems(
        candidates;
        max_problems = MAX_PROBLEMS_PER_CATEGORY,
        max_var = MAX_NVAR,
        max_con = MAX_NCON,
        hard_timeout = true,
    )
    selected = String[]

    for name in candidates
        lowercase(name) in KNOWN_BAD_PROBLEMS && continue

        meta = problem_metadata(name; hard_timeout)

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
#
# With `hard_timeout = true`, warm-up uses the guarded worker path so hangs cannot stall
# the page and JIT lands in the process that will be measured.
function warmup_solvers(
        solvers; problem = nothing, maxtime = SOLVE_TIMEOUT_SECONDS, hard_timeout = false
    )
    rows = NamedTuple[]

    println()
    println("Warming up solvers (JIT), results are not benchmarked:")
    for solver_name in solvers
        sync_warmup_worker!()
        if hard_timeout && solver_name in WARMUP_FAILED_SOLVERS
            @printf("  %-18s %-24s skipped (prior warm-up timeout)\n", solver_name, "-")
            continue
        end
        if hard_timeout && solver_warmed_here(solver_name)
            @printf(
                "  %-18s %-24s already warm on worker %d\n", solver_name, "-",
                BENCHMARK_WORKER[],
            )
            continue
        end

        problem_name = problem === nothing ? warmup_problem_for(solver_name) : problem
        @printf("  %-18s %-24s", solver_name, problem_name)
        started = time_ns()
        row = if hard_timeout
            run_single_solve_guarded(problem_name, solver_name; maxtime)
        else
            run_single_solve(problem_name, solver_name; maxtime)
        end
        push!(rows, row)
        @printf(
            " %s %s first %.3fs compile %.3fs clean %.3fs (total %.3fs)\n", row.status,
            row.retcode, row.first_solve_secs, row.compile_secs, row.secs,
            elapsed_seconds(started)
        )

        if hard_timeout
            sync_warmup_worker!()
            if row.status == "TIMEOUT" || row.retcode in ("HARD_TIMEOUT", "WORKER_FAILED")
                push!(WARMUP_FAILED_SOLVERS, solver_name)
            elseif BENCHMARK_WORKER[] in workers()
                mark_solver_warmed!(solver_name)
            end
        else
            WARMED_ON_WORKER[] = -1
            push!(WARMED_SOLVERS, solver_name)
        end
    end
    return DataFrame(rows)
end

# ---------------------------------------------------------------------------
# Hard wall-clock guard
#
# Cooperative `maxtime` cannot interrupt a stalled callback or hung SIF decode.
# With `hard_timeout = true`, solves run on one long-lived Distributed worker per
# page; a row that exceeds `hard_timeout_seconds(maxtime)` kills that worker's
# process group, records TIMEOUT/HARD_TIMEOUT, and continues on a fresh worker.
# ---------------------------------------------------------------------------

const HARNESS_FILE = abspath(@__FILE__)
const BENCHMARK_WORKER = Ref{Int}(0)
# Worker PID after setpgid(0,0); used for process-group cleanup (public getpid).
const BENCHMARK_WORKER_PGID = Ref{Int}(0)
# Solvers successfully warmed on WARMED_ON_WORKER; cleared when the worker is replaced.
const WARMED_ON_WORKER = Ref{Int}(0)
const WARMED_SOLVERS = Set{String}()
# Solvers whose warm-up hung once; do not retry warm-up (avoids restart loops).
const WARMUP_FAILED_SOLVERS = Set{String}()

hard_timeout_seconds(maxtime = SOLVE_TIMEOUT_SECONDS) =
    2 * Float64(maxtime) + HARD_TIMEOUT_GRACE_SECONDS

function reset_warmup_state!()
    empty!(WARMED_SOLVERS)
    empty!(WARMUP_FAILED_SOLVERS)
    WARMED_ON_WORKER[] = 0
    return nothing
end

function sync_warmup_worker!()
    worker = BENCHMARK_WORKER[]
    if worker == 0 || !(worker in workers())
        # Drop marks from a dead worker; keep in-process marks (WARMED_ON_WORKER == -1).
        if WARMED_ON_WORKER[] > 0
            empty!(WARMED_SOLVERS)
            WARMED_ON_WORKER[] = 0
        end
        return nothing
    end
    if WARMED_ON_WORKER[] != worker
        empty!(WARMED_SOLVERS)
        WARMED_ON_WORKER[] = worker
    end
    return nothing
end

solver_warmed_here(solver_name) =
    solver_name in WARMED_SOLVERS && (
    WARMED_ON_WORKER[] == -1 || WARMED_ON_WORKER[] == BENCHMARK_WORKER[]
)

function mark_solver_warmed!(solver_name)
    sync_warmup_worker!()
    BENCHMARK_WORKER[] in workers() || return nothing
    WARMED_ON_WORKER[] = BENCHMARK_WORKER[]
    push!(WARMED_SOLVERS, solver_name)
    return nothing
end

function ensure_solver_warmed!(
        solver_name; problem = nothing, maxtime = SOLVE_TIMEOUT_SECONDS, hard_timeout = true
    )
    sync_warmup_worker!()
    solver_warmed_here(solver_name) && return nothing
    solver_name in WARMUP_FAILED_SOLVERS && return nothing
    warmup_solvers([solver_name]; problem, maxtime, hard_timeout)
    return nothing
end

function proc_pgrp(pid::Integer)
    stat_path = "/proc/$pid/stat"
    isfile(stat_path) || return -1
    s = read(stat_path, String)
    rparen = findlast(')', s)
    rparen === nothing && return -1
    fields = split(SubString(s, nextind(s, rparen)))
    length(fields) >= 3 || return -1
    return parse(Int, fields[3])
end

function process_group_alive(pgid::Integer)
    pgid <= 1 && return false
    isdir("/proc") || return false
    for name in readdir("/proc")
        all(isdigit, name) || continue
        try
            proc_pgrp(parse(Int, name)) == pgid && return true
        catch
        end
    end
    return false
end

function wait_process_group_exit(pgid::Integer; waitfor = 10.0)
    deadline = time() + waitfor
    while time() < deadline
        process_group_alive(pgid) || return true
        sleep(0.05)
    end
    return !process_group_alive(pgid)
end

function kill_process_group(pgid::Integer)
    pgid <= 1 && return false
    return ccall(:kill, Cint, (Cint, Cint), -Int32(pgid), Int32(9)) == 0
end

function configure_worker_process_group(worker)
    pgid = remotecall_fetch(worker) do
        ccall(:setpgid, Cint, (Cint, Cint), 0, 0)
        return Int(getpid())
    end
    BENCHMARK_WORKER_PGID[] = pgid
    return pgid
end

function start_benchmark_worker()
    started = time()
    blas_threads = BLAS.get_num_threads()
    worker = only(
        addprocs(
            1;
            exeflags = "--project=$(Base.active_project())",
            enable_threaded_blas = blas_threads > 1,
        )
    )
    remotecall_fetch(
        Core.eval, worker, Main, quote
            using CUTEst
            using LinearAlgebra
            BLAS.set_num_threads($blas_threads)
            include($HARNESS_FILE)
        end
    )
    BENCHMARK_WORKER[] = worker
    empty!(WARMED_SOLVERS)
    WARMED_ON_WORKER[] = worker
    configure_worker_process_group(worker)
    @printf("Benchmark worker %d ready (startup %.1fs)\n", worker, time() - started)
    return worker
end

function ensure_benchmark_worker()
    worker = BENCHMARK_WORKER[]
    worker in workers() && return worker
    return start_benchmark_worker()
end

function kill_benchmark_worker(worker)
    pgid = BENCHMARK_WORKER_PGID[]

    try
        rmprocs(worker; waitfor = 10)
    catch
    end

    if pgid > 1
        kill_process_group(pgid)
        wait_process_group_exit(pgid; waitfor = 10) || kill_process_group(pgid)
        wait_process_group_exit(pgid; waitfor = 5)
    end

    try
        worker in workers() && rmprocs(worker; waitfor = 5)
    catch err
        print_exception("Unable to remove benchmark worker $worker", err, catch_backtrace())
    end

    BENCHMARK_WORKER[] = 0
    BENCHMARK_WORKER_PGID[] = 0
    empty!(WARMED_SOLVERS)
    WARMED_ON_WORKER[] = 0
    return nothing
end

function synthetic_solve_row(
        problem_name, solver_name; secs, retcode, status, first_retcode = retcode
    )
    return (;
        problem = problem_name,
        solver = solver_name,
        solver_class = solver_class(solver_name),
        n_vars = -1,
        secs = secs,
        first_solve_secs = NaN,
        compile_secs = NaN,
        decode_secs = NaN,
        solver_reported_secs = NaN,
        retcode = retcode,
        first_retcode = first_retcode,
        status = status,
        NO_QUALITY...,
    )
end

function run_single_solve_guarded(
        problem_name, solver_name; maxtime = SOLVE_TIMEOUT_SECONDS
    )
    worker = ensure_benchmark_worker()
    deadline = hard_timeout_seconds(maxtime)
    started = time()

    waiter = @async remotecall_fetch(
        Core.eval, worker, Main,
        :(run_single_solve($problem_name, $solver_name; maxtime = $maxtime))
    )

    if timedwait(() -> istaskdone(waiter), deadline; pollint = 0.5) === :timed_out
        println(
            " HARD TIMEOUT: no result after $(deadline)s, killing worker $worker " *
                "(a fresh worker is started for the next solve)"
        )
        kill_benchmark_worker(worker)

        return synthetic_solve_row(
            problem_name, solver_name;
            secs = deadline, retcode = "HARD_TIMEOUT", status = "TIMEOUT",
        )
    end

    try
        return fetch(waiter)
    catch err
        bt = catch_backtrace()
        print_exception(
            "Benchmark worker $worker failed on $problem_name with $solver_name", err, bt
        )
        if !(worker in workers()) || BENCHMARK_WORKER[] == worker
            kill_benchmark_worker(worker)
        end

        return synthetic_solve_row(
            problem_name, solver_name;
            secs = time() - started, retcode = "WORKER_FAILED", status = "FAILED",
        )
    end
end

function run_benchmarks(
        category, problems, solvers;
        warmup = true, maxtime = SOLVE_TIMEOUT_SECONDS, hard_timeout = true,
    )
    rows = NamedTuple[]
    reset_warmup_state!()

    if hard_timeout
        println("Hard timeout: ", hard_timeout_seconds(maxtime), "s per solve (worker process)")
        ensure_benchmark_worker()
    end

    warmup && warmup_solvers(solvers; maxtime, hard_timeout)

    println()
    println("Running $category benchmarks")
    println("Problems: ", length(problems))
    println("Solvers: ", join(solvers, ", "))
    println("Per-solve time cap: ", maxtime, " s")

    for problem_name in problems
        for solver_name in solvers
            if hard_timeout
                ensure_benchmark_worker()
                sync_warmup_worker!()
            end
            if warmup
                ensure_solver_warmed!(solver_name; maxtime, hard_timeout)
            end
            @printf("  %-18s %-24s", solver_name, problem_name)
            row = if hard_timeout
                run_single_solve_guarded(problem_name, solver_name; maxtime)
            else
                run_single_solve(problem_name, solver_name; maxtime)
            end
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
