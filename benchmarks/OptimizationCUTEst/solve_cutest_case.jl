using Optimization
using OptimizationNLPModels
using CUTEst
using OptimizationOptimJL
using OptimizationOptimJL: LBFGS, ConjugateGradient, NelderMead
using OptimizationMOI
using OptimizationMOI: MOI
using Ipopt

function optimizer_from_name(name)
    if name == "LBFGS"
        return LBFGS()
    elseif name == "ConjugateGradient"
        return ConjugateGradient()
    elseif name == "NelderMead"
        return NelderMead()
    elseif name == "Ipopt"
        return MOI.OptimizerWithAttributes(
            Ipopt.Optimizer,
            "max_iter" => 1000,
            "tol" => 1.0e-6,
            "print_level" => 0,
        )
    else
        error("Unknown optimizer: $name")
    end
end

function solve_seconds(sol, fallback)
    try
        if hasfield(typeof(sol), :stats) && hasfield(typeof(sol.stats), :time)
            return Float64(sol.stats.time)
        end
    catch
    end
    return fallback
end

function print_result(problem, solver, n_vars, secs, retcode, status)
    return println("CUTEST_RESULT\t", join((problem, solver, n_vars, secs, retcode, status), '\t'))
end

function main()
    problem_name = ARGS[1]
    optimizer_name = ARGS[2]

    nlp = nothing
    started = time()

    return try
        nlp = CUTEstModel(problem_name)
        prob = OptimizationNLPModels.OptimizationProblem(nlp, Optimization.AutoFiniteDiff())
        sol = solve(prob, optimizer_from_name(optimizer_name); maxiters = 1000, maxtime = 30.0)

        print_result(
            problem_name,
            optimizer_name,
            length(sol.u),
            solve_seconds(sol, time() - started),
            string(sol.retcode),
            "OK",
        )
    catch
        nvar = nlp === nothing ? -1 : nlp.meta.nvar
        print_result(problem_name, optimizer_name, nvar, time() - started, "FAILED", "FAILED")
    finally
        if nlp !== nothing
            try
                finalize(nlp)
            catch
            end
        end
    end
end

main()
