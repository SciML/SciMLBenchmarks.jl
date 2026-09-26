#=
Worker process for the Turing.jl AD benchmark (`TuringADTests.jmd`).

    julia --project=. turing_ad_worker.jl <model_name> <backend> [<backend>...]

Loads a single Turing model and runs `DynamicPPL.TestUtils.AD.run_ad` on it once
per requested backend, printing one tab-separated record per backend to stdout.

Each backend is run in a subprocess rather than in the parent because AD on some
of these models can take the whole process down (Enzyme in particular segfaults
sporadically), and a crash must not lose the results already collected. The
parent restarts the worker with whatever is left. This mirrors the Python/Julia
split used upstream in https://github.com/TuringLang/ADTests.

Protocol (all lines tab-separated, flushed as they are produced):

    META    <category>       <dimension>
    BEGIN   <backend>
    RESULT  <backend>        <status>  <relative>  <gradient_seconds>  <primal_seconds>
    MSG     <backend>        <one-line error description>

`status` is one of `ok`, `wrong`, `NaN`, or `error`; the numeric fields are
`NaN` unless the status is `ok`. A backend with a `BEGIN` but no `RESULT` was
in flight when the process died, and is recorded as a crash by the parent.
=#

using ADTypes:
    AutoEnzyme, AutoFiniteDifferences, AutoForwardDiff, AutoMooncake,
    AutoMooncakeForward, AutoReverseDiff
using DynamicPPL: DynamicPPL, LogDensityFunction, getlogjoint_internal, LinkAll
using DynamicPPL.TestUtils.AD: run_ad, ADIncorrectException, WithBackend
using LogDensityProblems: LogDensityProblems
using Random: Random, Xoshiro

import FiniteDifferences: central_fdm
import ForwardDiff
import ReverseDiff
import Mooncake
import Enzyme: set_runtime_activity, Forward, Reverse, Const

include(joinpath(@__DIR__, "turing_ad_models.jl"))

const ADTYPES = Dict(
    "FiniteDifferences" => AutoFiniteDifferences(; fdm = central_fdm(5, 1)),
    "ForwardDiff" => AutoForwardDiff(),
    "ReverseDiff" => AutoReverseDiff(; compile = false),
    "ReverseDiffCompiled" => AutoReverseDiff(; compile = true),
    "MooncakeRvs" => AutoMooncake(),
    "MooncakeFwd" => AutoMooncakeForward(),
    "EnzymeFwd" => AutoEnzyme(;
        mode = set_runtime_activity(Forward, true),
        function_annotation = Const,
    ),
    "EnzymeRvs" => AutoEnzyme(;
        mode = set_runtime_activity(Reverse, true),
        function_annotation = Const,
    ),
)

"""
    load_model(model_name)

Evaluate `models/\$(model_name).jl` inside a fresh anonymous module and return
the `model` object it defines.

Each model lives in its own module so that models can declare data and helper
functions without worrying about name clashes, and so that the model files
themselves stay free of the boilerplate (`using Turing`, `module ... end`) that
would obscure the model definition for a reader.
"""
function load_model(model_name::AbstractString)
    path = joinpath(@__DIR__, "models", model_name * ".jl")
    isfile(path) || error("no model file at $(path)")
    # Several models synthesise their data at load time from the global RNG.
    # Left unseeded, the reported status of a handful of models flips between
    # runs: `control_flow` takes a different branch, and `dppl_hmm_semisup`
    # draws a state sequence that some backends can and some cannot handle.
    Random.seed!(468)
    modname = gensym(model_name)
    # `module` is only valid at top level, hence the `Expr(:toplevel, ...)`.
    Core.eval(
        Main,
        Expr(
            :toplevel,
            :(
                module $modname
                using Turing
                include($path)
                end
            ),
        ),
    )
    mod = Base.invokelatest(getfield, Main, modname)
    return Base.invokelatest(getfield, mod, :model)
end

# Some models are more numerically sensitive than the default `sqrt(eps())`
# tolerance allows for; the reference gradient itself is only a finite
# difference approximation.
function relative_tolerance(model_name)
    return if model_name == "dppl_logistic_regression"
        1.0e-1
    elseif model_name == "lux_nn"
        1.0e-2
    elseif model_name == "ordinarydiffeq" || model_name == "delaydiffeq"
        1.0e-3
    else
        sqrt(eps())
    end
end

# FiniteDifferences errors on `dppl_hmm_semisup`, which would mark every backend
# on that model as wrong. See https://github.com/TuringLang/ADTests/issues/40.
function reference_backend(model_name)
    return model_name == "dppl_hmm_semisup" ? ADTYPES["ForwardDiff"] :
        ADTYPES[REFERENCE_BACKEND]
end

emit(args...) = (println(join(args, '\t')); flush(stdout))

function benchmark_backends(model, model_name, backend_names)
    ldf = LogDensityFunction(model, getlogjoint_internal, LinkAll())
    emit("META", CATEGORY_OF[model_name], LogDensityProblems.dimension(ldf))

    rtol = relative_tolerance(model_name)
    ref = reference_backend(model_name)

    for backend_name in backend_names
        emit("BEGIN", backend_name)
        try
            result = run_ad(
                model,
                ADTYPES[backend_name];
                rng = Xoshiro(468),
                test = WithBackend(ref),
                benchmark = true,
                rtol = rtol,
                verbose = false,
            )
            emit(
                "RESULT", backend_name, "ok",
                result.grad_time / result.primal_time,
                result.grad_time, result.primal_time,
            )
        catch e
            showerror(stderr, e)
            println(stderr)
            status = if e isa ADIncorrectException
                any(isnan, e.grad_expected) || any(isnan, e.grad_actual) ? "NaN" : "wrong"
            else
                "error"
            end
            # Tabs would break the record format; the first line of an error is
            # the informative part and the rest is often hundreds of lines.
            description = replace(first(split(sprint(showerror, e), '\n')), '\t' => ' ')
            emit("MSG", backend_name, description)
            emit("RESULT", backend_name, status, NaN, NaN, NaN)
        end
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error(
        "usage: julia --project=. turing_ad_worker.jl <model_name> <backend> [<backend>...]"
    )
    model = load_model(ARGS[1])
    # The model's methods were defined by the `eval` above, so they are newer
    # than any world age visible from here without `invokelatest`.
    Base.invokelatest(benchmark_backends, model, ARGS[1], ARGS[2:end])
end
