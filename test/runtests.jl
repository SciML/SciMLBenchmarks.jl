using SciMLTesting

run_tests(
    core = joinpath(@__DIR__, "core.jl"),
    qa = (
        env = joinpath(@__DIR__, "qa"),
        body = joinpath(@__DIR__, "qa", "qa.jl"),
    ),
)
