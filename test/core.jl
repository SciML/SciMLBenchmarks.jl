using SciMLBenchmarks
using Test

@testset "weave_file" begin
    benchmarks_dir = joinpath(dirname(@__DIR__), "benchmarks")
    SciMLBenchmarks.weave_file(joinpath(benchmarks_dir, "Testing"), "test.jmd")

    @test isfile(joinpath(dirname(@__DIR__), "markdown", "Testing", "test.md"))
end
