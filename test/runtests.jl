using SciMLBenchmarks, Test

@testset "weave_file" begin
    benchmarks_dir = joinpath(dirname(@__DIR__), "benchmarks")
    SciMLBenchmarks.weave_file(benchmarks_dir, "test.jmd")

    @test isfile(joinpath(dirname(@__DIR__), "script", "benchmarks", "test.jl"))
    @test isfile(joinpath(dirname(@__DIR__), "html", "benchmarks", "test.html"))
    @test isfile(joinpath(dirname(@__DIR__), "markdown", "benchmarks", "test.md"))
    @test isfile(joinpath(dirname(@__DIR__), "pdf", "benchmarks", "test.pdf"))
end
