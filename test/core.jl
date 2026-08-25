using SciMLBenchmarks
using Test

@testset "IJulia extension" begin
    fallback_method = which(SciMLBenchmarks.open_notebooks, Tuple{})
    @test fallback_method.module === SciMLBenchmarks

    @eval using IJulia
    extension = Base.get_extension(SciMLBenchmarks, :SciMLBenchmarksIJuliaExt)
    @test extension !== nothing
    @test fallback_method in methods(SciMLBenchmarks.open_notebooks)
    @test which(SciMLBenchmarks.open_notebooks, Tuple{}).module === extension
end

@testset "weave_file" begin
    benchmarks_dir = joinpath(dirname(@__DIR__), "benchmarks")
    SciMLBenchmarks.weave_file(joinpath(benchmarks_dir, "Testing"), "test.jmd")

    @test isfile(joinpath(dirname(@__DIR__), "markdown", "Testing", "test.md"))
end
